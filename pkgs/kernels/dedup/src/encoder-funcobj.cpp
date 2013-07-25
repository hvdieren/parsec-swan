/*
 * Decoder for dedup files
 *
 * Copyright 2010 Princeton University.
 * All rights reserved.
 *
 * Originally written by Minlan Yu.
 * Largely rewritten by Christian Bienia.
 */

/*
 * The pipeline model for Encode is Fragment->FragmentRefine->Deduplicate->Compress->Reorder
 * Each stage has basically three steps:
 * 1. fetch a group of items from the queue
 * 2. process the items
 * 3. put them in the queue for the next stage
 */

#include <assert.h>
#include <strings.h>
#include <math.h>
#include <limits.h>
#include <fcntl.h>
#include <errno.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <fcntl.h>
#include <string.h>

#include <list>

#include "tbb/pipeline.h"
#include "tbb/task_scheduler_init.h"

#define OUTER_PIPELINE_NUM_TOKENS 1024 
#define INNER_PIPELINE_NUM_TOKENS 1024

extern "C" {
#include "util.h"
#include "dedupdef.h"
#include "debug.h"
#include "hashtable.h"
#include "config.h"
#include "rabin.h"
#include "mbuffer.h"
}

#include "encoder.h"

#ifdef ENABLE_GZIP_COMPRESSION
#include <zlib.h>
#endif //ENABLE_GZIP_COMPRESSION

#ifdef ENABLE_BZIP2_COMPRESSION
#include <bzlib.h>
#endif //ENABLE_BZIP2_COMPRESSION

#ifdef ENABLE_PARSEC_HOOKS
#include <hooks.h>
#endif //ENABLE_PARSEC_HOOKS


#define INITIAL_SEARCH_TREE_SIZE 4096


//The configuration block defined in main
extern config_t * conf;

//Hash table data structure & utility functions
struct hashtable *cache;

static unsigned int hash_from_key_fn( void *k ) {
  //NOTE: sha1 sum is integer-aligned
  return ((unsigned int *)k)[0];
}

static int keys_equal_fn ( void *key1, void *key2 ) {
  return (memcmp(key1, key2, SHA1_LEN) == 0);
}

//Arguments to pass to each thread
struct thread_args {
  //thread id, unique within a thread pool (i.e. unique for a pipeline stage)
  int tid;
  //number of queues available, first and last pipeline stage only
  int nqueues;
  //file descriptor, first pipeline stage only
  int fd;
  //input file buffer, first pipeline stage & preloading only
  struct {
    void *buffer;
    size_t size;
  } input_file;
};


#ifdef ENABLE_STATISTICS
//Keep track of block granularity with 2^CHUNK_GRANULARITY_POW resolution (for statistics)
#define CHUNK_GRANULARITY_POW (7)
//Number of blocks to distinguish, CHUNK_MAX_NUM * 2^CHUNK_GRANULARITY_POW is biggest block being recognized (for statistics)
#define CHUNK_MAX_NUM (8*32)
//Map a chunk size to a statistics array slot
#define CHUNK_SIZE_TO_SLOT(s) ( ((s)>>(CHUNK_GRANULARITY_POW)) >= (CHUNK_MAX_NUM) ? (CHUNK_MAX_NUM)-1 : ((s)>>(CHUNK_GRANULARITY_POW)) )
//Get the average size of a chunk from a statistics array slot
#define SLOT_TO_CHUNK_SIZE(s) ( (s)*(1<<(CHUNK_GRANULARITY_POW)) + (1<<((CHUNK_GRANULARITY_POW)-1)) )
//Deduplication statistics (only used if ENABLE_STATISTICS is defined)
typedef struct {
  /* Cumulative sizes */
  size_t total_input; //Total size of input in bytes
  size_t total_dedup; //Total size of input without duplicate blocks (after global compression) in bytes
  size_t total_compressed; //Total size of input stream after local compression in bytes
  size_t total_output; //Total size of output in bytes (with overhead) in bytes

  /* Size distribution & other properties */
  unsigned int nChunks[CHUNK_MAX_NUM]; //Coarse-granular size distribution of data chunks
  unsigned int nDuplicates; //Total number of duplicate blocks
} stats_t;

//Initialize a statistics record
static void init_stats(stats_t *s) {
  int i;

  assert(s!=NULL);
  s->total_input = 0;
  s->total_dedup = 0;
  s->total_compressed = 0;
  s->total_output = 0;

  for(i=0; i<CHUNK_MAX_NUM; i++) {
    s->nChunks[i] = 0;
  }
  s->nDuplicates = 0;
}

//Print statistics
static void print_stats(stats_t *s) {
  const unsigned int unit_str_size = 7; //elements in unit_str array
  const char *unit_str[] = {"Bytes", "KB", "MB", "GB", "TB", "PB", "EB"};
  unsigned int unit_idx = 0;
  size_t unit_div = 1;

  assert(s!=NULL);

  //determine most suitable unit to use
  for(unit_idx=0; unit_idx<unit_str_size; unit_idx++) {
    unsigned int unit_div_next = unit_div * 1024;

    if(s->total_input / unit_div_next <= 0) break;
    if(s->total_dedup / unit_div_next <= 0) break;
    if(s->total_compressed / unit_div_next <= 0) break;
    if(s->total_output / unit_div_next <= 0) break;

    unit_div = unit_div_next;
  }

  printf("Total input size:              %14.2f %s\n", (float)(s->total_input)/(float)(unit_div), unit_str[unit_idx]);
  printf("Total output size:             %14.2f %s\n", (float)(s->total_output)/(float)(unit_div), unit_str[unit_idx]);
  printf("Effective compression factor:  %14.2fx\n", (float)(s->total_input)/(float)(s->total_output));
  printf("\n");

  //Total number of chunks
  unsigned int i;
  unsigned int nTotalChunks=0;
  for(i=0; i<CHUNK_MAX_NUM; i++) nTotalChunks+= s->nChunks[i];

  //Average size of chunks
  float mean_size = 0.0;
  for(i=0; i<CHUNK_MAX_NUM; i++) mean_size += (float)(SLOT_TO_CHUNK_SIZE(i)) * (float)(s->nChunks[i]);
  mean_size = mean_size / (float)nTotalChunks;

  //Variance of chunk size
  float var_size = 0.0;
  for(i=0; i<CHUNK_MAX_NUM; i++) var_size += (mean_size - (float)(SLOT_TO_CHUNK_SIZE(i))) *
                                             (mean_size - (float)(SLOT_TO_CHUNK_SIZE(i))) *
                                             (float)(s->nChunks[i]);

  printf("Mean data chunk size:          %14.2f %s (stddev: %.2f %s)\n", mean_size / 1024.0, "KB", sqrtf(var_size) / 1024.0, "KB");
  printf("Amount of duplicate chunks:    %14.2f%%\n", 100.0*(float)(s->nDuplicates)/(float)(nTotalChunks));
  printf("Data size after deduplication: %14.2f %s (compression factor: %.2fx)\n", (float)(s->total_dedup)/(float)(unit_div), unit_str[unit_idx], (float)(s->total_input)/(float)(s->total_dedup));
  printf("Data size after compression:   %14.2f %s (compression factor: %.2fx)\n", (float)(s->total_compressed)/(float)(unit_div), unit_str[unit_idx], (float)(s->total_dedup)/(float)(s->total_compressed));
  printf("Output overhead:               %14.2f%%\n", 100.0*(float)(s->total_output-s->total_compressed)/(float)(s->total_output));
}

//variable with global statistics
stats_t stats;
#endif //ENABLE_STATISTICS


//Simple write utility function
static int write_file(int fd, u_char type, u_long len, u_char * content) {
  if (xwrite(fd, &type, sizeof(type)) < 0){
    perror("xwrite:");
    EXIT_TRACE("xwrite type fails\n");
    return -1;
  }
  if (xwrite(fd, &len, sizeof(len)) < 0){
    EXIT_TRACE("xwrite content fails\n");
  }
  if (xwrite(fd, content, len) < 0){
    EXIT_TRACE("xwrite content fails\n");
  }
  return 0;
}

/*
 * Helper function that creates and initializes the output file
 * Takes the file name to use as input and returns the file handle
 * The output file can be used to write chunks without any further steps
 */
static int create_output_file(char *outfile) {
  int fd;

  //Create output file
  fd = open(outfile, O_CREAT|O_TRUNC|O_WRONLY|O_TRUNC, S_IRGRP | S_IWUSR | S_IRUSR | S_IROTH);
  if (fd < 0) {
    EXIT_TRACE("Cannot open output file.");
  }

  //Write header
  if (write_header(fd, conf->compress_type)) {
    EXIT_TRACE("Cannot write output file header.\n");
  }

  return fd;
}



/*
 * Helper function that writes a chunk to an output file depending on
 * its state. The function will write the SHA1 sum if the chunk has
 * already been written before, or it will write the compressed data
 * of the chunk if it has not been written yet.
 *
 * This function will block if the compressed data is not available yet.
 * This function might update the state of the chunk if there are any changes.
 */
//NOTE: The serial version relies on the fact that chunks are processed in-order,
//      which means if it reaches the function it is guaranteed all data is ready.
static void write_chunk_to_file(int fd, chunk_t *chunk) {
  assert(chunk!=NULL);

  if(!chunk->header.isDuplicate) {
    //Unique chunk, data has not been written yet, do so now
    write_file(fd, TYPE_COMPRESS, chunk->compressed_data.n, (uint8_t*)chunk->compressed_data.ptr);
    mbuffer_free(&chunk->compressed_data);
  } else {
    //Duplicate chunk, data has been written to file before, just write SHA1
    write_file(fd, TYPE_FINGERPRINT, SHA1_LEN, (unsigned char *)(chunk->sha1));
  }

  if(chunk->header.isDuplicate)
      free((void*)chunk);
} 

int rf_win;
int rf_win_dataprocess;

/*
 * Computational kernel of compression stage
 *
 * Actions performed:
 *  - Compress a data chunk
 */
void sub_Compress(chunk_t *chunk) {
    size_t n;
    int r;

    assert(chunk!=NULL);

    // Stage is empty if this is a duplicate.
    if(chunk->header.isDuplicate)
        return;

    //compress the item and add it to the database
    switch (conf->compress_type) {
      case COMPRESS_NONE:
        //Simply duplicate the data
        n = chunk->uncompressed_data.n;
        r = mbuffer_create(&chunk->compressed_data, n);
        if(r != 0) {
          EXIT_TRACE("Creation of compression buffer failed.\n");
        }
        //copy the block
        memcpy(chunk->compressed_data.ptr, chunk->uncompressed_data.ptr, chunk->uncompressed_data.n);
        break;
#ifdef ENABLE_GZIP_COMPRESSION
      case COMPRESS_GZIP:
        //Gzip compression buffer must be at least 0.1% larger than source buffer plus 12 bytes
        n = chunk->uncompressed_data.n + (chunk->uncompressed_data.n >> 9) + 12;
        r = mbuffer_create(&chunk->compressed_data, n);
        if(r != 0) {
          EXIT_TRACE("Creation of compression buffer failed.\n");
        }
        //compress the block
        r = compress( (Bytef*)chunk->compressed_data.ptr, (uLongf*)&n, (const Bytef*)chunk->uncompressed_data.ptr, (uLong)chunk->uncompressed_data.n);
        if (r != Z_OK) {
          EXIT_TRACE("Compression failed\n");
        }
        //Shrink buffer to actual size
        if(n < chunk->compressed_data.n) {
          r = mbuffer_realloc(&chunk->compressed_data, n);
          assert(r == 0);
        }
        break;
#endif //ENABLE_GZIP_COMPRESSION
#ifdef ENABLE_BZIP2_COMPRESSION
      case COMPRESS_BZIP2:
        //Bzip compression buffer must be at least 1% larger than source buffer plus 600 bytes
        n = chunk->uncompressed_data.n + (chunk->uncompressed_data.n >> 6) + 600;
        r = mbuffer_create(&chunk->compressed_data, n);
        if(r != 0) {
          EXIT_TRACE("Creation of compression buffer failed.\n");
        }
        //compress the block
        unsigned int int_n = n;
        r = BZ2_bzBuffToBuffCompress(chunk->compressed_data.ptr, &int_n, chunk->uncompressed_data.ptr, chunk->uncompressed_data.n, 9, 0, 30);
        n = int_n;
        if (r != BZ_OK) {
          EXIT_TRACE("Compression failed\n");
        }
        //Shrink buffer to actual size
        if(n < chunk->compressed_data.n) {
          r = mbuffer_realloc(&chunk->compressed_data, n);
          assert(r == 0);
        }
        break;
#endif //ENABLE_BZIP2_COMPRESSION
      default:
        EXIT_TRACE("Compression type not implemented.\n");
        break;
    }
    mbuffer_free(&chunk->uncompressed_data);

#ifdef ENABLE_STATISTICS
	// TODO: reduce
        stats.total_compressed += chunk->compressed_data.n;
#endif //ENABLE_STATISTICS

     return;
}


/*
 * Computational kernel of deduplication stage
 *
 * Actions performed:
 *  - Calculate SHA1 signature for each incoming data chunk
 *  - Perform database lookup to determine chunk redundancy status
 *  - On miss add chunk to database
 *  - Returns chunk redundancy status
 */
void sub_Deduplicate(chunk_t *chunk) {
  int isDuplicate;
  chunk_t *entry;

  assert(chunk!=NULL);
  assert(chunk->uncompressed_data.ptr!=NULL);

  // TODO: split off SHA1_Digest as separate stage?
  SHA1_Digest( (const void*)chunk->uncompressed_data.ptr, chunk->uncompressed_data.n, (unsigned char *)(chunk->sha1));

  //Query database to determine whether we've seen the data chunk before
  entry = (chunk_t *)hashtable_search(cache, (void *)(chunk->sha1));
  isDuplicate = (entry != NULL);
  chunk->header.isDuplicate = isDuplicate;
  if (!isDuplicate) {
    // Cache miss: Create entry in hash table and forward data to compression stage
    //NOTE: chunk->compressed_data.buffer will be computed in compression stage
    if (hashtable_insert(cache, (void *)(chunk->sha1), (void *)chunk) == 0) {
      EXIT_TRACE("hashtable_insert failed");
    }
  } else {
    // Cache hit: Skipping compression stage
    chunk->compressed_data_ref = entry;
    mbuffer_free(&chunk->uncompressed_data);
  }

#ifdef ENABLE_STATISTICS
      // HV: TODO: use reducer on stats ...
      if(isDuplicate) {
        stats.nDuplicates++;
      } else {
        stats.total_dedup += chunk->uncompressed_data.n;
      }
#endif //ENABLE_STATISTICS

  return;
}

class InnerPipeline : public tbb::filter {
    class FragmentRefine : public tbb::filter {
	u32int rabintab[256];
	u32int rabinwintab[256];
	chunk_t * chunk;
	int split;
	static const int rf_win = 0;

    public:
	FragmentRefine( chunk_t * chunk_in ) : filter( true ), chunk( chunk_in ), split( 1 ) {
	    rabininit(rf_win, rabintab, rabinwintab);
	}

	void * operator() ( void * ) { // no inputs
	    chunk_t * temp;

	    if( !split )
		return NULL;

	    //partition input block into fine-granular chunks
	    split = 0;
	    //Try to split the buffer at least ANCHOR_JUMP bytes away from its beginning
	    int offset = rabinseg((uchar*)chunk->uncompressed_data.ptr, (int)chunk->uncompressed_data.n, rf_win, rabintab, rabinwintab);
	    //Did we find a split location?
	    if(offset < (int)chunk->uncompressed_data.n) {
		//Split found somewhere in the middle of the buffer
		//Allocate a new chunk and create a new memory buffer
		temp = (chunk_t *)malloc(sizeof(chunk_t));
		if(temp==NULL) EXIT_TRACE("Memory allocation failed.\n");

		//split it into two pieces
		int r = mbuffer_split(&chunk->uncompressed_data, &temp->uncompressed_data, (size_t)offset);
		if(r!=0) EXIT_TRACE("Unable to split memory buffer.\n");
		temp->header.state = CHUNK_STATE_UNCOMPRESSED;

#ifdef ENABLE_STATISTICS
		//update statistics
		stats.nChunks[CHUNK_SIZE_TO_SLOT(chunk->uncompressed_data.n)]++;
#endif //ENABLE_STATISTICS

		split = 1;
	    } else {
		// Quit loop. Last piece of fragment.
		split = 0;
	    }

	    chunk_t * ret = chunk;

	    //prepare for next iteration
	    chunk = temp;
	    temp = NULL;
	    return (void *)ret;
	}
    };

    class Deduplicate : public tbb::filter {
    public:
	Deduplicate() : filter( false ) { }
	void * operator() ( void * item ) {
	    chunk_t * chunk = (chunk_t *)item;
	    sub_Deduplicate( chunk );
	    return item;
	}
    };

    class Compress : public tbb::filter {
    public:
	Compress() : filter( false ) { }
	void * operator() ( void * item ) {
	    chunk_t * chunk = (chunk_t *)item;
	    sub_Compress( chunk );
	    return item;
	}
    };

    class Reassemble : public tbb::filter {
	std::list<chunk_t *> & list_out;

    public:
	Reassemble( std::list<chunk_t *> & l ) : filter( true ), list_out( l ) { }
	void * operator() ( void * item ) {
	    chunk_t * chunk = (chunk_t *)item;
	    list_out.push_back( chunk );
	    return NULL;
	}
    };

public:
    InnerPipeline() : filter( false ) { }

    void * operator() ( void * item ) {
	chunk_t * chunk = (chunk_t *)item;
	std::list<chunk_t *> * list_out = new std::list<chunk_t *>();

	tbb::pipeline inner;

	FragmentRefine refine( chunk );
	Deduplicate deduplicate;
	Compress compress;
	Reassemble reassemble( *list_out );

	inner.add_filter( refine );
	inner.add_filter( deduplicate );
	inner.add_filter( compress );
	inner.add_filter( reassemble );

	inner.run( INNER_PIPELINE_NUM_TOKENS );
	inner.clear();

	return (void *)list_out;
    }
};

class OuterPipeline {
    int fd_out;

    class Fragment : public tbb::filter {
	u32int rabintab[256];
	u32int rabinwintab[256];
	struct thread_args * args;
	int fd;
	int split;
	size_t preloading_buffer_seek;
	chunk_t *temp;
	chunk_t *chunk;
	static const int rf_win_dataprocess = 0;

    public:
	Fragment( struct thread_args * args ) : tbb::filter( true ), args( args ) {
	    rabininit(rf_win_dataprocess, rabintab, rabinwintab);
	    fd = args->fd;
	    split = 0;
	    preloading_buffer_seek = 0;
	    temp = NULL;
	    chunk = NULL;

	    //Sanity check
	    if(MAXBUF < 8 * ANCHOR_JUMP) {
		printf("WARNING: I/O buffer size is very small. Performance degraded.\n");
		fflush(NULL);
	    }
	}

	void * operator()( void * ) {
	    int r;

	    if( split < 0 )
		return NULL;

	    //read from input file / buffer
	outer_loop:
	    if( !split ) {
		size_t bytes_left; //amount of data left over in last_mbuffer from previous iteration

		//Check how much data left over from previous iteration resp. create an initial chunk
		if(temp != NULL) {
		    bytes_left = temp->uncompressed_data.n;
		} else {
		    bytes_left = 0;
		}

		//Make sure that system supports new buffer size
		if(MAXBUF+bytes_left > SSIZE_MAX) {
		    EXIT_TRACE("Input buffer size exceeds system maximum.\n");
		}
		//Allocate a new chunk and create a new memory buffer
		chunk = (chunk_t *)malloc(sizeof(chunk_t));
		if(chunk==NULL) EXIT_TRACE("Memory allocation failed.\n");
		r = mbuffer_create(&chunk->uncompressed_data, MAXBUF+bytes_left);
		if(r!=0) {
		    EXIT_TRACE("Unable to initialize memory buffer.\n");
		}
		chunk->header.state = CHUNK_STATE_UNCOMPRESSED;
		if(bytes_left > 0) {
		    //FIXME: Short-circuit this if no more data available

		    //"Extension" of existing buffer, copy sequence number and left over data to beginning of new buffer
		    //NOTE: We cannot safely extend the current memory region because it has already been given to another thread
		    memcpy(chunk->uncompressed_data.ptr, temp->uncompressed_data.ptr, temp->uncompressed_data.n);
		    mbuffer_free(&temp->uncompressed_data);
		    free((void*)temp);
		    temp = NULL;
		}
		//Read data until buffer full
		size_t bytes_read=0;
		if(conf->preloading) {
		    size_t max_read = MIN(MAXBUF, args->input_file.size-preloading_buffer_seek);
		    memcpy((char*)chunk->uncompressed_data.ptr+bytes_left, (char*)args->input_file.buffer+preloading_buffer_seek, max_read);
		    bytes_read = max_read;
		    preloading_buffer_seek += max_read;
		} else {
		    while(bytes_read < MAXBUF) {
			r = read( fd, (void*)((char*)chunk->uncompressed_data.ptr+bytes_left+bytes_read), (size_t)MAXBUF-bytes_read);
			if(r<0) switch(errno) {
			    case EAGAIN:
				EXIT_TRACE("I/O error: No data available\n");break;
			    case EBADF:
				EXIT_TRACE("I/O error: Invalid file descriptor\n");break;
			    case EFAULT:
				EXIT_TRACE("I/O error: Buffer out of range\n");break;
			    case EINTR:
				EXIT_TRACE("I/O error: Interruption\n");break;
			    case EINVAL:
				EXIT_TRACE("I/O error: Unable to read from file descriptor\n");break;
			    case EIO:
				EXIT_TRACE("I/O error: Generic I/O error\n");break;
			    case EISDIR:
				EXIT_TRACE("I/O error: Cannot read from a directory\n");break;
			    default:
				EXIT_TRACE("I/O error: Unrecognized error\n");break;
			    }
			if(r==0) break;
			bytes_read += r;
		    }
		}
		//No data left over from last iteration and also nothing new read in, simply clean up and quit
		if(bytes_left + bytes_read == 0) {
		    mbuffer_free( &chunk->uncompressed_data);
		    free( (void*)chunk);
		    chunk = NULL;
		    split = -1;
		    return NULL;
		}
		//Shrink buffer to actual size
		if(bytes_left+bytes_read < chunk->uncompressed_data.n) {
		    r = mbuffer_realloc( &chunk->uncompressed_data, bytes_left+bytes_read);
		    assert(r == 0);
		}

		//Check whether any new data was read in, process last chunk if not
		if(bytes_read == 0) {
#ifdef ENABLE_STATISTICS
		    //update statistics
		    stats.nChunks[CHUNK_SIZE_TO_SLOT(chunk->uncompressed_data.n)]++;
#endif //ENABLE_STATISTICS

		    return (void *)chunk;
		    split = -1;
		}
	    }

	    //partition input block into fine-granular chunks
	    split = 0;
	    //Try to split the buffer at least ANCHOR_JUMP bytes away from its beginning
	    if(ANCHOR_JUMP < chunk->uncompressed_data.n) {
		int offset = rabinseg((uchar*)chunk->uncompressed_data.ptr + ANCHOR_JUMP, (int)chunk->uncompressed_data.n - ANCHOR_JUMP, rf_win_dataprocess, rabintab, rabinwintab);
		//Did we find a split location?
		if(offset == 0) {
		    //Split found at the very beginning of the buffer (should never happen due to technical limitations)
		    assert(0);
		    split = 0;
		} else if(offset < (int)chunk->uncompressed_data.n) {
		    //Split found somewhere in the middle of the buffer
		    //Allocate a new chunk and create a new memory buffer
		    temp = (chunk_t *)malloc(sizeof(chunk_t));
		    if(temp==NULL) EXIT_TRACE("Memory allocation failed.\n");

		    //split it into two pieces
		    r = mbuffer_split(&chunk->uncompressed_data, &temp->uncompressed_data, (size_t)offset);
		    if(r!=0) EXIT_TRACE("Unable to split memory buffer.\n");
		    temp->header.state = CHUNK_STATE_UNCOMPRESSED;

#ifdef ENABLE_STATISTICS
		    //update statistics
		    stats.nChunks[CHUNK_SIZE_TO_SLOT(chunk->uncompressed_data.n)]++;
#endif //ENABLE_STATISTICS

		    //prepare for next iteration
		    chunk_t * ret = chunk;
		    chunk = temp;
		    split = 1;
		    return (void *)ret;
		} else {
		    //Due to technical limitations we can't distinguish the cases "no split" and "split at end of buffer"
		    //This will result in some unnecessary (and unlikely) work but yields the correct result eventually.
		    split = 0;
		}
	    } else {
		split = 0;
	    }

	    temp = (chunk_t*)chunk;
	    chunk = 0;

	    goto outer_loop;
	}
    };

    class Reorder : public tbb::filter {
	int fd_out;
    public:
	Reorder( int fd_out ) : tbb::filter( true ), fd_out( fd_out ) {
	}

	void * operator()( void * l ) {
	    std::list<chunk_t *> & ll = *(std::list<chunk_t *>*)l;
	    for( std::list<chunk_t *>::const_iterator I=ll.begin(), E=ll.end(); I != E; ++I ) {
		write_chunk_to_file( fd_out, *I );
	    }
	    delete &ll;
	    return NULL;
	}
    };
public:
    OuterPipeline( struct thread_args * args ) {
	fd_out = create_output_file(conf->outfile);

	tbb::pipeline outer;

	Fragment fragment( args );
	InnerPipeline inner;
	Reorder reorder( fd_out );

	outer.add_filter( fragment );
	outer.add_filter( inner );
	outer.add_filter( reorder );

	outer.run( OUTER_PIPELINE_NUM_TOKENS );
	outer.clear();

	close( (int)fd_out);
    }
};

#if 0
void TBBIntegratedPipeline(struct thread_args * args) {
  size_t preloading_buffer_seek = 0;
  int fd = args->fd;
  int fd_out = create_output_file(conf->outfile);
  int r;

  chunk_t *temp = NULL;
  chunk_t *chunk = NULL;
  u32int * rabintab = (uint32_t*)malloc(256*sizeof rabintab[0]);
  u32int * rabinwintab = (uint32_t*)malloc(256*sizeof rabintab[0]);
  if(rabintab == NULL || rabinwintab == NULL) {
    EXIT_TRACE("Memory allocation failed.\n");
  }

  rf_win_dataprocess = 0;
  rabininit(rf_win_dataprocess, rabintab, rabinwintab);

  //Sanity check
  if(MAXBUF < 8 * ANCHOR_JUMP) {
    printf("WARNING: I/O buffer size is very small. Performance degraded.\n");
    fflush(NULL);
  }

  //read from input file / buffer
  while (1) {
    size_t bytes_left; //amount of data left over in last_mbuffer from previous iteration

    //Check how much data left over from previous iteration resp. create an initial chunk
    if(temp != NULL) {
      bytes_left = temp->uncompressed_data.n;
    } else {
      bytes_left = 0;
    }

    //Make sure that system supports new buffer size
    if(MAXBUF+bytes_left > SSIZE_MAX) {
      EXIT_TRACE("Input buffer size exceeds system maximum.\n");
    }
    //Allocate a new chunk and create a new memory buffer
    chunk = (chunk_t *)malloc(sizeof(chunk_t));
    if(chunk==NULL) EXIT_TRACE("Memory allocation failed.\n");
    r = mbuffer_create(&chunk->uncompressed_data, MAXBUF+bytes_left);
    if(r!=0) {
      EXIT_TRACE("Unable to initialize memory buffer.\n");
    }
    chunk->header.state = CHUNK_STATE_UNCOMPRESSED;
    if(bytes_left > 0) {
      //FIXME: Short-circuit this if no more data available

      //"Extension" of existing buffer, copy sequence number and left over data to beginning of new buffer
      //NOTE: We cannot safely extend the current memory region because it has already been given to another thread
      memcpy(chunk->uncompressed_data.ptr, temp->uncompressed_data.ptr, temp->uncompressed_data.n);
      mbuffer_free(&temp->uncompressed_data);
      free((void*)temp);
      temp = NULL;
    }
    //Read data until buffer full
    size_t bytes_read=0;
    if(conf->preloading) {
      size_t max_read = MIN(MAXBUF, args->input_file.size-preloading_buffer_seek);
      memcpy((char*)chunk->uncompressed_data.ptr+bytes_left, (char*)args->input_file.buffer+preloading_buffer_seek, max_read);
      bytes_read = max_read;
      preloading_buffer_seek += max_read;
    } else {
      while(bytes_read < MAXBUF) {
        r = read( fd, (void*)((char*)chunk->uncompressed_data.ptr+bytes_left+bytes_read), (size_t)MAXBUF-bytes_read);
        if(r<0) switch(errno) {
          case EAGAIN:
            EXIT_TRACE("I/O error: No data available\n");break;
          case EBADF:
            EXIT_TRACE("I/O error: Invalid file descriptor\n");break;
          case EFAULT:
            EXIT_TRACE("I/O error: Buffer out of range\n");break;
          case EINTR:
            EXIT_TRACE("I/O error: Interruption\n");break;
          case EINVAL:
            EXIT_TRACE("I/O error: Unable to read from file descriptor\n");break;
          case EIO:
            EXIT_TRACE("I/O error: Generic I/O error\n");break;
          case EISDIR:
            EXIT_TRACE("I/O error: Cannot read from a directory\n");break;
          default:
            EXIT_TRACE("I/O error: Unrecognized error\n");break;
        }
        if(r==0) break;
        bytes_read += r;
      }
    }
    //No data left over from last iteration and also nothing new read in, simply clean up and quit
    if(bytes_left + bytes_read == 0) {
      mbuffer_free( &chunk->uncompressed_data);
      free( (void*)chunk);
      chunk = NULL;
      break;
    }
    //Shrink buffer to actual size
    if(bytes_left+bytes_read < chunk->uncompressed_data.n) {
      r = mbuffer_realloc( &chunk->uncompressed_data, bytes_left+bytes_read);
      assert(r == 0);
    }

    //Check whether any new data was read in, process last chunk if not
    if(bytes_read == 0) {
#ifdef ENABLE_STATISTICS
      //update statistics
      stats.nChunks[CHUNK_SIZE_TO_SLOT(chunk->uncompressed_data.n)]++;
#endif //ENABLE_STATISTICS

      // Pthreads code will send to FragmentRefine here...
      {
	  InnerPipeline inner( chunk );
	  std::list<chunk_t *> & ll = inner.get_list();
	  for( std::list<chunk_t *>::const_iterator I=ll.begin(), E=ll.end(); I != E; ++I ) {
	      write_chunk_to_file( fd_out, *I );
	  }
      }

      //stop fetching from input buffer, terminate processing
      break;
    }

    // "FragmentRefine"

    //partition input block into fine-granular chunks
    int split;
    do {
      split = 0;
      //Try to split the buffer at least ANCHOR_JUMP bytes away from its beginning
      if(ANCHOR_JUMP < chunk->uncompressed_data.n) {
	  int offset = rabinseg((uchar*)chunk->uncompressed_data.ptr + ANCHOR_JUMP, (int)chunk->uncompressed_data.n - ANCHOR_JUMP, rf_win_dataprocess, rabintab, rabinwintab);
      //Did we find a split location?
	  if(offset == 0) {
	    //Split found at the very beginning of the buffer (should never happen due to technical limitations)
	    assert(0);
	    split = 0;
	  } else if(offset < (int)chunk->uncompressed_data.n) {
	    //Split found somewhere in the middle of the buffer
	    //Allocate a new chunk and create a new memory buffer
	    temp = (chunk_t *)malloc(sizeof(chunk_t));
	    if(temp==NULL) EXIT_TRACE("Memory allocation failed.\n");

	    //split it into two pieces
	    r = mbuffer_split(&chunk->uncompressed_data, &temp->uncompressed_data, (size_t)offset);
	    if(r!=0) EXIT_TRACE("Unable to split memory buffer.\n");
	    temp->header.state = CHUNK_STATE_UNCOMPRESSED;

    #ifdef ENABLE_STATISTICS
	    //update statistics
	    stats.nChunks[CHUNK_SIZE_TO_SLOT(chunk->uncompressed_data.n)]++;
    #endif //ENABLE_STATISTICS

	    {
		InnerPipeline inner( chunk );
		std::list<chunk_t *> & ll = inner.get_list();
		for( std::list<chunk_t *>::const_iterator I=ll.begin(), E=ll.end(); I != E; ++I ) {
		    write_chunk_to_file( fd_out, *I );
		}
	    }

	    //prepare for next iteration
	    chunk = temp;
	    temp = NULL;
	    split = 1;
	  } else {
	    //Due to technical limitations we can't distinguish the cases "no split" and "split at end of buffer"
	    //This will result in some unnecessary (and unlikely) work but yields the correct result eventually.
	    split = 0;
	  }
      } else {
	  split = 0;
      }
    } while(split);
    temp = (chunk_t*)chunk;
    chunk = 0;
  }

  free( (void*)rabintab);
  free( (void*)rabinwintab);

  close( (int)fd_out);

  return;
}
#endif

void *SerialIntegratedPipeline(void * targs) {
  struct thread_args *args = (struct thread_args *)targs;
  size_t preloading_buffer_seek = 0;
  int fd = args->fd;
  int fd_out = create_output_file(conf->outfile);
  int r;

  chunk_t *temp = NULL;
  chunk_t *chunk = NULL;
  u32int * rabintab = (u32int*)malloc(256*sizeof rabintab[0]);
  u32int * rabinwintab = (u32int*)malloc(256*sizeof rabintab[0]);
  if(rabintab == NULL || rabinwintab == NULL) {
    EXIT_TRACE("Memory allocation failed.\n");
  }

  rf_win_dataprocess = 0;
  rabininit(rf_win_dataprocess, rabintab, rabinwintab);

  //Sanity check
  if(MAXBUF < 8 * ANCHOR_JUMP) {
    printf("WARNING: I/O buffer size is very small. Performance degraded.\n");
    fflush(NULL);
  }

  //read from input file / buffer
  while (1) {
    size_t bytes_left; //amount of data left over in last_mbuffer from previous iteration

    //Check how much data left over from previous iteration resp. create an initial chunk
    if(temp != NULL) {
      bytes_left = temp->uncompressed_data.n;
    } else {
      bytes_left = 0;
    }

    //Make sure that system supports new buffer size
    if(MAXBUF+bytes_left > SSIZE_MAX) {
      EXIT_TRACE("Input buffer size exceeds system maximum.\n");
    }
    //Allocate a new chunk and create a new memory buffer
    chunk = (chunk_t *)malloc(sizeof(chunk_t));
    if(chunk==NULL) EXIT_TRACE("Memory allocation failed.\n");
    r = mbuffer_create(&chunk->uncompressed_data, MAXBUF+bytes_left);
    if(r!=0) {
      EXIT_TRACE("Unable to initialize memory buffer.\n");
    }
    chunk->header.state = CHUNK_STATE_UNCOMPRESSED;
    if(bytes_left > 0) {
      //FIXME: Short-circuit this if no more data available

      //"Extension" of existing buffer, copy sequence number and left over data to beginning of new buffer
      //NOTE: We cannot safely extend the current memory region because it has already been given to another thread
      memcpy(chunk->uncompressed_data.ptr, temp->uncompressed_data.ptr, temp->uncompressed_data.n);
      mbuffer_free(&temp->uncompressed_data);
      free(temp);
      temp = NULL;
    }
    //Read data until buffer full
    size_t bytes_read=0;
    if(conf->preloading) {
      size_t max_read = MIN(MAXBUF, args->input_file.size-preloading_buffer_seek);
      memcpy((char*)chunk->uncompressed_data.ptr+bytes_left, (char*)args->input_file.buffer+preloading_buffer_seek, max_read);
      bytes_read = max_read;
      preloading_buffer_seek += max_read;
    } else {
      while(bytes_read < MAXBUF) {
        r = read(fd, (char*)chunk->uncompressed_data.ptr+bytes_left+bytes_read, (size_t)MAXBUF-bytes_read);
        if(r<0) switch(errno) {
          case EAGAIN:
            EXIT_TRACE("I/O error: No data available\n");break;
          case EBADF:
            EXIT_TRACE("I/O error: Invalid file descriptor\n");break;
          case EFAULT:
            EXIT_TRACE("I/O error: Buffer out of range\n");break;
          case EINTR:
            EXIT_TRACE("I/O error: Interruption\n");break;
          case EINVAL:
            EXIT_TRACE("I/O error: Unable to read from file descriptor\n");break;
          case EIO:
            EXIT_TRACE("I/O error: Generic I/O error\n");break;
          case EISDIR:
            EXIT_TRACE("I/O error: Cannot read from a directory\n");break;
          default:
            EXIT_TRACE("I/O error: Unrecognized error\n");break;
        }
        if(r==0) break;
        bytes_read += r;
      }
    }
    //No data left over from last iteration and also nothing new read in, simply clean up and quit
    if(bytes_left + bytes_read == 0) {
      mbuffer_free(&chunk->uncompressed_data);
      free(chunk);
      chunk = NULL;
      break;
    }
    //Shrink buffer to actual size
    if(bytes_left+bytes_read < chunk->uncompressed_data.n) {
      r = mbuffer_realloc(&chunk->uncompressed_data, bytes_left+bytes_read);
      assert(r == 0);
    }

    //Check whether any new data was read in, process last chunk if not
    if(bytes_read == 0) {
#ifdef ENABLE_STATISTICS
      //update statistics
      stats.nChunks[CHUNK_SIZE_TO_SLOT(chunk->uncompressed_data.n)]++;
#endif //ENABLE_STATISTICS

      //Deduplicate
      sub_Deduplicate(chunk);
      int isDuplicate = chunk->header.isDuplicate;
#ifdef ENABLE_STATISTICS
      if(isDuplicate) {
        stats.nDuplicates++;
      } else {
        stats.total_dedup += chunk->uncompressed_data.n;
      }
#endif //ENABLE_STATISTICS

      //If chunk is unique compress & archive it.
      if(!isDuplicate) {
        sub_Compress(chunk);
#ifdef ENABLE_STATISTICS
        stats.total_compressed += chunk->compressed_data.n;
#endif //ENABLE_STATISTICS
      }

      write_chunk_to_file(fd_out, chunk);

      //stop fetching from input buffer, terminate processing
      break;
    }

    //partition input block into fine-granular chunks
    int split;
    do {
      split = 0;
      //Try to split the buffer
      int offset = rabinseg( (uchar*)chunk->uncompressed_data.ptr, (int)chunk->uncompressed_data.n, rf_win_dataprocess, rabintab, rabinwintab);
      //Did we find a split location?
      if(offset == 0) {
        //Split found at the very beginning of the buffer (should never happen due to technical limitations)
        assert(0);
        split = 0;
      } else if(offset < (int)chunk->uncompressed_data.n) {
        //Split found somewhere in the middle of the buffer
        //Allocate a new chunk and create a new memory buffer
        temp = (chunk_t *)malloc(sizeof(chunk_t));
        if(temp==NULL) EXIT_TRACE("Memory allocation failed.\n");

        //split it into two pieces
        r = mbuffer_split(&chunk->uncompressed_data, &temp->uncompressed_data, offset);
        if(r!=0) EXIT_TRACE("Unable to split memory buffer.\n");
        temp->header.state = CHUNK_STATE_UNCOMPRESSED;

#ifdef ENABLE_STATISTICS
        //update statistics
        stats.nChunks[CHUNK_SIZE_TO_SLOT(chunk->uncompressed_data.n)]++;
#endif //ENABLE_STATISTICS

        //Deduplicate
	sub_Deduplicate(chunk);
        int isDuplicate = chunk->header.isDuplicate;

        //If chunk is unique compress & archive it.
        if(!isDuplicate) {
          sub_Compress(chunk);
        }

        write_chunk_to_file(fd_out, chunk);

        //prepare for next iteration
        chunk = temp;
        temp = NULL;
        split = 1;
      } else {
        //Due to technical limitations we can't distinguish the cases "no split" and "split at end of buffer"
        //This will result in some unnecessary (and unlikely) work but yields the correct result eventually.
        temp = chunk;
        chunk = NULL;
        split = 0;
      }
    } while(split);
  }

  free(rabintab);
  free(rabinwintab);

  close(fd_out);

  return NULL;
}

/*--------------------------------------------------------------------------*/
/* Encode
 * Compress an input stream
 *
 * Arguments:
 *   conf:    Configuration parameters
 *
 */
void Encode(config_t * _conf) {
  struct stat filestat;
  int32 fd;

  conf = _conf;

#ifdef ENABLE_STATISTICS
  init_stats(&stats);
#endif

  //Create chunk cache
  cache = hashtable_create(65536, hash_from_key_fn, keys_equal_fn, FALSE);
  if(cache == NULL) {
    printf("ERROR: Out of memory\n");
    exit(1);
  }

  struct thread_args generic_args;

  assert(!mbuffer_system_init());

  /* src file stat */
  if (stat(conf->infile, &filestat) < 0) 
      EXIT_TRACE("stat() %s failed: %s\n", conf->infile, strerror(errno));

  if (!S_ISREG(filestat.st_mode)) 
    EXIT_TRACE("not a normal file: %s\n", conf->infile);
#ifdef ENABLE_STATISTICS
  stats.total_input = filestat.st_size;
#endif //ENABLE_STATISTICS

  /* src file open */
  if((fd = open(conf->infile, O_RDONLY | O_LARGEFILE)) < 0) 
    EXIT_TRACE("%s file open error %s\n", conf->infile, strerror(errno));

  //Load entire file into memory if requested by user
  void *preloading_buffer = NULL;
  if(conf->preloading) {
    size_t bytes_read=0;
    int r;

    preloading_buffer = malloc(filestat.st_size);
    if(preloading_buffer == NULL)
      EXIT_TRACE("Error allocating memory for input buffer.\n");

    //Read data until buffer full
    while(bytes_read < (size_t)filestat.st_size) {
      r = read(fd, (char*)preloading_buffer+bytes_read, filestat.st_size-bytes_read);
      if(r<0) switch(errno) {
        case EAGAIN:
          EXIT_TRACE("I/O error: No data available\n");break;
        case EBADF:
          EXIT_TRACE("I/O error: Invalid file descriptor\n");break;
        case EFAULT:
          EXIT_TRACE("I/O error: Buffer out of range\n");break;
        case EINTR:
          EXIT_TRACE("I/O error: Interruption\n");break;
        case EINVAL:
          EXIT_TRACE("I/O error: Unable to read from file descriptor\n");break;
        case EIO:
          EXIT_TRACE("I/O error: Generic I/O error\n");break;
        case EISDIR:
          EXIT_TRACE("I/O error: Cannot read from a directory\n");break;
        default:
          EXIT_TRACE("I/O error: Unrecognized error\n");break;
      }
      if(r==0) break;
      bytes_read += r;
    }
    generic_args.input_file.size = filestat.st_size;
    generic_args.input_file.buffer = preloading_buffer;
  }

#ifdef TBB_FUNCOBJ

  generic_args.tid = 0;
  generic_args.nqueues = -1;
  generic_args.fd = fd;

#ifdef ENABLE_PARSEC_HOOKS
  __parsec_roi_begin();
#endif

  //Do the processing
  // TBBIntegratedPipeline(&generic_args);
  tbb::task_scheduler_init init( conf->nthreads );
  OuterPipeline outer( &generic_args );

#ifdef ENABLE_PARSEC_HOOKS
  __parsec_roi_end();
#endif

#else //serial version

  generic_args.tid = 0;
  generic_args.nqueues = -1;
  generic_args.fd = fd;

#ifdef ENABLE_PARSEC_HOOKS
  __parsec_roi_begin();
#endif

  //Do the processing
  SerialIntegratedPipeline(&generic_args);

#ifdef ENABLE_PARSEC_HOOKS
  __parsec_roi_end();
#endif

#endif //TBB_FUNCOBJ

  //clean up after preloading
  if(conf->preloading) {
    free(preloading_buffer);
  }

  /* clean up with the src file */
  if (conf->infile != NULL)
    close(fd);

  assert(!mbuffer_system_destroy());

  hashtable_destroy(cache, TRUE);

#ifdef ENABLE_STATISTICS
  /* dest file stat */
  if (stat(conf->outfile, &filestat) < 0) 
      EXIT_TRACE("stat() %s failed: %s\n", conf->outfile, strerror(errno));
  stats.total_output = filestat.st_size;

  //Analyze and print statistics
  if(conf->verbose) print_stats(&stats);
#endif //ENABLE_STATISTICS
}

