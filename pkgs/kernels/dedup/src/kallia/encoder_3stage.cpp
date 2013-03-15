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
#include <unistd.h>
#include <string.h>
#include "/home/kchronaki/scheduler/scheduler/wf_interface.h"
//#define fprintf(stderr, ...) 
extern "C"
{
#include "util.h"
#include "dedupdef.h"
#include "debug.h"

#include "dconfig.h"
#include "binheap.h"
#include "tree.h"
#include "queue.h"
#include "mbuffer.h"
#include "rabin.h"
#include "hashtable.h"
}
#include "encoder.h"

#define MAX_NUM_ITERATIONS 500
#define MAX_SPLITTED_CHUNKS 100000

#ifdef ENABLE_PTHREADS
#include "queue.h"
#include "binheap.h"
#include "tree.h"
#endif //ENABLE_PTHREADS

#ifdef ENABLE_GZIP_COMPRESSION
#include <zlib.h>
#endif //ENABLE_GZIP_COMPRESSION

#ifdef ENABLE_BZIP2_COMPRESSION
#include <bzlib.h>
#endif //ENABLE_BZIP2_COMPRESSION

#ifdef ENABLE_PTHREADS
#include <pthread.h>
#endif //ENABLE_PTHREADS

#ifdef ENABLE_PARSEC_HOOKS
#include <hooks.h>
#endif //ENABLE_PARSEC_HOOKS


#define INITIAL_SEARCH_TREE_SIZE 4096
#include <iostream>
#include <list>
using namespace std;
using namespace obj;

const int compress_type = COMPRESS_GZIP;
const int preloading = 1;

//The configuration block defined in main
config_t * confi;
int maxchcount[8];
//Hash table data structure & utility functions
struct hashtable *cache;
//-----------------TIMERS---------------------

#include <sys/time.h>

typedef struct {
        struct  timeval start, end;
        float   diff;
} stimer_t;

void stimer_tick(stimer_t *timer);
float stimer_tuck(stimer_t *timer, const char *msg);
void stimer_tick(stimer_t *timer)
{
        gettimeofday(&(timer->start), 0);
}

float stimer_tuck(stimer_t *timer, const char *msg)
{
        gettimeofday(&(timer->end), 0);

        timer->diff = (timer->end.tv_sec - timer->start.tv_sec)
                                + (timer->end.tv_usec - timer->start.tv_usec) * 0.000001;

        if (msg)
                printf("%s: %.3f seconds\n", msg, timer->diff);

        return timer->diff;
}


//--------------END-TIMERS---------------------
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
queue_t *deduplicate_que, *refine_que, *reorder_que, *compress_que;
#ifdef ENABLE_PTHREADS
//The queues between the pipeline stages
queue_t *deduplicate_que, *refine_que, *reorder_que, *compress_que;

//Merge two statistics records: s1=s1+s2
static void merge_stats(stats_t *s1, stats_t *s2) {
  int i;

  assert(s1!=NULL);
  assert(s2!=NULL);
  s1->total_input += s2->total_input;
  s1->total_dedup += s2->total_dedup;
  s1->total_compressed += s2->total_compressed;
  s1->total_output += s2->total_output;

  for(i=0; i<CHUNK_MAX_NUM; i++) {
    s1->nChunks[i] += s2->nChunks[i];
  }
  s1->nDuplicates += s2->nDuplicates;
}
#endif //ENABLE_PTHREADS
stats_t *thread_stats;
//Print statistics
static void print_stats(stats_t *s) {
  const unsigned int unit_str_size = 7; //elements in unit_str array
  const char *unit_str[] = {"Bytes", "KB", "MB", "GB", "TB", "PB", "EB"};
  unsigned int unit_idx = 0;
  size_t unit_div = 1;
	s = thread_stats;
  assert(s!=NULL);
	printf("IN THREAD STATS thread_stats->total input = %14.2f\n", (float)thread_stats->total_input/(float)(unit_div));
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
	printf("Mean size = %f\n", mean_size);
  printf("Mean data chunk size:          %14.2f %s (stddev: %.2f %s)\n", mean_size / 1024.0, "KB", sqrtf(var_size) / 1024.0, "KB");
  printf("Amount of duplicate chunks:    %14.2f%%\n", 100.0*(float)(s->nDuplicates)/(float)(nTotalChunks));
  printf("Data size after deduplication: %14.2f %s (compression factor: %.2fx)\n", (float)(s->total_dedup)/(float)(unit_div), unit_str[unit_idx], (float)(s->total_input)/(float)(s->total_dedup));
  printf("Data size after compression:   %14.2f %s (compression factor: %.2fx)\n", (float)(s->total_compressed)/(float)(unit_div), unit_str[unit_idx], (float)(s->total_dedup)/(float)(s->total_compressed));
  printf("Output overhead:               %14.2f%%\n", 100.0*(float)(s->total_output-s->total_compressed)/(float)(s->total_output));
}

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

#include <sys/stat.h>
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
    //EXIT_TRACE("Cannot open output file.");
	fd = open("output.dat.ddp", O_CREAT|O_TRUNC|O_WRONLY|O_TRUNC, S_IRGRP | S_IWUSR | S_IRUSR | S_IROTH);
	if (fd < 0) {
    	EXIT_TRACE("Cannot open output file.");
  }
  }
	
  //Write header
  if (write_header(fd, compress_type)){//confi->compress_type)) {
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
#ifdef ENABLE_PTHREADS
//NOTE: The parallel version checks the state of each chunk to make sure the
//      relevant data is available. If it is not then the function waits.
static void write_chunk_to_file(int fd, chunk_t *chunk) {
  assert(chunk!=NULL);

  //Find original chunk
  if(chunk->header.isDuplicate) chunk = chunk->compressed_data_ref;

  pthread_mutex_lock(&chunk->header.lock);
  while(chunk->header.state == CHUNK_STATE_UNCOMPRESSED) {
    pthread_cond_wait(&chunk->header.update, &chunk->header.lock);
  }

  //state is now guaranteed to be either COMPRESSED or FLUSHED
  if(chunk->header.state == CHUNK_STATE_COMPRESSED) {
    //Chunk data has not been written yet, do so now
    write_file(fd, TYPE_COMPRESS, chunk->compressed_data.n, chunk->compressed_data.ptr);
    mbuffer_free(&chunk->compressed_data);
    chunk->header.state = CHUNK_STATE_FLUSHED;
  } else {
    //Chunk data has been written to file before, just write SHA1
    write_file(fd, TYPE_FINGERPRINT, SHA1_LEN, (unsigned char *)(chunk->sha1));
  }
  pthread_mutex_unlock(&chunk->header.lock);
}
#else
int numwrites = 0;
//NOTE: The serial version relies on the fact that chunks are processed in-order,
//      which means if it reaches the function it is guaranteed all data is ready.
static void write_chunk_to_file(int fd, chunk_t *chunk) {
  assert(chunk!=NULL);
	
  if(!chunk->header.isDuplicate) {
    //Unique chunk, data has not been written yet, do so now
    write_file(fd, TYPE_COMPRESS, chunk->compressed_data.n, chunk->compressed_data.ptr);
    mbuffer_free(&chunk->compressed_data);
  } else {
		
    //Duplicate chunk, data has been written to file before, just write SHA1
    write_file(fd, TYPE_FINGERPRINT, SHA1_LEN, (unsigned char *)(chunk->sha1));
  }
}
#endif //ENABLE_PTHREADS

int rf_win;
int rf_win_dataprocess;

/*
 * Computational kernel of compression stage
 *
 * Actions performed:
 *  - Compress a data chunk
 */
void sub_Compress(chunk_t **chunk_ptr) {
    size_t n;
    int r;
	chunk_t *chunk = *chunk_ptr;
    assert(chunk!=NULL);
    //compress the item and add it to the database
#ifdef ENABLE_PTHREADS
    pthread_mutex_lock(&chunk->header.lock);
    assert(chunk->header.state == CHUNK_STATE_UNCOMPRESSED);
#endif //ENABLE_PTHREADS
    switch (compress_type) {
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
//#ifdef ENABLE_GZIP_COMPRESSION
      case COMPRESS_GZIP:
        //Gzip compression buffer must be at least 0.1% larger than source buffer plus 12 bytes
        n = chunk->uncompressed_data.n + (chunk->uncompressed_data.n >> 9) + 12;
        r = mbuffer_create(&chunk->compressed_data, n);
        if(r != 0) {
          EXIT_TRACE("Creation of compression buffer failed.\n");
        }
        //compress the block
        r = compress(chunk->compressed_data.ptr, &n, chunk->uncompressed_data.ptr, chunk->uncompressed_data.n);
        if (r != Z_OK) {
          EXIT_TRACE("Compression failed\n");
        }
        //Shrink buffer to actual size
        if(n < chunk->compressed_data.n) {
          r = mbuffer_realloc(&chunk->compressed_data, n);
          assert(r == 0);
        }
        break;
//#endif //ENABLE_GZIP_COMPRESSION
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
	(*chunk_ptr)->header.state = CHUNK_STATE_COMPRESSED;

#ifdef ENABLE_PTHREADS
    
    pthread_cond_broadcast(&chunk->header.update);
    pthread_mutex_unlock(&chunk->header.lock);
#endif //ENABLE_PTHREADS

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
int sub_Deduplicate(chunk_t **chunk_ptr) {
  int isDuplicate;
  chunk_t *entry;
	chunk_t *chunk = *chunk_ptr;

  assert(chunk!=NULL);
  assert(chunk->uncompressed_data.ptr!=NULL);

  SHA1_Digest((*chunk_ptr)->uncompressed_data.ptr, (*chunk_ptr)->uncompressed_data.n, (unsigned char *)((*chunk_ptr)->sha1));

  //Query database to determine whether we've seen the data chunk before
#ifdef ENABLE_PTHREADS
  pthread_mutex_t *ht_lock = hashtable_getlock(cache, (void *)(chunk->sha1));
  pthread_mutex_lock(ht_lock);
#endif
	//chunk->sha1 = (void *)(*chunk_ptr)->sha1;
  entry = (chunk_t *)hashtable_search(cache, (void *)((*chunk_ptr)->sha1));//chunk->sha1));
  isDuplicate = (entry != NULL);
  (*chunk_ptr)->header.isDuplicate = isDuplicate;
  if (!isDuplicate) {
    // Cache miss: Create entry in hash table and forward data to compression stage
#ifdef ENABLE_PTHREADS
    pthread_mutex_init(&chunk->header.lock, NULL);
    pthread_cond_init(&chunk->header.update, NULL);
#endif
    //NOTE: chunk->compressed_data.buffer will be computed in compression stage
    if (hashtable_insert(cache, (void *)((*chunk_ptr)->sha1), (void *)(*chunk_ptr)) == 0) {
      EXIT_TRACE("hashtable_insert failed");
    }
  } else {
    // Cache hit: Skipping compression stage
    (*chunk_ptr)->compressed_data_ref = entry;
    mbuffer_free(&(*chunk_ptr)->uncompressed_data);
  }
#ifdef ENABLE_PTHREADS
  pthread_mutex_unlock(ht_lock);
#endif

  return isDuplicate;
}



/* 
 * Integrate all computationally intensive pipeline
 * stages to improve cache efficiency.
 */
void *SerialIntegratedPipeline(void * targs) {
  struct thread_args *args = (struct thread_args *)targs;
  size_t preloading_buffer_seek = 0;
  int fd = args->fd;
  int fd_out = create_output_file(confi->outfile);
  int r;

  chunk_t *temp = NULL;
  chunk_t *chunk = NULL;
  u32int * rabintab = malloc(256*sizeof rabintab[0]);
  u32int * rabinwintab = malloc(256*sizeof rabintab[0]);
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
    if(confi->preloading) {
      size_t max_read = MIN(MAXBUF, args->input_file.size-preloading_buffer_seek);
      memcpy(chunk->uncompressed_data.ptr+bytes_left, args->input_file.buffer+preloading_buffer_seek, max_read);
      bytes_read = max_read;
      preloading_buffer_seek += max_read;
    } else {
      while(bytes_read < MAXBUF) {
        r = read(fd, chunk->uncompressed_data.ptr+bytes_left+bytes_read, MAXBUF-bytes_read);
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
      thread_stats->nChunks[CHUNK_SIZE_TO_SLOT(chunk->uncompressed_data.n)]++;
#endif //ENABLE_STATISTICS

      //Deduplicate
      int isDuplicate = sub_Deduplicate(&chunk);
#ifdef ENABLE_STATISTICS
      if(isDuplicate) {
        thread_stats->nDuplicates++;
      } else {
        thread_stats->total_dedup += chunk->uncompressed_data.n;
      }
#endif //ENABLE_STATISTICS

      //If chunk is unique compress & archive it.
      if(!isDuplicate) {
        sub_Compress(&chunk);
#ifdef ENABLE_STATISTICS
        thread_stats->total_compressed += chunk->compressed_data.n;
#endif //ENABLE_STATISTICS
      }
		//if this is the last chunk write it and return
      write_chunk_to_file(fd_out, chunk);
      if(chunk->header.isDuplicate) {
        free(chunk);
        chunk=NULL;
      }

      //stop fetching from input buffer, terminate processing
      break;
    }

    //partition input block into fine-granular chunks
    int split;
    do {		//while(split)
      split = 0;
      //Try to split the buffer
      int offset = rabinseg(chunk->uncompressed_data.ptr, chunk->uncompressed_data.n, rf_win_dataprocess, rabintab, rabinwintab);
      //Did we find a split location?
      if(offset == 0) {
        //Split found at the very beginning of the buffer (should never happen due to technical limitations)
        assert(0);
        split = 0;
      } else if(offset < chunk->uncompressed_data.n) {
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
        thread_stats->nChunks[CHUNK_SIZE_TO_SLOT(chunk->uncompressed_data.n)]++;
#endif //ENABLE_STATISTICS

        //Deduplicate
        int isDuplicate = sub_Deduplicate(&chunk);
#ifdef ENABLE_STATISTICS
        if(isDuplicate) {
          thread_stats->nDuplicates++;
        } else {
          thread_stats->total_dedup += chunk->uncompressed_data.n;
        }
#endif //ENABLE_STATISTICS

        //If chunk is unique compress & archive it.
        if(!isDuplicate) {
          sub_Compress(&chunk);
#ifdef ENABLE_STATISTICS
          thread_stats->total_compressed += chunk->compressed_data.n;
#endif //ENABLE_STATISTICS
        }

        write_chunk_to_file(fd_out, chunk);
        if(chunk->header.isDuplicate){
          free(chunk);
          chunk=NULL;
        }

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
  }	//while(1)

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

  confi = _conf;

#ifdef ENABLE_STATISTICS
//  init_stats(&stats);
#endif

  //Create chunk cache
  cache = hashtable_create(65536, hash_from_key_fn, keys_equal_fn, FALSE);
  if(cache == NULL) {
    printf("ERROR: Out of memory\n");
    exit(1);
  }

#ifdef ENABLE_PTHREADS
  struct thread_args data_process_args;
  int i;

  //queue allocation & initialization
  const int nqueues = (confi->nthreads / MAX_THREADS_PER_QUEUE) +
                      ((confi->nthreads % MAX_THREADS_PER_QUEUE != 0) ? 1 : 0);
  deduplicate_que = malloc(sizeof(queue_t) * nqueues);
  refine_que = malloc(sizeof(queue_t) * nqueues);
  reorder_que = malloc(sizeof(queue_t) * nqueues);
  compress_que = malloc(sizeof(queue_t) * nqueues);
  if( (deduplicate_que == NULL) || (refine_que == NULL) || (reorder_que == NULL) || (compress_que == NULL)) {
    printf("Out of memory\n");
    exit(1);
  }
  int threads_per_queue;
  for(i=0; i<nqueues; i++) {
    if (i < nqueues -1 || confi->nthreads %MAX_THREADS_PER_QUEUE == 0) {
      //all but last queue
      threads_per_queue = MAX_THREADS_PER_QUEUE;
    } else {
      //remaining threads work on last queue
      threads_per_queue = confi->nthreads %MAX_THREADS_PER_QUEUE;
    }

    //call queue_init with threads_per_queue
    queue_init(&deduplicate_que[i], QUEUE_SIZE, threads_per_queue);
    queue_init(&refine_que[i], QUEUE_SIZE, 1);
    queue_init(&reorder_que[i], QUEUE_SIZE, threads_per_queue);
    queue_init(&compress_que[i], QUEUE_SIZE, threads_per_queue);
  }
#else
  struct thread_args generic_args;
#endif //ENABLE_PTHREADS

  assert(!mbuffer_system_init());

  /* src file stat */
  if (stat(confi->infile, &filestat) < 0) 
      EXIT_TRACE("stat() %s failed: %s\n", confi->infile, strerror(errno));

  if (!S_ISREG(filestat.st_mode)) 
    EXIT_TRACE("not a normal file: %s\n", confi->infile);
#ifdef ENABLE_STATISTICS
  thread_stats->total_input = filestat.st_size;
#endif //ENABLE_STATISTICS

  /* src file open */
  if((fd = open(confi->infile, O_RDONLY | O_LARGEFILE)) < 0) 
    EXIT_TRACE("%s file open error %s\n", confi->infile, strerror(errno));

  //Load entire file into memory if requested by user
  void *preloading_buffer = NULL;
  if(confi->preloading) {
    size_t bytes_read=0;
    int r;

    preloading_buffer = malloc(filestat.st_size);
    if(preloading_buffer == NULL)
      EXIT_TRACE("Error allocating memory for input buffer.\n");

    //Read data until buffer full
    while(bytes_read < filestat.st_size) {
      r = read(fd, preloading_buffer+bytes_read, filestat.st_size-bytes_read);
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
#ifdef ENABLE_PTHREADS
    data_process_args.input_file.size = filestat.st_size;
    data_process_args.input_file.buffer = preloading_buffer;
#else
    generic_args.input_file.size = filestat.st_size;
    generic_args.input_file.buffer = preloading_buffer;
#endif //ENABLE_PTHREADS
  }

#ifdef ENABLE_PTHREADS
  /* Variables for 3 thread pools and 2 pipeline stage threads.
   * The first and the last stage are serial (mostly I/O).
   */
  pthread_t threads_anchor[MAX_THREADS],
    threads_chunk[MAX_THREADS],
    threads_compress[MAX_THREADS],
    threads_send, threads_process;

  data_process_args.tid = 0;
  data_process_args.nqueues = nqueues;
  data_process_args.fd = fd;

#ifdef ENABLE_PARSEC_HOOKS
    __parsec_roi_begin();
#endif

  //thread for first pipeline stage (input)
  pthread_create(&threads_process, NULL, Fragment, &data_process_args);

  //Create 3 thread pools for the intermediate pipeline stages
  struct thread_args anchor_thread_args[confi->nthreads];
  for (i = 0; i < confi->nthreads; i ++) {
     anchor_thread_args[i].tid = i;
     pthread_create(&threads_anchor[i], NULL, FragmentRefine, &anchor_thread_args[i]);
  }

  struct thread_args chunk_thread_args[confi->nthreads];
  for (i = 0; i < confi->nthreads; i ++) {
    chunk_thread_args[i].tid = i;
    pthread_create(&threads_chunk[i], NULL, Deduplicate, &chunk_thread_args[i]);
  }

  struct thread_args compress_thread_args[confi->nthreads];
  for (i = 0; i < confi->nthreads; i ++) {
    compress_thread_args[i].tid = i;
    pthread_create(&threads_compress[i], NULL, Compress, &compress_thread_args[i]);
  }

  //thread for last pipeline stage (output)
  struct thread_args send_block_args;
  send_block_args.tid = 0;
  send_block_args.nqueues = nqueues;
  pthread_create(&threads_send, NULL, Reorder, &send_block_args);

  /*** parallel phase ***/

  //Return values of threads
  stats_t *threads_anchor_rv[confi->nthreads];
  stats_t *threads_chunk_rv[confi->nthreads];
  stats_t *threads_compress_rv[confi->nthreads];

  //join all threads 
  pthread_join(threads_process, NULL);
  for (i = 0; i < confi->nthreads; i ++)
    pthread_join(threads_anchor[i], (void **)&threads_anchor_rv[i]);
  for (i = 0; i < confi->nthreads; i ++)
    pthread_join(threads_chunk[i], (void **)&threads_chunk_rv[i]);
  for (i = 0; i < confi->nthreads; i ++)
    pthread_join(threads_compress[i], (void **)&threads_compress_rv[i]);
  pthread_join(threads_send, NULL);

#ifdef ENABLE_PARSEC_HOOKS
  __parsec_roi_end();
#endif

  /* free queues */
  for(i=0; i<nqueues; i++) {
    queue_destroy(&deduplicate_que[i]);
    queue_destroy(&refine_que[i]);
    queue_destroy(&reorder_que[i]);
    queue_destroy(&compress_que[i]);
  }
  free(deduplicate_que);
  free(refine_que);
  free(reorder_que);
  free(compress_que);

#ifdef ENABLE_STATISTICS
  //Merge everything into global `stats' structure
  for(i=0; i<confi->nthreads; i++) {
    merge_stats(&stats, threads_anchor_rv[i]);
    free(threads_anchor_rv[i]);
  }
  for(i=0; i<confi->nthreads; i++) {
    merge_stats(&stats, threads_chunk_rv[i]);
    free(threads_chunk_rv[i]);
  }
  for(i=0; i<confi->nthreads; i++) {
    merge_stats(&stats, threads_compress_rv[i]);
    free(threads_compress_rv[i]);
  }
#endif //ENABLE_STATISTICS

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

#endif //ENABLE_PTHREADS

  //clean up after preloading
  if(confi->preloading) {
    free(preloading_buffer);
  }

  /* clean up with the src file */
  if (confi->infile != NULL)
    close(fd);

  assert(!mbuffer_system_destroy());

  hashtable_destroy(cache, TRUE);

#ifdef ENABLE_STATISTICS
  /* dest file stat */
  if (stat(confi->outfile, &filestat) < 0)
      EXIT_TRACE("stat() %s failed: %s\n", confi->outfile, strerror(errno));
  thread_stats->total_output = filestat.st_size;

  //Analyze and print statistics
  if(confi->verbose) print_stats(thread_stats);
#endif //ENABLE_STATISTICS
}


inline void errorNumber(int r)
{
	if(r<0) switch(errno) 
	{
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
	//if(r==0) break;

}

//-------------------------------------------SWAN IMPLEMENTATION-------------------------------------------


/*
 * Pipeline stage function of reorder stage
 *
 * Actions performed:
 *  - Receive chunks from compression and deduplication stage
 *  - Check sequence number of each chunk to determine correct order
 *  - Cache chunks that arrive out-of-order until predecessors are available
 *  - Write chunks in-order to file (or preloading buffer)
 *
 * Notes:
 *  - This function blocks if the compression stage has not finished supplying
 *    the compressed data for a duplicate chunk.
 */
void ReorderSwan(struct thread_args *args, indep<chunk_t*>chunk_obj, int fd, inoutdep<int> wait)
{
	chunk_t *chunk;
	chunk = (chunk_t*)chunk_obj;
	// if(chunk == NULL)
	// {
		// fprintf(stderr, "NULL item in the LIST!!!!!!\n");
		// exit(0);					
	// }
	//fprintf(stderr, "In Reorder: chunk->l1num = %d, chunk->l2num = %d\n", chunk->sequence.l1num, chunk->sequence.l2num);
	leaf_call(write_chunk_to_file, fd, chunk);
	return;
}


//stage 4 function
void CompressSwan(struct thread_args *args, inoutdep<chunk_t*> chunk_obj)
{
	chunk_t * chunk = (chunk_t*)chunk_obj;
	assert(chunk!=NULL);
	if(chunk->header.isDuplicate == 0 )
		leaf_call(sub_Compress, &chunk);

	chunk_obj = chunk;
#ifdef ENABLE_STATISTICS
	thread_stats->total_compressed += (chunk)->compressed_data.n;
#endif //ENABLE_STATISTICS
}


//stage 3 function
void DeduplicateSwan(struct thread_args *args, chunk_t ** chunk_ptr)
{
	chunk_t *chunk = *chunk_ptr;
	assert(chunk!=NULL);

		//Do the processing
		int isDuplicate = sub_Deduplicate(chunk_ptr);

	#ifdef ENABLE_STATISTICS
		if(isDuplicate)  	thread_stats->nDuplicates++;
		else   				thread_stats->total_dedup += chunk->uncompressed_data.n;
	#endif //ENABLE_STATISTICS		
}

//stage 2 function
void FragmentRefineSwan(struct thread_args * targs, chunk_t *chunk, u32int * rabintab, u32int * rabinwintab, 
							int it, list<chunk_t*>*chunklist_ptr)
{
	struct thread_args *args = (struct thread_args *)targs;
	int r;
	list<chunk_t*>chunklist = *chunklist_ptr;
	chunk_t *temp, *newchunk;
	assert(chunk!=NULL);
	rabininit(rf_win, rabintab, rabinwintab);
	int split;
	sequence_number_t chcount = 0;
	do
	{
		newchunk = (chunk_t*)malloc(sizeof(chunk_t));
		assert(chunk!=NULL);
	
		//Find next anchor with Rabin fingerprint
		int offset = rabinseg((&chunk->uncompressed_data)->ptr, (&chunk->uncompressed_data)->n, rf_win, rabintab, rabinwintab);
		//Can we split the buffer?
		if(offset < chunk->uncompressed_data.n)
		{
			//Allocate a new chunk and create a new memory buffer
			temp = (chunk_t *)malloc(sizeof(chunk_t));
			if(temp==NULL) 			EXIT_TRACE("Memory allocation failed.\n");
			temp->header.state = chunk->header.state;
			temp->sequence.l1num = chunk->sequence.l1num;
			//split it into two pieces
			r = mbuffer_split(&chunk->uncompressed_data, &temp->uncompressed_data, offset);
			if(r!=0) 					EXIT_TRACE("Unable to split memory buffer.\n");

			//Set correct state and sequence numbers
			chunk->sequence.l2num = chcount;
			chunk->isLastL2Chunk = FALSE;
			newchunk->sequence.l2num = chcount;
			newchunk->isLastL2Chunk = FALSE;
			newchunk->header.state = chunk->header.state;
			newchunk->sequence.l1num = chunk->sequence.l1num;// = it;
			r = mbuffer_create(&newchunk->uncompressed_data, chunk->uncompressed_data.n);
			if(r!=0) 					EXIT_TRACE("Unable to allocate memory buffer...\n");
			memcpy(newchunk->uncompressed_data.ptr, chunk->uncompressed_data.ptr, chunk->uncompressed_data.n);
			newchunk->uncompressed_data.n = chunk->uncompressed_data.n;		//n is now offset
  			leaf_call(DeduplicateSwan, args, &newchunk);
			chcount++;
	#ifdef ENABLE_STATISTICS
			  //update statistics
			thread_stats->nChunks[CHUNK_SIZE_TO_SLOT(chunk->uncompressed_data.n)]++;
	#endif //ENABLE_STATISTICS
			
			//prepare for next iteration
			chunk = temp;
			split = 1;
		}
		else
		{
			//End of buffer reached, don't split but simply enqueue it
			//Set correct state and sequence numbers
			chunk->sequence.l2num = chcount;
			chunk->isLastL2Chunk = TRUE;
	
	#ifdef ENABLE_STATISTICS
			  //update statistics
			thread_stats->nChunks[CHUNK_SIZE_TO_SLOT(chunk->uncompressed_data.n)]++;
	#endif //ENABLE_STATISTICS
		
			newchunk->header.state = chunk->header.state;
			newchunk->sequence.l1num = chunk->sequence.l1num;// = it;
			r = mbuffer_create(&newchunk->uncompressed_data, chunk->uncompressed_data.n);
			if(r!=0) 					EXIT_TRACE("Unable to allocate memory buffer...\n");
			memcpy(newchunk->uncompressed_data.ptr, chunk->uncompressed_data.ptr, chunk->uncompressed_data.n);
			newchunk->uncompressed_data.n = chunk->uncompressed_data.n;		//n is not offset now
			newchunk->sequence.l2num = chcount;
			newchunk->isLastL2Chunk = TRUE;			
			leaf_call(DeduplicateSwan,args, &newchunk);
			chcount++;
			//prepare for next iteration
			chunk = NULL;
			split = 0;
		}
		chunklist.push_back(newchunk);
	} while(split);
	*chunklist_ptr = chunklist;
}


inline chunk_t* FragmentSplit(int split, int *split_ptr, sequence_number_t *anchorcount, chunk_t **chunk_ptr, chunk_t **temp_ptr, u32int * rabintab, u32int * rabinwintab, int *condition, int mode)
{
	chunk_t *chunk = *chunk_ptr;
	chunk_t *temp = *temp_ptr;
	int r;
	if(!(split == 0 && mode == 1))
	{
		do
		{
			split = 0;
			*split_ptr = split;
			//Try to split the buffer at least ANCHOR_JUMP bytes away from its beginning
			if(ANCHOR_JUMP < chunk->uncompressed_data.n)
			{
				int offset = rabinseg(chunk->uncompressed_data.ptr+ANCHOR_JUMP, chunk->uncompressed_data.n-ANCHOR_JUMP, rf_win_dataprocess, rabintab, rabinwintab);
				//Did we find a split location?
				if(offset == 0)
				{
					//Split found at the very beginning of the buffer (should never happen due to technical limitations)
					assert(0);
					exit(1);
					split = 0;
					*split_ptr = split;
				}
				else if(offset + ANCHOR_JUMP < chunk->uncompressed_data.n)
				{
					//Split found somewhere in the middle of the buffer
					//Allocate a new chunk and create a new memory buffer
					temp = (chunk_t *)malloc(sizeof(chunk_t));
					if(temp==NULL) 		EXIT_TRACE("Memory allocation failed.\n");
					//split it into two pieces
					r = mbuffer_split(&chunk->uncompressed_data, &temp->uncompressed_data, offset + ANCHOR_JUMP);
					if(r!=0) 				EXIT_TRACE("Unable to split memory buffer.\n");

					(temp)->header.state = CHUNK_STATE_UNCOMPRESSED;
					(temp)->sequence.l1num = *anchorcount;
					(*anchorcount)++;
					chunk_t *result = chunk;
					//prepare for next iteration
					*chunk_ptr = temp;
					*temp_ptr = NULL;
					split = 1;
					*split_ptr = split;
					return result;
				}
				else 
				{
					//Due to technical limitations we can't distinguish the cases "no split" and "split at end of buffer"
					//This will result in some unnecessary (and unlikely) work but yields the correct result eventually.
					*temp_ptr = chunk;
					*chunk_ptr = NULL;
					split = 0;
					*split_ptr = split;
					continue;
				}
			} //if(ANCHOR_JUMP < chunk->uncompressed_data.n) 
			else 
			{
				//NOTE: We don't process the stub, instead we try to read in more data so we might be able to find a proper split.
				//Only once the end of the file is reached do we get a genuine stub which will be enqueued right after the read operation.
				*temp_ptr = chunk;
				if(chunk == NULL)		fprintf(stderr, "chunk is NULL what the fuck???\n");
				*chunk_ptr = NULL;
				split = 0;
				*split_ptr = split;
				return NULL;
				if(mode == 0)		*condition = 1;
			}
		}while(split);
	}

	return chunk;
}

//Stage 1 function
chunk_t * FragmentSwan(struct thread_args  *targs, chunk_t **chunk_ptr, chunk_t **temp_ptr, 
						int *read_done, size_t *bytes_left_ptr, size_t *preloading_buffer_seek_ptr, 
						int *split_ptr, u32int * rabintab, u32int * rabinwintab)
						//amount of data left over in last_mbuffer from previous iteration
{
	struct thread_args *args = (struct thread_args *)targs;
	size_t preloading_buffer_seek = *preloading_buffer_seek_ptr;
	int qid = 0;
	int fd = args->fd;
	int i;
	size_t bytes_left = *bytes_left_ptr;
	int condition = 0;
	sequence_number_t anchorcount = 0;
	int r, split = *split_ptr;
	chunk_t *result, *chunk = *chunk_ptr, *temp = *temp_ptr;

	if(split)
	{
		result = FragmentSplit(split, split_ptr, &anchorcount, chunk_ptr, temp_ptr, rabintab, rabinwintab, &condition, 1);
		if(result != NULL)		return result;
	}
	temp = *temp_ptr;

	//read from input file / buffer
	while (1)
	{
		//Check how much data left over from previous iteration resp. create an initial chunk
		if(temp != NULL)
		{
			bytes_left = temp->uncompressed_data.n;
			*bytes_left_ptr = bytes_left;
		}
		else 
		{
			bytes_left = 0;
			*bytes_left_ptr = bytes_left;
		}
		 
		//Allocate a new chunk and create a new memory buffer
		*chunk_ptr = (chunk_t *)malloc(sizeof(chunk_t));
		if(*chunk_ptr==NULL) 			EXIT_TRACE("Memory allocation failed.\n");
		r = mbuffer_create(&(*chunk_ptr)->uncompressed_data, MAXBUF+bytes_left);
		if(r!=0)   				EXIT_TRACE("Unable to initialize memory buffer.\n");
		chunk = *chunk_ptr;
		if(bytes_left > 0)
		{
		   //FIXME: Short-circuit this if no more data available

		   //"Extension" of existing buffer, copy sequence number and left over data to beginning of new buffer
		   chunk->header.state = CHUNK_STATE_UNCOMPRESSED;
		   chunk->sequence.l1num = temp->sequence.l1num;

		   //NOTE: We cannot safely extend the current memory region because it has already been given to another thread
		   memcpy(chunk->uncompressed_data.ptr, temp->uncompressed_data.ptr, temp->uncompressed_data.n);
		   mbuffer_free(&temp->uncompressed_data);
		   free(temp);
		   *temp_ptr = NULL;
			temp = *temp_ptr;
		}
		else 
		{
		   //brand new mbuffer, increment sequence number
		   chunk->header.state = CHUNK_STATE_UNCOMPRESSED;
		   chunk->sequence.l1num = anchorcount;
		   anchorcount++;
		}
		 //Read data until buffer full
		size_t bytes_read=0;
		size_t max_read;
		if(preloading)
		{
			max_read = MIN(MAXBUF, args->input_file.size-preloading_buffer_seek);
		   memcpy(chunk->uncompressed_data.ptr+bytes_left, args->input_file.buffer+preloading_buffer_seek, max_read);
		   bytes_read = max_read;
		   preloading_buffer_seek += max_read;
			*preloading_buffer_seek_ptr = preloading_buffer_seek;
			result = chunk;
		}
		else 
		{
			while(bytes_read < MAXBUF)
			{
				r = read(fd, chunk->uncompressed_data.ptr+bytes_left+bytes_read, MAXBUF-bytes_read);
				errorNumber(r);
				if(r == 0) 		break;
				bytes_read += r;
		   }
		}
		//No data left over from last iteration and also nothing new read in, simply clean up and quit
		if(bytes_left + bytes_read == 0)
		{
			mbuffer_free(&chunk->uncompressed_data);
			free(chunk);
			*chunk_ptr = NULL;
			*read_done = 1;//se auto to shmeio mporw na valw mia metavlhth read_done kai na thn 8etw =1
			break;
		}
		 //Shrink buffer to actual size
		if(bytes_left+bytes_read < chunk->uncompressed_data.n)
		{
			r = mbuffer_realloc(&chunk->uncompressed_data, bytes_left+bytes_read);
			assert(r == 0);
		}
		 //Check whether any new data was read in, enqueue last chunk if not
		if(bytes_read == 0)
		{
			//put it into send buffer return this chunk
			//NOTE: No need to empty a full send_buf, we will break now and pass everything on to the queue
			*chunk_ptr = NULL;
			result = chunk;
			break;
		}

		result  = FragmentSplit(/*split:*/1, split_ptr, &anchorcount, chunk_ptr, temp_ptr, rabintab, rabinwintab, &condition, 0);
		//partition input block into large, coarse-granular chunks
		assert(chunk!=NULL);
		if(condition){	fprintf(stderr, "continuing....\n");		continue;}

		break;
	}		//while(1)
	

	return result;
	free(rabintab);
	free(rabinwintab);
	return;
}

inline int stage1(struct thread_args *args, chunk_t **chunk_ptr, chunk_t **temp_ptr, size_t *pr_buff_seek_ptr, int *split_ptr,
					u32int ** rabintab_ptr, u32int ** rabinwintab_ptr, u32int ** rabintabRef_ptr, u32int ** rabinwintabRef_ptr,
					int iteration, outdep<list<chunk_t*>> chunklist_obj)
{
	int read_done = 0;
	size_t bytes_left;
	size_t pr_buff_seek = *pr_buff_seek_ptr;
	int split = *split_ptr;
	u32int 	* rabintab = *rabintab_ptr,
			* rabinwintab = *rabinwintab_ptr, 
			* rabintabRef = *rabintabRef_ptr,  
			* rabinwintabRef = *rabinwintabRef_ptr;
	chunk_t *result = NULL;
	list<chunk_t*>chunklist = (list<chunk_t*>)chunklist_obj;

	result = leaf_call(FragmentSwan, args, chunk_ptr, temp_ptr, &read_done, &bytes_left, &pr_buff_seek, &split, rabintab, rabinwintab);	
	*split_ptr = split;
	*pr_buff_seek_ptr = pr_buff_seek;
	*split_ptr = split;
	
	if(!read_done || iteration == 0)
	{
		call(FragmentRefineSwan, args, result, rabintabRef, rabinwintabRef, iteration, &chunklist);
	}
	chunklist_obj = chunklist;
	return read_done;
}
void nestedPipeline(struct thread_args *args, indep<list<chunk_t*>>chunks, int fd, int iteration, inoutdep<int>align)
{
	list<chunk_t*>			chunklist;
	chunklist = (list<chunk_t*>)chunks;
	list<chunk_t*>::iterator iter;
	chunk_t * list_item;
	
	unversioned<int> wait;
	unversioned<int> wait2;
	for(iter = chunklist.begin(); iter != chunklist.end(); iter++)
	{
		object_t<chunk_t*> chunk_obj;
		list_item = *iter;
		chunk_obj = (chunk_t*)list_item;
		spawn(CompressSwan, args, (inoutdep<chunk_t*>) chunk_obj);
		spawn(ReorderSwan, args, (indep<chunk_t*>)chunk_obj, fd, (inoutdep<int>) wait2);
	}
	ssync();

}
void SwanPipeline(struct thread_args *args)
{
	int read_done = 0, split = 0, iteration = -1;
	size_t bytes_left;// = (size_t *)malloc(sizeof(size_t));
	chunk_t ** temp_ptr = (chunk_t**)malloc(sizeof(chunk_t*));
	chunk_t *temp = *temp_ptr;
	size_t pr_buff_seek = 0;
	temp = NULL;	

	u32int * rabintab = malloc(256*sizeof rabintab[0]);
	u32int * rabinwintab = malloc(256*sizeof rabintab[0]);
	u32int * rabintabRef = malloc(256*sizeof rabintabRef[0]);
	u32int * rabinwintabRef = malloc(256*sizeof rabintabRef[0]);
	if(rabintab == NULL || rabinwintab == NULL || rabintabRef == NULL || rabinwintabRef == NULL)
		EXIT_TRACE("Memory allocation failed.\n");
	u32int ** rabintab_ptr = &rabintab;
	u32int ** rabinwintab_ptr = &rabintab;
	u32int ** rabintabRef_ptr = &rabintabRef;
	u32int ** rabinwintabRef_ptr = &rabintabRef;

	chunk_t *chunk = NULL, *result = NULL, *curr_chunk = NULL;
	split = 0;
	rf_win_dataprocess = 0;
	rabininit(rf_win_dataprocess, rabintab, rabinwintab);
	int fd = create_output_file(confi->outfile);
	object_t<int> wait; 
	wait = 0;
	while(1)
	{
		iteration++;
		list<chunk_t*>	chunklist;
		list<chunk_t*>::iterator iter;
		object_t<list<chunk_t*>> chunklist_obj;
		chunklist_obj = (list<chunk_t*>)chunklist;
		read_done =  call(stage1, 	args, 
									&chunk, 
									&temp, 
									&pr_buff_seek, 
									&split,
									rabintab_ptr, 
									rabinwintab_ptr, 
									rabintabRef_ptr, 
									rabinwintabRef_ptr,
									iteration, 
									(outdep<list<chunk_t*>>)chunklist_obj);
					
		//chunklist = (list<chunk_t*>)chunklist_obj;
		if(!read_done || iteration == 0)
		{
			spawn(nestedPipeline, args, (indep<list<chunk_t*>>)chunklist_obj, fd, iteration, (inoutdep<int>)wait);
		}
		else
		{
			break;
		}
	}
	ssync();
	close(fd);
}



/*--------------------------------------------------------------------------*/
/* Encode
 * Compress an input stream
 *
 * Arguments:
 *   conf:    Configuration parameters
 *
 */
void EncodeSwan(config_t * _conf)
{
	fprintf(stderr, "Encode Starting.... conf-> inputfile = %s\n", _conf->infile);
	struct stat filestat;
	int32 fd;
	confi = _conf;
	//compress_type = _conf->compress_type;
	struct thread_args data_process_args;
	//preloading = confi->preloading;
	#ifdef ENABLE_STATISTICS
	thread_stats = (stats_t *)malloc(sizeof(stats_t));
	init_stats(thread_stats);
	#endif

	//Create chunk cache
	cache = hashtable_create(65536, hash_from_key_fn, keys_equal_fn, FALSE);
	if(cache == NULL)
	{
		printf("ERROR: Out of memory\n");
		exit(1);
	}

	assert(!mbuffer_system_init());

	/* src file stat */
	if (stat(_conf->infile, &filestat) < 0)	EXIT_TRACE("stat() %s failed: %s\n", _conf->infile, strerror(errno));
	//if (stat("kallia.txt", &filestat) < 0)	EXIT_TRACE("stat() %s failed: %s\n", _conf->infile, strerror(errno));

	if (!S_ISREG(filestat.st_mode))		EXIT_TRACE("not a normal file: %s\n", _conf->infile);

	#ifdef ENABLE_STATISTICS
	thread_stats->total_input = filestat.st_size;
	fprintf(stderr, "file size = %lf \n", thread_stats->total_input);
	fprintf(stderr, "confi->infile = %s\n", _conf->infile);
	#endif //ENABLE_STATISTICS

	/* src file open */
	if((fd = open(_conf->infile, O_RDONLY | O_LARGEFILE)) < 0)	EXIT_TRACE("%s file open error %s\n", _conf->infile, strerror(errno));

	//Load entire file into memory if requested by user
	void *preloading_buffer = NULL;
	if(_conf->preloading)
	{
		 size_t bytes_read=0;
		 int r;

		 preloading_buffer = malloc(filestat.st_size);
		 if(preloading_buffer == NULL)	EXIT_TRACE("Error allocating memory for input buffer.\n");

		 //Read data until buffer full
		 while(bytes_read < filestat.st_size) 
		 {
				r = read(fd, preloading_buffer+bytes_read, filestat.st_size-bytes_read);
				errorNumber(r);
				if(r == 0) 		break;
				bytes_read += r;
	  }
	  data_process_args.input_file.size = filestat.st_size;
	  data_process_args.input_file.buffer = preloading_buffer;

	}

	//------------------SWAN-PIPELINE---------------------
	fprintf(stderr, "before SwanPipeline!!\n");
	struct thread_args *dpargs = &data_process_args;
	stimer_t tmr_pipeline;
	stimer_tick(&tmr_pipeline);
	run(SwanPipeline, dpargs);
	stimer_tuck(&tmr_pipeline, "Pipeline took: ");
	//SerialIntegratedPipeline(&data_process_args);

	//------------------END SWAN---------------------------
	
	//clean up after preloading
	if(_conf->preloading)					 free(preloading_buffer);

	/* clean up with the src file */
	//	  if (_conf->infile != NULL)			 
	close(fd);
	fprintf(stderr, "file size = %lf \n", thread_stats->total_input);
	fprintf(stderr, "confi->infile = %s\n", _conf->infile);
	assert(!mbuffer_system_destroy());
	hashtable_destroy(cache, TRUE);

	#ifdef ENABLE_STATISTICS
	/* dest file stat */
	//fprintf(stderr, "confi->outfile = %s\n", confi->outfile);
	//if (stat(_conf->outfile, &filestat) < 0) 
	if (stat("output.dat.ddp", &filestat) < 0) 
		EXIT_TRACE("stat() %s failed: %s\n", _conf->outfile, strerror(errno));
	thread_stats->total_output = filestat.st_size;

	//Analyze and print statistics
	//if(confi->verbose) 
	print_stats(thread_stats);
	#endif //ENABLE_STATISTICS
}




