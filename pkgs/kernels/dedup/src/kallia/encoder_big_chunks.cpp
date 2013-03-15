/*
 * BIG CHUNKS, NO FRAGMENT REFINE
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

#include "swan/wf_interface.h"
#include "swan/swan_config.h"

//#define fprintf(stderr, ...) 
extern "C"
{
#include "util.h"
#include "dedupdef.h"
#include "debug.h"

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

#include <iostream>
#include <list>

using namespace obj;
using namespace std;
#define INITIAL_SEARCH_TREE_SIZE 4096

const int compress_type = COMPRESS_GZIP;
const int preloading = 1;

//The configuration block defined in main
config_t * confi;
//int maxchcount[MAX_NUM_ITERATIONS];
//Hash table data structure & utility functions
struct hashtable *cache;
#undef ENABLE_STATISTICS
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
// queue_t *deduplicate_que, *refine_que, *reorder_que, *compress_que;
#ifdef ENABLE_PTHREADS
//The queues between the pipeline stages
// queue_t *deduplicate_que, *refine_que, *reorder_que, *compress_que;

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
	//chunk_t *chunk = *chunk_ptr;
    assert(chunk!=NULL);
    //compress the item and add it to the database
//#ifdef ENABLE_PTHREADS
    //pthread_mutex_lock(&chunk->header.lock);
    assert(chunk->header.state == CHUNK_STATE_UNCOMPRESSED);
//#endif //ENABLE_PTHREADS
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
	(chunk)->header.state = CHUNK_STATE_COMPRESSED;

//#ifdef ENABLE_PTHREADS
    
    //pthread_cond_broadcast(&chunk->header.update);
   // pthread_mutex_unlock(&chunk->header.lock);
//#endif //ENABLE_PTHREADS

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
int sub_Deduplicate(chunk_t *chunk) {
  int isDuplicate;
  chunk_t *entry;
	//chunk_t *chunk = *chunk_ptr;

  assert(chunk!=NULL);
  assert(chunk->uncompressed_data.ptr!=NULL);

  SHA1_Digest((chunk)->uncompressed_data.ptr, (chunk)->uncompressed_data.n, (unsigned char *)((chunk)->sha1));

  //Query database to determine whether we've seen the data chunk before
#ifdef ENABLE_PTHREADS
  pthread_mutex_t *ht_lock = hashtable_getlock(cache, (void *)(chunk->sha1));
  pthread_mutex_lock(ht_lock);
#endif
	//chunk->sha1 = (void *)(*chunk_ptr)->sha1;
  entry = (chunk_t *)hashtable_search(cache, (void *)((chunk)->sha1));//chunk->sha1));
  isDuplicate = (entry != NULL);
  (chunk)->header.isDuplicate = isDuplicate;
  if (!isDuplicate) {
    // Cache miss: Create entry in hash table and forward data to compression stage
//#ifdef ENABLE_PTHREADS
    //pthread_mutex_init(&chunk->header.lock, NULL);
    //pthread_cond_init(&chunk->header.update, NULL);
//#endif
    //NOTE: chunk->compressed_data.buffer will be computed in compression stage
    if (hashtable_insert(cache, (void *)((chunk)->sha1), (void *)(chunk)) == 0) {
      EXIT_TRACE("hashtable_insert failed");
    }
  } else {
    // Cache hit: Skipping compression stage
    (chunk)->compressed_data_ref = entry;
    mbuffer_free(&chunk->uncompressed_data);
  }
#ifdef ENABLE_PTHREADS
  pthread_mutex_unlock(ht_lock);
#endif

  return isDuplicate;

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
cas_mutex output_mutex;
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
//void ReorderSwan(struct thread_args * targs, indep<list<chunk_t *>> chunks_obj, int fd, int priority)//indep<list<chunk_t *>> chunks_obj, int fd)

void ReorderSwan(struct thread_args * targs, indep<chunk_t*> chunk_obj, int fd, int iteration, inoutdep<int> wait)
{
	//fprintf(stderr, "priority = %d\n", priority);
	struct thread_args *args = (struct thread_args *)targs;
	chunk_t *chunk = chunk_obj;
	//list<chunk_t*> chunks = chunks_obj;
	//sequence_t next;
	//sequence_reset(&next);
	int r;
	int i;

	//list<chunk_t*>::iterator iter;
	wait = 0;
	//for(iter = chunks.begin(); iter != chunks.end(); iter++)
	{
		//chunk = *iter;
		if(chunk == NULL)
		{
			fprintf(stderr, "NULL item in the LIST!!!!!!\n");
			//output_mutex.unlock();
			exit(0);					
			//continue;
		}
		//fprintf(stderr, "In Reorder: chunk->l1num = %d, chunk->l2num = %d\n", chunk->sequence.l1num, chunk->sequence.l2num);
		
		leaf_call(write_chunk_to_file, fd, chunk);
	}
	return;
}


//stage 4 function
void CompressSwan(struct thread_args * targs, inoutdep<chunk_t *> chunk_obj)
{
	struct thread_args *args = (struct thread_args *)targs;
	//list<chunk_t *> chunks = chunk_obj;
	chunk_t *chunk = (chunk_t*) chunk_obj;
	//list<chunk_t*>::iterator iter;

	//for(iter = chunks.begin(); iter != chunks.end(); iter++)
	//{
		//chunk = *iter;
		if(!(chunk->header.isDuplicate))
		{
			leaf_call(sub_Compress, chunk);
			#ifdef ENABLE_STATISTICS
			thread_stats->total_compressed += (chunk)->compressed_data.n;
			#endif //ENABLE_STATISTICS
		}
		//else
			//continue;//return;
//	}
	chunk_obj = chunk;
	//put the item in the next queue for the write thread

}

//stage 3 function
void DeduplicateSwan(struct thread_args * targs, inoutdep<chunk_t*> chunk_obj, inoutdep<int> wait)
{
	struct thread_args *args = (struct thread_args *)targs;
	const int qid = args->tid / MAX_THREADS_PER_QUEUE;
	chunk_t *chunk = (chunk_t*)chunk_obj;// = 
	//chunk_t **chunk_ptr;
	int r;
	//list<chunk_t*> chunklist = chunklist_obj;
	assert(chunk!=NULL);
	//list<chunk_t*>::iterator iter;

	wait = 0;
			//for(iter = chunklist.begin(); iter != chunklist.end(); iter++)
			{
				//chunk = *iter;
				//chunk_ptr = &chunk;
				//Do the processing
				int isDuplicate =  leaf_call(sub_Deduplicate, chunk);
			#ifdef ENABLE_STATISTICS
				if(isDuplicate)  	thread_stats->nDuplicates++;
				else   				thread_stats->total_dedup += chunk->uncompressed_data.n;
			#endif //ENABLE_STATISTICS
				//Enqueue chunk either into compression queue or into send queue
			/*	if(!isDuplicate)
				{
					//compress = 1;
					//call compressSwan then return to call ReorderSwan
					//leaf_call(CompressSwan, args, chunk_ptr);
				}
				else 
				{
					//return to call ReorderSwan
					continue;
				}
*/			}
	//chunklist_obj = chunklist;
	chunk_obj = chunk;

}

//stage 2 function
void FragmentRefineSwan(struct thread_args * targs, 	indep<chunk_t*> chunk_obj, 
																		inoutdep<u32int*> rabintab_obj, 
																		inoutdep<u32int*> rabinwintab_obj, 
																		//chunk_t** chunks, 
																		outdep<std::list<chunk_t*>> chunkslist_obj,
																		int it)
{
	struct thread_args *args = (struct thread_args *)targs;
	const int qid = args->tid / MAX_THREADS_PER_QUEUE;
	int r = 0;

	chunk_t 	*temp;
	chunk_t	*chunk			=(chunk_t*)chunk_obj;
	u32int	*rabintab		=(u32int*)rabintab_obj,
			 	*rabinwintab	=(u32int*)rabinwintab_obj;
	chunk_t	*newchunk;
	list<chunk_t*> chunkslist = chunkslist_obj;

	rabininit(rf_win, rabintab, rabinwintab);
	int split;
	sequence_number_t chcount = 0;
	do
	{
		//Find next anchor with Rabin fingerprint
		int offset = rabinseg((&chunk->uncompressed_data)->ptr, (&chunk->uncompressed_data)->n, rf_win, rabintab, rabinwintab);
		//Can we split the buffer?
		newchunk = (chunk_t*)malloc(sizeof(chunk_t));
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

			newchunk->sequence.l2num = chcount;
			newchunk->isLastL2Chunk = FALSE;
			newchunk->header.state = chunk->header.state;
			newchunk->sequence.l1num = chunk->sequence.l1num;// = it;
			r = mbuffer_create(&newchunk->uncompressed_data, chunk->uncompressed_data.n);
			if(r!=0) 					EXIT_TRACE("Unable to allocate memory buffer...\n");
			memcpy(newchunk->uncompressed_data.ptr, chunk->uncompressed_data.ptr, chunk->uncompressed_data.n);
			newchunk->uncompressed_data.n = chunk->uncompressed_data.n;		//n is now offset
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
			chcount++;
			//prepare for next iteration
			chunk = NULL;
			split = 0;
		}
		chunkslist.push_back(newchunk);
	} while(split);
	chunkslist_obj = chunkslist;
	//maxchcount[it] = chcount;
	return;
	  //drain buffer
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
					exit(1);
					split = 0;
					*split_ptr = 0;
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
					*split_ptr = 1;
					return result;
				}
				else 
				{
					//Due to technical limitations we can't distinguish the cases "no split" and "split at end of buffer"
					//This will result in some unnecessary (and unlikely) work but yields the correct result eventually.
					*temp_ptr = chunk;
					*chunk_ptr = NULL;
					split = 0;
					*split_ptr = 0;
					continue;
				}
			} //if(ANCHOR_JUMP < chunk->uncompressed_data.n) 
			else 
			{
				//NOTE: We don't process the stub, instead we try to read in more data so we might be able to find a proper split.
				//Only once the end of the file is reached do we get a genuine stub which will be enqueued right after the read operation.
				*temp_ptr = chunk;
				*chunk_ptr = NULL;
				split = 0;
				*split_ptr = 0;
				return NULL;
				if(mode == 0)		*condition = 1;
			}
		}while(split);
	}

	return chunk;
}

//Stage 1 function
int FragmentSwan(outdep<chunk_t*> 		result_obj, 
						struct thread_args  *targs, 
						//inoutdep<chunk_t*> 	chunk_obj, 
						//inoutdep<chunk_t*> 	temp_obj, 
						chunk_t** chunk_ptr, chunk_t** temp_ptr,
						//outdep<int>			read_done_obj, 
						//int						*read_done_ptr,
						//inoutdep<size_t*>		bytes_left_obj, 
//						inoutdep<size_t*>		preloading_buffer_seek_obj, 
						size_t *					preloading_buffer_seek_ptr,
						//inoutdep<int*>			split_obj, 
						int *						split_ptr,
						//inoutdep<u32int*>		rabintab_obj, 
						u32int	*rabintab,
						//inoutdep<u32int*> 	rabinwintab_obj,
						u32int	*rabinwintab,
						int it)//amount of data left over in last_mbuffer from previous iteration
{
	struct thread_args *args = (struct thread_args *)targs;
	
	chunk_t 	*result 								= (chunk_t*)result_obj,//,		
				*chunk 								= (chunk_t*)*chunk_ptr,//chunk_obj, 			
				*temp 								= (chunk_t*)*temp_ptr;//temp_obj;
	//u32int	*rabintab 							= (u32int*)rabintab_obj,	
	//			*rabinwintab 						= (u32int*)rabinwintab_obj;
	//size_t 	//*bytes_left_ptr 					= (size_t*)bytes_left_obj,
			 	//*preloading_buffer_seek_ptr 	= (size_t*)preloading_buffer_seek_obj;
	//int 		*split_ptr 							= (int*)split_obj;

	size_t 	preloading_buffer_seek 			= *preloading_buffer_seek_ptr;
	int 		fd 									= args->fd;
	int 		i;	
	size_t 	bytes_left;// = *bytes_left_ptr;
	int 		condition 							= 0;
	sequence_number_t anchorcount 			= 0;
	//chunk_t 	**temp_ptr 							= &temp;
	//chunk_t 	**chunk_ptr = &chunk;
	int 		r, 
				split 								= *split_ptr;
	int read_done = 0;
	if(split)
	{
		result = FragmentSplit(split, split_ptr, &anchorcount, chunk_ptr, temp_ptr, rabintab, rabinwintab, &condition, 1);
		if(result != NULL)
		{
			result->sequence.l1num = it;
			result_obj = result;
			//chunk_obj = *chunk_ptr;
			//temp_obj = *temp_ptr;
			//rabintab_obj = rabintab;
			//rabinwintab_obj = rabinwintab;
			//split_obj = split_ptr;
			//bytes_left_obj = bytes_left_ptr;
			//preloading_buffer_seek_obj = preloading_buffer_seek_ptr;
			return read_done;
		}
	}
	temp = *temp_ptr;
	//read from input file / buffer
	while (1)
	{
		//Check how much data left over from previous iteration resp. create an initial chunk
		if(temp != NULL)
		{
			bytes_left = temp->uncompressed_data.n;
			//*bytes_left_ptr = bytes_left;
		}
		else 
		{
			bytes_left = 0;
			//*bytes_left_ptr = bytes_left;
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
		   //chunk->sequence.l1num = temp->sequence.l1num;

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
		   //chunk->sequence.l1num = anchorcount;
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
			result_obj = (chunk_t*)result;
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
			read_done = 1;//se auto to shmeio mporw na valw mia metavlhth read_done kai na thn 8etw =1
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
			result = chunk;
			break;
		}
		result  = FragmentSplit(/*split:*/1, split_ptr, &anchorcount, chunk_ptr, temp_ptr, rabintab, rabinwintab, &condition, 0);
		//partition input block into large, coarse-granular chunks
		if(condition)		continue;

		break;
	}		//while(1)
	result->sequence.l1num = it;
	result_obj = result;
	//chunk_obj = *chunk_ptr;
	//temp_obj = *temp_ptr;
	//rabintab_obj = rabintab;
	//rabinwintab_obj = rabinwintab;
	//split_obj = split_ptr;
	//bytes_left_obj = bytes_left_ptr;
	//preloading_buffer_seek_obj = preloading_buffer_seek_ptr;
	
	return read_done;
}

void terminate_function(indep<int*>	read_done_obj, int *read_done)
{
	int rd = *read_done_obj;
	*read_done = rd;
	return;
}

void storeChunks(indep<list<chunk_t*>> chunks, list<chunk_t*>allchunks[MAX_NUM_ITERATIONS], int iteration)
{

	allchunks[iteration] = chunks;
	return;
}

void SwanPipeline(struct thread_args *args)
{

	int 			split = 0, iteration = -1;
	//size_t 		bytes_left;// = (size_t *)malloc(sizeof(size_t));
	//chunk_t 		** temp_ptr = (chunk_t**)malloc(sizeof(chunk_t*));
	chunk_t 		*temp;// = *temp_ptr;
	size_t 		pr_buff_seek = 0;
	temp = NULL;	

	u32int * rabintab = (u32int*)malloc(256*sizeof rabintab[0]);
	u32int * rabinwintab = (u32int*)malloc(256*sizeof rabintab[0]);
	u32int * rabintabRef = (u32int*)malloc(256*sizeof rabintabRef[0]);
	u32int * rabinwintabRef = (u32int*)malloc(256*sizeof rabintabRef[0]);
	if(rabintab == NULL || rabinwintab == NULL || rabintabRef == NULL || rabinwintabRef == NULL)
		EXIT_TRACE("Memory allocation failed.\n");
	//thread_stats have been initiallized!

	//struct thread_args *args;
	//int i;
	chunk_t *chunk = NULL;//, *result = NULL;
	chunk_t **chunk_ptr = &chunk;
	//size_t *bytes_left_ptr = &bytes_left;
	size_t *pr_buff_seek_ptr = (size_t*)malloc(sizeof(size_t));
	*pr_buff_seek_ptr = pr_buff_seek;
	int *split_ptr = &split;
	split = 0;
	rf_win_dataprocess = 0;
	leaf_call(rabininit, rf_win_dataprocess, rabintab, rabinwintab);

	object_t<chunk_t*>	chunk_obj, temp_obj;//, result_obj;
	
	object_t<int*>		 	split_obj;
	object_t<size_t*>		bytes_left_obj, pr_buff_seek_obj;
	object_t<u32int*>		rabintab_obj, rabinwintab_obj,
								rabintabRef_obj, rabinwintabRef_obj;
	object_t<int>			wait;
	object_t<int>			wait2;
	wait2						= 0;
	wait 						= 0;

	
	//chunk_obj 		 		=(chunk_t*)chunk;
	//temp_obj 		 		=(chunk_t*)temp;
	//split_obj				=(int*)split_ptr;
	//bytes_left_obj			=(size_t*)bytes_left_ptr;
	pr_buff_seek_obj		=(size_t*)pr_buff_seek_ptr;
	//rabintab_obj			=(u32int*)rabintab;
	//rabinwintab_obj		=(u32int*)rabinwintab;
	rabintabRef_obj		=(u32int*)rabintabRef;
	rabinwintabRef_obj	=(u32int*)rabinwintabRef;

	int fd 					= create_output_file(confi->outfile);

	while(1)
	{
		int						read_done = 0;
		list<chunk_t*>			chunklist;
		//object_t<std::list<chunk_t*>>	chunks;
		//chunks = chunklist;

		chunk_t *result = NULL;
		object_t<chunk_t*>	result_obj;
		result_obj 		 		=(chunk_t*)result;

		iteration++;

		read_done = call(FragmentSwan, (outdep<chunk_t*>)result_obj, args, 
										//(inoutdep<chunk_t*>)chunk_obj, 
										&chunk, 
										&temp,
										//(inoutdep<chunk_t*>)temp_obj, 
										//(inoutdep<size_t*>)bytes_left_obj, 
										//(inoutdep<size_t*>)pr_buff_seek_obj, 
										pr_buff_seek_ptr,
										split_ptr,
										//(inoutdep<int*>)split_obj, 
										//(inoutdep<u32int*>)rabintab_obj, 
										rabintab,
										rabinwintab, iteration);
										//(inoutdep<u32int*>)rabinwintab_obj, iteration);

		if(!(read_done) || iteration == 0)
		{
			/*spawn(FragmentRefineSwan, args, 	(indep<chunk_t*>)result_obj, 
															(inoutdep<u32int*>)rabintabRef_obj, 
															(inoutdep<u32int*>)rabinwintabRef_obj, 
															(outdep<list<chunk_t*>>)chunks,
															iteration);									*/
			spawn(DeduplicateSwan, args, 	(inoutdep<chunk_t*>)result_obj, (inoutdep<int>) wait);
			spawn(CompressSwan, args, (inoutdep<chunk_t*>) result_obj);

			spawn(ReorderSwan, args, 	(indep<chunk_t*>)result_obj, fd, iteration, (inoutdep<int>) wait2);
		}
		else
			break;
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
void Encode(config_t * _conf)
{
	//fprintf(stderr, "Encode Starting.... conf-> inputfile = %s\n", _conf->infile);
	struct stat filestat;
	int32 fd;
	confi = _conf;
	//compress_type = _conf->compress_type;
	struct thread_args data_process_args;
	struct thread_args *data_process_args_ptr;
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
	//fprintf(stderr, "file size = %lf \n", thread_stats->total_input);
	//fprintf(stderr, "confi->infile = %s\n", _conf->infile);
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
	//fprintf(stderr, "before SwanPipeline!!\n");
	data_process_args_ptr = &data_process_args;
	stimer_t tmr_pipeline;
	stimer_tick(&tmr_pipeline);
	run(SwanPipeline, data_process_args_ptr);
	stimer_tuck(&tmr_pipeline, "Pipeline took: ");
	//SerialIntegratedPipeline(&data_process_args);

	//------------------END SWAN---------------------------
	
	//clean up after preloading
	if(_conf->preloading)					 free(preloading_buffer);

	/* clean up with the src file */
	//	  if (_conf->infile != NULL)			 
	close(fd);
	//fprintf(stderr, "file size = %lf \n", thread_stats->total_input);
	//fprintf(stderr, "confi->infile = %s\n", _conf->infile);
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

