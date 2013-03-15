/* AUTORIGHTS
Copyright (C) 2007 Princeton University
      
This file is part of Ferret Toolkit.

Ferret Toolkit is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2, or (at your option)
any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software Foundation,
Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
*/
#include <stdio.h>
#include <math.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <dirent.h>
#include <unistd.h>
#include <pthread.h>
#include <iostream>
#include <list>

extern "C"
{
#include <cass.h>
#include <cass_timer.h>
#include <../image/image.h>
// #include "tpool.h"
// #include "queue.h"
}

#include "wf_interface.h"
using namespace obj;

#ifdef ENABLE_PARSEC_HOOKS
#include <hooks.h>
#endif

#define DEFAULT_DEPTH	25
#define MAXR	100
#define IMAGE_DIM	14

const char *db_dir = NULL;
const char *table_name = NULL;
const char *query_dir = NULL;
const char *output_path = NULL;

FILE *fout;

int NTHREAD_LOAD = 1;
int NTHREAD_SEG	= 1;
int NTHREAD_EXTRACT = 1;
int NTHREAD_VEC	= 1;
int NTHREAD_RANK = 1;
int NTHREAD_OUT	= 1;
int DEPTH = DEFAULT_DEPTH;

int top_K = 10;

char *extra_params = "-L 8 - T 20";

int input_end, output_end;
pthread_cond_t done;
pthread_mutex_t done_mutex;

cass_env_t *env;
cass_table_t *table;
cass_table_t *query_table;

int vec_dist_id = 0;
int vecset_dist_id = 0;

struct load_data
{
	/*int*/ long width, height;
	char *name;
	unsigned char *HSV, *RGB;
};

namespace platform_x86_64 {

template<size_t ireg, size_t freg, size_t loff>
struct arg_passing<ireg, freg, loff, load_data >
    : arg_passing_struct5<ireg, freg, loff, load_data, int, int, char *, unsigned char *, unsigned char *> {
};

}


struct seg_data
{
	int width, height, nrgn;
	char *name;
	unsigned char *mask;
	unsigned char *HSV;
};

struct extract_data
{
	cass_dataset_t ds;
	char *name;
};

struct vec_query_data
{
	char *name;
	cass_dataset_t *ds;
	cass_result_t result;
};

struct rank_data
{
	char *name;
	cass_dataset_t *ds;
	cass_result_t result;
};

/* ------- The Helper Functions ------- */
int cnt_enqueue;
int cnt_dequeue;
char path[BUFSIZ];

int scan_dir (const char *dir, char *head, std::list<load_data> * load_list);//__attribute__((always_inline));
void file_helper (const char *file, std::list<load_data>* load_list)__attribute__((always_inline));
struct load_data* load_array[5000];
unsigned int array_size = -1;

int dir_helper (char *dir, char *head, std::list<load_data> * load_list)
{
	DIR *pd = NULL;
	struct dirent *ent = NULL;
	volatile struct load_data ** retdata;
	int result = 0;
	pd = opendir(dir);
	if (pd == NULL) goto except;
	for (;;)
	{
		ent = readdir(pd);
		if (ent == NULL) break;
		if (scan_dir( ent->d_name, head, load_list) != 0) return -1;
	}
	goto final;

except:
	result = -1;
	perror("Error:");
final:
	if (pd != NULL) closedir(pd);
	return result;
}

/* the whole path to the file */
void file_helper (const char *file, std::list<load_data> * load_list)
{
	int r;
	struct load_data data;

	//data = (struct load_data *)malloc(sizeof(struct load_data));
	//assert(data != NULL);

	data.name = strdup(file);

	// r = image_read_rgb_hsv(file, &data.width, &data.height, &data.RGB, &data.HSV);
	int w, h;
	r = image_read_rgb_hsv(file, &w, &h, &data.RGB, &data.HSV);
	data.width = w;
	data.height = h;
	assert(r == 0);

	/*
		r = image_read_rgb(file, &data->width, &data->height, &data->RGB);
		r = image_read_hsv(file, &data->width, &data->height, &data->HSV);
		*/

	cnt_enqueue++;
	// enqueue(&q_load_seg, data);
	load_list->push_back( data );
}

int scan_dir (const char *dir, char *head, std::list<load_data> * load_list)
{
	struct stat st;
	int ret;
	/* test for . and .. */
	if (dir[0] == '.')
	{
		if (dir[1] == 0) return 0;
		else if (dir[1] == '.')
		{
			if (dir[2] == 0) return 0;
		}
	}

	/* append the name to the path */
	strcat(head, dir);
	ret = stat(path, &st);
	if (ret != 0)
	{
		perror("Error:");
		return -1;
	}
	if (S_ISREG(st.st_mode)) file_helper(path, load_list);
	else if (S_ISDIR(st.st_mode))
	{
		strcat(head, "/");
		dir_helper(path, head + strlen(head), load_list);
	}
	/* removed the appended part */
	head[0] = 0;
	return 0;
}


/* ------ The Stages ------ */
void t_load (const char * dir, std::list<load_data> * load_list)
{
	path[0] = 0;

	if (strcmp(dir, ".") == 0)
	{
	    dir_helper(".", path, load_list);
	}
	else
	{
		scan_dir(dir, path, load_list);
	}
}


void t_seg (load_data * load, seg_data* seg)
{
	seg->name = load->name;
	seg->width = load->width;
	seg->height = load->height;
	seg->HSV = load->HSV;
	//image segment is inline
	if(image_segment( (void**)&seg->mask, &seg->nrgn, (uchar*)load->RGB, load->width, load->height)!=0)
		fprintf(fout, "wrong image_segment return!!\n");

	free(load->RGB);
	
	return;
}

// void t_seg_wrapper(popdep<load_data> load_queue, outdep<seg_data> segarg)
void t_seg_wrapper(load_data load, outdep<seg_data> segarg)
{
	// load_data load = load_queue.pop();
	// fprintf(stderr, "pop: %s\n", load.name );
	seg_data * d = (seg_data*)segarg;
	leaf_call(t_seg, &load, d);
	return;
}

void t_extract (seg_data* seg, extract_data* extract)
{
	assert(seg != NULL);
	extract->name = seg->name;							
	image_extract_helper( seg->HSV, seg->mask, seg->width, seg->height, seg->nrgn, &extract->ds);

	free(seg->mask);
	free(seg->HSV);
}

void t_extract_wrapper(indep<seg_data> segarg, outdep<extract_data> extractarg)
{
	seg_data* seg = (seg_data*)segarg;
	extract_data * x = (extract_data*)extractarg;
	leaf_call(t_extract, seg, x);
}

void t_vec (struct extract_data* extract, struct vec_query_data* vec)
{
	cass_query_t query;
	assert(extract != NULL);
	vec->name = extract->name;
	memset(&query, 0, sizeof query);
	query.flags = CASS_RESULT_LISTS | CASS_RESULT_USERMEM;
	vec->ds = query.dataset = &extract->ds;
	query.vecset_id = 0;
	query.vec_dist_id = vec_dist_id;
	query.vecset_dist_id = vecset_dist_id;
	query.topk = 2*top_K;
	query.extra_params = extra_params;
	//cass_result_alloc_list, cass_table_query are inline
	cass_result_alloc_list(&vec->result, vec->ds->vecset[0].num_regions, query.topk);
	cass_table_query( table, &query, &vec->result);
	//vecarg = vec;
	

	return;
}

void t_vec_wrapper(indep<extract_data> extractarg, outdep<struct vec_query_data> vecarg)
{
	extract_data* extract = (extract_data*)extractarg;
	vec_query_data * query = (vec_query_data*)vecarg;
	leaf_call(t_vec, extract, query);
	return;
}

void t_rank (struct vec_query_data* vec, struct rank_data* rank)
{
	cass_result_t *candidate;
	cass_query_t query;
	rank->name = vec->name;
	query.flags = CASS_RESULT_LIST | CASS_RESULT_USERMEM | CASS_RESULT_SORT;
	query.dataset = vec->ds; // this points to extract_data structure
	query.vecset_id = 0;
	query.vec_dist_id = vec_dist_id;
	query.vecset_dist_id = vecset_dist_id;
	query.topk = top_K;
	query.extra_params = NULL;
	//cass_result_merge_lists is inline
	candidate = cass_result_merge_lists(&vec->result, (cass_dataset_t *)query_table->__private, 0);
	query.candidate = candidate;

	//cass_result_alloc_list, cass_result_free, cass_dataset_release, cass_table_query are inline
	cass_result_alloc_list((cass_result_t *)&rank->result, (unsigned int)0, (int)top_K);	
	cass_table_query(query_table, &query, &rank->result);
	cass_result_free(&vec->result);
	cass_result_free(candidate);
	free(candidate);
	cass_dataset_release(vec->ds);
	// free(vec->ds); -- this is an error and only works by accident/bad design in pthreads version!
	
	return;
}

// Add extractarg here because vecarg stores a pointer to extractarg
void t_rank_wrapper(indep<vec_query_data> vecarg, indep<extract_data> extractarg, outdep<rank_data> rankarg)
{
	vec_query_data* vec = (vec_query_data*)vecarg;
	rank_data * rank = (rank_data*)rankarg;
	leaf_call(t_rank, vec, rank);
	return;
}

void t_sxvrs(load_data ld, outdep<rank_data> rank_obj) {
	load_data * load = &ld;
	seg_data seg;
	rank_data * rank = (rank_data*)rank_obj;
	extract_data x;
	vec_query_data vq;

	leaf_call( t_seg, load, &seg );
	leaf_call( t_extract, &seg, &x );
	leaf_call( t_vec, &x, &vq );
	leaf_call( t_rank, &vq, rank );
}

void t_out (struct rank_data* rank);

void t_out_wrapper(indep<struct rank_data> rankarg, inoutdep<int> wait)
{
	leaf_call( t_out, (rank_data*)rankarg );
}

void t_sxvr(std::list<load_data> * load_list) {
	object_t<rank_data> rank_obj;
	object_t<int> serialize;
	object_t<seg_data> seg;
	object_t<extract_data> x;
	object_t<vec_query_data> vq;
	while( !load_list->empty() ) {
		load_data ld = load_list->front();
		load_list->pop_front();
		// spawn(t_sxvrs, ld, (outdep<rank_data>)rank_obj);
		spawn( t_seg_wrapper, ld, (outdep<seg_data>)seg );
		spawn( t_extract_wrapper, (indep<seg_data>)seg, (outdep<extract_data>)x );
		spawn( t_vec_wrapper, (indep<extract_data>)x, (outdep<vec_query_data>)vq );
		spawn( t_rank_wrapper, (indep<vec_query_data>)vq, (indep<extract_data>)x, (outdep<rank_data>)rank_obj );
		spawn(t_out_wrapper, (indep<rank_data>)rank_obj, (inoutdep<int>)serialize);
	}
	ssync();
}

void t_xvr(indep<seg_data> segdataobj, outdep<rank_data> rankarg) {
	seg_data * seg = (seg_data *)segdataobj;
	rank_data * rank = (rank_data*)rankarg;
	extract_data x;
	vec_query_data vq;
	leaf_call( t_extract, seg, &x );
	leaf_call( t_vec, &x, &vq );
	leaf_call( t_rank, &vq, rank );
}

void t_out (struct rank_data* rank)
{
	//struct rank_data *rank;
	//rank = rankarg;
	assert(rank != NULL);
	fprintf(fout, "%s", rank->name);
	ARRAY_BEGIN_FOREACH(rank->result.u.list, cass_list_entry_t p)
	{
		char *obj = NULL;
		if (p.dist == HUGE) continue;
		//cass_map_id_to_dataobj is inline
		cass_map_id_to_dataobj(query_table->map, p.id, &obj);
		assert(obj != NULL);
		fprintf(fout, "\t%s:%g", obj, p.dist);
	} ARRAY_END_FOREACH;

	fprintf(fout, "\n");

	//cass_result_free is inline
	cass_result_free(&rank->result);
	free(rank->name);
	cnt_dequeue++;
		
	// fprintf(stderr, "(%d,%d)\n", rank->count, cnt_dequeue);
	fprintf(stderr, "(%d,%d)\n", cnt_enqueue, cnt_dequeue);

/*	if (input_end && (cnt_enqueue == cnt_dequeue))
	{
		pthread_mutex_lock(&done_mutex);
		output_end = 1;
		pthread_cond_signal(&done);
		pthread_mutex_unlock(&done_mutex);
	}
*/	return;
}

using namespace std;

void pipeline()
{
//	const char* dummy;
	object_t<seg_data> segdataobj;
	object_t<extract_data> extractarg;
	object_t<vec_query_data> vecarg;
	object_t<rank_data> rankarg;
	unversioned<int>	wait;
/*	unversioned<seg_data> segdataobj;
	unversioned<extract_data> extractarg;
	unversioned<vec_query_data> vecarg;
	unversioned<rank_data> rankarg;
*/
	std::list<load_data> load_list;

	call(t_load, query_dir, &load_list);
	call(t_sxvr, &load_list );

	assert(cnt_enqueue == cnt_dequeue);
}

int main (int argc, char *argv[])
{
	stimer_t tmr;

	int ret, i;
	// fprintf(stderr, "main starting...\n");
#ifdef PARSEC_VERSION
#define __PARSEC_STRING(x) #x
#define __PARSEC_XSTRING(x) __PARSEC_STRING(x)
        printf("PARSEC Benchmark Suite Version "__PARSEC_XSTRING(PARSEC_VERSION)"\n");
        fflush(NULL);
#else
        printf("PARSEC Benchmark Suite\n");
        fflush(NULL);
#endif //PARSEC_VERSION
#ifdef ENABLE_PARSEC_HOOKS
	__parsec_bench_begin(__parsec_ferret);
#endif

	cnt_enqueue = cnt_dequeue = 0;

	fprintf(stderr, "Running with pipeline model\n");

	if (argc < 8)
	{
		printf("%s <database> <table> <query dir> <top K> <depth> <n> <out>\n", argv[0]); 
		return 0;
	}

	db_dir = argv[1];
	table_name = argv[2];
	query_dir = argv[3];
	top_K = atoi(argv[4]);

	DEPTH = atoi(argv[5]);
	NTHREAD_SEG = atoi(argv[6]);
	NTHREAD_EXTRACT = atoi(argv[6]);
	NTHREAD_VEC = atoi(argv[6]);
	NTHREAD_RANK = atoi(argv[6]);

	output_path = argv[7];

	fout = fopen(output_path, "w+");
	assert(fout != NULL);

	cass_init();

	ret = cass_env_open(&env, db_dir, 0);
	if (ret != 0) { printf("ERROR: %s\n", cass_strerror(ret)); return 0; }

	vec_dist_id = cass_reg_lookup(&env->vec_dist, "L2_float");
	assert(vec_dist_id >= 0);

	vecset_dist_id = cass_reg_lookup(&env->vecset_dist, "emd");
	assert(vecset_dist_id >= 0);

	i = cass_reg_lookup(&env->table, table_name);


	table = query_table = cass_reg_get(&env->table, i);

	i = table->parent_id;

	if (i >= 0)
	{
		query_table = cass_reg_get(&env->table, i);
	}

	if (query_table != table) cass_table_load(query_table);
	
	cass_map_load(query_table->map);

	cass_table_load(table);

	image_init(argv[0]);
	
	
	stimer_tick(&tmr);
	run(pipeline);
	stimer_tuck(&tmr, "QUERY TIME");
	ret = cass_env_close(env, 0);
	if (ret != 0) { printf("ERROR: %s\n", cass_strerror(ret)); return 0; }

	cass_cleanup();

	image_cleanup();

	fclose(fout);

#ifdef ENABLE_PARSEC_HOOKS
	__parsec_bench_end();
#endif
	return 0;
}

