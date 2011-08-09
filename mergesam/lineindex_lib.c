#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <omp.h>
#include "lineindex_lib.h"
#include <assert.h>

#define DEF_LINEINDEX_SIZE 1000


void lineindex_resize(lineindex_table * lineindex, size_t size) {
	assert(lineindex->size<size);
	char ** ptr = (char**)malloc(sizeof(char**)*size);
	if (lineindex->size>0) {
		size_t start_mod = lineindex->start%lineindex->size;
		size_t end_mod = lineindex->end%lineindex->size;
		if (start_mod>=end_mod && lineindex->start!=lineindex->end) {
			size_t used=0;
			size_t entries_to_move=lineindex->size-start_mod;
			memmove(ptr,lineindex->table+start_mod,entries_to_move*sizeof(char*));
			used+=entries_to_move;
			entries_to_move=end_mod;
			memmove(ptr+used,lineindex->table,entries_to_move*sizeof(char*));
			used+=entries_to_move;
			lineindex->end=used;
			lineindex->start=0;
		} else {
			//end_mod > start_mod
			size_t entries_to_move=end_mod-start_mod;
			memmove(ptr,lineindex->table+start_mod,entries_to_move*sizeof(char*));
			lineindex->start=0;
			lineindex->end=entries_to_move;
		}
		free(lineindex->table);
	}
	lineindex->size=size;
	lineindex->table=ptr;
	return;
}

static inline void lineindex_resize_if_needed(lineindex_table * lineindex) {
	if (lineindex->start+lineindex->size==lineindex->end) {
		//fprintf(stderr,"Resizing!\n");
		lineindex_resize(lineindex,lineindex->size*2);
	}
	return;
}
lineindex_table * lineindex_init(size_t size) {
	lineindex_table * lineindex=(lineindex_table*)malloc(sizeof(lineindex_table));
	if (lineindex==NULL) {
		fprintf(stderr,"ran out of memory!\n");	
		exit(1);
	}
	memset(lineindex,0,sizeof(lineindex_table));
	lineindex_resize(lineindex,size);
	return lineindex;	
}

void lineindex_destroy(lineindex_table * lt) {
	assert(lt!=NULL);
	if (lt->size>0) {
		free(lt->table);
	}
	free(lt);
}

//does not add first as a line!\n
lineindex_table * add_lineindex_from_memory(lineindex_table * lineindex, char * data, size_t size, char comment_start,bool has_extra_byte) {
	size_t i;
	for (i=0; i<size; i++) {
		if (data[i]=='\n') {
			if (i+1<size || has_extra_byte) {
				if (data[i+1]!=comment_start && data[i+1]!='\0') {
					lineindex_resize_if_needed(lineindex);
					//fprintf(stderr,"ENTRY : %d | %d\n",(lineindex->end)%lineindex->size,lineindex->size);
					lineindex->table[(lineindex->end++)%lineindex->size]=i+data+1;
				}
			}
			data[i]='\0';
		}	
	}
	return lineindex;	
}

size_t add_lineindex_from_memory_threaded(lineindex_table * lineindex, lineindex_table ** thread_lineindexes,char *data, size_t data_size, int threads, char comment_start) {
	if (data_size==0) {
		return 0;
	}
	if (threads>data_size) {
		threads=(int)data_size;
	}
	//set up the first entry in the file	
	if (data_size>0 && data[0]!='\0' && data[0]!=comment_start) {
		lineindex_resize_if_needed(lineindex);
		lineindex->table[(lineindex->end++)%lineindex->size]=data;
	}
	//offsets of where to copy into
	size_t offsets[threads];
	//populate
	//fprintf(stderr,"%d LINE INDEX IS %d\n",lineindex->size,lineindex->end);
	omp_set_num_threads(threads);
	#pragma omp parallel
	{
		size_t i = (size_t)omp_get_thread_num();
		//fprintf(stderr,"THIS THREAD %d\n",i);
		size_t thread_chunk_size = data_size/threads + ((data_size%threads)>i ? 1 : 0);
		size_t offset = (data_size/threads)*i + ((data_size%threads)>i ?  i : (data_size%threads) );
		//fprintf(stderr,"%lu: Starting at %lu, size %lu\n",i,offset,thread_chunk_size);
		thread_lineindexes[i]->start=0; thread_lineindexes[i]->end=0;
		add_lineindex_from_memory(thread_lineindexes[i],data+offset,thread_chunk_size,comment_start,(offset+thread_chunk_size)!=data_size);
		//int j;
		//fprintf(stderr,"START %d end %d\n",thread_lineindexes[i]->start,thread_lineindexes[i]->end);
		//for (j=thread_lineindexes[i]->start; thread_lineindexes[i]->end>j; j++ ) {
		//	fprintf(stderr,"RETTT %d |%s|\n",j,thread_lineindexes[i]->table[j]);
		//}
	}
	//fprintf(stderr,"%d LINE INDEX IS %d,%d\n",lineindex->size,lineindex->end,lineindex->start);
	//merge
	size_t newly_added=0;
	int i;
	for (i=0; i<threads; i++) {
		offsets[i]=(i>0 ? offsets[i-1]+thread_lineindexes[i-1]->end-thread_lineindexes[i-1]->start : 0);
		newly_added+=thread_lineindexes[i]->end-thread_lineindexes[i]->start;
	}
	//fprintf(stderr,"%d LINE INDEX IS %d\n",lineindex->size,lineindex->end);
	size_t new_used = newly_added+lineindex->end-lineindex->start;
	while (new_used>lineindex->size) {
		lineindex_resize(lineindex,lineindex->size*2);				
		//fprintf(stderr,"X %d LINE INDEX IS %d\n",lineindex->size,lineindex->end);
	}
	//fprintf(stderr,"%d LINE INDEX IS %d\n",lineindex->size,lineindex->end);
	#pragma omp parallel
	{
		int i = omp_get_thread_num();
		assert(thread_lineindexes[i]->start==0);
		size_t thread_index_size=thread_lineindexes[i]->end-thread_lineindexes[i]->start;
		if (thread_index_size>0) {
			//fprintf(stderr,"lineindex->end %d, offsets[i] %d\n",lineindex->end,offsets[i]);
			size_t dest_start_mod_inclusive=(lineindex->end+offsets[i])%lineindex->size;
			size_t dest_end_mod_inclusive=(lineindex->end+offsets[i]+thread_index_size-1)%lineindex->size;	
			//fprintf(stderr,"THREAD MOVE : %d %lu, %lu %lu\n",i,thread_index_size,dest_start_mod_inclusive,dest_end_mod_inclusive);
			if (dest_end_mod_inclusive<dest_start_mod_inclusive) {
				//need to do two move operations
				size_t to_move_to_end=lineindex->size-dest_start_mod_inclusive;
				memmove(lineindex->table+dest_start_mod_inclusive,thread_lineindexes[i]->table,sizeof(char*)*to_move_to_end);
				memmove(lineindex->table,thread_lineindexes[i]->table+(lineindex->size-dest_start_mod_inclusive),sizeof(char*)*(thread_index_size-to_move_to_end));
			} else {
				//just one move operation
				memmove(lineindex->table+dest_start_mod_inclusive,thread_lineindexes[i]->table,sizeof(char*)*thread_index_size);
			}
		}
	}
	lineindex->end=new_used+lineindex->start;
	return newly_added;
}

