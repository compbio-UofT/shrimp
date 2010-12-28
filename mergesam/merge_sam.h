#ifndef __MERGE_SAM_H__
#define __MERGE_SAM_H__
#include <stdio.h>
#include <stdint.h>
#include <stdbool.h>
#include <zlib.h>

#define	DEF_FILE_BUFFER_SIZE	1024*1024*40

typedef struct file_buffer file_buffer;
struct file_buffer {
        size_t size;
        size_t filled;
        size_t left_over;
	char * buffer;
	gzFile file;
	uint64_t last_read;
	int eof;
};

typedef struct output_filter output_filter;
struct output_filter {
	int number_outputs;
	int max_alignments;
	bool strata;
	bool unaligned;
	bool half_paired;	
};

typedef struct mapq_info mapq_info;
struct mapq_info {
	bool calculate;
	double score_matrix_lambda;
	double sw_scale_constant;
	int top_hits;
	bool strata;
};

void merge_sam(char* reads_filename, bool reads_fastq, bool colour_space, int read_threads, int compute_threads, int number_of_files, char ** filenames, FILE * output_file,output_filter of,mapq_info mqi,int buffer_size);
#endif 
