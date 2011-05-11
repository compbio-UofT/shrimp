#ifndef __MERGESAM_H__
#define __MERGESAM_H__
#define DEF_MAX_ALIGNMENTS	-1
#define DEF_MAX_OUTPUTS	20
#define DEF_INSERT_SIZE -1
#define GROWTH_FACTOR 1.3
#define SIZE_READ_NAME 255
#define MIN(a,b) (((a) < (b)) ? (a) : (b))
#define MAX_INT32	2147483647

#define PAIRED		0
#define UNPAIRED	1
#define FIRST_LEG	2
#define SECOND_LEG	3
#define UNMAPPED	4
#define LL_ALL		5

//IO SETTINGS
#define DEF_READ_SIZE	1024*10
#define DEF_BUFFER_SIZE	1024*1024*50
#define DEF_READ_RATE 3000
#define DEF_ALIGNMENTS_STACK_SIZE 21
typedef struct runtime_options {
	//input options
	size_t buffer_size;
	size_t read_size;
	size_t alignments_stack_size;
	int read_rate;
	bool fastq;
	bool colour_space;
	bool letter_space;
	int expected_insert_size;
	int threads;
	int number_of_sam_files;
	//output options
	bool strata;
	bool half_paired;
	bool sam_unaligned;
	bool sam_format;
	int32_t max_alignments;
	int32_t max_outputs;
	//determined at runtime options
	bool paired;
	bool unpaired; 
};
extern runtime_options options;
#endif
