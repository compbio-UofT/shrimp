#ifndef __MERGESAM_H__
#define __MERGESAM_H__


#define DEF_MAX_OUTPUTS	10

#define DEF_INSERT_SIZE -1

#define GROWTH_FACTOR 1.7
#define SIZE_READ_NAME 255
//#define MIN(a,b) (((a) < (b)) ? (a) : (b))
#define MAX_INT32	2147483647

#define PAIRED	0
#define PAIRED_SECOND	1
#define UNPAIRED	2
#define FIRST_LEG	3
#define SECOND_LEG	4
#define UNMAPPED	5
#define LL_ALL		6

//IO SETTINGS
#define DEF_READ_SIZE	1024*10
#define DEF_BUFFER_SIZE	1024*1024*100
#define DEF_READ_RATE 40000
#define DEF_ALIGNMENTS_STACK_SIZE DEF_READ_RATE*2
typedef struct runtime_options {
	//efficiency options
	size_t buffer_size;
	size_t read_size;
	size_t alignments_stack_size;
	int read_rate;
	int threads;
	int min_mapq;
	//gmapper supported options
	bool fastq;
	bool fastq_set;
	bool mode_set;
	bool colour_space;
	bool letter_space;
	double insert_size_mean;
	double insert_size_stddev;
	bool strata;
	bool half_paired;
	bool sam_unaligned;
	bool sam_format;
	int32_t max_alignments;
	int32_t max_outputs;
	FILE * unaligned_reads_file;
	FILE * aligned_reads_file;
	
	bool no_mapping_qualities;
	bool no_improper_mappings;
	bool no_autodetect_input;
	//mergesam specific
	int number_of_sam_files;
	bool single_best;
	bool all_contigs;
	//determined at runtime options
	bool paired;
	bool unpaired; 
	bool leave_mapq;
};
extern runtime_options options;

typedef struct output_buffer {
	size_t size;
	size_t used;
	char * base;
};

extern int64_t genome_length;
#endif
