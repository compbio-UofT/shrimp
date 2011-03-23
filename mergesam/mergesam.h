#ifndef __MERGESAM_H__
#define __MERGESAM_H__
typedef struct runtime_options {
	//input options
	size_t buffer_size;
	size_t read_size;
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
