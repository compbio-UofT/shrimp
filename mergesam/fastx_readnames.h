#ifndef __FILE_BUFFER_PARSERS__
#define __FILE_BUFFER_PARSERS__
#include "file_buffer.h"
#include "sam2pretty_lib.h"
#include "mergesam.h"
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
#define DEF_READ_SIZE	1024*1024*5
#define DEF_BUFFER_SIZE	1024*1024*15
#define DEF_READ_RATE 10
typedef struct fastx_readnames {
	size_t reads_seen;
	size_t reads_unseen;
	size_t reads_filled; //reads_seen+reads_unseen = reads_filled
	
	size_t reads_inmem; //how many reads is space allocated for

	char * read_names;
	
	bool reads_exhausted;	
	file_buffer * fb;	
};
void reconsolidate_reads(fastx_readnames  * fxrn);
void parse_reads(fastx_readnames * fxrn);
#endif
