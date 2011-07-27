#ifndef __FILE_BUFFER_PARSERS__
#define __FILE_BUFFER_PARSERS__

#include "file_buffer.h"
#include "sam2pretty_lib.h"
#include "mergesam.h"

typedef struct fastx_readnames {
	size_t reads_seen;
	size_t reads_unseen;
	size_t reads_filled; //reads_seen+reads_unseen = reads_filled
	
	size_t reads_inmem; //how many reads is space allocated for

	char * read_names;
	
	bool reads_exhausted;	
	file_buffer * fb;	
	bool fastq_seen_name;//=false;
	bool fastq_seen_plus;//=false;
	int64_t fastq_seen_seq;//=0;
	int64_t fastq_seen_qual;//=0;  
};


void reconsolidate_reads(fastx_readnames  * fxrn);
void parse_reads(fastx_readnames * fxrn);
#endif
