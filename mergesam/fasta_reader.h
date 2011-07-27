#ifndef __FASTA_READER_PARSERS__
#define __FASTA_READER_PARSERS__

#include "file_buffer.h"

typedef struct fasta_reader {
	size_t reads_seen;
	size_t reads_unseen;
	size_t reads_filled; //reads_seen+reads_unseen = reads_filled
	size_t reads_inmem; //how many reads is space allocated for

	bool exhausted;
	file_buffer * fb;
	int i;
};
void fasta_last_entry(char * base,size_t used, char ** start,char ** end);
void parse_fasta(fasta_reader * fard);
size_t fasta_move(file_buffer * fb, char * dest);
#endif
