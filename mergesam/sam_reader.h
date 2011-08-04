#include "file_buffer.h"
#include "fastx_readnames.h"
#include "sam2pretty_lib.h"
#include "mergesam_heap.h"
#include "mergesam.h"

#ifndef __SAM_READER__
#define __SAM_READER__

typedef struct pp_ll pp_ll;
struct pp_ll {
	pretty * head;
	pretty * tail;
	size_t length;
	double * file_z1s;
};
struct sam_reader {
	file_buffer * fb;
	pretty * pretty_stack;
	size_t pretty_stack_size;
	size_t pretty_stack_end;
	size_t pretty_stack_start;
	pp_ll * pp_lls;
	pp_ll * sam_headers;
	size_t last_tested;
	size_t * inter_offsets;
	size_t * pretty_stack_ends;
	int fileno;
};
void parse_sam(sam_reader * sr,fastx_readnames * fxrn);
void pp_ll_combine_and_check(pp_ll * m_ll,pp_ll ** ll,heap_pa * h,output_buffer * ob);
void grow_sam_pretty(sam_reader * sr);
int sam_header_sort(const void * a, const void *b);
void sam_close(sam_reader * sr);
sam_reader * sam_open(char * sam_filename,fastx_readnames * fxrn);
extern bool found_sam_headers;
#endif
