#include <stdint.h>
#include <stdbool.h>
#include <stdlib.h>
#include <stdio.h>
#ifndef __SAM2PRETTY_H__
#define __SAM2PRETTY_H__

#define SIZE_NEWLINE 1
#define SIZE_TAB 1
#define SIZE_NULL 1
#define SIZE_SAM_AUX 5
#define SIZE_FLAG 5
#define SIZE_POS 9
#define SIZE_MAPQ 3
#define SIZE_MPOS 9
#define SIZE_ISIZE 10
#define SIZE_32bit 10
#define SIZE_DOUBLE 12

#define SAM2PRETTY_NUM_ZS	7

#define HAS_Z0	1<<0
#define HAS_Z1	1<<1
#define HAS_Z2	1<<2
#define HAS_Z3	1<<3
#define HAS_Z4	1<<4
#define HAS_Z5	1<<5
#define HAS_Z6	1<<6

#define HAS_ZPAIRED	(HAS_Z2 | HAS_Z3 | HAS_Z4 | HAS_Z6)
#define HAS_ZHALF	(HAS_Z0 | HAS_Z1 | HAS_Z4 | HAS_Z5)
#define HAS_ZUNPAIRED	(HAS_Z0 | HAS_Z1)

typedef struct pretty pretty;
struct pretty {
	//char* genome_string;
	char* read_string;
	char* read_qualities;
	char* pretty_genome_string;
	char* pretty_read_string;
	char* pretty_read_qualities;
	char* pretty_cs_string;
	char* pretty_cs_qualities;
	char* cigar;
	char* reference_name;
	char* mate_reference_name;
	char* pretty_match_string;
	int32_t num_cigar;
	char* cigar_ops;
	uint32_t * cigar_lengths;
	//for colour space
	bool colour_space;
	uint32_t genome_start_padded;
	uint32_t genome_end_padded;
	uint32_t genome_start_unpadded;
	uint32_t mate_genome_start_unpadded;
	uint32_t genome_end_unpadded;
	int32_t clipped_read_start;
	int32_t clipped_read_end;
	int32_t pretty_clipped_read_start;
	int32_t pretty_clipped_read_end;

	char* genome_sequence;
	int32_t genome_length;

	uint32_t read_id;
	char * read_name;
	size_t read_name_length;
	int32_t isize;

	int fileno;
	int has_zs;
	double z[SAM2PRETTY_NUM_ZS];

	//int32_t strand; //0 is postive, 1 is reverse
	int32_t mapq;
	double temp_mapq;
	int32_t clipped;
	int32_t mismatches;
	int32_t deletions;
	int32_t insertions;
	int32_t matches;
	int32_t skipped;

	//flags
	int32_t flags;
	//0x0001 the read is paired in sequencing, no matter whether it is mapped in a pair
	bool paired_sequencing;
	//0x0002 the read is mapped in a proper pair (depends on the protocol, normally inferred during alignment) 1
	bool proper_pair;
	//0x0004 the query sequence itself is unmapped
	bool mapped;
	//0x0008 the mate is unmapped 1
	bool mp_mapped;
	//0x0010 strand of the query (0 for forward; 1 for reverse strand)
	bool reverse;
	//0x0020 strand of the mate 1
	bool mp_reverse;
	//0x0040 the read is the first read in a pair 1,2
	bool first_in_pair;
	//0x0080 the read is the second read in a pair 1,2
	bool second_in_pair;
	//0x0100 the alignment is not primary (a read having split hits may have multiple primary alignment records)
	bool primary_alignment;
	//0x0200 the read fails platform/vendor quality checks
	bool platform_quality_fail;
	//0x0400 the read is either a PCR duplicate or an optical duplicate
	bool pcr_duplicate;

	int32_t unclipped_read_length;
	int32_t pretty_clipped_read_length;
	int32_t clipped_read_length;
	int32_t pretty_length;

	//aux fields
	bool has_cs_string;
	char* cs_string;
	bool has_cs_qualities;
	char* cs_qualities;
	bool has_edit_distance;
	int32_t edit_distance;
	bool has_cs_mismatches;
	int32_t cs_mismatches;
	bool has_cs_edit_string;
	char * cs_edit_string;
	bool has_read_group;
	char * read_group;
	bool has_score;
	int32_t score;
	bool has_ih;
	int32_t ih;
	bool has_hi;
	int32_t hi;
	bool has_r2;
	char * r2;
	bool has_h0;
	int32_t h0;
	bool has_h1;
	int32_t h1;
	bool has_h2;
	int32_t h2;
	//relative to read orientation
//	int32_t alignemtn_begins_with_insert;
	char * sam_string;
	size_t sam_string_length;
	pretty * mate_pair;
	pretty * next;
	bool mark;
	char * aux;

	bool sam_header;
};

extern bool sam2pretty_lib_verbose;

int32_t pretty_remap_header(FILE* f,char* contig_name,uint32_t offset_start, uint32_t offset_end);
int32_t pretty_header(FILE* f,char* contig_name,uint32_t sequence_size);
pretty * pretty_sub_prettys(pretty * parent);
pretty * pretty_sub_pretty(pretty * parent, int32_t from, int32_t to);
void pretty_genome(pretty* pa);
void pretty_ls(pretty* pa);
void pretty_cs(pretty* pa);
void pretty_qualities(pretty* pa);
void pretty_stats(pretty* pa);
//pretty *  pretty_from_sam(samfile_t * sam_file, bam1_t * entry);
pretty * pretty_from_string(char *sam_string);
pretty * pretty_from_string_inplace(char * sam_string,size_t length_of_string,pretty * pa);
void pretty_print(pretty* pa);
void pretty_free(pretty* pa);
void pretty_free_fast(pretty* pa);
void pretty_match(pretty* pa);
int32_t pretty_remap(pretty * pa, uint32_t offset_start, uint32_t offset_end);
void pretty_print_sam(FILE * f, pretty * pa);
void pretty_print_sam_force_ls(FILE * f, pretty * pa);
pretty * pretty_new();
void pretty_print_sam_update(pretty * pa,bool inplace);
void pretty_print_sam_unaligned(pretty * pa,bool inplace);
void pretty_from_aux_inplace(pretty * pa);
void pretty_print_sam_fastx(pretty* pa, bool inplace );
int pretty_get_flag(pretty * pa );

void calculate_genome_end(pretty * pa);
void calculate_insert_size(pretty * pa, pretty * pa_mp);
void pretty_cigar_parse(pretty * pa);
void fill_cigar_len(pretty *pa);

double inv_tnlog(int x);

int tnlog(double x);


#endif
