#ifndef __SAM2PRETTY_H__
#define __SAM2PRETTY_H__

#include <stdint.h>
#include <stdbool.h>

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
	int num_cigar;
	char* cigar_ops;
	uint32_t* cigar_lengths;
	//for colour space
	bool colour_space;
	unsigned long genome_start_padded;
	unsigned long genome_end_padded;
	unsigned long genome_start_unpadded;
	unsigned long mate_genome_start_unpadded;
	unsigned long genome_end_unpadded;
	int clipped_read_start;
	int clipped_read_end;
	int pretty_clipped_read_start;
	int pretty_clipped_read_end;


	int id;


	char* genome_sequence;
	int genome_length;

	char * read_name;
	int isize;

	//int strand; //0 is postive, 1 is reverse
	int mapq;
	double temp_mapq;
	int clipped;
	int mismatches;
	int deletions;
	int insertions;
	int matches;
	int skipped;

	//flags
	int flags;
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

	int unclipped_read_length;
	int pretty_clipped_read_length;
	int clipped_read_length;
	int pretty_length;

	//aux fields
	bool has_cs_string;
	char* cs_string;
	bool has_cs_qualities;
	char* cs_qualities;
	bool has_edit_distance;
	int edit_distance;
	bool has_cs_mismatches;
	int cs_mismatches;
	bool has_cs_edit_string;
	char * cs_edit_string;
	bool has_read_group;
	char * read_group;
	bool has_score;
	int score;
	//relative to read orientation
//	int alignemtn_begins_with_insert;
	char * sam_string;
	pretty * mate_pair;
	pretty * next;
	bool mark;
};

int pretty_remap_header(FILE* f,char* contig_name,unsigned long offset_start, unsigned long offset_end);
int pretty_header(FILE* f,char* contig_name,unsigned long sequence_size);
pretty * pretty_sub_prettys(pretty * parent);
pretty * pretty_sub_pretty(pretty * parent, int from, int to);
void pretty_genome(pretty* pa);
void pretty_ls(pretty* pa);
void pretty_cs(pretty* pa);
void pretty_qualities(pretty* pa);
void pretty_stats(pretty* pa);
//pretty *  pretty_from_sam(samfile_t * sam_file, bam1_t * entry);
pretty * pretty_from_string(char *sam_string);
void pretty_print(pretty* pa);
void pretty_free(pretty* pa);
void pretty_match(pretty* pa);
int pretty_remap(pretty * pa, unsigned long offset_start, unsigned long offset_end);
void pretty_print_sam(FILE * f, pretty * pa);
pretty * pretty_from_string_fast(char* sam_string) ;
pretty * pretty_new();
void pretty_print_sam_update(pretty * pa);
#endif
