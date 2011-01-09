#include "util.h"

/*	$Id: fasta.h,v 1.12 2009/06/16 23:26:21 rumble Exp $	*/

#ifndef _FASTA_H_
#define _FASTA_H_

#define LETTER_SPACE	1
#define COLOUR_SPACE	2

/*
 * We're presently using 4 bits, so we max out at 16 different bases. This just
 * works out if we treat N and X as the same.
 *
 * from: http://genome.ucsc.edu/FAQ/FAQdownloads#download5
 */

#define BASE_A		0		/* Adenine */
#define BASE_C		1		/* Cytosine */
#define BASE_G		2		/* Guanine */
#define BASE_T		3		/* Thymine */
#define BASE_U		4		/* Uracil */
#define BASE_M		5		/* A or C */
#define BASE_R		6		/* A or G (Purine) */
#define BASE_W		7		/* A or T */
#define BASE_S		8		/* C or G */
#define BASE_Y		9		/* C or T (Pyrimidine) */
#define BASE_K		10		/* G or T */
#define BASE_V		11		/* A or C or G (not T) */
#define BASE_H		12		/* A or C or T (not G) */
#define BASE_D		13		/* A or G or T (not C) */
#define BASE_B		14		/* C or G or T (not A) */
#define BASE_X		15		/* G or A or T or C (any base) */
#define BASE_N		15		/* G or A or T or C (any base) */

#define BASE_LS_MIN	BASE_A
#define BASE_LS_MAX	BASE_B

#define BASE_0		0
#define BASE_1		1
#define BASE_2		2
#define BASE_3		3

#define BASE_CS_MIN	BASE_0
#define BASE_CS_MAX	BASE_3

#define FASTA_PER_LINE 80

typedef struct _fasta_t {
	gzFile fp;
	char  *file;
	shrimp_mode_t space;
	char   buffer[8*1024*1024];
	char   translate[256];
	bool   leftover;
	bool	fastq;
	char *	parse_buffer;
	uint32_t parse_buffer_size;
	//for fast_gzgets 
	char *save_buf;
	int   save_len;
	int   save_bytes;
	int   save_skip;
	bool header;
} * fasta_t;

typedef struct _fasta_stats_t {
	uint64_t	total_ticks;
} * fasta_stats_t;

struct read_entry {
  char *        name;
  char *        seq;
  char *	orig_seq;
  char *        qual;
  char *	orig_qual;
  char *        plus_line; //The '+' line in fastq
  uint32_t *    read[2];        /* the read as a bitstring */

  uint32_t *    mapidx[2];      /* per-seed list of mapidxs in read */

  struct anchor *       anchors[2];     /* list of anchors */

  struct read_hit *     hits[2];        /* list of hits */

  struct range_restriction * ranges;
  char *        range_string;

  int           n_anchors[2];
  int           n_hits[2];
  int           n_ranges;

  int           final_matches;
  int           final_dup_matches;

  int           initbp[2];              /* colour space init letter */
  int           read_len;
  int           window_len;

  int           max_n_kmers;    /* = read_len - min_seed_span + 1 */
  int           min_kmer_pos;   /* = 0 in LS; = 1 in CS */
  int           input_strand;
  bool          is_rna;
  bool		ignore;
  bool		paired;
  bool		first_in_pair;
  struct read_entry * mate_pair;
};


fasta_t	  fasta_open(const char *, shrimp_mode_t, bool);
void	  fasta_close(fasta_t);
//bool	  fasta_get_next_with_range(fasta_t, char **, char **, bool *, char **, char **);
bool	  fasta_get_next_read_with_range(fasta_t, read_entry * re);
int	  fasta_get_initial_base(int, char *);
uint32_t *fasta_bitfield_to_colourspace(fasta_t, uint32_t *, uint32_t, bool);
uint32_t *fasta_sequence_to_bitfield(fasta_t, char *);
fasta_stats_t fasta_stats(void);
char      base_translate(int, bool);
void	fasta_write_read(FILE* file, read_entry * re);
void	fasta_write_fasta(FILE* file, char* seq);

static inline bool
fasta_get_next_contig(fasta_t file, char **name, char **seq, bool *is_rna) {
	read_entry re;
  //return fasta_get_next_with_range(file, name, seq, is_rna, NULL, NULL);
  	bool ret = fasta_get_next_read_with_range(file, &re);
	if (ret) {
		*seq=re.seq;
		*name=re.name;
		if (is_rna!=NULL) {
			*is_rna=re.is_rna;
		}
	}
	return ret;	
	
}


#endif
