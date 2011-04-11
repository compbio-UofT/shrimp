/*
 * This file should contain definitions only. This file should be safe to include by itself.
 *
 * DO include here:
 *   - #defines (excluding default settings)
 *   - type definitions
 *
 * DO NOT include here:
 *   - #defines with default settings
 *   - constants which require memory allocation
 */

#ifndef _GMAPPER_DEFINITIONS_H
#define _GMAPPER_DEFINITIONS_H

#include <stdlib.h>
#include "../common/bitmap.h"
#include "../common/heap.h"
#include "../common/stats.h"

typedef long long int llint;


typedef struct {
  char ** argv;
  int argc;
} shrimp_args_t;

typedef enum {
  MODE_LETTER_SPACE = 1,
  MODE_COLOUR_SPACE = 2,
  MODE_HELICOS_SPACE = 3
} shrimp_mode_t;


/* pair mode */
#define PAIR_NONE	0
#define PAIR_OPP_IN	1
#define PAIR_OPP_OUT	2
#define PAIR_COL_FW	3
#define PAIR_COL_BW	4


/* seeds */
#define MAX_SEED_WEIGHT		14
#define MAX_SEED_SPAN		64
#define MAX_HASH_SEED_WEIGHT	64
#define MAX_HASH_SEED_SPAN	64
#define HASH_TABLE_POWER	12	/* 4^HASH_POWER entries in table */

#define MAX_N_DEFAULT_SEEDS 5
typedef char const * const default_seed_array_t[MAX_N_DEFAULT_SEEDS];

typedef struct seed_type {
  bitmap_type	mask[1];	/* a bitmask, least significant bit = rightmost match */
  int		span;		/* max 64 (could be uint8_t) */
  int		weight;		/* max 64 */
} seed_type;


/* mapping structs */
struct anchor {
  llint		x;
  llint		y;
  int		length;
  int		width;
  int		weight;
  int		cn;
};

typedef struct read_entry {
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

  int		delta_g_off_min[2];
  int		delta_g_off_max[2];
  int		delta_region_min[2];
  int		delta_region_max[2];	// where to find the mp; only used for first read in a pair

  int           max_n_kmers;    /* = read_len - min_seed_span + 1 */
  int           min_kmer_pos;   /* = 0 in LS; = 1 in CS */
  int           input_strand;
  bool          is_rna;
  bool		ignore;
  bool		paired;
  bool		first_in_pair;
  struct read_entry * mate_pair;
} read_entry;


typedef struct read_hit {
  struct sw_full_results *      sfrp;
  struct anchor anchor;
  llint         g_off;
  int           score_window_gen;
  int           score_vector;
  int           pct_score_vector;
  int           score_full;
  double	pct_score_full;
  int		pass1_key;
  int		pass2_key;
  int           score_max;
  int           matches;
  int           cn;
  int           pair_min; // -1 means none
  int           pair_max;
  int           w_len;
  int           st;
  int           gen_st;
} read_hit;

typedef struct read_hit_pair {
  struct read_hit *	rh[2];
  int			score_max;
  int			score;
  int			pct_score;
  int			key;
  int			insert_size;
} read_hit_pair;

struct read_hit_pair_holder {
  struct read_hit *     hit[2];
  int                   insert_size;
};
 
struct read_hit_holder {
  struct read_hit *     hit;
};


/* other */
struct range_restriction {
  int		cn;
  int		st;
  llint		g_start;
  llint		g_end;
};


/* cigar string */
typedef struct {
	uint32_t * lengths;
	char * ops;
	int size;
} cigar_t;


typedef struct regions_options {
  bool		recompute;
  int		min_seed;
  int		max_seed;
} regions_options;

typedef struct anchor_list_options {
  bool		recompute;			// whether to recompute anchor list for each read
  bool		collapse;
  bool		use_region_counts;		// whether to use region counts for each read
  int		min_count[2];
  int		max_count[2];			// min/max[0]: min/max count for this read; min/max[1]: for mp
  //int		min_seed;
  //int		max_seed;			// which seeds to use in creating the anchor list
} anchor_list_options;

typedef struct hit_list_options {
  bool		recompute;
  bool		gapless;
  int		match_mode;
  double	threshold;
} hit_list_options;

typedef struct pass1_options {
  bool		recompute;
  bool		gapless;
  bool		only_paired;
  int		num_outputs;
  double	threshold;
  double	window_overlap;
} pass1_options;

typedef struct pass2_options {
  bool		strata;
  int		num_outputs;
  int		stop_count;
  double	threshold;
  double	stop_threshold;
} pass2_options;


typedef struct read_mapping_options_t {
  // region handling
  struct regions_options regions;

  // anchor list
  struct anchor_list_options anchor_list;

  // hit list
  struct hit_list_options hit_list;

  // vector SW
  struct pass1_options pass1;

  // scalar/full SW
  struct pass2_options pass2;

} read_mapping_options_t;

typedef struct pairing_options {
  int		pair_mode;
  int		min_insert_size;
  int		max_insert_size;		// for read 1 relative to read 0
  int		min_num_matches;
  int		pass1_num_outputs;
  int		pass2_num_outputs;
  int		stop_count;

  double	pass1_threshold;
  double	pass2_threshold;		// thresholds for the pair
  double	stop_threshold;

  bool		strata;
  bool		pair_up_hits;
} pairing_options;

typedef struct readpair_mapping_options_t {
  // initial computation of region counts controlled by global flag

  // handle readpair or each read
  struct pairing_options pairing;

  struct read_mapping_options_t read[2];

} readpair_mapping_options_t;


static inline void
read_free_anchor_list(struct read_entry * re)
{
  if (re->anchors[0] != NULL) {
    free(re->anchors[0]);
    re->anchors[0] = NULL;
    re->n_anchors[0] = 0;
  }
  if (re->anchors[1] != NULL) {
    free(re->anchors[1]);
    re->anchors[1] = NULL;
    re->n_anchors[1] = 0;
  }
}


static inline void
read_free_hit_list(struct read_entry * re)
{
  if (re->hits[0] != NULL) {
    free(re->hits[0]);
    re->hits[0] = NULL;
    re->n_hits[0] = 0;
  }
  if (re->hits[1] != NULL) {
    free(re->hits[1]);
    re->hits[1] = NULL;
    re->n_hits[1] = 0;
  }
}


#endif
