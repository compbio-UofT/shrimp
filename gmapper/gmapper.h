#ifndef _GMAPPER_H
#define _GMAPPER_H

#include <stdlib.h>
#include "../common/bitmap.h"
#include "../common/sw-full-common.h"
#include "../common/debug.h"
#include "../common/anchors.h"
#include "../common/heap.h"

/* pair mode definitions */
#define PAIR_NONE	0
#define PAIR_OPP_IN	1
#define PAIR_OPP_OUT	2
#define PAIR_COL_FW	3
#define PAIR_COL_BW	4

static char const * const pair_mode_string[5] =
  { "none", "opposing strands; inwards", "opposing strands; outwards",
    "same strand; second is forward", "same strand; second is backward" };

static bool const pair_reverse[5][2] =
  { { 0, 0 }, // PAIR_NONE
    { 0, 0 }, // PAIR_OPP_IN
    { 1, 1 }, // PAIR_OPP_OUT
    { 0, 1 }, // PAIR_COL_FW
    { 1, 0 }  // PAIR_COL_BW
  };

/* defaults */
#define DEF_NUM_THREADS 1
#define DEF_CHUNK_SIZE 1000

#define DEF_HASH_FILTER_CALLS	true
#define DEF_GAPLESS_SW		false
#define DEF_LIST_CUTOFF		4294967295u // 2^32 - 1

#define DEF_PAIR_MODE		PAIR_NONE
#define DEF_MIN_INSERT_SIZE	50
#define DEF_MAX_INSERT_SIZE	2000

#define DEF_WINDOW_LEN		140.0
#define DEF_WINDOW_OVERLAP	20.0
#define DEF_NUM_MATCHES		2
#define DEF_NUM_OUTPUTS		10
#define DEF_ANCHOR_WIDTH	8	/* width around anchors in full SW */

/* SW Scores */
#define DEF_MATCH_VALUE		10
#define DEF_MISMATCH_VALUE	-15
#define DEF_A_GAP_OPEN		-40
#define DEF_B_GAP_OPEN		 DEF_A_GAP_OPEN
#define DEF_A_GAP_EXTEND	-7
#define DEF_B_GAP_EXTEND	 DEF_A_GAP_EXTEND
#define DEF_XOVER_PENALTY	-14	/* CS only */

/* Score Thresholds */
#define DEF_WINDOW_GEN_THRESHOLD	55.0	/* Min required to generate match window */
#define DEF_SW_VECT_THRESHOLD	60.0	/* == DEF_SW_FULL_THRESHOLD in lspace */
#define DEF_SW_FULL_THRESHOLD	68.0	/* read_length x match_value x .68 */

/*
 * The maximum seed weight (maximum number of 1's in the seed) sets an
 * upper limit on our lookup table allocation size. The memory usage of
 * rmapper corresponds strongly to 4^MAX_SEED_WEIGHT * (sizeof(void *) +
 * sizeof(uint32_t)). At 16, this is 32GB on 32-bit and 48GB on 64-bit
 * architectures.
 */
#define MAX_SEED_WEIGHT		14
#define MAX_SEED_SPAN		64

/*
 * For larger seeds we'll just use a hash table. Presently, we're restricted to
 * 128 bytes in kmer_to_mapidx, but it's trivially extended.
 */
#define MAX_HASH_SEED_WEIGHT	64
#define MAX_HASH_SEED_SPAN	64
#define HASH_TABLE_POWER	12	/* 4^HASH_POWER entries in table */


#define MAX_ANCHOR_LIST		10000000


/*
 * If window_len, sw_vect_threshold, sw_full_threshold are absolute values,
 * we'll set them negative to distinguish.
 */
#define IS_ABSOLUTE(x)	((x) < 0)

static inline double abs_or_pct(double x, double base) {
  return IS_ABSOLUTE(x) ? -x : base * (x / 100.0);
}

static int const default_spaced_seeds_cs_cnt = 4;
static char const * const default_spaced_seeds_cs[] =
  //{ "111110011111", "111100110001111", "111100100100100111", "111001000100001001111" };
  //{ "1111001111111", "1111100110001111", "11110010010001001111", "11100110010000100100111" };
{ "111110001111111", "111100111001001111", "111001001000111001111", "1111001000010001001001111" };

static int const default_spaced_seeds_ls_cnt = 4;
static char const * const default_spaced_seeds_ls[] =
  //{ "111110011111", "111100110001111", "111100100100100111", "111001000100001001111" };
  //{ "1111001111111", "1111100110001111", "11110010010001001111", "11100110010000100100111" };
{ "111101101011111", "111010110011001111", "1110110001011010111", "11110010100000100110111" };

struct seed_type {
  bitmap_type	mask[1];	/* a bitmask, least significant bit = rightmost match */
  int		span;		/* max 64 (could be uint8_t) */
  int		weight;		/* max 64 */
};

struct range_restriction {
  int		cn;
  int		st;
  llint		g_start;
  llint		g_end;
};

struct read_entry {
  char *	name;
  char *	seq;
  uint32_t *	read[2];	/* the read as a bitstring */

  uint32_t *	mapidx[2];	/* per-seed list of mapidxs in read */

  struct anchor *	anchors[2];	/* list of anchors */

  struct read_hit *	hits[2];	/* list of hits */

  struct range_restriction * ranges;
  char *	range_string;

  int		n_anchors[2];
  int		n_hits[2];
  int		n_ranges;

  int		final_matches;
  int		final_dup_matches;

  int		initbp[2];		/* colour space init letter */
  int		read_len;
  int		window_len;

  int		max_n_kmers;	/* = read_len - min_seed_span + 1 */
  int		min_kmer_pos;	/* = 0 in LS; = 1 in CS */
  int		input_strand;
  bool		is_rna;
};

struct read_hit {
  struct sw_full_results *	sfrp;
  struct anchor	anchor;
  llint		g_off;
  int		score_window_gen;
  int		score_vector;
  int		pct_score_vector;
  int		score_full;
  int		pct_score_full;
  int		score_max;
  int		matches;
  int		cn;
  int		pair_min; // -1 means none
  int		pair_max;
  int		w_len;
  int		st;
  int		gen_st;
};

struct read_hit_pair_holder {
  struct read_hit *	hit[2];
  int			insert_size;
};

struct read_hit_holder {
  struct read_hit *	hit;
};


#endif
