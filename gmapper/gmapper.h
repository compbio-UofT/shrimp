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
#define DEF_WINDOW_OVERLAP	80.0
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

#define HASH_TABLE_POWER	12	/* 4^HASH_POWER entries in table */

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
#define MAX_HASH_SEED_WEIGHT	128
#define MAX_HASH_SEED_SPAN	128
#define HASH_TABLE_POWER	12	/* 4^HASH_POWER entries in table */

/*
 * If window_len, sw_vect_threshold, sw_full_threshold are absolute values,
 * we'll set them negative to distinguish.
 */
#define IS_ABSOLUTE(x)	((x) < 0)

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
  uint32_t	span;		/* max 64 (could be uint8_t) */
  uint32_t	weight;		/* max 64 */
  // aligned to 8B
};

struct re_score {
  struct sw_full_results * sfrp;	/* alignment results (final pass) */
  struct anchor	anchor;

  uint		contig_num;	/* contig index (for filename)*/
  union {
    int		score;		/* doubles as heap cnt in [0] */
    uint	heap_elems;
  };
  union {
    uint	g_idx;		/* doubles as heap alloc in [0]*/
    uint	heap_capacity;
  };
  bool		rev_cmpl;	/* from contig's reverse cmpl */
};

struct range_restriction {
  uint	cn;
  uint	st;
  uint	g_start;
  uint	g_end;
};

struct read_entry {
  char *	name;
  char *	seq;
  re_score *	scores;
  uint32_t *	read[2];	/* the read as a bitstring */

  uint32_t *	mapidx[2];	/* per-seed list of mapidxs in read */
  bool *	mapidx_pos[2];	/* per-seed list of validity of mapidx positions in read; only if read has Ns */

  struct uw_anchor *	anchors[2];	/* list of anchors */

  struct read_hit *	hits[2];	/* list of hits */

  struct range_restriction * ranges;
  char *		range_string;

  uint		n_anchors[2];
  uint		n_hits[2];
  uint		n_ranges;

  uint32_t	sw_hits;
  uint32_t	final_matches;

  int8_t	initbp[2];		/* colour space init letter */
  uint8_t	read_len;
  uint8_t	window_len;

  uint8_t	max_n_kmers;	/* = read_len - min_seed_span + 1 */
  uint8_t	min_kmer_pos;	/* = 0 in LS; = 1 in CS */
  uint8_t	input_strand;
  bool		has_Ns;
  bool		is_rna;
};

struct read_hit {
  struct sw_full_results *	sfrp;
  struct anchor	anchor;
  uint		g_off;
  int		score_window_gen;
  int		score_vector;
  int		pct_score_vector;
  int		score_full;
  int		pct_score_full;
  int		score_max;
  uint		matches;
  uint		cn;
  int		pair_min;
  int		pair_max;
  uint16_t	w_len;
  uint8_t	st;
  uint8_t	gen_st;
};

struct read_hit_pair_holder {
  struct read_hit *	hit[2];
  int			insert_size;
};

struct read_hit_holder {
  struct read_hit *	hit;
};


#endif
