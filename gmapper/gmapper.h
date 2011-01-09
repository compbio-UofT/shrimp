#ifndef _GMAPPER_H
#define _GMAPPER_H

#include <getopt.h>
#include <stdlib.h>
#include "../common/bitmap.h"
#include "../common/sw-full-common.h"
#include "../common/debug.h"
#include "../common/anchors.h"
#include "../common/heap.h"

#define DEF_MAX_ALIGNMENTS 0

/* pair mode definitions */
#define PAIR_NONE	0
#define PAIR_OPP_IN	1
#define PAIR_OPP_OUT	2
#define PAIR_COL_FW	3
#define PAIR_COL_BW	4

#define DEF_LONGEST_READ_LENGTH	1000

struct option standard_options[] =
{
	{"un",1,0,10},
	{"al",1,0,11},
	{"upstream",1,0,'1'},
	{"downstream",1,0,'2'},
	{"sam-unaligned",0,0,12},
	{"longest-read",1,0,13},	
	{"seeds",1,0,'s'}, 
	{"report",1,0,'o'},
	{"match-window",1,0,'w'},
	{"cmw-mode",1,0,'n'},
	{"cmw-overlap",1,0,'l'},
	{"anchor-width",1,0,'a'},
	{"save",1,0,'S'},
	{"load",1,0,'L'},
	{"cutoff",1,0,'z'},
	{"match",1,0,'m'},
	{"mismatch",1,0,'i'},
	{"open-r",1,0,'g'},
	{"open-q",1,0,'q'},	
	{"ext-r",1,0,'e'},
	{"ext-q",1,0,'f'},
	{"cmv-threshold",1,0,'r'},
	{"full-threshold",1,0,'h'},
	{"threads",1,0,'N'},
	{"thread-chunk",1,0,'K'},
	{"pair-mode",1,0,'p'},
	{"isize",1,0,'I'},
	{"ungapped",0,0,'U'},
	{"negative",0,0,'C'},
	{"positive",0,0,'F'},
	{"pretty",0,0,'P'},
	{"sam",0,0,'E'},
	{"fastq",0,0,'Q'},
	{"print-reads",0,0,'R'},
	{"rev-tiebreak",0,0,'T'},
	{"tiebreak-off",0,0,'t'},
	{"isize-histogram",0,0,'X'},
	{"proj-histogram",0,0,'Y'},
	{"cachebypass-off",0,0,'Z'},
	{"help",0,0,'?'},
	{"spaced-kmers",0,0,'H'},
	{"thread-stats",0,0,'D'},
	{"trim-off",0,0,'V'},
	{"strata",0,0,9},
	{"max-alignments",1,0,14},
	{"global",0,0,15},
	{"read-group",1,0,17},
	{"sam-header",1,0,18},
	{"half-paired",0,0,19},
	{"sam-r2",0,0,20},
	{"mode",1,0,'M'},
	{"trim-front",1,0,200},
	{"trim-end",1,0,201},
	{"trim-first",1,0,202},
	{"trim-second",1,0,203}
};
struct option colour_space_options[] = {
	{"crossover",1,0,'x'},
	{"vec-threshold",1,0,'v'},
	{"bfast",0,0,16},
	{0,0,0,0}
};
struct option letter_space_options[] = {
		{0,0,0,0}
	}; 
size_t standard_entries = sizeof(standard_options)/sizeof(struct option);
size_t letter_entries = sizeof(letter_space_options)/sizeof(struct option);
size_t colour_entries = sizeof(colour_space_options)/sizeof(struct option);
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
//SHRiMP v 2.0.1
//#define DEF_SW_VECT_THRESHOLD	60.0	/* == DEF_SW_FULL_THRESHOLD in lspace */
//#define DEF_SW_FULL_THRESHOLD	68.0	/* read_length x match_value x .68 */
//SHRiMP v 2.0.2
#define DEF_SW_VECT_THRESHOLD	50.0	/* == DEF_SW_FULL_THRESHOLD in lspace */
#define DEF_SW_FULL_THRESHOLD	55.0	/* read_length x match_value x .55 */

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

#define DEF_MAX_N_DEFAULT_SEEDS 5
typedef char const * const default_seed_array_t[DEF_MAX_N_DEFAULT_SEEDS];

static int const default_min_spaced_seed_weight_cs = 10;
static int const default_max_spaced_seed_weight_cs = 18;
static int const default_spaced_seed_weight_cs = 12;
static int const default_spaced_seeds_cs_cnt[9] = { 4, 4, 4, 0, 0, 0, 4, 0, 4 };
static default_seed_array_t const default_spaced_seeds_cs[9] = {
  { "111110011111", "111100110001111", "111100100100100111", "111001000100001001111" },
  { "1111001111111", "1111100110001111", "11110010010001001111", "11100110010000100100111" },
  { "111110001111111", "111100111001001111", "111001001000111001111", "1111001000010001001001111" },
  { },
  { },
  { },
  { "111111101110111111", "1111100101101101011111", "11110011001010100011011111", "111101001100000100110011010111" },
  { },
  { "11111011111110111111", "11110111011010111011111", "11111100110101101001011111", "11111010101100100010011101111" } };

static int const default_min_spaced_seed_weight_ls = 10;
static int const default_max_spaced_seed_weight_ls = 18;
static int const default_spaced_seed_weight_ls = 12;
static int const default_spaced_seeds_ls_cnt[9] = { 4, 4, 4, 0, 0, 0, 4, 0, 4 };
static default_seed_array_t const default_spaced_seeds_ls[9] = {
  { "111110011111", "111100110001111", "111100100100100111", "111001000100001001111" },
  { "1111001111111", "1111100110001111", "11110010010001001111", "11100110010000100100111" },
  { "111101101011111", "111010110011001111", "1110110001011010111", "11110010100000100110111" },
  { },
  { },
  { },
  { "111111101110111111", "1111100101101101011111", "11110011001010100011011111", "111101001100000100110011010111" },
  { },
  { "11111011111110111111", "11110111011010111011111", "11111100110101101001011111", "11111010101100100010011101111" } };

static int const default_spaced_seeds_mirna_cnt = 5;
static default_seed_array_t const default_spaced_seeds_mirna =
  { "00111111001111111100", "00111111110011111100", "00111111111100111100", "00111111111111001100", "00111111111111110000" } ;

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

typedef struct {
	uint32_t * lengths;
	char * ops;
	int size;
} cigar_t;

struct read_hit_pair_holder {
  struct read_hit *	hit[2];
  int			insert_size;
};

struct read_hit_holder {
  struct read_hit *	hit;
};


#endif
