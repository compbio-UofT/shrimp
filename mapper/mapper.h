/*
 * mapper.h
 *
 *  Created on: 2009-10-02
 *      Author: dlister
 */

#ifndef MAPPER_H_
#define MAPPER_H_

#include "../common/bitmap.h"
#include <stdlib.h>
#include "../common/debug.h"


/* Default parameters */
#define DEF_WINDOW_LEN		135.0
#define DEF_WINDOW_OVERLAP	80.0
#define DEF_NUM_MATCHES		4
#define DEF_NUM_OUTPUTS		100

/* SW Scores */
#define DEF_MATCH_VALUE		10
#define DEF_MISMATCH_VALUE	-15
#define DEF_A_GAP_OPEN		-40
#define DEF_B_GAP_OPEN		 DEF_A_GAP_OPEN
#define DEF_A_GAP_EXTEND	-7
#define DEF_B_GAP_EXTEND	 DEF_A_GAP_EXTEND
#define DEF_XOVER_PENALTY	-14	/* CS only */

#define DEF_SW_VECT_THRESHOLD	60.0	/* == DEF_SW_FULL_THRESHOLD in lspace */
#define DEF_SW_FULL_THRESHOLD	68.0	/* read_length x match_value x .68 */

#define DEF_ANCHOR_WIDTH	8	/* width around anchors in full SW */


#define HASH_TABLE_POWER	12	/* 4^HASH_POWER entries in table */

/*
 * The maximum seed weight (maximum number of 1's in the seed) sets an
 * upper limit on our lookup table allocation size. The memory usage of
 * rmapper corresponds strongly to 4^MAX_SEED_WEIGHT * (sizeof(void *) +
 * sizeof(uint32_t)). At 16, this is 32GB on 32-bit and 48GB on 64-bit
 * architectures.
 */
#ifndef MAX_SEED_WEIGHT
#define MAX_SEED_WEIGHT		14
#endif

/*
 * We hold seeds as bitmaps to reduce cache footprint.
 */
#define MAX_SEED_SPAN		64

/*
 * If window_len, sw_vect_threshold, sw_full_threshold are absolute values,
 * we'll set them negative to distinguish.
 */
#define IS_ABSOLUTE(x)	((x) < 0)

/*
 * For larger seeds we'll just use a hash table. Presently, we're restricted to
 * 128 bytes in kmer_to_mapidx, but it's trivially extended.
 */
#define MAX_HASH_SEED_WEIGHT	128
#define MAX_HASH_SEED_SPAN	128
#define HASH_TABLE_POWER	12	/* 4^HASH_POWER entries in table */

#ifndef DEBUG_SEEDS
static int const default_spaced_seeds_cs_cnt = 4;
static char const * const default_spaced_seeds_cs[] =
  //{ "111110011111", "111100110001111", "111100100100100111", "111001000100001001111" };
  //{ "1111001111111", "1111100110001111", "11110010010001001111", "11100110010000100100111" };
{ "111110001111111", "111100111001001111", "111001001000111001111", "1111001000010001001001111" };
#else
static int const default_spaced_seeds_cs_cnt = 1;
static char const * const default_spaced_seeds_cs[] =
  { "101"};
#endif

#ifndef DEBUG_SEEDS
static int const default_spaced_seeds_ls_cnt = 4;
static char const * const default_spaced_seeds_ls[] =
  //{ "111110011111", "111100110001111", "111100100100100111", "111001000100001001111" };
  //{ "1111001111111", "1111100110001111", "11110010010001001111", "11100110010000100100111" };
{ "111101101011111", "111010110011001111", "1110110001011010111", "11110010100000100110111" };
#else
static int const default_spaced_seeds_ls_cnt = 1;
static char const * const default_spaced_seeds_ls[] =
  {"101"};
#endif

static int const default_spaced_seeds_hs_cnt = 4;
static char const * const default_spaced_seeds_hs[] =
  //{ "111110011111", "111100110001111", "111100100100100111", "111001000100001001111" };
  { "1111001111111", "1111100110001111", "11110010010001001111", "11100110010000100100111" };


extern struct seed_type *seed;
extern uint32_t **seed_hash_mask;
extern uint max_seed_span;
extern uint avg_seed_span;
extern uint n_seeds;
extern u_int	nkmers;				/* total kmers of reads loaded*/


struct seed_type {
  bitmap_type	mask[1];	/* a bitmask, least significant bit = rightmost match */
  uint32_t	span;		/* max 64 (could be uint8_t) */
  uint32_t	weight;		/* max 64 */
  // aligned to 8B
};

extern size_t
power(size_t base, size_t exp);

extern uint32_t
kmer_to_mapidx_hash(uint32_t *kmerWindow, u_int sn);

extern uint32_t
kmer_to_mapidx_orig(uint32_t *kmerWindow, u_int sn);

extern bool
add_spaced_seed(const char *seedStr);

extern void
load_default_seeds();

extern void
init_seed_hash_mask(void);

char *
seed_to_string(uint sn);

#endif /* MAPPER_H_ */
