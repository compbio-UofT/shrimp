#ifndef _F1_WRAPPER_H
#define _F1_WRAPPER_H

#include <stdlib.h>
#include <stdint.h>
#include <inttypes.h>
#include <sys/types.h>
#include <stdbool.h>
#include <assert.h>
#include <string.h>

#include "../common/util.h"
#include "../common/sw-vector.h"
#include "../common/sw-gapless.h"

static int const f1_window_cache_size = 1048576;

static uint64_t f1_calls_bypassed;


/* Thread-private */
static uint32_t f1_hash_tag;
struct f1_window_cache_entry {
  uint32_t tag;
  uint32_t score;
};
static struct f1_window_cache_entry * f1_window_cache;

#pragma omp threadprivate(f1_hash_tag, f1_window_cache)


/*
 * Set up SW filter.
 * Must be called by each thread prior to using the filter.
 */
static inline int
f1_setup(int _dblen, int _qrlen,
	 int _a_gap_open, int _a_gap_ext, int _b_gap_open, int _b_gap_ext, int _match, int _mismatch,
	 int _use_colours, bool reset_stats)
{
  f1_hash_tag = 0;
  f1_window_cache = (struct f1_window_cache_entry *)xcalloc(f1_window_cache_size * sizeof(f1_window_cache[0]));

  if (!gapless_sw)
    return sw_vector_setup(_dblen, _qrlen,
			   _a_gap_open, _a_gap_ext, _b_gap_open, _b_gap_ext, _match, _mismatch,
			   _use_colours, reset_stats);
  else
    return sw_gapless_setup(_match, _mismatch, reset_stats);
}


/*
 * Get statistics of this filter.
 * Called independently by various threads.
 */
static inline void
f1_stats(uint64_t *invocs, uint64_t *cells, uint64_t *ticks, uint64_t *calls_bypassed)
{
  if (calls_bypassed != NULL)
    *calls_bypassed = f1_calls_bypassed;

  if (!gapless_sw)
    return sw_vector_stats(invocs, cells, ticks);
  else
    return sw_gapless_stats(invocs, cells, ticks);
}


/*
 * Run SW filter. Called independently by different threads.
 * If tag != 0, look up score in hash table.
 */
static inline int
f1_run(uint32_t * genome, int glen, int goff, int wlen, uint32_t * read, int rlen, int g_idx, int r_idx,
       uint32_t * genome_ls, int init_bp, bool is_rna, uint tag)
{
  uint32_t hash_val = 0;
  int score;
  
  /* Look-up */
  if (hash_filter_calls && tag != 0) {
    hash_val = hash_genome_window(genome, goff, wlen) % f1_window_cache_size;

    if (f1_window_cache[hash_val].tag == tag) { // Cache hit
#pragma omp atomic
      f1_calls_bypassed++;

      return f1_window_cache[hash_val].score;
    }
  }

  /* Compute */
  if (!gapless_sw) {
    score = sw_vector(genome, goff, wlen, read, rlen,
		      genome_ls, init_bp, is_rna);
  } else {
    score = sw_gapless(genome, glen, read, rlen, g_idx, r_idx,
		       genome_ls, init_bp, is_rna);
  }

  /* Save */
  if (hash_filter_calls && tag != 0) {
    f1_window_cache[hash_val].tag = tag;
    f1_window_cache[hash_val].score = score;
  }

  return score;
}


#endif
