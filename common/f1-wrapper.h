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

#undef EXTERN
#undef STATIC
#ifdef _MODULE_GMAPPER
#include "../gmapper/gmapper-defaults.h"
#define EXTERN(_type, _id, _init_val) _type _id = _init_val
#define STATIC(_type, _id, _init_val) static _type _id = _init_val
#else
#define EXTERN(_type, _id, _init_val) extern _type _id
#define STATIC(_type, _id, _init_val)
#endif

static int const f1_window_cache_size = 1048576;

EXTERN(uint64_t, f1_calls_bypassed, 0);

/* Thread-private */
EXTERN(uint32_t, f1_hash_tag, 0);

typedef struct f1_window_cache_entry {
  uint32_t tag;
  uint32_t score;
} f1_window_cache_entry;
EXTERN(struct f1_window_cache_entry *, f1_window_cache, NULL);

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

  if (gapless_sw && shrimp_mode == MODE_LETTER_SPACE)
    return sw_gapless_setup(_match, _mismatch, reset_stats)
      || sw_vector_setup(_dblen, _qrlen,
			   _a_gap_open, _a_gap_ext, _b_gap_open, _b_gap_ext, _match, _mismatch,
			   _use_colours, reset_stats);

  else if (gapless_sw) // and colour space
    return sw_gapless_setup(_match, _mismatch, reset_stats);

  else // gapped alignment
    return sw_vector_setup(_dblen, _qrlen,
			   _a_gap_open, _a_gap_ext, _b_gap_open, _b_gap_ext, _match, _mismatch,
			   _use_colours, reset_stats);
}
static inline int
f1_free(void) 
{
	free(f1_window_cache);
	return 0;
}

/*
 * Get statistics of this filter.
 * Called independently by various threads.
 */
static inline void
f1_stats(uint64_t *invocs, uint64_t *cells, double *secs, uint64_t *calls_bypassed)
{
  if (calls_bypassed != NULL)
    *calls_bypassed = f1_calls_bypassed;

  if (!gapless_sw)
    return sw_vector_stats(invocs, cells, secs);
  else
    return sw_gapless_stats(invocs, cells, NULL);
}


/*
 * Run SW filter. Called independently by different threads.
 * If tag != 0, look up score in hash table.
 */
static inline int
f1_run(uint32_t * genome, int glen, int goff, int wlen, uint32_t * read, int rlen, int g_idx, int r_idx,
       uint32_t * genome_ls, int init_bp, bool is_rna, uint tag, bool gapless)
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
  if (!gapless) {
	//fprintf(stderr,"Doing a vector call, %d, %d, %d\n", goff, wlen, rlen);
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
