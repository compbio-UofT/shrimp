#include <assert.h>
#include <ctype.h>
#include <errno.h>
#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <unistd.h>
#include <zlib.h>
#include <sys/time.h>

#include "../common/util.h"
#include "../common/sw-gapless.h"
#include "../common/stats.h"


static int	initialised;
static int	match, mismatch;

/* statistics */
static count_t	ticks, cells, invocs;

#pragma omp threadprivate(initialised, match, mismatch, ticks, cells, invocs)


int
sw_gapless_setup(int _match, int _mismatch, bool reset_stats)
{
  match = _match;
  mismatch = _mismatch;

  if (reset_stats) {
    count_init(&ticks);
    count_init(&cells);
    count_init(&invocs);
  }

  initialised = 1;

  return 0;
}


void
sw_gapless_stats(uint64_t * _invocs, uint64_t * _cells, uint64_t * _ticks)
{
  if (_invocs != NULL)
    *_invocs = (uint64_t)count_get_count(&invocs);
  if (_cells != NULL)
    *_cells = (uint64_t)count_get_count(&cells);
  if (_ticks != NULL)
    *_ticks = (uint64_t)count_get_count(&ticks);
}

int
sw_gapless(uint32_t * genome, int glen, uint32_t * read, int rlen, int g_idx, int r_idx,
	   uint32_t * genome_ls, int init_bp, bool is_rna)
{
  int score;
  int g_left, g_right, r_left, r_right;
  int max_score;

  llint before = rdtsc(), after;

  if (!initialised)
    abort();

  count_increment(&invocs);

  if (g_idx < r_idx) {
    g_left = 0;
    r_left = r_idx - g_idx;
  } else {
    g_left = g_idx - r_idx;
    r_left = 0;
  }
  g_right = g_left;
  r_right = r_left;

  score = 0;
  if (genome_ls != NULL && r_left == 0) { // forcefully match first colour in read
    int real_colour = lstocs(EXTRACT(genome_ls, g_right), init_bp, is_rna);
    if (real_colour == (int)EXTRACT(read, 0)) {
      score = match;
    } else {
      r_left++;
      g_left++;
    }
    r_right++;
    g_right++;
  }

  max_score = score;

  while (g_right < glen && r_right < rlen) {
    score += (EXTRACT(genome, g_right) == EXTRACT(read, r_right)? match : mismatch);

    if (score > max_score)
      max_score = score;

    g_right++;
    r_right++;
    if (score < 0) {
      g_left = g_right;
      r_left = r_right;
      score = 0;
    }
  }

  count_add(&cells, rlen);
  after = rdtsc();
  count_add(&ticks, MAX(after - before, 0));

  return max_score;
}
