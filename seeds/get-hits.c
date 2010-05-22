#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../common/bitmap.h"
#include "../common/util.h"
#include "../common/stats.h"

#define MAX_SEED_SPAN 64

struct seed_type {
  bitmap_type	mask[1];	/* a bitmask, least significant bit = rightmost match */
  uint	span;		/* max 64 (could be uint8_t) */
  uint	weight;		/* max 64 */
  uint step;
  // aligned to 8B
};

struct seed_type * seed;
uint max_seed_span;
uint min_seed_span = MAX_SEED_SPAN;
uint n_seeds;

int target_length, target_mismatches;
int d[10];

int min_total_hits;
int min_d[10];
stat_t total_hits;

bool
add_spaced_seed(char *seedStr)
{
  uint i;
  char * c;

  seed = (struct seed_type *)xrealloc(seed, sizeof(struct seed_type) * (n_seeds + 1));
  seed[n_seeds].mask[0] = 0x0;
  c = strtok(seedStr, ",");

  seed[n_seeds].span = strlen(c);
  seed[n_seeds].weight = strchrcnt(c, '1');

  if (seed[n_seeds].span < 1
      || seed[n_seeds].span > MAX_SEED_SPAN
      || seed[n_seeds].weight < 1
      || strchrcnt(c, '0') != seed[n_seeds].span - seed[n_seeds].weight)
    return false;

  for (i = 0; i < seed[n_seeds].span; i++)
    bitmap_prepend(seed[n_seeds].mask, 1, (c[i] == '1' ? 1 : 0));

  if (seed[n_seeds].span > max_seed_span)
    max_seed_span = seed[n_seeds].span;

  if (seed[n_seeds].span < min_seed_span)
    min_seed_span = seed[n_seeds].span;

  c = strtok(NULL, ",");
  if (c == NULL) {
    seed[n_seeds].step = 1;
  } else {
    seed[n_seeds].step = atoi(c);
    if (seed[n_seeds].step < 1 || seed[n_seeds].step > 16)
      return false;
  }

  n_seeds++;

  return true;
}


int count_hits(int sn, int first_pos)
{
  int i, r = 0;

  for ( ; first_pos + seed[sn].span - 1 < target_length; first_pos += seed[sn].step) {
    for (i = 0; i < target_mismatches; i++) {
      if (d[i] >= first_pos
	  && d[i] <= first_pos + seed[sn].span - 1
	  && bitmap_extract(seed[sn].mask, 1, d[i] - first_pos) == 1)
	break;
    }
    if (i == target_mismatches)
      r++;
  }

  return r;
}


void process()
{
  int i;
  for (i = 0; i < target_mismatches; i++) {
    fprintf(stderr, " %d", d[i]);
  }
  fprintf(stderr, "\n");

  // count number of seed hits
  int sn, sum_hits = 0;
  for (sn = 0; sn < n_seeds; sn++) {
    int first_pos, min_hits = target_length;
    for (first_pos = 0; first_pos < seed[sn].step; first_pos++) {
      int hits = count_hits(sn, first_pos);
      if (hits < min_hits)
	min_hits = hits;
    }

    fprintf(stderr, " (sn:%d,hits:%d)", sn, min_hits);
    sum_hits += min_hits;
  }
  fprintf(stderr, "\n");

  if (sum_hits < min_total_hits) {
    for (i = 0; i < target_mismatches; i++) {
      min_d[i] = d[i];
    }
    min_total_hits = sum_hits;
  }
  stat_add(&total_hits, sum_hits);
}


void rec_gen(int k)
{
  if (k == target_mismatches) {
    process();
  } else {
    int i;
    for (i = (k == 0? 0 : d[k-1] + 1); i < target_length; i++) {
      d[k] = i;
      rec_gen(k+1);
    }
  }
}


int main(int argc, char * argv[]) {
  int n, i;

  if (argc < 5) {
    fprintf(stderr, "use: %s <n_seeds> <seed_1> .. <target_length> <target_mismatches>\n", argv[0]);
    exit(1);
  }

  n = atoi(argv[1]);
  if (n <= 0 || n > 127) {
    fprintf(stderr, "error: invalid number of seeds [%s]\n", argv[1]);
    exit(1);
  }

  for (i = 0; i < n; i++) {
    if (!add_spaced_seed(argv[2+i])) {
      fprintf(stderr, "error: invalid seed [%s]\n", argv[2+i]);
      exit(1);
    }
  }

  target_length = atoi(argv[2+n]);
  if (target_length <= 0 || target_length > 127) {
    fprintf(stderr, "error: invalid target_length [%s]\n", argv[2+n]);
    exit(1);
  }
  min_total_hits = target_length;

  target_mismatches = atoi(argv[2+n+1]);
  if(target_mismatches < 0 || target_mismatches > 9) {
    fprintf(stderr, "error: invalid target_mismatches [%s]\n", argv[2+n+1]);
    exit(1);
  }

  rec_gen(0);
  fprintf(stdout, "avg_hits:%.2f std_dev:%.2f min_hits:%d\n",
	  stat_get_mean(&total_hits), stat_get_sample_stddev(&total_hits), min_total_hits);
  for (i = 0; i < target_mismatches; i++) {
    fprintf(stdout, " %d", min_d[i]);
  }
  fprintf(stdout, "\n");

  return 0;
}
