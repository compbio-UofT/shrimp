#include <stdlib.h>
#include <math.h>
#include <limits.h>
#include <string.h>
#include <assert.h>

#include "gen-st.h"


/*
 * Fill in the subtree with the root at the given depth and level index
 * with the keys from the given a. The subtree will have the given h.
 */
static void
gen_st_fill(gen_st * t, int d, int lev_idx, int h, uint32_t * a, int n)
{
  int abs_idx, delta, k, i, j, prev_j, left_child_lev_idx;
  uint32_t * node;

  if (n == 0) return;

  abs_idx = (t->pow[d] - 1) / (t->b - 1) + lev_idx;
  assert(0 <= abs_idx && abs_idx < t->n_nodes);

  node = &t->a[abs_idx * (t->b - 1)];

  assert(t->pow[h - 1] - 1 < n && n <= t->pow[h] - 1);

  if (h == 1)
    { // no children
      assert(n <= t->b - 1);
      for (i = 0; i < n; i++)
	node[i] = a[i];
    }
  else
    { // has children
      left_child_lev_idx = t->b * lev_idx;
      if (n == t->pow[h] - 1)
	{ // all children are full
	  for (i = 0, prev_j = -1; i < t->b - 1; i++, prev_j = j) {
	    j = (i + 1) * t->pow[h - 1] - 1;
	    assert(j - prev_j - 1 == t->pow[h - 1] - 1);
	    node[i] = a[j];
	    gen_st_fill(t, d + 1, left_child_lev_idx + i, h - 1, &a[prev_j + 1], j - prev_j - 1);
	  }
	  assert(n - prev_j - 1 == t->pow[h - 1] - 1);
	  gen_st_fill(t, d + 1, left_child_lev_idx + i, h - 1, &a[prev_j + 1], n - prev_j - 1);
	}
      else
	{ // some children are not full
	  k = (n - (t->pow[h - 1] - 1)) / ((t->b - 1) * t->pow[h - 2]);
	  delta = n - (t->b - 1) - k * (t->pow[h - 1] - 1) - (t->b - k - 1) * (t->pow[h - 2] - 1);
	  for (i = 0, prev_j = -1; i < t->b - 1; i++, prev_j = j) {
	    if (i < k) {
	      j = (i + 1) * t->pow[h - 1] - 1;
	    } else {
	      j = k * t->pow[h - 1] + (delta + 1) + (i - k) * t->pow[h - 2] - 1;
	    }
	    assert((i < k && j - prev_j - 1 == t->pow[h - 1] - 1)
		   || (i == k && j - prev_j - 1 == delta)
		   || (i > k && j - prev_j - 1 == t->pow[h - 2] - 1));
	    node[i] = a[j];
	    gen_st_fill(t, d + 1, left_child_lev_idx + i,
			(i > k || (i == k && delta == t->pow[h - 2] - 1)? h - 2 : h - 1),
			&a[prev_j + 1], j - prev_j - 1);
	  }
	  assert((i < k && n - prev_j - 1 == t->pow[h - 1] - 1)
		 || (i == k && n - prev_j - 1 == delta)
		 || (i > k && n - prev_j - 1 == t->pow[h - 2] - 1));
	  gen_st_fill(t, d + 1, left_child_lev_idx + i,
		      (i > k || (i == k && delta == t->pow[h - 2] - 1)? h - 2 : h - 1),
		      &a[prev_j + 1], n - prev_j - 1);
	}
    }
}


void
gen_st_init(gen_st * t, int b, uint32_t * a, int n)
{
  int tmp;
  uint32_t * a_aux;

  assert(t != NULL);
  assert(b >= 2);

  t->b = b;
  t->b = GEN_ST_BASE; // hard-coded to be equal to 17 during searching
  t->n_keys = (n == 0? 0 : ((n - 1) / (t->b - 1) + 1) * (t->b - 1));
  a_aux = (uint32_t *)malloc(t->n_keys * sizeof(uint32_t));
  memcpy(a_aux, a, n * sizeof(uint32_t));
  for (int i = n; i < t->n_keys; i++)
    a_aux[i] = UINT_MAX;
  t->n_nodes = t->n_keys / (t->b - 1);

  for (t->h = 0, tmp = 1; n > tmp - 1; t->h++, tmp *= t->b);
  //t->h = (int)ceil(log(t->n_keys + 1) / log(t->b));

  t->pow = (int *)malloc((t->h + 2) * sizeof(int));
  t->pow[0] = 1;
  for (int i = 1; i <= t->h + 1; i++)
    t->pow[i] = t->pow[i - 1] * t->b;

  // finally, set up a
  t->a = (uint32_t *)malloc(t->n_keys * sizeof(uint32_t));
  gen_st_fill(t, 0, 0, t->h, a_aux, t->n_keys);

  free(a_aux);
}


void
gen_st_delete(gen_st * t)
{
  assert(t != NULL);

  free(t->a);
  free(t->pow);
}
