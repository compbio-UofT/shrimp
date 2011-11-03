#ifndef __GEN_ST_H
#define __GEN_ST_H

#include <stdint.h>
#include <assert.h>

//#define GEN_ST_BASE t->b
#define GEN_ST_BASE 17

typedef struct {
  uint32_t *	a;
  int *	pow;
  int	b;
  int	h;
  int	n_keys;
  int	n_nodes;
} gen_st;


void gen_st_init(gen_st *, int, uint32_t *, int);
void gen_st_delete(gen_st *);


static inline int
gen_st_search_node(uint32_t * node, int load, uint32_t val)
{
  assert(node != NULL);

  
  int i;

  for (i = 0; i < load; i++)
    if (node[i] > val)
      return i;

  return i;
  
  /*
  int l, r, m;
  l = 0;
  r = load;
  while (l < r) {
    m = (r + l)/2;
    if (node[m] <= val)
      l = m + 1;
    else
      r = m;
  }
  return l;
  */
}


static inline int
gen_st_search(gen_st * t, uint32_t val)
{
  int node_depth, node_lev_idx, node_abs_idx, nodes_above;
  uint32_t * node;
  int range_start, range_end;
  int k, h, idx, delta;

  assert(t != NULL);

  node_depth = 0;
  node_lev_idx = 0;
  node_abs_idx = 0;
  nodes_above = 0;
  node = &t->a[0];
  range_start = 0;
  range_end = t->n_keys;
  h = t->h;

  while (h > 1) {
    assert(range_end - range_start > GEN_ST_BASE - 1);
    assert(t->pow[h - 1] - 1 < range_end - range_start && range_end - range_start <= t->pow[h] - 1);

    idx = gen_st_search_node(node, GEN_ST_BASE - 1, val);
    idx--;

    k = ((range_end - range_start) - (t->pow[h - 1] - 1)) / ((GEN_ST_BASE - 1) * t->pow[h - 2]);
    assert(0 <= k && k <= GEN_ST_BASE);

    if (k == GEN_ST_BASE) {
      range_start += (idx + 1) * t->pow[h - 1];
      range_end = range_start + t->pow[h - 1] - 1;
    } else {
      delta = (range_end - range_start) - (GEN_ST_BASE - 1) - k * (t->pow[h - 1] - 1) - (GEN_ST_BASE - k - 1) * (t->pow[h - 2] - 1);
      if (idx < k) {
	range_start += (idx + 1) * t->pow[h - 1];
	if (idx + 1 < k) {
	  range_end = range_start + t->pow[h - 1] - 1;
	} else {
	  assert(idx + 1 == k);
	  range_end = range_start + delta;
	}
      } else { // idx >= k
	range_start += k * t->pow[h - 1] + (delta + 1) + (idx - k) * t->pow[h - 2];
	range_end = range_start + t->pow[h - 2] - 1;
      }
    }
    node_lev_idx = node_lev_idx * GEN_ST_BASE + idx + 1;
    nodes_above += t->pow[node_depth];
    node_depth++;
    node_abs_idx = nodes_above + node_lev_idx;
    node = &t->a[node_abs_idx * (GEN_ST_BASE - 1)];

    assert(0 <= node_lev_idx && node_lev_idx < t->pow[node_depth]);
    assert(node_depth < t->h);

    if (range_end - range_start == t->pow[h - 2] - 1)
      h -= 2;
    else
      h--;

    assert(h == 0 || node_abs_idx < t->n_nodes);
  }

  assert(range_end - range_start <= GEN_ST_BASE - 1);

  if (h == 1) {
    idx = gen_st_search_node(node, GEN_ST_BASE - 1, val);
    range_start += idx;
  }

  return range_start - 1;
}


#endif
