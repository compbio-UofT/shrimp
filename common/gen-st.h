#ifndef __GEN_ST_H
#define __GEN_ST_H

#include <assert.h>


typedef struct {
  int *	a;
  int *	pow;
  int	b;
  int	h;
  int	n_keys;
  int	n_nodes;
} gen_st;


void gen_st_init(gen_st *, int, int *, int);
void gen_st_delete(gen_st *);


static inline int
gen_st_search_node(int * node, int load, int val)
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
gen_st_search(gen_st * t, int val)
{
  int node_depth, node_lev_idx, node_abs_idx, nodes_above;
  int * node;
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
    assert(range_end - range_start > t->b - 1);
    assert(t->pow[h - 1] - 1 < range_end - range_start && range_end - range_start <= t->pow[h] - 1);

    idx = gen_st_search_node(node, t->b - 1, val);
    idx--;

    k = ((range_end - range_start) - (t->pow[h - 1] - 1)) / ((t->b - 1) * t->pow[h - 2]);
    assert(0 <= k && k <= t->b);

    if (k == t->b) {
      range_start += (idx + 1) * t->pow[h - 1];
      range_end = range_start + t->pow[h - 1] - 1;
    } else {
      delta = (range_end - range_start) - (t->b - 1) - k * (t->pow[h - 1] - 1) - (t->b - k - 1) * (t->pow[h - 2] - 1);
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
    node_lev_idx = node_lev_idx * t->b + idx + 1;
    nodes_above += t->pow[node_depth];
    node_depth++;
    node_abs_idx = nodes_above + node_lev_idx;
    node = &t->a[node_abs_idx * (t->b - 1)];

    assert(0 <= node_lev_idx && node_lev_idx < t->pow[node_depth]);
    assert(node_depth < t->h);

    if (range_end - range_start == t->pow[h - 2] - 1)
      h -= 2;
    else
      h--;

    assert(h == 0 || node_abs_idx < t->n_nodes);
  }

  assert(range_end - range_start <= t->b - 1);

  if (h == 1) {
    idx = gen_st_search_node(node, range_end - range_start, val);
    range_start += idx;
  }

  return range_start - 1;
}


#endif
