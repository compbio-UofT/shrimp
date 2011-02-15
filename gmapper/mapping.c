#define _MODULE_MAPPING

#include <omp.h>

#include "mapping.h"
#include "../common/sw-full-common.h"
#include "../common/sw-full-cs.h"
#include "../common/sw-full-ls.h"
#include "../common/sw-vector.h"
#include "../common/read_hit_heap.h"


DEF_HEAP(uint32_t, uint, uu)
DEF_HEAP(double, struct read_hit_holder, unpaired)
DEF_HEAP(double, struct read_hit_pair_holder, paired)


/*
 * Mapping routines
 */
static void
read_get_mapidxs_per_strand(struct read_entry * re, int st)
{
  int i, sn, load, base, r_idx;
  uint32_t * kmerWindow = (uint32_t *)xcalloc(sizeof(kmerWindow[0]) * BPTO32BW(max_seed_span));
  
  re->mapidx[st] = (uint32_t *)xmalloc(n_seeds * re->max_n_kmers * sizeof(re->mapidx[0][0]));

  load = 0;
  for (i = 0; i < re->read_len; i++) {
    base = EXTRACT(re->read[st], i);
    bitfield_prepend(kmerWindow, max_seed_span, base);

    if (load < max_seed_span)
      load++;

    for (sn = 0; sn < n_seeds; sn++) {
      if (i < re->min_kmer_pos + seed[sn].span - 1)
	continue;

      r_idx = i - seed[sn].span + 1;
      re->mapidx[st][sn*re->max_n_kmers + (r_idx - re->min_kmer_pos)] = KMER_TO_MAPIDX(kmerWindow, sn);
    }
  }

  free(kmerWindow);
}


/*
 * Extract spaced kmers from read, save them in re->mapidx.
 */
static inline void
read_get_mapidxs(struct read_entry * re)
{
  read_get_mapidxs_per_strand(re, 0);
  read_get_mapidxs_per_strand(re, 1);
}


static int
bin_search(uint32_t * array, int l, int r, uint32_t value)
{
  int m;
  while (l + 1 < r) {
    m = (l + r - 1)/2;
    if (array[m] < value)
      l = m + 1;
    else
      r = m + 1;
  }
  if (l < r && array[l] < value)
    return l + 1;
  else
    return l;
}


static void
read_get_restricted_anchor_list_per_strand(struct read_entry * re, int st, bool collapse)
{
  int i, j, offset, sn;
  llint g_start, g_end;
  int idx_start, idx_end, k;
  uint32_t mapidx;
  int anchor_cache[re->read_len];
  uint diag;
  int l;

  assert(0); // unmaintained!!

  assert(re->mapidx[st] != NULL);

  re->n_anchors[st] = 0;

  if ((st == 0 && !Fflag) || (st == 1 && !Cflag))
    return;

  for (j = 0; j < re->n_ranges; j++) {
    if (re->ranges[j].st != st)
      continue;

    g_start = (llint)contig_offsets[re->ranges[j].cn] + re->ranges[j].g_start;
    g_end = (llint)contig_offsets[re->ranges[j].cn] + re->ranges[j].g_end;

    for (sn = 0; sn < n_seeds; sn++) {
      for (i = 0; re->min_kmer_pos + i + seed[sn].span - 1 < re->read_len; i++) {
	offset = sn*re->max_n_kmers + i;
	mapidx = re->mapidx[st][offset];

	idx_start = bin_search(genomemap[sn][mapidx], 0, (int)genomemap_len[sn][mapidx], g_start);
	idx_end = bin_search(genomemap[sn][mapidx], idx_start, (int)genomemap_len[sn][mapidx], g_end + 1);

	if (idx_start >= idx_end)
	  continue;

	re->anchors[st] = (struct anchor *)xrealloc(re->anchors[st],
						    sizeof(re->anchors[0][0])
						    * (re->n_anchors[st] + (idx_end - idx_start)));
	for (k = 0; idx_start + k < idx_end; k++) {
	  re->anchors[st][re->n_anchors[st] + k].cn = re->ranges[j].cn;
	  re->anchors[st][re->n_anchors[st] + k].x =
	    genomemap[sn][mapidx][idx_start + k] - contig_offsets[re->ranges[j].cn];
	  re->anchors[st][re->n_anchors[st] + k].y = re->min_kmer_pos + i;
	  re->anchors[st][re->n_anchors[st] + k].length = seed[sn].span;
	  re->anchors[st][re->n_anchors[st] + k].weight = 1;
	}
	re->n_anchors[st] += (int)(idx_end - idx_start);
      }
    }
  }

  qsort(re->anchors[st], re->n_anchors[st], sizeof(re->anchors[0][0]), anchor_uw_cmp);

  if (collapse) {
    for (i = 0; i < re->read_len; i++)
      anchor_cache[i] = -1;

    for (k = 0, i = 0; i < re->n_anchors[st]; i++) {
      re->anchors[st][k] = re->anchors[st][i];
      diag = (re->anchors[st][k].x + re->read_len - re->anchors[st][k].y) % re->read_len;
      l = anchor_cache[diag];
      if (l >= 0
	  && re->anchors[st][l].cn == re->anchors[st][k].cn
	  && anchor_uw_intersect(&re->anchors[st][l], &re->anchors[st][k])) {
	anchor_uw_join(&re->anchors[st][l], &re->anchors[st][k]);
      } else {
	anchor_cache[diag] = k;
	k++;
      }
    }
    re->n_anchors[st] = k;
  }
}


static void
read_get_anchor_list_per_strand(struct read_entry * re, int st, bool collapse)
{
  uint list_sz;
  uint offset;
  int i, sn;
  struct heap_uu h;
  uint * idx;
  struct heap_uu_elem tmp;
  int anchor_cache[re->read_len];

  assert(re->mapidx[st] != NULL);

  re->n_anchors[st] = 0;

  if ((st == 0 && !Fflag) || (st == 1 && !Cflag))
    return;

  // compute size of anchor list
  list_sz = 0;
  for (sn = 0; sn < n_seeds; sn++) {
    for (i = 0; re->min_kmer_pos + i + seed[sn].span - 1 < re->read_len; i++) {
      offset = sn*re->max_n_kmers + i;
      list_sz += genomemap_len[sn][re->mapidx[st][offset]];
    }
  }

  // init anchor list
  re->anchors[st] = (struct anchor *)xmalloc(list_sz * sizeof(re->anchors[0][0]));

  // init min heap, indices in genomemap lists, and anchor_cache
  heap_uu_init(&h, n_seeds * re->max_n_kmers);
  idx = (uint *)xcalloc(n_seeds * re->max_n_kmers * sizeof(idx[0]));
  for (i = 0; i < re->read_len; i++)
    anchor_cache[i] = -1;

  // load inital anchors in min heap
  for (sn = 0; sn < n_seeds; sn++) {
    for (i = 0; re->min_kmer_pos + i + seed[sn].span - 1 < re->read_len; i++) {
      offset = sn*re->max_n_kmers + i;

      if (genomemap_len[sn][re->mapidx[st][offset]] > list_cutoff) {
	idx[offset] = genomemap_len[sn][re->mapidx[st][offset]];
      }

      if (idx[offset] < genomemap_len[sn][re->mapidx[st][offset]]) {
	tmp.key = genomemap[sn][re->mapidx[st][offset]][idx[offset]];
	tmp.rest = offset;
	heap_uu_insert(&h, &tmp);
	idx[offset]++;
      }
    }
  }

  while (h.load > 0) {
    // extract min
    heap_uu_get_min(&h, &tmp);

    // add to anchor list
    offset = tmp.rest;
    sn = offset / re->max_n_kmers;
    i = offset % re->max_n_kmers;
    re->anchors[st][re->n_anchors[st]].x = tmp.key;
    re->anchors[st][re->n_anchors[st]].y = re->min_kmer_pos + i;
    re->anchors[st][re->n_anchors[st]].length = seed[sn].span;
    re->anchors[st][re->n_anchors[st]].width = 1;
    re->anchors[st][re->n_anchors[st]].weight = 1;
    get_contig_num(re->anchors[st][re->n_anchors[st]].x, &re->anchors[st][re->n_anchors[st]].cn);
    re->n_anchors[st]++;

    if (collapse) {
      // check if current anchor intersects the cached one on the same diagonal
      uint diag = (re->anchors[st][re->n_anchors[st]-1].x + re->read_len - re->anchors[st][re->n_anchors[st]-1].y) % re->read_len;
      int j = anchor_cache[diag];

      if (j >= 0
	  && re->anchors[st][j].cn == re->anchors[st][re->n_anchors[st]-1].cn
	  && anchor_uw_intersect(&re->anchors[st][j], &re->anchors[st][re->n_anchors[st]-1])) {
	anchor_uw_join(&re->anchors[st][j], &re->anchors[st][re->n_anchors[st]-1]);
	re->n_anchors[st]--;
      } else {
	anchor_cache[diag] = re->n_anchors[st]-1;
      }
    }

    // load next anchor for that seed/mapidx
    if (idx[offset] < genomemap_len[sn][re->mapidx[st][offset]]) {
      tmp.key = genomemap[sn][re->mapidx[st][offset]][idx[offset]];
      tmp.rest = offset;
      heap_uu_replace_min(&h, &tmp);
      idx[offset]++;
    } else {
      heap_uu_extract_min(&h, &tmp);
    }
  }

  heap_uu_destroy(&h);
  free(idx);
}


/*
 * Given the kmer lists in mapidx, lookup matching kmers in genomemap.
 * Create a list of unit-width anchors, possibly collapsing intersecting kmers.
 * Save anchor lists in re->anchors[][] and their sizes in re->n_anchors[]
 */
static inline void
read_get_anchor_list(struct read_entry * re, bool collapse) {
  if (re->n_ranges == 0) {
    read_get_anchor_list_per_strand(re, 0, collapse);
    read_get_anchor_list_per_strand(re, 1, collapse);
  } else {
    read_get_restricted_anchor_list_per_strand(re, 0, collapse);
    read_get_restricted_anchor_list_per_strand(re, 1, collapse);    
  }

#ifdef DEBUG_ANCHOR_LIST
  {
#warning Dumping anchor list.
    uint i, st;

    fprintf(stderr,"Dumping anchors for read:[%s]\n", re->name);
    for (st = 0; st < 2; st++) {
      fprintf(stderr, "st:%u ", st);
      for(i = 0; i < re->n_anchors[st]; i++){
	fprintf(stderr,"(%u,%u,%u,%u)%s", re->anchors[st][i].x, re->anchors[st][i].y,
		re->anchors[st][i].length, re->anchors[st][i].weight,
		i < re->n_anchors[st]-1? "," : "\n");
      }
    }
  }
#endif
}


#if defined (DEBUG_HIT_LIST_CREATION) || defined (DEBUG_HIT_LIST_PAIR_UP) \
  || defined (DEBUG_HIT_LIST_PASS1) || defined (DEBUG_HIT_LIST_PAIRED_HITS)
static void
dump_hit(struct read_hit * hit) {
  fprintf(stderr, "(cn:%d,st:%d,gen_st:%d,g_off:%lld,w_len:%d,scores:(%d,%d,%d),"
	          "matches:%d,pair_min:%d,pair_max:%d,anchor:(%lld,%lld,%d,%d))\n",
	  hit->cn, hit->st, hit->gen_st, hit->g_off, hit->w_len,
	  hit->score_window_gen, hit->score_vector, hit->score_full,
	  hit->matches, hit->pair_min, hit->pair_max,
	  hit->anchor.x, hit->anchor.y, hit->anchor.length, hit->anchor.width);
}
#endif


#if defined (DEBUG_HIT_LIST_CREATION) || defined (DEBUG_HIT_LIST_PAIR_UP) \
  || defined (DEBUG_HIT_LIST_PASS1)
static void
dump_hit_list(struct read_entry * re, int st, bool only_paired, bool only_after_vector)
{
  int i;
  for (i = 0; i < (int)re->n_hits[st]; i++) {
    if (only_paired && re->hits[st][i].pair_min < 0)
      continue;
    if (only_after_vector && re->hits[st][i].score_vector < 0)
      continue;

    dump_hit(&re->hits[st][i]);
  }
}
#endif


static void
read_get_hit_list_per_strand(struct read_entry * re, int match_mode, int st)
{
  llint goff, gstart, gend;
  int max_score, tmp_score = 0;
  int i, j, cn, max_idx;
  int w_len;
  int short_len = 0, long_len = 0;
  struct anchor a[3];

  re->hits[st] = (struct read_hit *)xcalloc(re->n_anchors[st] * sizeof(re->hits[0][0]));
  re->n_hits[st] = 0;

  for (i = 0; i < re->n_anchors[st]; i++) {
    // contig num of crt anchor
    cn = re->anchors[st][i].cn;

    // w_len
    w_len = re->window_len;
    if ((uint32_t)w_len > genome_len[cn])
      w_len = (int)genome_len[cn];

    // set gstart and gend
    gend = (re->anchors[st][i].x - contig_offsets[cn]) + re->read_len - 1 - re->anchors[st][i].y;
    if (gend > genome_len[cn] - 1) {
      gend = genome_len[cn] - 1;
    }

    if (gend >= re->window_len) {
      gstart = gend - re->window_len;
    } else {
      gstart = 0;
    }
    /*
     * Modes of matching:
     * 1. gapless: only extend current anchor; no window_gen_threshold check
     * 2. n=1: one kmer match & no window_gen_threshold check; or at least two matches
     * 3. n=2: at least two kmer matches & window_gen_threshold check
     */
    max_idx = i;
    max_score = re->anchors[st][i].length * match_score;

    if (!gapless_sw) {
      // avoid single matches when n=2
      if (match_mode > 1 && re->anchors[st][i].weight == 1)
	max_score = 0;

      for (j = i - 1;
	   j >= 0
	     && re->anchors[st][j].x >= (llint)contig_offsets[cn] + gstart;
	   j--) {
	if (re->anchors[st][j].y >= re->anchors[st][i].y) {
	  continue;
	}
	if (re->anchors[st][i].x - (llint)contig_offsets[cn] - re->anchors[st][i].y
	    > re->anchors[st][j].x - (llint)contig_offsets[cn] - re->anchors[st][j].y)
	  { // deletion in read
	    short_len = (int)(re->anchors[st][i].y - re->anchors[st][j].y) + re->anchors[st][i].length;
	    long_len = (int)(re->anchors[st][i].x - re->anchors[st][j].x) + re->anchors[st][i].length;
	  }
	else
	  { // insertion in read
	    short_len = (int)(re->anchors[st][i].x - re->anchors[st][j].x) + re->anchors[st][i].length;
	    long_len = (int)(re->anchors[st][i].y - re->anchors[st][j].y) + re->anchors[st][i].length;
	  }

	if (long_len > short_len) {
	  tmp_score = short_len * match_score + b_gap_open_score
	    + (long_len - short_len) * b_gap_extend_score;
	} else {
	  tmp_score = short_len * match_score;
	}

	if (tmp_score > max_score) {
	  max_idx = j;
	  max_score = tmp_score;
	}
      }
    }

    if (gapless_sw
	|| match_mode == 1
	|| max_score >= (int)abs_or_pct(window_gen_threshold,
					(re->read_len < w_len ? re->read_len : w_len) * match_score))
      {
	// set goff
	int x_len = (int)(re->anchors[st][i].x - re->anchors[st][max_idx].x) + re->anchors[st][i].length;

	if ((re->window_len - x_len)/2 < re->anchors[st][max_idx].x - contig_offsets[cn]) {
	  goff = (re->anchors[st][max_idx].x - contig_offsets[cn]) - (re->window_len - x_len)/2;
	} else {
	  goff = 0;
	}
	if (goff + w_len > genome_len[cn]) {
	  goff = genome_len[cn] - w_len;
	}

	// compute anchor
	if (max_idx < i) {
	  a[0] = re->anchors[st][i];
	  anchor_to_relative(&a[0], contig_offsets[cn] + goff);
	  a[1] = re->anchors[st][max_idx];
	  anchor_to_relative(&a[1], contig_offsets[cn] + goff);
	  anchor_join(a, 2, &a[2]);
	} else {
	  a[2] = re->anchors[st][i];
	  anchor_to_relative(&a[2], contig_offsets[cn] + goff);
	}

	// add hit
	re->hits[st][re->n_hits[st]].g_off = goff;
	re->hits[st][re->n_hits[st]].w_len = w_len;
	re->hits[st][re->n_hits[st]].cn = cn;
	re->hits[st][re->n_hits[st]].st = st;
	re->hits[st][re->n_hits[st]].gen_st = 0;
	re->hits[st][re->n_hits[st]].anchor = a[2];
	re->hits[st][re->n_hits[st]].score_window_gen = max_score;
	re->hits[st][re->n_hits[st]].matches = (gapless_sw || max_idx == i?
						re->anchors[st][i].weight :
						re->anchors[st][i].weight + re->anchors[st][max_idx].weight);
	re->hits[st][re->n_hits[st]].score_vector = -1;
	re->hits[st][re->n_hits[st]].score_full = -1;
	re->hits[st][re->n_hits[st]].score_max = (re->read_len < w_len? re->read_len : w_len) * match_score;
	re->hits[st][re->n_hits[st]].pair_min = -1;
	re->hits[st][re->n_hits[st]].pair_max = -1;
	re->n_hits[st]++;
      }
  }

  // sort list (there might be few misordered pairs because of the goff computation
  for (i = 1; i < re->n_hits[st]; i++) {
    j = i;
    while (j >= 1
	   && re->hits[st][j-1].cn == re->hits[st][i].cn
	   && re->hits[st][j-1].g_off > re->hits[st][i].g_off)
      j--;
    if (j < i) { // shift elements at indexes j..i-1 higher
      struct read_hit tmp = re->hits[st][i];
      int k;
      for (k = i - 1; k >= j; k--)
	re->hits[st][k+1] = re->hits[st][k];
      re->hits[st][j] = tmp;
    }
  }
}


/*
 * Given the list of unit-width anchors, create a set of potential hits that might be later
 * investigated by the SW vector filter.
 */
static inline void
read_get_hit_list(struct read_entry * re, int match_mode)
{
  read_get_hit_list_per_strand(re, match_mode, 0);
  read_get_hit_list_per_strand(re, match_mode, 1);

#ifdef DEBUG_HIT_LIST_CREATION
  fprintf(stderr, "Dumping hit list after creation for read:[%s]\n", re->name);
  dump_hit_list(re, 0, false, false);
  dump_hit_list(re, 1, false, false);
#endif
}


/*
 * Reverse read hit.
 *
 * The 'st' strand of the read matches the 'gen_st' strand of the genome. Negate both.
 */
static inline void
reverse_hit(struct read_entry * re, struct read_hit * rh)
{
  assert(re != NULL && rh != NULL);

  rh->g_off = genome_len[rh->cn] - rh->g_off - rh->w_len;
  anchor_reverse(&rh->anchor, rh->w_len, re->read_len);
  rh->gen_st = 1 - rh->gen_st;
  rh->st = 1 - rh->st;
}


static void
read_pass1_per_strand(struct read_entry * re, bool only_paired, int st)
{
  int i, j;

  f1_hash_tag++;
  j = -1; // last good hit
  //fprintf(stderr,"read %s, found %d hits on strand %d\n",re->name,re->n_hits[st],st);
  //int vsws=0;
  for (i = 0; i < re->n_hits[st]; i++) {
    if (only_paired && re->hits[st][i].pair_min < 0) {
      continue;
    }
    //TODO make this tunable?
    if (sam_half_paired && !only_paired && pair_mode!=PAIR_NONE) {
      if (re->hits[st][i].gen_st!=0) {
	reverse_hit(re,&re->hits[st][i]);	
      }
      if (re->hits[st][i].matches<2) {
	re->hits[st][i].score_vector=0;
	re->hits[st][i].pct_score_vector=0;
	continue;
      }
      int32_t w_len = re->window_len;
      if ((uint32_t)w_len > genome_len[re->hits[st][i].cn])
	w_len = (int)genome_len[re->hits[st][i].cn];
      if (re->hits[st][i].score_window_gen < (int)abs_or_pct(window_gen_threshold,
							     (re->read_len < w_len ? re->read_len : w_len) * match_score)) {
	re->hits[st][i].score_vector=0;
	re->hits[st][i].pct_score_vector=0;
	continue;
      }
    }

    // check window overlap
    if (j >= 0
	&& re->hits[st][i].cn == re->hits[st][j].cn
	&& re->hits[st][i].g_off <= re->hits[st][j].g_off + re->window_len
	- (int)abs_or_pct(window_overlap, re->window_len)) {
      re->hits[st][i].score_vector = 0;
      re->hits[st][i].pct_score_vector = 0;
      continue;
    }
    //fprintf(stderr,"running vector with %d matches, %lld \n",re->hits[st][i].matches,re->hits[st][i].g_off);
    //vsws++;
    if (shrimp_mode == MODE_COLOUR_SPACE)
      re->hits[st][i].score_vector = f1_run(genome_cs_contigs[re->hits[st][i].cn], genome_len[re->hits[st][i].cn],
					    re->hits[st][i].g_off, re->hits[st][i].w_len,
					    re->read[st], re->read_len,
					    re->hits[st][i].g_off + re->hits[st][i].anchor.x, re->hits[st][i].anchor.y,
					    genome_contigs[re->hits[st][i].cn], re->initbp[st], genome_is_rna, f1_hash_tag);
    else
      re->hits[st][i].score_vector = f1_run(genome_contigs[re->hits[st][i].cn], genome_len[re->hits[st][i].cn],
					    re->hits[st][i].g_off, re->hits[st][i].w_len,
					    re->read[st], re->read_len,
					    re->hits[st][i].g_off + re->hits[st][i].anchor.x, re->hits[st][i].anchor.y,
					    NULL, -1, genome_is_rna, f1_hash_tag);

    re->hits[st][i].pct_score_vector = (100 * re->hits[st][i].score_vector)/re->hits[st][i].score_max;
    if (re->hits[st][i].score_vector >= (int)abs_or_pct(sw_vect_threshold, re->hits[st][i].score_max)) {
      //save_score(re, re->hits[st][i].score, re->hits[st][i].g_off, re->hits[st][i].cn,
      //	   re->hits[st][i].st > 0, &re->hits[st][i].anchor);
      j = i;
    }
  }
  //fprintf(stderr,"Ran %d vec sw's\n",vsws);
}


/*
 * Go through hit list, apply vector filter, and save top scores.
 */
static inline void
read_pass1(struct read_entry * re, bool only_paired)
{
  read_pass1_per_strand(re, only_paired, 0);
  read_pass1_per_strand(re, only_paired, 1);

#ifdef DEBUG_HIT_LIST_PASS1
  fprintf(stderr, "Dumping hit list after pass1 for read:[%s]\n", re->name);
  dump_hit_list(re, 0, only_paired, false);
  dump_hit_list(re, 1, only_paired, false);
#endif
}


static void
readpair_pair_up_hits(struct read_entry * re1, struct read_entry * re2)
{
  int st1, st2, i, j, k, l;

  for (st1 = 0; st1 < 2; st1++) {
    st2 = 1 - st1; // opposite strand
    int re1_strand = pair_reverse[pair_mode][0] ? 1-st1 : st1;
    int re2_strand = pair_reverse[pair_mode][1] ? 1-st2 : st2;
    int re1_correction = re1->window_len-re1->read_len;
    int re2_correction = re2->window_len-re2->read_len;

    int max_correction = 0;
    int min_correction = 0;
    if (pair_mode==PAIR_OPP_IN) {
	max_correction = re1_correction + re2_correction;
    } else if (pair_mode==PAIR_OPP_OUT) {
	min_correction = - re1_correction - re2_correction;
    } else if (pair_mode==PAIR_COL_FW) {
	max_correction = re2_correction;
	min_correction = -re1_correction;
    }
   
   // printf("%d / %d i, %d , and %d / %d, i , %d\n",re1->n_hits[0],re1->n_hits[1],re1->hits[0][0].matches,re2->n_hits[0],re2->n_hits[1],re2->hits[0][0].matches  );  
    j = 0; // invariant: matching hit at index j or larger
    for (i = 0; i < re1->n_hits[st1]; i++) {
      int64_t fivep = re1->hits[st1][i].g_off + ((re1_strand==1) ? re1->window_len : 0);
      int64_t min_fivep_mp=0;
      int64_t max_fivep_mp=0;
      if (re1_strand==0) {
          min_fivep_mp = (int64_t)(fivep+min_insert_size+min_correction);
          max_fivep_mp = (int64_t)(fivep+max_insert_size+max_correction);
      } else {
          max_fivep_mp = (int64_t)(fivep-min_insert_size-min_correction);
          min_fivep_mp = (int64_t)(fivep-max_insert_size-max_correction);
      }
      //printf("MAX %ld, MIN %ld\n",max_fivep_mp,min_fivep_mp);
      // find matching hit, if any
      while (j < re2->n_hits[st2]
	     && (re2->hits[st2][j].cn < re1->hits[st1][i].cn // prev contig
		 || (re2->hits[st2][j].cn == re1->hits[st1][i].cn // same contig, but too close
		     && (int64_t)(re2->hits[st2][j].g_off + ((re2_strand==1) ? re2->window_len : 0)) < min_fivep_mp //five_mp < min_fivep_mp
		     )
		 )
	     ) {
	//printf("skipping g_off %lld, %d, translated %ld, min_fivep_mp %ld\n",re2->hits[st2][j].g_off,(re2_strand==1) ? re2->window_len : 0,(int64_t)(re2->hits[st2][j].g_off + ((re2_strand==1) ? re2->window_len : 0)),min_fivep_mp);
	j++;
      }

      k = j;
      while (k < re2->n_hits[st2]
	     && re2->hits[st2][k].cn == re1->hits[st1][i].cn
	     && (int64_t)(re2->hits[st2][k].g_off + ((re2_strand==1) ? re2->window_len : 0)) <= max_fivep_mp) {
	//printf("taking g_off %d, translated %d, max_fivep_mp %d\n",re2->hits[st2][k].g_off,(int64_t)(re2->hits[st2][k].g_off + ((re2_strand==1) ? re2->window_len : 0)),max_fivep_mp);
	k++;
      }
      if (j == k) // no paired hit
	continue;

      re1->hits[st1][i].pair_min = j;
      re1->hits[st1][i].pair_max = k-1;
      for (l = j; l < k; l++) {
	if (re2->hits[st2][l].pair_min < 0) {
	  re2->hits[st2][l].pair_min = i;
	}
	re2->hits[st2][l].pair_max = i;
      }
    }
  }

#ifdef DEBUG_HIT_LIST_PAIR_UP
  fprintf(stderr, "Dumping hit list after pair-up for read:[%s]\n", re1->name);
  dump_hit_list(re1, 0, false, false);
  dump_hit_list(re1, 1, false, false);
  fprintf(stderr, ".. and read:[%s]\n", re2->name);
  dump_hit_list(re2, 0, false, false);
  dump_hit_list(re2, 1, false, false);
#endif
}


/*
 * Go through the list adding hits to the heap.
 */
static void
read_get_vector_hits(struct read_entry * re, struct heap_unpaired * h)
{
  int st, i;
  heap_unpaired_elem tmp;

  assert(re != NULL && h != NULL);
  assert(pair_mode == PAIR_NONE || sam_half_paired);

  for (st = 0; st < 2; st++) {
    assert(re->n_hits[st] == 0 || re->hits[st] != NULL);

    for (i = 0; i < re->n_hits[st]; i++) {
      if (re->hits[st][i].score_vector >= (int)abs_or_pct(sw_vect_threshold, re->hits[st][i].score_max)
	  && (h->load < h->capacity
	      || ( (IS_ABSOLUTE(sw_vect_threshold)
		    && re->hits[st][i].score_vector > (int)h->array[0].key)
		   || (!IS_ABSOLUTE(sw_vect_threshold)
		       && re->hits[st][i].pct_score_vector > (int)h->array[0].key)))) {
        //fprintf(stderr,"%d matches\n",re->hits[st][i].matches);
	tmp.key = (IS_ABSOLUTE(sw_vect_threshold)? re->hits[st][i].score_vector : re->hits[st][i].pct_score_vector);
	tmp.rest.hit = &re->hits[st][i];
        assert(re->hits[st][i].gen_st==0);
	if (h->load < h->capacity)
	  heap_unpaired_insert(h, &tmp);
	else
	  heap_unpaired_replace_min(h, &tmp);
      }
    }
  }
}


/*
 * Go through the hit lists, constructing paired hits.
 */
static void
readpair_get_vector_hits(struct read_entry * re1, struct read_entry * re2, struct heap_paired * h)
{
  int st1, st2, i, j;
  heap_paired_elem tmp;

  assert(re1 != NULL && re2 != NULL && h != NULL);
  assert(pair_mode != PAIR_NONE);

  for (st1 = 0; st1 < 2; st1++) {
    st2 = 1 - st1; // opposite strand

    for (i = 0; i < re1->n_hits[st1]; i++) {
      if (re1->hits[st1][i].pair_min < 0)
	continue;

      for (j = re1->hits[st1][i].pair_min; j <= re1->hits[st1][i].pair_max; j++) {
	if (num_matches == 3 // require at least two matches on one of the feet
	    && re1->hits[st1][i].matches == 1
	    && re2->hits[st2][j].matches == 1)
	  continue;

	if (re1->hits[st1][i].score_vector + re2->hits[st2][j].score_vector
	    >= (int)abs_or_pct(sw_vect_threshold, re1->hits[st1][i].score_max + re2->hits[st2][j].score_max)
	    && (h->load < h->capacity
		|| ( (IS_ABSOLUTE(sw_vect_threshold)
		      && re1->hits[st1][i].score_vector + re2->hits[st2][j].score_vector > (int)h->array[0].key)
		     || (!IS_ABSOLUTE(sw_vect_threshold)
			 && (re1->hits[st1][i].pct_score_vector + re2->hits[st2][j].pct_score_vector)/2 > (int)h->array[0].key)))) {
	  tmp.key = (IS_ABSOLUTE(sw_vect_threshold)?
		     re1->hits[st1][i].score_vector + re2->hits[st2][j].score_vector
		     : (re1->hits[st1][i].pct_score_vector + re2->hits[st2][j].pct_score_vector)/2);
	  tmp.rest.hit[0] = &re1->hits[st1][i];
	  tmp.rest.hit[1] = &re2->hits[st2][j];
	//TODO HISTORGRAM IS OFF! NOT SAM
	  tmp.rest.insert_size = (int)(st1 == 0?
				       re2->hits[st2][j].g_off - (re1->hits[st1][i].g_off + re1->hits[st1][i].w_len) :
				       re1->hits[st1][i].g_off - (re2->hits[st2][j].g_off + re2->hits[st2][j].w_len));

	  if (h->load < h->capacity)
	    heap_paired_insert(h, &tmp);
	  else
	    heap_paired_replace_min(h, &tmp);
	}
      }
    }
  }
}


/*
 * Run full SW filter on this hit.
 */
static void
hit_run_full_sw(struct read_entry * re, struct read_hit * rh, int thresh)
{
  uint32_t * gen = NULL;

  assert(re != NULL && rh != NULL);
  if (rh->gen_st!=0) {
	fprintf(stderr,"rh->gen_st is %d\n%d and %d\n",rh->gen_st,min_insert_size,max_insert_size);
  }
  assert(rh->gen_st == 0);

  if (rh->st == re->input_strand) {
    gen = genome_contigs[rh->cn];
  } else {
    reverse_hit(re, rh);
    gen = genome_contigs_rc[rh->cn];
  }

  assert(rh->st == re->input_strand);

  rh->sfrp = (struct sw_full_results *)xcalloc(sizeof(rh->sfrp[0]));

#ifdef DEBUG_SW_FULL_CALLS
  fprintf(stderr, "SW full call: (name:[%s],cn:%d,st:%d,gen_st:%d,g_off:%lld,w_len:%d,anchor:(%lld,%lld,%d,%d))\n",
	  re->name, rh->cn, rh->st, rh->gen_st, rh->g_off, rh->w_len,
	  rh->anchor.x, rh->anchor.y, rh->anchor.length, rh->anchor.width);
#endif

  if (shrimp_mode == MODE_COLOUR_SPACE) {
    sw_full_cs(gen, rh->g_off, rh->w_len,
	       re->read[rh->st], re->read_len, re->initbp[rh->st],
	       thresh, rh->sfrp, rh->gen_st && Tflag, genome_is_rna,
	       &rh->anchor, 1,Gflag ? 0 : 1);
  } else {
    /*
     * The full SW in letter space assumes it's given the correct max score.
     * This might not be true just yet if we're using hashing&caching because
     * of possible hash collosions.
     */
    rh->score_vector = sw_vector(gen, rh->g_off, rh->w_len,
				 re->read[rh->st], re->read_len,
				 NULL, -1, genome_is_rna);

    if (rh->score_vector >= thresh) {
      sw_full_ls(gen, rh->g_off, rh->w_len,
		 re->read[rh->st], re->read_len,
		 thresh, rh->score_vector, rh->sfrp, rh->gen_st && Tflag,
		 &rh->anchor, 1, Gflag ? 0 : 1);
      //assert(rh->sfrp->score == rh->score_vector);
    } else { // this wouldn't have passed the filter
      rh->sfrp->score = 0;
    }
  }
  rh->score_full = rh->sfrp->score;
  rh->pct_score_full = (100 * rh->score_full)/(double)rh->score_max;
}


static inline int32_t
get_isize(struct read_hit *rh, struct read_hit * rh_mp)
{
  if (rh_mp==NULL || rh==NULL) {
    return 0;
  }
  //get the mate pair info
  int read_start_mp = rh_mp->sfrp->read_start+1; //1based
  int read_end_mp = read_start_mp + rh_mp->sfrp->rmapped -1; //1base
  int genome_length_mp = genome_len[rh_mp->cn];
  int genome_start_mp;
  bool reverse_strand_mp = (rh_mp->gen_st ==1);
  if (!reverse_strand_mp) {
    genome_start_mp = rh_mp->sfrp->genome_start+1; // 0 based -> 1 based
  } else {
    int genome_right_most_coordinate = genome_length_mp - rh_mp->sfrp->genome_start;
    //rh->sfrp->deletions is deletions in the reference
    // This is when the read has extra characters that dont match into ref
    genome_start_mp = genome_right_most_coordinate - (read_end_mp - read_start_mp - rh_mp->sfrp->deletions + rh_mp->sfrp->insertions);
  }
  int genome_end_mp=genome_start_mp+rh_mp->sfrp->gmapped-1;
  //get the other hit info
  int genome_start;
  bool reverse_strand=(rh->gen_st == 1);
  if (!reverse_strand) {
    genome_start = rh->sfrp->genome_start+1; // 0 based -> 1 based
  } else {
    int read_start = rh->sfrp->read_start+1; //1based
    int read_end = read_start + rh->sfrp->rmapped -1; //1base
    int genome_length = genome_len[rh->cn];
    int genome_right_most_coordinate = genome_length - rh->sfrp->genome_start;
    //rh->sfrp->deletions is deletions in the reference
    // This is when the read has extra characters that dont match into ref
    genome_start = genome_right_most_coordinate - (read_end - read_start - rh->sfrp->deletions + rh->sfrp->insertions);
  }
  int genome_end=genome_start+rh->sfrp->gmapped-1;
	
  int fivep = 0;
  int fivep_mp = 0;
  if (reverse_strand){
    fivep = genome_end;
  } else {
    fivep = genome_start - 1;
  }

  if (reverse_strand_mp){
    fivep_mp = genome_end_mp;
  } else {
    fivep_mp = genome_start_mp-1;
  }
  return (fivep_mp - fivep);
}


/*
 * Do a final pass for given read.
 * Highest scoring matches are in scores heap.
 */
static int
read_pass2(struct read_entry * re, struct heap_unpaired * h) {
  int i;
  assert(re != NULL && h != NULL);
  assert(pair_mode == PAIR_NONE || sam_half_paired);

  read_hit_heap rhh;
  read_hit_heap_init(&rhh,num_tmp_outputs,true);
  read_hit_heap_e tmp;

  /* compute full alignment scores */
  for (i = 0; i <(int)h->load; i++) {
    struct read_hit * rh = h->array[i].rest.hit;
    hit_run_full_sw(re, rh, (int)abs_or_pct(sw_full_threshold, rh->score_max));
    if ( (IS_ABSOLUTE(sw_full_threshold)
		&& rh->score_full >= abs_or_pct(sw_full_threshold, rh->score_max))
		|| (!IS_ABSOLUTE(sw_full_threshold) && rh->score_full*100/(double)rh->score_max >= sw_full_threshold) ) {
	    tmp.score=rh->score_full;
	    tmp.isize_score=0;
	    tmp.isize=0;
	    tmp.max_score=rh->score_max;
	    tmp.rh[0]=rh;
	    tmp.rh[1]=NULL;
	    rhh.array[rhh.load++]=tmp;
    } 
  }

  if (strata_flag) {
	read_hit_heap_remove_dups_and_sort_strata(&rhh);
  } else {
	read_hit_heap_remove_dups_and_sort(&rhh);
  }

  int hits[] = {0,0,0};
  for (i=0; i<rhh.load; i++) {
	struct sw_full_results * sfrp = rhh.array[i].rh[0]->sfrp;
	int edit_distance=sfrp->mismatches+sfrp->insertions+sfrp->deletions;
	if (0<=edit_distance && edit_distance<3) {
		hits[edit_distance]++;
	}
  }


  int satisfying_alignments=rhh.load;
  if (rhh.load>0) {
	  if (max_alignments==0 || satisfying_alignments<=max_alignments) {
	#pragma omp atomic
	    total_reads_matched++;
	  } else {
	#pragma omp atomic
	    total_reads_dropped++;
	    rhh.load=0;
	  }
  }
  rhh.load=MIN(num_outputs,rhh.load);
  int outputted_alignments=rhh.load;
  re->final_matches+=outputted_alignments;
	  /* Output sorted list, removing any duplicates. */
	  for (i = 0; i<rhh.load; i++) {
	    struct read_hit * rh = rhh.array[i].rh[0];
	    char * output1 = NULL, * output2 = NULL;
	    if (sam_half_paired) {
		int other_hits[]={0,0,0};
		if (re->first_in_pair) {
	      		hit_output(re, rh, NULL, &output1, NULL ,true,hits,satisfying_alignments);
	      		hit_output(re->mate_pair, NULL, rh, &output2, NULL,false,other_hits,0);
		} else {
	      		hit_output(re->mate_pair, NULL, rh, &output1, NULL,true,other_hits,0);
	      		hit_output(re, rh, NULL, &output2, NULL,false,hits,satisfying_alignments);
		}
	    } else {
	      	hit_output(re, rh, NULL,  &output1, &output2, false, hits,satisfying_alignments);
	    }
	    //legacy support for SHRiMP output
	    if (!Eflag) {
		if ( !Pflag) {
		#pragma omp critical (stdout)
			{
			  fprintf(stdout, "%s\n", output1);
			}
		      } else {
		#pragma omp critical (stdout)
			{
			  fprintf(stdout, "%s\n\n%s\n", output1, output2);
			}
		      }
		      free(output1);
		      if (Pflag || sam_half_paired) {
			free(output2);
		      }
	    }


	  }
#pragma omp atomic
  total_single_matches += re->final_matches;

  read_hit_heap_destroy(&rhh);
  return outputted_alignments;
}


/*
 * Do a final pass for given pair of reads.
 * Highest scoring matches are in scores heap.
 */
static int
readpair_pass2(struct read_entry * re1, struct read_entry * re2, struct heap_paired * h) {
  int i, j;

  assert(re1 != NULL && re2 != NULL && h != NULL);
  read_hit_heap rhh;
  read_hit_heap_init(&rhh,num_tmp_outputs,true);
  read_hit_heap_e tmp;
  /* compute full alignment scores */
  for (i = 0; i < (int)h->load; i++) {
    for (j = 0; j < 2; j++) {
      struct read_hit * rh = h->array[i].rest.hit[j];
      struct read_entry * re = (j == 0? re1 : re2);

      if (rh->score_full >= 0) // previously insepcted foot
	continue;

      hit_run_full_sw(re, rh, (int)abs_or_pct(sw_full_threshold, h->array[i].rest.hit[0]->score_max + h->array[i].rest.hit[1]->score_max)/3);
    }
    if (h->array[i].rest.hit[0]->score_full == 0 || h->array[i].rest.hit[1]->score_full == 0) {
      h->array[i].key = 0;
    } else {
      int32_t score=h->array[i].rest.hit[0]->score_full+h->array[i].rest.hit[1]->score_full; 
      int32_t max_score=h->array[i].rest.hit[0]->score_max+h->array[i].rest.hit[1]->score_max;
      if  ( (IS_ABSOLUTE(sw_full_threshold) && score >= abs_or_pct(sw_full_threshold,h->array[i].rest.hit[0]->score_max + h->array[i].rest.hit[1]->score_max))
	      || (!IS_ABSOLUTE(sw_full_threshold) && score*100/(double)max_score >= sw_full_threshold) ) {
		tmp.score=score;
		tmp.isize=abs(get_isize(h->array[i].rest.hit[0],h->array[i].rest.hit[1]));
		tmp.isize_score=expected_isize==-1 ? 0 : abs(tmp.isize-expected_isize);	
		tmp.max_score=max_score;
		tmp.rh[0]=h->array[i].rest.hit[0];
		tmp.rh[1]=h->array[i].rest.hit[1];
		rhh.array[rhh.load++]=tmp;
	}
    }
  }


  if (strata_flag) {
        read_hit_heap_remove_dups_and_sort_strata(&rhh);
  } else {
        read_hit_heap_remove_dups_and_sort(&rhh);
  }


  int hits1[] = {0,0,0};
  int hits2[] = {0,0,0};
  
  for (i = 0; i<rhh.load; i++) {
	struct sw_full_results * sfrp1 = rhh.array[i].rh[0]->sfrp;
	struct sw_full_results * sfrp2 = rhh.array[i].rh[1]->sfrp;
	int edit_distance1=sfrp1->mismatches+sfrp1->insertions+sfrp1->deletions;
	int edit_distance2=sfrp2->mismatches+sfrp2->insertions+sfrp1->deletions;
	if (0<=edit_distance1 && edit_distance1<3) {
		hits1[edit_distance1]++;
	}
	if (0<=edit_distance2 && edit_distance2<3) {
		hits2[edit_distance2]++;
	}
  }
	
  int satisfying_alignments=rhh.load;
  if (rhh.load>0) {
    if (max_alignments==0 || satisfying_alignments<=max_alignments) {
#pragma omp atomic
    	total_pairs_matched++;
    } else {
#pragma omp atomic
	total_pairs_dropped++;
	rhh.load=0;
    }
  }

  rhh.load=MIN(num_outputs,rhh.load);
  int outputted_alignments=rhh.load;

  /* Output sorted list, removing any duplicates. */
	  for (i = 0; i<rhh.load; i++) {
	    struct read_hit * rh1 = rhh.array[i].rh[0];
	    struct read_hit * rh2 = rhh.array[i].rh[1];
	    uint bucket;

	      char * output1 = NULL, * output2 = NULL, * output3 = NULL, * output4 = NULL;

	      re1->final_matches++;
	      assert(rh1!=NULL && rh2!=NULL);
	      hit_output(re1, rh1,  rh2, &output1, &output2,true,hits1,satisfying_alignments);
	      hit_output(re2, rh2,  rh1, &output3, &output4,false,hits2,satisfying_alignments);
	      if (!Eflag) {
		      if (!Pflag) {
		#pragma omp critical (stdout)
			{
			  fprintf(stdout, "%s\n%s\n", output1, output3);
			}
		      } else {
		#pragma omp critical (stdout)
			{
			  fprintf(stdout, "%s\n%s\n%s\n%s\n", output1, output2, output3, output4);
			}
		      }
	      }

	      if (Xflag) {
		if (rhh.array[i].isize < min_insert_size)
		  bucket = 0;
		else if (rhh.array[i].isize > max_insert_size)
		  bucket = 99;
		else
		  bucket = (uint)((rhh.array[i].isize - min_insert_size) / insert_histogram_bucket_size);

		if (bucket >= 100) {
		  fprintf(stderr, "error: insert_size:%d re1->name:(%s)\n", rhh.array[i].isize, re1->name);
		}
		assert(bucket < 100);

	#pragma omp atomic
		insert_histogram[bucket]++;
	      }

	      if (Pflag) {
	        free(output1);
	        free(output3);
	      	free(output2);
	      	free(output4);
	      }
	  }

	#pragma omp atomic
	  total_paired_matches += re1->final_matches;
   read_hit_heap_destroy(&rhh);
   return outputted_alignments;
}


void
handle_read(read_entry *re)
{
  uint i;
  uint64_t before = rdtsc();

  read_get_mapidxs(re);
#ifdef DEBUG_KMERS
  {
    uint sn, i, j;
    fprintf(stderr, "max_n_kmers:%u, min_kmer_pos:%u\n",
	    re->max_n_kmers, re->min_kmer_pos);
    fprintf(stderr, "collapsed kmers from read:\n");
    for (sn = 0; sn < n_seeds; sn++) {
      fprintf(stderr, "sn:%u\n", sn);
      for (i = 0; re->min_kmer_pos + i + seed[sn].span <= re->read_len; i++) {
	fprintf(stderr, "\tpos:%u kmer:", re->min_kmer_pos + i);
	for (j = 0; j < seed[sn].weight; j++) {
	  fprintf(stderr, "%c%s",
		  base_translate((re->mapidx[0][sn*re->max_n_kmers + i] >> 2*(seed[sn].weight - 1 - j)) & 0x3,
				 shrimp_mode == MODE_COLOUR_SPACE),
		  j < seed[sn].weight - 1? "," : "\n");
	}
      }
    }
    fprintf(stderr, "collapsed kmers from read_rc:\n");
    for (sn = 0; sn < n_seeds; sn++) {
      fprintf(stderr, "sn:%u\n", sn);
      for (i = 0; re->min_kmer_pos + i + seed[sn].span <= re->read_len; i++) {
	fprintf(stderr, "\tpos:%u kmer:", re->min_kmer_pos + i);
	for (j = 0; j < seed[sn].weight; j++) {
	  fprintf(stderr, "%c%s",
		  base_translate((re->mapidx[1][sn*re->max_n_kmers + i] >> 2*(seed[sn].weight - 1 - j)) & 0x3,
				 shrimp_mode == MODE_COLOUR_SPACE),
		  j < seed[sn].weight - 1? "," : "\n");
	}
      }
    }
  }
#endif

  read_get_anchor_list(re, true);

  read_get_hit_list(re, (num_matches >= 2? 2 : 1));
  read_pass1(re, false);

  // initialize heap of best hits for this read
  heap_unpaired h;
  heap_unpaired_init(&h, num_tmp_outputs);
  read_get_vector_hits(re, &h);

  DEBUG("second pass");
  int printed_alignments=0;
  if (h.load > 0) {
    printed_alignments = read_pass2(re, &h);
    if (aligned_reads_file!=NULL) {
#pragma omp critical (stdout) 
      {
	fasta_write_read(aligned_reads_file,re);
      }
    }
  }
  if (printed_alignments==0) {
    if (Eflag && sam_unaligned) {
      //no alignments, print to sam empty record
      hit_output(re, NULL, NULL, NULL, NULL, false, NULL,0);
    }
    if (unaligned_reads_file!=NULL) {
#pragma omp critical (stdout) 
      {
	fasta_write_read(unaligned_reads_file,re);
      }
    }
  }


  // Done with this read; deallocate memory.
  for (i = 0; i < h.load; i++)
    hit_free_sfrp(h.array[i].rest.hit);
  heap_unpaired_destroy(&h);

  read_free_full(re);

  scan_ticks[omp_get_thread_num()] += rdtsc() - before;
}


void
handle_readpair(struct read_entry * re1, struct read_entry * re2)
{
  heap_paired h;
  uint i;
  uint64_t before = rdtsc();

  read_get_mapidxs(re1);
  read_get_mapidxs(re2);

  read_get_anchor_list(re1, true);
  read_get_anchor_list(re2, true);

  read_get_hit_list(re1, (num_matches >= 4? 2 : 1));
  read_get_hit_list(re2, (num_matches >= 4? 2 : 1));

  readpair_pair_up_hits(re1, re2);

  read_pass1(re1, true); // check pairing
  read_pass1(re2, true);

  heap_paired_init(&h, num_tmp_outputs);

  readpair_get_vector_hits(re1, re2, &h);

  DEBUG("second pass");
  int printed_alignments=0;
  if (h.load > 0) {
    printed_alignments+=readpair_pass2(re1, re2, &h);
    if (aligned_reads_file!=NULL) {
#pragma omp critical (stdout)
      {
	fasta_write_read(aligned_reads_file,re1);
	fasta_write_read(aligned_reads_file,re2);
      }
    }
  }
  if (printed_alignments==0 || (max_alignments>0 && printed_alignments>max_alignments)) {
    if (sam_half_paired) {
      //fprintf(stderr,"in half paired mode\n");
      heap_unpaired h;
      unsigned int i;

      //free(re1->hits[0]);
      //free(re1->hits[1]);
      //read_get_hit_list(re1, (num_matches>1? 2 : 1));
      read_pass1(re1, false);
      heap_unpaired_init(&h, num_tmp_outputs);
      read_get_vector_hits(re1, &h);
      if (h.load>0) {
	printed_alignments+=read_pass2(re1, &h);	
      }
      // Done with this read; deallocate memory.
      for (i = 0; i < h.load; i++)
	hit_free_sfrp(h.array[i].rest.hit);
      heap_unpaired_destroy(&h);
 
      //free(re2->hits[0]);
      //free(re2->hits[1]);
      //read_get_hit_list(re2, (num_matches>1? 2 : 1));
      read_pass1(re2, false);
      heap_unpaired_init(&h, num_tmp_outputs);
      read_get_vector_hits(re2, &h);
      if (h.load>0) {
	printed_alignments+=read_pass2(re2, &h);	
      }
      // Done with this read; deallocate memory.
      for (i = 0; i < h.load; i++)
	hit_free_sfrp(h.array[i].rest.hit);
      heap_unpaired_destroy(&h);
    }
    if (Eflag && sam_unaligned && printed_alignments==0) {
      //no alignments, print to sam empty record
      hit_output(re1, NULL, NULL, NULL, NULL, true, NULL,0);
      hit_output(re2, NULL, NULL, NULL, NULL, false, NULL,0);

    }
    if (unaligned_reads_file!=NULL) {
#pragma omp critical (stdout)
      {
	fasta_write_read(unaligned_reads_file,re1);
	fasta_write_read(unaligned_reads_file,re2);
      }
    }
  }

  /* Done; free read entry */
  for (i = 0; i < h.load; i++) {
    hit_free_sfrp(h.array[i].rest.hit[0]);
    hit_free_sfrp(h.array[i].rest.hit[1]);
  }

  heap_paired_destroy(&h);

  read_free_full(re1);
  read_free_full(re2);

  scan_ticks[omp_get_thread_num()] += rdtsc() - before;
}
