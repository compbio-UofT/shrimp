#define _MODULE_MAPPING

#include <omp.h>
#ifdef USE_PREFETCH
#include <xmmintrin.h>
#endif

#include "mapping.h"
#include "output.h"
#include "../common/sw-full-common.h"
#include "../common/sw-full-cs.h"
#include "../common/sw-full-ls.h"
#include "../common/sw-vector.h"
#include "../common/read_hit_heap.h"
#include "../common/sw-post.h"

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

  assert(re != NULL);
  assert(re->mapidx[st] == NULL);
  
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
  fprintf(stderr, "(cn:%d,st:%d,gen_st:%d,g_off:%lld,w_len:%d,scores:(wg=%d,vc=%d,fl=%d,poster=%.5g),"
	          "matches:%d,pair_min:%d,pair_max:%d,anchor:(x=%lld,y=%lld,ln=%d,wd=%d))\n",
	  hit->cn, hit->st, hit->gen_st, hit->g_off, hit->w_len,
	  hit->score_window_gen, hit->score_vector, hit->score_full, hit->sfrp != NULL? hit->sfrp->posterior : -1.0,
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
    if (re->hits[st][i].score_vector >= 0) {
      //fprintf(stderr, "found precomputed sw score for hit [%s]; skipping\n", re->name);
      continue;
    }
    //TODO make this tunable?
    if (sam_half_paired && !only_paired && pair_mode!=PAIR_NONE) {
      //fprintf(stderr, "got here, [%s]\n", re->name);
      if (re->hits[st][i].gen_st!=0) {
	reverse_hit(re,&re->hits[st][i]);
	fprintf(stderr, "reverse_hit??? XXX\n");
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
					    genome_contigs[re->hits[st][i].cn], re->initbp[st], genome_is_rna, f1_hash_tag,
					    gapless_sw);
    else
      re->hits[st][i].score_vector = f1_run(genome_contigs[re->hits[st][i].cn], genome_len[re->hits[st][i].cn],
					    re->hits[st][i].g_off, re->hits[st][i].w_len,
					    re->read[st], re->read_len,
					    re->hits[st][i].g_off + re->hits[st][i].anchor.x, re->hits[st][i].anchor.y,
					    NULL, -1, genome_is_rna, f1_hash_tag,
					    gapless_sw);

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

    j = 0; // invariant: matching hit at index j or larger
    for (i = 0; i < re1->n_hits[st1]; i++) {

      // find matching hit, if any
      while (j < re2->n_hits[st2]
	     && (re2->hits[st2][j].cn < re1->hits[st1][i].cn // prev contig
		 || (re2->hits[st2][j].cn == re1->hits[st1][i].cn // same contig, but too close
		     && (int64_t)(re2->hits[st2][j].g_off ) < (int64_t)(re1->hits[st1][i].g_off) + (int64_t)re1->delta_g_off_min[st1]
		     )
		 )
	     )
	{
	  j++;
	}

      k = j;
      while (k < re2->n_hits[st2]
	     && re2->hits[st2][k].cn == re1->hits[st1][i].cn
	     && (int64_t)(re2->hits[st2][k].g_off) <= (int64_t)(re1->hits[st1][i].g_off) + (int64_t)re1->delta_g_off_max[st1]
	     )
	{
	  k++;
	}
      //fprintf(stderr,"DONE\n");
      if (j == k) {
	//fprintf(stderr,"no paired hit\n");
	continue;
      }

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
        //assert(re->hits[st][i].gen_st==0);
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

  /*
  if (rh->gen_st!=0) {
    fprintf(stderr,"[%s] rh->gen_st is %d\n%d and %d\n",re->name, rh->gen_st,min_insert_size,max_insert_size);
  }
  assert(rh->gen_st == 0);

  if (rh->st == re->input_strand) {
    gen = genome_contigs[rh->cn];
  } else {
    reverse_hit(re, rh);
    //fprintf(stderr, "reverse_hit from hit_run_full_sw [%s]\n", re->name);
    gen = genome_contigs_rc[rh->cn];
  }
  */

  if (rh->st != re->input_strand) {
    reverse_hit(re, rh);
  }

  if (rh->gen_st == 0) {
    gen = genome_contigs[rh->cn];
  } else {
    gen = genome_contigs_rc[rh->cn];
  }

  assert(rh->sfrp == NULL);
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
  rh->pct_score_full = (1000 * 100 * rh->score_full)/rh->score_max;
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
    if (rh->score_full < 0) {
      hit_run_full_sw(re, rh, (int)abs_or_pct(sw_full_threshold, rh->score_max));
    }
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
    if (aligned_reads_file!=NULL) { // ???!!! JUST BECAUSE READ HAS VECTOR HITS DOESN'T MEAN IT HAS FINAL ALIGNMENTS
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
    free_sfrp(&h.array[i].rest.hit->sfrp);
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
    if (aligned_reads_file!=NULL) { // ???!!! SAME PROBLEM EXPLAINED ABOVE
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
	free_sfrp(&h.array[i].rest.hit->sfrp);
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
	free_sfrp(&h.array[i].rest.hit->sfrp);
      heap_unpaired_destroy(&h);
    }
    if (Eflag && sam_unaligned && printed_alignments==0) {
      //no alignments, print to sam empty record
      hit_output(re1, NULL, NULL, NULL, NULL, true, NULL,0);
      hit_output(re2, NULL, NULL, NULL, NULL, false, NULL,0);

    }
    if (unaligned_reads_file!=NULL) { // ???!!! THE READS MIGHT BOTH HAVE UNPAIRED ALIGNMENTS. SHOULD THEY BE CONSIDERED UNALIGNED?
#pragma omp critical (stdout)
      {
	fasta_write_read(unaligned_reads_file,re1);
	fasta_write_read(unaligned_reads_file,re2);
      }
    }
  }

  /* Done; free read entry */
  for (i = 0; i < h.load; i++) {
    free_sfrp(&h.array[i].rest.hit[0]->sfrp);
    free_sfrp(&h.array[i].rest.hit[1]->sfrp);
  }

  heap_paired_destroy(&h);

  read_free_full(re1);
  read_free_full(re2);

  scan_ticks[omp_get_thread_num()] += rdtsc() - before;
}


static void
read_get_region_counts(struct read_entry * re, int st, struct regions_options * options)
{
  int sn, i, offset, region;
  uint j;
  llint before = rdtsc();
  int number_in_pair = re->first_in_pair? 0 : 1;

  assert(use_regions);

  if (region_map_id == 0) {
    region_map_id = 1;
    for (int _nip = 0; _nip < 2; _nip++) {
      for (int _st = 0; _st < 2; _st++) {
	free(region_map[_nip][_st]);
	region_map[_nip][_st] = (int32_t *)xcalloc(n_regions * sizeof(region_map[0][0][0]));
      }
    }
  }

  assert(region_map[0][0] != NULL);

  for (sn = (options->min_seed >= 0? options->min_seed : 0);
       sn <= (options->max_seed >= 0? options->max_seed: n_seeds - 1);
       sn++) {
    for (i = 0; re->min_kmer_pos + i + seed[sn].span - 1 < re->read_len; i++) {
      offset = sn*re->max_n_kmers + i;

      if (genomemap_len[sn][re->mapidx[st][offset]] > list_cutoff)
        continue;

      for (j = 0; j < genomemap_len[sn][re->mapidx[st][offset]]; j++) {
#ifdef USE_PREFETCH
	if (j + 4 < genomemap_len[sn][re->mapidx[st][offset]]) {
	  int region_ahead = (int)(genomemap[sn][re->mapidx[st][offset]][j + 4] >> region_bits);
	  _mm_prefetch((char *)&region_map[number_in_pair][st][region_ahead], _MM_HINT_T0);
	}
#endif

        region = (int)(genomemap[sn][re->mapidx[st][offset]][j] >> region_bits);

        if (((region_map[number_in_pair][st][region] >> 16) & ((1 << region_map_id_bits) - 1)) == region_map_id) {
	  if ((region_map[number_in_pair][st][region] & ((1 << 8) - 1)) < region_map_max_count) // don't overflow
	    region_map[number_in_pair][st][region]++;
	} else {
	  region_map[number_in_pair][st][region] = 0; // clear old entry
	  region_map[number_in_pair][st][region] |= (region_map_id << 16); // set new id
	  region_map[number_in_pair][st][region]++; // set new count
	}

	// extend regions by region_overlap
	if ((genomemap[sn][re->mapidx[st][offset]][j] & ((1 << region_bits) - 1)) < (uint)region_overlap && region > 0) {
	  region--;

	  // BEGIN copy-paste from above
	  if (((region_map[number_in_pair][st][region] >> 16) & ((1 << region_map_id_bits) - 1)) == region_map_id) {
	    if ((region_map[number_in_pair][st][region] & ((1 << 8) - 1)) < region_map_max_count) // don't overflow
	      region_map[number_in_pair][st][region]++;
	  } else {
	    region_map[number_in_pair][st][region] = 0; // clear old entry
	    region_map[number_in_pair][st][region] |= (region_map_id << 16); // set new id
	    region_map[number_in_pair][st][region]++; // set new count
	  }
	  // END
	}
      }
    }
  }

  region_counts_ticks[omp_get_thread_num()] += rdtsc() - before;
}


/*
static void
read_get_mp_region_counts_per_strand(struct read_entry * re, int st, struct read_entry * re_mp,
				     struct regions_options * options)
{
  int i, first, last, j, max;

  for (i = 0; i < n_regions; i++) {
    max = 0;
    first = i + re->delta_region_min[st];
    if (first < 0)
      first = 0;
    last = i + re->delta_region_max[st];
    if (last > n_regions - 1)
      last = n_regions - 1;
    for (j = first; j <= last; j++) {
      if (re_mp->region_map[1-st][0][j] == re_mp->region_map_id
	  && re_mp->region_map[1-st][1][j] > max)
	max = re_mp->region_map[1-st][1][j];
    }
    re->region_map[st][2][i] = max;
  }
}


static inline void
read_get_mp_region_counts(struct read_entry * re, struct read_entry * re_mp,
			  struct regions_options * options)
{
  read_get_mp_region_counts_per_strand(re, 0, re_mp, options);
  read_get_mp_region_counts_per_strand(re, 1, re_mp, options);
}
*/


static inline void
advance_index_in_genomemap(struct read_entry * re, int st,
			   struct anchor_list_options * options,
			   uint * idx, uint max_idx, uint32_t * map, int * anchors_discarded)
{
  int first, last, max, k;
  int nip = re->first_in_pair? 0 : 1;

  while (*idx < max_idx) {
#ifdef USE_PREFETCH
    if (*idx + 2 < max_idx) {
      int region_ahead = (int)(map[*idx + 2] >> region_bits);
      _mm_prefetch((char *)&region_map[number_in_pair][st][region_ahead], _MM_HINT_T0);
    }
#endif
    int region = (int)(map[*idx] >> region_bits);

    assert(((region_map[nip][st][region] >> 16) & ((1 << region_map_id_bits) - 1)) == region_map_id);

    if ((options->min_count[0] == 0 || options->min_count[0] <= (region_map[nip][st][region] & ((1 << 8) - 1)))
	&& (options->max_count[0] == 0 || options->max_count[0] >= (region_map[nip][st][region] & ((1 << 8) - 1))))
      {
	if (options->min_count[1] != 0 && options->max_count[1] != 0)
	  break;
	if ((region_map[nip][st][region] >> 31) == 0) {
	  // compute mp counts
	  first = MAX(0, region + re->delta_region_min[st]);
	  last = MIN(n_regions - 1, region + re->delta_region_max[st]);
	  max = 0;
	  for (k = first; k <= last; k++) {
	    if (((region_map[1-nip][1-st][k] >> 16) & ((1 << region_map_id_bits) - 1)) == region_map_id
		&& (region_map[1-nip][1-st][k] & ((1 << 8) - 1)) > max) {
	      max = (region_map[1-nip][1-st][k] & ((1 << 8) - 1));
	    }
	  }
	  region_map[nip][st][region] |= (1 << 31);
	  region_map[nip][st][region] |= (max << 8);
	}
	if ((options->min_count[1] == 0 || options->min_count[1] <= ((region_map[nip][st][region] >> 8) & ((1 << 8) - 1)))
	    && (options->max_count[1] == 0 || options->max_count[1] >= ((region_map[nip][st][region] >> 8) & ((1 << 8) - 1)))
	    )
	  break;
      }

    if ((map[*idx] & ((1 << region_bits) - 1)) < (uint)region_overlap && region > 0) {
      region--;

      // copy-paste from above -- SERIOUSLY
      if ((options->min_count[0] == 0 || options->min_count[0] <= (region_map[nip][st][region] & ((1 << 8) - 1)))
	  && (options->max_count[0] == 0 || options->max_count[0] >= (region_map[nip][st][region] & ((1 << 8) - 1))))
	{
	  if (options->min_count[1] != 0 && options->max_count[1] != 0)
	    break;
	  if ((region_map[nip][st][region] >> 31) == 0) {
	    // compute mp counts
	    first = MAX(0, region + re->delta_region_min[st]);
	    last = MIN(n_regions - 1, region + re->delta_region_max[st]);
	    max = 0;
	    for (k = first; k <= last; k++) {
	      if (((region_map[1-nip][1-st][k] >> 16) & ((1 << region_map_id_bits) - 1)) == region_map_id
		  && (region_map[1-nip][1-st][k] & ((1 << 8) - 1)) > max) {
		max = (region_map[1-nip][1-st][k] & ((1 << 8) - 1));
	      }
	    }
	    region_map[nip][st][region] |= (1 << 31);
	    region_map[nip][st][region] |= (max << 8);
	  }
	  if ((options->min_count[1] == 0 || options->min_count[1] <= ((region_map[nip][st][region] >> 8) & ((1 << 8) - 1)))
	      && (options->max_count[1] == 0 || options->max_count[1] >= ((region_map[nip][st][region] >> 8) & ((1 << 8) - 1)))
	      )
	    break;
	}
    }

    (*anchors_discarded)++;
    (*idx)++;
  }
}


static void
new_read_get_anchor_list_per_strand(struct read_entry * re, int st,
				    struct anchor_list_options * options)
{
  uint list_sz;
  uint offset;
  int i, sn;
  uint * idx;
  struct heap_uu h;
  struct heap_uu_elem tmp;
  int anchor_cache[re->read_len];
  int anchors_discarded = 0;
  int big_gaps = 0;

  assert(re != NULL && options != NULL);
  assert(re->anchors[st] == NULL && re->n_anchors[st] == 0);

  if (re->mapidx[st] == NULL)
    return;

  if ((st == 0 && !Fflag) || (st == 1 && !Cflag))
    return;

  // compute estimate size of anchor list
  list_sz = 0;
  for (sn = 0; sn < n_seeds; sn++) {
    for (i = 0; re->min_kmer_pos + i + seed[sn].span - 1 < re->read_len; i++) {
      offset = sn*re->max_n_kmers + i;
      if (genomemap_len[sn][re->mapidx[st][offset]] > list_cutoff)
        continue;
      list_sz += genomemap_len[sn][re->mapidx[st][offset]];
    }
  }
  stat_add(&anchor_list_init_size[omp_get_thread_num()], list_sz);

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

      if (options->use_region_counts) {
	advance_index_in_genomemap(re, st, options,
				   &idx[offset], genomemap_len[sn][re->mapidx[st][offset]],
				   genomemap[sn][re->mapidx[st][offset]],
				   &anchors_discarded);
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

    if (re->n_anchors[st] > 0 && tmp.key - re->anchors[st][re->n_anchors[st] - 1].x >= anchor_list_big_gap)
      big_gaps++;

    re->n_anchors[st]++;

    if (options->collapse) {
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

    if (options->use_region_counts) {
      advance_index_in_genomemap(re, st, options,
				 &idx[offset], genomemap_len[sn][re->mapidx[st][offset]],
				 genomemap[sn][re->mapidx[st][offset]],
				 &anchors_discarded);
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

  stat_add(&n_anchors_discarded[omp_get_thread_num()], anchors_discarded);
  stat_add(&n_big_gaps_anchor_list[omp_get_thread_num()], big_gaps);
}

static inline void
new_read_get_anchor_list(struct read_entry * re, struct anchor_list_options * options)
{
  llint before = rdtsc();

  new_read_get_anchor_list_per_strand(re, 0, options);
  new_read_get_anchor_list_per_strand(re, 1, options);

  anchor_list_ticks[omp_get_thread_num()] += rdtsc() - before;
}


static void
new_read_get_hit_list_per_strand(struct read_entry * re, int st, struct hit_list_options * options)
{
  llint goff, gstart, gend;
  int max_score, tmp_score = 0;
  int i, j, cn, max_idx;
  int w_len;
  int short_len = 0, long_len = 0;
  struct anchor a[3];

  assert(re != NULL && options != NULL);
  assert(re->hits[st] == NULL && re->n_hits[st] == 0);

  if (re->n_anchors[st] == 0)
    return;

  re->hits[st] = (struct read_hit *)xcalloc(re->n_anchors[st] * sizeof(re->hits[0][0]));

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

    if (!options->gapless) {
      // avoid single matches when n=2
      if (options->match_mode > 1 && re->anchors[st][i].weight == 1)
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

    if (options->gapless
	|| options->match_mode == 1
	|| max_score >= (int)abs_or_pct(options->threshold,
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
	re->hits[st][re->n_hits[st]].matches = (options->gapless || max_idx == i?
						re->anchors[st][i].weight :
						re->anchors[st][i].weight + re->anchors[st][max_idx].weight);
	re->hits[st][re->n_hits[st]].score_vector = -1;
	re->hits[st][re->n_hits[st]].score_full = -1;
	re->hits[st][re->n_hits[st]].score_max = (re->read_len < w_len? re->read_len : w_len) * match_score;
	re->hits[st][re->n_hits[st]].pair_min = -1;
	re->hits[st][re->n_hits[st]].pair_max = -1;
	re->hits[st][re->n_hits[st]].mapping_quality = 255;
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


static inline void
new_read_get_hit_list(struct read_entry * re, struct hit_list_options * options)
{
  llint before = rdtsc();

  new_read_get_hit_list_per_strand(re, 0, options);
  new_read_get_hit_list_per_strand(re, 1, options);

  hit_list_ticks[omp_get_thread_num()] += rdtsc() - before;

#ifdef DEBUG_HIT_LIST_CREATION
  fprintf(stderr, "Dumping hit list after creation for read:[%s]\n", re->name);
  dump_hit_list(re, 0, false, false);
  dump_hit_list(re, 1, false, false);
#endif
}


static void
new_read_pass1_per_strand(struct read_entry * re, int st, struct pass1_options * options)
{
  int i, j;

  f1_hash_tag++;
  j = -1; // last good hit

  for (i = 0; i < re->n_hits[st]; i++) {
    if (options->only_paired && re->hits[st][i].pair_min < 0) {
      continue;
    }

    // check window overlap
    if (j >= 0
	&& re->hits[st][i].cn == re->hits[st][j].cn
	&& re->hits[st][i].g_off <= re->hits[st][j].g_off + re->window_len
	- (int)abs_or_pct(options->window_overlap, re->window_len)) {
      re->hits[st][i].score_vector = 0;
      re->hits[st][i].pct_score_vector = 0;
      continue;
    }

    if (re->hits[st][i].score_vector <= 0) {
      if (re->hits[st][i].gen_st != 0)
	reverse_hit(re, &re->hits[st][i]);

      if (shrimp_mode == MODE_COLOUR_SPACE)
	re->hits[st][i].score_vector = f1_run(genome_cs_contigs[re->hits[st][i].cn], genome_len[re->hits[st][i].cn],
					      re->hits[st][i].g_off, re->hits[st][i].w_len,
					      re->read[st], re->read_len,
					      re->hits[st][i].g_off + re->hits[st][i].anchor.x, re->hits[st][i].anchor.y,
					      genome_contigs[re->hits[st][i].cn], re->initbp[st], genome_is_rna, f1_hash_tag,
					      options->gapless);
      else
	re->hits[st][i].score_vector = f1_run(genome_contigs[re->hits[st][i].cn], genome_len[re->hits[st][i].cn],
					      re->hits[st][i].g_off, re->hits[st][i].w_len,
					      re->read[st], re->read_len,
					      re->hits[st][i].g_off + re->hits[st][i].anchor.x, re->hits[st][i].anchor.y,
					      NULL, -1, genome_is_rna, f1_hash_tag,
					      options->gapless);

      re->hits[st][i].pct_score_vector = (1000 * 100 * re->hits[st][i].score_vector)/re->hits[st][i].score_max;
    }

    if (re->hits[st][i].score_vector >= (int)abs_or_pct(options->threshold, re->hits[st][i].score_max)) {
      j = i;
    }
  }
}


/*
 * Go through hit list, apply vector filter, and save top scores.
 */
static inline void
new_read_pass1(struct read_entry * re, struct pass1_options * options)
{
  new_read_pass1_per_strand(re, 0, options);
  new_read_pass1_per_strand(re, 1, options);

#ifdef DEBUG_HIT_LIST_PASS1
  fprintf(stderr, "Dumping hit list after pass1 for read:[%s]\n", re->name);
  dump_hit_list(re, 0, options->only_paired, false);
  dump_hit_list(re, 1, options->only_paired, false);
#endif
}


#define EXTHEAP_unpaired_pass1_CMP(a, b) ((a)->pass1_key < (b)->pass1_key)
DEF_EXTHEAP(struct read_hit *,unpaired_pass1)


/*
 * Go through the list adding hits to the heap.
 */
static void
new_read_get_vector_hits(struct read_entry * re, struct read_hit * * a, int * load, struct pass1_options * options)
{
  int st, i;

  assert(re != NULL && a != NULL && load != NULL && *load == 0);
  assert(pair_mode == PAIR_NONE || sam_half_paired);

  for (st = 0; st < 2; st++) {
    assert(re->n_hits[st] == 0 || re->hits[st] != NULL);

    for (i = 0; i < re->n_hits[st]; i++) {
      if (re->hits[st][i].score_vector >= (int)abs_or_pct(options->threshold, re->hits[st][i].score_max)
	  && (*load < options->num_outputs
	      || ( (IS_ABSOLUTE(options->threshold)
		    && re->hits[st][i].score_vector > a[0]->pass1_key)
		   || (!IS_ABSOLUTE(options->threshold)
		       && re->hits[st][i].pct_score_vector > a[0]->pass1_key)))) {
        //fprintf(stderr,"%d matches\n",re->hits[st][i].matches);
	re->hits[st][i].pass1_key = (IS_ABSOLUTE(options->threshold)? re->hits[st][i].score_vector : re->hits[st][i].pct_score_vector);
        //assert(re->hits[st][i].gen_st==0);
	if (*load < options->num_outputs)
	  extheap_unpaired_pass1_insert(a, load, &re->hits[st][i]);
	else
	  extheap_unpaired_pass1_replace_min(a, load, &re->hits[st][i]);
      }
    }
  }
}


/*
// don't use this for qsort directly;
static inline int
read_hit_purealign_cmp(struct read_hit * rh1, struct read_hit * rh2)
{
  if (rh1->cn < rh2->cn) {
    return -1;
  } else if (rh1->cn == rh2->cn) {
    if (rh1->st < rh2->st) {
      return -1;
    } else if (rh1->st == rh2->st) {
      return rh1->g_off - rh2->g_off;
    }
  }
  return 1;
}


// don't use this for qsort directly;
static inline int
read_hit_overlap_cmp(struct read_hit * rh1, struct read_hit * rh2)
{
  if (rh1->cn < rh2->cn) {
    return -1;
  } else if (rh1->cn == rh2->cn) {
    if (rh1->st < rh2->st) {
      return -1;
    } else if (rh1->st == rh2->st) {
      if (rh1->g_off + rh1->w_len < rh2->g_off) {
	return -1;
      } else if (rh2->g_off + rh2->w_len < rh1->g_off) {
	return 1;
      } else {
	return 0;
      }
    }
  }
  return 1;
}


// sort by: conting number; then strand; then g_off
static int
pass2_read_hit_align_cmp(void const * e1, void const * e2)
{
  struct read_hit * rh1 = *(struct read_hit * *)e1;
  struct read_hit * rh2 = *(struct read_hit * *)e2;

  return read_hit_purealign_cmp(rh1, rh2);
}


// sort by: conting number; then strand; then g_off; but return 0 if overlapping
static int
pass2_read_hit_overlap_cmp(void const * e1, void const * e2)
{
  struct read_hit * rh1 = *(struct read_hit * *)e1;
  struct read_hit * rh2 = *(struct read_hit * *)e2;

  return read_hit_overlap_cmp(rh1, rh2);
}
*/


// sort by score
static int
pass2_read_hit_score_cmp(void const * e1, void const * e2) {
  return (*(struct read_hit * *)e2)->pass2_key - (*(struct read_hit * *)e1)->pass2_key;
}


static inline int
pass2_read_hit_sfrp_gen_start_cmp_base(struct read_hit const * rh1, struct read_hit const * rh2) {
  if (rh1->cn != rh2->cn)
    return rh1->cn - rh2->cn;

  if (rh1->gen_st != rh2->gen_st)
    return rh1->gen_st - rh2->gen_st;

  return rh1->sfrp->genome_start - rh2->sfrp->genome_start;
}

static inline int
pass2_read_hit_sfrp_gen_end_cmp_base(struct read_hit const * rh1, struct read_hit const * rh2) {
  if (rh1->cn != rh2->cn)
    return rh1->cn - rh2->cn;

  if (rh1->gen_st != rh2->gen_st)
    return rh1->gen_st - rh2->gen_st;

  return (- rh1->sfrp->genome_start - rh1->sfrp->rmapped + rh1->sfrp->deletions - rh1->sfrp->insertions)
    - (- rh2->sfrp->genome_start - rh2->sfrp->rmapped + rh2->sfrp->deletions - rh2->sfrp->insertions);
}

static int
pass2_read_hit_sfrp_gen_start_cmp(void const * e1, void const * e2) {
  return pass2_read_hit_sfrp_gen_start_cmp_base(*(struct read_hit * *)e1, *(struct read_hit * *)e2);
}

static int
pass2_read_hit_sfrp_gen_end_cmp(void const * e1, void const * e2) {
  return pass2_read_hit_sfrp_gen_end_cmp_base(*(struct read_hit * *)e1, *(struct read_hit * *)e2);
}


// remove duplicate hits
static void
read_remove_duplicate_hits(struct read_hit * * hits_pass2, int * n_hits_pass2)
{
  int i, j, k, max, max_idx;
  llint before = rdtsc();

  /*
  qsort(hits_pass2, *n_hits_pass2, sizeof(hits_pass2[0]), pass2_read_hit_align_cmp);
  i = 0;
  k = 0;
  while (i < *n_hits_pass2) {
    max = hits_pass2[i]->pass2_key;
    max_idx = i;
    j = i + 1;
    while (j < *n_hits_pass2 && !pass2_read_hit_overlap_cmp(&hits_pass2[i], &hits_pass2[j])) {
      if (hits_pass2[j]->pass2_key > max) {
	max = hits_pass2[j]->pass2_key;
	max_idx = j;
      }
      j++;
    }
    if (max_idx != k) {
      hits_pass2[k] = hits_pass2[max_idx];
    }
    k++;
    i = j;
  }
  return k;
  */

  qsort(hits_pass2, *n_hits_pass2, sizeof(hits_pass2[0]), pass2_read_hit_sfrp_gen_start_cmp);
  i = 0;
  k = 0;
  while (i < *n_hits_pass2) {
    max = hits_pass2[i]->pass2_key;
    max_idx = i;
    j = i + 1;
    while (j < *n_hits_pass2 && !pass2_read_hit_sfrp_gen_start_cmp(&hits_pass2[i], &hits_pass2[j])) {
      if (hits_pass2[j]->pass2_key > max) {
	max = hits_pass2[j]->pass2_key;
	max_idx = j;
      }
      j++;
    }
    if (max_idx != k) {
      hits_pass2[k] = hits_pass2[max_idx];
    }
    k++;
    i = j;
  }
#pragma omp atomic
  total_dup_single_matches += (*n_hits_pass2) - k;

  *n_hits_pass2 = k;

  qsort(hits_pass2, *n_hits_pass2, sizeof(hits_pass2[0]), pass2_read_hit_sfrp_gen_end_cmp);
  i = 0;
  k = 0;
  while (i < *n_hits_pass2) {
    max = hits_pass2[i]->pass2_key;
    max_idx = i;
    j = i + 1;
    while (j < *n_hits_pass2 && !pass2_read_hit_sfrp_gen_end_cmp(&hits_pass2[i], &hits_pass2[j])) {
      if (hits_pass2[j]->pass2_key > max) {
	max = hits_pass2[j]->pass2_key;
	max_idx = j;
      }
      j++;
    }
    if (max_idx != k) {
      hits_pass2[k] = hits_pass2[max_idx];
    }
    k++;
    i = j;
  }
#pragma omp atomic
  total_dup_single_matches += (*n_hits_pass2) - k;

  *n_hits_pass2 = k;

  duplicate_removal_ticks[omp_get_thread_num()] += rdtsc() - before;
}


static void
hit_run_post_sw(struct read_entry * re, struct read_hit * rh)
{
  if (shrimp_mode == MODE_COLOUR_SPACE) {
    post_sw(re->read[0], re->initbp[0], re->qual, rh->sfrp);
  } else { // LS: cheat; reuse SW score to get posterior
    rh->sfrp->posterior = pow(2.0, ((double)rh->sfrp->score - (double)rh->sfrp->rmapped * (2.0 * score_alpha + score_beta))/score_alpha);
  }

  rh->sfrp->posterior_score = (int)rint(score_alpha * log(rh->sfrp->posterior) / log(2.0) + (double)rh->sfrp->rmapped * (2.0 * score_alpha + score_beta));
  rh->sfrp->pct_posterior_score = (1000 * 100 * rh->sfrp->posterior_score)/rh->score_max;
  rh->score_full = rh->sfrp->posterior_score;
  rh->pct_score_full = rh->sfrp->pct_posterior_score;
}


/*
 * Do a final pass for given read.
 */
static bool
new_read_pass2(struct read_entry * re,
	       struct read_hit * * hits_pass1, int * n_hits_pass1,
	       struct read_hit * * hits_pass2, int * n_hits_pass2,
	       struct pass2_options * options)
{
  int i, cnt;
  assert(re != NULL && hits_pass1 != NULL && hits_pass2 != NULL);

  /* compute full alignment scores */
  for (i = 0; i < *n_hits_pass1; i++) {
    struct read_hit * rh = hits_pass1[i];
    if (rh->score_full < 0 || rh->sfrp == NULL) {
      hit_run_full_sw(re, rh, (int)abs_or_pct(options->threshold, rh->score_max));
      if (compute_mapping_qualities && rh->score_full > 0) {
	hit_run_post_sw(re, rh);
	/*
	fprintf(stderr, "read:%s\tSW-prob:%g\tposterior:%g\n", re->name,
		pow(2.0, ((double)rh->sfrp->score - (double)rh->sfrp->rmapped * (2.0 * score_alpha + score_beta))/score_alpha)
		* pow(1.0 - pr_xover, rh->sfrp->rmapped - rh->sfrp->crossovers),
		rh->sfrp->posterior);
	*/
      }

      rh->pass2_key = (IS_ABSOLUTE(options->threshold)? rh->score_full : (int)rh->pct_score_full);
    }

    if (rh->score_full >= abs_or_pct(options->threshold, rh->score_max)) {
      hits_pass2[*n_hits_pass2] = rh;
      (*n_hits_pass2)++;
    } 
  }

#ifdef DEBUG_HIT_LIST_PASS2
  fprintf(stderr, "Dumping hit list after pass2 (before duplicates removal and sorting) for read:[%s]\n", re->name);
  for (i = 0; i < *n_hits_pass1; i++) {
    dump_hit(hits_pass1[i]);
  }
#endif

  // remove duplicate hits
  read_remove_duplicate_hits(hits_pass2, n_hits_pass2);

  // sort by non-increasing score
  qsort(hits_pass2, *n_hits_pass2, sizeof(hits_pass2[0]), pass2_read_hit_score_cmp);

  if (compute_mapping_qualities && *n_hits_pass2 > 0) {
    // compute Z
    re->mq_denominator = hits_pass2[0]->sfrp->posterior;
    for (i = 1; i < *n_hits_pass2 && hits_pass2[i]->score_full > hits_pass2[0]->score_full - score_difference_mq_cutoff; i++) {
      re->mq_denominator += hits_pass2[i]->sfrp->posterior;
    }

    // compute mapping qualities
    for (i = 0; i < *n_hits_pass2 && hits_pass2[i]->score_full > hits_pass2[0]->score_full - score_difference_mq_cutoff; i++) {
      hits_pass2[i]->mapping_quality = qv_from_pr_corr(hits_pass2[i]->sfrp->posterior / re->mq_denominator);
      if (hits_pass2[i]->mapping_quality >= 10) {
	#pragma omp atomic
	total_reads_matched_conf++;
      }
    }
    for ( ; i < *n_hits_pass2; i++) {
      hits_pass2[i]->mapping_quality = 0;
    }

  }

  // trim excess mappings
  if (*n_hits_pass2 > options->num_outputs)
    *n_hits_pass2 = options->num_outputs;

  // if strata is set, keep only top scoring hits
  if (options->strata && *n_hits_pass2 > 0) {
    for (i = 1; i < *n_hits_pass2 && hits_pass2[0]->score_full == hits_pass2[i]->score_full; i++);
    *n_hits_pass2 = i;
  }

  // drop reads with too many mappings
  if (*n_hits_pass2 > 0) {
    if (max_alignments == 0 || *n_hits_pass2 <= max_alignments) {
#pragma omp atomic
      total_reads_matched++;
    } else {
#pragma omp atomic
      total_reads_dropped++;
      *n_hits_pass2 = 0;
    }
  }

  // update counts
  re->final_matches += *n_hits_pass2;
#pragma omp atomic
  total_single_matches += re->final_matches;

  // check stop condition
  if (options->stop_count == 0)
    return true;

  for (i = 0, cnt = 0; i < *n_hits_pass2; i++) {
    if (hits_pass2[i]->score_full >= (int)abs_or_pct(options->stop_threshold, hits_pass2[i]->score_max)) {
      cnt++;
    }
  }

  return cnt >= options->stop_count;
}


void
new_handle_read(struct read_entry * re, struct read_mapping_options_t * options, int n_options)
{
  bool done;
  int option_index = 0;
  int i;
  struct read_hit * * hits_pass1;
  struct read_hit * * hits_pass2;
  int n_hits_pass1;
  int n_hits_pass2;

  uint64_t before = rdtsc();

  if (re->mapidx[0] == NULL) {
    read_get_mapidxs(re);
  }

  do {

    if (options[option_index].regions.recompute) {
      region_map_id++;
      region_map_id &= ((1 << region_map_id_bits) - 1);
      read_get_region_counts(re, 0, &options[option_index].regions);
      read_get_region_counts(re, 1, &options[option_index].regions);
    }

    if (options[option_index].anchor_list.recompute) {
      read_free_anchor_list(re);
      new_read_get_anchor_list(re, &options[option_index].anchor_list);
    }

    if (options[option_index].hit_list.recompute) {
      read_free_hit_list(re);
      new_read_get_hit_list(re, &options[option_index].hit_list);
    }

    if (options[option_index].pass1.recompute) {
      new_read_pass1(re, &options[option_index].pass1);
    }

    hits_pass1 = (struct read_hit * *)xmalloc(options[option_index].pass1.num_outputs * sizeof(hits_pass1[0]));
    n_hits_pass1 = 0;
    new_read_get_vector_hits(re, hits_pass1, &n_hits_pass1, &options[option_index].pass1);

    hits_pass2 = (struct read_hit * *)xmalloc(options[option_index].pass1.num_outputs * sizeof(hits_pass2[0]));
    n_hits_pass2 = 0;
    done = new_read_pass2(re, hits_pass1, &n_hits_pass1, hits_pass2, &n_hits_pass2, &options[option_index].pass2);

    new_read_output(re, hits_pass2, &n_hits_pass2);

    for (i = 0; i < n_hits_pass1; i++) {
      free_sfrp(&hits_pass1[i]->sfrp);
    }

    free(hits_pass1);
    free(hits_pass2);

  } while (!done && ++option_index < n_options);

  //if (options_index >= n_options) {
    // this read fell through all the option sets
  //}

  scan_ticks[omp_get_thread_num()] += rdtsc() - before;
}


#define EXTHEAP_paired_pass1_CMP(a, b) ((a).key < (b).key)
DEF_EXTHEAP(struct read_hit_pair, paired_pass1)

/*
 * Go through the hit lists, constructing paired hits.
 */
static void
new_readpair_get_vector_hits(struct read_entry * re1, struct read_entry * re2,
			     struct read_hit_pair * a, int * load,
			     struct pairing_options * options)
{
  int st1, st2, i, j;
  read_hit_pair tmp;

  assert(re1 != NULL && re2 != NULL && a != NULL);

  for (st1 = 0; st1 < 2; st1++) {
    st2 = 1 - st1; // opposite strand

    for (i = 0; i < re1->n_hits[st1]; i++) {
      if (re1->hits[st1][i].pair_min < 0)
	continue;

      for (j = re1->hits[st1][i].pair_min; j <= re1->hits[st1][i].pair_max; j++) {
	if (re1->hits[st1][i].matches + re2->hits[st2][j].matches < options->min_num_matches)
	  continue;

	tmp.score = re1->hits[st1][i].score_vector + re2->hits[st2][j].score_vector;
	tmp.score_max = re1->hits[st1][i].score_max + re2->hits[st2][j].score_max;
	tmp.pct_score = (1000 * 100 * tmp.score)/tmp.score_max;
	tmp.key = (IS_ABSOLUTE(options->pass1_threshold)? tmp.score : tmp.pct_score);

	if (tmp.score >= (int)abs_or_pct(options->pass1_threshold, tmp.score_max)
	    && (*load < options->pass1_num_outputs || tmp.key > a[0].key)) {
	  tmp.rh[0] = &re1->hits[st1][i];
	  tmp.rh[1] = &re2->hits[st2][j];
	  //TODO HISTORGRAM IS OFF! NOT SAM
	  tmp.insert_size = (int)(st1 == 0?
				  re2->hits[st2][j].g_off - (re1->hits[st1][i].g_off + re1->hits[st1][i].w_len) :
				  re1->hits[st1][i].g_off - (re2->hits[st2][j].g_off + re2->hits[st2][j].w_len));

	  if (*load < options->pass1_num_outputs)
	    extheap_paired_pass1_insert(a, load, tmp);
	  else
	    extheap_paired_pass1_replace_min(a, load, tmp);
	}
      }
    }
  }
}


/*
// sort by: contig number; then strand; then g_off of hit[0]
static int
pass2_readpair_hit0_align_cmp(void const * e1, void const * e2) {
  struct read_hit_pair * rhp1 = (struct read_hit_pair *)e1;
  struct read_hit_pair * rhp2 = (struct read_hit_pair *)e2;

  return read_hit_purealign_cmp(rhp1->rh[0], rhp2->rh[0]);
}


// sort by: contig number; then strand; then g_off of hit[1]
static int
pass2_readpair_hit1_align_cmp(void const * e1, void const * e2) {
  struct read_hit_pair * rhp1 = (struct read_hit_pair *)e1;
  struct read_hit_pair * rhp2 = (struct read_hit_pair *)e2;

  return read_hit_purealign_cmp(rhp1->rh[1], rhp2->rh[1]);
}


// sort by: contig number; then strand; then g_off of hit[0]; but return 0 if hit[0] overlapping
static int
pass2_readpair_hit0_overlap_cmp(void const * e1, void const * e2) {
  struct read_hit_pair * rhp1 = (struct read_hit_pair *)e1;
  struct read_hit_pair * rhp2 = (struct read_hit_pair *)e2;

  return read_hit_overlap_cmp(rhp1->rh[0], rhp2->rh[0]);
}


// sort by: contig number; then strand; then g_off of hit[1]; but return 0 if hit[1] overlapping
static int
pass2_readpair_hit1_overlap_cmp(void const * e1, void const * e2) {
  struct read_hit_pair * rhp1 = (struct read_hit_pair *)e1;
  struct read_hit_pair * rhp2 = (struct read_hit_pair *)e2;

  return read_hit_overlap_cmp(rhp1->rh[1], rhp2->rh[1]);
}
*/


// sort by score
static int
pass2_read_hit_pair_score_cmp(void const * e1, void const * e2)
{
  return ((struct read_hit_pair *)e2)->key - ((struct read_hit_pair *)e1)->key;
}

static int
pass2_readpair_hit0_sfrp_gen_start_cmp(void const * e1, void const * e2)
{
  return pass2_read_hit_sfrp_gen_start_cmp_base(((read_hit_pair *)e1)->rh[0], ((read_hit_pair *)e2)->rh[0]);
}

static int
pass2_readpair_hit1_sfrp_gen_start_cmp(void const * e1, void const * e2)
{
  return pass2_read_hit_sfrp_gen_start_cmp_base(((read_hit_pair *)e1)->rh[1], ((read_hit_pair *)e2)->rh[1]);
}

static int
pass2_readpair_hit0_sfrp_gen_end_cmp(void const * e1, void const * e2)
{
  return pass2_read_hit_sfrp_gen_end_cmp_base(((read_hit_pair *)e1)->rh[0], ((read_hit_pair *)e2)->rh[0]);
}

static int
pass2_readpair_hit1_sfrp_gen_end_cmp(void const * e1, void const * e2)
{
  return pass2_read_hit_sfrp_gen_end_cmp_base(((read_hit_pair *)e1)->rh[1], ((read_hit_pair *)e2)->rh[1]);
}

static int
pass2_readpair_pointer_cmp(void const * e1, void const * e2)
{
  struct read_hit_pair * rhpp1 = (struct read_hit_pair *)e1;
  struct read_hit_pair * rhpp2 = (struct read_hit_pair *)e2;

  if (rhpp1->rh[0] < rhpp2->rh[0])
    return -1;
  else if (rhpp1->rh[0] > rhpp2->rh[0])
    return 1;
  else { // equal rh[0]
    if (rhpp1->rh[1] < rhpp2->rh[1])
      return -1;
    else if (rhpp1->rh[1] > rhpp2->rh[1])
      return 1;
  }
  return 0;
}


static inline void
readpair_compute_paired_hit(struct read_hit * rh1, struct read_hit * rh2, bool threshold_is_absolute,
			    struct read_hit_pair * dest)
{ 
  dest->rh[0] = rh1;
  dest->rh[1] = rh2;
  dest->score_max = rh1->score_max + rh2->score_max;
  dest->score = rh1->score_full + rh2->score_full;
  dest->pct_score = (1000 * 100 * dest->score)/dest->score_max;
  dest->key = threshold_is_absolute? dest->score : dest->pct_score;
  dest->insert_size = abs(get_isize(rh1, rh2));
  //tmp.isize_score=expected_isize==-1 ? 0 : abs(tmp.isize-expected_isize);
}


static void
readpair_push_dominant_single_hits(struct read_hit_pair * hits_pass2, int * n_hits_pass2, bool threshold_is_absolute,
				   int num_in_pair, int (*cmp)(void const *, void const *))
{
  int i, j, k, max, max_idx;

  qsort(hits_pass2, *n_hits_pass2, sizeof(hits_pass2[0]), cmp);
  i = 0;
  while (i < *n_hits_pass2) {
    max = hits_pass2[i].rh[num_in_pair]->score_full;
    max_idx = i;
    j = i + 1;
    while (j < *n_hits_pass2 && !cmp((void *)&hits_pass2[i], (void *)&hits_pass2[j])) {
      if (hits_pass2[j].rh[num_in_pair]->score_full > max) {
	max = hits_pass2[j].rh[num_in_pair]->score_full;
	max_idx = j;
      }
      j++;
    }
    for (k = i; k < j; k++) {
      if (k != max_idx) {
	hits_pass2[k].rh[num_in_pair] = hits_pass2[max_idx].rh[num_in_pair];
	readpair_compute_paired_hit(hits_pass2[k].rh[0], hits_pass2[k].rh[1], threshold_is_absolute, &hits_pass2[k]);
      }
    }
    i = j;
  }
}


// remove duplicate hits
static void
readpair_remove_duplicate_hits(struct read_hit_pair * hits_pass2, int * n_hits_pass2, bool threshold_is_absolute)
{
  /*
  int i, j, k, l, max, max_idx;

  qsort(hits_pass2, *n_hits_pass2, sizeof(hits_pass2[0]), pass2_readpair_hit0_align_cmp);
  i = 0;
  k = 0;
  while (i < *n_hits_pass2) {
    j = i + 1;
    while (j < *n_hits_pass2 && !pass2_readpair_hit0_overlap_cmp(&hits_pass2[j-1], &hits_pass2[j])) {
      j++;
    }
    if (j > i + 1) {
      qsort(&hits_pass2[i], j - i, sizeof(hits_pass2[0]), pass2_readpair_hit1_align_cmp);
    }

    while (i < j) {
      max = hits_pass2[i].key;
      max_idx = i;
      l = i + 1;
      while (l < j && !pass2_readpair_hit1_overlap_cmp(&hits_pass2[max_idx], &hits_pass2[l])) {
	if (hits_pass2[l].key > max) {
	  max = hits_pass2[l].key;
	  max_idx = l;
	}
	l++;
      }
      if (max_idx != k) {
	hits_pass2[k] = hits_pass2[max_idx];
      }
      k++;
      i = l;
    }
  }
  return k;
  */

  int tmp;
  llint before = rdtsc();

  readpair_push_dominant_single_hits(hits_pass2, n_hits_pass2, threshold_is_absolute, 0, pass2_readpair_hit0_sfrp_gen_start_cmp);
  readpair_push_dominant_single_hits(hits_pass2, n_hits_pass2, threshold_is_absolute, 0, pass2_readpair_hit0_sfrp_gen_end_cmp);
  readpair_push_dominant_single_hits(hits_pass2, n_hits_pass2, threshold_is_absolute, 1, pass2_readpair_hit1_sfrp_gen_start_cmp);
  readpair_push_dominant_single_hits(hits_pass2, n_hits_pass2, threshold_is_absolute, 1, pass2_readpair_hit1_sfrp_gen_end_cmp);

  qsort(hits_pass2, *n_hits_pass2, sizeof(hits_pass2[0]), pass2_readpair_pointer_cmp);
  tmp = removedups(hits_pass2, *n_hits_pass2, sizeof(hits_pass2[0]), pass2_readpair_pointer_cmp);

#pragma omp atomic
  total_dup_paired_matches += (*n_hits_pass2) - tmp;

  *n_hits_pass2 = tmp;

  duplicate_removal_ticks[omp_get_thread_num()] += rdtsc() - before;
}


/*
 * Do a final pass for given read.
 */
static bool
new_readpair_pass2(struct read_entry * re1, struct read_entry * re2,
		   struct read_hit_pair * hits_pass1, int * n_hits_pass1,
		   struct read_hit_pair * hits_pass2, int * n_hits_pass2,
		   struct pairing_options * options)
{
  int i, j, cnt;

  /* compute full alignment scores */
  for (i = 0; i < *n_hits_pass1; i++) {
    for (j = 0; j < 2; j++) {
      struct read_hit * rh = hits_pass1[i].rh[j];
      struct read_entry * re = (j == 0? re1 : re2);

      if (rh->score_full < 0 || rh->sfrp == NULL) {
	hit_run_full_sw(re, rh, (int)abs_or_pct(options->pass2_threshold, hits_pass1[i].score_max)/3);
      }
    }

    //hitpair_run_post_sw(re1, re2, hits_pass1[i].rh[0], hits_pass1[i].rh[1]);

    if (hits_pass1[i].rh[0]->score_full == 0 || hits_pass1[i].rh[1]->score_full == 0) {
      continue;
    }

    if (hits_pass1[i].rh[0]->score_full + hits_pass1[i].rh[1]->score_full
	>= (int)abs_or_pct(options->pass2_threshold, hits_pass1[i].score_max)) {
      readpair_compute_paired_hit(hits_pass1[i].rh[0], hits_pass1[i].rh[1], IS_ABSOLUTE(options->pass2_threshold),
				  &hits_pass2[*n_hits_pass2]);
      (*n_hits_pass2)++;
    }
  }

#ifdef DEBUG_HIT_LIST_PASS2
  fprintf(stderr, "Dumping paired hits after pass2 (before duplicates removal and sorting) for reads:[%s,%s]\n",
	  re1->name, re2->name);
  for (i = 0; i < *n_hits_pass1; i++) {
    dump_hit(hits_pass1[i].rh[0]);
    dump_hit(hits_pass1[i].rh[1]);
  }
#endif

  // remove duplicates
  readpair_remove_duplicate_hits(hits_pass2, n_hits_pass2, IS_ABSOLUTE(options->pass2_threshold));

  // sort by score
  qsort(hits_pass2, *n_hits_pass2, sizeof(hits_pass2[0]), pass2_read_hit_pair_score_cmp);

  // trim excess mappings
  if (*n_hits_pass2 > options->pass2_num_outputs)
    *n_hits_pass2 = options->pass2_num_outputs;

  // if strata is set, keep only top scoring hits
  if (options->strata && *n_hits_pass2 > 0) {
    for (i = 1; i < *n_hits_pass2 && hits_pass2[0].score == hits_pass2[i].score; i++);
    *n_hits_pass2 = i;
  }

  // drop pairs with too many mappings
  if (*n_hits_pass2 > 0) {
    if (max_alignments == 0 || *n_hits_pass2 <= max_alignments) {
#pragma omp atomic
      total_pairs_matched++;
    } else {
#pragma omp atomic
      total_pairs_dropped++;
      *n_hits_pass2 = 0;
    }
  }

  // update counts
  re1->final_matches += *n_hits_pass2;
#pragma omp atomic
  total_paired_matches += re1->final_matches;

  // check stop condition
  if (options->stop_count == 0)
    return true;

  for (i = 0, cnt = 0; i < *n_hits_pass2; i++) {
    if (hits_pass2[i].score >= (int)abs_or_pct(options->stop_threshold, hits_pass2[i].score_max)) {
      cnt++;
    }
  }

  return cnt >= options->stop_count;
}


static void
readpair_compute_mp_ranges(struct read_entry * re1, struct read_entry * re2,
			   struct pairing_options * options)
{
  switch (pair_mode) {
  case PAIR_OPP_IN:
    re1->delta_g_off_min[0] =   options->min_insert_size						- re2->window_len;
    re1->delta_g_off_max[0] =   options->max_insert_size + (re1->window_len - re1->read_len)	- re2->read_len;
    re1->delta_g_off_min[1] = - options->max_insert_size + re1->read_len				+ (re2->read_len - re2->window_len);
    re1->delta_g_off_max[1] = - options->min_insert_size + re1->window_len;

    re2->delta_g_off_min[0] = - re1->delta_g_off_max[1];
    re2->delta_g_off_max[0] = - re1->delta_g_off_min[1];
    re2->delta_g_off_min[1] = - re1->delta_g_off_max[0];
    re2->delta_g_off_max[1] = - re1->delta_g_off_min[0];
    break;

  case PAIR_OPP_OUT:
    re1->delta_g_off_min[0] = - options->max_insert_size						- re2->window_len;
    re1->delta_g_off_max[0] = - options->min_insert_size + (re1->window_len - re1->read_len)	- re2->read_len;
    re1->delta_g_off_min[1] =   options->min_insert_size + re1->read_len				+ (re2->read_len - re2->window_len);
    re1->delta_g_off_max[1] =   options->max_insert_size + re1->window_len;

    re2->delta_g_off_min[0] = - re1->delta_g_off_max[1];
    re2->delta_g_off_max[0] = - re1->delta_g_off_min[1];
    re2->delta_g_off_min[1] = - re1->delta_g_off_max[0];
    re2->delta_g_off_max[1] = - re1->delta_g_off_min[0];
    break;

  case PAIR_COL_FW:
    re1->delta_g_off_min[0] =   options->min_insert_size						+ (re2->read_len - re2->window_len);
    re1->delta_g_off_max[0] =   options->max_insert_size + (re1->window_len - re1->read_len);
    re1->delta_g_off_min[1] = - options->max_insert_size + re1->read_len				- re2->window_len;
    re1->delta_g_off_max[1] = - options->min_insert_size + re1->window_len			- re2->read_len;

    re2->delta_g_off_min[0] = - re1->delta_g_off_max[0];
    re2->delta_g_off_max[0] = - re1->delta_g_off_min[0];
    re2->delta_g_off_min[1] = - re1->delta_g_off_max[1];
    re2->delta_g_off_max[1] = - re1->delta_g_off_min[1];
    break;

  case PAIR_COL_BW:
    re1->delta_g_off_min[0] = - options->max_insert_size						+ (re2->read_len - re2->window_len);
    re1->delta_g_off_max[0] = - options->min_insert_size + (re1->window_len - re1->read_len);
    re1->delta_g_off_min[1] =   options->min_insert_size + re1->read_len				- re2->window_len;
    re1->delta_g_off_max[1] =   options->max_insert_size + re1->window_len			- re2->read_len;

    re2->delta_g_off_min[0] = - re1->delta_g_off_max[0];
    re2->delta_g_off_max[0] = - re1->delta_g_off_min[0];
    re2->delta_g_off_min[1] = - re1->delta_g_off_max[1];
    re2->delta_g_off_max[1] = - re1->delta_g_off_min[1];
    break;

  default:
    assert(0);
    break;
  }

  re1->delta_region_min[0] = re1->delta_g_off_min[0] >= 0? re1->delta_g_off_min[0]/(1 << region_bits) : - 1 - (- re1->delta_g_off_min[0] - 1)/(1 << region_bits);
  re1->delta_region_max[0] = re1->delta_g_off_max[0] > 0? 1 + (re1->delta_g_off_max[0] - 1)/(1 << region_bits) : - (- re1->delta_g_off_max[0]/(1 << region_bits));
  re1->delta_region_min[1] = re1->delta_g_off_min[1] >= 0? re1->delta_g_off_min[1]/(1 << region_bits) : - 1 - (- re1->delta_g_off_min[1] - 1)/(1 << region_bits);
  re1->delta_region_max[1] = re1->delta_g_off_max[1] > 0? 1 + (re1->delta_g_off_max[1] - 1)/(1 << region_bits) : - (- re1->delta_g_off_max[1]/(1 << region_bits));

  re2->delta_region_min[0] = re2->delta_g_off_min[0] >= 0? re2->delta_g_off_min[0]/(1 << region_bits) : - 1 - (- re2->delta_g_off_min[0] - 1)/(1 << region_bits);
  re2->delta_region_max[0] = re2->delta_g_off_max[0] > 0? 1 + (re2->delta_g_off_max[0] - 1)/(1 << region_bits) : - (- re2->delta_g_off_max[0]/(1 << region_bits));
  re2->delta_region_min[1] = re2->delta_g_off_min[1] >= 0? re2->delta_g_off_min[1]/(1 << region_bits) : - 1 - (- re2->delta_g_off_min[1] - 1)/(1 << region_bits);
  re2->delta_region_max[1] = re2->delta_g_off_max[1] > 0? 1 + (re2->delta_g_off_max[1] - 1)/(1 << region_bits) : - (- re2->delta_g_off_max[1]/(1 << region_bits));

#ifdef DEBUG_DUMP_MP_RANGES
  fprintf(stderr, "mp_ranges read[%s]: goff_min[0]:%d goff_max[0]:%d goff_min[1]:%d goff_max[1]:%d reg_min[0]:%d reg_max[0]:%d reg_min[1]:%d reg_max[1]:%d\n",
	  re1->name,
	  re1->delta_g_off_min[0], re1->delta_g_off_max[0], re1->delta_g_off_min[1], re1->delta_g_off_max[1],
	  re1->delta_region_min[0], re1->delta_region_max[0], re1->delta_region_min[1], re1->delta_region_max[1]);
  fprintf(stderr, "mp_ranges read[%s]: goff_min[0]:%d goff_max[0]:%d goff_min[1]:%d goff_max[1]:%d reg_min[0]:%d reg_max[0]:%d reg_min[1]:%d reg_max[1]:%d\n",
	  re2->name,
	  re2->delta_g_off_min[0], re2->delta_g_off_max[0], re2->delta_g_off_min[1], re2->delta_g_off_max[1],
	  re2->delta_region_min[0], re2->delta_region_max[0], re2->delta_region_min[1], re2->delta_region_max[1]);
#endif
}


void
new_handle_readpair(struct read_entry * re1, struct read_entry * re2,
		    struct readpair_mapping_options_t * options, int n_options)
{
  bool done;
  int option_index = 0;
  int i;
  struct read_hit_pair * hits_pass1;
  struct read_hit_pair * hits_pass2;
  int n_hits_pass1;
  int n_hits_pass2;

  uint64_t before = rdtsc();

  read_get_mapidxs(re1);
  read_get_mapidxs(re2);

  do {
    readpair_compute_mp_ranges(re1, re2, &options[option_index].pairing);

    if (options[option_index].read[0].regions.recompute || options[option_index].read[1].regions.recompute) {
      region_map_id++;
      region_map_id &= ((1 << region_map_id_bits) - 1);
      read_get_region_counts(re1, 0, &options[option_index].read[0].regions);
      read_get_region_counts(re1, 1, &options[option_index].read[0].regions);
      read_get_region_counts(re2, 0, &options[option_index].read[1].regions);
      read_get_region_counts(re2, 1, &options[option_index].read[1].regions);
    }

    if (options[option_index].read[0].anchor_list.recompute) {
      read_free_anchor_list(re1);
      new_read_get_anchor_list(re1, &options[option_index].read[0].anchor_list);
    }
    if (options[option_index].read[1].anchor_list.recompute) {
      read_free_anchor_list(re2);
      new_read_get_anchor_list(re2, &options[option_index].read[1].anchor_list);
    }

    if (options[option_index].read[0].hit_list.recompute) {
      read_free_hit_list(re1);
      new_read_get_hit_list(re1, &options[option_index].read[0].hit_list);
    }
    if (options[option_index].read[1].hit_list.recompute) {
      read_free_hit_list(re2);
      new_read_get_hit_list(re2, &options[option_index].read[1].hit_list);
    }

    readpair_pair_up_hits(re1, re2);

    if (options[option_index].read[0].pass1.recompute) {
      new_read_pass1(re1, &options[option_index].read[0].pass1);
    }
    if (options[option_index].read[1].pass1.recompute) {
      new_read_pass1(re2, &options[option_index].read[1].pass1);
    }

    hits_pass1 = (struct read_hit_pair *)xmalloc(options[option_index].pairing.pass1_num_outputs * sizeof(hits_pass1[0]));
    n_hits_pass1 = 0;
    new_readpair_get_vector_hits(re1, re2, hits_pass1, &n_hits_pass1, &options[option_index].pairing);

    hits_pass2 = (struct read_hit_pair *)xmalloc(options[option_index].pairing.pass1_num_outputs * sizeof(hits_pass2[0]));
    n_hits_pass2 = 0;
    done = new_readpair_pass2(re1, re2, hits_pass1, &n_hits_pass1, hits_pass2, &n_hits_pass2, &options[option_index].pairing);

    new_readpair_output(re1, re2, hits_pass2, &n_hits_pass2);

    for (i = 0; i < n_hits_pass1; i++) {
      free_sfrp(&hits_pass1[i].rh[0]->sfrp);
      free_sfrp(&hits_pass1[i].rh[1]->sfrp);
    }

    free(hits_pass1);
    free(hits_pass2);

  } while (!done && ++option_index < n_options);

  scan_ticks[omp_get_thread_num()] += rdtsc() - before;

  if (option_index >= n_options && sam_half_paired) {
    // this read pair fell through all the option sets; try unpaired mapping
    new_handle_read(re1, unpaired_mapping_options[0], n_unpaired_mapping_options[0]);
    new_handle_read(re2, unpaired_mapping_options[1], n_unpaired_mapping_options[1]);
  }
}
