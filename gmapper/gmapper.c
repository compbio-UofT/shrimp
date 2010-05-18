#include <assert.h>
#include <ctype.h>
#include <errno.h>
#include <fcntl.h>
#include <inttypes.h>
#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <unistd.h>
#include <zlib.h>
#include <omp.h>	// OMP multithreading
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/time.h>

#include "../common/hash.h"
#include "../common/fasta.h"
#include "../common/util.h"
#include "../gmapper/gmapper.h"
#include "../common/version.h"
#include "../common/sw-full-common.h"
#include "../common/sw-full-cs.h"
#include "../common/sw-full-ls.h"
#include "../common/sw-vector.h"
#include "../common/output.h"
#include "../common/input.h"

#include "../common/heap.h"
DEF_HEAP(uint32_t,uint,uu)
DEF_HEAP(uint, struct read_hit_holder, unpaired)
DEF_HEAP(uint, struct read_hit_pair_holder, paired)

/* Parameters */
static double	window_len		= DEF_WINDOW_LEN;
static double	window_overlap		= DEF_WINDOW_OVERLAP;
static int	num_matches		= DEF_NUM_MATCHES;
static int	num_outputs		= DEF_NUM_OUTPUTS;
static int	num_tmp_outputs		= 20 + num_outputs;
static int	anchor_width		= DEF_ANCHOR_WIDTH;
static uint32_t	list_cutoff		= DEF_LIST_CUTOFF;
static bool	gapless_sw		= DEF_GAPLESS_SW;
static bool	hash_filter_calls	= DEF_HASH_FILTER_CALLS;

/* contains inlined calls; uses gapless_sw and hash_filter_calls vars */
#include "../common/f1-wrapper.h"

/* Scores */
static int	match_score		= DEF_MATCH_VALUE;
static int	mismatch_score		= DEF_MISMATCH_VALUE;
static int	a_gap_open_score	= DEF_A_GAP_OPEN;
static int	a_gap_extend_score	= DEF_A_GAP_EXTEND;
static int	b_gap_open_score	= DEF_B_GAP_OPEN;
static int	b_gap_extend_score	= DEF_B_GAP_EXTEND;
static int	crossover_score		= DEF_XOVER_PENALTY;

/* Score Thresholds */
static double	window_gen_threshold	= DEF_WINDOW_GEN_THRESHOLD;
static double	sw_vect_threshold	= DEF_SW_VECT_THRESHOLD;
static double	sw_full_threshold	= DEF_SW_FULL_THRESHOLD;

/* Flags */
static bool Cflag = false;			/* do complement only */
static bool Fflag = false;			/* do positive (forward) only */
static bool Hflag = false;			/* use hash table, not lookup */
static bool Pflag = false;			/* pretty print results */
static bool Rflag = false;			/* add read sequence to output*/
static bool Tflag = false;			/* reverse sw full tie breaks */
static bool Dflag = false;			/* print statistics for each thread */
static bool Eflag = false;			/* output sam format */
static bool Xflag = false;			/* print insert histogram */
static bool Yflag = false;			/* print genome projection histogram */
static bool Vflag = true;			/* automatic genome index trimming */

/* Mate Pairs */
static int	pair_mode		= DEF_PAIR_MODE;
static int	min_insert_size		= DEF_MIN_INSERT_SIZE;
static int	max_insert_size		= DEF_MAX_INSERT_SIZE;
static llint	insert_histogram[100];
static int	insert_histogram_bucket_size = 1;

/* Statistics */
static llint	nreads;
static llint	total_reads_matched;
static llint	total_pairs_matched;
static llint	total_single_matches;
static llint	total_paired_matches;
static llint	total_dup_single_matches;			/* number of duplicate hits */
static llint	total_dup_paired_matches;
static llint	total_work_usecs;
static llint	map_usecs;

static count_t mem_genomemap;

/* files to use when saving and loading genome maps */
char *		save_file = NULL;
char *		load_file = NULL;

/* Kmer to genome index */
static uint32_t ***genomemap;
static uint32_t **genomemap_len;

/* offset info for genome contigs */
static uint32_t *contig_offsets;
static char **contig_names = NULL;
static int	num_contigs;

/* Genomic sequence, stored in 32-bit words, first is in the LSB */
static uint32_t **genome_contigs;			/* genome -- always in letter */
static uint32_t **genome_contigs_rc;			/* reverse complemets */
static uint32_t **genome_cs_contigs;
static int *	genome_initbp;
static uint32_t	 *genome_len;

static bool      genome_is_rna = false;		/* is genome RNA (has uracil)?*/

/* constants for thread control */
static int	num_threads = DEF_NUM_THREADS;
static int	chunk_size = DEF_CHUNK_SIZE;

static llint	scan_ticks[50];
static llint	wait_ticks[50];


/* kmer_to_mapidx function */
static uint32_t (*kmer_to_mapidx)(uint32_t *, int) = NULL;


/* seed management */
int			n_seeds = 0;
struct seed_type *	seed = NULL;
uint32_t * *		seed_hash_mask = NULL;
int			max_seed_span = 0;
int			min_seed_span = MAX_SEED_SPAN;
int			avg_seed_span = 0;


/* pulled off the web; this may or may not be any good */
static uint32_t
hash(uint32_t a)
{
  a = (a+0x7ed55d16) + (a<<12);
  a = (a^0xc761c23c) ^ (a>>19);
  a = (a+0x165667b1) + (a<<5);
  a = (a+0xd3a2646c) ^ (a<<9);
  a = (a+0xfd7046c5) + (a<<3);
  a = (a^0xb55a4f09) ^ (a>>16);
  return a;
}

/* hash-based version or kmer -> map index function for larger seeds */
uint32_t
kmer_to_mapidx_hash(uint32_t *kmerWindow, int sn)
{
  static uint32_t maxidx = ((uint32_t)1 << 2*HASH_TABLE_POWER) - 1;
  uint32_t mapidx = 0;
  int i;

  assert(seed_hash_mask != NULL);

  for (i = 0; i < BPTO32BW(max_seed_span); i++)
    mapidx = hash((kmerWindow[i] & seed_hash_mask[sn][i]) ^ mapidx);

  return mapidx & maxidx;
}

/*
 * Compress the given kmer into an index in 'readmap' according to the seed.
 * While not optimal, this is only about 20% of the spaced seed scan time.
 *
 * This is the original version for smaller seeds.
 *
 * XXX- This algorithm only considers bases 0-3, which implies overlap
 *      when we have other bases (mainly uracil, but also wobble codes).
 *      This won't affect sensitivity, but may cause extra S-W calls.
 */
uint32_t
kmer_to_mapidx_orig(uint32_t *kmerWindow, int sn)
{
  bitmap_type a = seed[sn].mask[0];
  uint32_t mapidx = 0;
  int i = 0;

  do {
    if ((a & 0x1) == 0x1) {
      mapidx <<= 2;
      mapidx |= ((kmerWindow[i/8] >> (i%8)*4) & 0x3);
    }
    a >>= 1;
    i++;
  } while (a != 0x0);

  assert(mapidx < power(4, seed[sn].weight));

  return mapidx;
}


/* Seed management */
static bool
add_spaced_seed(char const * seed_string)
{
  int i;

  seed = (struct seed_type *)xrealloc(seed, sizeof(struct seed_type) * (n_seeds + 1));
  seed[n_seeds].mask[0] = 0x0;
  seed[n_seeds].span = strlen(seed_string);
  seed[n_seeds].weight = strchrcnt(seed_string, '1');

  if (seed[n_seeds].span < 1
      || seed[n_seeds].span > MAX_SEED_SPAN
      || seed[n_seeds].weight < 1
      || (int)strchrcnt(seed_string, '0') != seed[n_seeds].span - seed[n_seeds].weight)
    return false;

  for (i = 0; i < seed[n_seeds].span; i++)
    bitmap_prepend(seed[n_seeds].mask, 1, (seed_string[i] == '1' ? 1 : 0));

  max_seed_span = MAX(max_seed_span, seed[n_seeds].span);
  min_seed_span = MIN(min_seed_span, seed[n_seeds].span);

  n_seeds++;

  avg_seed_span = 0;
  for(i =0; i < n_seeds;i++){
    avg_seed_span += seed[i].span;
  }
  avg_seed_span = avg_seed_span/n_seeds;

  return true;
}

static void
load_default_seeds() {
  int i;

  n_seeds = 0;
  switch(shrimp_mode) {
  case MODE_COLOUR_SPACE:
    for (i = 0; i < default_spaced_seeds_cs_cnt; i++)
      add_spaced_seed(default_spaced_seeds_cs[i]);
    break;
  case MODE_LETTER_SPACE:
    for (i = 0; i < default_spaced_seeds_ls_cnt; i++)
      add_spaced_seed(default_spaced_seeds_ls[i]);
    break;
  case MODE_HELICOS_SPACE:
    fprintf(stderr, "error: helicos mode not implemented\n");
    exit(1);
    break;
  }
}

static void
init_seed_hash_mask(void)
{
  int i, sn;

  seed_hash_mask = (uint32_t **)xmalloc(sizeof(seed_hash_mask[0]) * n_seeds);
  for (sn = 0; sn < n_seeds; sn++) {
    seed_hash_mask[sn] = (uint32_t *)xcalloc(sizeof(seed_hash_mask[sn][0]) * BPTO32BW(max_seed_span));

    for (i = seed[sn].span - 1; i >= 0; i--)
      bitfield_prepend(seed_hash_mask[sn], max_seed_span,
		       bitmap_extract(seed[sn].mask, 1, i) == 1? 0xf : 0x0);
  }
}

static char *
seed_to_string(int sn)
{
  static char buffer[100];
  bitmap_type tmp;
  int i;

  buffer[seed[sn].span] = 0;
  for (i = seed[sn].span - 1, tmp = seed[sn].mask[0];
       i >= 0;
       i--, tmp >>= 1) {
    if (bitmap_extract(&tmp, 1, 0) == 1)
      buffer[i] = '1';
    else
      buffer[i] = '0';
  }

  return buffer;
}

static bool
valid_spaced_seeds()
{
  int sn;

  for (sn = 0; sn < n_seeds; sn++) {
    if (!Hflag && seed[sn].weight > MAX_SEED_WEIGHT)
      return false;

    if (Hflag && (seed[sn].span > MAX_HASH_SEED_SPAN ||
		  seed[sn].weight > MAX_HASH_SEED_WEIGHT))
      return false;
  }

  return true;
}


/* get contig number from absolute index */
static inline void
get_contig_num(uint32_t idx, int * cn) {
  *cn = 0;
  while (*cn < num_contigs - 1
	 && idx >= contig_offsets[*cn + 1])
    (*cn)++;

  assert(contig_offsets[*cn] <= idx && idx < contig_offsets[*cn] + genome_len[*cn]);
}


/*
 * Loading and saving the genome projection.
 */

static bool save_genome_map_seed(const char *file, int sn) {
  gzFile fp = gzopen(file, "wb");
  if (fp == NULL){
    return false;
  }

  // shrimp_mode
  uint32_t m;
  m = (uint32_t)shrimp_mode;
  xgzwrite(fp, &m, sizeof(uint32_t));

  // Hflag
  uint32_t h = (uint32_t)Hflag;
  xgzwrite(fp, &h, sizeof(uint32_t));

  // Seed
  xgzwrite(fp, &seed[sn], sizeof(seed_type));

  // genomemap_len
  uint32_t capacity;
  if (Hflag) {
    capacity = (uint32_t)power4(HASH_TABLE_POWER);
  } else {
    capacity = (uint32_t)power4(seed[sn].weight);
  }
  xgzwrite(fp, genomemap_len[sn], sizeof(genomemap_len[0][0]) * capacity);

  // total
  uint32_t total = 0;
  uint32_t j;
  for (j = 0; j < capacity; j++){
    total += genomemap_len[sn][j];
  }
  xgzwrite(fp, &total, sizeof(uint32_t)); //TODO do not need to write this but makes things simpler

  // genome_map
  // TODO if memory usage to high change this
  /*
    uint32_t *map,*ptr;
    map = (uint32_t *)xmalloc(sizeof(uint32_t)*total);
    ptr = map;
    for (i = 0;i<capacity;i++){
    memcpy(ptr,genomemap[sn][i],sizeof(uint32_t)*genomemap_len[sn][i]);
    ptr += genomemap_len[sn][i];
    }
    xgzwrite(fp,map,sizeof(uint32_t)*total);
    free(map);
  */

  for (j = 0; j < capacity; j++) {
    xgzwrite(fp, (void *)genomemap[sn][j], sizeof(genomemap[0][0][0]) * genomemap_len[sn][j]);
  }

  gzclose(fp);
  return true;
}

static bool load_genome_map_seed(const char *file){
  int i;
  uint32_t j;
  uint32_t total;

  gzFile fp = gzopen(file, "rb");
  if (fp == NULL){
    fprintf(stderr,"Could not open file [%s]\n",file);
    return false;
  }

  // shrimp_mode
  uint32_t m;
  xgzread(fp, &m, sizeof(uint32_t));
  if(m != (uint32_t)shrimp_mode) {
    fprintf(stderr,"Shrimp mode in file %s does not match\n",file);
  }

  // Hflag
  uint32_t h;
  xgzread(fp, &h, sizeof(uint32_t));
  if (h != (uint32_t)Hflag){
    fprintf(stderr,"Hash settings do not match in file %s\n",file);
  }

  // Seed
  int sn = n_seeds;
  n_seeds++;
  seed = (seed_type *)xrealloc(seed, sizeof(seed_type) * n_seeds);
  genomemap_len = (uint32_t **)xrealloc_c(genomemap_len,
					  sizeof(genomemap_len[0]) * n_seeds,
					  sizeof(genomemap_len[0]) * (n_seeds - 1),
					  &mem_genomemap);
  genomemap = (uint32_t ***)xrealloc_c(genomemap,
				       sizeof(genomemap[0]) * n_seeds,
				       sizeof(genomemap[0]) * (n_seeds - 1),
				       &mem_genomemap);
  xgzread(fp,seed + sn,sizeof(seed_type));

  max_seed_span = MAX(max_seed_span, seed[sn].span);
  min_seed_span = MIN(min_seed_span, seed[sn].span);

  avg_seed_span = 0;
  for(i = 0; i < n_seeds; i++) {
    avg_seed_span += seed[i].span;
  }
  avg_seed_span = avg_seed_span/n_seeds;

  // genomemap_len
  uint32_t capacity;
  if(Hflag) {
    capacity = (uint32_t)power4(HASH_TABLE_POWER);
  } else {
    capacity = (uint32_t)power4(seed[sn].weight);
  }
  genomemap_len[sn] = (uint32_t *)xmalloc_c(sizeof(genomemap_len[0][0]) * capacity, &mem_genomemap);
  genomemap[sn] = (uint32_t **)xmalloc_c(sizeof(genomemap[0][0]) * capacity, &mem_genomemap);
  xgzread(fp, genomemap_len[sn], sizeof(uint32_t) * capacity);

  // total
  xgzread(fp, &total, sizeof(uint32_t)); //TODO do not need to write this but makes things simpler

  // genome_map
  uint32_t * map;
  map = (uint32_t *)xmalloc_c(sizeof(uint32_t) * total, &mem_genomemap);
  gzread(fp,map,sizeof(uint32_t) * total);
  uint32_t * ptr;
  ptr = map;

  for (j = 0; j < capacity; j++) {
    genomemap[sn][j] = ptr;
    ptr += genomemap_len[sn][j];
  }

  gzclose(fp);
  return true;
}

static bool save_genome_map(const char *prefix) {
  char * name;
  name = (char *)xmalloc(strlen(prefix) + n_seeds + 10);

  int sn;
  for(sn = 0; sn < n_seeds; sn++) {
    sprintf(name,"%s.seed.%d", prefix, sn);
    save_genome_map_seed(name, sn);
  }

  sprintf(name, "%s.genome", prefix);
  gzFile fp = gzopen(name, "wb");
  if (fp == NULL){
    return false;
  }

  //shrimp mode
  uint32_t m;
  m = (uint32_t)(shrimp_mode);
  xgzwrite(fp,&m,sizeof(uint32_t));

  //Hflag
  uint32_t h = (uint32_t)Hflag;
  xgzwrite(fp,&h,sizeof(uint32_t));

  // num contigs
  xgzwrite(fp,&num_contigs,sizeof(uint32_t));

  // genome_len
  xgzwrite(fp,genome_len,sizeof(uint32_t)*num_contigs);

  // contig_offsets
  xgzwrite(fp,contig_offsets,sizeof(uint32_t)*num_contigs);

  //names / total
  int i;
  uint32_t total = 0;
  for(i = 0; i < num_contigs; i++){
    uint32_t len = (uint32_t)strlen(contig_names[i]);
    xgzwrite(fp, &len, sizeof(uint32_t));
    xgzwrite(fp, contig_names[i], len + 1);
    total += BPTO32BW(genome_len[i]);
  }
  xgzwrite(fp,&total,sizeof(uint32_t));

  //genome_contigs / genome_contigs_rc / genome_cs_contigs / genome_initbp
  /*
    uint32_t *gen, *gen_rc, *gen_cs = NULL, *ptr1, *ptr2, *ptr3 = NULL;
    gen = (uint32_t *)xmalloc(sizeof(uint32_t)*total);
    ptr1 = gen;
    gen_rc = (uint32_t *)xmalloc(sizeof(uint32_t)*total);
    ptr2 = gen_rc;
    if(shrimp_mode == MODE_COLOUR_SPACE){
    gen_cs = (uint32_t *)xmalloc(sizeof(uint32_t)*total);
    ptr3 = gen_cs;
    }

    for(i = 0; i < num_contigs; i++){
    memcpy(ptr1,genome_contigs[i],sizeof(uint32_t)*BPTO32BW(genome_len[i]));
    ptr1 += BPTO32BW(genome_len[i]);
    memcpy(ptr2,genome_contigs_rc[i],sizeof(uint32_t)*BPTO32BW(genome_len[i]));
    ptr2 += BPTO32BW(genome_len[i]);
    if (shrimp_mode == MODE_COLOUR_SPACE){
    memcpy(ptr3,genome_cs_contigs[i],sizeof(uint32_t)*BPTO32BW(genome_len[i]));
    ptr3 += BPTO32BW(genome_len[i]);
    }
    }

    xgzwrite(fp,gen,sizeof(uint32_t)*total);
    xgzwrite(fp,gen_rc,sizeof(uint32_t)*total);
    if(shrimp_mode == MODE_COLOUR_SPACE){
    xgzwrite(fp,gen_cs,sizeof(uint32_t)*total);
    xgzwrite(fp,genome_initbp,sizeof(uint32_t)*num_contigs);
    free(gen_cs);
    }
    free(gen);
    free(gen_rc);
  */

  for (i = 0; i < num_contigs; i++) {
    xgzwrite(fp, (void *)genome_contigs[i], BPTO32BW(genome_len[i]) * sizeof(uint32_t));
  }
  for (i = 0; i < num_contigs; i++) {
    xgzwrite(fp, (void *)genome_contigs_rc[i], BPTO32BW(genome_len[i]) * sizeof(uint32_t));
  }
  if (shrimp_mode == MODE_COLOUR_SPACE) {
    for (i = 0; i < num_contigs; i++) {
      xgzwrite(fp, (void *)genome_cs_contigs[i], BPTO32BW(genome_len[i]) * sizeof(uint32_t));
    }
    xgzwrite(fp, (void *)genome_initbp, num_contigs * sizeof(uint32_t));
  }

  gzclose(fp);
  return true;
}


static bool load_genome_map(const char *file){
  int i;

  gzFile fp = gzopen(file, "rb");
  if (fp == NULL){
    return false;
  }

  //shrimp mode
  uint32_t m;
  xgzread(fp, &m, sizeof(uint32_t));
  if (shrimp_mode != (shrimp_mode_t)m) {
    fprintf(stderr, "error: shrimp mode does not match genome file (%s)\n", file);
    exit(1);
  }

  //Hflag
  uint32_t h;
  xgzread(fp, &h, sizeof(uint32_t));
  Hflag = h;

  // num_contigs
  xgzread(fp, &num_contigs, sizeof(uint32_t));

  genome_len = (uint32_t *)xmalloc(sizeof(uint32_t) * num_contigs);
  contig_offsets = (uint32_t *)xmalloc(sizeof(uint32_t) * num_contigs);
  contig_names = (char **)xmalloc(sizeof(char *) * num_contigs);

  genome_contigs = (uint32_t **)xmalloc(sizeof(uint32_t *) * num_contigs);
  genome_contigs_rc = (uint32_t **)xmalloc(sizeof(uint32_t *) * num_contigs);
  if (shrimp_mode == MODE_COLOUR_SPACE) {
    genome_cs_contigs = (uint32_t **)xmalloc(sizeof(genome_cs_contigs[0]) * num_contigs);
    genome_initbp = (int *)xmalloc(sizeof(genome_initbp[0]) * num_contigs);
  }

  //genome_len
  xgzread(fp, genome_len, sizeof(uint32_t) * num_contigs);

  // contig_offfsets
  xgzread(fp, contig_offsets, sizeof(uint32_t) * num_contigs);

  // names / total
  for (i = 0; i < num_contigs; i++) {
    uint32_t len;
    xgzread(fp, &len,sizeof(uint32_t));
    contig_names[i] = (char *)xmalloc(sizeof(char) * (len + 1));
    xgzread(fp, contig_names[i], len + 1);
  }

  uint32_t total;
  xgzread(fp, &total, sizeof(uint32_t));

  //genome_contigs / genome_contigs_rc / genome_cs_contigs / genome_initbp
  uint32_t *gen, *gen_rc, *gen_cs, *ptr1, *ptr2, *ptr3 = NULL;
  gen = (uint32_t *)xmalloc(sizeof(uint32_t) * total);
  xgzread(fp, gen, sizeof(uint32_t) * total);
  ptr1 = gen;
  gen_rc = (uint32_t *)xmalloc(sizeof(uint32_t) * total);
  xgzread(fp, gen_rc, sizeof(uint32_t) * total);
  ptr2 = gen_rc;
  if (shrimp_mode == MODE_COLOUR_SPACE) {
    gen_cs = (uint32_t *)xmalloc(sizeof(uint32_t) * total);
    xgzread(fp, gen_cs, sizeof(uint32_t) * total);
    ptr3 = gen_cs;
    xgzread(fp, genome_initbp, sizeof(uint32_t) * num_contigs);
  }
  for (i = 0; i < num_contigs; i++) {
    genome_contigs[i] = ptr1;
    ptr1 += BPTO32BW(genome_len[i]);
    genome_contigs_rc[i] = ptr2;
    ptr2 += BPTO32BW(genome_len[i]);
    if (shrimp_mode == MODE_COLOUR_SPACE) {
      genome_cs_contigs[i] = ptr3;
      ptr3 += BPTO32BW(genome_len[i]);
    }
  }

  gzclose(fp);
  return true;
}


/*
 * Mapping routines
 */
static void
read_get_mapidxs_per_strand(struct read_entry * re, int st) {
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
      re->mapidx[st][sn*re->max_n_kmers + (r_idx - re->min_kmer_pos)] = kmer_to_mapidx(kmerWindow, sn);
    }
  }

  free(kmerWindow);
}


/*
 * Extract spaced kmers from read, save them in re->mapidx.
 */
static inline void
read_get_mapidxs(struct read_entry * re) {
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
read_get_restricted_anchor_list_per_strand(struct read_entry * re, int st, bool collapse) {
  int i, j, offset, sn;
  llint g_start, g_end;
  int idx_start, idx_end, k;
  uint32_t mapidx;
  int anchor_cache[re->read_len];
  uint diag;
  int l;

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
read_get_anchor_list_per_strand(struct read_entry * re, int st, bool collapse) {
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
  fprintf(stderr, "(cn:%d,st:%d,gen_st:%d,g_off:%lld,w_len:%d,scores:(%d,%d,%d),matches:%d,pair_min:%d,pair_max:%d,anchor:(%lld,%lld,%d,%d))\n",
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
read_get_hit_list_per_strand(struct read_entry * re, int match_mode, int st) {
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
    if (gend > genome_len[cn] - 1)
      gend = genome_len[cn] - 1;

    if (gend >= re->window_len)
      gstart = gend - re->window_len;
    else
      gstart = 0;

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
	if (re->anchors[st][j].y >= re->anchors[st][i].y)
	  continue;

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

	if (long_len > short_len)
	  tmp_score = short_len * match_score + b_gap_open_score
	    + (long_len - short_len) * b_gap_extend_score;
	else
	  tmp_score = short_len * match_score;

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

	if ((re->window_len - x_len)/2 < re->anchors[st][max_idx].x - contig_offsets[cn])
	  goff = (re->anchors[st][max_idx].x - contig_offsets[cn]) - (re->window_len - x_len)/2;
	else
	  goff = 0;

	if (goff + w_len > genome_len[cn])
	  goff = genome_len[cn] - w_len;

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
read_get_hit_list(struct read_entry * re, int match_mode) {
  read_get_hit_list_per_strand(re, match_mode, 0);
  read_get_hit_list_per_strand(re, match_mode, 1);

#ifdef DEBUG_HIT_LIST_CREATION
  fprintf(stderr, "Dumping hit list after creation for read:[%s]\n", re->name);
  dump_hit_list(re, 0, false, false);
  dump_hit_list(re, 1, false, false);
#endif
}


static void
read_pass1_per_strand(struct read_entry * re, bool only_paired, int st) {
  int i, j;

  f1_hash_tag++;
  j = -1; // last good hit

  for (i = 0; i < re->n_hits[st]; i++) {
    if (only_paired && re->hits[st][i].pair_min < 0)
      continue;

    // check window overlap
    if (j >= 0
	&& re->hits[st][i].cn == re->hits[st][j].cn
	&& re->hits[st][i].g_off <= re->hits[st][j].g_off + re->window_len
	  - (int)abs_or_pct(window_overlap, re->window_len)) {
      re->hits[st][i].score_vector = 0;
      continue;
    }

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
}


/*
 * Go through hit list, apply vector filter, and save top scores.
 */
static inline void
read_pass1(struct read_entry * re, bool only_paired) {
  read_pass1_per_strand(re, only_paired, 0);
  read_pass1_per_strand(re, only_paired, 1);

#ifdef DEBUG_HIT_LIST_PASS1
  fprintf(stderr, "Dumping hit list after pass1 for read:[%s]\n", re->name);
  dump_hit_list(re, 0, only_paired, false);
  dump_hit_list(re, 1, only_paired, false);
#endif
}


static void
readpair_pair_up_hits(struct read_entry * re1, struct read_entry * re2) {
  int st1, st2, i, j, k, l;
  int64_t delta_min, delta_max; // add to re1->g_off

  for (st1 = 0; st1 < 2; st1++) {
    st2 = 1 - st1; // opposite strand
    if (st1 == 0) {
      delta_min = (int64_t)re1->window_len + (int64_t)min_insert_size;
      delta_max = (int64_t)re1->window_len + (int64_t)max_insert_size;
    } else {
      delta_min = - (int64_t)max_insert_size - (int64_t)re2->window_len;
      delta_max = - (int64_t)min_insert_size - (int64_t)re2->window_len;
    }

    j = 0; // invariant: matching hit at index j or larger
    for (i = 0; i < re1->n_hits[st1]; i++) {
      // find matching hit, if any
      while (j < re2->n_hits[st2]
	     && (re2->hits[st2][j].cn < re1->hits[st1][i].cn // prev contig
		 || (re2->hits[st2][j].cn == re1->hits[st1][i].cn // same contig, but too close
		     && (int64_t)re2->hits[st2][j].g_off < (int64_t)re1->hits[st1][i].g_off + delta_min
		     )
		 )
	     )
	j++;

      k = j;
      while (k < re2->n_hits[st2]
	     && re2->hits[st2][k].cn == re1->hits[st1][i].cn
	     && (int64_t)re2->hits[st2][k].g_off <= (int64_t)re1->hits[st1][i].g_off + delta_max)
	k++;

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
  assert(pair_mode == PAIR_NONE);

  for (st = 0; st < 2; st++) {
    assert(re->n_hits[st] == 0 || re->hits[st] != NULL);

    for (i = 0; i < re->n_hits[st]; i++) {
      if (re->hits[st][i].score_vector >= (int)abs_or_pct(sw_vect_threshold, re->hits[st][i].score_max)
	  && (h->load < h->capacity
	      || ( (IS_ABSOLUTE(sw_vect_threshold)
		    && re->hits[st][i].score_vector > (int)h->array[0].key)
		   || (~IS_ABSOLUTE(sw_vect_threshold)
		       && re->hits[st][i].pct_score_vector > (int)h->array[0].key)))) {
	tmp.key = (IS_ABSOLUTE(sw_vect_threshold)? re->hits[st][i].score_vector : re->hits[st][i].pct_score_vector);
	tmp.rest.hit = &re->hits[st][i];

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
		     || (~IS_ABSOLUTE(sw_vect_threshold)
			 && (re1->hits[st1][i].pct_score_vector + re2->hits[st2][j].pct_score_vector)/2 > (int)h->array[0].key)))) {
	  tmp.key = (IS_ABSOLUTE(sw_vect_threshold)?
		     re1->hits[st1][i].score_vector + re2->hits[st2][j].score_vector
		     : (re1->hits[st1][i].pct_score_vector + re2->hits[st2][j].pct_score_vector)/2);
	  tmp.rest.hit[0] = &re1->hits[st1][i];
	  tmp.rest.hit[1] = &re2->hits[st2][j];
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


/*
 * Run full SW filter on this hit.
 */
static void
hit_run_full_sw(struct read_entry * re, struct read_hit * rh, int thresh)
{
  uint32_t * gen = NULL;

  assert(re != NULL && rh != NULL);
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
	       &rh->anchor, 1);
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
		 &rh->anchor, 1);
      assert(rh->sfrp->score == rh->score_vector);
    } else { // this wouldn't have passed the filter
      rh->sfrp->score = 0;
    }
  }
  rh->score_full = rh->sfrp->score;
  rh->pct_score_full = (100 * rh->score_full)/rh->score_max;
}
//TODO move this to utils
void reverse(char* s, char* t) {
       int l=strlen(s);
       int i;
       for (i=0; i<l; i++) {
               switch (s[i]) {
                       case 'A':
                               t[l-i-1]='T';
                               break;
                       case 'a':
                               t[l-i-1]='t';
                               break;
                       case 'C':
                               t[l-i-1]='G';
                               break;
                       case 'c':
                               t[l-i-1]='g';
                               break;
                       case 'G':
                               t[l-i-1]='C';
                               break;
                       case 'g':
                               t[l-i-1]='c';
                               break;
                       case 'T':
                               t[l-i-1]='A';
                               break;
                       case 't':
                               t[l-i-1]='a';
                               break;
                       case '-':
                               t[l-i-1]='-';
                               break;
                       default:
                               fprintf(stderr,"There has been a error in getting reverse complement of %s\n",s);
                               exit(1);
               }
       }
       t[l]='\0';
       //printf("%s vs %s\n", s , t);
       strcpy(s,t);
}
/*
 * Print given hit.
 *
 */
static inline void
hit_output(struct read_entry * re, struct read_hit * rh,struct read_entry * re_mp, struct read_hit * rh_mp, char ** output1, char ** output2, bool paired, bool first)
{
  assert(re != NULL && rh != NULL);
  assert(rh->sfrp != NULL);

  *output1 = output_normal(re->name, contig_names[rh->cn], rh->sfrp,
			   genome_len[rh->cn], shrimp_mode == MODE_COLOUR_SPACE, re->read[rh->st],
			   re->read_len, re->initbp[rh->st], rh->gen_st, Rflag);

  if(Pflag) {
    *output2 = output_pretty(re->name, contig_names[rh->cn], rh->sfrp,
			     genome_contigs[rh->cn], genome_len[rh->cn],
			     (shrimp_mode == MODE_COLOUR_SPACE), re->read[rh->st],
			     re->read_len, re->initbp[rh->st], rh->gen_st);
  } else if(Eflag){
    struct format_spec *fsp;
    char * output = output_format_line(Rflag);
    fsp = format_get_from_string(output);
    //fprintf(stdout,"5\n");
    free(output);

    struct input inp, inp_mp;
    memset(&inp, 0, sizeof(inp));
    memset(&inp_mp, 0, sizeof(inp_mp));
    //fprintf(stdout,*output1);
    //fprintf(stdout,"\n");
    input_parse_string(*output1,fsp,&inp);
    //fprintf(stdout,"%d\n",inp.flags);
    if(re_mp != NULL){
    	//fprintf(stdout,"4\n");
    	free(*output1);
		*output1 = output_normal(re_mp->name, contig_names[rh_mp->cn], rh_mp->sfrp,
    			   genome_len[rh_mp->cn], shrimp_mode == MODE_COLOUR_SPACE, re_mp->read[rh_mp->st],
    			   re_mp->read_len, re_mp->initbp[rh_mp->st], rh_mp->gen_st, Rflag);
		//fprintf(stdout,*output1);
		//fprintf(stdout,"\n");
		input_parse_string(*output1,fsp,&inp_mp);
		//fprintf(stdout,"%d\n",inp_mp.flags);
    }
    char * cigar, *read;
    uint32_t * read_bitstring; //, *read_ls_bitstring;
    cigar = (char *)xmalloc(sizeof(char)*200);
	int first_bp = 0;
    edit2cigar(inp.edit,inp.read_start,inp.read_end,inp.read_length,cigar);

    int contigstart = inp.genome_start + 1;

    if(inp.flags & INPUT_FLAG_IS_REVCMPL){
    	char * cigar_reverse = (char *)xmalloc(sizeof(char)*strlen(cigar)+1);
    	*cigar_reverse = '\0';
    	char * tmp = (char *)xmalloc(sizeof(char)*strlen(cigar)+1);
    	*tmp = '\0';
    	//char * tmp2 = (char *)xmalloc(sizeof(char)*strlen(cigar)+1);
    	char * ptr = cigar;
    	char * last = tmp;
    	for (ptr = cigar; *ptr != '\0'; ptr++){
    		*last = *ptr;
    		last ++;
    		*last = '\0';
    		if (!(*ptr <= '9' && *ptr >= '0')){
    			strcat(tmp,cigar_reverse);
    			strcpy(cigar_reverse,tmp);
    			last = tmp;
    			*tmp = '\0';
    		}
    	}
    	//fprintf(stdout,"3\n");
    	free(cigar);
    	cigar = cigar_reverse;
    	//fprintf(stdout,"2\n");
    	free(tmp);
    }
    

	read_bitstring = re->read[rh->gen_st];

    /*if (shrimp_mode == COLOUR_SPACE){
		first_bp = re->initbp[rh->gen_st];
    	uint index;
    	read_ls_bitstring = (uint32_t *)xmalloc(sizeof(uint32_t)*BPTO32BW(re->read_len + 1));
    	int current_bp = first_bp;
    	for (index = 0; index < re->read_len;index++){
    		bitfield_append(read_ls_bitstring, index, current_bp);
    		current_bp = cstols(current_bp,EXTRACT(read_bitstring,index),re->is_rna);
    	}
    	bitfield_append(read_ls_bitstring, index, current_bp);
    } else {
    	read_ls_bitstring = read_bitstring;
    }*/
    if(shrimp_mode == COLOUR_SPACE){
    	first_bp = re->initbp[rh->gen_st];
    	char * tmp_read = rh->sfrp->qralign;
    	int len = strlen(tmp_read);
    	read = (char *)xmalloc(sizeof(char)*len+1);
    	int i;
    	int j = 0;
    	for (i = 0; i < len;i++){
    		if(tmp_read[i] != '-'){
    			read[j] = tmp_read[i];
    			j++;
    		}
    	}
    	read[j] = '\0';
    	
    	tmp_read = (char *)xmalloc(sizeof(char)*len+1);
    	if (inp.flags & INPUT_FLAG_IS_REVCMPL){
    		reverse(read,tmp_read);
    	}
    	
    	//fprintf(stdout,"1\n");
    	free(tmp_read); 
    	
    	len = strlen(cigar);
    	for (i = 0; i < len; i++){
    		if (cigar[i] == 'S'){
    			cigar[i] = 'H';
    		}
    	}
	} else {
		read = readtostr(read_bitstring,re->read_len,false,0);
    }
    char *name;
    bool second = !first && paired;
    if (re_mp == NULL){
    	name = inp.read;
    } else {
    	int len = strlen(inp.read);
    	name = (char *)xmalloc(sizeof(char *)*len+1);
    	strncpy(name,inp.read,strlen(inp.read)+1);
		int i = 0;
		for (i = 0; i < len && *(inp.read + i) == *(inp_mp.read + i); i++);
		i--;
		name[i] = '\0';
    }

    int fivep = 0;
    int fivep_mp = 0;
    if (inp.flags & INPUT_FLAG_IS_REVCMPL){
    	//fprintf(stdout,"1\n");
    	fivep = inp.genome_end + 1;
    } else {
    	//fprintf(stdout,"2\n");
    	fivep = inp.genome_start;
    }

    if (inp_mp.flags & INPUT_FLAG_IS_REVCMPL){
    	//fprintf(stdout,"3\n");
    	fivep_mp = inp_mp.genome_end + 1;
    } else {
    	//fprintf(stdout,"4\n");
    	fivep_mp = inp_mp.genome_start;
    }
    int ins_size = (re_mp ==NULL)?0:(fivep_mp - fivep);

    free(*output1);
    *output1 = (char *)xmalloc(sizeof(char *)*2000);
    char *extra = *output1 + sprintf(*output1,"%s\t%i\t%s\t%u\t%i\t%s\t%s\t%u\t%i\t%s\t%s\tAS:i:%i",
	    name,
	    ((inp.flags & INPUT_FLAG_IS_REVCMPL) ? 16 : 0) | ((re_mp != NULL) && (inp.flags & INPUT_FLAG_IS_REVCMPL) ? 32 : 0) | ((paired) ? 1 : 0) | (first ? 64 : 0) | (second ? 128 : 0),
	    inp.genome,
	    contigstart,
	    255,
	    cigar,
	    ((re_mp == NULL)?"*":(strcmp(inp.genome,inp_mp.genome)== 0) ? "=": inp_mp.genome),
	    ((re_mp == NULL) ? 0:(inp_mp.genome_start + 1)),
	    ins_size,
	    read,
	    "*",
	    inp.score);
    if (shrimp_mode == COLOUR_SPACE){
    	sprintf(extra,"\tCS:Z:%s\tXX:Z:%s",readtostr(re->read[0],re->read_len,true,first_bp),rh->sfrp->qralign);
    	//free(read_ls_bitstring);
    }
    if(re_mp != NULL){
    	free(name);
    }
    free(read);
    free(cigar);
    format_free(fsp);
    //input_free(&inp);
    //input_free(&inp_mp);
  }
}


/*
 * Free sfrp for given hit.
 */
static void
hit_free_sfrp(struct read_hit * rh)
{
  assert(rh != NULL);

  if (rh->sfrp != NULL) {
    free(rh->sfrp->dbalign);
    free(rh->sfrp->qralign);
  }
  free(rh->sfrp);
  rh->sfrp = NULL;
}


/*
 * Do a final pass for given read.
 * Highest scoring matches are in scores heap.
 */
static void
read_pass2(struct read_entry * re, struct heap_unpaired * h) {
  int i;

  assert(re != NULL && h != NULL);
  assert(pair_mode == PAIR_NONE);

  /* compute full alignment scores */
  for (i = 0; i < (int)h->load; i++) {
    struct read_hit * rh = h->array[i].rest.hit;

    hit_run_full_sw(re, rh, (int)abs_or_pct(sw_full_threshold, rh->score_max));
    h->array[i].key = (IS_ABSOLUTE(sw_full_threshold)? rh->score_full : rh->pct_score_full);
  }

  /* sort scores */
  //qsort(&re->scores[1], re->scores[0].heap_elems, sizeof(re->scores[1]), score_cmp);
  
  //heap_unpaired_heapify(h);
  //heap_unpaired_heapsort(h);
  heap_unpaired_qsort(h);

  if ( (IS_ABSOLUTE(sw_full_threshold)
	&& (int)h->array[0].key >= (int)abs_or_pct(sw_full_threshold, h->array[0].rest.hit->score_max))
       || (~IS_ABSOLUTE(sw_full_threshold)
	   && (int)h->array[0].key >= (int)sw_full_threshold) ) {
#pragma omp atomic
    total_reads_matched++;
  }

  /* Output sorted list, removing any duplicates. */
  for (i = 0;
       i < (int)h->load && i < num_outputs
	 && ( (IS_ABSOLUTE(sw_full_threshold)
	       && (int)h->array[i].key >= (int)abs_or_pct(sw_full_threshold, h->array[i].rest.hit->score_max))
	      || (~IS_ABSOLUTE(sw_full_threshold)
		  && (int)h->array[i].key >= (int)sw_full_threshold) );
       i++) {
    struct read_hit * rh = h->array[i].rest.hit;
    bool dup;

    if (i == 0)
      dup = false;
    else
      dup = sw_full_results_equal(h->array[i-1].rest.hit->sfrp, rh->sfrp);

    if (!dup) {
      char * output1 = NULL, * output2 = NULL;

      re->final_matches++;

      hit_output(re, rh, NULL, NULL, &output1, &output2, false, false);

      if (!Pflag) {
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
      free(output2);

    } else {
      re->final_dup_matches++;
    }
  }

#pragma omp atomic
  total_single_matches += re->final_matches;

#pragma omp atomic
  total_dup_single_matches += re->final_dup_matches;

}


/*
 * Do a final pass for given pair of reads.
 * Highest scoring matches are in scores heap.
 */
static void
readpair_pass2(struct read_entry * re1, struct read_entry * re2, struct heap_paired * h) {
  int i, j;

  assert(re1 != NULL && re2 != NULL && h != NULL);

  /* compute full alignment scores */
  for (i = 0; i < (int)h->load; i++) {
    for (j = 0; j < 2; j++) {
      struct read_hit * rh = h->array[i].rest.hit[j];
      struct read_entry * re = (j == 0? re1 : re2);

      if (rh->score_full >= 0) // previously insepcted foot
	continue;

      hit_run_full_sw(re, rh, (int)abs_or_pct(sw_full_threshold, h->array[i].rest.hit[0]->score_max + h->array[i].rest.hit[1]->score_max)/3);
    }
    if (h->array[i].rest.hit[0]->score_full == 0 || h->array[i].rest.hit[1]->score_full == 0)
      h->array[i].key = 0;
    else
      h->array[i].key = (IS_ABSOLUTE(sw_full_threshold)?
			 h->array[i].rest.hit[0]->score_full + h->array[i].rest.hit[1]->score_full
			 : (h->array[i].rest.hit[0]->pct_score_full + h->array[i].rest.hit[1]->pct_score_full)/2);
  }

  /* sort scores */
  //qsort(&re->scores[1], re->scores[0].heap_elems, sizeof(re->scores[1]), score_cmp);
  
  //heap_paired_heapify(h);
  //heap_paired_heapsort(h);
  heap_paired_qsort(h);

  if ( (IS_ABSOLUTE(sw_full_threshold)
	&& (int)h->array[0].key >= (int)abs_or_pct(sw_full_threshold, h->array[0].rest.hit[0]->score_max + h->array[0].rest.hit[1]->score_max))
       || (~IS_ABSOLUTE(sw_full_threshold)
	   && (int)h->array[0].key >= (int)sw_full_threshold) ) {
#pragma omp atomic
    total_pairs_matched++;
  }

  /* Output sorted list, removing any duplicates. */
  for (i = 0;
       i < (int)h->load && i < num_outputs
	 && ( (IS_ABSOLUTE(sw_full_threshold)
	       && (int)h->array[i].key >= (int)abs_or_pct(sw_full_threshold,
							  h->array[i].rest.hit[0]->score_max + h->array[i].rest.hit[1]->score_max))
	      || (~IS_ABSOLUTE(sw_full_threshold)
		  && (int)h->array[i].key >= (int)sw_full_threshold) );
       i++) {
    struct read_hit * rh1 = h->array[i].rest.hit[0];
    struct read_hit * rh2 = h->array[i].rest.hit[1];
    bool dup;
    uint bucket;

    if (i == 0)
      dup = false;
    else
      dup = sw_full_results_equal(h->array[i-1].rest.hit[0]->sfrp, rh1->sfrp)
	&& sw_full_results_equal(h->array[i-1].rest.hit[1]->sfrp, rh2->sfrp);

    if (!dup) {
      char * output1 = NULL, * output2 = NULL, * output3 = NULL, * output4 = NULL;

      re1->final_matches++;

      hit_output(re1, rh1, re2, rh2, &output1, &output2,true,true);
      hit_output(re2, rh2, re1, rh1, &output3, &output4,true,false);

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

      if (Xflag) {
	if (h->array[i].rest.insert_size < min_insert_size)
	  bucket = 0;
	else if (h->array[i].rest.insert_size > max_insert_size)
	  bucket = 99;
	else
	  bucket = (uint)((h->array[i].rest.insert_size - min_insert_size) / insert_histogram_bucket_size);

	if (bucket >= 100) {
	  fprintf(stderr, "error: insert_size:%d re1->name:(%s)\n", h->array[i].rest.insert_size, re1->name);
	}
	assert(bucket < 100);

#pragma omp atomic
	insert_histogram[bucket]++;
      }

      free(output1);
      free(output2);
      free(output3);
      free(output4);
    } else {
      re1->final_dup_matches++;
    }
  }

#pragma omp atomic
  total_paired_matches += re1->final_matches;
#pragma omp atomic
  total_dup_paired_matches += re1->final_dup_matches;

}


/*
 * Free memory allocated by this read.
 */
static void
read_free(struct read_entry * re)
{
  free(re->name);
  free(re->read[0]);
  free(re->read[1]);

  free(re->mapidx[0]);
  free(re->mapidx[1]);
  free(re->anchors[0]);
  free(re->anchors[1]);
  free(re->hits[0]);
  free(re->hits[1]);

  if (re->n_ranges > 0) {
    free(re->ranges);
  }
}


static void
handle_read(read_entry *re){
  heap_unpaired h;
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
  heap_unpaired_init(&h, num_tmp_outputs);

  read_get_vector_hits(re, &h);

  DEBUG("second pass");
  if (h.load > 0)
    read_pass2(re, &h);

  // Done with this read; deallocate memory.
  for (i = 0; i < h.load; i++)
    hit_free_sfrp(h.array[i].rest.hit);
  heap_unpaired_destroy(&h);

  read_free(re);

  scan_ticks[omp_get_thread_num()] += rdtsc() - before;
}


static void
handle_readpair(struct read_entry * re1, struct read_entry * re2) {
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

  if (h.load > 0)
    readpair_pass2(re1, re2, &h);

  /* Done; free read entry */
  for (i = 0; i < h.load; i++) {
    hit_free_sfrp(h.array[i].rest.hit[0]);
    hit_free_sfrp(h.array[i].rest.hit[1]);
  }

  heap_paired_destroy(&h);

  read_free(re1);
  read_free(re2);

  scan_ticks[omp_get_thread_num()] += rdtsc() - before;
}

/*
 * Reverse read.
 *
 * This is useful for dealing with the various possibilities for matepair orientation in a unified way.
 */
static void
read_reverse(struct read_entry * re) {
  uint32_t * tmp1 = re->read[0];
  re->read[0] = re->read[1];
  re->read[1] = tmp1;
  
  int tmp2 = re->initbp[0];
  re->initbp[0] = re->initbp[1];
  re->initbp[1] = tmp2;

  re->input_strand = 1 - re->input_strand;
}

static uint
get_contig_number_from_name(char const * c)
{
  int cn;
  // hack: accept '>' for cn=0
  if (*c == '>')
    return 0;
  for (cn = 0; cn < num_contigs; cn++) {
    if (!strcmp(c, contig_names[cn]))
      break;
  }
  return cn;
}

/*
 * Compute range limitations for this read.
 */
static void
read_compute_ranges(struct read_entry * re)
{
  char * r, * r_save, * c, * c_save;
  uint g_start, g_end;
  int cn, st;

  assert(re->range_string != NULL);

  for (r = strtok_r(re->range_string, " ", &r_save); r != NULL; r = strtok_r(NULL, " ", &r_save))
    {
      c = strtok_r(r, ",", &c_save); // contig name
      if (c == NULL)
	continue;
      cn = get_contig_number_from_name(c);
      if (cn >= num_contigs)
	continue;

      c = strtok_r(NULL, ",", &c_save); // strand
      if (c == NULL)
	continue;
      if (*c == '+')
	st = 0;
      else if (*c == '-')
	st = 1;
      else
	continue;

      c = strtok_r(NULL, ",", &c_save); // g_start
      if (c == NULL)
	continue;
      g_start = (uint)atoll(c);
      if (g_start == 0)
	continue;
      g_start--;

      c = strtok_r(NULL, ",", &c_save); // g_end
      if (c == NULL)
	continue;
      g_end = (uint)atoll(c);
      if (g_end == 0)
	continue;
      g_end--;

      re->n_ranges++;
      re->ranges = (struct range_restriction *)xrealloc(re->ranges, re->n_ranges * sizeof(re->ranges[0]));
      re->ranges[re->n_ranges - 1].cn = cn;
      re->ranges[re->n_ranges - 1].st = st;
      re->ranges[re->n_ranges - 1].g_start = g_start;
      re->ranges[re->n_ranges - 1].g_end = g_end;
    }
}


/*
 * Launch the threads that will scan the reads
 */
static bool
launch_scan_threads(const char *file){
  fasta_t fasta;
  int space;

  //open the fasta file and check for errors
  if (shrimp_mode == MODE_LETTER_SPACE)
    space = LETTER_SPACE;
  else
    space = COLOUR_SPACE;
  fasta = fasta_open(file, space);
  if (fasta == NULL) {
    fprintf(stderr, "error: failded to open read file [%s]\n", file);
    return (false);
  } else {
    fprintf(stderr, "- Processing read file [%s]\n", file);
  }

  if (pair_mode != PAIR_NONE)
    assert(chunk_size % 2 == 0); // read an even number of reads

  bool more = true;

#pragma omp parallel shared(more, fasta) num_threads(num_threads)
  {
    struct read_entry * re_buffer;
    int load, i;
    uint64_t before;

    re_buffer = (struct read_entry *)xcalloc(chunk_size * sizeof(re_buffer[0]));

    while (more){
      before = rdtsc();

#pragma omp critical (fill_reads_buffer)
      {
	wait_ticks[omp_get_thread_num()] += rdtsc() - before;

	for (load = 0; load < chunk_size; load++) {
	  if (!fasta_get_next_with_range(fasta, &re_buffer[load].name, &re_buffer[load].seq, &re_buffer[load].is_rna,
					 &re_buffer[load].range_string)) {
	    more = false;
	    break;
	  }
	}

	nreads += load;
      } // end critical section

      if (pair_mode != PAIR_NONE)
	assert(load % 2 == 0); // read even number of reads

      for (i = 0; i < load; i++) {

	// if running in paired mode and first foot is ignored, ignore this one, too
	if (pair_mode != PAIR_NONE && i % 2 == 1 && re_buffer[i-1].name == NULL) {
	  free(re_buffer[i].name);
	  re_buffer[i].name = NULL;
	  free(re_buffer[i].seq);
	  continue;
	}

	if (!(strcspn(re_buffer[i].seq, "nNxX.") == strlen(re_buffer[i].seq))) {
	  // ignore this read
	  if (pair_mode != PAIR_NONE && i % 2 == 1) {
	    free(re_buffer[i-1].name);
	  }
	  free(re_buffer[i].name);
	  re_buffer[i].name = NULL;
	  free(re_buffer[i].seq);
	  continue;
	}

	re_buffer[i].read[0] = fasta_sequence_to_bitfield(fasta, re_buffer[i].seq);
	re_buffer[i].read_len = strlen(re_buffer[i].seq);
	re_buffer[i].max_n_kmers = re_buffer[i].read_len - min_seed_span + 1;
	re_buffer[i].min_kmer_pos = 0;
	if (shrimp_mode == MODE_COLOUR_SPACE){
	  re_buffer[i].read_len--;
	  re_buffer[i].max_n_kmers -= 2; // 1st color always discarded from kmers
	  re_buffer[i].min_kmer_pos = 1;
	  re_buffer[i].initbp[0] = fasta_get_initial_base(fasta,re_buffer[i].seq);
	  re_buffer[i].initbp[1] = re_buffer[i].initbp[0];
	  re_buffer[i].read[1] = reverse_complement_read_cs(re_buffer[i].read[0], (int8_t)re_buffer[i].initbp[0], (int8_t)re_buffer[i].initbp[1],
							    re_buffer[i].read_len, re_buffer[i].is_rna);
	} else {
	  re_buffer[i].read[1] = reverse_complement_read_ls(re_buffer[i].read[0], re_buffer[i].read_len, re_buffer[i].is_rna);
	}
	re_buffer[i].window_len = (uint16_t)abs_or_pct(window_len,re_buffer[i].read_len);
	re_buffer[i].input_strand = 0;

	re_buffer[i].mapidx[0] = NULL;
	re_buffer[i].mapidx[1] = NULL;
	re_buffer[i].anchors[0] = NULL;
	re_buffer[i].anchors[1] = NULL;
	re_buffer[i].hits[0] = NULL;
	re_buffer[i].hits[1] = NULL;
	re_buffer[i].ranges = NULL;
	re_buffer[i].n_ranges = 0;
	re_buffer[i].final_matches = 0;
	re_buffer[i].final_dup_matches = 0;

	if (re_buffer[i].range_string != NULL) {
	  read_compute_ranges(&re_buffer[i]);
	  free(re_buffer[i].range_string);
	  re_buffer[i].range_string = NULL;
	}

	free(re_buffer[i].seq);

	if (pair_mode == PAIR_NONE)
	  {
	    handle_read(&re_buffer[i]);
	  }
	else if (i % 2 == 1) 
	  {
	    if (pair_reverse[pair_mode][0])
	      read_reverse(&re_buffer[i-1]);
	    if (pair_reverse[pair_mode][1])
	      read_reverse(&re_buffer[i]);

	    handle_readpair(&re_buffer[i-1], &re_buffer[i]);
	  }
      }
    }

    free(re_buffer);
  } // end parallel section

  fasta_close(fasta);
  return true;
}


static void
print_genomemap_stats() {
  stat_t list_size, list_size_non0;
  int sn;
  uint64_t capacity, mapidx;
  uint max;

  uint64_t histogram[100];
  uint64_t cummulative_histogram[100];
  int bucket_size;
  int i, bucket;


  fprintf(stderr, "Genome Map stats:\n");

  for (sn = 0; sn < n_seeds; sn++) {
    if (Hflag)
      capacity = power(4, HASH_TABLE_POWER);
    else
      capacity = power(4, seed[sn].weight);

    stat_init(&list_size);
    stat_init(&list_size_non0);
    max = 0;
    for (mapidx = 0; mapidx < capacity; mapidx++) {
      if (genomemap_len[sn][mapidx] > list_cutoff) {
	stat_add(&list_size, 0);
	continue;
      }

      stat_add(&list_size, genomemap_len[sn][mapidx]);
      if (genomemap_len[sn][mapidx] > 0)
	stat_add(&list_size_non0, genomemap_len[sn][mapidx]);

      if (genomemap_len[sn][mapidx] > max)
	max = genomemap_len[sn][mapidx];
    }

    fprintf(stderr, "sn:%d weight:%d total_kmers:%llu lists:%llu (non-zero:%llu) list_sz_avg:%.2f (%.2f) list_sz_stddev:%.2f (%.2f) max:%u\n",
	    sn, seed[sn].weight, (long long unsigned int)stat_get_sum(&list_size),
	    (long long unsigned int)capacity, (long long unsigned int)stat_get_count(&list_size_non0),
	    stat_get_mean(&list_size), stat_get_mean(&list_size_non0),
	    stat_get_sample_stddev(&list_size), stat_get_sample_stddev(&list_size_non0), max);

    for (i = 0; i < 100; i++) {
      histogram[i] = 0;
    }

    bucket_size = ceil_div((max+1), 100); // values in [0..max]
    for (mapidx = 0; mapidx < capacity; mapidx++) {
      if (genomemap_len[sn][mapidx] > list_cutoff) {
	bucket = 0;
      } else {
	bucket = genomemap_len[sn][mapidx] / bucket_size;
	if (bucket >= 100)
	  bucket = 99;
      }
      histogram[bucket]++;
    }

    cummulative_histogram[0] = histogram[0];
    for (i = 1; i < 100; i++) {
      cummulative_histogram[i] = histogram[i] + cummulative_histogram[i-1];
    }

    for (i = 0; i < 100; i++) {
      fprintf(stderr, "[%d-%d]: %llu (cummulative: %.4f%%)\n", i*bucket_size, (i+1)*bucket_size,
	      (long long unsigned int)histogram[i], (((double)cummulative_histogram[i])/capacity)*100.0);
    }

  }
}


/*
 * index the kmers in the genome contained in the file.
 * This can then be used to align reads against.
 */
static bool
load_genome(char **files, int nfiles)
{
	fasta_t fasta;
	size_t seqlen, capacity;
	uint32_t *read;
	char *seq, *name;
	uint32_t *kmerWindow;
	int sn;
	char *file;
	bool is_rna;

	//allocate memory for the genome map
	genomemap = (uint32_t ***) xmalloc_c(n_seeds * sizeof(genomemap[0]),
			&mem_genomemap);
	genomemap_len = (uint32_t **)xmalloc_c(n_seeds * sizeof(genomemap_len[0]),
			&mem_genomemap);


	for (sn = 0; sn < n_seeds; sn++) {
	  if(Hflag){
	    capacity = power(4, HASH_TABLE_POWER);
	  } else {
	    capacity = power(4, seed[sn].weight);
	  }

	  genomemap[sn] = (uint32_t **)xcalloc_c(sizeof(uint32_t *) * capacity, &mem_genomemap);
	  genomemap_len[sn] = (uint32_t *)xcalloc_c(sizeof(uint32_t) * capacity, &mem_genomemap);

	}
	num_contigs = 0;
	u_int i = 0;
	int cfile;
	for(cfile = 0; cfile < nfiles; cfile++){
		file = files[cfile];
		//open the fasta file and check for errors
		fasta = fasta_open(file,LETTER_SPACE);
		if (fasta == NULL) {
			fprintf(stderr,"error: failded to open genome file [%s]\n",file);
			return (false);
		} else {
			fprintf(stderr,"- Processing genome file [%s]\n",file);
		}

		//Read the contigs and record their sizes

		while(fasta_get_next(fasta, &name, &seq, &is_rna)){
			genome_is_rna = is_rna;
			num_contigs++;
			contig_offsets = (uint32_t *)xrealloc(contig_offsets,sizeof(uint32_t)*num_contigs);
			contig_offsets[num_contigs - 1] = i;
			contig_names = (char **)xrealloc(contig_names,sizeof(char *)*num_contigs);
			contig_names[num_contigs - 1] = name;

			fprintf(stderr,"- Processing contig %s\n",name);

			if (strchr(name, '\t') != NULL || strchr(seq, '\t') != NULL) {
				fprintf(stderr, "error: tabs are not permitted in fasta names "
						"or sequences. Tag: [%s].\n", name);
				return false;
			}

			seqlen = strlen(seq);
			if (seqlen == 0) {
				fprintf(stderr, "error: genome [%s] had no sequence!\n",
						name);
				return false;
			}

			read = fasta_sequence_to_bitfield(fasta,seq);

			if (read == NULL) {
				fprintf(stderr, "error: invalid sequence; tag: [%s]\n", name);
				return false;
			}
			genome_contigs = (uint32_t **)xrealloc(genome_contigs,sizeof(uint32_t *)*num_contigs);
			genome_contigs[num_contigs-1] = read;
			genome_contigs_rc = (uint32_t **)xrealloc(genome_contigs_rc,sizeof(uint32_t *)*num_contigs);
			genome_contigs_rc[num_contigs-1] = reverse_complement_read_ls(read,seqlen,is_rna);
			if (shrimp_mode == MODE_COLOUR_SPACE){
				genome_cs_contigs = (uint32_t **)xrealloc(genome_cs_contigs,sizeof(uint32_t *)*num_contigs);
				genome_cs_contigs[num_contigs-1] = fasta_bitfield_to_colourspace(fasta,genome_contigs[num_contigs -1],seqlen,is_rna);
			}
			genome_len = (uint32_t *)xrealloc(genome_len,sizeof(uint32_t)*num_contigs);
			genome_len[num_contigs -1] = seqlen;

			if (shrimp_mode == MODE_COLOUR_SPACE){
			  genome_initbp = (int *)xrealloc(genome_initbp,sizeof(genome_initbp[0])*num_contigs);
			  genome_initbp[num_contigs - 1] = EXTRACT(read,0);
			  read = fasta_bitfield_to_colourspace(fasta,read,seqlen,is_rna);
			}
			kmerWindow = (uint32_t *)xcalloc(sizeof(kmerWindow[0])*BPTO32BW(max_seed_span));
			int load = 0;
			for ( ; i < seqlen + contig_offsets[num_contigs-1]; i++) {
				int base;

				base = EXTRACT(read, i - contig_offsets[num_contigs-1]);
				bitfield_prepend(kmerWindow, max_seed_span, base);

				//skip past any Ns or Xs
				if (base == BASE_N || base == BASE_X)
				  load = 0;
				else if (load < max_seed_span)
				  load++;;
				for (sn = 0; sn < n_seeds; sn++) {
					if (load < seed[sn].span)
						continue;


					uint32_t mapidx = kmer_to_mapidx(kmerWindow, sn);
					//increase the match count and store the location of the match
					genomemap_len[sn][mapidx]++;
					genomemap[sn][mapidx] = (uint32_t *)xrealloc_c(genomemap[sn][mapidx],
							sizeof(uint32_t) * (genomemap_len[sn][mapidx]),
							sizeof(uint32_t) * (genomemap_len[sn][mapidx] - 1),
							&mem_genomemap);
					genomemap[sn][mapidx][genomemap_len[sn][mapidx] - 1] = i - seed[sn].span + 1;

				}
			}
			free(seq);
			seq = name = NULL;

			free(kmerWindow);
		}
		fasta_close(fasta);
	}

	fprintf(stderr,"Loaded Genome\n");

	return (true);

}

/*
 * Trim long genome lists
 */
static void
trim_genome()
{
  uint capacity;
  int sn;
  uint32_t mapidx;

  for (sn = 0; sn < n_seeds; sn++) {
    if(Hflag) {
      capacity = power(4, HASH_TABLE_POWER);
    } else {
      capacity = power(4, seed[sn].weight);
    }

    for (mapidx = 0; mapidx < capacity; mapidx++) {
      if (genomemap_len[sn][mapidx] > list_cutoff) {
	genomemap_len[sn][mapidx] = 0;
	if (load_file == NULL) {
	  free(genomemap[sn][mapidx]);
	} // otherwise, this memory is block-allocated
	genomemap[sn][mapidx] = NULL;
      }
    }
  }
}


static void
print_insert_histogram()
{
  int i;
  for (i = 0; i < 100; i++) {
    fprintf(stderr, "[%d-%d]: %.2f%%\n",
	    min_insert_size + i * insert_histogram_bucket_size,
	    min_insert_size + (i + 1) * insert_histogram_bucket_size - 1,
	    total_paired_matches == 0? 0 : ((double)insert_histogram[i] / (double)total_paired_matches) * 100);
  }
}


static void
print_statistics()
{
  static char const my_tab[] = "    ";

	uint64_t f1_invocs[num_threads], f1_cells[num_threads], f1_ticks[num_threads];
	double f1_secs[num_threads], f1_cellspersec[num_threads];
	uint64_t f1_total_invocs = 0, f1_total_cells = 0;
	double f1_total_secs = 0, f1_total_cellspersec = 0;
	uint64_t f1_calls_bypassed = 0;

	uint64_t f2_invocs[num_threads], f2_cells[num_threads], f2_ticks[num_threads];
	double f2_secs[num_threads], f2_cellspersec[num_threads];
	uint64_t f2_total_invocs = 0, f2_total_cells = 0;
	double f2_total_secs = 0, f2_total_cellspersec = 0;

	double scan_secs[num_threads], readload_secs[num_threads];
	double total_scan_secs = 0, total_wait_secs = 0, total_readload_secs = 0;

	double hz;
	uint64_t fasta_load_ticks;
	fasta_stats_t fs;

	fs = fasta_stats();
	fasta_load_ticks = fs->total_ticks;
	free(fs);
	hz = cpuhz();

#pragma omp parallel num_threads(num_threads) shared(hz)
	{
	  int tid = omp_get_thread_num();

	  f1_stats(&f1_invocs[tid], &f1_cells[tid], &f1_ticks[tid], NULL);

	  f1_secs[tid] = (double)f1_ticks[tid] / hz;
	  f1_cellspersec[tid] = (double)f1_cells[tid] / f1_secs[tid];
	  if (isnan(f1_cellspersec[tid]))
	    f1_cellspersec[tid] = 0;

	  if (shrimp_mode == MODE_COLOUR_SPACE)
	    sw_full_cs_stats(&f2_invocs[tid], &f2_cells[tid], &f2_ticks[tid]);
	  else
	    sw_full_ls_stats(&f2_invocs[tid], &f2_cells[tid], &f2_ticks[tid]);

	  f2_secs[tid] = (double)f2_ticks[tid] / hz;
	  f2_cellspersec[tid] = (double)f2_cells[tid] / f2_secs[tid];
	  if (isnan(f2_cellspersec[tid]))
	    f2_cellspersec[tid] = 0;

	  scan_secs[tid] = ((double)scan_ticks[tid] / hz) - f1_secs[tid] - f2_secs[tid];
	  scan_secs[tid] = MAX(0, scan_secs[tid]);
	  readload_secs[tid] = ((double)total_work_usecs / 1.0e6) - ((double)scan_ticks[tid] / hz) - ((double)wait_ticks[tid] / hz);
	}
	f1_stats(NULL, NULL, NULL, &f1_calls_bypassed);

	fprintf(stderr, "\nStatistics:\n");

	fprintf(stderr, "%sOverall:\n", my_tab);
	fprintf(stderr, "%s%s%-24s" "%.2f seconds\n", my_tab, my_tab,
		"Load Genome Time:", (double)map_usecs / 1.0e6);
	fprintf(stderr, "%s%s%-24s" "%.2f seconds\n", my_tab, my_tab,
		"Read Mapping Time:", (double)total_work_usecs / 1.0e6);

	fprintf(stderr, "\n");

	int i;
	for(i = 0; i < num_threads; i++){
	  total_scan_secs += scan_secs[i];
	  total_readload_secs += readload_secs[i];
	  total_wait_secs += (double)wait_ticks[i] / hz;

	  f1_total_secs += f1_secs[i];
	  f1_total_invocs += f1_invocs[i];
	  f1_total_cells += f1_cells[i];

	  f2_total_secs += f2_secs[i];
	  f2_total_invocs += f2_invocs[i];
	  f2_total_cells += f2_cells[i];
	}
	f1_total_cellspersec = f1_total_secs == 0? 0 : (double)f1_total_cells / f1_total_secs;
	f2_total_cellspersec = f2_total_secs == 0? 0 : (double)f2_total_cells / f2_total_secs;

	if (Dflag) {
	  fprintf(stderr, "%sPer-Thread Stats:\n", my_tab);
	  fprintf(stderr, "%s%s" "%11s %9s %9s %25s %25s %9s\n", my_tab, my_tab,
		  "", "Read Load", "Scan", "Vector SW", "Scalar SW", "Wait");
	  fprintf(stderr, "%s%s" "%11s %9s %9s %15s %9s %15s %9s %9s\n", my_tab, my_tab,
		  "", "Time", "Time", "Invocs", "Time", "Invocs", "Time", "Time");
	  fprintf(stderr, "\n");
	  for(i = 0; i < num_threads; i++) {
	    fprintf(stderr, "%s%s" "Thread %-4d %9.2f %9.2f %15s %9.2f %15s %9.2f %9.2f\n", my_tab, my_tab,
		    i, readload_secs[i], scan_secs[i], comma_integer(f1_invocs[i]), f1_secs[i],
		    comma_integer(f2_invocs[i]), f2_secs[i], (double)wait_ticks[i] / hz);
	  }
	  fprintf(stderr, "\n");
	}

	fprintf(stderr, "%sSpaced Seed Scan:\n", my_tab);
	fprintf(stderr, "%s%s%-24s" "%.2f seconds\n", my_tab, my_tab,
		"Run-time:", total_scan_secs);

	fprintf(stderr, "\n");

	fprintf(stderr, "%sVector Smith-Waterman:\n", my_tab);
	fprintf(stderr, "%s%s%-24s" "%.2f seconds\n", my_tab, my_tab,
		"Run-time:", f1_total_secs);
	fprintf(stderr, "%s%s%-24s" "%s\n", my_tab, my_tab,
		"Invocations:", comma_integer(f1_total_invocs));
	fprintf(stderr, "%s%s%-24s" "%s\n", my_tab, my_tab,
		"Bypassed Calls:", comma_integer(f1_calls_bypassed));
	fprintf(stderr, "%s%s%-24s" "%.2f million\n", my_tab, my_tab,
		"Cells Computed:", (double)f1_total_cells / 1.0e6);
	fprintf(stderr, "%s%s%-24s" "%.2f million\n", my_tab, my_tab,
		"Cells per Second:", f1_total_cellspersec / 1.0e6);

	fprintf(stderr, "\n");

	fprintf(stderr, "%sScalar Smith-Waterman:\n", my_tab);
	fprintf(stderr, "%s%s%-24s" "%.2f seconds\n", my_tab, my_tab,
		"Run-time:", f2_total_secs);
	fprintf(stderr, "%s%s%-24s" "%s\n", my_tab, my_tab,
		"Invocations:", comma_integer(f2_total_invocs));
	fprintf(stderr, "%s%s%-24s" "%.2f million\n", my_tab, my_tab,
		"Cells Computed:", (double)f2_total_cells / 1.0e6);
	fprintf(stderr, "%s%s%-24s" "%.2f million\n", my_tab, my_tab,
		"Cells per Second:", f2_total_cellspersec / 1.0e6);

	fprintf(stderr, "\n");

	fprintf(stderr, "%sMiscellaneous Totals:\n", my_tab);
	fprintf(stderr, "%s%s%-24s" "%.2f seconds\n", my_tab, my_tab,
		"Fasta Lib Time:", (double)fasta_load_ticks / hz);
	fprintf(stderr, "%s%s%-24s" "%.2f seconds\n", my_tab, my_tab,
		"Read Load Time:", total_readload_secs);
	fprintf(stderr, "%s%s%-24s" "%.2f seconds\n", my_tab, my_tab,
		"Wait Time:", total_wait_secs);

	fprintf(stderr, "\n");

	fprintf(stderr, "%sGeneral:\n", my_tab);
	if (pair_mode == PAIR_NONE)
	  {
	    fprintf(stderr, "%s%s%-24s" "%s    (%.4f%%)\n", my_tab, my_tab,
		    "Reads Matched:",
		    comma_integer(total_reads_matched),
		    (nreads == 0) ? 0 : ((double)total_reads_matched / (double)nreads) * 100);
	    fprintf(stderr, "%s%s%-24s" "%s\n", my_tab, my_tab,
		    "Total Matches:",
		    comma_integer(total_single_matches));
	    fprintf(stderr, "%s%s%-24s" "%.2f\n", my_tab, my_tab,
		    "Avg Hits/Matched Read:",
		    (total_reads_matched == 0) ? 0 : ((double)total_single_matches / (double)total_reads_matched));
	    fprintf(stderr, "%s%s%-24s" "%s\n", my_tab, my_tab,
		    "Duplicate Hits Pruned:",
		    comma_integer(total_dup_single_matches));
	  }
	else // paired hits
	  {
	    fprintf(stderr, "%s%s%-40s" "%s    (%.4f%%)\n", my_tab, my_tab,
		    "Pairs Matched:",
		    comma_integer(total_pairs_matched),
		    (nreads == 0) ? 0 : ((double)total_pairs_matched / (double)(nreads/2)) * 100);
	    fprintf(stderr, "%s%s%-40s" "%s\n", my_tab, my_tab,
		    "Total Paired Matches:",
		    comma_integer(total_paired_matches));
	    fprintf(stderr, "%s%s%-40s" "%.2f\n", my_tab, my_tab,
		    "Avg Matches/Pair Matched:",
		    (total_pairs_matched == 0) ? 0 : ((double)total_paired_matches / (double)total_pairs_matched));
	    fprintf(stderr, "%s%s%-40s" "%s\n", my_tab, my_tab,
		    "Duplicate Paired Matches Pruned:",
		    comma_integer(total_dup_paired_matches));

	    /*
	    fprintf(stderr, "\n");

	    fprintf(stderr, "%s%s%-40s" "%s    (%.4f%%)\n", my_tab, my_tab,
		    "Additional Reads Matched Unpaired:",
		    comma_integer(total_reads_matched),
		    (nreads == 0) ? 0 : ((double)total_reads_matched / (double)nreads) * 100);
	    fprintf(stderr, "%s%s%-40s" "%s\n", my_tab, my_tab,
		    "Total Single Matches:",
		    comma_integer(total_single_matches));
	    fprintf(stderr, "%s%s%-40s" "%.2f\n", my_tab, my_tab,
		    "Avg Matches/Unpaired Matched Read:",
		    (total_reads_matched == 0) ? 0 : ((double)total_single_matches / (double)total_reads_matched));
	    fprintf(stderr, "%s%s%-40s" "%s\n", my_tab, my_tab,
		    "Duplicate Unpaired Matches Pruned:",
		    comma_integer(total_dup_single_matches));
	    */
	  }

	fprintf(stderr, "\n");

	fprintf(stderr, "%sMemory usage:\n", my_tab);
	fprintf(stderr, "%s%s%-24s" "%s\n", my_tab, my_tab,
		"Genomemap:",
		comma_integer(count_get_count(&mem_genomemap)));

	if (Xflag) {
	  print_insert_histogram();
	}
}

static void
usage(char * progname, bool full_usage){
  char *slash;
  int sn;

  load_default_seeds();

  slash = strrchr(progname, '/');
  if (slash != NULL)
    progname = slash + 1;

  fprintf(stderr, "usage: %s [parameters] [options] reads_file genome_file1 genome_file2...\n", progname);

  fprintf(stderr, "Parameters:\n");

  fprintf(stderr,
	  "    -s    Spaced Seed(s)                          (default: ");
  for (sn = 0; sn < n_seeds; sn++) {
    if (sn > 0)
      fprintf(stderr, "                                                            ");
    fprintf(stderr, "%s%s\n", seed_to_string(sn), (sn == n_seeds - 1? ")" : ","));
  }

  fprintf(stderr,
	  "    -o    Maximum Hits per Read                   (default: %d)\n",
	  DEF_NUM_OUTPUTS);
  fprintf(stderr,
	  "    -w    Match Window Length                     (default: %.02f%%)\n",
	  DEF_WINDOW_LEN);
  fprintf(stderr,
	  "    -n    Seed Matches per Window                 (default: %d)\n",
	  DEF_NUM_MATCHES);
  if (full_usage) {
  fprintf(stderr,
	  "    -l    Match Window Overlap Length             (default: %.02f%%)\n",
	  DEF_WINDOW_OVERLAP);
  fprintf(stderr,
	  "    -a    Anchor Width Limiting Full SW           (default: %d; disable: -1)\n",
	  DEF_ANCHOR_WIDTH);

  fprintf(stderr, "\n");
  fprintf(stderr,
	  "    -S    Save Genome Projection in File          (default: no)\n");
  fprintf(stderr,
	  "    -L    Load Genome Projection from File        (default: no)\n");
  fprintf(stderr,
	  "    -z    Projection List Cut-off Length          (default: %u)\n",
	  DEF_LIST_CUTOFF);
  }

  fprintf(stderr, "\n");
  fprintf(stderr,
	  "    -m    SW Match Score                          (default: %d)\n",
	  DEF_MATCH_VALUE);
  fprintf(stderr,
	  "    -i    SW Mismatch Score                       (default: %d)\n",
	  DEF_MISMATCH_VALUE);
  fprintf(stderr,
	  "    -g    SW Gap Open Score (Reference)           (default: %d)\n",
	  DEF_A_GAP_OPEN);
  fprintf(stderr,
	  "    -q    SW Gap Open Score (Query)               (default: %d)\n",
	  DEF_B_GAP_OPEN);
  fprintf(stderr,
	  "    -e    SW Gap Extend Score (Reference)         (default: %d)\n",
	  DEF_A_GAP_EXTEND);
  fprintf(stderr,
	  "    -f    SW Gap Extend Score (Query)             (default: %d)\n",
	  DEF_B_GAP_EXTEND);
  if (shrimp_mode == MODE_COLOUR_SPACE) {
  fprintf(stderr,
	  "    -x    SW Crossover Score                      (default: %d)\n",
	  DEF_XOVER_PENALTY);
  }
  fprintf(stderr,
	  "    -r    Window Generation Threshold             (default: %.02f%%)\n",
	  DEF_WINDOW_GEN_THRESHOLD);
  if (shrimp_mode == MODE_COLOUR_SPACE) {
  fprintf(stderr,
	  "    -v    SW Vector Hit Threshold                 (default: %.02f%%)\n",
	  DEF_SW_VECT_THRESHOLD);
  }
  fprintf(stderr,
	  "    -h    SW Full Hit Threshold                   (default: %.02f%%)\n",
	  DEF_SW_FULL_THRESHOLD);

  fprintf(stderr, "\n");

  fprintf(stderr,
	  "    -N    Number of Threads                       (default: %d)\n",
	  DEF_NUM_THREADS);
  if (full_usage) {
  fprintf(stderr,
	  "    -K    Thread Chunk Size                       (default: %d)\n",
	  DEF_CHUNK_SIZE);
  }

  fprintf(stderr, "\n");
  fprintf(stderr,
	  "    -p    Paired Mode                             (default: %s)\n",
	  pair_mode_string[pair_mode]);
  fprintf(stderr,
	  "    -I    Min and Max Insert Size                 (default: %d,%d)\n",
	  DEF_MIN_INSERT_SIZE, DEF_MAX_INSERT_SIZE);

  fprintf(stderr, "\n");
  fprintf(stderr, "Options:\n");

  fprintf(stderr,
	  "    -U    Perform Ungapped Alignment                    (default: disabled)\n");
  fprintf(stderr,
	  "    -C    Only Process Negative Strand (Rev. Compl.)    (default: disabled)\n");
  fprintf(stderr,
	  "    -F    Only Process Positive Strand                  (default: disabled)\n");
  fprintf(stderr,
	  "    -P    Pretty Print Alignments                       (default: disabled)\n");
  fprintf(stderr,
	  "    -E    Output SAM Format                             (default: disabled)\n");
  if (full_usage) {
  fprintf(stderr,
	  "    -R    Print Reads in Output                         (default: disabled)\n");
  fprintf(stderr,
	  "    -T    Reverse Tie-break on Negative Strand          (default: disabled)\n");
  fprintf(stderr,
	  "    -X    Print Insert Size Histogram                   (default: disabled)\n");
  fprintf(stderr,
	  "    -Y    Print Genome Projection Histogram             (default: disabled)\n");
  fprintf(stderr,
	  "    -Z    Disable Cache Bypass for SW Vector Calls      (default: enabled)\n");
  fprintf(stderr,
	  "    -H    Hash Spaced Kmers in Genome Projection        (default: disabled)\n");
  fprintf(stderr,
	  "    -D    Individual Thread Statistics                  (default: disabled)\n");
  fprintf(stderr,
	  "    -V    Disable Automatic Genome Index Trimming       (default: enabled)\n");
  }
  fprintf(stderr,
	  "    -?    Full List of Parameters and Options\n");

  exit(1);
}

static void
print_settings() {
  static char const my_tab[] = "    ";
  int sn;

  fprintf(stderr, "Settings:\n");
  fprintf(stderr, "%s%-40s%s (%d/%d)\n", my_tab,
	  (n_seeds == 1) ? "Spaced Seed (weight/span)" : "Spaced Seeds (weight/span)",
	  seed_to_string(0), seed[0].weight, seed[0].span);
  for (sn = 1; sn < n_seeds; sn++) {
    fprintf(stderr, "%s%-40s%s (%d/%d)\n", my_tab, "",
	    seed_to_string(sn), seed[sn].weight, seed[sn].span);
  }

  fprintf(stderr, "%s%-40s%d\n", my_tab, "Number of Outputs per Read:", num_outputs);
  fprintf(stderr, "%s%-40s%d\n", my_tab, "Window Generation Mode:", num_matches);

  if (IS_ABSOLUTE(window_len)) {
    fprintf(stderr, "%s%-40s%u\n", my_tab, "Window Length:", (uint)-window_len);
  } else {
    fprintf(stderr, "%s%-40s%.02f%%\n", my_tab, "Window Length:", window_len);
  }

  if (IS_ABSOLUTE(window_overlap)) {
    fprintf(stderr, "%s%-40s%u\n", my_tab, "Window Overlap Length:", (uint)-window_overlap);
  } else {
    fprintf(stderr, "%s%-40s%.02f%%\n", my_tab, "Window Overlap Length:", window_overlap);
  }

  fprintf(stderr, "\n");
  fprintf(stderr, "%s%-40s%d\n", my_tab, "SW Match Score:", match_score);
  fprintf(stderr, "%s%-40s%d\n", my_tab, "SW Mismatch Score:", mismatch_score);
  fprintf(stderr, "%s%-40s%d\n", my_tab, "SW Gap Open Score (Ref):", a_gap_open_score);
  fprintf(stderr, "%s%-40s%d\n", my_tab, "SW Gap Open Score (Qry):", b_gap_open_score);
  fprintf(stderr, "%s%-40s%d\n", my_tab, "SW Gap Extend Score (Ref):", a_gap_extend_score);
  fprintf(stderr, "%s%-40s%d\n", my_tab, "SW Gap Extend Score (Qry):", b_gap_extend_score);
  if (shrimp_mode == MODE_COLOUR_SPACE) {
    fprintf(stderr, "%s%-40s%d\n", my_tab, "SW Crossover Score:", crossover_score);
  }

  fprintf(stderr, "\n");

  if (IS_ABSOLUTE(window_gen_threshold)) {
    fprintf(stderr, "%s%-40s%u\n", my_tab,
	    "Window Generation Threshold:", (uint)-window_gen_threshold);
  } else {
    fprintf(stderr, "%s%-40s%.02f%%\n", my_tab,
	    "Window Generation Threshold:", window_gen_threshold);
  }
  if (shrimp_mode == MODE_COLOUR_SPACE) {
    if (IS_ABSOLUTE(sw_vect_threshold)) {
      fprintf(stderr, "%s%-40s%u\n", my_tab, "SW Vector Hit Threshold:", (uint)-sw_vect_threshold);
    } else {
      fprintf(stderr, "%s%-40s%.02f%%\n", my_tab, "SW Vector Hit Threshold:", sw_vect_threshold);
    }
  }
  if (IS_ABSOLUTE(sw_full_threshold)) {
    fprintf(stderr, "%s%-40s%u\n", my_tab,
	    shrimp_mode == MODE_COLOUR_SPACE? "SW Full Hit Threshold:" : "SW Hit Threshold",
	    (uint)-sw_full_threshold);
  } else {
    fprintf(stderr, "%s%-40s%.02f%%\n", my_tab,
	    shrimp_mode == MODE_COLOUR_SPACE? "SW Full Hit Threshold:" : "SW Hit Threshold",
	    sw_full_threshold);
  }

  fprintf(stderr, "\n");

  fprintf(stderr, "%s%-40s%s\n", my_tab, "Paired mode:", pair_mode_string[pair_mode]);
  if (pair_mode != PAIR_NONE) {
    fprintf(stderr, "%s%-40smin:%d max:%d\n", my_tab, "Insert sizes:", min_insert_size, max_insert_size);
    if (Xflag) {
      fprintf(stderr, "%s%-40s%d\n", my_tab, "Bucket size:", insert_histogram_bucket_size);
    }
  }

  fprintf(stderr, "\n");

  fprintf(stderr, "%s%-40s%s\n", my_tab, "Gapless mode:", gapless_sw? "yes" : "no");
  fprintf(stderr, "%s%-40s%d\n", my_tab, "Number of threads:", num_threads);
  fprintf(stderr, "%s%-40s%d\n", my_tab, "Thread chuck size:", chunk_size);
  fprintf(stderr, "%s%-40s%s\n", my_tab, "Hash Filter Calls:", hash_filter_calls? "yes" : "no");
  fprintf(stderr, "%s%-40s%d%s\n", my_tab, "Anchor Width:", anchor_width,
	  anchor_width == -1? " (disabled)" : "");
  if (list_cutoff < DEF_LIST_CUTOFF) {
  fprintf(stderr, "%s%-40s%u\n", my_tab, "Index List Cutoff Length:", list_cutoff);
  }

}


int main(int argc, char **argv){
	char **genome_files = NULL;
	int ngenome_files = 0;
	char *reads_file = NULL;

	char *progname = argv[0];
	char const * optstr = NULL;
	char *c;
	int ch;

	bool a_gap_open_set, b_gap_open_set;
	bool a_gap_extend_set, b_gap_extend_set;
	bool num_matches_set = false;

	set_mode_from_argv(argv);

	a_gap_open_set = b_gap_open_set = a_gap_extend_set = b_gap_extend_set = false;

	fprintf(stderr, "--------------------------------------------------"
			"------------------------------\n");
	fprintf(stderr, "gmapper: %s.\nSHRiMP %s\n[%s]\n", get_mode_string(),
			SHRIMP_VERSION_STRING, get_compiler());
	fprintf(stderr, "--------------------------------------------------"
			"------------------------------\n");

	//TODO -t -9 -d -Z -D -Y
	switch(shrimp_mode){
	case MODE_COLOUR_SPACE:
		optstr = "?s:n:w:l:o:p:m:i:g:q:e:f:h:r:a:z:DCEFHI:K:L:N:PRS:TUVXYZx:v:";
		break;
	case MODE_LETTER_SPACE:
		optstr = "?s:n:w:l:o:p:m:i:g:q:e:f:h:r:a:z:DCEFHI:K:L:N:PRS:TUVXYZ";
		break;
	case MODE_HELICOS_SPACE:
		fprintf(stderr,"Helicose currently unsuported\n");
		exit(1);
		break;
	}

	while ((ch = getopt(argc,argv,optstr)) != -1){
		switch (ch) {
		case 's':
			if (strchr(optarg, ',') == NULL) { // allow comma-separated seeds
				if (!add_spaced_seed(optarg)) {
					fprintf(stderr, "error: invalid spaced seed \"%s\"\n", optarg);
					exit (1);
				}
			} else {
				c = strtok(optarg, ",");
				do {
					if (!add_spaced_seed(c)) {
						fprintf(stderr, "error: invalid spaced seed \"%s\"\n", c);
						exit (1);
					}
					c = strtok(NULL, ",");
				} while (c != NULL);
			}
			break;
		case 'n':
			num_matches = atoi(optarg);
			num_matches_set = true;
			break;
		case 'w':
			window_len = atof(optarg);
			if (window_len <= 0.0) {
				fprintf(stderr, "error: invalid window "
						"length\n");
				exit(1);
			}
			if (strcspn(optarg, "%.") == strlen(optarg))
				window_len = -window_len;		//absol.
			break;
		case 'o':
			num_outputs = atoi(optarg);
			break;
		case 'm':
			match_score = atoi(optarg);
			break;
		case 'i':
			mismatch_score = atoi(optarg);
			break;
		case 'g':
			a_gap_open_score = atoi(optarg);
			a_gap_open_set = true;
			break;
		case 'q':
			b_gap_open_score = atoi(optarg);
			b_gap_open_set = true;
			break;
		case 'e':
			a_gap_extend_score = atoi(optarg);
			a_gap_extend_set = true;
			break;
		case 'f':
			b_gap_extend_score = atoi(optarg);
			b_gap_extend_set = true;
			break;
		case 'x':
			assert(shrimp_mode == MODE_COLOUR_SPACE);
			crossover_score = atoi(optarg);
			break;
		case 'h':
			sw_full_threshold = atof(optarg);
			if (sw_full_threshold < 0.0) {
				fprintf(stderr, "error: invalid s-w full "
						"hit threshold\n");
				exit(1);
			}
			if (strcspn(optarg, "%.") == strlen(optarg))
				sw_full_threshold = -sw_full_threshold;	//absol.
			break;
		case 'v':
		  assert(shrimp_mode == MODE_COLOUR_SPACE);
			sw_vect_threshold = atof(optarg);
			if (sw_vect_threshold < 0.0) {
				fprintf(stderr, "error: invalid s-w vector "
						"hit threshold\n");
				exit(1);
			}
			if (strcspn(optarg, "%.") == strlen(optarg))
				sw_vect_threshold = -sw_vect_threshold;	//absol.
			break;
		case 'r':
		  window_gen_threshold = atof(optarg);
		  if (window_gen_threshold < 0.0) {
		    fprintf(stderr, "error: invalid window generation threshold [%s]\n", optarg);
		    exit(1);
		  }
		  if (strcspn(optarg, "%.") == strlen(optarg))
		    window_gen_threshold = -window_gen_threshold;	//absol.
		  break;
		case 'C':
			if (Fflag) {
				fprintf(stderr, "error: -C and -F are mutually "
						"exclusive\n");
				exit(1);
			}
			Cflag = true;
			break;
		case 'F':
			if (Cflag) {
				fprintf(stderr, "error: -C and -F are mutually "
						"exclusive\n");
				exit(1);
			}
			Fflag = true;
			break;
		case 'H':
			Hflag = true;
			break;
		case 'P':
			Pflag = true;
			break;
		case 'R':
			Rflag = true;
			break;
		case 'T':
			Tflag = true;
			break;
			/*
			 * New options/parameters since SHRiMP 1.2.1
			 */
		case 'a':
			anchor_width = atoi(optarg);
			if (anchor_width < -1 || anchor_width >= 100) {
				fprintf(stderr, "error: anchor_width requested is invalid (%s)\n",
						optarg);
				exit(1);
			}
			break;
		case 'X':
			Xflag = true;
			break;
		case 'Y':
			Yflag = true;
			break;
		case 'l':
			window_overlap = atof(optarg);
			if (window_overlap <= 0.0) {
				fprintf(stderr, "error: invalid window overlap\n");
				exit(1);
			}
			if (strcspn(optarg, "%.") == strlen(optarg))
				window_overlap = -window_overlap;		//absol.
			break;
		case 'N':
			num_threads = atoi(optarg);
			break;
		case 'K':
			chunk_size = atoi(optarg);
			break;
		case 'S':
			save_file = optarg;
			break;
		case 'L':
			load_file = optarg;
			break;
		case 'D':
			Dflag = true;
			break;
		case '?':
			usage(progname, true);
			break;
		case 'Z':
		  hash_filter_calls = false;
		  break;
		case 'U':
		  gapless_sw = true;
		  anchor_width = 0;
		  a_gap_open_score = -255;
		  b_gap_open_score = -255;
		  hash_filter_calls = false;
		  break;
		case 'p':
		  if (!strcmp(optarg, "none")) {
		    pair_mode = PAIR_NONE;
		  } else if (!strcmp(optarg, "opp-in")) {
		    pair_mode = PAIR_OPP_IN;
		  } else if (!strcmp(optarg, "opp-out")) {
		    pair_mode = PAIR_OPP_OUT;
		  } else if (!strcmp(optarg, "col-fw")) {
		    pair_mode = PAIR_COL_FW;
		  } else if (!strcmp(optarg, "col-bw")) {
		    pair_mode = PAIR_COL_BW;
		  } else {
		    fprintf(stderr, "error: unrecognized pair mode (%s)\n", optarg);
		    exit(1);
		  }
		  break;
		case 'I':
		  c = strtok(optarg, ",");
		  if (c == NULL) {
		    fprintf(stderr, "error: format for insert sizes is \"-I 200,1000\"\n");
		    exit(1);
		  }
		  min_insert_size = atoi(c);
		  c = strtok(NULL, ",");
		  if (c == NULL) {
		    fprintf(stderr, "error: format for insert sizes is \"-I 200,1000\"\n");
		    exit(1);
		  }
		  max_insert_size = atoi(c);
		  if (min_insert_size > max_insert_size) {
		    fprintf(stderr, "error: invalid insert sizes (min:%d,max:%d)\n",
			    min_insert_size, max_insert_size);
		    exit(1);
		  }
		  break;
		case 'E':
			Eflag = true;
			break;
		case 'z':
		  list_cutoff = atoi(optarg);
		  if (list_cutoff == 0) {
		    fprintf(stderr, "error: invalid list cutoff (%s)\n", optarg);
		    exit(1);
		  }
		  break;
		case 'V':
		  Vflag = false;
		  break;
		default:
			usage(progname, false);
		}
	}

	argc -= optind;
	argv += optind;

	if (pair_mode != PAIR_NONE && !num_matches_set) {
	  num_matches = 4;
	}

	if (Xflag) {
	  if (pair_mode == PAIR_NONE) {
	    fprintf(stderr, "warning: insert histogram not available in unpaired mode; ignoring\n");
	    Xflag = false;
	  } else {
	    insert_histogram_bucket_size = ceil_div(max_insert_size - min_insert_size + 1, 100);
	  }
	}

	if(load_file != NULL && n_seeds != 0){
	  fprintf(stderr,"error: cannot specify seeds when loading genome map\n");
	  usage(progname,false);
	}

	if (n_seeds == 0 && load_file == NULL) {
	  load_default_seeds();
	}

	kmer_to_mapidx = kmer_to_mapidx_orig;
	if (Hflag){
	  kmer_to_mapidx = kmer_to_mapidx_hash;
	  init_seed_hash_mask();
	}

	if (save_file != NULL && load_file != NULL && list_cutoff == DEF_LIST_CUTOFF){
	  fprintf(stderr,"error: -L and -S allowed together only if list_cutoff is specified\n");
	  exit(1);
	}

	if (load_file != NULL && save_file != NULL)
	  { // args: none
	    if (argc != 0) {
	      fprintf(stderr, "error: when using both -L and -S, no extra files can be given\n");
	      usage(progname, false);
	    }
	  } 
	else if (load_file != NULL)
	  { // args: reads file
	    if (argc == 0) {
	      fprintf(stderr,"error: read_file not specified\n");
	      usage(progname, false);
	    } else if (argc == 1) {
	      reads_file    = argv[0];
	    } else {
	      fprintf(stderr,"error: too many arguments with -L\n");
	      usage(progname, false);
	    }
	  }
	else if (save_file != NULL)
	  { // args: genome file(s)
	    if (argc == 0){
	      fprintf(stderr, "error: genome_file(s) not specified\n");
	      usage(progname,false);
	    }
	    genome_files  = &argv[0];
	    ngenome_files = argc;
	  }
	else
	  { // args: reads file, genome file(s)
	    if (argc < 2) {
	      fprintf(stderr, "error: %sgenome_file(s) not specified\n",
		      (argc == 0) ? "reads_file, " : "");
	      usage(progname, false);
	    }
	    reads_file    = argv[0];
	    genome_files  = &argv[1];
	    ngenome_files = argc - 1;
	  }

	if (!Cflag && !Fflag) {
	  Cflag = Fflag = true;
	}

	if (pair_mode != PAIR_NONE && (!Cflag || !Fflag)) {
	  fprintf(stderr, "warning: in paired mode, both strands must be inspected; ignoring -C and -F\n");
	  Cflag = Fflag = true;
	}

	if (shrimp_mode == MODE_LETTER_SPACE) {
	  sw_vect_threshold = sw_full_threshold;
	}

	if (Eflag && Pflag) {
	  fprintf(stderr,"-E and -P are incompatable\n");
	  exit(1);
	}

	if (Eflag && Rflag) {
	  fprintf(stderr,"-E and -R are incompatable\n");
	  exit(1);
	}

	if (!valid_spaced_seeds()) {
	  fprintf(stderr, "error: invalid spaced seed\n");
	  if (!Hflag)
	    fprintf(stderr, "       for longer seeds, try using the -H flag\n");
	  exit(1);
	}

	if (!IS_ABSOLUTE(window_len) && window_len < 100.0) {
	  fprintf(stderr, "error: window length < 100%% of read length\n");
	  exit(1);
	}

	if (!IS_ABSOLUTE(window_overlap) && window_overlap > 100.0) {
	  fprintf(stderr, "warning: window overlap length > 100%% of window_length; resetting to 100%%\n");
	  window_overlap = 100.0;
	}

	if (num_matches < 1) {
	  fprintf(stderr, "error: invalid number of matches\n");
	  exit(1);
	}

	if (num_outputs < 1) {
	  fprintf(stderr, "error: invalid maximum hits per read\n");
	  exit(1);
	}

	if (a_gap_open_score > 0 || b_gap_open_score > 0) {
	  fprintf(stderr, "error: invalid gap open penalty\n");
	  exit(1);
	}

	if (a_gap_extend_score > 0 || b_gap_extend_score > 0) {
	  fprintf(stderr, "error: invalid gap extend penalty\n");
	  exit(1);
	}

	if (!IS_ABSOLUTE(sw_full_threshold) && sw_full_threshold > 100.0) {
	  fprintf(stderr, "error: invalid s-w full hit threshold\n");
	  exit(1);
	}

	if (shrimp_mode == MODE_COLOUR_SPACE && !IS_ABSOLUTE(sw_vect_threshold)
	    && sw_vect_threshold > 100.0) {
	  fprintf(stderr, "error: invalid s-w vector threshold\n");
	  exit(1);
	}

	if (!IS_ABSOLUTE(window_gen_threshold) && window_gen_threshold > 100.0) {
	  fprintf(stderr, "error: invalid window generation threshold\n");
	  exit(1);
	}

	if ((IS_ABSOLUTE(window_gen_threshold) && IS_ABSOLUTE(sw_full_threshold)
	     && -window_gen_threshold > -sw_full_threshold)
	    ||
	    (!IS_ABSOLUTE(window_gen_threshold) && !IS_ABSOLUTE(sw_full_threshold)
	     && window_gen_threshold > sw_full_threshold)) {
	  fprintf(stderr, "warning: window generation threshold is larger than sw threshold\n");
	}

	if ((a_gap_open_set && !b_gap_open_set) || (a_gap_extend_set && !b_gap_extend_set)) {
	  fputc('\n', stderr);
	}

	if (a_gap_open_set && !b_gap_open_set) {
	  fprintf(stderr, "Notice: Gap open penalty set for reference but not query; assuming symmetry.\n");
	  b_gap_open_score = a_gap_open_score;
	}
	if (a_gap_extend_set && !b_gap_extend_set) {
	  fprintf(stderr, "Notice: Gap extend penalty set for reference but not query; assuming symmetry.\n");
	  b_gap_extend_score = a_gap_extend_score;
	}

	if ((a_gap_open_set && !b_gap_open_set) || (a_gap_extend_set && !b_gap_extend_set)) {
	  fputc('\n', stderr);
	}

	if(load_file == NULL){
	  print_settings();
	}

	uint64_t before;
	before = gettimeinusecs();
	if (load_file != NULL){
		if (strchr(load_file, ',') == NULL) {
			//use prefix
			int buf_size = strlen(load_file) + 20;
			char * genome_name = (char *)xmalloc(sizeof(char)*buf_size);
			strncpy(genome_name,load_file,buf_size);
			strncat(genome_name,".genome",buf_size);
			fprintf(stderr,"Loading genome from %s\n",genome_name);
			if (!load_genome_map(genome_name)){
				fprintf(stderr, "error: loading from genome file \"%s\"\n", genome_name);
				exit (1);
			}
			free(genome_name);
			int seed_count = 0;
			char * seed_name = (char *)xmalloc(sizeof(char)*buf_size);
			char * buff = (char *)xmalloc(sizeof(char)*buf_size);
			strncpy(seed_name,load_file,buf_size);
			strncat(seed_name,".seed.",buf_size);
			sprintf(buff,"%d",seed_count);
			strncat(seed_name,buff,buf_size);
			FILE *f = fopen(seed_name,"r");
			while(f != NULL){
				fclose(f);
				fprintf(stderr,"Loading seed from %s\n",seed_name);
				if (!load_genome_map_seed(seed_name)) {
					fprintf(stderr, "error: loading from map file \"%s\"\n", seed_name);
					exit (1);
				}
				seed_count++;
				strncpy(seed_name,load_file,buf_size);
				strncat(seed_name,".seed.",buf_size);
				sprintf(buff,"%d",seed_count);
				strncat(seed_name,buff,buf_size);
				f = fopen(seed_name,"r");
			}
			free(seed_name);
			free(buff);

		} else {
			c = strtok(load_file, ",");
			fprintf(stderr,"Loading genome from %s\n",c);
			if (!load_genome_map(c)){
				fprintf(stderr, "error: loading from genome file \"%s\"\n", c);
				exit (1);
			}
			c = strtok(NULL, ",");
			do {
				fprintf(stderr,"Loading seed from %s\n",c);
				if (!load_genome_map_seed(c)) {
					fprintf(stderr, "error: loading from map file \"%s\"\n", c);
					exit (1);
				}
				c = strtok(NULL, ",");
			} while (c != NULL);
		}

		if (Hflag) {
		  kmer_to_mapidx = kmer_to_mapidx_hash;
		  init_seed_hash_mask();
		}

		print_settings();
	} else {
		if (!load_genome(genome_files,ngenome_files)){
			exit(1);
		}
	}

	map_usecs += (gettimeinusecs() - before);

	//
	// Automatic genome index trimming
	//
	if (Vflag && save_file == NULL && list_cutoff == DEF_LIST_CUTOFF) {
	  // this will be a mapping job; enable automatic trimming
	  int i, sn;
	  long long unsigned int total_genome_len = 0;
	  int max_seed_weight = 0;

	  for (i = 0; i < num_contigs; i++) {
	    total_genome_len += (long long unsigned int)genome_len[i];
	  }

	  if (Hflag) {
	    max_seed_weight = HASH_TABLE_POWER;
	  } else {
	    for (sn = 0; sn < n_seeds; sn++) {
	      if (seed[sn].weight > max_seed_weight) {
		max_seed_weight = seed[sn].weight;
	      }
	    }
	  }

	  // cutoff := max (1000, 100*(total_genome_len/4^max_seed_weight))
	  list_cutoff = 1000;
	  if ((uint32_t)((100llu * total_genome_len)/power(4, max_seed_weight)) > list_cutoff) {
	    list_cutoff = (uint32_t)((100llu * total_genome_len)/power(4, max_seed_weight));
	  }
	  fprintf(stderr, "Automatically trimming genome index lists longer than: %u\n", list_cutoff);
	}

	if (Yflag)
	  print_genomemap_stats();

	if (save_file != NULL) {
	  if (list_cutoff != DEF_LIST_CUTOFF) {
	    fprintf(stderr, "Trimming genome map lists longer than %u\n", list_cutoff);
	    trim_genome();
	  }

	  fprintf(stderr,"Saving genome map to %s\n",save_file);
	  if(save_genome_map(save_file)){
	    exit(0);
	  }
	  exit(1);
	}

	//TODO setup need max window and max read len
	int longest_read_len = 1000;
	int max_window_len = (int)abs_or_pct(window_len,longest_read_len);
#pragma omp parallel shared(longest_read_len,max_window_len,a_gap_open_score, a_gap_extend_score, b_gap_open_score, b_gap_extend_score,\
		match_score, mismatch_score,shrimp_mode,crossover_score,anchor_width) num_threads(num_threads)
	{
	  //hash_mark = 0;
	  //window_cache = (struct window_cache_entry *)xcalloc(1048576 * sizeof(window_cache[0]));

		if (f1_setup(max_window_len, longest_read_len,
			     a_gap_open_score, a_gap_extend_score, b_gap_open_score, b_gap_extend_score,
			     match_score, mismatch_score,
			     shrimp_mode == MODE_COLOUR_SPACE, false)) {
			fprintf(stderr, "failed to initialise vector "
					"Smith-Waterman (%s)\n", strerror(errno));
			exit(1);
		}

		int ret;
		if (shrimp_mode == MODE_COLOUR_SPACE) {
			/* XXX - a vs. b gap */
			ret = sw_full_cs_setup(max_window_len, longest_read_len,
					a_gap_open_score, a_gap_extend_score, match_score, mismatch_score,
					crossover_score, false, anchor_width);
		} else {
			ret = sw_full_ls_setup(max_window_len, longest_read_len,
					a_gap_open_score, a_gap_extend_score, b_gap_open_score, b_gap_extend_score,
					match_score, mismatch_score, false, anchor_width);
		}
		if (ret) {
			fprintf(stderr, "failed to initialise scalar "
					"Smith-Waterman (%s)\n", strerror(errno));
			exit(1);
		}
	}


	char * output;
	if (Eflag){
		//Print sam header
		fprintf(stdout,"@HD\tVN:%i\tSO:%s\n",1,"unsorted");

		int s;
		for(s = 0; s < num_contigs; s++){
			fprintf(stdout,"@SQ\tSN:%s\tLN:%u\n",contig_names[s],genome_len[s]);
		}
		fprintf(stdout,"@PG\tID:%s\tVN:%s\n","gmapper",SHRIMP_VERSION_STRING);
	} else {
		output = output_format_line(Rflag);
		puts(output);
		free(output);
	}

	before = gettimeinusecs();
	launch_scan_threads(reads_file);
	total_work_usecs += (gettimeinusecs() - before);

	print_statistics();
}
