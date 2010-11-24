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
#include <getopt.h>

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

/* Mode */
static int	mode_mirna		= false;

/* Parameters */
static double	window_len		= DEF_WINDOW_LEN;
static double	window_overlap		= DEF_WINDOW_OVERLAP;
static int	num_matches		= DEF_NUM_MATCHES;
static int	num_outputs		= DEF_NUM_OUTPUTS;
static int	max_alignments		= DEF_MAX_ALIGNMENTS;
static int	num_tmp_outputs		= 20 + num_outputs;
static int	anchor_width		= DEF_ANCHOR_WIDTH;
static uint32_t	list_cutoff		= DEF_LIST_CUTOFF;
static bool	gapless_sw		= DEF_GAPLESS_SW;
static bool	hash_filter_calls	= DEF_HASH_FILTER_CALLS;
static int	longest_read_len	= DEF_LONGEST_READ_LENGTH;


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
static bool strata_flag = false;		/* get only top scoring hits */
static bool Cflag = false;			/* do complement only */
static bool Fflag = false;			/* do positive (forward) only */
static bool Hflag = false;			/* use hash table, not lookup */
static bool Pflag = false;			/* pretty print results */
static bool Rflag = false;			/* add read sequence to output*/
static bool Tflag = true;			/* reverse sw full tie breaks */
static bool Dflag = false;			/* print statistics for each thread */
static bool Eflag = false;			/* output sam format */
static bool Xflag = false;			/* print insert histogram */
static bool Yflag = false;			/* print genome projection histogram */
static bool Vflag = true;			/* automatic genome index trimming */
static bool Qflag = false;			/* use fastq reads */
static bool Gflag = false;			/* global alignment flag ! */
static bool Bflag = false;			/* be like bfast - cs only! */


/* Mate Pairs */
static int	pair_mode		= DEF_PAIR_MODE;
static int	min_insert_size		= DEF_MIN_INSERT_SIZE;
static int	max_insert_size		= DEF_MAX_INSERT_SIZE;
static llint	insert_histogram[100];
static int	insert_histogram_bucket_size = 1;
static char *	reads_filename		= NULL;
static char * 	left_reads_filename	= NULL;
static char *	right_reads_filename	= NULL;
bool		single_reads_file	= true;

/* SAM stuff */
FILE*	unaligned_reads_file =NULL;
FILE*	aligned_reads_file =NULL;
bool	sam_unaligned =false;
bool	sam_half_paired = false; //output reads in paired mode that only have one mapping
bool	sam_r2 = false;
static char * sam_header_filename = NULL;
static char * sam_read_group_name = NULL;
static char * sam_sample_name = NULL;


/* Statistics */
static llint	nreads;
static llint	total_reads_matched;
static llint	total_pairs_matched;
static llint	total_reads_dropped;
static llint	total_pairs_dropped;
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
static uint32_t *contig_offsets=NULL;
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
load_default_mirna_seeds() {
  for (int i = 0; i < default_spaced_seeds_mirna_cnt; i++)
    add_spaced_seed(default_spaced_seeds_mirna[i]);
}

static bool
load_default_seeds(int weight) {
  int i;

  //n_seeds = 0;
  switch(shrimp_mode) {
  case MODE_COLOUR_SPACE:
    if (weight == 0)
      weight = default_spaced_seed_weight_cs;
    else if (weight < default_min_spaced_seed_weight_cs || weight > default_max_spaced_seed_weight_cs)
      return false;
    for (i = 0; i < default_spaced_seeds_cs_cnt[weight - default_min_spaced_seed_weight_cs]; i++)
      add_spaced_seed(default_spaced_seeds_cs[weight - default_min_spaced_seed_weight_cs][i]);
    break;
  case MODE_LETTER_SPACE:
    if (weight == 0)
      weight = default_spaced_seed_weight_ls;
    else if (weight < default_min_spaced_seed_weight_ls || weight > default_max_spaced_seed_weight_ls)
      return false;
    for (i = 0; i < default_spaced_seeds_ls_cnt[weight - default_min_spaced_seed_weight_ls]; i++)
      add_spaced_seed(default_spaced_seeds_ls[weight - default_min_spaced_seed_weight_ls][i]);
    break;
  case MODE_HELICOS_SPACE:
    assert(0);
    break;
  }
  return true;
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
	/*
	 * This function writes the .seed.X file appropriate for loading
	 * file is the file name to write to and sn is the seed number to write
	 *
	 * The file format is a gziped binary format as follows
	 *
	 * uint32_t				: shrimp_mode
	 * uint32_t				: Hflag
	 * seed_type			: Seed
	 * uint32_t				: capacity
	 * uint32_t * capacity	: genomemap_len
	 * uint32_t				: total (= sum from 0 to capacity - 1 of genomemap_len)
	 * uint32_t * total		: genomemap (each entry of length genomemap_len)
	 *
	 */
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
  xgzwrite(fp, &total, sizeof(uint32_t));

  // genome_map
  for (j = 0; j < capacity; j++) {
    xgzwrite(fp, (void *)genomemap[sn][j], sizeof(genomemap[0][0][0]) * genomemap_len[sn][j]);
  }

  gzclose(fp);
  return true;
}

static bool load_genome_map_seed(const char *file){
	/*
	 * This function reads the .seed.X file
	 * file is the file name to read
	 *
	 * The file format is a gziped binary format as follows
	 *
	 * uint32_t				: shrimp_mode
	 * uint32_t				: Hflag
	 * seed_type			: Seed
	 * uint32_t				: capacity
	 * uint32_t * capacity	: genomemap_len
	 * uint32_t				: total (= sum from 0 to capacity - 1 of genomemap_len)
	 * uint32_t * total		: genomemap (each entry of length genomemap_len)
	 *
	 */
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
  xgzread(fp, &total, sizeof(uint32_t));

  // genome_map
  uint32_t * map;
  map = (uint32_t *)xmalloc_c(sizeof(uint32_t) * total, &mem_genomemap);
  xgzread(fp,map,sizeof(uint32_t) * total);
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
	/*
	 * This function writes the .genome and .seed.X files appropriate for loading
	 * prefix will be used to name the files eg. prefix.genome, prefix.seed.0, ...
	 *
	 * The file format for .seed.X files is descriped in save_genome_map_seed
	 *
	 * The file format for the .genome file is a gziped binary format as follows
	 *
	 * uint32_t					: shrimp_mode
	 * uint32_t					: Hflag
	 * uint32_t 				: num_contigs
	 * uint32_t * num_contigs	: genome_len (the length of each contig)
	 * uint32_t * num_contigs	: contig_offsets
	 * per contig
	 * 		uint32_t					: name_length
	 * 		char * (name_length + 1)	: name including null termination
	 * uint32_t					: total (= sum of BPTO32BW(genome_len)
	 * per contig
	 * 		uint32_t * BPTO32BW(contig_len)	: genome_contigs
	 * per contig
	 * 		uint32_t * BPTO32BW(contig_len)	: genome_contigs_rc
	 * if colour space
	 * 		per contig
	 * 			uint32_t * BPTO32BW(contig_len) : genome_cs_cntigs
	 * 		uint32_t * num_contigs	: genome_initbp
	 */
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
	/*
	 * This function loads the .genome
	 * file is the name of the file to load from
	 *
	 * The file format for the .genome file is a gziped binary format as follows
	 *
	 * uint32_t					: shrimp_mode
	 * uint32_t					: Hflag
	 * uint32_t 				: num_contigs
	 * uint32_t * num_contigs	: genome_len (the length of each contig)
	 * uint32_t * num_contigs	: contig_offsets
	 * per contig
	 * 		uint32_t					: name_length
	 * 		char * (name_length + 1)	: name including null termination
	 * uint32_t					: total (= sum of BPTO32BW(genome_len)
	 * per contig
	 * 		uint32_t * BPTO32BW(contig_len)	: genome_contigs
	 * per contig
	 * 		uint32_t * BPTO32BW(contig_len)	: genome_contigs_rc
	 * if colour space
	 * 		per contig
	 * 			uint32_t * BPTO32BW(contig_len) : genome_cs_cntigs
	 * 		uint32_t * num_contigs	: genome_initbp
	 */
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
		     || (!IS_ABSOLUTE(sw_vect_threshold)
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
	Get hit stats
*/
int
hit_edit_distance(struct read_hit * rh) {
  //find how many perfect matches, off by 1 matches and off by 2 matches
  //       int matches;                            /* # of matches */
  //        int mismatches;                         /* # of substitutions */
  //      int insertions;                         /* # of insertions */
  //      int deletions;                          /* # of deletions */
	int edit_distance=0;
	edit_distance+=rh->sfrp->mismatches;
	edit_distance+=rh->sfrp->insertions;
	edit_distance+=rh->sfrp->deletions;
	assert(edit_distance>=0);	
	return edit_distance;	
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
  rh->pct_score_full = (100 * rh->score_full)/rh->score_max;
}
//TODO move this to utils
void reverse(char* s, char* t) {
       int l=strlen(s);
       int i;
       for (i=0; i<l; i++) {
               switch (s[i]) {
                       case 'A': t[l-i-1]='T'; break;
                       case 'a': t[l-i-1]='t'; break;
                       case 'T': t[l-i-1]='A'; break;
                       case 't': t[l-i-1]='a'; break;

                       case 'C': t[l-i-1]='G'; break;
                       case 'c': t[l-i-1]='g'; break;
                       case 'G': t[l-i-1]='C'; break;
                       case 'g': t[l-i-1]='c'; break;

                       case '-': t[l-i-1]='-'; break;

                       case 'N': t[l-i-1]='N'; break;
                       case 'n': t[l-i-1]='n'; break;
		       case '.': t[l-i-1]='.'; break;

                       case 'R': t[l-i-1]='Y'; break;
                       case 'r': t[l-i-1]='y'; break;
                       case 'Y': t[l-i-1]='R'; break;
                       case 'y': t[l-i-1]='r'; break;

                       case 'S': t[l-i-1]='S'; break;
                       case 's': t[l-i-1]='s'; break;
                       case 'W': t[l-i-1]='W'; break;
                       case 'w': t[l-i-1]='w'; break;

                       case 'K': t[l-i-1]='M'; break;
                       case 'k': t[l-i-1]='m'; break;
                       case 'M': t[l-i-1]='K'; break;
                       case 'm': t[l-i-1]='k'; break;

                       case 'B': t[l-i-1]='V'; break;
                       case 'b': t[l-i-1]='v'; break;
                       case 'V': t[l-i-1]='B'; break;
                       case 'v': t[l-i-1]='b'; break;

                       case 'D': t[l-i-1]='H'; break;
                       case 'd': t[l-i-1]='h'; break;
                       case 'H': t[l-i-1]='D'; break;
                       case 'h': t[l-i-1]='d'; break;
	
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
	read_start and read_end are 1 based
*/
cigar_t * make_cigar(int read_start, int read_end , int read_length, char* qralign,char* dbalign) {
	cigar_t * cigar = (cigar_t*)xmalloc(sizeof(cigar_t));
	cigar->size=2; 
	assert(cigar->size>=2); //for start and end without loop
	cigar->ops=(char*)xmalloc(sizeof(char)*cigar->size);
	cigar->lengths=(uint32_t*)xmalloc(sizeof(uint32_t)*cigar->size);
	int used=0;
	if (read_start>1) {
		assert(cigar->size-used>0);
		cigar->ops[used]='S';
		cigar->lengths[used]=read_start-1;
		used++;
	}
	int i=0; int qralign_length=strlen(qralign);
	while (i<qralign_length) {
		int length; char op;
		if (qralign[i]=='-') {
			for (length=0; qralign[i+length]=='-' && i+length<qralign_length; length++);
			op='D';
		} else if (dbalign[i]=='-') {
			for (length=0; dbalign[i+length]=='-' && i+length<qralign_length; length++);
			op='I';
		} else {
			for (length=0; dbalign[i+length]!='-' && qralign[i+length]!='-' && i+length<qralign_length; length++);
			op='M';
		}
		while ((used+1)>=cigar->size) { //make it bigger, want to make sure have enough for one more after loop!
			assert(cigar->size!=0);
			cigar->size*=2;
			cigar->ops=(char*)xrealloc(cigar->ops,sizeof(char)*cigar->size);
			cigar->lengths=(uint32_t*)xrealloc(cigar->lengths,sizeof(uint32_t)*cigar->size);			
		}
		cigar->ops[used]=op;
		cigar->lengths[used]=length;
		i+=length;
		used++;		
	}
	if (read_end!=read_length) {
		assert(used<cigar->size); //by loop invariant and initial size >=2
		cigar->ops[used]='S';
		cigar->lengths[used]=read_length-read_end; 
		used++;
	}
	cigar->ops=(char*)xrealloc(cigar->ops,sizeof(char)*used);
	cigar->lengths=(uint32_t*)xrealloc(cigar->lengths,sizeof(uint32_t)*used);
	cigar->size=used;
	return cigar;
} 


void reverse_cigar(cigar_t * cigar) {
	char ops[cigar->size];
	uint32_t lengths[cigar->size];
	memcpy(ops,cigar->ops,sizeof(char)*cigar->size);
	memcpy(lengths,cigar->lengths,sizeof(uint32_t)*cigar->size);
	int i;	
	for (i=0; i<cigar->size; i++) {
		cigar->ops[i]=ops[cigar->size-i-1];
		cigar->lengths[i]=lengths[cigar->size-i-1];
	}
	return;
}

char* make_cigar_string(cigar_t * cigar) {
	int string_length=cigar->size; //1 char for each op
	int i;
	for (i=0; i<cigar->size; i++) {
		int j=0,length;
		for (length=cigar->lengths[i]; length>0; j++)
			length/=10;
		string_length+=j;
	}
	string_length++; // for null term
	char * ret = (char*)xmalloc(sizeof(char)*(string_length));
	int used=0;
	for (i=0; i<cigar->size; i++) {
		//printf("%d%c\n",cigar->lengths[i],cigar->ops[i]);
		used+=sprintf(ret+used,"%d%c",cigar->lengths[i],cigar->ops[i]);
	}
	//printf("%d vs %d, %d\n",used,string_length,cigar->size);
	assert(used+1==string_length);
	ret[used]='\0';
	return ret;
}

void free_cigar(cigar_t * cigar) {
	if (cigar->size>0) {
		assert(cigar->ops!=NULL);
		assert(cigar->lengths!=NULL);
		free(cigar->ops);
		free(cigar->lengths);
	}
	free(cigar);
}

/*
 * Print given hit.
 *
 */
static inline void
hit_output(struct read_entry * re, struct read_hit * rh, struct read_hit * rh_mp, char ** output1, char ** output2, bool first_in_pair, int* hits, int found_alignments)
/*
 * This function sets the strings output1 and output2 to be the output for the current read and if in sam mode its matepair
 * It is capable of outputting regular shrimp output, pretty print output, and sam output
 *
 * re is the read_entry for the current read
 * rh is the read_hit for the current read
 * rh_mp is the read_hit for the current reads mate pair
 *
 * output1 is a pointer to a string to be used for the firest line of output
 * output1 is a pointer to a string to be used for remaining output lines (used for pretty print)
 *
 * paired is true if this read is paired
 * first is true if this is the first read in the pair
 *
 */
{
  assert(re !=NULL);
  if(!Eflag) {
  	assert(rh != NULL);
  	assert(rh->sfrp != NULL);
  	*output1 = output_normal(re->name, contig_names[rh->cn], rh->sfrp,
			   genome_len[rh->cn], shrimp_mode == MODE_COLOUR_SPACE, re->read[rh->st],
			   re->read_len, re->initbp[rh->st], rh->gen_st, Rflag);
	if (Pflag) {
		//pretty print output
		*output2 = output_pretty(re->name, contig_names[rh->cn], rh->sfrp,
				     genome_contigs[rh->cn], genome_len[rh->cn],
				     (shrimp_mode == MODE_COLOUR_SPACE), re->read[rh->st],
				     re->read_len, re->initbp[rh->st], rh->gen_st);
	}
  } else {
	//TODO change this size?
	int buffer_size=MAX(longest_read_len,1000)*8;
	*output1 = (char *)xmalloc(sizeof(char *)*buffer_size);
	//qname
	char * read_name = re->name;
	char qname[strlen(read_name)+1];
	strcpy(qname,read_name);
	//flag
	int flag;
	//rname
	char * rname = "*";
	//pos
	int pos=0;
	//mapq
	int mapq=255;
	//cigar
	char * cigar="*";
	cigar_t * cigar_binary=NULL;
	//mrnm
	const char * mrnm = "*"; //mate reference name
	//mpos
	int mpos=0;
	//isize
	int isize=0;
	//seq
	assert(shrimp_mode==MODE_COLOUR_SPACE || (signed int)strlen(re->seq)==re->read_len);
	assert(shrimp_mode==MODE_LETTER_SPACE || (signed int)strlen(re->seq)==re->read_len+1);
	char seq[re->read_len+1];
	if (shrimp_mode == MODE_LETTER_SPACE) {
		int i; 
		for (i=0; i<re->read_len; i++) {
			char c = re->seq[i];				
			switch(c) {
				case 'R':
				case 'Y':
				case 'S':
				case 'W':
				case 'K':
				case 'M':
				case 'B':
				case 'D':
				case 'H':
				case 'V':
					seq[i]='N';
					break;
				default:
					if (c>='a') {
						c-=32;
					}
					seq[i]=c;
					break;	
			}
		} 
		assert(i==re->read_len);
		seq[re->read_len]='\0';
	} else {
		seq[0]='*'; seq[1]='\0';
	}
	//qual
	int read_length = re->read_len;
	char qual[read_length+10];
	strcpy(qual,"*");
	//initialize flags	
	bool paired_read = re->paired;
	struct read_entry * re_mp = re->mate_pair;
	assert(!paired_read || re_mp!=NULL);
	bool paired_alignment = paired_read && (rh!=NULL && rh_mp!=NULL); //paired mapping, not paired read!
	//bool proper_pair = (paired_read && !query_unmapped && !mate_unmapped);
	//bool query_unmapped = (re->n_hits[0] + re->n_hits[1])>0 ? false : true;
	bool query_unmapped = (rh==NULL);
	bool mate_unmapped=false;
	bool reverse_strand = false;
	bool reverse_strand_mp = false;
	int genome_end_mp=0; int genome_start_mp=0;
	if (paired_read) {
		int min_read_name_length=MIN(strlen(read_name),strlen(re_mp->name));
		int i;
		for (i=0; i<min_read_name_length; i++) {
			if (read_name[i]==re_mp->name[i]) {
				qname[i]=read_name[i];
			} else {
				break;
			}
		}
		if (i>0 && (qname[i-1]==':' || qname[i-1]=='/')) {
			i--;
		}
		qname[i]='\0';
		//mate_unmapped=(re_mp->n_hits[0]+re_mp->n_hits[1])>0 ? false : true;
		mate_unmapped= (rh_mp==NULL);
		if (!mate_unmapped) {
			//char * read_name_mp = re->name;
			//char * rname_mp = contig_names[rh_mp->cn];
			int read_start_mp = rh_mp->sfrp->read_start+1; //1based
			//int read_length_mp = re_mp->read_len;
			int read_end_mp = read_start_mp + rh_mp->sfrp->rmapped -1; //1base
			int genome_length_mp = genome_len[rh_mp->cn];
			reverse_strand_mp = (rh_mp->gen_st ==1);
			if (!reverse_strand_mp) {
				genome_start_mp = rh_mp->sfrp->genome_start+1; // 0 based -> 1 based
			} else {
				int genome_right_most_coordinate = genome_length_mp - rh_mp->sfrp->genome_start;
				//rh->sfrp->deletions is deletions in the reference
				// This is when the read has extra characters that dont match into ref
				genome_start_mp = genome_right_most_coordinate - (read_end_mp - read_start_mp - rh_mp->sfrp->deletions + rh_mp->sfrp->insertions);
			}
			genome_end_mp=genome_start_mp+rh_mp->sfrp->gmapped-1;
			mpos=genome_start_mp;
			mrnm = contig_names[rh_mp->cn];
		}
	}
	bool second_in_pair = (paired_read && !first_in_pair);
	bool primary_alignment = false;
	bool platform_quality_fail = false;
	bool pcr_duplicate = false;
	int stored_alignments = MIN(num_outputs,found_alignments); //IH
	//if the read has no mapping or if not in half_paired mode and the mate has no mapping
	if (query_unmapped || (!sam_half_paired && paired_read && mate_unmapped)) {
		mapq=0;
		if (Qflag && shrimp_mode == MODE_LETTER_SPACE ){
			strcpy(qual,re->qual);
		}
		flag = 
			( paired_read ?  0x0001 : 0) |
			( paired_alignment ? 0x0002 : 0) |
			( query_unmapped ? 0x0004 : 0) |
			( mate_unmapped ? 0x0008 : 0) |
			( reverse_strand ? 0x0010 : 0) |
			( reverse_strand_mp ? 0x0020 : 0) |
			( first_in_pair ? 0x0040 : 0) |
			( second_in_pair ? 0x0080 : 0) |
			( primary_alignment ? 0x0100 : 0) |
			( platform_quality_fail ? 0x0200 : 0) |
			( pcr_duplicate ? 0x0400 : 0);
		char *extra = *output1 + sprintf(*output1,"%s\t%i\t%s\t%u\t%i\t%s\t%s\t%u\t%i\t%s\t%s",
			qname,flag,rname,pos,mapq,cigar,mrnm,mpos,
			isize,seq,qual);
		if (shrimp_mode == MODE_COLOUR_SPACE) {
			if (Qflag) {
				extra = extra + sprintf(extra,"\tCQ:Z:%s",re->qual);
			} else {
				extra = extra + sprintf(extra,"\tCQ:Z:%s",qual);
			}
			extra = extra + sprintf(extra, "\tCS:Z:%s",re->seq);
		}
		if (sam_r2) {
			if (shrimp_mode == MODE_COLOUR_SPACE) {
				extra = extra + sprintf(extra, "\tX2:Z:%s",re_mp->seq);
			} else {
				extra = extra + sprintf(extra, "\tR2:Z:%s",re_mp->seq);
			}
		}
		if (sam_read_group_name!=NULL ){
			extra+=sprintf(extra,"\tRG:Z:%s",sam_read_group_name);
		}	
		assert(found_alignments==0);
		assert(stored_alignments==0);
		//assert(hits[0]==0);
		//assert(hits[1]==0);
		//assert(hits[2]==0);
		//assert(hits[3]==0);
		//extra = extra + sprintf(extra,"\tH0:i:%d\tH1:i:%d\tH2:i:%d\tNH:i:%d\tIH:i:%d",hits[0],hits[1],hits[2],found_alignments,stored_alignments);
		return;
	}
	assert(rh!=NULL);
	assert( !paired_read || (!query_unmapped && !mate_unmapped) || sam_half_paired);
	//start filling in the fields
	rname = contig_names[rh->cn];
	reverse_strand = (rh->gen_st == 1);
	int read_start = rh->sfrp->read_start+1; //1based
	int read_end = read_start + rh->sfrp->rmapped -1; //1base
	int genome_length = genome_len[rh->cn];
	cigar_binary = make_cigar(read_start,read_end,read_length,rh->sfrp->qralign,rh->sfrp->dbalign);

	int qralign_length=strlen(rh->sfrp->qralign);
	int i,j=0;
	int seq_length=0;
	if (shrimp_mode == MODE_LETTER_SPACE ) {
		j=read_start-1;
		seq_length=re->read_len;
	} else if (shrimp_mode == MODE_COLOUR_SPACE ) {
		seq_length=read_end-read_start+1;
	}
	assert(seq_length<=re->read_len);
	for(i=0;i<qralign_length;i++) {
		char c=rh->sfrp->qralign[i];
		if (c!='-') { 
			if (c>='a') {
				c-=32;
			}
			if (c!='A' && c!='a' &&
				c!='G' && c!='g' &&
				c!='C' && c!='c' &&
				c!='T' && c!='t' &&
				c!='N' && c!='n') {
				//see if we can figure out what its suppose to be
				c='N';
				if (rh->sfrp->dbalign[i]!='-') {
					char r = rh->sfrp->dbalign[i];
					if (r>='a') {
						r-=32;
					}
					switch (r) {
						case 'A':
							if (c=='R' || c=='W' || c=='M' || c=='D' || c=='H' || c=='V' )
								c='A';
							break;
						case 'C':
							if (c=='Y' || c=='S' || c=='M' || c=='B' || c=='H' || c=='V')
								c='C';
							break;
						case 'G':
							if (c=='R' || c=='S' || c=='K' || c=='B' || c=='D' || c=='V')
								c='G';
							break;
						case 'T':
							if (c=='Y' || c=='W' || c=='K' || c=='B' || c=='D' || c=='H') 
								c='T';
							break;
						default: 
							fprintf(stderr,"There has been an error in printing an alignment, %c\n",r);
							break;
					}
				}
			}
			seq[j++]=c;
		}
	}
	if (shrimp_mode == MODE_LETTER_SPACE ) {
		assert(j+(re->read_len-read_end)==seq_length);
		seq[j+(re->read_len-read_end)]='\0';
	} else if (shrimp_mode == MODE_COLOUR_SPACE) {
		assert(j==seq_length);
		seq[seq_length]='\0';
	}

	//if its letter space need to reverse the qual string if its backwards
	if (shrimp_mode == MODE_LETTER_SPACE) {
		//seq=re->seq;
		if (Qflag) {
			if (!reverse_strand) {
				strcpy(qual,re->qual);
			} else {
				int qual_len = strlen(re->qual); //not same as read_len, for color space reads... apperently.....
				assert((qual_len+1)<(re->read_len+10));
				int i;
				for (i=0; i<qual_len; i++) {
					qual[(qual_len-1)-i]=re->qual[i];
				}
				qual[qual_len]='\0';
			}
		}
	//else in colour space dont print a qual string
	//but get the seq differently and change 'S' to 'H' in cigar
	} else if (shrimp_mode == MODE_COLOUR_SPACE) {
		//also change 'S' in cigar to 'H'
		//clip the qual values
		for (i=0; i<cigar_binary->size; i++) {
			if (cigar_binary->ops[i]=='S') {
				cigar_binary->ops[i]='H';
			}
		}
		if (Bflag && Qflag) {
			int read_length=(read_end-read_start+1);
			for (i=0; i<read_length; i++) {
				qual[i]=re->qual[i+read_start-1];
			}
			qual[i]='\0';
			for (i=0; i<read_length-1; i++) {
				//this is different from bfast
				//qralign is already clipped! i.e. doesn't have clipped stuff and is
				//read orientation (not always on positive reference strand!)
				int first_position_mismatch = rh->sfrp->qralign[i] > 96;
				int second_position_mismatch = rh->sfrp->qralign[i+1] > 96;
				int base_qual=0;
				if (first_position_mismatch && second_position_mismatch ) {
					base_qual+=0;
				} else if (first_position_mismatch) {
					base_qual+=qual[i+1]-qual[i];
				} else if (second_position_mismatch) {
					base_qual+=qual[i]-qual[i+1]+33;
				} else {
					base_qual+=qual[i]+qual[i+1]+10-33;
				}
				base_qual=MIN('`',MAX(base_qual,'"'));
				qual[i]=base_qual;
			}
			if (reverse_strand) {
				for (i=0; i<read_length/2; i++) {
					char temp = qual[i];
					qual[i]=qual[read_length-i-1];
					qual[read_length-i-1]=temp;
				}
			}
		}
	}
	//get the pos
	int genome_start;
	if (!reverse_strand) {
		genome_start = rh->sfrp->genome_start+1; // 0 based -> 1 based
	} else {
		int genome_right_most_coordinate = genome_length - rh->sfrp->genome_start;
		//rh->sfrp->deletions is deletions in the reference
		// This is when the read has extra characters that dont match into ref
		genome_start = genome_right_most_coordinate - (read_end - read_start - rh->sfrp->deletions + rh->sfrp->insertions);
		char * tmp = (char*)xmalloc(sizeof(char)*(strlen(seq)+1));
		reverse(seq,tmp);
		free(tmp);
		reverse_cigar(cigar_binary);
	}
	int genome_end=genome_start+rh->sfrp->gmapped-1;
	pos=genome_start;
	//make the cigar string
	cigar = make_cigar_string(cigar_binary); 

	//do some stats using matepair
	if (paired_read && !mate_unmapped) {
		assert(rh_mp!=NULL && re_mp!=NULL);

		mrnm = (strcmp(rname,mrnm)==0) ? "=" : mrnm;
		//printf("%d %d %c, %d %d %c\n",genome_start, genome_end,reverse_strand ? '-' : '+' , genome_start_mp, genome_end_mp, reverse_strand_mp ? '-' : '+');
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
		isize = (fivep_mp - fivep);
	}
	flag = 
		( paired_read ?  0x0001 : 0) |
		( paired_alignment ? 0x0002 : 0) |
		( query_unmapped ? 0x0004 : 0) |
		( mate_unmapped ? 0x0008 : 0) |
		( reverse_strand ? 0x0010 : 0) |
		( reverse_strand_mp ? 0x0020 : 0) |
		( first_in_pair ? 0x0040 : 0) |
		( second_in_pair ? 0x0080 : 0) |
		( primary_alignment ? 0x0100 : 0) |
		( platform_quality_fail ? 0x0200 : 0) |
		( pcr_duplicate ? 0x0400 : 0);
	char *extra = *output1 + sprintf(*output1,"%s\t%i\t%s\t%u\t%i\t%s\t%s\t%u\t%i\t%s\t%s",
		qname,flag,rname,pos,mapq,cigar,mrnm,mpos,
		isize,seq,qual);
	extra = extra + sprintf(extra,"\tAS:i:%d\tH0:i:%d\tH1:i:%d\tH2:i:%d\tNM:i:%d\tNH:i:%d\tIH:i:%d",rh->sfrp->score,hits[0],hits[1],hits[2],rh->sfrp->mismatches+rh->sfrp->deletions+rh->sfrp->insertions,found_alignments,stored_alignments);
	if (shrimp_mode == COLOUR_SPACE){
		//TODO
		//int first_bp = re->initbp[0];
		//printf("%s vs %s\n",readtostr(re->read[0],re->read_len,true,first_bp),re->seq);
		//assert(strcmp(readtostr(re->read[0],re->read_len,true,first_bp),re->seq)==0);
		if (Qflag) {
			extra = extra + sprintf(extra,"\tCQ:Z:%s",re->qual);
		}
		extra = extra + sprintf(extra, "\tCS:Z:%s\tCM:i:%d\tXX:Z:%s",re->seq,rh->sfrp->crossovers,rh->sfrp->qralign);
	} 
	if (sam_r2) {
		if (shrimp_mode == MODE_COLOUR_SPACE) {
			extra = extra + sprintf(extra, "\tX2:Z:%s",re_mp->seq);
		} else {
			extra = extra + sprintf(extra, "\tR2:Z:%s",re_mp->seq);
		}
	}
	
	if (sam_read_group_name!=NULL) {
			extra+=sprintf(extra,"\tRG:Z:%s",sam_read_group_name);
	}
	if (cigar_binary!=NULL) {
		free_cigar(cigar_binary);
		free(cigar);
	}

    //to calculate the insert size we need to find the five' end of the reads
/*
	    inp.score);
    if(re_mp != NULL){
    	free(name);
    }
    free(read);
    free(cigar);
    format_free(fsp);*/

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
static int
read_pass2(struct read_entry * re, struct heap_unpaired * h) {
  int i;

  assert(re != NULL && h != NULL);
  assert(pair_mode == PAIR_NONE || sam_half_paired);

  /* compute full alignment scores */
  for (i = 0; i < (int)h->load; i++) {
    struct read_hit * rh = h->array[i].rest.hit;
    //fprintf(stderr,"g_off is %llu\n",rh->g_off);
    hit_run_full_sw(re, rh, (int)abs_or_pct(sw_full_threshold, rh->score_max));
    h->array[i].key = (IS_ABSOLUTE(sw_full_threshold)? rh->score_full : rh->pct_score_full);
  }

  /* sort scores */
  //qsort(&re->scores[1], re->scores[0].heap_elems, sizeof(re->scores[1]), score_cmp);
  
  //heap_unpaired_heapify(h);
  //heap_unpaired_heapsort(h);
  heap_unpaired_qsort(h);


  int found_alignments=0;
  int hits[] = {0,0,0};
  for (i=0; i<(int)h->load; i++) {
	struct sw_full_results * sfrp = h->array[i].rest.hit->sfrp;
	if (i==0) {
		sfrp->dup=false;
	} else {
		sfrp->dup=sw_full_results_equal(h->array[i-1].rest.hit->sfrp, sfrp);
	}
	if (!sfrp->dup) {
		if ( (IS_ABSOLUTE(sw_full_threshold)
			&& (int)h->array[i].key >= (int)abs_or_pct(sw_full_threshold, h->array[i].rest.hit->score_max))
		   || (!IS_ABSOLUTE(sw_full_threshold)
			   && (int)h->array[i].key >= (int)sw_full_threshold) ) {
			if (!strata_flag) {
				found_alignments++;
				//int edit_distance=re->read_len-sfrp->matches+sfrp->insertions;
				int edit_distance=sfrp->mismatches+sfrp->insertions+sfrp->deletions;
				if (0<=edit_distance && edit_distance<3) {
					hits[edit_distance]++;
				}
			} else if ((int)h->array[i].key==(int)h->array[0].key) {
				found_alignments++;
				//int edit_distance=re->read_len-sfrp->matches+sfrp->insertions;
				int edit_distance=sfrp->mismatches+sfrp->insertions+sfrp->deletions;
				if (0<=edit_distance && edit_distance<3) {
					hits[edit_distance]++;
				}
			}
		}
	}
	//fprintf(stderr,"%s\n%s\n",sfrp->dbalign,sfrp->qralign);
  }

  if ( (IS_ABSOLUTE(sw_full_threshold)
	&& (int)h->array[0].key >= (int)abs_or_pct(sw_full_threshold, h->array[0].rest.hit->score_max))
   || (!IS_ABSOLUTE(sw_full_threshold)
	   && (int)h->array[0].key >= (int)sw_full_threshold) ) {
	  if (max_alignments==0 || found_alignments<=max_alignments) {
	#pragma omp atomic
	    total_reads_matched++;
	  } else {
	#pragma omp atomic
	    total_reads_dropped++;
	  }
  }

  if (max_alignments==0 || found_alignments<=max_alignments) {
	  int outputted=0;
	  /* Output sorted list, removing any duplicates. */
	  for (i = 0;
	       i < (int)h->load && outputted < num_outputs
		 && ( (IS_ABSOLUTE(sw_full_threshold)
		       && (int)h->array[i].key >= (int)abs_or_pct(sw_full_threshold, h->array[i].rest.hit->score_max))
		      || (!IS_ABSOLUTE(sw_full_threshold)
			  && (int)h->array[i].key >= (int)sw_full_threshold) );
	       i++) {
	    struct read_hit * rh = h->array[i].rest.hit;
	    if (strata_flag && h->array[i].key != h->array[0].key) {
		break;
	    }
	    /*bool dup;

	    if (i == 0)
	      dup = false;
	    else
	      dup = sw_full_results_equal(h->array[i-1].rest.hit->sfrp, rh->sfrp);
		*/
	    if (!rh->sfrp->dup) {
	      outputted++;
	      char * output1 = NULL, * output2 = NULL;
	      re->final_matches++;
	      if (sam_half_paired) {
		int other_hits[]={0,0,0};
		if (re->first_in_pair) {
	      		hit_output(re, rh, NULL, &output1, NULL ,true,hits,found_alignments);
	      		hit_output(re->mate_pair, NULL, rh, &output2, NULL,false,other_hits,0);
		} else {
	      		hit_output(re->mate_pair, NULL, rh, &output1, NULL,true,other_hits,0);
	      		hit_output(re, rh, NULL, &output2, NULL,false,hits,found_alignments);
		}
	      } else {
	      	hit_output(re, rh, NULL,  &output1, &output2, false, hits,found_alignments);
	      }
	//TODO - maybe buffer printf output here ? might make it go faster?
	      if (!Pflag) {
	#pragma omp critical (stdout)
		{
		  fprintf(stdout, "%s\n", output1);
		  if (sam_half_paired) {
		  	fprintf(stdout, "%s\n", output2);
		  }
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

	    } else {
	      re->final_dup_matches++;
	    }
	  }
//TODO maybe have each thread own counter and then sum?
//This might be pretty heavy, maybe not because of 
//printf? 
#pragma omp atomic
  total_single_matches += re->final_matches;

#pragma omp atomic
  total_dup_single_matches += re->final_dup_matches;
   }
	return found_alignments;
}


/*
 * Do a final pass for given pair of reads.
 * Highest scoring matches are in scores heap.
 */
static int
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

  int found_alignments=0;
  int hits1[] = {0,0,0};
  int hits2[] = {0,0,0};
  for (i=0; i<(int)h->load; i++) {
	struct sw_full_results * sfrp1 = h->array[i].rest.hit[0]->sfrp;
	struct sw_full_results * sfrp2 = h->array[i].rest.hit[1]->sfrp;
	if (i==0) {
		sfrp1->dup=false;	
		sfrp2->dup=false;
	} else {
		sfrp1->dup=sw_full_results_equal(h->array[i-1].rest.hit[0]->sfrp, sfrp1);
		sfrp2->dup=sw_full_results_equal(h->array[i-1].rest.hit[1]->sfrp, sfrp2);
	}
	if ((!sfrp1->dup) || (!sfrp2->dup)) {
		//int edit_distance1=re1->read_len-sfrp1->matches+sfrp1->insertions;
		//int edit_distance2=re2->read_len-sfrp2->matches+sfrp2->insertions;
		int edit_distance1=sfrp1->mismatches+sfrp1->insertions+sfrp1->deletions;
		int edit_distance2=sfrp2->mismatches+sfrp2->insertions+sfrp1->deletions;
		if ( (IS_ABSOLUTE(sw_full_threshold)
			&& (int)h->array[i].key >= (int)abs_or_pct(sw_full_threshold, h->array[i].rest.hit[0]->score_max + h->array[i].rest.hit[1]->score_max))
		       || (!IS_ABSOLUTE(sw_full_threshold)
			   && (int)h->array[i].key >= (int)sw_full_threshold) ) {
			if (!strata_flag) {
				found_alignments++;
				if (0<=edit_distance1 && edit_distance1<3) {
					hits1[edit_distance1]++;
				}
				if (0<=edit_distance2 && edit_distance2<3) {
					hits2[edit_distance2]++;
				}
			} else if ((int)h->array[i].key==(int)h->array[0].key) {
				found_alignments++;
				if (0<=edit_distance1 && edit_distance1<3) {
					hits1[edit_distance1]++;
				}
				if (0<=edit_distance2 && edit_distance2<3) {
					hits2[edit_distance2]++;
				}
			}
		}
	}
  }

  if ( (IS_ABSOLUTE(sw_full_threshold)
	&& (int)h->array[0].key >= (int)abs_or_pct(sw_full_threshold, h->array[0].rest.hit[0]->score_max + h->array[0].rest.hit[1]->score_max))
       || (!IS_ABSOLUTE(sw_full_threshold)
	   && (int)h->array[0].key >= (int)sw_full_threshold) ) {
    if (max_alignments==0 || found_alignments<=max_alignments) {
#pragma omp atomic
    	total_pairs_matched++;
    } else {
#pragma omp atomic
	total_pairs_dropped++;
    }
  }

  int outputted=0;
  /* Output sorted list, removing any duplicates. */
  if (max_alignments==0 || found_alignments<=max_alignments) {
	  for (i = 0;
	       i < (int)h->load && outputted < num_outputs
		 && ( (IS_ABSOLUTE(sw_full_threshold)
		       && (int)h->array[i].key >= (int)abs_or_pct(sw_full_threshold,
								  h->array[i].rest.hit[0]->score_max + h->array[i].rest.hit[1]->score_max))
		      || (!IS_ABSOLUTE(sw_full_threshold)
			  && (int)h->array[i].key >= (int)sw_full_threshold) );
	       i++) {
	    struct read_hit * rh1 = h->array[i].rest.hit[0];
	    struct read_hit * rh2 = h->array[i].rest.hit[1];
	    if (strata_flag && h->array[i].key != h->array[0].key) {
		break;
	    }
	    uint bucket;

		/*
	    bool dup;
	    if (i == 0)
	      dup = false;
	    else
			dup= sw_full_results_equal(h->array[i-1].rest.hit[0]->sfrp, rh1->sfrp)
			&& sw_full_results_equal(h->array[i-1].rest.hit[1]->sfrp, rh2->sfrp);
		*/

	    if ((!rh1->sfrp->dup) || (!rh2->sfrp->dup)) {
	      outputted++;
	      char * output1 = NULL, * output2 = NULL, * output3 = NULL, * output4 = NULL;

	      re1->final_matches++;
	      assert(rh1!=NULL && rh2!=NULL);
	      hit_output(re1, rh1,  rh2, &output1, &output2,true,hits1,found_alignments);
	      hit_output(re2, rh2,  rh1, &output3, &output4,false,hits2,found_alignments);

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
	      free(output3);
	      if (Pflag) {
	      	free(output2);
	      	free(output4);
	      }
	    } else {
	      re1->final_dup_matches++;
	    }
	  }

	#pragma omp atomic
	  total_paired_matches += re1->final_matches;
	#pragma omp atomic
	  total_dup_paired_matches += re1->final_dup_matches;
   }
   return found_alignments;
}


/*
 * Free memory allocated by this read.
 */
static void 
read_free(read_entry * re) {
  free(re->name);
  free(re->seq);
  if (Qflag) {
        assert(re->qual!=NULL);
        free(re->qual);
        assert(re->plus_line!=NULL);
        free(re->plus_line);
  }
}
static void
read_free_full(struct read_entry * re)
{
  read_free(re);
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
  int found_alignments=0;
  if (h.load > 0) {
    found_alignments = read_pass2(re, &h);
    if (aligned_reads_file!=NULL) {
	#pragma omp critical (stdout) 
	{
		fasta_write_read(aligned_reads_file,re);
	}
    }
  }
  if (found_alignments==0) {
	if (Eflag && sam_unaligned) {
		//no alignments, print to sam empty record
		char* output1=NULL; 
      		hit_output(re, NULL, NULL, &output1, NULL, false, NULL,0);
		#pragma omp critical (stdout)
		{
			fprintf(stdout, "%s\n", output1);
		}
		free(output1);
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

  DEBUG("second pass");
  int found_alignments=0;
  if (h.load > 0) {
    found_alignments+=readpair_pass2(re1, re2, &h);
    if (aligned_reads_file!=NULL) {
	#pragma omp critical (stdout)
	{
		fasta_write_read(aligned_reads_file,re1);
		fasta_write_read(aligned_reads_file,re2);
	}
    }
  }
  if (found_alignments==0) {
	if (sam_half_paired) {
  		heap_unpaired h;
		unsigned int i;
		
  		read_pass1(re1, false);
  		heap_unpaired_init(&h, num_tmp_outputs);
  		read_get_vector_hits(re1, &h);
		if (h.load>0) {
			found_alignments+=read_pass2(re1, &h);	
		}
		// Done with this read; deallocate memory.
		for (i = 0; i < h.load; i++)
		  hit_free_sfrp(h.array[i].rest.hit);
		heap_unpaired_destroy(&h);
 
 		read_pass1(re2, false);
  		heap_unpaired_init(&h, num_tmp_outputs);
  		read_get_vector_hits(re2, &h);
		if (h.load>0) {
			found_alignments+=read_pass2(re2, &h);	
		}
		// Done with this read; deallocate memory.
		for (i = 0; i < h.load; i++)
		  hit_free_sfrp(h.array[i].rest.hit);
		heap_unpaired_destroy(&h);
	}
	if (Eflag && sam_unaligned && found_alignments==0) {
		//no alignments, print to sam empty record
		char* output1=NULL; char * output2=NULL;
      		hit_output(re1, NULL, NULL, &output1, NULL, true, NULL,0);
      		hit_output(re2, NULL, NULL, &output2, NULL, false, NULL,0);
		#pragma omp critical (stdout)
		{
			fprintf(stdout, "%s\n", output1);
			fprintf(stdout, "%s\n", output2);
		}
		free(output1);
		free(output2);		

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
launch_scan_threads(){
  fasta_t fasta=NULL;
  fasta_t left_fasta=NULL;
  fasta_t right_fasta=NULL;

  //open the fasta file and check for errors
  if (single_reads_file) {  
  	fasta = fasta_open(reads_filename, shrimp_mode, Qflag);
	if (fasta == NULL) {
		fprintf(stderr, "error: failed to open read file [%s]\n", reads_filename);
		return (false);
	} else {
		fprintf(stderr, "- Processing read file [%s]\n", reads_filename);
	}
  } else {
	left_fasta = fasta_open(left_reads_filename,shrimp_mode,Qflag);
	if (left_fasta == NULL) {
		fprintf(stderr, "error: failed to open read file [%s]\n", left_reads_filename);
		return (false);
	}
	right_fasta = fasta_open(right_reads_filename,shrimp_mode,Qflag);
	if (right_fasta == NULL) {
		fprintf(stderr, "error: failed to open read file [%s]\n", right_reads_filename);
		return (false);
	}
	if (right_fasta->space != left_fasta->space) {
		fprintf(stderr,"error: when using -1 and -2, both files must be either only colour space or only letter space!\n");
		return (false);
	}
	fasta=left_fasta;
	fprintf(stderr, "- Processing read files [%s , %s]\n", left_reads_filename,right_reads_filename);
  }

  if ((pair_mode != PAIR_NONE || !single_reads_file) && (chunk_size%2)!=0) {
	fprintf(stderr,"error: when in paired mode or using options -1 and -2, then thread chunk size must be even\n"); 
  }

  bool read_more = true;
  bool more_in_left_file=true; bool more_in_right_file=true;

#pragma omp parallel shared(read_more,more_in_left_file,more_in_right_file, fasta) num_threads(num_threads)
  {
    struct read_entry * re_buffer;
    int load, i;
    uint64_t before;

    re_buffer = (struct read_entry *)xcalloc(chunk_size * sizeof(re_buffer[0]));
    memset(re_buffer,0,chunk_size * sizeof(re_buffer[0]));

    while (read_more){
      before = rdtsc();
//Maybe just have one thread read in everything and then split the work ? TODO
#pragma omp critical (fill_reads_buffer)
      {
	wait_ticks[omp_get_thread_num()] += rdtsc() - before;

	load = 0;
	assert(chunk_size>2);
	while( read_more && ((single_reads_file && load < chunk_size) || (!single_reads_file && load < chunk_size-1))) {
	  //if (!fasta_get_next_with_range(fasta, &re_buffer[load].name, &re_buffer[load].seq, &re_buffer[load].is_rna,
	//				 &re_buffer[load].range_string, &re_buffer[load].qual))
	  if (single_reads_file) { 
		  if (fasta_get_next_read_with_range(fasta, &re_buffer[load])) {
		    load++;
		  } else { 
		    read_more = false;
		  }
	  } else {
		  //read from the left file
		  if (fasta_get_next_read_with_range(left_fasta, &re_buffer[load])) {
		  	load++;
		  } else {
		  	more_in_left_file = false;
		  }
		  //read from the right file
		  if (fasta_get_next_read_with_range(right_fasta, &re_buffer[load])) {
		  	load++;
		  } else {
		  	more_in_right_file = false;
		  }
		  //make sure that one is not smaller then the other
		  if (more_in_left_file!=more_in_right_file) {
			fprintf(stderr,"error: when using options -1 and -2, both files specified must have the same number of entries\n");
			exit(1);
		  }
		  //keep reading?
		  read_more=more_in_left_file && more_in_right_file; 
	  }
	}
	nreads += load;
      } // end critical section
      if (pair_mode != PAIR_NONE)
	assert(load % 2 == 0); // read even number of reads

      for (i = 0; i < load; i++) {
	// if running in paired mode and first foot is ignored, ignore this one, too
	if (pair_mode != PAIR_NONE && i % 2 == 1 && re_buffer[i-1].ignore) {
	  read_free(re_buffer+i-1);
	  read_free(re_buffer+i);
	  continue;
	}

	//if (!(strcspn(re_buffer[i].seq, "nNxX.") == strlen(re_buffer[i].seq))) {
	//  if (pair_mode != PAIR_NONE && i % 2 == 1) {
	//    read_free(re_buffer+i-1);
	//    read_free(re_buffer+i);
	//  }
	//  continue;
	//}

	re_buffer[i].ignore = false;
	re_buffer[i].read[0] = fasta_sequence_to_bitfield(fasta, re_buffer[i].seq);
	re_buffer[i].read_len = strlen(re_buffer[i].seq);
	re_buffer[i].max_n_kmers = re_buffer[i].read_len - min_seed_span + 1;
	re_buffer[i].min_kmer_pos = 0;
	if (shrimp_mode == MODE_COLOUR_SPACE){
	  re_buffer[i].read_len--;
	  re_buffer[i].max_n_kmers -= 2; // 1st color always discarded from kmers
	  re_buffer[i].min_kmer_pos = 1;
	  re_buffer[i].initbp[0] = fasta_get_initial_base(shrimp_mode,re_buffer[i].seq);
	  re_buffer[i].initbp[1] = re_buffer[i].initbp[0];
	  re_buffer[i].read[1] = reverse_complement_read_cs(re_buffer[i].read[0], (int8_t)re_buffer[i].initbp[0], (int8_t)re_buffer[i].initbp[1],
							    re_buffer[i].read_len, re_buffer[i].is_rna);
	} else {
	  re_buffer[i].read[1] = reverse_complement_read_ls(re_buffer[i].read[0], re_buffer[i].read_len, re_buffer[i].is_rna);
	}
	//Check if we can actually use this read
	if (re_buffer[i].max_n_kmers<0) {
		fprintf(stderr,"warning: Read smaller then any seed, skipping read!\n");
		read_free_full(&re_buffer[i]);
		continue;	
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

	//free(re_buffer[i].seq);
	if (pair_mode == PAIR_NONE)
	  {
	    if (re_buffer[i].read_len>longest_read_len) {
		fprintf(stderr,"warning: read \"%s\" has length %d, maximum length allowed is %d. Use --longest-read ?\n",re_buffer[i].name,re_buffer[i].read_len,longest_read_len);
	    } else {
	    	handle_read(&re_buffer[i]);
	    }
	  }
	else if (i % 2 == 1) 
	  {
	    if (re_buffer[i].read_len>longest_read_len) {
		fprintf(stderr,"warning: read \"%s\" has length %d, maximum length allowed is %d. Use --longest-read ?\n",re_buffer[i].name,re_buffer[i].read_len,longest_read_len);
	    } else if (re_buffer[i-1].read_len>longest_read_len) {
		fprintf(stderr,"warning: read \"%s\" has length %d, maximum length allowed is %d. Use --longest-read ?\n",re_buffer[i-1].name,re_buffer[i-1].read_len,longest_read_len);
	    } else {
	
		    if (pair_reverse[pair_mode][0])
		      read_reverse(&re_buffer[i-1]);
		    if (pair_reverse[pair_mode][1])
		      read_reverse(&re_buffer[i]);
		    re_buffer[i-1].paired=true;
		    re_buffer[i-1].first_in_pair=true;
		    re_buffer[i-1].mate_pair=&re_buffer[i];
		    re_buffer[i].paired=true;
		    re_buffer[i].first_in_pair=false;
		    re_buffer[i].mate_pair=&re_buffer[i-1];
		    handle_readpair(&re_buffer[i-1], &re_buffer[i]);
	   }
	  }
      }
    }

    free(re_buffer);
  } // end parallel section

  if (single_reads_file) {
    fasta_close(fasta);
  } else {
    fasta_close(left_fasta);
    fasta_close(right_fasta);
  }
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


void free_genome(void) {
	int sn,capacity;
	for (sn=0; sn<n_seeds; sn++){
	  	if(Hflag){
	 	   capacity = power(4, HASH_TABLE_POWER);
	 	 } else {
	    	capacity = power(4, seed[sn].weight);
	  	}
		//uint32_t mapidx = kmer_to_mapidx(kmerWindow, sn);
		int i;
		for(i=0; i<capacity; i++) {
			if (genomemap[sn][i]!=NULL && (i==0 || load_file==NULL)) {
				free(genomemap[sn][i]);
			}
		}	
		free(genomemap[sn]);
		free(genomemap_len[sn]);
	}
	free(genomemap);
	free(genomemap_len);
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
	  memset(genomemap[sn],0,sizeof(uint32_t *) * capacity);
	  genomemap_len[sn] = (uint32_t *)xcalloc_c(sizeof(uint32_t) * capacity, &mem_genomemap);

	}
	num_contigs = 0;
	u_int i = 0;
	int cfile;
	for(cfile = 0; cfile < nfiles; cfile++){
		file = files[cfile];
		//open the fasta file and check for errors
		fasta = fasta_open(file, MODE_LETTER_SPACE, false);
		if (fasta == NULL) {
			fprintf(stderr,"error: failded to open genome file [%s]\n",file);
			return (false);
		} else {
			fprintf(stderr,"- Processing genome file [%s]\n",file);
		}

		//Read the contigs and record their sizes

		while(fasta_get_next_contig(fasta, &name, &seq, &is_rna)){
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
			if (shrimp_mode==MODE_COLOUR_SPACE) {
				free(read);
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
	    fprintf(stderr, "%s%s%-24s" "%s    (%.4f%%)\n", my_tab, my_tab,
		    "Reads Dropped:",
		    comma_integer(total_reads_dropped),
		    (nreads == 0) ? 0 : ((double)total_reads_dropped / (double)nreads) * 100);
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
	    fprintf(stderr, "%s%s%-40s" "%s    (%.4f%%)\n", my_tab, my_tab,
		    "Pairs Dropped:",
		    comma_integer(total_pairs_dropped),
		    (nreads == 0) ? 0 : ((double)total_pairs_dropped / (double)(nreads/2)) * 100);
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

  if (n_seeds == 0)
    load_default_seeds(0);

  slash = strrchr(progname, '/');
  if (slash != NULL)
    progname = slash + 1;

  fprintf(stderr, 
	  "usage: %s [options/parameters] { <r> | -1 <r1> -2 <r2> } <g1> <g2>...\n", progname);
  fprintf(stderr,
	  "   <r>                  Reads filename, paired or unpaired\n");
  fprintf(stderr,
	  "   <r1>                 Upstream reads filename\n");
  fprintf(stderr,
	  "   <r2>                 Downstream reads filename\n");
  fprintf(stderr,
	  "   <g1> <g2>...         Space seperated list of genome filenames\n");
  fprintf(stderr,
	  "Parameters:\n");
  fprintf(stderr,
	  "   -s/--seeds           Spaced Seed(s)                (default: ");
  for (sn = 0; sn < n_seeds; sn++) {
    if (sn > 0)
      fprintf(stderr,
	  "                                                       ");
    fprintf(stderr, "%s%s\n", seed_to_string(sn), (sn == n_seeds - 1? ")" : ","));
  }
  fprintf(stderr,
	  "   -o/--report          Maximum Hits per Read         (default: %d)\n",
	  DEF_NUM_OUTPUTS);
  fprintf(stderr,
          "      --max-alignments  Max. align. per read  (0=all) (default: %d)\n",
	  DEF_MAX_ALIGNMENTS);
  fprintf(stderr,
	  "   -w/--match-window    Match Window Length           (default: %.02f%%)\n",
	  DEF_WINDOW_LEN);
  fprintf(stderr,
	  "   -n/--cmw-mode        Seed Matches per Window       (default: %d)\n",
	  DEF_NUM_MATCHES);
  if (full_usage) {
  fprintf(stderr,
	  "   -l/--cmw-overlap     Match Window Overlap Length   (default: %.02f%%)\n",
	  DEF_WINDOW_OVERLAP);
  fprintf(stderr,
	  "   -a/--anchor-width    Anchor Width Limiting Full SW (default: %d; disable: -1)\n",
	  DEF_ANCHOR_WIDTH);

  fprintf(stderr, "\n");
  fprintf(stderr,
	  "   -S/--save            Save Genome Proj. in File     (default: no)\n");
  fprintf(stderr,
	  "   -L/--load            Load Genome Proj. from File   (default: no)\n");
  fprintf(stderr,
	  "   -z/--cutoff          Projection List Cut-off Len.  (default: %u)\n",
	  DEF_LIST_CUTOFF);
  }

  fprintf(stderr, "\n");
  fprintf(stderr,
	  "   -m/--match           SW Match Score                (default: %d)\n",
	  DEF_MATCH_VALUE);
  fprintf(stderr,
	  "   -i/--mismatch        SW Mismatch Score             (default: %d)\n",
	  DEF_MISMATCH_VALUE);
  fprintf(stderr,
	  "   -g/--open-r          SW Gap Open Score (Reference) (default: %d)\n",
	  DEF_A_GAP_OPEN);
  fprintf(stderr,
	  "   -q/--open-q          SW Gap Open Score (Query)     (default: %d)\n",
	  DEF_B_GAP_OPEN);
  fprintf(stderr,
	  "   -e/--ext-r           SW Gap Extend Score(Reference)(default: %d)\n",
	  DEF_A_GAP_EXTEND);
  fprintf(stderr,
	  "   -f/--ext-q           SW Gap Extend Score (Query)   (default: %d)\n",
	  DEF_B_GAP_EXTEND);
  if (shrimp_mode == MODE_COLOUR_SPACE) {
  fprintf(stderr,
	  "   -x/--crossover       SW Crossover Score            (default: %d)\n",
	  DEF_XOVER_PENALTY);
  }
  fprintf(stderr,
	  "   -r/--cmw-threshold   Window Generation Threshold   (default: %.02f%%)\n",
	  DEF_WINDOW_GEN_THRESHOLD);
  if (shrimp_mode == MODE_COLOUR_SPACE) {
  fprintf(stderr,
	  "   -v/--vec-threshold   SW Vector Hit Threshold       (default: %.02f%%)\n",
	  DEF_SW_VECT_THRESHOLD);
  }
  fprintf(stderr,
	  "   -h/--full-threshold  SW Full Hit Threshold         (default: %.02f%%)\n",
	  DEF_SW_FULL_THRESHOLD);

  fprintf(stderr, "\n");

  fprintf(stderr,
	  "   -N/--threads         Number of Threads             (default: %d)\n",
	  DEF_NUM_THREADS);
  if (full_usage) {
  fprintf(stderr,
	  "   -K/--thread-chunk    Thread Chunk Size             (default: %d)\n",
	  DEF_CHUNK_SIZE);
  }

  fprintf(stderr, "\n");
  fprintf(stderr,
	  "   -p/--pair-mode       Paired Mode                   (default: %s)\n",
	  pair_mode_string[pair_mode]);
  fprintf(stderr,
	  "   -I/--isize           Min and Max Insert Size       (default: %d,%d)\n",
	  DEF_MIN_INSERT_SIZE, DEF_MAX_INSERT_SIZE);
  fprintf(stderr,
	  "      --longest-read    Maximum read length           (default: %d)\n",
	  DEF_LONGEST_READ_LENGTH);
  fprintf(stderr,
	  "   -1/--upstream        Upstream read pair file\n");
  fprintf(stderr,
	  "   -2/--downstream      Downstream read pair file\n");
  fprintf(stderr,
	  "      --un              Dump unaligned reads to file\n");
  fprintf(stderr,
	  "      --al              Dump aligned reads to file\n");
  fprintf(stderr,
	  "      --read-group      Attach SAM Read Group name\n");
  fprintf(stderr,
          "      --sam-header      Use file as SAM header\n");

  fprintf(stderr, "\n");
  fprintf(stderr, "Options:\n");

  fprintf(stderr,
	  "   -U/--ungapped        Perform Ungapped Alignment    (default: disabled)\n");
  fprintf(stderr,
          "      --global          Perform full global alignment (default: disabled)\n");
  if (shrimp_mode == MODE_COLOUR_SPACE) {
  fprintf(stderr,
          "      --bfast           Try to align like bfast       (default: disabled)\n");
  }
  fprintf(stderr,
	  "   -C/--negative        Negative Strand Aln. Only     (default: disabled)\n");
  fprintf(stderr,
	  "   -F/--positive        Positive Strand Aln. Only     (default: disabled)\n");
  fprintf(stderr,
	  "   -P/--pretty          Pretty Print Alignments       (default: disabled)\n");
  fprintf(stderr,
	  "   -E/--sam             Output SAM Format             (default: disabled)\n");
  fprintf(stderr,
  	  "   -Q/--fastq           Reads are in fastq format     (default: disabled)\n");
  if (full_usage) {
  fprintf(stderr,
	  "   -R/--print-reads     Print Reads in Output         (default: disabled)\n");
 // fprintf(stderr,
//	  "    -T    (does nothing since default) Reverse Tie-break on Negative Strand          (default: enabled)\n");
  fprintf(stderr,
	  "   -t/--tiebreak-off    Disable Reverse Tie-break\n");
  fprintf(stderr,
	  "                                  on Negative Strand  (default: enabled)\n");
  fprintf(stderr,
	  "   -X/--isize-hist      Print Insert Size Histogram   (default: disabled)\n");
  fprintf(stderr,
	  "   -Y/--proj-hist       Print Genome Proj. Histogram  (default: disabled)\n");
  fprintf(stderr,
	  "   -Z/--bypass-off      Disable Cache Bypass for SW\n");
  fprintf(stderr,
	  "                                    Vector Calls      (default: enabled)\n");
  fprintf(stderr,
	  "   -H/--spaced-kmers    Hash Spaced Kmers in Genome\n");
  fprintf(stderr,
	  "                                    Projection        (default: disabled)\n");
  fprintf(stderr,
	  "   -D/--thread-stats    Individual Thread Statistics  (default: disabled)\n");
  fprintf(stderr,
	  "   -V/--trim-off        Disable Automatic Genome\n");
  fprintf(stderr,
	  "                                 Index Trimming       (default: enabled)\n");
  }
  fprintf(stderr,
	  "      --sam-unaligned   Unaligned reads in SAM output (default: disabled)\n");
  fprintf(stderr,
	  "      --half-paired     Output half mapped read pairs (default: disabled)\n");
  fprintf(stderr,
	  "      --strata          Print only the best scoring hits\n");
  fprintf(stderr,
	  "   -?/--help            Full List of Parameters and Options\n");

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


static int
set_mode_from_string(char const * s) {
  if (!strcmp(s, "mirna")) {
    mode_mirna = true;

    //load_default_mirna_seeds();

    Hflag = true;
    gapless_sw = true;
    anchor_width = 0;
    a_gap_open_score = -255;
    b_gap_open_score = -255;
    hash_filter_calls = false;
    num_matches = 1;
    window_len = 100.0;

    return 1;
  } else {
    return 0;
  }
}


int main(int argc, char **argv){
	char **genome_files = NULL;
	int ngenome_files = 0;

	char *progname = argv[0];
	char const * optstr = NULL;
	char *c;
	int ch;

	bool a_gap_open_set, b_gap_open_set;
	bool a_gap_extend_set, b_gap_extend_set;
	bool num_matches_set = false;

	shrimp_args.argc=argc;
	shrimp_args.argv=argv;
	set_mode_from_argv(argv);

	a_gap_open_set = b_gap_open_set = a_gap_extend_set = b_gap_extend_set = false;

	fprintf(stderr, "--------------------------------------------------"
			"------------------------------\n");
	fprintf(stderr, "gmapper: %s.\nSHRiMP %s\n[%s]\n", get_mode_string(),
			SHRIMP_VERSION_STRING, get_compiler());
	fprintf(stderr, "--------------------------------------------------"
			"------------------------------\n");

	struct option getopt_long_options[standard_entries+MAX(colour_entries,letter_entries)];
	memcpy(getopt_long_options,standard_options,sizeof(standard_options));
	//TODO -t -9 -d -Z -D -Y
	switch(shrimp_mode){
	case MODE_COLOUR_SPACE:
		optstr = "?1:2:s:n:w:l:o:p:m:i:g:q:e:f:h:r:a:z:DCEFHI:K:L:M:N:PRS:TtUVXYZQx:v:";
		memcpy(getopt_long_options+standard_entries,colour_space_options,sizeof(colour_space_options));
		break;
	case MODE_LETTER_SPACE:
		optstr = "?1:2:s:n:w:l:o:p:m:i:g:q:e:f:h:r:a:z:DCEFHI:K:L:M:N:PRS:TtUVXYZQ";
		memcpy(getopt_long_options+standard_entries,letter_space_options,sizeof(letter_space_options));
		break;
	case MODE_HELICOS_SPACE:
		fprintf(stderr,"Helicose currently unsupported\n");
		exit(1);
		break;
	}

	

	while ((ch = getopt_long(argc,argv,optstr,getopt_long_options,NULL)) != -1){
		switch (ch) {
		case 9:
			strata_flag = true;
			break;
		case 10:
			unaligned_reads_file=fopen(optarg,"w");
			if (unaligned_reads_file==NULL) {
				fprintf(stderr,"error: cannot open file \"%s\" for writting\n",optarg);	
			}
			break;
		case 11:
			aligned_reads_file=fopen(optarg,"w");
			if (aligned_reads_file==NULL) {
				fprintf(stderr,"error: cannot open file \"%s\" for writting\n",optarg);	
			}
			break;
		case 12:
			sam_unaligned=true;
			break;
		case 13:
			longest_read_len=atoi(optarg);
			if (longest_read_len<200) {
				fprintf(stderr,"error: longest read length must be at least 200\n");
				exit(1);
			}
			break;
		case 14:
			max_alignments=atoi(optarg);
			break;
		case 15:
			Gflag = true;
			break;
		case 16:
			Bflag = true;
			Gflag = true;
			break;
		case 17:
			sam_read_group_name=optarg;
			sam_sample_name = strchr(optarg,',');
			if (sam_sample_name==NULL) {
				fprintf(stderr,"error: sam read group needs to be two values, delimited by commas.\n");
				fprintf(stderr,"       the first value is unique read group identifier\n");
				fprintf(stderr,"       the second is the sample (use pool name where a pool is being sequence\n");
				exit(1);	
			}
			sam_sample_name[0]='\0';
			sam_sample_name++;
			break;
		case 18:
			sam_header_filename=optarg;
			break;	
		case 19:
			sam_half_paired=true;
			break;
		case 20:
			sam_r2=true;
			break;
		case '1':
			left_reads_filename = optarg;
			break;
		case '2':
			right_reads_filename = optarg;
			break;
		case 's':
			if (strchr(optarg, ',') == NULL) { // allow comma-separated seeds
				if (optarg[0] == 'w') {
					int weight = (int)atoi(&optarg[1]);
					if (!load_default_seeds(weight)) {
						fprintf(stderr, "error: invalid spaced seed weight (%d)\n", weight);
						exit(1);
					}
				} else {
					if (!add_spaced_seed(optarg)) {
						fprintf(stderr, "error: invalid spaced seed \"%s\"\n", optarg);
						exit (1);
					}
				}
			} else {
				c = strtok(optarg, ",");
				do {
                                	if (c[0] == 'w') {
	                                        int weight = (int)atoi(&c[1]);
        	                                if (!load_default_seeds(weight)) {
                	                                fprintf(stderr, "error: invalid spaced seed weight (%d)\n", weight);
                        	                        exit(1);
                                	        }
	                                } else {
        	                                if (!add_spaced_seed(c)) {
                	                                fprintf(stderr, "error: invalid spaced seed \"%s\"\n", c);
                        	                        exit (1);
                                	        }
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
			num_tmp_outputs = 30 + num_outputs;
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
		case 't':
			Tflag= false;
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
		case 'Q':
			Qflag = true;
			break;
		case 'M':
			c = strtok(optarg, ",");
			do {
			  if (!set_mode_from_string(c)) {
			    fprintf(stderr, "error: unrecognized mode (%s)\n", c);
			    exit(1);
			  }
			  c = strtok(NULL, ",");
			} while (c != NULL);		  
			break;
		default:
			usage(progname, false);
		}
	}

	argc -= optind;
	argv += optind;

	if (Gflag && gapless_sw) {
		fprintf(stderr,"error: cannot use global (or bfast) and ungapped mode at the same time!\n");
		exit(1);
	}
	if (sam_unaligned && !Eflag) {
		fprintf(stderr,"error: when using flag --sam-unaligned must also use -E/--sam\n");
		usage(progname,false);
	}
	if (right_reads_filename != NULL || left_reads_filename !=NULL) {
		if (right_reads_filename == NULL || left_reads_filename == NULL ){
			fprintf(stderr,"error: when using \"%s\" must also specify \"%s\"\n",
				(left_reads_filename != NULL) ? "-1" : "-2",
				(left_reads_filename != NULL) ? "-2" : "-1");
			usage(progname,false);
		}
		single_reads_file=false;
		if (strcmp(right_reads_filename,"-")==0 && strcmp(left_reads_filename,"-")==0) {
			fprintf(stderr,"error: both -1 and -2 arguments cannot be stdin (\"-\")\n");
			usage(progname,false);
		}
	}
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
	  if (mode_mirna)
	    load_default_mirna_seeds();
	  else
	    load_default_seeds(0);
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
	    if (argc == 0 && single_reads_file) {
	      fprintf(stderr,"error: read_file not specified\n");
	      usage(progname, false);
	    } else if (argc == 1) {
	      if (single_reads_file) {
	      	reads_filename    = argv[0];
	      } else {
		fprintf(stderr,"error: cannot specify a reads file when using -L, -1 and -2\n");
		usage(progname,false);
	      }
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
	else if (single_reads_file)
	  { // args: reads file, genome file(s)
	    if (argc < 2) {
	      fprintf(stderr, "error: %sgenome_file(s) not specified\n",
		      (argc == 0) ? "reads_file, " : "");
	      usage(progname, false);
	    }
	    reads_filename    = argv[0];
	    genome_files  = &argv[1];
	    ngenome_files = argc - 1;
	  }
	else 
	  {
	   if( argc < 1) {
	      fprintf(stderr, "error: genome_file(s) not specified\n");
	      usage(progname, false);
	   }
	    genome_files  = &argv[0];
	    ngenome_files = argc;
	  }

	if (!Cflag && !Fflag) {
	  Cflag = Fflag = true;
	}

	if (pair_mode != PAIR_NONE && (!Cflag || !Fflag)) {
	  fprintf(stderr, "warning: in paired mode, both strands must be inspected; ignoring -C and -F\n");
	  Cflag = Fflag = true;
	}
	if (pair_mode == PAIR_NONE && sam_half_paired) {
	  fprintf(stderr, "error: cannot use option half-paired in non-paired mode!\n");
	  exit(1);
	}
	if (pair_mode == PAIR_NONE && sam_r2) {
	  fprintf(stderr, "error: cannot use option sam-r2 in non-paired mode!\n");
	  exit(1);
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
	//int longest_read_len = 2000;
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
		int i;
		if (sam_header_filename!=NULL) {
			FILE * sam_header_file = fopen(sam_header_filename,"r");
			if (sam_header_file==NULL) {
				perror("Failed to open sam header file ");
				exit(1);
			}
			size_t buffer_size=2046;
			char buffer[buffer_size];
			size_t read; bool ends_in_newline=true;
			while ((read=fread(buffer,1,buffer_size-1,sam_header_file))) {
				buffer[read]='\0';
				fprintf(stdout,buffer);
				if (buffer[read-1]=='\n') {
					ends_in_newline=true;
				} else {
					ends_in_newline=false;
				}
			}
			if (!ends_in_newline) {
				fprintf(stdout,"\n");
			}
		} else {
			//Print sam header
			fprintf(stdout,"@HD\tVN:%s\tSO:%s\n","1.0","unsorted");

			for(i = 0; i < num_contigs; i++){
				fprintf(stdout,"@SQ\tSN:%s\tLN:%u\n",contig_names[i],genome_len[i]);
			}
		}
		//read group
		if (sam_read_group_name!=NULL) {
			fprintf(stdout,"@RG\tID:%s\tSM:%s\n",sam_read_group_name,sam_sample_name);
		}
		//print command line args used to invoke SHRiMP
		fprintf(stdout,"@PG\tID:%s\tVN:%s\tCL:","gmapper",SHRIMP_VERSION_STRING);
		for (i=0; i<(shrimp_args.argc-1); i++) {
			fprintf(stdout,"%s ",shrimp_args.argv[i]);
		}
		fprintf(stdout,"%s\n",shrimp_args.argv[i]);
	} else {
		output = output_format_line(Rflag);
		puts(output);
		free(output);
	}
	before = gettimeinusecs();
	bool launched = launch_scan_threads();
	if (!launched) {
		fprintf(stderr,"error: a fatal error occured while launching scan thread(s)!\n");
		exit(1);
	}
	total_work_usecs += (gettimeinusecs() - before);
	
	//if (load_file==NULL) {
		free_genome();
	//}
	print_statistics();
#pragma omp parallel shared(longest_read_len,max_window_len,a_gap_open_score, a_gap_extend_score, b_gap_open_score, b_gap_extend_score,\
		match_score, mismatch_score,shrimp_mode,crossover_score,anchor_width) num_threads(num_threads)
	{
		sw_vector_cleanup();
		if (shrimp_mode==MODE_COLOUR_SPACE) {
			sw_full_cs_cleanup();
		}
		sw_full_ls_cleanup();
		f1_free();	
	}
	int i;
	for (i=0; i<num_contigs; i++){
		free(contig_names[i]);
		if (i==0 || load_file==NULL) {
			free(genome_contigs[i]);
			free(genome_contigs_rc[i]);
		}
	}
	if (shrimp_mode==MODE_COLOUR_SPACE) {
		for (i=0; i<num_contigs && (i==0 || load_file==NULL); i++){
			free(genome_cs_contigs[i]);
		}
		free(genome_cs_contigs);
		free(genome_initbp);
	}
	free(genome_len);
	free(genome_contigs_rc);
	free(genome_contigs);
	free(contig_names);
	free(contig_offsets);
	free(seed);
}
