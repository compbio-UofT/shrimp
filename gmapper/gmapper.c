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
#include <omp.h>

#include <xmmintrin.h>	// for _mm_prefetch

#include <sys/types.h>
#include <sys/stat.h>
#include <sys/time.h>

#include "../common/hash.h"
#include "../common/fasta.h"
#include "../common/util.h"
#include "../gmapper/gmapper.h"
#include "../mapper/mapper.h" // why? to reduce code duplication contains common functions from gmapper and rmapper
#include "../common/version.h"

#include "../common/sw-full-common.h"
#include "../common/sw-full-cs.h"
#include "../common/sw-full-ls.h"
#include "../common/sw-vector.h"
#include "../common/output.h"

/* Parameters */
static double	window_len		= DEF_WINDOW_LEN;
static double	window_overlap		= DEF_WINDOW_OVERLAP;
static uint	num_matches		= DEF_NUM_MATCHES;
static uint	num_outputs		= DEF_NUM_OUTPUTS;
static double	sw_vect_threshold	= DEF_SW_VECT_THRESHOLD;
static double	sw_full_threshold	= DEF_SW_FULL_THRESHOLD;

/* Scores */
static int	match_score		= DEF_MATCH_VALUE;
static int	mismatch_score		= DEF_MISMATCH_VALUE;
static int	a_gap_open_score	= DEF_A_GAP_OPEN;
static int	a_gap_extend_score	= DEF_A_GAP_EXTEND;
static int	b_gap_open_score	= DEF_B_GAP_OPEN;
static int	b_gap_extend_score	= DEF_B_GAP_EXTEND;
static int	crossover_score		= DEF_XOVER_PENALTY;

/* Flags */
//static int Bflag = false;			/* print a progress bar */
//static int Cflag = false;			/* do complement only */
//static int Fflag = false;			/* do positive (forward) only */
//static int Hflag = false;			/* use hash table, not lookup */
static int Pflag = true;			/* pretty print results */
static int Rflag = false;			/* add read sequence to output*/
static int Tflag = false;			/* reverse sw full tie breaks */
//static int Uflag = false;			/* output unmapped reads, too */

/* Statistics */
static uint	nduphits;			/* number of duplicate hits */


/* Kmer to genome index */
static uint32_t ***genomemap;
static uint32_t **genomemap_len;

/* offset info for genome contigs */
static uint32_t *contig_offsets;
static char **contig_names = NULL;
static uint32_t num_contigs;

/* Genomic sequence, stored in 32-bit words, first is in the LSB */
static uint32_t **genome_contigs;			/* genome -- always in letter */
static uint32_t **genome_contigs_rc;			/* reverse complemets */
static uint32_t **genome_cs_contigs;
static uint32_t  *genome_initbp;
static uint32_t	 *genome_len;

static bool      genome_is_rna = false;		/* is genome RNA (has uracil)?*/

static count_t mem_genomemap;

static uint numThreads = omp_get_num_procs();
static uint chunkSize = 1000;

extern size_t
power(size_t base, size_t exp);

extern uint32_t
kmer_to_mapidx_hash(uint32_t *kmerWindow, u_int sn);

/* If x is negative, return its absolute value; else return base*x% */
static inline double
abs_or_pct(double x, double base) {
	if (IS_ABSOLUTE(x))
		return -x;
	else
		return base * (x / 100.0);
}

/* percolate down in our min-heap */
static void
reheap(struct re_score *scores, uint node)
{
  struct re_score tmp;
  uint left, right, max;

  assert(node >= 1 && node <= (uint)scores[0].heap_capacity);

  left  = node * 2;
  right = left + 1;
  max   = node;
  
  if (left <= scores[0].heap_elems &&
      scores[left].score < scores[node].score)
    max = left;

  if (right <= scores[0].heap_elems &&
      scores[right].score < scores[max].score)
    max = right;

  if (max != node) {
    tmp = scores[node];
    scores[node] = scores[max];
    scores[max] = tmp;
    reheap(scores, max);
  }
}

/* percolate up in our min-heap */
static void
percolate_up(struct re_score *scores, uint node)
{
  struct re_score tmp;
  int parent;

  assert(node >= 1 && node <= (uint)scores[0].heap_capacity);

  if (node == 1)
    return;

  parent = node / 2;

  if (scores[parent].score > scores[node].score) {
    tmp = scores[node];
    scores[node] = scores[parent];
    scores[parent] = tmp;
    percolate_up(scores, parent);
  }
}

/* update our score min-heap */
static void
save_score(read_entry * re, int score, uint g_idx, int contig_num, bool rev_cmpl)
{
  struct re_score * scores = re->scores;

  if (scores[0].heap_elems == num_outputs) {
    if (score < scores[1].score)
      return;

    /* replace min score */
    scores[1].score = score;
    scores[1].g_idx = g_idx;
    scores[1].contig_num = contig_num;
    scores[1].rev_cmpl = rev_cmpl;

    reheap(scores, 1);
  } else {
    uint idx = 1 + scores[0].heap_elems++;

    assert(idx <= scores[0].heap_capacity);

    scores[idx].score = score;
    scores[idx].g_idx = g_idx;
    scores[idx].contig_num = contig_num;
    scores[idx].rev_cmpl = rev_cmpl;
    
    percolate_up(scores, idx);
  }
}

static int
score_cmp(const void *arg1, const void *arg2)
{
  const struct re_score *one = (const struct re_score *)arg1;
  const struct re_score *two = (const struct re_score *)arg2;

  //assert(one->revcmpl == two->revcmpl);

  if (one->score > two->score)
    return (-1);
  if (one->score < two->score)
    return (1);

  if (one->sfrp->genome_start > two->sfrp->genome_start)
    return (1);
  if (one->sfrp->genome_start < two->sfrp->genome_start)
    return (-1);

  if (one->sfrp->matches > two->sfrp->matches)
    return (-1);
  if (one->sfrp->matches < two->sfrp->matches)
    return (1);

  return (0);
}

void print_info(){
	fprintf(stderr,"num_contigs = %u\n",num_contigs);
	uint i;
	for (i = 0; i < num_contigs; i++){
		fprintf(stderr,"contig %u: name=%s offset = %u\n",i,contig_names[i],contig_offsets[i]);
	}
	fprintf(stderr,"n_seeds = %u\n",n_seeds);
	for (i=0; i < n_seeds;i++){
		fprintf(stderr,"seed %u: %s\n",i,seed_to_string(i));
	}
#ifdef DEBUGGING
	for(i=0; i < n_seeds;i++){
		fprintf(stderr,"map for seed %u\n",i);
		uint j;
		for(j = 0; j< power(4, seed[i].weight);j++){
			//fprintf(stderr,"%u\n",genomemap_len[i][j]);
			if(genomemap_len[i][j] != 0){
				fprintf(stderr,"entry at %u\n",j);
				uint k;
				for(k=0;k<genomemap_len[i][j];k++){
					fprintf(stderr,"%u, ",genomemap[i][j][k]);
				}
				fprintf(stderr,"\n");
			}

		}
	}
#endif

}

static bool save_genome_map(const char *file) {
	gzFile fp = gzopen(file, "wb");
	if (fp == NULL){
		return false;
	}

	//write the header
	//fprintf(stderr,"saving num_contgs: %u\n",num_contigs);
	uint32_t i;
	gzwrite(fp,&num_contigs,sizeof(uint32_t));
	//fprintf(stderr,"saved num_contigs\n");
	for (i = 0; i < num_contigs; i++) {
		gzwrite(fp,&contig_offsets[i],sizeof(uint32_t));
		uint32_t len = strlen(contig_names[i]);
		gzwrite(fp,&len,sizeof(uint32_t));
		gzwrite(fp,contig_names[i],len +1);
	}
	//fprintf(stderr,"saving seeds\n");
	//write the seeds and genome_maps
	gzwrite(fp,&n_seeds,sizeof(uint32_t));
	for (i = 0; i < n_seeds; i++) {
		gzwrite(fp,&seed[i], sizeof(seed_type));

		//write the genome_map for this seed
		uint32_t j;
		uint32_t p = power(4, seed[i].weight);
		//fprintf(stderr,"saving index\n");
		gzwrite(fp,&p, sizeof(uint32_t));
		for (j = 0; j < p; j++) {
			uint32_t len = genomemap_len[i][j];
			gzwrite(fp, &len, sizeof(uint32_t));
			gzwrite(fp, genomemap[i][j], sizeof(uint32_t) * len);
		}
	}
	gzclose(fp);
	return true;
}

static bool load_genome_map(const char *file){
	gzFile fp = gzopen(file,"rb");
	if (fp == NULL){
		return false;
	}

	uint32_t i;
	gzread(fp,&num_contigs,sizeof(uint32_t));
	//fprintf(stderr,"num_contigs = %u\n",num_contigs);

	contig_names = (char **)xrealloc(contig_names,sizeof(char *)*num_contigs);
	contig_offsets = (uint32_t *)xrealloc(contig_offsets,sizeof(uint32_t)*num_contigs);

	for (i = 0; i < num_contigs; i++){
		gzread(fp,&contig_offsets[i],sizeof(uint32_t));

		uint32_t len;
		gzread(fp,&len,sizeof(uint32_t));
		contig_names[i] = (char *)xrealloc(contig_names[i],sizeof(char)*len);
		gzread(fp,contig_names[i],len+1);
		//fprintf(stderr,"contig %u: name=%s offset = %u\n",i,contig_names[i],contig_offsets[i]);
	}
	gzread(fp,&n_seeds,sizeof(uint32_t));
	//fprintf(stderr,"n_seeds = %u\n",n_seeds);
	seed =(seed_type *)xrealloc(seed,sizeof(seed_type)*n_seeds);
	genomemap_len = (uint32_t **)xrealloc(genomemap_len,sizeof(uint32_t *)*n_seeds);
	genomemap = (uint32_t ***)xrealloc(genomemap,sizeof(uint32_t **) * n_seeds);
	for (i = 0; i < n_seeds; i++){
		gzread(fp,&seed[i], sizeof(seed_type));
		//fprintf(stderr,"seed %u: span=%u\n",i,seed[i].span);
		uint32_t j;
		uint32_t p;
		gzread(fp,&p,sizeof(uint32_t));
		genomemap_len[i] = (uint32_t *)xrealloc(genomemap_len[i],sizeof(uint32_t)*p);
		genomemap[i] = (uint32_t **)xrealloc(genomemap[i],sizeof(uint32_t *)*p);
		for (j = 0; j < p; j++){
			gzread(fp,&genomemap_len[i][j],sizeof(uint32_t));
			genomemap[i][j] = (uint32_t *)xrealloc(genomemap[i][j],
					sizeof(uint32_t)*genomemap_len[i][j]);
			gzread(fp,genomemap[i][j],sizeof(uint32_t) * genomemap_len[i][j]);
		}
	}
	gzclose(fp);
	return true;
}

static int comp(const void *a, const void *b){
	return *(uint32_t *)a - *(uint32_t *)b;
}

static void
read_locations(read_entry *re, int *len,int *len_rc,uint32_t **list, uint32_t **list_rc){
	*len = 0;
	*len_rc =0;
	uint sn,i;
	uint load;
	uint32_t *kmerWindow = (uint32_t *)xcalloc(sizeof(kmerWindow[0])*BPTO32BW(max_seed_span));
	uint32_t *kmerWindow_rc = (uint32_t *)xcalloc(sizeof(kmerWindow_rc[0])*BPTO32BW(max_seed_span));
	for (load = 0, i = 0 ; i < re->read_len; i++){
		uint base;
		base = EXTRACT(re->read, i);
		bitfield_prepend(kmerWindow, max_seed_span, base);
		base = EXTRACT(re->read_rc,i);
		bitfield_prepend(kmerWindow_rc, max_seed_span, base);
		//skip past any Ns or Xs
		if (base == BASE_N || base == BASE_X)
			load = 0;
		else if (load < max_seed_span)
			load++;
		for (sn = 0; sn < n_seeds; sn++){
			if (load < seed[sn].span)
				continue;
			/*
			 * For simplicity we throw out the first kmer when in colour space. If
			 * we did not do so, we'd run into a ton of colour-letter space
			 * headaches. For instance, how should the first read kmer match
			 * against a kmer from the genome? The first colour of the genome kmer
			 * depends on the previous letter in the genome, so we may have a
			 * matching read, but the colour representation doesn't agree due to
			 * different initialising bases.
			 *
			 * If we wanted to be complete, we could compute the four permutations
			 * and add them, but I'm not so sure that'd be a good idea. Perhaps
			 * this should be investigated in the future.
			 */
			if (shrimp_mode == MODE_COLOUR_SPACE && i == seed[sn].span - 1)
				continue;

			uint32_t mapidx = kmer_to_mapidx_orig(kmerWindow, sn);
			*len += genomemap_len[sn][mapidx];
			if(*len){
				*list = (uint32_t *)xrealloc(*list,sizeof(uint32_t)* *len);
				memcpy(*list + *len - genomemap_len[sn][mapidx],genomemap[sn][mapidx],genomemap_len[sn][mapidx]*sizeof(uint32_t));
			}

			mapidx = kmer_to_mapidx_orig(kmerWindow_rc, sn);
			*len_rc += genomemap_len[sn][mapidx];
			if(*len_rc){
				*list_rc = (uint32_t *)xrealloc(*list_rc,sizeof(uint32_t)* *len_rc);
				memcpy(*list_rc + *len_rc - genomemap_len[sn][mapidx],genomemap[sn][mapidx],genomemap_len[sn][mapidx]*sizeof(uint32_t));
			}

		}
	}
	qsort(*list,*len,sizeof(uint32_t),&comp);
	qsort(*list_rc,*len_rc,sizeof(uint32_t),&comp);
}

/*
 * Find matches for the given read.
 *
 * list[0..len-1] is a sorted list containing start locations
 * of matching spaced kmers between this read and the genome.
 *
 * Does not assume list of locations is tagged with actual kmer and originating seed.
 * Instead, use avg_seed_span and center genome window around matching locations.
 *
 * Assumes that all contigs (num_contigs) that were indexed
 * are loaded as letter space bitmaps in genome[].
 * If the reads are colourspace, we also expect contigs in colour space in genome_cs[].
 */
static void
scan_read_lscs_pass1(read_entry *re, uint32_t *list, uint len, bool rev_cmpl){
  uint cn;
  uint32_t goff, glen;
  uint i, j;
  int score, thresh;
  score = 0; //Shut up compliler TODO don't do this

  cn = 0;
  for (i = num_matches - 1, j = 0; i < len; i++, j++) {
    assert(i == j + num_matches - 1);

    /* adjust current contig number */
    while (cn < num_contigs - 1 && list[i] >= contig_offsets[cn + 1])
      cn++;
    assert (contig_offsets[cn] <= list[i] && list[i] < contig_offsets[cn] + genome_len[cn]);

    /* test if last num_matches are in same contig and fit in a window of size window_len */
    if (contig_offsets[cn] <= list[j]
	&& list[j] + re->window_len >= list[i] + avg_seed_span) {

      /* set goff, glen */
      if ((re->window_len - (list[i] + avg_seed_span - list[j]))/2 <= list[j])
	goff = list[j] - (re->window_len - (list[i] + avg_seed_span - list[j]))/2;
      else
	goff = 0;
      glen = re->window_len;

      if (goff + glen > genome_len[cn])
	continue;

      // optionally use cache for this genomic window

      if (shrimp_mode == MODE_COLOUR_SPACE) {
	score = sw_vector(genome_cs_contigs[cn], goff, glen,
			  rev_cmpl? re->read_rc : re->read, re->read_len,
			  genome_contigs[cn], rev_cmpl? re->initbp_rc : re->initbp, genome_is_rna);
      } else if (shrimp_mode == MODE_LETTER_SPACE) {
	score = sw_vector(genome_contigs[cn], goff, glen,
			  rev_cmpl? re->read_rc : re->read, re->read_len,
			  NULL, -1, genome_is_rna);
      }

      thresh = (int)abs_or_pct(sw_vect_threshold, match_score * re->read_len);
      if (score >= thresh) {
	/* save hit */
	//ES_count_increment(&es_total_filter_passes);
	//ES_count32_increment(&re->es_filter_passes);

	save_score(re, score, goff, cn, rev_cmpl);
	re->sw_hits++;

	/* advance i&j appropriately */
	j++;
	while (j < len
	       && !(cn < num_contigs - 1 && list[j] >= contig_offsets[cn + 1])
	       && !(list[j] >= contig_offsets[cn] + goff + re->window_len - (int)abs_or_pct(window_overlap, re->read_len))
	       )
	  j++;
	j--;
	i = j + num_matches - 1;
      }
    }
  }
}

/*
 * Do a final pass for given read.
 * Highest scoring matches are in scores heap.
 */
static void
scan_read_lscs_pass2(read_entry * re) {
  uint i;
  int thresh = (int)abs_or_pct(sw_full_threshold, match_score * re->read_len);
  struct sw_full_results last_sfr;

  assert(re != NULL && re->scores != NULL);

  /* compute full alignment scores */
  for (i = 1; i <= re->scores[0].heap_elems; i++) {
    struct re_score * rs = &re->scores[i];
    uint32_t * gen;
    uint goff;

    rs->sfrp = (struct sw_full_results *)xmalloc(sizeof(struct sw_full_results));

    if (rs->rev_cmpl) {
      gen = genome_contigs_rc[rs->contig_num];
      goff = rs->g_idx + re->window_len - 1;
    } else {
      gen = genome_contigs[rs->contig_num];
      goff = rs->g_idx;
    }

    if (shrimp_mode == MODE_COLOUR_SPACE) {
	sw_full_cs(gen, goff, re->window_len,
		   re->read, re->read_len, re->initbp,
		   thresh, rs->sfrp, rs->rev_cmpl && Tflag, genome_is_rna, NULL, 0);
    } else {
      /*
       * The full SW in letter space assumes it's given the correct max score.
       * This might not be true just yet if we're using hashing&caching because
       * of possible hash collosions.
       */
      rs->score = sw_vector(gen, goff, re->window_len,
			    re->read, re->read_len,
			    NULL, -1, genome_is_rna);
      if (rs->score >= thresh) {
	sw_full_ls(gen, goff, re->window_len,
		   re->read, re->read_len,
		   thresh, rs->score, rs->sfrp, rs->rev_cmpl && Tflag,  NULL, 0);
	assert(rs->sfrp->score == rs->score);
      } else { // this wouldn't have passed the filter; eliminated in loop below
	rs->sfrp->score = rs->score;
	rs->sfrp->dbalign = NULL;
	rs->sfrp->qralign = NULL;
      }
    }
  }

  /* sort scores */
  qsort(&re->scores[1], re->scores[0].heap_elems, sizeof(re->scores[1]), score_cmp);

  /* Output sorted list, removing any duplicates. */
  memset(&last_sfr, 0, sizeof(last_sfr));
  for (i = 1; i <= re->scores[0].heap_elems; i++) {
    struct re_score * rs = &re->scores[i];
    bool dup;

    dup = sw_full_results_equal(&last_sfr, rs->sfrp);
    if (dup)
      nduphits++;

    if (rs->score >= thresh && dup == false) {
    	char *output;

    	re->final_matches++;

    	output = output_normal(re->name, contig_names[rs->contig_num], rs->sfrp,
    			genome_len[rs->contig_num], (shrimp_mode == MODE_COLOUR_SPACE), re->read,
    			re->read_len, re->initbp, rs->rev_cmpl, Rflag);
    	puts(output);
    	free(output);

    	if (Pflag) {
    		output = output_pretty(re->name, contig_names[rs->contig_num], rs->sfrp,
    				genome_contigs[rs->contig_num], genome_len[rs->contig_num],
    				(shrimp_mode == MODE_COLOUR_SPACE), re->read,
    				re->read_len, re->initbp, rs->rev_cmpl);
    		puts("");
    		puts(output);
    		free(output);
    	}
    }

    last_sfr = *rs->sfrp;
    if (rs->sfrp->dbalign != NULL)
      free(rs->sfrp->dbalign);
    if (rs->sfrp->qralign != NULL)
      free(rs->sfrp->qralign);
    free(rs->sfrp);
  }

}


static void
handle_read(read_entry *re){
	int len,len_rc = 0;
	uint32_t *list = NULL, *list_rc = NULL;
	DEBUG("computing read locations for read %s",re->name);
	read_locations(re,&len,&len_rc,&list,&list_rc);
#ifdef DEBUGGING
	int i;
	fprintf(stderr,"read locations:");
	for(i = 0; i < len; i++){
		fprintf(stderr,"%u,",list[i]);
	}
	fprintf(stderr,"\n");
#endif

  // initialize heap of best hits for this read
  re->scores = (struct re_score *)xcalloc((num_outputs + 1) * sizeof(struct re_score));
  re->scores[0].heap_elems = 0;
  re->scores[0].heap_capacity = num_outputs;

	DEBUG("fisrt pass on list");
	scan_read_lscs_pass1(re,list,len,false);
	DEBUG("fisrt pass on list_rc");
	scan_read_lscs_pass1(re,list_rc,len_rc,true);
	DEBUG("second pass");
	scan_read_lscs_pass2(re);

  // TODO This v
  // done with this read; deallocate memory.
  // deallocate name/bitmap ???
  free(re->scores);
}

/*
 * Launch the threads that will scan the reads
 */

static bool
launch_scan_threads(const char *file){
	fasta_t fasta;
	int space;
	int s = chunkSize*numThreads;
	read_entry *res;
	char *seq;

	bool is_rna;

	res = (read_entry *)xmalloc(sizeof(read_entry)*s);

	//open the fasta file and check for errors
	if (shrimp_mode == MODE_LETTER_SPACE)
			space = LETTER_SPACE;
		else
			space = COLOUR_SPACE;
	fasta = fasta_open(file,space);
	if (fasta == NULL) {
		fprintf(stderr,"error: failded to open read file [%s]\n",file);
		return (false);
	} else {
		fprintf(stderr,"- Processing read file [%s]\n",file);
	}

	//read the fasta file, s sequences at a time, and process in threads.
	bool more = true;
	while (more){
		int i;
		for(i = 0; i < s; i++){
			if(!fasta_get_next(fasta, &res[i].name, &seq, &is_rna)){
				more = false;
				break;
			}
			res[i].read = fasta_sequence_to_bitfield(fasta, seq);
			res[i].read_len = strlen(seq);
			if (shrimp_mode == MODE_COLOUR_SPACE){
				res[i].initbp = fasta_get_initial_base(fasta,seq);
				res[i].initbp_rc = res[i].initbp;
				res[i].read_rc = reverse_complement_read_cs(res[i].read,&res[i].initbp_rc,res[i].read_len,is_rna);
			} else {
				res[i].read_rc = reverse_complement_read_ls(res[i].read,res[i].read_len,is_rna);
			}

			res[i].window_len = (uint16_t)abs_or_pct(window_len,res[i].read_len);
		}
#pragma omp parallel shared(res,i) num_threads(numThreads)
		{
			int j;
#pragma omp for
			for (j = 0; j < i; j++){
				handle_read(&res[j]);
			}
		}

	}
	return true;
}

/*
 * index the kmers in the genome contained in the file.
 * This can then be used to align reads against.
 */
static bool
load_genome(char **files, int nfiles)
{
	fasta_t fasta;
	size_t seqlen,  bytes;
	uint32_t *read;
	char *seq, *name;
	uint32_t *kmerWindow;
	uint sn;
	char *file;
	bool is_rna;

	//allocate memory for the genome map
	genomemap = (uint32_t ***) xmalloc_c(n_seeds * sizeof(genomemap[0]),
			&mem_genomemap);
	genomemap_len = (uint32_t **)xmalloc_c(n_seeds * sizeof(genomemap_len[0]),
			&mem_genomemap);


	for (sn = 0; sn < n_seeds; sn++) {

		//TODO impliment Hflag
		//bytes = sizeof(uint32_t *) * power(4, HASH_TABLE_POWER);
		bytes = sizeof(uint32_t *) * power(4, seed[sn].weight);

		genomemap[sn] = (uint32_t **)xcalloc_c(bytes, &mem_genomemap);
		genomemap_len[sn] = (uint32_t *)xcalloc_c(bytes,&mem_genomemap);

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
				exit(1);
			}

			seqlen = strlen(seq);
			if (seqlen == 0) {
				fprintf(stderr, "error: genome [%s] had no sequence!\n",
						name);
				exit(1);
			}

			read = fasta_sequence_to_bitfield(fasta,seq);

			if (read == NULL) {
				fprintf(stderr, "error: invalid sequence; tag: [%s]\n", name);
				exit(1);
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
				genome_initbp[num_contigs -1] = EXTRACT(read,0);
				read = fasta_bitfield_to_colourspace(fasta,read,seqlen,is_rna);
				genome_initbp = (uint32_t *)xrealloc(genome_initbp,sizeof(uint32_t *)*num_contigs);
			}
			kmerWindow = (uint32_t *)xcalloc(sizeof(kmerWindow[0])*BPTO32BW(max_seed_span));
			u_int load;
			DEBUG("looping seq");
			for (load = 0; i < seqlen + contig_offsets[num_contigs-1]; i++) {
				uint base, sn;

				base = EXTRACT(read, i - contig_offsets[num_contigs-1]);
				bitfield_prepend(kmerWindow, max_seed_span, base);

				//skip past any Ns or Xs
				if (base == BASE_N || base == BASE_X)
					load = 0;
				else if (load < max_seed_span)
					load++;
				DEBUG("looping seeds");
				for (sn = 0; sn < n_seeds; sn++) {
					if (load < seed[sn].span)
						continue;


					DEBUG("hashing");
					uint32_t mapidx = kmer_to_mapidx_orig(kmerWindow, sn);
					//increase the match count and store the location of the match
					DEBUG("updating len");
					genomemap_len[sn][mapidx]++;
					DEBUG("reallocing map");
					genomemap[sn][mapidx] = (uint32_t *)xrealloc_c(genomemap[sn][mapidx],
							sizeof(uint32_t) * (genomemap_len[sn][mapidx]),
							sizeof(uint32_t) * (genomemap_len[sn][mapidx] - 1),
							&mem_genomemap);
					if (genomemap[sn][mapidx] == NULL){
						DEBUG("realloc returned null");
					}
					DEBUG("updateing map");
					genomemap[sn][mapidx][genomemap_len[sn][mapidx] - 1] = i - seed[sn].span + 1;
					nkmers++;

				}
			}
			DEBUG("done indexing sequence");
			free(seq);
			seq = name = NULL;

			free(kmerWindow);
		}
		fasta_close(fasta);
	}
	fprintf(stderr,"Loaded Genome\n");

	return (true);

}

void usage(char *progname,bool full_usage){
	char *slash;
	uint sn;

	load_default_seeds();
	slash = strrchr(progname, '/');
	if (slash != NULL)
		progname = slash + 1;
	fprintf(stderr, "usage: %s [parameters] [options] "
			"genome_file1\n", progname);

	fprintf(stderr, "Parameters:\n");

	fprintf(stderr,
			"    -s    Spaced Seed(s)                          (default: ");
	for (sn = 0; sn < n_seeds; sn++) {
		if (sn > 0)
			fprintf(stderr, "                                                            ");
		fprintf(stderr, "%s%s\n", seed_to_string(sn), (sn == n_seeds - 1? ")" : ","));
	}
	fprintf(stderr,
			"    -?    Full list of parameters and options\n");

}

void print_settings() {
  static char const my_tab[] = "    ";
  uint i;

  fprintf(stderr, "--------------------------------------------------------------------------------\n");
  fprintf(stderr, "gmapper: %s.\nSHRiMP %s\n[%s]\n", get_mode_string(), SHRIMP_VERSION_STRING, get_compiler());
  fprintf(stderr, "--------------------------------------------------------------------------------\n");

  fprintf(stderr, "Settings:\n");
  fprintf(stderr, "%s%-40s%s (%u/%u)\n", my_tab,
	  (n_seeds == 1) ? "Spaced Seed (weight/span)" : "Spaced Seeds (weight/span)",
	  seed_to_string(0), seed[0].weight, seed[0].span);
  for (i = 1; i < n_seeds; i++) {
    fprintf(stderr, "                                          %s (%u/%u)\n",
	    seed_to_string(i), seed[i].weight, seed[i].span);
  }

  fprintf(stderr, "%s%-40s%u\n", my_tab, "Seed Matches per Window:", num_matches);

  if (IS_ABSOLUTE(window_len)) {
    fprintf(stderr, "%s%-40s%u\n", my_tab, "Seed Window Length:", (uint)-window_len);
  } else {
    fprintf(stderr, "%s%-40s%.02f%%\n", my_tab, "Seed Window Length:", window_len);
  }

  if (IS_ABSOLUTE(window_overlap)) {
    fprintf(stderr, "%s%-40s%u\n", my_tab, "Seed Window Overlap Length:", (uint)-window_overlap);
  } else {
    fprintf(stderr, "%s%-40s%.02f%%\n", my_tab, "Seed Window Overlap Length:", window_overlap);
  }

  fprintf(stderr, "\n");
  fprintf(stderr, "%s%-40s%d\n", my_tab, "S-W Match Score:", match_score);
  fprintf(stderr, "%s%-40s%d\n", my_tab, "S-W Mismatch Score:", mismatch_score);
  fprintf(stderr, "%s%-40s%d\n", my_tab, "S-W Gap Open Score (Ref):", a_gap_open_score);
  fprintf(stderr, "%s%-40s%d\n", my_tab, "S-W Gap Open Score (Qry):", b_gap_open_score);
  fprintf(stderr, "%s%-40s%d\n", my_tab, "S-W Gap Extend Score (Ref):", a_gap_extend_score);
  fprintf(stderr, "%s%-40s%d\n", my_tab, "S-W Gap Extend Score (Qry):", b_gap_extend_score);
  if (shrimp_mode == MODE_COLOUR_SPACE) {
    fprintf(stderr, "%s%-40s%d\n", my_tab, "S-W Crossover Score:", crossover_score);
  }
  if (shrimp_mode == MODE_COLOUR_SPACE) {
    if (IS_ABSOLUTE(sw_vect_threshold)) {
      fprintf(stderr, "%s%-40s%u\n", my_tab, "S-W Vector Hit Threshold:", (uint)-sw_vect_threshold);
    } else {
      fprintf(stderr, "%s%-40s%.02f%%\n", my_tab, "S-W Vector Hit Threshold:", sw_vect_threshold);
    }
  }
  if (IS_ABSOLUTE(sw_full_threshold)) {
    fprintf(stderr, "%s%-40s%u\n", my_tab,
	    shrimp_mode == MODE_COLOUR_SPACE? "S-W Full Hit Threshold:" : "S-W Hit Threshold",
	    (uint)-sw_full_threshold);
  } else {
    fprintf(stderr, "%s%-40s%.02f%%\n", my_tab,
	    shrimp_mode == MODE_COLOUR_SPACE? "S-W Full Hit Threshold:" : "S-W Hit Threshold",
	    sw_full_threshold);
  }

}


#ifdef DEBUGMAIN
int main(int argc, char **argv){
	if (argc > 1){
		set_mode_from_argv(argv);
		if (n_seeds == 0)
				load_default_seeds();
		init_seed_hash_mask();

  print_settings();

  int max_window_len = 200;
  int longest_read_len = 100;
  if (sw_vector_setup(max_window_len, longest_read_len,
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
			   crossover_score, false, 0);
  } else {
    ret = sw_full_ls_setup(max_window_len, longest_read_len,
			   a_gap_open_score, a_gap_extend_score, b_gap_open_score, b_gap_extend_score,
			   match_score, mismatch_score, false, 0);
  }
  if (ret) {
    fprintf(stderr, "failed to initialise scalar "
	    "Smith-Waterman (%s)\n", strerror(errno));
    exit(1);
  }


		if (1){
			fprintf(stderr,"loading gneomfile\n");
			load_genome(argv+1,2);
		}
		if (0){
			fprintf(stderr,"saving compressed index\n");
			save_genome_map("testfile.gz");
		}
		if (0){
			print_info();
		}
		if(1){
			char * output;

			output = output_format_line(Rflag);
			puts(output);
			free(output);
		  launch_scan_threads(argv[3]);
		}
		if (0){
			fprintf(stderr,"loading compressed index\n");
			load_genome_map("testfile.gz");
		}
		fprintf(stderr,"done\n");
	}
}

#else
int main(int argc, char **argv){
	char *genome_file;
	char *progname = argv[0];
	char const * optstr;
	char *c;
	int ch;

	set_mode_from_argv(argv);

	switch(shrimp_mode){
	case MODE_COLOUR_SPACE:
		optstr = "?s:";
		break;
	case MODE_LETTER_SPACE:
		optstr = "?s:";
		break;
	case MODE_HELICOS_SPACE:
		fprintf(stderr,"Helicose currently unsuported\n");
		exit(1);
		break;
	default:
		assert(0);
	}

	fprintf(stderr, "--------------------------------------------------"
			"------------------------------\n");
	fprintf(stderr, "gmapper: %s.\nSHRiMP %s\n[%s]\n", get_mode_string(),
			SHRIMP_VERSION_STRING, get_compiler());
	fprintf(stderr, "--------------------------------------------------"
			"------------------------------\n");

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
		case '?':
			usage(progname, true);
			break;
		default:
			usage(progname, false);
		}
	}

	argc -= optind;
	argv += optind;

	if (n_seeds == 0)
			load_default_seeds();

	init_seed_hash_mask();

	if (argc < 1){
		fprintf(stderr, "Genome File not specified\n");
		exit(1);
	}

	genome_file = argv[0];

	fprintf(stderr,"loading gneomfile\n");
	load_genome(&genome_file,1);
	fprintf(stderr, "        Genomemap:                %s\n",
					comma_integer(count_get_count(&mem_genomemap)));
	fprintf(stderr,"saving compressed index\n");
	save_genome_map("testfile.gz");
	print_info();
	fprintf(stderr,"loading compressed index\n");
	load_genome_map("testfile.gz");

	launch_scan_threads(genome_file);

}
#endif
