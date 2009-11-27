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
static int	anchor_width		= DEF_ANCHOR_WIDTH;

/* Flags */
static int Bflag = false;			/* print a progress bar */
static int Cflag = false;			/* do complement only */
static int Fflag = false;			/* do positive (forward) only */
static int Hflag = false;			/* use hash table, not lookup */
static int Pflag = false;			/* pretty print results */
static int Rflag = false;			/* add read sequence to output*/
static int Tflag = false;			/* reverse sw full tie breaks */
static int Uflag = false;			/* output unmapped reads, too */
static int Mflag = true;			/* print memory usage stats */
static int Dflag = false;			/* print statistics for each thread */


/* Statistics */
static count_t	dup_hits_c;			/* number of duplicate hits */
static count_t	reads_c;
static stat_t	matches_s;

static uint64_t reads_matched = 0;
static uint64_t total_matches = 0;

static uint64_t nreads = 0;

static uint64_t scan_usecs;
static uint64_t map_usecs;

static count_t mem_genomemap;

size_t capacity;

/* files to use when saving and loading genome maps */
char *save_file = NULL;
char *load_file = NULL;

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


static uint num_threads = DEF_NUM_THREADS;
static uint chunk_size = DEF_CHUNK_SIZE;

/* kmer_to_mapidx function */
static uint32_t (*kmer_to_mapidx)(uint32_t *, u_int) = NULL;

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

static bool
valid_spaced_seeds()
{
	uint sn;

	for (sn = 0; sn < n_seeds; sn++) {
		if (seed[sn].weight > MAX_SEED_WEIGHT && !Hflag)
			return false;

		if (Hflag && (seed[sn].span > MAX_HASH_SEED_SPAN ||
				seed[sn].weight > MAX_HASH_SEED_WEIGHT))
			return false;
	}

	return true;
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
	gzwrite(fp,&shrimp_mode,sizeof(shrimp_mode_t));
	uint32_t h = (uint32_t)Hflag;
	gzwrite(fp,&h,sizeof(int32_t));
	uint32_t i;
	gzwrite(fp,&num_contigs,sizeof(uint32_t));
	//fprintf(stderr,"saved num_contigs\n");
	for (i = 0; i < num_contigs; i++) {
		gzwrite(fp,&contig_offsets[i],sizeof(uint32_t));
		uint32_t len = strlen(contig_names[i]);
		gzwrite(fp,&len,sizeof(uint32_t));
		gzwrite(fp,contig_names[i],len +1);

		gzwrite(fp,genome_len + i,sizeof(uint32_t));
		gzwrite(fp,genome_contigs[i],sizeof(uint32_t)*BPTO32BW(genome_len[i]));
		gzwrite(fp,genome_contigs_rc[i],sizeof(uint32_t)*BPTO32BW(genome_len[i]));
		if(shrimp_mode == MODE_COLOUR_SPACE){
			gzwrite(fp,genome_cs_contigs[i],sizeof(uint32_t)*BPTO32BW(genome_len[i]));
			gzwrite(fp,genome_initbp+i,sizeof(uint32_t));
		}
	}
	//fprintf(stderr,"saving seeds\n");
	//write the seeds and genome_maps
	gzwrite(fp,&n_seeds,sizeof(uint32_t));
	for (i = 0; i < n_seeds; i++) {
		gzwrite(fp,&seed[i], sizeof(seed_type));

		//write the genome_map for this seed
		uint32_t j;
		//uint32_t p = power(4, seed[i].weight);
		//fprintf(stderr,"saving index\n");
		gzwrite(fp,&capacity, sizeof(size_t));
		for (j = 0; j < capacity; j++) {
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
	gzread(fp,&shrimp_mode,sizeof(shrimp_mode_t)); //TODO make sure no conflict
	int32_t f;
	gzread(fp,&f,sizeof(int32_t));
	Hflag = f;
	gzread(fp,&num_contigs,sizeof(uint32_t));
	//fprintf(stderr,"num_contigs = %u\n",num_contigs);

	contig_names = (char **)xrealloc(contig_names,sizeof(char *)*num_contigs);
	contig_offsets = (uint32_t *)xrealloc(contig_offsets,sizeof(uint32_t)*num_contigs);

	genome_len = (uint32_t *)xrealloc(genome_len,sizeof(uint32_t)*num_contigs);
	genome_contigs = (uint32_t **)xrealloc(genome_contigs,sizeof(uint32_t *)*num_contigs);
	genome_contigs_rc = (uint32_t **)xrealloc(genome_contigs,sizeof(uint32_t *)*num_contigs);
	if(shrimp_mode == MODE_COLOUR_SPACE){
		genome_cs_contigs = (uint32_t **)xrealloc(genome_cs_contigs,sizeof(uint32_t *)*num_contigs);
		genome_initbp = (uint32_t *)xrealloc(genome_initbp,sizeof(uint32_t)*num_contigs);
	}

	for (i = 0; i < num_contigs; i++){
		gzread(fp,&contig_offsets[i],sizeof(uint32_t));

		uint32_t len;
		gzread(fp,&len,sizeof(uint32_t));
		contig_names[i] = (char *)xrealloc(contig_names[i],sizeof(char)*len);
		gzread(fp,contig_names[i],len+1);

		gzread(fp,genome_len + i,sizeof(uint32_t));

		genome_contigs[i] = (uint32_t *)xrealloc(genome_contigs[i],sizeof(uint32_t)*BPTO32BW(genome_len[i]));
		gzread(fp,genome_contigs[i],sizeof(uint32_t)*BPTO32BW(genome_len[i]));

		genome_contigs_rc[i] = (uint32_t *)xrealloc(genome_contigs_rc[i],sizeof(uint32_t)*BPTO32BW(genome_len[i]));
		gzread(fp,genome_contigs_rc[i],sizeof(uint32_t)*BPTO32BW(genome_len[i]));

		if(shrimp_mode == MODE_COLOUR_SPACE){
			genome_cs_contigs[i] = (uint32_t *)xrealloc(genome_cs_contigs[i],sizeof(uint32_t)*BPTO32BW(genome_len[i]));
			gzread(fp,genome_cs_contigs[i],sizeof(uint32_t)*BPTO32BW(genome_len[i]));
			gzread(fp,genome_initbp+i,sizeof(uint32_t));
		}
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
		gzread(fp,&capacity,sizeof(size_t));
		genomemap_len[i] = (uint32_t *)xrealloc(genomemap_len[i],sizeof(uint32_t)*capacity);
		genomemap[i] = (uint32_t **)xrealloc(genomemap[i],sizeof(uint32_t *)*capacity);
		for (j = 0; j < capacity; j++){
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

			uint32_t mapidx = kmer_to_mapidx(kmerWindow, sn);
			*len += genomemap_len[sn][mapidx];
			if(*len){
				*list = (uint32_t *)xrealloc(*list,sizeof(uint32_t)* *len);
				memcpy(*list + *len - genomemap_len[sn][mapidx],genomemap[sn][mapidx],genomemap_len[sn][mapidx]*sizeof(uint32_t));
			}

			mapidx = kmer_to_mapidx(kmerWindow_rc, sn);
			*len_rc += genomemap_len[sn][mapidx];
			if(*len_rc){
				*list_rc = (uint32_t *)xrealloc(*list_rc,sizeof(uint32_t)* *len_rc);
				memcpy(*list_rc + *len_rc - genomemap_len[sn][mapidx],genomemap[sn][mapidx],genomemap_len[sn][mapidx]*sizeof(uint32_t));
			}

		}
	}
	qsort(*list,*len,sizeof(uint32_t),&comp);
	qsort(*list_rc,*len_rc,sizeof(uint32_t),&comp);

	free(kmerWindow);
	free(kmerWindow_rc);
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
	uint first, last;
	int score = 0; //Shut up compliler TODO don't do this

	/*
  uint i, j;
  cn = 0;
  for (i = num_matches - 1, j = 0; i < len; i++, j++) {
    assert(i == j + num_matches - 1);

    // adjust current contig number
    while (cn < num_contigs - 1 && list[i] >= contig_offsets[cn + 1])
      cn++;
    assert (contig_offsets[cn] <= list[i] && list[i] < contig_offsets[cn] + genome_len[cn]);

    // test if last num_matches are in same contig and fit in a window of size window_len
    if (contig_offsets[cn] <= list[j]
	&& list[j] + re->window_len >= list[i] + avg_seed_span) {

      // set goff, glen
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

      if (score >= (int)abs_or_pct(sw_vect_threshold, match_score * re->read_len)) {
	// save hit
	//ES_count_increment(&es_total_filter_passes);
	//ES_count32_increment(&re->es_filter_passes);

	save_score(re, score, goff, cn, rev_cmpl);
	re->sw_hits++;

	// advance i&j appropriately
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
	 */

	cn = 0;
	first = 0;
	last = 0;
	do {
		assert(0 <= first && first < len);
		assert(first <= last && last < len);

		// fix contig num
		while (cn < num_contigs - 1 && list[first] >= contig_offsets[cn + 1])
			cn++;
		assert(contig_offsets[cn] <= list[first] && list[last] < contig_offsets[cn] + genome_len[cn]);

		// advance last as far as possible
		while (last < len - 1
				&& list[last + 1] < contig_offsets[cn] + genome_len[cn]
				                                                    && list[last + 1] + avg_seed_span <= list[first] + re->window_len)
			last++;

		if (last - first + 1 >= num_matches) {
			// fit window around first...last
			if ((re->window_len - (list[last] + avg_seed_span - list[first]))/2 <= list[first] - contig_offsets[cn])
				goff = list[first] - contig_offsets[cn] - (re->window_len - (list[last] + avg_seed_span - list[first]))/2;
			else
				goff = 0;
			glen = re->window_len;

			if (goff + glen > genome_len[cn])
				glen = genome_len[cn] - goff;

			//DEBUG( "trying read:[%s] cn:%u goff:%u glen:%u", re->name, cn, goff, glen);

			// run sw filter on window goff, glen
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

			if (score >= (int)abs_or_pct(sw_vect_threshold, match_score * re->read_len)) {
				// save hit and continue
				//ES_count_increment(&es_total_filter_passes);
				//ES_count32_increment(&re->es_filter_passes);

				save_score(re, score, goff, cn, rev_cmpl);
				re->sw_hits++;

				while (first <= last
						&& list[first] - contig_offsets[cn] < goff + glen - (int)abs_or_pct(window_overlap, re->read_len))
					first++;
			} else {
				// filter call didn't pass threshold
				if (last == len - 1
						|| (cn < num_contigs - 1 && list[last + 1] >= contig_offsets[cn + 1])) {
					first = last + 1;
				} else {
					last++;
					first++;
					assert(re->window_len > avg_seed_span);
					while (list[first] + re->window_len < list[last] + avg_seed_span)
						first++;
				}
			}
		} else {
			// not enough matches to run filter
			first++;
		}

		if (last < first)
			last = first;
	} while (first < len);

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
		uint goff, glen;

		rs->sfrp = (struct sw_full_results *)xmalloc(sizeof(struct sw_full_results));

		if (rs->rev_cmpl) {
			gen = genome_contigs_rc[rs->contig_num];
			goff = rs->g_idx + re->window_len - 1;
			if (goff > genome_len[rs->contig_num] - 1) {
			  goff = 0;
			} else {
			  goff = genome_len[rs->contig_num] - 1 - goff;
			}
		} else {
			gen = genome_contigs[rs->contig_num];
			goff = rs->g_idx;
		}

		glen = re->window_len;
		if (goff + glen > genome_len[rs->contig_num])
		  glen = genome_len[rs->contig_num] - goff;

		if (shrimp_mode == MODE_COLOUR_SPACE) {
			sw_full_cs(gen, goff, glen,
					re->read, re->read_len, re->initbp,
					thresh, rs->sfrp, rs->rev_cmpl && Tflag, genome_is_rna, NULL, 0);
		} else {
			/*
			 * The full SW in letter space assumes it's given the correct max score.
			 * This might not be true just yet if we're using hashing&caching because
			 * of possible hash collosions.
			 */
			rs->score = sw_vector(gen, goff, glen,
					re->read, re->read_len,
					NULL, -1, genome_is_rna);
			if (rs->score >= thresh) {
				sw_full_ls(gen, goff, glen,
						re->read, re->read_len,
						thresh, rs->score, rs->sfrp, rs->rev_cmpl && Tflag,  NULL, 0);
				assert(rs->sfrp->score == rs->score);
			} else { // this wouldn't have passed the filter; eliminated in loop below
				rs->sfrp->score = rs->score;
				rs->sfrp->dbalign = NULL;
				rs->sfrp->qralign = NULL;
			}
		}
		rs->score = rs->sfrp->score;
	}

	/* sort scores */
	qsort(&re->scores[1], re->scores[0].heap_elems, sizeof(re->scores[1]), score_cmp);

	/* Output sorted list, removing any duplicates. */
	if (re->scores[0].heap_elems > 0){
		reads_matched++;
	}
	memset(&last_sfr, 0, sizeof(last_sfr));
	for (i = 1; i <= re->scores[0].heap_elems; i++) {
		struct re_score * rs = &re->scores[i];
		bool dup;

		dup = sw_full_results_equal(&last_sfr, rs->sfrp);
		if (dup)
			count_increment(&dup_hits_c);

		if (rs->score >= thresh && dup == false) {
			char *output1 = NULL,*output2 = NULL,*format;

			re->final_matches++;

			output1 = output_normal(re->name, contig_names[rs->contig_num], rs->sfrp,
					genome_len[rs->contig_num], (shrimp_mode == MODE_COLOUR_SPACE), re->read,
					re->read_len, re->initbp, rs->rev_cmpl, Rflag);

			if (Pflag) {
				output2 = output_pretty(re->name, contig_names[rs->contig_num], rs->sfrp,
						genome_contigs[rs->contig_num], genome_len[rs->contig_num],
						(shrimp_mode == MODE_COLOUR_SPACE), re->read,
						re->read_len, re->initbp, rs->rev_cmpl);
			}
			if(output2 != NULL){
				format = (char *)"%s\n\n%s\n";
			} else{
				format = (char *)"%s\n";
			}
#pragma omp critical (stdout)
			{
				fprintf(stdout,format,output1,output2);
			}
			free(output1);
			free(output2);

			total_matches++;
		}

		last_sfr = *rs->sfrp;
		if (rs->sfrp->dbalign != NULL)
			free(rs->sfrp->dbalign);
		if (rs->sfrp->qralign != NULL)
			free(rs->sfrp->qralign);
		free(rs->sfrp);
	}

	stat_add(&matches_s, re->final_matches);
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
	fprintf(stderr,"read locations revcmpl:");
	for(i = 0; i < len_rc; i++){
		fprintf(stderr,"%u,",list_rc[i]);
	}
	fprintf(stderr,"\n");
#endif

	// initialize heap of best hits for this read
	re->scores = (struct re_score *)xcalloc((num_outputs + 1) * sizeof(struct re_score));
	re->scores[0].heap_elems = 0;
	re->scores[0].heap_capacity = num_outputs;

	if (len > 0 && Fflag) {
		DEBUG("first pass on list");
		scan_read_lscs_pass1(re,list,len,false);
	}
	if (len_rc > 0 && Cflag) {
		DEBUG("first pass on list_rc");
		scan_read_lscs_pass1(re,list_rc,len_rc,true);
	}
	if (re->scores[0].heap_elems > 0) {
		DEBUG("second pass");
		scan_read_lscs_pass2(re);
	}

	// done with this read; deallocate memory.
	free(re->scores);
	free(re->name);
	free(re->read);
	free(re->read_rc);
	free(list);
	free(list_rc);

}

/*
 * Launch the threads that will scan the reads
 */

static bool
launch_scan_threads(const char *file){
	fasta_t fasta;
	int space;
	int s = chunk_size*num_threads;
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
			count_increment(&reads_c);
			if (shrimp_mode == MODE_COLOUR_SPACE){
			  res[i].read_len--;
				res[i].initbp = fasta_get_initial_base(fasta,seq);
				res[i].initbp_rc = res[i].initbp;
				res[i].read_rc = reverse_complement_read_cs(res[i].read, res[i].initbp, res[i].initbp_rc,
									    res[i].read_len, is_rna);
			} else {
				res[i].read_rc = reverse_complement_read_ls(res[i].read,res[i].read_len,is_rna);
			}

			res[i].window_len = (uint16_t)abs_or_pct(window_len,res[i].read_len);
			free(seq);
		}
		nreads += i;
		#pragma omp parallel shared(res,i) num_threads(num_threads)
		{
			int j;
			#pragma omp for
			for (j = 0; j < i; j++){
				handle_read(&res[j]);
			}
		}

	}
	free(res);
	return true;
}


void print_genomemap_stats() {
	stat_t list_size, list_size_non0;
	uint sn;
	uint64_t capacity, mapidx;

	fprintf(stderr, "Genome Map stats:\n");

	for (sn = 0; sn < n_seeds; sn++) {
		capacity = power(4, seed[sn].weight);

		stat_init(&list_size);
		stat_init(&list_size_non0);
		for (mapidx = 0; mapidx < capacity; mapidx++) {
			stat_add(&list_size, genomemap_len[sn][mapidx]);
			if (genomemap_len[sn][mapidx] > 0)
				stat_add(&list_size_non0, genomemap_len[sn][mapidx]);
		}

		fprintf(stderr, "sn:%u weight:%u total_kmers:%llu lists:%llu (non-zero:%llu) list_sz_avg:%.2f (%.2f) list_sz_stddev:%.2f (%.2f)\n",
				sn, seed[sn].weight, (long long unsigned int)stat_get_sum(&list_size),
				(long long unsigned int)capacity, (long long unsigned int)stat_get_count(&list_size_non0),
				stat_get_mean(&list_size), stat_get_mean(&list_size_non0),
				stat_get_sample_stddev(&list_size), stat_get_sample_stddev(&list_size_non0));
	}
}


void print_run_stats() {
	static char const my_tab[] = "    ";

	fprintf(stderr, "%sGeneral:\n", my_tab);
	fprintf(stderr, "%s%s%-24s" "%s    (%.4f%%)\n", my_tab, my_tab, "Reads Matched:",
			comma_integer(stat_get_count(&matches_s)),
			count_get_count(&reads_c) == 0? 0 : ((double)stat_get_count(&matches_s) / (double)count_get_count(&reads_c)) * 100);
	fprintf(stderr, "%s%s%-24s" "%s\n", my_tab, my_tab, "Total Matches:",
			comma_integer(stat_get_sum(&matches_s)));
	fprintf(stderr, "%s%s%-24s" "%.2f\n", my_tab, my_tab, "Avg Hits/Matched Read:",
			stat_get_sum(&matches_s) == 0 ? 0 : ((double)stat_get_sum(&matches_s) / (double)stat_get_count(&matches_s)));
	fprintf(stderr, "%s%s%-24s" "%s\n", my_tab, my_tab, "Duplicate Hits Pruned:",
			comma_integer(count_get_count(&dup_hits_c)));
}  


/*
 * index the kmers in the genome contained in the file.
 * This can then be used to align reads against.
 */
static bool
load_genome(char **files, int nfiles)
{
	fasta_t fasta;
	size_t seqlen;
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

		if(Hflag){
			capacity = sizeof(uint32_t *) * power(4, HASH_TABLE_POWER);
		} else {
			capacity = power(4, seed[sn].weight);
		}

		genomemap[sn] = (uint32_t **)xcalloc_c(sizeof(uint32_t *) * capacity, &mem_genomemap);
		genomemap_len[sn] = (uint32_t *)xcalloc_c(sizeof(uint32_t) * capacity,&mem_genomemap);

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
				genome_initbp = (uint32_t *)xrealloc(genome_initbp,sizeof(uint32_t *)*num_contigs);
				genome_initbp[num_contigs -1] = EXTRACT(read,0);
				read = fasta_bitfield_to_colourspace(fasta,read,seqlen,is_rna);
			}
			kmerWindow = (uint32_t *)xcalloc(sizeof(kmerWindow[0])*BPTO32BW(max_seed_span));
			u_int load;
			//DEBUG("indexing sequence");
			for (load = 0; i < seqlen + contig_offsets[num_contigs-1]; i++) {
				uint base;
				int sn;

				base = EXTRACT(read, i - contig_offsets[num_contigs-1]);
				bitfield_prepend(kmerWindow, max_seed_span, base);

				//skip past any Ns or Xs
				if (base == BASE_N || base == BASE_X)
					load = 0;
				else if (load < max_seed_span)
					load++;
				//DEBUG("looping seeds");
#pragma omp parallel for shared(kmerWindow,genomemap_len,genomemap,nkmers,mem_genomemap,load,seed) num_threads(num_threads)
				for (sn = 0; sn < (int)n_seeds; sn++) {
					if (load < seed[sn].span)
						continue;


					//DEBUG("hashing");
					uint32_t mapidx = kmer_to_mapidx(kmerWindow, sn);
					//increase the match count and store the location of the match
					//DEBUG("updating len");
					genomemap_len[sn][mapidx]++;
					//DEBUG("reallocing map");
					genomemap[sn][mapidx] = (uint32_t *)xrealloc_c(genomemap[sn][mapidx],
							sizeof(uint32_t) * (genomemap_len[sn][mapidx]),
							sizeof(uint32_t) * (genomemap_len[sn][mapidx] - 1),
							&mem_genomemap);
					//DEBUG("updateing map");
					genomemap[sn][mapidx][genomemap_len[sn][mapidx] - 1] = i - seed[sn].span + 1;
					nkmers++;

				}
			}
			//DEBUG("done indexing sequence");
			free(seq);
			seq = name = NULL;

			free(kmerWindow);
		}
		fasta_close(fasta);
	}
	fprintf(stderr,"Loaded Genome\n");

	return (true);

}

static void
print_statistics()
{
	double totalseedscantime = 0;
	double totalruntime = 0;
	uint64_t totalinvocs = 0;
	uint64_t totalcells = 0;
	double cellspersec[num_threads], hz, seedscantime[num_threads];
	uint64_t invocs[num_threads], cells[num_threads], fasta_load_ticks, vticks[num_threads];
	fasta_stats_t fs;

	fs = fasta_stats();
	fasta_load_ticks = fs->total_ticks;
	free(fs);
	hz = cpuhz();

#pragma omp parallel num_threads(num_threads) shared(hz)
	{
		int tid = omp_get_thread_num();
		sw_vector_stats(&invocs[tid], &cells[tid], &vticks[tid], &cellspersec[tid]);

		seedscantime[tid] = ((double)scan_usecs / 1000000.0) - (vticks[tid] / hz);
		seedscantime[tid] = MAX(0, seedscantime[tid]);
	}

	fprintf(stderr, "\nStatistics:\n");

	int i;
	for(i = 0; i < (int)num_threads; i++){
		if (Dflag){
			fprintf(stderr,"    Spaced Seed Scan[Thread %i]:\n",i);
			fprintf(stderr, "        Run-time:               %.2f seconds\n",
					seedscantime[i]);
			fprintf(stderr,"\n");
		}
		totalseedscantime += seedscantime[i];
	}



	fprintf(stderr, "    Spaced Seed Scan[Total]:\n");
	fprintf(stderr, "        Run-time:               %.2f seconds\n",
			totalseedscantime);
	fprintf(stderr, "        Total Kmers:            %s\n",
			comma_integer(nkmers));



	fprintf(stderr, "\n");

	for(i = 0; i < (int)num_threads; i++){
		if (Dflag){
			fprintf(stderr, "    Vector Smith-Waterman[Thread %i]:\n",i);
			fprintf(stderr, "        Run-time:               %.2f seconds\n",
					(cellspersec[i] == 0) ? 0 : cells[i] / cellspersec[i]);
			fprintf(stderr, "        Invocations:            %s\n",
					comma_integer(invocs[i]));
			fprintf(stderr, "        Cells Computed:         %.2f million\n",
					(double)cells[i] / 1.0e6);
			fprintf(stderr, "        Cells per Second:       %.2f million\n",
					cellspersec[i] / 1.0e6);
			fprintf(stderr, "\n");
		}
		totalruntime += (cellspersec[i] == 0) ? 0 : cells[i] / cellspersec[i];
		totalinvocs += invocs[i];
		totalcells += cells[i];

	}
	fprintf(stderr, "    Vector Smith-Waterman[Total]:\n");
	fprintf(stderr, "        Run-time:               %.2f seconds\n",
			totalruntime);
	fprintf(stderr, "        Invocations:            %s\n",
			comma_integer(totalinvocs));
	fprintf(stderr, "        Cells Computed:         %.2f million\n",
			(double)totalcells / 1.0e6);
	fprintf(stderr, "\n");

#pragma omp parallel num_threads(num_threads)
	{
		int tid = omp_get_thread_num();
		if (shrimp_mode == MODE_COLOUR_SPACE)
			sw_full_cs_stats(&invocs[tid], &cells[tid], NULL, &cellspersec[tid]);
		else
			sw_full_ls_stats(&invocs[tid], &cells[tid], NULL, &cellspersec[tid]);
	}


	totalruntime = 0;
	totalinvocs = 0;
	totalcells = 0;
	for(i = 0; i < (int)num_threads; i++){
		if (Dflag){
			fprintf(stderr, "    Scalar Smith-Waterman[Thread %i]:\n",i);
			fprintf(stderr, "        Run-time:               %.2f seconds\n",
					(cellspersec[i] == 0) ? 0 : cells[i] / cellspersec[i]);
			fprintf(stderr, "        Invocations:            %s\n",
					comma_integer(invocs[i]));
			fprintf(stderr, "        Cells Computed:         %.2f million\n",
					(double)cells[i] / 1.0e6);
			fprintf(stderr, "        Cells per Second:       %.2f million\n",
					cellspersec[i] / 1.0e6);

			fprintf(stderr, "\n");
		}
		totalruntime += (cellspersec[i] == 0) ? 0 : cells[i] / cellspersec[i];
		totalinvocs += invocs[i];
		totalcells += cells[i];

	}
	fprintf(stderr, "    Scalar Smith-Waterman[Total]:\n");
	fprintf(stderr, "        Run-time:               %.2f seconds\n",
			totalruntime);
	fprintf(stderr, "        Invocations:            %s\n",
			comma_integer(totalinvocs));
	fprintf(stderr, "        Cells Computed:         %.2f million\n",
			(double)totalcells / 1.0e6);

	fprintf(stderr, "\n");

	fprintf(stderr, "    Miscellaneous:\n");
	fprintf(stderr, "        Load Time:              %.2f seconds\n",
			(double)fasta_load_ticks / hz);
	fprintf(stderr, "        Map time:               %.2f seconds\n",
			(double)map_usecs / 1000000.0);
	fprintf(stderr, "        Scan time:              %.2f seconds\n",
			(double)scan_usecs / 1000000.0);
	fprintf(stderr, "\n");

	fprintf(stderr, "    General:\n");
	fprintf(stderr, "        Reads Matched:          %s    "
			"(%.4f%%)\n", comma_integer(reads_matched),
			(nreads == 0) ? 0 : ((double)reads_matched / (double)nreads) * 100);
	fprintf(stderr, "        Total Matches:          %s\n",
			comma_integer(total_matches));
	fprintf(stderr, "        Avg Hits/Matched Read:  %.2f\n",
			(total_matches == 0) ? 0 : ((double)total_matches /
					(double)reads_matched));
	fprintf(stderr, "        Duplicate Hits Pruned:  %s\n",
			comma_integer(dup_hits_c));

	if (Mflag) {
		fprintf(stderr, "\n    Memory usage:\n");
		fprintf(stderr, "        Genomemap:                %s\n",
				comma_integer(count_get_count(&mem_genomemap)));
	}
}

void usage(char *progname,bool full_usage){
	char *slash;
	double vect_sw_threshold = -1;
	uint sn;

	switch (shrimp_mode) {
	case MODE_LETTER_SPACE:
		vect_sw_threshold = DEF_SW_VECT_THRESHOLD;
		break;
	case MODE_COLOUR_SPACE:
		vect_sw_threshold = DEF_SW_VECT_THRESHOLD;
		break;
	case MODE_HELICOS_SPACE:
		fprintf(stderr,"Helicos mode not supported");
		exit(1);
		break;
	default:
		assert(0);
	}

	load_default_seeds();

	slash = strrchr(progname, '/');
	if (slash != NULL)
		progname = slash + 1;

	fprintf(stderr, "usage: %s [parameters] [options] "
			"reads_file genome_file1 genome_file2...\n", progname);

	fprintf(stderr, "Parameters:\n");

	fprintf(stderr,
			"    -s    Spaced Seed(s)                          (default: ");
	for (sn = 0; sn < n_seeds; sn++) {
		if (sn > 0)
			fprintf(stderr, "                                                            ");
		fprintf(stderr, "%s%s\n", seed_to_string(sn), (sn == n_seeds - 1? ")" : ","));
	}

	fprintf(stderr,
			"    -n    Seed Matches per Window                 (default: %d)\n",
			DEF_NUM_MATCHES);

	fprintf(stderr,
			"    -w    Seed Window Length                      (default: %.02f%%)\n",
			DEF_WINDOW_LEN);

	fprintf(stderr,
			"    -W    Seed Window Overlap Length              (default: %.02f%%)\n",
			DEF_WINDOW_OVERLAP);

	fprintf(stderr,
			"    -o    Maximum Hits per Read                   (default: %d)\n",
			DEF_NUM_OUTPUTS);

	fprintf(stderr, "\n");
	fprintf(stderr,
			"    -m    S-W Match Value                         (default: %d)\n",
			DEF_MATCH_VALUE);

	fprintf(stderr,
			"    -i    S-W Mismatch Value                      (default: %d)\n",
			DEF_MISMATCH_VALUE);

	fprintf(stderr,
			"    -g    S-W Gap Open Penalty (Reference)        (default: %d)\n",
			DEF_A_GAP_OPEN);

	fprintf(stderr,
			"    -q    S-W Gap Open Penalty (Query)            (default: %d)\n",
			DEF_B_GAP_OPEN);

	fprintf(stderr,
			"    -e    S-W Gap Extend Penalty (Reference)      (default: %d)\n",
			DEF_A_GAP_EXTEND);

	fprintf(stderr,
			"    -f    S-W Gap Extend Penalty (Query)          (default: %d)\n",
			DEF_B_GAP_EXTEND);

	if (shrimp_mode == MODE_COLOUR_SPACE) {
		fprintf(stderr,
				"    -x    S-W Crossover Penalty                   ("
				"default: %d)\n", DEF_XOVER_PENALTY);
	}

	if (shrimp_mode == MODE_COLOUR_SPACE) {
		fprintf(stderr,
				"    -h    S-W Full Hit Threshold                  "
				"(default: %.02f%%)\n", DEF_SW_FULL_THRESHOLD);
	}
	if (shrimp_mode == MODE_COLOUR_SPACE || shrimp_mode == MODE_HELICOS_SPACE) {
		fprintf(stderr,
				"    -v    S-W Vector Hit Threshold                "
				"(default: %.02f%%)\n", vect_sw_threshold);
	} else {
		fprintf(stderr,
				"    -h    S-W Hit Threshold                       "
				"(default: %.02f%%)\n", DEF_SW_FULL_THRESHOLD);
	}

	if (full_usage) {
		fprintf(stderr, "\n");
		if (shrimp_mode == MODE_COLOUR_SPACE || shrimp_mode == MODE_LETTER_SPACE) {
			fprintf(stderr,
					"    -A    Anchor width limiting full SW           (default: %d; disable: -1)\n",
					DEF_ANCHOR_WIDTH);
		}

	}
	fprintf(stderr,"\n");
	fprintf(stderr,"    -N    Set the number of threads               (default: %u)\n",DEF_NUM_THREADS);
	fprintf(stderr,"    -K    Set the thread chunk size               (default: %u)\n",DEF_CHUNK_SIZE);

	fprintf(stderr, "\n");
	fprintf(stderr, "Options:\n");

	fprintf(stderr,
			"    -B    Print Scan Progress Bar                       (default: "
			"disabled)\n");

	fprintf(stderr,
			"    -C    Only Process Negative Strand (Rev. Compl.)    (default: "
			"disabled)\n");

	fprintf(stderr,
			"    -F    Only Process Positive Strand                  (default: "
			"disabled)\n");

	fprintf(stderr,
			"    -P    Pretty Print Alignments                       (default: "
			"disabled)\n");

	fprintf(stderr,
			"    -R    Print Reads in Output                         (default: "
			"disabled)\n");

	if (shrimp_mode != MODE_HELICOS_SPACE) {
		fprintf(stderr,
				"    -T    Reverse tie-break choice on negative strand   (default: "
				"disabled)\n");
	}

	fprintf(stderr,
			"    -U    Print Unmapped Read Names in Output           (default: "
			"disabled)\n");
	fprintf(stderr,
			"    -M    Toggle brief memory usage statistics          (default: %s)\n",
			Mflag? "enabled" : "disabled");

	fprintf(stderr,
				"    -M    Toggle thread statistics                  (default: %s)\n",
				Dflag? "enabled" : "disabled");
	fprintf(stderr,
			"    -?    Full list of parameters and options\n");


	exit(1);
}

void print_settings() {
	static char const my_tab[] = "    ";
	uint i;

	fprintf(stderr, "Settings:\n");
	fprintf(stderr, "%s%-40s%s (%u/%u)\n", my_tab,
			(n_seeds == 1) ? "Spaced Seed (weight/span)" : "Spaced Seeds (weight/span)",
					seed_to_string(0), seed[0].weight, seed[0].span);
	for (i = 1; i < n_seeds; i++) {
		fprintf(stderr, "%s%-40s%s (%u/%u)\n", my_tab, "",
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
	fprintf(stderr, "\n");
	fprintf(stderr,"%s%-40s%u\n",my_tab,"Number of threads:",num_threads);
	fprintf(stderr,"%s%-40s%u\n",my_tab,"Thread chuck size:",chunk_size);

}


#ifdef DEBUG_MAIN
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
			load_genome(argv+1,1);
			print_genomemap_stats();
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
			launch_scan_threads(argv[2]);
		}
		if (0){
			fprintf(stderr,"loading compressed index\n");
			load_genome_map("testfile.gz");
		}
		fprintf(stderr,"done\n");
	}
	print_run_stats();
}

#else
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

	set_mode_from_argv(argv);

	a_gap_open_set = b_gap_open_set = a_gap_extend_set = b_gap_extend_set = false;

	set_mode_from_argv(argv);

	fprintf(stderr, "--------------------------------------------------"
			"------------------------------\n");
	fprintf(stderr, "gmapper: %s.\nSHRiMP %s\n[%s]\n", get_mode_string(),
			SHRIMP_VERSION_STRING, get_compiler());
	fprintf(stderr, "--------------------------------------------------"
			"------------------------------\n");

	//TODO -t -9 -d -Z -D -Y
	switch(shrimp_mode){
	case MODE_COLOUR_SPACE:
		optstr = "?s:n:w:o:r:m:i:g:q:e:f:x:h:v:N:K:BCFHPRTUA:MW:S:L:D";
		break;
	case MODE_LETTER_SPACE:
		optstr = "?s:n:w:o:r:m:i:g:q:e:f:h:X:N:K:BCFHPRTUA:MW:S:L:D";
		break;
	case MODE_HELICOS_SPACE:
		fprintf(stderr,"Helicose currently unsuported\n");
		exit(1);
		break;
	default:
		assert(0);
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
		case 'r':
			/* max_read_len is gone - silently ignore */
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
			assert(shrimp_mode != MODE_HELICOS_SPACE);
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
			assert(shrimp_mode == MODE_COLOUR_SPACE ||
					shrimp_mode == MODE_HELICOS_SPACE);
			sw_vect_threshold = atof(optarg);
			if (sw_vect_threshold < 0.0) {
				fprintf(stderr, "error: invalid s-w vector "
						"hit threshold\n");
				exit(1);
			}
			if (strcspn(optarg, "%.") == strlen(optarg))
				sw_vect_threshold = -sw_vect_threshold;	//absol.
			break;
		case 'B':
			Bflag = true;
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
		case 'U':
			Uflag = true;
			break;
			/*
			 * New options/parameters since SHRiMP 1.2.1
			 */
		case 'A':
			anchor_width = atoi(optarg);
			if (anchor_width < -1 || anchor_width >= 100) {
				fprintf(stderr, "error: anchor_width requested is invalid (%s)\n",
						optarg);
				exit(1);
			}
			break;
		case 'M':
			Mflag = !Mflag;
			break;
		case 'W':
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
		default:
			usage(progname, false);
		}
	}

	argc -= optind;
	argv += optind;

	if(load_file != NULL && n_seeds != 0){
		fprintf(stderr,"error: cannot specify seeds when loading genome map\n");
		usage(progname,false);
	}

	if (n_seeds == 0)
		load_default_seeds();

	kmer_to_mapidx = kmer_to_mapidx_orig;
	if (Hflag){
		kmer_to_mapidx = kmer_to_mapidx_hash;
		init_seed_hash_mask();
	}

	if (save_file != NULL && load_file != NULL){
		fprintf(stderr,"error: -L and -S incompatible\n");
		exit(1);
	}

	if(load_file != NULL){
		if (argc == 0){
			fprintf(stderr,"error: read_file not specified\n");
			usage(progname,false);
		}
		if (argc == 1){
			reads_file    = argv[0];
		} else {
			fprintf(stderr,"error: too many arguments with -L\n");
			usage(progname,false);
		}
	} else if (save_file != NULL){
		if (argc == 0){
			fprintf(stderr, "error: genome_file(s) not specified\n");
			usage(progname,false);
		}
		genome_files  = &argv[0];
		ngenome_files = argc;
	} else {
		if (argc < 2) {
			fprintf(stderr, "error: %sgenome_file(s) not specified\n",
					(argc == 0) ? "reads_file, " : "");
			usage(progname, false);
		}

		reads_file    = argv[0];
		genome_files  = &argv[1];
		ngenome_files = argc - 1;
	}
	if (!Cflag && !Fflag)
		Cflag = Fflag = true;

	if (shrimp_mode == MODE_LETTER_SPACE)
		sw_vect_threshold = sw_full_threshold;

	if (!valid_spaced_seeds()) {
		fprintf(stderr, "error: invalid spaced seed\n");
		if (!Hflag)
			fprintf(stderr, "       for longer seeds, try using the -H flag\n");
		exit(1);
	}

	if (!IS_ABSOLUTE(window_len) && window_len < 100.0) {
		fprintf(stderr, "error: window length < 100%% "
				"of read length\n");
		exit(1);
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

	if (shrimp_mode == MODE_COLOUR_SPACE && !IS_ABSOLUTE(sw_vect_threshold) &&
			sw_vect_threshold > 100.0) {
		fprintf(stderr, "error: invalid s-w full hit threshold\n");
		exit(1);
	}

	if ((a_gap_open_set && !b_gap_open_set)
			|| (a_gap_extend_set && !b_gap_extend_set))
		fputc('\n', stderr);
	if (a_gap_open_set && !b_gap_open_set) {
		fprintf(stderr, "Notice: Gap open penalty set for reference but not query; assuming symmetry.\n");
		b_gap_open_score = a_gap_open_score;
	}
	if (a_gap_extend_set && !b_gap_extend_set) {
		fprintf(stderr, "Notice: Gap extend penalty set for reference but not query; assuming symmetry.\n");
		b_gap_extend_score = a_gap_extend_score;
	}
	if ((a_gap_open_set && !b_gap_open_set)
			|| (a_gap_extend_set && !b_gap_extend_set))
		fputc('\n', stderr);

	print_settings();

	uint64_t before;
	before = gettimeinusecs();
	if (load_file != NULL){
		if (!load_genome_map(load_file)){
			exit(1);
		}
	} else {
		if (!load_genome(genome_files,ngenome_files)){
			exit(1);
		}
	}
	map_usecs += (gettimeinusecs() - before);

	if (save_file != NULL){
		if(save_genome_map(save_file)){
			exit(0);
		}
		exit(1);
	}

		//TODO setup need max window and max read len
		int longest_read_len = 100;
		int max_window_len = (int)abs_or_pct(window_len,longest_read_len);
#pragma omp parallel shared(longest_read_len,max_window_len,a_gap_open_score, a_gap_extend_score, b_gap_open_score, b_gap_extend_score,\
		match_score, mismatch_score,shrimp_mode,crossover_score,anchor_width) num_threads(num_threads)
		{
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

		output = output_format_line(Rflag);
		puts(output);
		free(output);

		before = gettimeinusecs();
		launch_scan_threads(reads_file);
		scan_usecs += (gettimeinusecs() - before);

		print_statistics();
}
#endif
//TODO ask matei about Fflag Cflag, setting up sw,
