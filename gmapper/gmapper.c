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
#include "../common/input.h"
#include "../common/heap.h"

DEF_HEAP(uint32_t,uint,uu)

/* Parameters */
static double	window_len		= DEF_WINDOW_LEN;
static double	window_overlap		= DEF_WINDOW_OVERLAP;
static uint	num_matches		= DEF_NUM_MATCHES;
static uint	num_outputs		= DEF_NUM_OUTPUTS;
static uint	num_tmp_outputs		= 20 + num_outputs;
static bool	hash_filter_calls	= DEF_HASH_FILTER_CALLS;
static int	anchor_width		= DEF_ANCHOR_WIDTH;
static bool	gapless_sw		= DEF_GAPLESS_SW;
static uint32_t	list_cutoff		= DEF_LIST_CUTOFF;

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
static int Bflag = false;			/* print a progress bar */
static int Cflag = false;			/* do complement only */
static int Fflag = false;			/* do positive (forward) only */
static int Hflag = false;			/* use hash table, not lookup */
static int Pflag = false;			/* pretty print results */
static int Rflag = false;			/* add read sequence to output*/
static int Tflag = false;			/* reverse sw full tie breaks */
static int Uflag = false;			/* output unmapped reads, too */
static int Mflag = false;			/* print insert histogram */
static int Dflag = false;			/* print statistics for each thread */
static int Eflag = false;			/* output sam format */

/* Mate Pairs */
static int	pair_mode		= DEF_PAIR_MODE;
static int	min_insert_size		= DEF_MIN_INSERT_SIZE;
static int	max_insert_size		= DEF_MAX_INSERT_SIZE;
static uint64_t	insert_histogram[100];
static int	insert_histogram_bucket_size;

/* Statistics */
static uint64_t nreads;

static uint64_t total_reads_matched;
static uint64_t total_pairs_matched;

static uint64_t total_single_matches;
static uint64_t total_paired_matches;

static uint64_t	total_dup_single_matches;			/* number of duplicate hits */
static uint64_t total_dup_paired_matches;


static uint64_t total_work_usecs;
static uint64_t map_usecs;

static count_t mem_genomemap;

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

/* constants for thread control */
static uint num_threads = DEF_NUM_THREADS;
static uint chunk_size = DEF_CHUNK_SIZE;

static uint64_t scan_ticks[50];
static uint64_t wait_ticks[50];


/* Thread-private */
//static uint32_t hash_mark;
//struct window_cache_entry {
//  uint32_t mark;
//  uint32_t score;
//};
//static struct window_cache_entry * window_cache;

//#pragma omp threadprivate(hash_mark, window_cache)


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

/* get contig number from absolute index */
static inline void
get_contig_num(uint idx, uint * cn) {
	*cn = 0;
	while (*cn < num_contigs - 1
			&& idx >= contig_offsets[*cn + 1])
		(*cn)++;

	assert(contig_offsets[*cn] <= idx && idx < contig_offsets[*cn] + genome_len[*cn]);
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


/*
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
*/

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


#ifdef DEBUGGING
static void print_everything(){
	uint sn;
	fprintf(stderr,"----------------------------------------\n");
	fprintf(stderr,"%s\n",get_mode_string());
	for (sn = 0; sn < n_seeds; sn++){
		fprintf(stderr,"seed %u: %s (%u/%u)\n",sn,seed_to_string(sn),seed[sn].weight,seed[sn].span);
		uint capacity;
		if(Hflag){
			capacity = power(4, HASH_TABLE_POWER);
		} else {
			capacity = power(4, seed[sn].weight);
		}
		uint mapidx;
		for(mapidx = 0; mapidx < capacity; mapidx ++){
			fprintf(stderr,"genomemap_len[%u][%u] = %u\n",sn,mapidx,genomemap_len[sn][mapidx]);
			uint j;
			for(j = 0; j < genomemap_len[sn][mapidx]; j++){
				fprintf(stderr,"%u ",genomemap[sn][mapidx][j]);
			}
			fprintf(stderr,"\n");
		}
	}
	fprintf(stderr,"----------------------------------------\n");
	fprintf(stderr,"nkemers = %u\n",nkmers);
	fprintf(stderr,"is rna? %s\n",genome_is_rna ? "yes":"no");
	uint c;
	for(c = 0; c < num_contigs; c++){
		fprintf(stderr,"%s\n",contig_names[c]);
		fprintf(stderr,"offset = %u, len = %u\n",contig_offsets[c],genome_len[c]);
		uint j;
		for(j = 0; j < genome_len[c];j++){
			fprintf(stderr,"%c",base_translate(EXTRACT(genome_contigs[c],j),false));
		}
		fprintf(stderr,"\n");
		for(j = 0; j < genome_len[c];j++){
			fprintf(stderr,"%c",base_translate(EXTRACT(genome_contigs_rc[c],j),false));
		}
		fprintf(stderr,"\n");
	}
	fprintf(stderr,"----------------------------------------\n");
}
#endif


static bool save_genome_map_seed(const char *file,uint sn){
	uint32_t total = 0;
	uint32_t i;

	gzFile fp = gzopen(file, "wb");
	if (fp == NULL){
		return false;
	}
	// shrimp_mode
	uint32_t m;
	m = (uint32_t)(shrimp_mode);
	xgzwrite(fp,&m,sizeof(uint32_t));
	// Hflag
	uint32_t h = (uint32_t)Hflag;
	xgzwrite(fp,&h,sizeof(uint32_t));
	// Seed
	xgzwrite(fp,&seed[sn], sizeof(seed_type));
	// genomemap_len
	size_t capacity;
	if(Hflag){
		capacity = power(4, HASH_TABLE_POWER);
	} else {
		capacity = power(4, seed[sn].weight);
	}
	xgzwrite(fp,genomemap_len[sn],sizeof(uint32_t) * capacity);

	//total
	for (i = 0; i < capacity;i++){
		total += genomemap_len[sn][i];
	}
	xgzwrite(fp,&total,sizeof(uint32_t)); //TODO do not need to write this but makes things simpler

	// genome_map
	// TODO if memory usage to high change this
	uint32_t *map,*ptr;
	map = (uint32_t *)xmalloc(sizeof(uint32_t)*total);
	ptr = map;
	for (i = 0;i<capacity;i++){
		memcpy(ptr,genomemap[sn][i],sizeof(uint32_t)*genomemap_len[sn][i]);
		ptr += genomemap_len[sn][i];
	}
	xgzwrite(fp,map,sizeof(uint32_t)*total);

	gzclose(fp);
	free(map);
	return true;

}

static bool load_genome_map_seed(const char *file){
	uint32_t i;
	uint32_t total;

	gzFile fp = gzopen(file, "rb");
	if (fp == NULL){
		fprintf(stderr,"Could not open file [%s]\n",file);
		return false;
	}

	// shrimp_mode
	uint32_t m;
	xgzread(fp,&m,sizeof(uint32_t));
	if(m != (uint32_t)shrimp_mode){
		fprintf(stderr,"Shrimp mode in file %s does not match\n",file);
	}
	// Hflag
	uint32_t h;
	xgzread(fp,&h,sizeof(uint32_t));
	if (h != (uint32_t)Hflag){
		fprintf(stderr,"Hash settings do not match in file %s\n",file);
	}
	// Seed
	uint sn = n_seeds;
	n_seeds++;
	seed = (seed_type *)xrealloc(seed,sizeof(seed_type)*n_seeds);
	genomemap_len = (uint32_t **)xrealloc(genomemap_len,sizeof(uint32_t *)*n_seeds);
	genomemap = (uint32_t ***)xrealloc(genomemap,sizeof(uint32_t **)*n_seeds);
	xgzread(fp,seed + sn,sizeof(seed_type));

	if (seed[sn].span > max_seed_span)
		max_seed_span = seed[sn].span;

	if (seed[sn].span < min_seed_span)
		min_seed_span = seed[sn].span;

	avg_seed_span = 0;
	for(i =0; i < n_seeds;i++){
		avg_seed_span += seed[i].span;
	}
	avg_seed_span = avg_seed_span/n_seeds;

	// genomemap_len
	size_t capacity;
	if(Hflag){
		capacity = power(4, HASH_TABLE_POWER);
	} else {
		capacity = power(4, seed[sn].weight);
	}
	genomemap_len[sn] = (uint32_t*)xmalloc(sizeof(uint32_t)*capacity);
	genomemap[sn] = (uint32_t **)xmalloc(sizeof(uint32_t *)*capacity);
	xgzread(fp,genomemap_len[sn],sizeof(uint32_t) * capacity);

	// total
	xgzread(fp,&total,sizeof(uint32_t)); //TODO do not need to write this but makes things simpler
	nkmers += total;
	// genome_map
	uint32_t * map;
	map = (uint32_t *)xmalloc(sizeof(uint32_t)*total);
	gzread(fp,map,sizeof(uint32_t)*total);
	uint32_t * ptr;
	ptr = map;

	for (i = 0; i < capacity;i++){
		genomemap[sn][i] = ptr;
		ptr += genomemap_len[sn][i];
	}

	gzclose(fp);
	return true;
}

static bool save_genome_map(const char *prefix) {
	char *name;
	name = (char *)xmalloc(strlen(prefix)+n_seeds+10);
	uint sn;
	for(sn = 0;sn < n_seeds;sn++){
		sprintf(name,"%s.seed.%u",prefix,sn);
		save_genome_map_seed(name,sn);
	}

	sprintf(name,"%s.genome",prefix);
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

	uint i;
	uint32_t total = 0;
	for(i = 0; i < num_contigs; i++){
		uint32_t len = (uint32_t)strlen(contig_names[i]);
		xgzwrite(fp,&len,sizeof(uint32_t));
		xgzwrite(fp,contig_names[i],len +1);
		total += BPTO32BW(genome_len[i]);
	}

	xgzwrite(fp,&total,sizeof(uint32_t));

	//genome_contigs / genome_contigs_rc / genome_cs_contigs / genome_initbp
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

	gzclose(fp);
	return true;
	//	gzFile fp = gzopen(file, "wb0");
	//	if (fp == NULL){
	//		return false;
	//	}
	//
	//	//write the header
	//	//fprintf(stderr,"saving num_contgs: %u\n",num_contigs);
	//	gzwrite(fp,&shrimp_mode,sizeof(shrimp_mode_t));
	//	fprintf(stderr,"shrimp_mode:%u\n",shrimp_mode);
	//	uint32_t h = (uint32_t)Hflag;
	//	fprintf(stderr,"h:%u\n",h);
	//	gzwrite(fp,&h,sizeof(uint32_t));
	//	uint32_t i;
	//	gzwrite(fp,&num_contigs,sizeof(uint32_t));
	//	//fprintf(stderr,"saved num_contigs\n");
	//	for (i = 0; i < num_contigs; i++) {
	//		fprintf(stderr,"offset:%u\n",contig_offsets[i]);
	//		gzwrite(fp,&contig_offsets[i],sizeof(uint32_t));
	//		size_t len = strlen(contig_names[i]);
	//		fprintf(stderr,"len:%u\n",len);
	//		gzwrite(fp,&len,sizeof(size_t));
	//		gzwrite(fp,contig_names[i],len +1);
	//
	//		gzwrite(fp,genome_len + i,sizeof(uint32_t));
	//		gzwrite(fp,genome_contigs[i],sizeof(uint32_t)*BPTO32BW(genome_len[i]));
	//		gzwrite(fp,genome_contigs_rc[i],sizeof(uint32_t)*BPTO32BW(genome_len[i]));
	//		if(shrimp_mode == MODE_COLOUR_SPACE){
	//			gzwrite(fp,genome_cs_contigs[i],sizeof(uint32_t)*BPTO32BW(genome_len[i]));
	//			gzwrite(fp,genome_initbp+i,sizeof(uint32_t));
	//		}
	//	}
	//	//fprintf(stderr,"saving seeds\n");
	//	//write the seeds and genome_maps
	//	gzwrite(fp,&n_seeds,sizeof(uint32_t));
	//	for (i = 0; i < n_seeds; i++) {
	//		gzwrite(fp,&seed[i], sizeof(seed_type));
	//
	//		//write the genome_map for this seed
	//		uint32_t j;
	//		//uint32_t p = power(4, seed[i].weight);
	//		//fprintf(stderr,"saving index\n");
	//		size_t capacity;
	//		if(Hflag){
	//			capacity = sizeof(uint32_t *) * power(4, HASH_TABLE_POWER);
	//		} else {
	//			capacity = power(4, seed[i].weight);
	//		}
	//		gzwrite(fp,&capacity, sizeof(size_t));
	//		for (j = 0; j < capacity; j++) {
	//			uint32_t len = genomemap_len[i][j];
	//			gzwrite(fp, &len, sizeof(uint32_t));
	//			gzwrite(fp, genomemap[i][j], sizeof(uint32_t) * len);
	//		}
	//	}
	//	gzclose(fp);
	//	return true;
}


static bool load_genome_map(const char *file){
	gzFile fp = gzopen(file,"rb");
	if (fp == NULL){
		return false;
	}

	//shrimp mode
	uint32_t m;
	xgzread(fp,&m,sizeof(uint32_t));
	if (shrimp_mode != (shrimp_mode_t)m) {
	  fprintf(stderr, "error: shrimp mode does not match genome file (%s)\n", file);
	  exit(1);
	}

	//Hflag
	uint32_t h;
	xgzread(fp,&h,sizeof(uint32_t));
	Hflag = h;

	// num_contigs
	xgzread(fp,&num_contigs,sizeof(uint32_t));

	genome_len = (uint32_t *)xmalloc(sizeof(uint32_t)*num_contigs);
	contig_offsets = (uint32_t *)xmalloc(sizeof(uint32_t)*num_contigs);
	contig_names = (char **)xmalloc(sizeof(char *)*num_contigs);

	genome_contigs = (uint32_t **)xmalloc(sizeof(uint32_t *)*num_contigs);
	genome_contigs_rc = (uint32_t **)xmalloc(sizeof(uint32_t *)*num_contigs);
	if(shrimp_mode == MODE_COLOUR_SPACE){
		genome_cs_contigs = (uint32_t **)xmalloc(sizeof(uint32_t *)*num_contigs);
		genome_initbp = (uint32_t *)xmalloc(sizeof(uint32_t)*num_contigs);
	}

	//genome_len
	xgzread(fp,genome_len,sizeof(uint32_t)*num_contigs);
	uint i;

	// contig_offfsets
	xgzread(fp,contig_offsets,sizeof(uint32_t)*num_contigs);

	// names / total

	for(i = 0; i < num_contigs; i++){
		uint32_t len;
		xgzread(fp,&len,sizeof(uint32_t));

		contig_names[i] = (char *)xmalloc(sizeof(char)*(len+1));
		xgzread(fp,contig_names[i],len+1);
	}

	uint32_t total;
	xgzread(fp,&total,sizeof(uint32_t));

	//genome_contigs / genome_contigs_rc / genome_cs_contigs / genome_initbp
	uint32_t *gen, *gen_rc, *gen_cs, *ptr1, *ptr2, *ptr3 = NULL;
	gen = (uint32_t *)xmalloc(sizeof(uint32_t)*total);
	xgzread(fp,gen,sizeof(uint32_t)*total);
	ptr1 = gen;
	gen_rc = (uint32_t *)xmalloc(sizeof(uint32_t)*total);
	xgzread(fp,gen_rc,sizeof(uint32_t)*total);
	ptr2 = gen_rc;
	if(shrimp_mode == MODE_COLOUR_SPACE){
		gen_cs = (uint32_t *)xmalloc(sizeof(uint32_t)*total);
		xgzread(fp,gen_cs,sizeof(uint32_t)*total);
		ptr3 = gen_cs;
		xgzread(fp,genome_initbp,sizeof(uint32_t)*num_contigs);
	}
	for(i = 0; i < num_contigs; i++){
		genome_contigs[i] = ptr1;
		ptr1 += BPTO32BW(genome_len[i]);
		genome_contigs_rc[i] = ptr2;
		ptr2 += BPTO32BW(genome_len[i]);
		if(shrimp_mode == MODE_COLOUR_SPACE){
			genome_cs_contigs[i] = ptr3;
			ptr3 += BPTO32BW(genome_len[i]);
		}
	}
	gzclose(fp);
	return true;
	//	gzFile fp = gzopen(file,"rb");
	//	if (fp == NULL){
	//		return false;
	//	}
	//
	//	uint32_t i;
	//	gzread(fp,&shrimp_mode,sizeof(shrimp_mode_t)); //TODO make sure no conflict
	//	uint32_t f;
	//	gzread(fp,&f,sizeof(uint32_t));
	//	Hflag = f;
	//	gzread(fp,&num_contigs,sizeof(uint32_t));
	//	//fprintf(stderr,"num_contigs = %u\n",num_contigs);
	//
	//	contig_names = (char **)xrealloc(contig_names,sizeof(char *)*num_contigs);
	//	contig_offsets = (uint32_t *)xrealloc(contig_offsets,sizeof(uint32_t)*num_contigs);
	//
	//	genome_len = (uint32_t *)xrealloc(genome_len,sizeof(uint32_t)*num_contigs);
	//	genome_contigs = (uint32_t **)xrealloc(genome_contigs,sizeof(uint32_t *)*num_contigs);
	//	genome_contigs_rc = (uint32_t **)xrealloc(genome_contigs,sizeof(uint32_t *)*num_contigs);
	//	if(shrimp_mode == MODE_COLOUR_SPACE){
	//		genome_cs_contigs = (uint32_t **)xrealloc(genome_cs_contigs,sizeof(uint32_t *)*num_contigs);
	//		genome_initbp = (uint32_t *)xrealloc(genome_initbp,sizeof(uint32_t)*num_contigs);
	//	}
	//
	//	for (i = 0; i < num_contigs; i++){
	//		gzread(fp,&contig_offsets[i],sizeof(uint32_t));
	//
	//		uint32_t len;
	//		gzread(fp,&len,sizeof(uint32_t));
	//		contig_names[i] = (char *)xrealloc(contig_names[i],sizeof(char)*len);
	//		gzread(fp,contig_names[i],len+1);
	//
	//		gzread(fp,genome_len + i,sizeof(uint32_t));
	//
	//		genome_contigs[i] = (uint32_t *)xrealloc(genome_contigs[i],sizeof(uint32_t)*BPTO32BW(genome_len[i]));
	//		gzread(fp,genome_contigs[i],sizeof(uint32_t)*BPTO32BW(genome_len[i]));
	//
	//		genome_contigs_rc[i] = (uint32_t *)xrealloc(genome_contigs_rc[i],sizeof(uint32_t)*BPTO32BW(genome_len[i]));
	//		gzread(fp,genome_contigs_rc[i],sizeof(uint32_t)*BPTO32BW(genome_len[i]));
	//
	//		if(shrimp_mode == MODE_COLOUR_SPACE){
	//			genome_cs_contigs[i] = (uint32_t *)xrealloc(genome_cs_contigs[i],sizeof(uint32_t)*BPTO32BW(genome_len[i]));
	//			gzread(fp,genome_cs_contigs[i],sizeof(uint32_t)*BPTO32BW(genome_len[i]));
	//			gzread(fp,genome_initbp+i,sizeof(uint32_t));
	//		}
	//	}
	//	gzread(fp,&n_seeds,sizeof(uint32_t));
	//	//fprintf(stderr,"n_seeds = %u\n",n_seeds);
	//	seed =(seed_type *)xrealloc(seed,sizeof(seed_type)*n_seeds);
	//	genomemap_len = (uint32_t **)xrealloc(genomemap_len,sizeof(uint32_t *)*n_seeds);
	//	genomemap = (uint32_t ***)xrealloc(genomemap,sizeof(uint32_t **) * n_seeds);
	//	for (i = 0; i < n_seeds; i++){
	//		gzread(fp,&seed[i], sizeof(seed_type));
	//		//fprintf(stderr,"seed %u: span=%u\n",i,seed[i].span);
	//		uint32_t j;
	//		size_t capacity;
	//		if(Hflag){
	//			capacity = sizeof(uint32_t *) * power(4, HASH_TABLE_POWER);
	//		} else {
	//			capacity = power(4, seed[i].weight);
	//		}
	//		gzread(fp,&capacity,sizeof(size_t));
	//		genomemap_len[i] = (uint32_t *)xrealloc(genomemap_len[i],sizeof(uint32_t)*capacity);
	//		genomemap[i] = (uint32_t **)xrealloc(genomemap[i],sizeof(uint32_t *)*capacity);
	//		for (j = 0; j < capacity; j++){
	//			gzread(fp,&genomemap_len[i][j],sizeof(uint32_t));
	//			genomemap[i][j] = (uint32_t *)xrealloc(genomemap[i][j],
	//					sizeof(uint32_t)*genomemap_len[i][j]);
	//			gzread(fp,genomemap[i][j],sizeof(uint32_t) * genomemap_len[i][j]);
	//		}
	//	}
	//	gzclose(fp);
	//	return true;
}


static void
read_get_mapidxs_per_strand(struct read_entry * re, uint st) {
  uint sn, i, load, base, r_idx;
  uint32_t * kmerWindow = (uint32_t *)xcalloc(sizeof(kmerWindow[0])*BPTO32BW(max_seed_span));

  re->mapidx[st] = (uint32_t *)xmalloc(n_seeds * re->max_n_kmers * sizeof(re->mapidx[0][0]));
  if (re->has_Ns) {
    re->mapidx_pos[st] = (bool *)xmalloc(n_seeds * re->max_n_kmers * sizeof(re->mapidx_pos[0][0]));
  }

  for (load = 0, i = 0; i < re->read_len; i++) {
    base = EXTRACT(re->read[st], i);
    bitfield_prepend(kmerWindow, max_seed_span, base);

    //skip past any Ns or Xs
    if (re->has_Ns && (base == BASE_N || base == BASE_X))
      load = 0;
    else if (load < max_seed_span)
      load++;

    for (sn = 0; sn < n_seeds; sn++) {
      if (i < re->min_kmer_pos + seed[sn].span - 1)
	continue;

      r_idx = i - seed[sn].span + 1;
      if (re->has_Ns && load < seed[sn].span) {
	re->mapidx_pos[st][sn*re->max_n_kmers + (r_idx - re->min_kmer_pos)] = false;
	continue;
      }

      re->mapidx[st][sn*re->max_n_kmers + (r_idx - re->min_kmer_pos)] = kmer_to_mapidx(kmerWindow, sn);
      if (re->has_Ns) {
	re->mapidx_pos[st][sn*re->max_n_kmers + (r_idx - re->min_kmer_pos)] = true;
      }
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
uw_anchor_cmp(void const * p1, void const * p2)
{
  if (((struct uw_anchor *)p1)->x < ((struct uw_anchor *)p2)->x)
    return -1;
  else if (((struct uw_anchor *)p1)->x > ((struct uw_anchor *)p2)->x)
    return 1;
  else
    return 0;
}


static uint
bin_search(uint32_t * array, uint l, uint r, uint32_t value)
{
  uint m;
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
read_get_restricted_anchor_list_per_strand(struct read_entry * re, uint st, bool collapse) {
  uint i, j, k, sn, offset;
  uint g_start, g_end;
  uint idx_start, idx_end;
  uint32_t mapidx;
  int anchor_cache[re->read_len];
  uint diag;
  int l;

  assert(re->mapidx[st] != NULL);

  re->n_anchors[st] = 0;

  for (j = 0; j < re->n_ranges; j++) {
    if (re->ranges[j].st != st)
      continue;

    g_start = contig_offsets[re->ranges[j].cn] + re->ranges[j].g_start;
    g_end = contig_offsets[re->ranges[j].cn] + re->ranges[j].g_end;

    for (sn = 0; sn < n_seeds; sn++) {
      for (i = 0; re->min_kmer_pos + i + seed[sn].span - 1 < re->read_len; i++) {
	offset = sn*re->max_n_kmers + i;
	mapidx = re->mapidx[st][offset];

	idx_start = bin_search(genomemap[sn][mapidx], 0, genomemap_len[sn][mapidx], g_start);
	idx_end = bin_search(genomemap[sn][mapidx], idx_start, genomemap_len[sn][mapidx], g_end + 1);

	if (idx_start >= idx_end)
	  continue;

	re->anchors[st] = (struct uw_anchor *)xrealloc(re->anchors[st],
						       (re->n_anchors[st] + (idx_end - idx_start))
						       * sizeof(re->anchors[0][0]));
	for (k = 0; idx_start + k < idx_end; k++) {
	  re->anchors[st][re->n_anchors[st] + k].cn = re->ranges[j].cn;
	  re->anchors[st][re->n_anchors[st] + k].x =
	    genomemap[sn][mapidx][idx_start + k] - contig_offsets[re->ranges[j].cn];
	  re->anchors[st][re->n_anchors[st] + k].y = re->min_kmer_pos + i;
	  re->anchors[st][re->n_anchors[st] + k].length = seed[sn].span;
	  re->anchors[st][re->n_anchors[st] + k].weight = 1;
	}
	re->n_anchors[st] += idx_end - idx_start;
      }
    }
  }

  qsort(re->anchors[st], re->n_anchors[st], sizeof(re->anchors[0][0]), uw_anchor_cmp);

  if (collapse) {
    for (i = 0; i < re->read_len; i++)
      anchor_cache[i] = -1;

    for (k = 0, i = 0; i < re->n_anchors[st]; i++) {
      re->anchors[st][k] = re->anchors[st][i];
      diag = (re->anchors[st][k].x + re->read_len - re->anchors[st][k].y) % re->read_len;
      l = anchor_cache[diag];
      if (l >= 0
	  && re->anchors[st][l].cn == re->anchors[st][k].cn
	  && uw_anchors_intersect(&re->anchors[st][l], &re->anchors[st][k])) {
	uw_anchors_join(&re->anchors[st][l], &re->anchors[st][k]);
      } else {
	anchor_cache[diag] = k;
	k++;
      }
    }
    re->n_anchors[st] = k;
  }
}


static void
read_get_anchor_list_per_strand(struct read_entry * re, uint st, bool collapse) {
  uint list_sz;
  uint i, sn, offset;
  struct heap_uu h;
  uint * idx;
  struct heap_uu_elem tmp;
  int anchor_cache[re->read_len];

  assert(re->mapidx[st] != NULL);

  // compute size of anchor list
  list_sz = 0;
  for (sn = 0; sn < n_seeds; sn++) {
    for (i = 0; re->min_kmer_pos + i + seed[sn].span - 1 < re->read_len; i++) {
      offset = sn*re->max_n_kmers + i;
      if (!re->has_Ns || re->mapidx_pos[st][offset]) {
	list_sz += genomemap_len[sn][re->mapidx[st][offset]];
      }
    }
  }

  // init anchor list
  re->anchors[st] = (struct uw_anchor *)xmalloc(list_sz * sizeof(re->anchors[0][0]));
  re->n_anchors[st] = 0;

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

      if ((!re->has_Ns || re->mapidx_pos[st][offset])
	  && idx[offset] < genomemap_len[sn][re->mapidx[st][offset]]
	  ) {
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
    re->anchors[st][re->n_anchors[st]].weight = 1;
    get_contig_num(re->anchors[st][re->n_anchors[st]].x, &re->anchors[st][re->n_anchors[st]].cn);
    re->n_anchors[st]++;

    if (collapse) {
      // check if current anchor intersects the cached one on the same diagonal
      uint diag = (re->anchors[st][re->n_anchors[st]-1].x + re->read_len - re->anchors[st][re->n_anchors[st]-1].y) % re->read_len;
      int j = anchor_cache[diag];

      if (j >= 0
	  && re->anchors[st][j].cn == re->anchors[st][re->n_anchors[st]-1].cn
	  && uw_anchors_intersect(&re->anchors[st][j], &re->anchors[st][re->n_anchors[st]-1])) {
	uw_anchors_join(&re->anchors[st][j], &re->anchors[st][re->n_anchors[st]-1]);
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
  fprintf(stderr, "(cn:%u,st:%u,gen_st:%u,g_off:%u,w_len:%u,scores:(%d,%d,%d),matches:%u,pair_min:%d,pair_max:%d,anchor:(%d,%d,%u,%u))\n",
	  hit->cn, hit->st, hit->gen_st, hit->g_off, hit->w_len,
	  hit->score_window_gen, hit->score_vector, hit->score_full,
	  hit->matches, hit->pair_min, hit->pair_max,
	  hit->anchor.x, hit->anchor.y, hit->anchor.length, hit->anchor.width);
}
#endif


#if defined (DEBUG_HIT_LIST_CREATION) || defined (DEBUG_HIT_LIST_PAIR_UP) \
  || defined (DEBUG_HIT_LIST_PASS1)
static void
dump_hit_list(struct read_entry * re, uint st, bool only_paired, bool only_after_vector)
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
read_get_hit_list_per_strand(struct read_entry * re, int match_mode, uint st) {
  uint cn, i, max_idx;
  uint32_t goff, w_len, gstart, gend;
  int max_score, tmp_score = 0;
  int j;
  int short_len = 0, long_len = 0;
  uint x_len;
  struct anchor a[3];

  re->hits[st] = (struct read_hit *)xcalloc(re->n_anchors[st] * sizeof(re->hits[0][0]));
  re->n_hits[st] = 0;

  for (i = 0; i < re->n_anchors[st]; i++) {
    // contig num of crt anchor
    cn = re->anchors[st][i].cn;

    // w_len
    w_len = re->window_len;
    if (w_len > genome_len[cn])
      w_len = genome_len[cn];

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

      for (j = (int)i - 1;
	   j >= 0
	     && re->anchors[st][j].x >= contig_offsets[cn] + gstart;
	   j--) {
	if (re->anchors[st][j].y >= re->anchors[st][i].y)
	  continue;

	if ((int64_t)(re->anchors[st][i].x - contig_offsets[cn]) - (int64_t)re->anchors[st][i].y
	    > (int64_t)(re->anchors[st][j].x - contig_offsets[cn]) - (int64_t)re->anchors[st][j].y)
	  { // deletion in read
	    short_len = (re->anchors[st][i].y - re->anchors[st][j].y) + re->anchors[st][i].length;
	    long_len = (re->anchors[st][i].x - re->anchors[st][j].x) + re->anchors[st][i].length;
	  }
	else
	  { // insertion in read
	    short_len = (re->anchors[st][i].x - re->anchors[st][j].x) + re->anchors[st][i].length;
	    long_len = (re->anchors[st][i].y - re->anchors[st][j].y) + re->anchors[st][i].length;
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
	x_len = (re->anchors[st][i].x - re->anchors[st][max_idx].x) + re->anchors[st][i].length;

	if ((re->window_len - x_len)/2 < re->anchors[st][max_idx].x - contig_offsets[cn])
	  goff = (re->anchors[st][max_idx].x - contig_offsets[cn]) - (re->window_len - x_len)/2;
	else
	  goff = 0;

	if (goff + w_len > genome_len[cn])
	  goff = genome_len[cn] - w_len;

	// compute anchor
	if (max_idx < i) {
	  uw_anchor_to_anchor(&re->anchors[st][i], &a[0], contig_offsets[cn] + goff);
	  uw_anchor_to_anchor(&re->anchors[st][max_idx], &a[1], contig_offsets[cn] + goff);
	  join_anchors(a, 2, &a[2]);
	} else {
	  uw_anchor_to_anchor(&re->anchors[st][i], &a[2], contig_offsets[cn] + goff);
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
    if (j < (int)i) { // shift elements at indexes j..i-1 higher
      struct read_hit tmp = re->hits[st][i];
      int k;
      for (k = (int)i-1; k >= j; k--)
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
read_pass1_per_strand(struct read_entry * re, bool only_paired, uint st) {
  uint i;
  int j;

  f1_hash_tag++;
  j = -1; // last good hit

  for (i = 0; i < re->n_hits[st]; i++) {
    if (only_paired && re->hits[st][i].pair_min < 0)
      continue;

    // check window overlap
    if (j >= 0
	&& re->hits[st][i].cn == re->hits[st][j].cn
	&& re->hits[st][i].g_off <= re->hits[st][j].g_off + (uint)abs_or_pct(window_overlap, re->window_len)) {
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
  uint st1, st2, i, j, k, l;
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
  uint st, i;
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
  uint st1, st2, i, j;
  heap_paired_elem tmp;

  assert(re1 != NULL && re2 != NULL && h != NULL);
  assert(pair_mode != PAIR_NONE);

  for (st1 = 0; st1 < 2; st1++) {
    st2 = 1 - st1; // opposite strand

    for (i = 0; i < re1->n_hits[st1]; i++) {
      if (re1->hits[st1][i].pair_min < 0)
	continue;

      for (j = re1->hits[st1][i].pair_min; (int)j <= re1->hits[st1][i].pair_max; j++) {
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
	  tmp.rest.insert_size = (st1 == 0?
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
  reverse_anchor(&rh->anchor, rh->w_len, re->read_len);
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
  fprintf(stderr, "SW full call: (name:[%s],cn:%u,st:%u,gen_st:%u,g_off:%u,w_len:%u,anchor:(%d,%d,%u,%u))\n",
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


/*
 * Print given hit.
 *
 */
static inline void
hit_output(struct read_entry * re, struct read_hit * rh,struct read_entry * re_mp, struct read_hit * rh_mp, char ** output1, char ** output2,bool paired)
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
    free(output);

    struct input inp, inp_mp;
    memset(&inp, 0, sizeof(inp));
    input_parse_string(*output1,fsp,&inp);
    if(re_mp != NULL){
    	free(*output1);
		*output1 = output_normal(re_mp->name, contig_names[rh_mp->cn], rh_mp->sfrp,
    			   genome_len[rh_mp->cn], shrimp_mode == MODE_COLOUR_SPACE, re_mp->read[rh_mp->st],
    			   re_mp->read_len, re_mp->initbp[rh_mp->st], rh_mp->gen_st, Rflag);
		input_parse_string(*output1,fsp,&inp_mp);
    }
    char * cigar, *read;
    uint32_t * read_bitstring, *read_ls_bitstring;
    cigar = (char *)xmalloc(sizeof(char)*200);
	int first_bp = 0;
    edit2cigar(inp.edit,inp.read_start,inp.read_end,inp.read_length,cigar);

	read_bitstring = re->read[rh->gen_st];

    if (shrimp_mode == COLOUR_SPACE){
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
    }
    read = readtostr(read_bitstring,re->read_len,false,0);
    char *name;
    bool first = false, second = false;
    if (re_mp == NULL){
    	name = inp.read;
    } else {
    	name = (char *)xmalloc(sizeof(char *)*strlen(inp.read)+1);
    	strncpy(name,inp.read,strlen(inp.read)+1);
    	char * end = strrchr(name,':');
    	*end = '\0';
    	end++;
    	if (*end == '1'){
    		first = true;
    	} else {
    		second = true;
    	}
    }

    int ins_size =0;
    if (re_mp != NULL){
    	if (pair_mode == PAIR_COL_FW || pair_mode == PAIR_COL_BW){
    		ins_size = inp_mp.genome_start - inp.genome_start;
    	} else if (pair_mode == PAIR_OPP_IN || pair_mode == PAIR_OPP_OUT) {
    		bool point_in = false;
    		if (( pair_mode == PAIR_OPP_IN && !((inp.flags & INPUT_FLAG_IS_REVCMPL) && first))
    				|| (pair_mode == PAIR_OPP_OUT && ((inp.flags & INPUT_FLAG_IS_REVCMPL) && first))){
    			point_in = true;
    		}
    		if (point_in && first){
    			ins_size = inp_mp.genome_end - inp.genome_start;
    		} else if(point_in && second){
    			ins_size = inp_mp.genome_start - inp.genome_end;
    		} else if(!point_in && first){
    			ins_size = inp_mp.genome_start - inp.genome_end;
    		} else if(!point_in && second){
    			ins_size = inp_mp.genome_end - inp.genome_start;
    		}
    	}
    }

    free(*output1);
    *output1 = (char *)xmalloc(sizeof(char *)*2000);
    char *extra = *output1 + sprintf(*output1,"%s\t%i\t%s\t%u\t%i\t%s\t%s\t%u\t%i\t%s\t%s\tAS:i:%i",
	    name,
	    ((inp.flags & INPUT_FLAG_IS_REVCMPL) ? 16 : 0) | ((re_mp != NULL) && (inp.flags & INPUT_FLAG_IS_REVCMPL) ? 32 : 0) | ((paired) ? 1 : 0) | (first ? 64 : 0) | (second ? 128 : 0),
	    inp.genome,
	    inp.genome_start + 1,
	    255,
	    cigar,
	    ((re_mp == NULL)?"*":(strcmp(inp.genome,inp_mp.genome)== 0) ? "=": inp_mp.genome),
	    ((re_mp == NULL) ? 0:(inp_mp.genome_start + 1)),
	    ins_size,
	    read,
	    "*",
	    inp.score);
    if (shrimp_mode == COLOUR_SPACE){
    	sprintf(extra,"\tCS:Z:%s",readtostr(read_bitstring,re->read_len,true,first_bp));
    	free(read_ls_bitstring);
    }
    if(re_mp != NULL){
    	free(name);
    }
    free(read);
    free(cigar);
    format_free(fsp);
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
  uint i;

  assert(re != NULL && h != NULL);
  assert(pair_mode == PAIR_NONE);

  /* compute full alignment scores */
  for (i = 0; i < h->load; i++) {
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
       i < h->load && i < num_outputs
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

      hit_output(re, rh, NULL, NULL, &output1, &output2,false);

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

#pragma omp atomic
      total_single_matches++;
    } else {
#pragma omp atomic
      total_dup_single_matches++;
    }
  }

}


/*
 * Do a final pass for given pair of reads.
 * Highest scoring matches are in scores heap.
 */
static void
readpair_pass2(struct read_entry * re1, struct read_entry * re2, struct heap_paired * h) {
  uint i, j;

  assert(re1 != NULL && re2 != NULL && h != NULL);

  /* compute full alignment scores */
  for (i = 0; i < h->load; i++) {
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
       i < h->load && i < num_outputs
	 && ( (IS_ABSOLUTE(sw_full_threshold)
	       && (int)h->array[i].key >= (int)abs_or_pct(sw_full_threshold,
							  h->array[i].rest.hit[0]->score_max + h->array[i].rest.hit[1]->score_max))
	      || (~IS_ABSOLUTE(sw_full_threshold)
		  && (int)h->array[0].key >= (int)sw_full_threshold) );
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

      //re1->paired_final_matches++;

      hit_output(re1, rh1, re2, rh2, &output1, &output2,true);
      hit_output(re2, rh2, re1, rh1, &output3, &output4,true);

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

      free(output1);
      free(output2);
      free(output3);
      free(output4);

#pragma omp atomic
      total_paired_matches++;
    } else {
#pragma omp atomic
      total_dup_paired_matches++;
    }
  }
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
  if (re->has_Ns) {
    free(re->mapidx_pos[0]);
    free(re->mapidx_pos[1]);
  }
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
    fprintf(stderr, "max_n_kmers:%u, min_kmer_pos:%u, has_Ns:%c\n",
	    re->max_n_kmers, re->min_kmer_pos, re->has_Ns ? 'Y' : 'N');
    fprintf(stderr, "collapsed kmers from read:\n");
    for (sn = 0; sn < n_seeds; sn++) {
      fprintf(stderr, "sn:%u\n", sn);
      for (i = 0; re->min_kmer_pos + i + seed[sn].span <= re->read_len; i++) {
	fprintf(stderr, "\tpos:%u kmer:", re->min_kmer_pos + i);
	if (re->has_Ns && !re->mapidx_pos[0][sn*re->max_n_kmers + i]) {
	  fprintf(stderr, "X\n");
	} else {
	  for (j = 0; j < seed[sn].weight; j++) {
	    fprintf(stderr, "%c%s",
		    base_translate((re->mapidx[0][sn*re->max_n_kmers + i] >> 2*(seed[sn].weight - 1 - j)) & 0x3,
				   shrimp_mode == MODE_COLOUR_SPACE),
		    j < seed[sn].weight - 1? "," : "\n");
	  }
	}
      }
    }
    fprintf(stderr, "collapsed kmers from read_rc:\n");
    for (sn = 0; sn < n_seeds; sn++) {
      fprintf(stderr, "sn:%u\n", sn);
      for (i = 0; re->min_kmer_pos + i + seed[sn].span <= re->read_len; i++) {
	fprintf(stderr, "\tpos:%u kmer:", re->min_kmer_pos + i);
	if (re->has_Ns && !re->mapidx_pos[1][sn*re->max_n_kmers + i]) {
	  fprintf(stderr, "X\n");
	} else {
	  for (j = 0; j < seed[sn].weight; j++) {
	    fprintf(stderr, "%c%s",
		    base_translate((re->mapidx[1][sn*re->max_n_kmers + i] >> 2*(seed[sn].weight - 1 - j)) & 0x3,
				   shrimp_mode == MODE_COLOUR_SPACE),
		    j < seed[sn].weight - 1? "," : "\n");
	  }
	}
      }
    }
  }
#endif

  read_get_anchor_list(re, true);

  read_get_hit_list(re, (num_matches >= 2? 2 : 1));
  read_pass1(re, false);

  // initialize heap of best hits for this read
  //re->scores = (struct re_score *)xcalloc((num_outputs + 1) * sizeof(struct re_score));
  //re->scores[0].heap_elems = 0;
  //re->scores[0].heap_capacity = num_outputs;
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
 * Launch the threads that will scan the reads
 *
static bool
launch_scan_threads_OLD(const char *file){
	fasta_t fasta;
	int space;
	int buffer_size = chunk_size*num_threads*64;
	read_entry *res;
	char *seq;

	bool is_rna;

	res = (read_entry *)xmalloc(sizeof(read_entry) * buffer_size);

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
	int i;
#pragma omp parallel shared(res,i,more) num_threads(num_threads)
	{
		while (more){
#pragma omp barrier
#pragma omp single
			{
				for(i = 0; i < buffer_size; i++){
					if(!fasta_get_next(fasta, &res[i].name, &seq, &is_rna)){
						more = false;
						break;
					}
					res[i].read[0] = fasta_sequence_to_bitfield(fasta, seq);
					res[i].read_len = strlen(seq);
					res[i].max_n_kmers = res[i].read_len - min_seed_span + 1;
					res[i].min_kmer_pos = 0;
					if (shrimp_mode == MODE_COLOUR_SPACE){
						res[i].read_len--;
						res[i].max_n_kmers -= 2; // 1st color always discarded from kmers
						res[i].min_kmer_pos = 1;
						res[i].initbp[0] = fasta_get_initial_base(fasta,seq);
						res[i].initbp[1] = res[i].initbp[0];
						res[i].read[1] = reverse_complement_read_cs(res[i].read[0], res[i].initbp[0], res[i].initbp[1],
								res[i].read_len, is_rna);
					} else {
						res[i].read[1] = reverse_complement_read_ls(res[i].read[0],
								res[i].read_len,is_rna);
					}
					res[i].window_len = (uint16_t)abs_or_pct(window_len,res[i].read_len);
					res[i].has_Ns = !(strcspn(seq, "nNxX.") == strlen(seq)); // from fasta.c
					free(seq);
				}
				nreads += i;
			}
			int j;
#pragma omp for schedule(dynamic,chunk_size)
			for (j = 0; j < i; j++){
				handle_read(&res[j]);
			}
		}

	}
	free(res);
	return true;
}
*/


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
  
  int8_t tmp2 = re->initbp[0];
  re->initbp[0] = re->initbp[1];
  re->initbp[1] = tmp2;

  re->input_strand = 1 - re->input_strand;
}


static uint
get_contig_number_from_name(char const * c)
{
  uint cn;
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
  uint cn, st, g_start, g_end;

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
    uint load, i;
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

	re_buffer[i].has_Ns = !(strcspn(re_buffer[i].seq, "nNxX.") == strlen(re_buffer[i].seq)); // from fasta.c
	if (re_buffer[i].has_Ns) { // ignore this read
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
	  re_buffer[i].read[1] = reverse_complement_read_cs(re_buffer[i].read[0], re_buffer[i].initbp[0], re_buffer[i].initbp[1],
							    re_buffer[i].read_len, re_buffer[i].is_rna);
	} else {
	  re_buffer[i].read[1] = reverse_complement_read_ls(re_buffer[i].read[0], re_buffer[i].read_len, re_buffer[i].is_rna);
	}
	re_buffer[i].window_len = (uint16_t)abs_or_pct(window_len,re_buffer[i].read_len);
	re_buffer[i].input_strand = 0;

	re_buffer[i].mapidx[0] = NULL;
	re_buffer[i].mapidx[1] = NULL;
	re_buffer[i].mapidx_pos[0] = NULL;
	re_buffer[i].mapidx_pos[1] = NULL;
	re_buffer[i].anchors[0] = NULL;
	re_buffer[i].anchors[1] = NULL;
	re_buffer[i].hits[0] = NULL;
	re_buffer[i].hits[1] = NULL;
	re_buffer[i].ranges = NULL;
	re_buffer[i].n_ranges = 0;

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


void print_genomemap_stats() {
  stat_t list_size, list_size_non0;
  uint sn;
  uint64_t capacity, mapidx;
  uint max;

  fprintf(stderr, "Genome Map stats:\n");

  for (sn = 0; sn < n_seeds; sn++) {
    if (Hflag)
      capacity = power(4, HASH_TABLE_POWER);
    else
      capacity = power(4, seed[sn].weight);

    stat_init(&list_size);
    stat_init(&list_size_non0);
    for (max = 0, mapidx = 0; mapidx < capacity; mapidx++) {
      if (genomemap_len[sn][mapidx] > list_cutoff)
	continue;

      stat_add(&list_size, genomemap_len[sn][mapidx]);
      if (genomemap_len[sn][mapidx] > 0)
	stat_add(&list_size_non0, genomemap_len[sn][mapidx]);

      if (genomemap_len[sn][mapidx] > max)
	max = genomemap_len[sn][mapidx];
    }

    fprintf(stderr, "sn:%u weight:%u total_kmers:%llu lists:%llu (non-zero:%llu) list_sz_avg:%.2f (%.2f) list_sz_stddev:%.2f (%.2f) max:%u\n",
	    sn, seed[sn].weight, (long long unsigned int)stat_get_sum(&list_size),
	    (long long unsigned int)capacity, (long long unsigned int)stat_get_count(&list_size_non0),
	    stat_get_mean(&list_size), stat_get_mean(&list_size_non0),
	    stat_get_sample_stddev(&list_size), stat_get_sample_stddev(&list_size_non0), max);

#ifdef DEBUG_DUMP_MAP_DIST
#warning Dumping genomemap length distribution
    {
      uint64_t histogram[100];
      uint bucket_size = ceil_div(max, 100);
      uint i, bucket;

      for (i = 0; i < 100; i++) {
	histogram[i] = 0;
      }

      for (mapidx = 0; mapidx < capacity; mapidx++) {
	if (genomemap_len[sn][mapidx] > list_cutoff)
	  continue;

	bucket = genomemap_len[sn][mapidx] / bucket_size;
	histogram[bucket]++;
      }

      for (i = 0; i < 100; i++) {
	fprintf(stderr, "[%d-%d]: %lld\n", i*bucket_size, (i+1)*bucket_size, histogram[i]);
      }
    }
#endif

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
			capacity = power(4, HASH_TABLE_POWER);
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
			for (load = 0; i < seqlen + contig_offsets[num_contigs-1]; i++) {
				uint base;
				int sn;

				base = EXTRACT(read, i - contig_offsets[num_contigs-1]);
				bitfield_prepend(kmerWindow, max_seed_span, base);

				//skip past any Ns or Xs
				if (base == BASE_N || base == BASE_X)
					load = 0;
				else if (load < max_seed_span)
					load++;;
				for (sn = 0; sn < (int)n_seeds; sn++) {
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
					nkmers++;

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


static void
print_insert_histogram()
{
  uint i;
  for (i = 0; i < 100; i++) {
    fprintf(stderr, "[%d-%d]: %.2f%%\n",
	    min_insert_size + i * insert_histogram_bucket_size,
	    min_insert_size + (i + 1) * insert_histogram_bucket_size - 1,
	    ((double)insert_histogram[i] / (double)total_paired_matches) * 100);
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
	for(i = 0; i < (int)num_threads; i++){
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
	f1_total_cellspersec = (double)f1_total_cells / f1_total_secs;
	f2_total_cellspersec = (double)f2_total_cells / f2_total_secs;

	if (Dflag) {
	  fprintf(stderr, "%sPer-Thread Stats:\n", my_tab);
	  fprintf(stderr, "%s%s" "%11s %9s %9s %25s %25s %9s\n", my_tab, my_tab,
		  "", "Read Load", "Scan", "Vector SW", "Scalar SW", "Wait");
	  fprintf(stderr, "%s%s" "%11s %9s %9s %15s %9s %15s %9s %9s\n", my_tab, my_tab,
		  "", "Time", "Time", "Invocs", "Time", "Invocs", "Time", "Time");
	  fprintf(stderr, "\n");
	  for(i = 0; i < (int)num_threads; i++) {
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
	  }

	fprintf(stderr, "\n");

	fprintf(stderr, "%sMemory usage:\n", my_tab);
	fprintf(stderr, "%s%s%-24s" "%s\n", my_tab, my_tab,
		"Genomemap:",
		comma_integer(count_get_count(&mem_genomemap)));

	if (Mflag) {
	  print_insert_histogram();
	}
}

void usage(char *progname,bool full_usage){
	char *slash;
	uint sn;

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
	fprintf(stderr,
		"    -N    Set the number of threads               (default: %u)\n",
		DEF_NUM_THREADS);

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

	fprintf(stderr, "\n");

	fprintf(stderr,
		"    -r    Window Generation Threshold:            "
		  "(default: %.02f%%)\n", DEF_WINDOW_GEN_THRESHOLD);
	if (shrimp_mode == MODE_COLOUR_SPACE) {
	  fprintf(stderr,
		  "    -v    S-W Vector Hit Threshold                "
		  "(default: %.02f%%)\n", DEF_SW_VECT_THRESHOLD);
	}
	fprintf(stderr,
		"    -h    S-W Full Hit Threshold                  "
		"(default: %.02f%%)\n", DEF_SW_FULL_THRESHOLD);

	if (full_usage) {
	  fprintf(stderr, "\n");
	  if (shrimp_mode == MODE_COLOUR_SPACE || shrimp_mode == MODE_LETTER_SPACE) {
	    fprintf(stderr,
		    "    -A    Anchor width limiting full SW           (default: %d; disable: -1)\n",
		    DEF_ANCHOR_WIDTH);
	    fprintf(stderr,
		    "    -K    Set the thread chunk size               (default: %u)\n",
		    DEF_CHUNK_SIZE);
	    fprintf(stderr,
		    "    -Y    Index list cutoff length                (default: %u)\n",
		    DEF_LIST_CUTOFF);
	  }
	}

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
			"    -E    Output in SAM format                          (default: "
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
		"    -M    Print insert size histogram                   (default: disabled)\n");

	fprintf(stderr,
			"    -D    Toggle thread statistics                  (default: %s)\n",
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

  fprintf(stderr, "%s%-40s%s\n", my_tab, "Pair mode:", pair_mode_string[pair_mode]);
  if (pair_mode != PAIR_NONE) {
    fprintf(stderr, "%s%-40smin:%d max:%d\n", my_tab, "Insert sizes:", min_insert_size, max_insert_size);
    if (Mflag) {
      fprintf(stderr, "%s%-40s%d\n", my_tab, "Bucket size:", insert_histogram_bucket_size);
    }
  }

  fprintf(stderr, "\n");

  fprintf(stderr, "%s%-40s%s\n", my_tab, "Gapless mode:", gapless_sw? "yes" : "no");
  fprintf(stderr, "%s%-40s%u\n", my_tab, "Number of threads:", num_threads);
  fprintf(stderr, "%s%-40s%u\n", my_tab, "Thread chuck size:", chunk_size);
  fprintf(stderr, "%s%-40s%s\n", my_tab, "Hash Filter Calls:", hash_filter_calls? "yes" : "no");
  fprintf(stderr, "%s%-40s%d%s\n", my_tab, "Anchor Width:", anchor_width,
	  anchor_width == -1? " (disabled)" : "");
  fprintf(stderr, "%s%-40s%u\n", my_tab, "Index List Cutoff Length:", list_cutoff);

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
		optstr = "?s:n:w:o:m:i:g:q:e:f:x:h:v:r:N:K:BCFHPRTUA:MW:S:L:DZGp:I:EY:";
		break;
	case MODE_LETTER_SPACE:
		optstr = "?s:n:w:o:m:i:g:q:e:f:h:r:X:N:K:BCFHPRTUA:MW:S:L:DZGp:I:EY:";
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
		case 'r':
		  window_gen_threshold = atof(optarg);
		  if (window_gen_threshold < 0.0) {
		    fprintf(stderr, "error: invalid window generation threshold [%s]\n", optarg);
		    exit(1);
		  }
		  if (strcspn(optarg, "%.") == strlen(optarg))
		    window_gen_threshold = -window_gen_threshold;	//absol.
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
		case 'Z':
		  hash_filter_calls = !hash_filter_calls;
		  break;
		case 'G':
		  gapless_sw = true;
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
		  min_insert_size = (uint)atoi(c);
		  c = strtok(NULL, ",");
		  if (c == NULL) {
		    fprintf(stderr, "error: format for insert sizes is \"-I 200,1000\"\n");
		    exit(1);
		  }
		  max_insert_size = (uint)atoi(c);
		  if (min_insert_size > max_insert_size) {
		    fprintf(stderr, "error: invalid insert sizes (min:%d,max:%d)\n",
			    min_insert_size, max_insert_size);
		    exit(1);
		  }
		  break;
		case 'E':
			Eflag = true;
			break;
		case 'Y':
		  list_cutoff = atoi(optarg);
		  if (list_cutoff == 0) {
		    fprintf(stderr, "error: invalid list cutoff (%s)\n", optarg);
		    exit(1);
		  }
		  break;
		default:
			usage(progname, false);
		}
	}

	argc -= optind;
	argv += optind;

	insert_histogram_bucket_size = ceil_div(max_insert_size - min_insert_size + 1, 100);

	if(load_file != NULL && n_seeds != 0){
		fprintf(stderr,"error: cannot specify seeds when loading genome map\n");
		usage(progname,false);
	}

	if (n_seeds == 0 && load_file == NULL)
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

	if (Eflag && Pflag){
		fprintf(stderr,"-E and -P are incompatable\n");
		exit(1);
	}
	if (Eflag && Rflag){
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
		fprintf(stderr, "error: invalid s-w vector threshold\n");
		exit(1);
	}

	if (!IS_ABSOLUTE(window_gen_threshold)
	    && window_gen_threshold > 100.0) {
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
	if(load_file == NULL){
		print_settings();
	}

	uint64_t before;
	before = gettimeinusecs();
	if (load_file != NULL){
		if (strchr(load_file, ',') == NULL) {
			fprintf(stderr,"Must specify genome save file and map save file\n");
			exit (1);
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
		print_settings();
	} else {
		if (!load_genome(genome_files,ngenome_files)){
			exit(1);
		}
	}
	map_usecs += (gettimeinusecs() - before);

	print_genomemap_stats();

	if (save_file != NULL){
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

		uint s;
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
#endif
