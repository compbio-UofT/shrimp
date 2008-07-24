/*	$Id$	*/

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

#include <xmmintrin.h>	// for _mm_prefetch

#include <sys/types.h>
#include <sys/stat.h>
#include <sys/time.h>

#include "../common/fasta.h"
#include "../common/dag_glue.h"
#include "../common/sw-full-common.h"
#include "../common/sw-full-cs.h"
#include "../common/sw-full-ls.h"
#include "../common/sw-vector.h"
#include "../common/output.h"
#include "../common/util.h"
#include "../common/version.h"

#include "rmapper.h"

/* External parameters */
static char    *spaced_seed		= "0";
static double	window_len		= DEF_WINDOW_LEN;
static u_int	num_matches		= DEF_NUM_MATCHES;
static u_int	hit_taboo_len		= DEF_HIT_TABOO_LEN;
static u_int	seed_taboo_len		= DEF_SEED_TABOO_LEN;
static u_int	num_outputs		= DEF_NUM_OUTPUTS;
static int	kmer_stddev_limit  	= DEF_KMER_STDDEV_LIMIT;

static u_int	dag_epsilon		= DEF_DAG_EPSILON;
static int	dag_read_match		= DEF_DAG_READ_MATCH;
static int	dag_read_mismatch	= DEF_DAG_READ_MISMATCH;
static int	dag_read_gap		= DEF_DAG_READ_GAP;

static int	dag_ref_match		= DEF_DAG_REF_MATCH;
static int 	dag_ref_mismatch	= DEF_DAG_REF_MISMATCH;
static int	dag_ref_half_match	= DEF_DAG_REF_HALF_MATCH;
static int	dag_ref_neither_match	= DEF_DAG_REF_NEITHER_MATCH;
static int	dag_ref_match_deletion	= DEF_DAG_REF_MATCH_DELETION;
static int	dag_ref_mismatch_deletion=DEF_DAG_REF_MISMATCH_DELETION;
static int	dag_ref_error_insertion	= DEF_DAG_REF_ERROR_INSERTION;
static double	dag_ref_weighted_thresh = DEF_DAG_REF_WEIGHTED_THRESHOLD;

static int	match_value		= DEF_MATCH_VALUE;
static int	mismatch_value		= DEF_MISMATCH_VALUE;
static int	a_gap_open		= DEF_A_GAP_OPEN;
static int	a_gap_extend		= DEF_A_GAP_EXTEND;
static int	b_gap_open		= DEF_B_GAP_OPEN;
static int	b_gap_extend		= DEF_B_GAP_EXTEND;
static int	xover_penalty		= DEF_XOVER_PENALTY;
static double	sw_vect_threshold	= DEF_SW_VECT_THRESHOLD;
static double	sw_full_threshold	= DEF_SW_FULL_THRESHOLD;

/*
 * If window_len, sw_vect_threshold, sw_full_threshold are absolute values,
 * we'll set them negative to distinguish.
 */
#define IS_ABSOLUTE(x)	((x) < 0)

/* Genomic sequence, stored in 32-bit words, first is in the LSB */
static uint32_t *genome;			/* genome -- always in letter*/
static uint32_t *genome_cs;			/* genome -- colourspace */
static uint32_t	 genome_len;
static char     *contig_name;			/* name of contig in genome */
static uint32_t  ncontigs;			/* number of contigs read */

/* Reads */
#define OFFSET_TO_READ(_o)	\
    (struct read_elem *)((char *)reads + (read_size * (_o)))
static struct read_elem      *reads;		/* list of reads */
static uint32_t		      reads_allocated;	/* slots allocated for reads */
static uint32_t		      read_size;	/* size of each read struct */

static u_int	nreads;				/* total reads loaded */
static u_int	nkmers;				/* total kmers of reads loaded*/
static u_int	longest_read_len;		/* longest read we've seen */
static u_int	nduphits;			/* number of duplicate hits */

/* Kmer to read index */
static struct readmap_entry **readmap;		/* kmer to read map */
static uint32_t		     *readmap_len;	/* entries in each readmap[n] */

/* Flags */
static int Bflag = false;			/* print a progress bar */
static int Cflag = false;			/* do complement only */
static int Fflag = false;			/* do positive (forward) only */
static int Pflag = false;			/* pretty print results */
static int Rflag = false;			/* add read sequence to output*/

/* Scan stats */
static uint64_t shortest_scanned_kmer_list;
static uint64_t longest_scanned_kmer_list;
static uint64_t kmer_lists_scanned;
static uint64_t kmer_list_entries_scanned;

/* Misc stats */
static uint64_t scan_ticks;
static uint64_t revcmpl_ticks;

#define PROGRESS_BAR(_a, _b, _c, _d)	\
    if (Bflag) progress_bar((_a), (_b), (_c), (_d))

#define FETCH_AHEAD 8
#define PF_HINT _MM_HINT_T2
#define USE_PREFETCH

static size_t
power(size_t base, size_t exp)
{
	size_t i, ret;

	if (exp == 0)
		return (1);

	ret = base;
	for (i = 1; i < exp; i++)
		ret *= base;

	return (ret);
}

/* percolate down in our min-heap */
static void
reheap(struct re_score *scores, int node)
{
	struct re_score tmp;
	int left, right, max;

	assert(node >= 1 && node <= (int)scores[0].index);

	left  = node * 2;
	right = left + 1;
	max   = node;
	
	if (left <= scores[0].score &&
	    scores[left].score < scores[node].score)
		max = left;

	if (right <= scores[0].score &&
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
percolate_up(struct re_score *scores, int node)
{
	struct re_score tmp;
	int parent;

	assert(node >= 1 && node <= (int)scores[0].index);

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
save_score(struct read_elem *re, int score, int index, int contig_num,
    bool revcmpl)
{
	struct re_score *scores;

	scores = re->ri->scores;

	if (scores[0].score == (int)num_outputs) {
		if (score < scores[1].score)
			return;

		scores[1].score = score;
		scores[1].index = index;
		scores[1].contig_num = contig_num;
		scores[1].revcmpl = revcmpl;
		reheap(scores, 1);
	} else {
		int idx = 1 + scores[0].score++;

		/* We do the array doubling trick for O(1) amortised alloc. */
		if (scores[0].score > (int)scores[0].index) {
			unsigned int alloc_slots;

			assert(scores[0].score > 0);

			if (((unsigned int)scores[0].score * 2) > num_outputs)
				alloc_slots = num_outputs;
			else
				alloc_slots = scores[0].score * 2;

			assert(alloc_slots <= num_outputs);

			scores[0].index = alloc_slots;
			scores = (struct re_score *)xrealloc(scores,
			    sizeof(struct re_score) * (alloc_slots + 1));
			memset(&scores[idx], 0,
			    sizeof(struct re_score) * (alloc_slots + 1 - idx));
			re->ri->scores = scores;
		}

		scores[idx].score = score;
		scores[idx].index = index;
		scores[idx].contig_num = contig_num;
		scores[idx].revcmpl = revcmpl;
		percolate_up(scores, idx);
	}
}

/*
 * Compress the given kmer into an index in 'readmap' according to the seed.
 * While not optimal, this is only about 20% of the spaced seed scan time.
 */
static uint32_t
kmer_to_mapidx(uint32_t *kmer)
{
	static u_int seed_span;

	uint32_t mapidx;
	int i;

	if (seed_span == 0)
		seed_span = strlen(spaced_seed);

	mapidx = 0;
	for (i = seed_span - 1; i >= 0; i--) {
		if (spaced_seed[i] == '1') {
			mapidx = mapidx << 2;
			mapidx |= (kmer[i / 8] & (0x3 << ((i % 8) * 4))) >>
			    ((i % 8) * 4);
		}
	}

	assert(mapidx < power(4, seed_span));

	return (mapidx);
}

/*
 * If 'kmer_stddev_limit' is set to a reasonable value, prune all kmers, whose
 * frequencies are more than 'kmer_stddev_limit' standard deviations greater
 * than the mean.
 */
static void
readmap_prune()
{
	double mean, stddev;
	uint64_t pruned = 0;
	uint32_t i, j;

	if (kmer_stddev_limit < 0)
		return;

	j = power(4, strchrcnt(spaced_seed, '1'));

	mean = 0;
	for (i = 0; i < j; i++)
		mean += readmap_len[i];
	mean /= j;

	stddev = 0;
	for (i = 0; i < j; i++)
		stddev += pow((double)readmap_len[i] - mean, 2);
	stddev = sqrt(stddev / j);

	fprintf(stderr, "- Pruning kmers; mu: %f, sigma: %f, "
	    "sigma limit: +/- %f\n", mean, stddev,
	    kmer_stddev_limit * stddev);
	if (mean < 1.0)
		fprintf(stderr, "WARNING: low mean - are you sure you want to "
		    "prune kmers?\n");

	for (i = 0; i < j; i++) {
		if (readmap_len[i] > mean + (kmer_stddev_limit * stddev)) {
			readmap[i] = NULL;
			pruned++;
		}
	}

	fprintf(stderr, "  - Pruned %" PRIu64 " kmer(s)\n", pruned);
}

/* reset the fields used by scan() for each read */
static void
reset_reads()
{
	uint32_t i, j;

	for (i = 0; i < nreads; i++) {
		struct read_elem *re = OFFSET_TO_READ(i);

		re->last_swhit_idx = UINT32_MAX;
		re->prev_hit = 0;
		re->next_hit = 0;
		for (j = 0; j < num_matches; j++)
			re->hits[j].g_idx = re->hits[j].r_idx = UINT32_MAX;
	}
}

static int
get_sw_threshold(struct read_elem *re, double which, int readnum)
{
	int thresh;

	assert(readnum == 1 || readnum == 2);

	if (IS_ABSOLUTE(which)) {
		thresh = (u_int)-which;
	} else {
		if (readnum == 1) {
			thresh = (u_int)floor(((double)
			    (match_value * re->ri->read1_len) *
			    (which / 100.0)));
		} else {
			assert(readnum == 2);
			assert(shrimp_mode == MODE_HELICOS_SPACE);
			thresh = (u_int)floor(((double)
			    (match_value * re->ri->read2_len) *
			    (which / 100.0)));
		}
	}

	return (thresh);
}


/*
 * Determine whether or not the kmer hits seen within this read
 * were colinear.
 */ 
#if 0
static int
are_hits_colinear(struct read_elem *re)
{
	int i, prev, this;

	this = prev = re->next_hit;
	for (i = 1; i < num_matches; i++) {
		this = (this + 1) % num_matches;

		if (re->hits[prev].r_idx > re->hits[this].r_idx &&
		    re->hits[prev].r_idx != UINT32_MAX && re->hits[this].r_idx != UINT32_MAX) {
			return (1);
		}

		prev = this;
	}

	return (1);
}
#endif

/* scan the genome by kmers, running S-W as needed, and updating scores */
static void
scan(int contig_num, bool revcmpl)
{
	struct read_elem *re;
	struct readmap_entry *rme1, *rme2;
	uint32_t *kmer, *scan_genome;
	uint32_t i, j, k, pf, idx, mapidx, next_mapidx = 0;
	uint32_t base, skip, score1, score2;
	uint32_t goff, glen, prevhit;
	u_int seed_span;
	u_int thresh;

	if (shrimp_mode == MODE_COLOUR_SPACE)
		scan_genome = genome_cs;
	else
		scan_genome = genome;

	seed_span = strlen(spaced_seed);

	PROGRESS_BAR(stderr, 0, 0, 10);

	kmer = (uint32_t *)xmalloc(sizeof(kmer[0]) * BPTO32BW(seed_span));
	memset(kmer, 0, sizeof(kmer[0]) * BPTO32BW(seed_span));

	skip = seed_span - 1;
	for (i = 0; i < genome_len; i++) {
		PROGRESS_BAR(stderr, i, genome_len, 10);

		base = EXTRACT(scan_genome, i);
		if (i == (seed_span - 1)) {
			bitfield_prepend(kmer, seed_span, base);
			mapidx = kmer_to_mapidx(kmer);
		} else
			mapidx = next_mapidx;

		/* compute next kmer early and prefetch the readmap entry */
		if ((i + 1) < genome_len) {
			int tmp = EXTRACT(scan_genome, i + 1);
			bitfield_prepend(kmer, seed_span, tmp);
			next_mapidx = kmer_to_mapidx(kmer);
#ifdef USE_PREFETCH
			_mm_prefetch((const char *)&readmap[next_mapidx], PF_HINT);
#endif
		}

		/*
		 * Ignore kmers that contain an 'N'.
		 */
		if (base > 3)
			skip = seed_span - 1;

		if (skip > 0) {
			skip--;
			continue;
		}

		idx = i - (seed_span - 1);

		rme1 = rme2 = readmap[mapidx];
		if (rme1 == NULL) {
			kmer_lists_scanned++;
			continue;
		}

		/* prefetch first FETCH_AHEAD struct read_elem's */
#ifdef USE_PREFETCH
		for (pf = 0; rme2->offset != UINT32_MAX && pf < FETCH_AHEAD; pf++)
			_mm_prefetch((const char *)OFFSET_TO_READ((rme2++)->offset), PF_HINT);
#endif

		for (re = OFFSET_TO_READ(rme1->offset), j = k = 0;
		     rme1->offset != UINT32_MAX;
		     re = OFFSET_TO_READ((++rme1)->offset), k++) {
			/* continue prefetching the struct read_elem's */
#ifdef USE_PREFETCH
			if (rme2->offset != UINT32_MAX)
				_mm_prefetch((const char *)OFFSET_TO_READ((rme2++)->offset), PF_HINT);
#endif

			prevhit = re->hits[re->prev_hit].g_idx;
			if ((idx - prevhit) <= hit_taboo_len && prevhit != UINT32_MAX)
				continue;

			re->hits[re->next_hit].g_idx = idx;
			re->hits[re->next_hit].r_idx = rme1->r_idx;
			re->prev_hit = re->next_hit;
			re->next_hit = (uint8_t)((re->next_hit + 1) % num_matches);

			if ((idx - re->hits[re->next_hit].g_idx) < re->window_len &&
			    (re->last_swhit_idx == UINT32_MAX ||
			     (i - re->last_swhit_idx) >= re->window_len) &&
			    re->hits[re->next_hit].g_idx != UINT32_MAX) { // && are_hits_colinear(re)) {
				bool meets_thresh = false;

				if (i < re->window_len)
					goff = 0;
				else
					goff = i - re->window_len;

				if (goff + (re->window_len * 2) >= genome_len)
					glen = genome_len - goff;
				else
					glen = re->window_len * 2;

				if (shrimp_mode == MODE_COLOUR_SPACE) {
					score1 = sw_vector(scan_genome, goff, glen,
					    re->ri->read1, re->ri->read1_len,
					    genome, re->ri->initbp);
					thresh = get_sw_threshold(re, sw_vect_threshold, 1);
					if (score1 >= thresh)
						meets_thresh = true;
				} else if (shrimp_mode == MODE_LETTER_SPACE) {
					score1 = sw_vector(scan_genome, goff, glen,
					    re->ri->read1, re->ri->read1_len, NULL, -1);
					thresh = get_sw_threshold(re, sw_vect_threshold, 1);
					if (score1 >= thresh)
						meets_thresh = true;
				} else if (shrimp_mode == MODE_HELICOS_SPACE) {
					score1 = sw_vector(scan_genome, goff, glen,
					    re->ri->read1, re->ri->read1_len, NULL, -1);
					thresh = get_sw_threshold(re, sw_vect_threshold, 1);
					if (score1 >= thresh) {
						score2 = sw_vector(scan_genome, goff, glen,
						    re->ri->read2, re->ri->read2_len, NULL, -1);
						thresh = get_sw_threshold(re, sw_vect_threshold, 2);
						if (score2 >= thresh) {
							score1 = score1 * re->ri->read1_len;
							score2 = score2 * re->ri->read2_len;
							score1 = (score1 + score2) /
							    (re->ri->read1_len + re->ri->read2_len);
							meets_thresh = true;
						}
					}
				}

				if (meets_thresh) {
					save_score(re, score1, i, contig_num,
					    revcmpl);
					re->last_swhit_idx = i;
					re->ri->swhits++;
				}
			}
			j++;
		}

		assert(k == readmap_len[mapidx]);

		kmer_list_entries_scanned += j;
		shortest_scanned_kmer_list = MIN(shortest_scanned_kmer_list, j);
		longest_scanned_kmer_list = MAX(longest_scanned_kmer_list, j);
		kmer_lists_scanned++;
	}

	if (Bflag) {
		PROGRESS_BAR(stderr, genome_len, genome_len, 10);
		putc('\n', stderr);
	}

	free(kmer);
}

static struct read_elem *
readalloc()
{
	struct read_elem *re;
	size_t bytes, prevbytes, i;

	assert(nreads <= reads_allocated);

	if (nreads == reads_allocated) {
		read_size = sizeof(*re) +
		    (sizeof(re->hits[0]) * (num_matches));

		prevbytes = read_size * reads_allocated;
		reads_allocated += 1000;
		bytes = read_size * reads_allocated;
		reads = (struct read_elem *)xrealloc(reads, bytes);
		memset((char *)reads + prevbytes, 0, bytes - prevbytes);
	}

	re = OFFSET_TO_READ(nreads);

	re->ri = (struct read_int *)xmalloc(sizeof(*re->ri));
	memset(re->ri, 0, sizeof(*re->ri));
	re->ri->offset = nreads;

	for (i = 0; i < num_matches; i++)
		re->hits[i].g_idx = re->hits[i].r_idx = UINT32_MAX;

	bytes = sizeof(struct re_score) * 1;
	re->ri->scores = (struct re_score *)xmalloc(bytes);
	memset(re->ri->scores, 0, bytes);

	nreads++;

	return (re);
}

static bool
load_reads_lscs(const char *file)
{
	struct read_elem *re;
	u_int seed_span;
	char *name, *seq;
	fasta_t fasta;
	size_t seqlen, skip = 0;
	ssize_t bases = 0;
	int space;

	seed_span = strlen(spaced_seed);

	if (shrimp_mode == MODE_LETTER_SPACE)
		space = LETTER_SPACE;
	else
		space = COLOUR_SPACE;
	fasta = fasta_open(file, space);
	if (fasta == NULL) {
		fprintf(stderr, "error: failed to open reads file [%s]\n",
		    file);
		return (false);
	}

	fprintf(stderr, "- Loading reads...");

	while (true) {
		u_int i;
		uint32_t *kmer;

		if (Bflag && (nreads % 17) == 0) 
			fprintf(stderr, "\r- Loading reads... %u", nreads);

		if (!fasta_get_next(fasta, &name, &seq))
			break;

		if (strchr(name, '\t') != NULL || strchr(seq, '\t') != NULL) {
			fprintf(stderr, "error: tabs are not permitted in fasta names "
			    "or sequences. Tag: [%s].\n", name);
			exit(1);
		}

		seqlen = strlen(seq);
		if (seqlen == 0) {
			fprintf(stderr, "error: read [%s] had no sequence!\n",
			    name);
			exit(1);
		}
		if (shrimp_mode == MODE_COLOUR_SPACE) {
			/* the sequence begins with the initial letter base */
			if (seqlen < 1) {
				fprintf(stderr, "error: read [%s] had sequence "
				    "with no colours!\n", name);
				exit(1);
			}
			seqlen--;
		}
		bases += seqlen;

		re = readalloc();
		re->last_swhit_idx = UINT32_MAX;
		re->ri->name = name;
		re->ri->read1 = fasta_sequence_to_bitfield(fasta, seq);
		if (re->ri->read1 == NULL) {
			fprintf(stderr, "error: invalid sequence; tag: [%s]\n", name); 
			exit(1);
		}
		re->ri->read1_len = seqlen;
		longest_read_len = MAX(re->ri->read1_len, longest_read_len);

		if (shrimp_mode == MODE_COLOUR_SPACE)
			re->ri->initbp = fasta_get_initial_base(fasta, seq);

		kmer = (uint32_t *)xmalloc(BPTO32BW(seed_span) * sizeof(uint32_t));
		memset(kmer, 0, BPTO32BW(seed_span) * sizeof(uint32_t));

		for (i = 0; i < seed_span - 1; i++)
			bitfield_prepend(kmer, seed_span, EXTRACT(re->ri->read1, i));

		for (; i < seqlen; i++) {
			int base;
			uint32_t mapidx;

			base = EXTRACT(re->ri->read1, i);
			bitfield_prepend(kmer, seed_span, base);

			/*
			 * Ignore kmers that contain an 'N'.
			 */
			if (base > 3)
				skip = MAX(seed_span - 1, skip);

			if (skip > 0) {
				skip--;
				continue;
			}

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
			if (shrimp_mode == MODE_COLOUR_SPACE && i == seed_span - 1)
				continue;

			mapidx = kmer_to_mapidx(kmer);

			/* permit only one kmer reference per read */
			if (readmap[mapidx] != NULL &&
			    readmap[mapidx][readmap_len[mapidx] - 1].offset == re->ri->offset)
				continue;
			
			readmap_len[mapidx]++;
			readmap[mapidx] = (struct readmap_entry *)xrealloc(
			    readmap[mapidx], (readmap_len[mapidx] + 1) * sizeof(*readmap[0]));
			readmap[mapidx][readmap_len[mapidx] - 1].offset	= re->ri->offset;
			readmap[mapidx][readmap_len[mapidx] - 1].r_idx	= UINT32_MAX;
			readmap[mapidx][readmap_len[mapidx]].offset	= UINT32_MAX;
			readmap[mapidx][readmap_len[mapidx]].r_idx	= UINT32_MAX;
			nkmers++;

			/* Taboo the seed generation as well, if requested. */
			if (seed_taboo_len)
				skip = seed_taboo_len;
		}

		free(seq);
		free(kmer);
	}
	if (Bflag)
		fprintf(stderr, "\r- Loading reads... %u\n", nreads);
	else
		fprintf(stderr, "\n");

	fprintf(stderr, "- Loaded %s %s in %s reads (%s kmers)\n",
	    comma_integer(bases), (shrimp_mode == MODE_COLOUR_SPACE) ?
	    "colours" : "letters", comma_integer(nreads),comma_integer(nkmers));

	fasta_close(fasta);

	return (true);
}

static bool
load_reads_dag(const char *file)
{
	struct read_elem *re;
	int seed_span;
	char *name1, *name2, *seq1, *seq2;
	fasta_t fasta;
	size_t seq1len, seq2len;
	ssize_t bases = 0;

	seed_span = strlen(spaced_seed);

	fasta = fasta_open(file, LETTER_SPACE);
	if (fasta == NULL) {
		fprintf(stderr, "error: failed to open reads file [%s]\n",
		    file);
		return (false);
	}

	fprintf(stderr, "- Loading read pairs...");

	while (true) {
		int i;
		bool b1, b2;
		char **kmers;

		if (Bflag && (nreads % 17) == 0) 
			fprintf(stderr, "\r- Loading read pairs... %u", nreads);

		b1 = fasta_get_next(fasta, &name1, &seq1);
		b2 = fasta_get_next(fasta, &name2, &seq2);

		if (b1 == false)
			break;

		if (b2 == false) {
			fprintf(stderr, "error: helicos fasta file [%s] does "
			    "not appear to have a pair for read [%s]\n",
			    file, name1);
			return (false);
		}

		if (strchr(name1, '\t') != NULL || strchr(seq1, '\t') != NULL) {
			fprintf(stderr, "error: tabs are not permitted in fasta names "
			    "or sequences. Tag: [%s].\n", name1);
			exit(1);
		}
		if (strchr(name2, '\t') != NULL || strchr(seq2, '\t') != NULL) {
			fprintf(stderr, "error: tabs are not permitted in fasta names "
			    "or sequences. Tag: [%s].\n", name2);
			exit(1);
		}

		seq1len = strlen(seq1);
		seq2len = strlen(seq2);
		bases += seq1len + seq2len;

		if (seq1len == 0 || seq2len == 0) {
			fprintf(stderr, "error: read [%s] had no sequence!\n",
			    name1);
			exit(1);
		}

		re = readalloc();
		re->last_swhit_idx = UINT32_MAX;
		re->ri->name = name1;
		re->ri->read1 = fasta_sequence_to_bitfield(fasta, seq1);
		if (re->ri->read1 == NULL) {
			fprintf(stderr, "error: invalid sequence; tag: [%s]\n", name1);
			exit(1);
		}
		re->ri->read2 = fasta_sequence_to_bitfield(fasta, seq2);
		if (re->ri->read2 == NULL) {
			fprintf(stderr, "error: invalid sequence; tag: [%s]\n", name2);
			exit(1);
		}
		re->ri->read1_len = seq1len;
		re->ri->read2_len = seq2len;

		longest_read_len = MAX(re->ri->read1_len, longest_read_len);
		longest_read_len = MAX(re->ri->read2_len, longest_read_len);

		re->ri->dag_cookie = dag_build_kmer_graph(seq1, seq2, dag_epsilon);
		kmers = dag_get_kmers(re->ri->dag_cookie, seed_span);
		assert(kmers != NULL);

		for (i = 0; kmers[i] != NULL; i++) {
			uint32_t *kmer;
			uint32_t mapidx;

			/* scan() builds kmers in reverse, so we must as well... */
			strrev(kmers[i]);

			kmer = fasta_sequence_to_bitfield(fasta, kmers[i]);
			if (kmer == NULL) {
				fprintf(stderr, "error: invalid hs kmer: [%s]\n", kmers[i]);
				exit(1);
			}
			mapidx = kmer_to_mapidx(kmer);
			free(kmer);

			/* permit only one kmer reference per read */
			if (readmap[mapidx] != NULL &&
			    readmap[mapidx][readmap_len[mapidx] - 1].offset == re->ri->offset)
				continue;

			readmap_len[mapidx]++;
			readmap[mapidx] = (struct readmap_entry *)xrealloc(
			    readmap[mapidx], (readmap_len[mapidx] + 1) * sizeof(*readmap[0]));
			readmap[mapidx][readmap_len[mapidx] - 1].offset	= re->ri->offset;
			readmap[mapidx][readmap_len[mapidx] - 1].r_idx	= UINT32_MAX;
			readmap[mapidx][readmap_len[mapidx]].offset	= UINT32_MAX;
			readmap[mapidx][readmap_len[mapidx]].r_idx	= UINT32_MAX;
			nkmers++;
		}
		if (i == 0)
			fprintf(stderr, "WARNING: Got %d kmers for [%s] and [%s]\n", i, name1, name2);

		dag_free_kmers(kmers);
		free(name2);
		free(seq1);
		free(seq2);
	}
	if (Bflag)
		fprintf(stderr, "\r- Loading read pairs... %u\n", nreads);
	else
		fprintf(stderr, "\n");

	fprintf(stderr, "- Loaded %s letters in %u read pairs (%u kmers)\n",
	    comma_integer(bases), nreads, nkmers);

	fasta_close(fasta);

	return (true);
}

static bool
load_reads(const char *file)
{
	int seed_weight;
	ssize_t ret, bytes;

	seed_weight = strchrcnt(spaced_seed, '1');

	/* allocate our read map */
	bytes = sizeof(readmap[0]) * power(4, seed_weight);
	readmap = (struct readmap_entry **)xmalloc(bytes);
	memset(readmap, 0, bytes);

	bytes = sizeof(readmap_len[0]) * power(4, seed_weight);
	readmap_len = (uint32_t *)xmalloc(bytes);
	memset(readmap_len, 0, bytes);

	if (shrimp_mode == MODE_HELICOS_SPACE)
		ret = load_reads_dag(file);
	else
		ret = load_reads_lscs(file);

	return (ret);
}

/*
 * XXX - this needs to be hooked into our regular output format.
 */
static void
generate_output_dag(struct re_score *rs_array, size_t rs_len, bool revcmpl)
{
	uint32_t goff, glen;
	size_t i;

	/* Obtain final results for each element of our array. */
	for (i = 0; i < rs_len; i++) {
		struct re_score *rs = &rs_array[i];
		struct read_elem *re = rs->parent;
		double score;
		bool meets_thresh = false;

		re = rs->parent;

		if (rs->index < re->window_len)
			goff = 0;
		else
			goff = rs->index - re->window_len;

		assert(goff < genome_len);

		if (goff + (re->window_len * 2) >= genome_len)
			glen = genome_len - goff;
		else
			glen = re->window_len * 2;

		u_int j;
		char *refseq = (char *)xmalloc(glen + 1);

		for (j = 0; j < glen; j++) {
			refseq[j] = base_translate(EXTRACT(genome, goff + j), false);
		}
		refseq[j] = '\0';

		struct dag_alignment *dap = dag_build_alignment(refseq, re->ri->dag_cookie);

		double avglen;

		avglen = (double)(re->ri->read1_len + re->ri->read2_len) / 2;
		if (((double)dap->score / avglen) >= dag_ref_weighted_thresh) {
			meets_thresh = true;
			score = ((double)dap->score / avglen);
		}

		if (meets_thresh) {
			uint32_t genome_start = goff + dap->start_index + 1;
			uint32_t genome_end = goff + dap->end_index;

			if (revcmpl) {
				uint32_t tmp = genome_start;
				genome_start = genome_len - genome_end - 1;
				genome_end = genome_len - tmp - 1;
			}

			printf("\n--------------------------------------------------------------------------------\n");
			printf("ALIGNMENT:\n");
			printf("  contig: [%s]\n", contig_name);
			printf("  read:   [%s]\n", re->ri->name);
			printf("  strand: %c, score %.3f, genome_start: %u, genome_end: %u\n",
			    revcmpl ? '-' : '+', score, genome_start, genome_end);
			printf("\n");
			printf("  READ 1:     %s\n", dap->read1);
			printf("  READ 2:     %s\n", dap->read2);
			printf("  SEQUENCE:   %s\n", dap->sequence);
			printf("--------------------------------------------------------------------------------\n");
			printf("\n");
			re->ri->final_matches++;
		}

		free(refseq);
		dag_free_alignment(dap);
	}
}

static int
score_cmp(const void *arg1, const void *arg2)
{
	const struct re_score *one = (const struct re_score *)arg1;
	const struct re_score *two = (const struct re_score *)arg2;

	assert(one->revcmpl == two->revcmpl);

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

static void
generate_output_lscs(struct re_score *rs_array, size_t rs_len, bool revcmpl)
{
	static bool firstcall = true;

	struct sw_full_results last_sfr;
	uint32_t goff, glen;
	int thresh;
	u_int i;

	assert(shrimp_mode == MODE_LETTER_SPACE ||
	       shrimp_mode == MODE_COLOUR_SPACE);

	if (rs_len == 0)
		return;

	/* shut up, gcc */
	thresh = 0;

	/* Obtain final results for each element of our array. */
	for (i = 0; i < rs_len; i++) {
		struct re_score *rs = &rs_array[i];
		struct read_elem *re = rs->parent;

		if (rs->index < re->window_len)
			goff = 0;
		else
			goff = rs->index - re->window_len;

		assert(goff < genome_len);

		if (goff + (re->window_len * 2) >= genome_len)
			glen = genome_len - goff;
		else
			glen = re->window_len * 2;

		if (IS_ABSOLUTE(sw_full_threshold)) {
			thresh = (u_int)-sw_full_threshold;
		} else {
			thresh = (u_int)floor(((double)
			    (match_value * re->ri->read1_len) *
			    (sw_full_threshold / 100.0)));
		}

		/*
		 * Hacky Optimisation:
		 *   In letter-space, our full alignment will be the same as
		 *   what our vector aligner detected. For this reason, we
		 *   can call the vector aligner a few more times to home in
		 *   on the alignment we want and reduce our genome length
		 *   in the much slower, full aligner, providing a net win.
		 */
		if (shrimp_mode == MODE_LETTER_SPACE) {
			uint32_t tryoff, trylen;
			int score;

			trylen = glen / 2;
			tryoff = goff + trylen / 2;
			assert(tryoff < (goff + glen));
			score = sw_vector(genome, tryoff, trylen,
			    re->ri->read1, re->ri->read1_len, NULL, -1);
			if (score == rs->score) {
				goff = tryoff;
				glen = trylen;
			}
		}

		rs->sfrp = (struct sw_full_results *)xmalloc(sizeof(*rs->sfrp));
		if (shrimp_mode == MODE_COLOUR_SPACE) {
			sw_full_cs(genome, goff, glen, re->ri->read1,
			    re->ri->read1_len, re->ri->initbp, thresh, rs->sfrp);
		} else {
			sw_full_ls(genome, goff, glen, re->ri->read1,
			    re->ri->read1_len, thresh, rs->score, rs->sfrp);
			assert(rs->sfrp->score == rs->score);
		}
		rs->score = rs->sfrp->score;
	}

	/* Sort the scores. */
	qsort(rs_array, rs_len, sizeof(rs_array[0]), score_cmp);

	/* Output sorted list, removing any duplicates. */
	memset(&last_sfr, 0, sizeof(last_sfr));
	for (i = 0; i < rs_len; i++) {
		struct re_score *rs = &rs_array[i];
		struct read_elem *re = rs->parent;
		bool dup;

		dup = sw_full_results_equal(&last_sfr, rs->sfrp);
		if (dup)
			nduphits++;

		if (rs->score >= thresh && dup == false) {
			char *output;

			re->ri->final_matches++;

			if (firstcall) {
				output = output_format_line(Rflag);
				puts(output);
				free(output);
			}
			firstcall = false;

			output = output_normal(re->ri->name, contig_name, rs->sfrp,
			    genome_len, (shrimp_mode == MODE_COLOUR_SPACE), re->ri->read1,
			    re->ri->read1_len, re->ri->initbp, revcmpl, Rflag);
			puts(output);
			free(output);

			if (Pflag) {
				output = output_pretty(re->ri->name, contig_name, rs->sfrp,
				    genome, genome_len,
				    (shrimp_mode == MODE_COLOUR_SPACE), re->ri->read1,
				    re->ri->read1_len, re->ri->initbp, revcmpl);
				puts("");
				puts(output);
				free(output);
			}
		}

		last_sfr = *rs->sfrp;
		free(rs->sfrp->dbalign);
		free(rs->sfrp->qralign);
		free(rs->sfrp);
	}
}

/*
 * Split our linked list of scores into segments based on read and send an
 * array of them to the appropriate space-specific output function.
 */
static void
generate_output(struct re_score *rs, bool revcmpl)
{
	struct re_score *byread;
	struct read_elem *re, *last_re;
	u_int i;

	byread = (struct re_score *)xmalloc(sizeof(*byread) * num_outputs);
	memset(byread, 0, sizeof(*byread) * num_outputs);

	i = 0;
	last_re = NULL;
	while (true) {
		re = (rs == NULL) ? NULL : rs->parent;

		if ((rs == NULL && byread[0].score > 0) ||
		    (re != NULL && last_re != NULL && re != last_re)) {
			if (shrimp_mode == MODE_HELICOS_SPACE)
				generate_output_dag(byread, i, revcmpl);
			else
				generate_output_lscs(byread, i, revcmpl);

			memset(byread, 0, sizeof(*byread) * num_outputs);
			i = 0;
		}

		if (rs == NULL)
			break;

		assert(re != NULL);
		assert(i < num_outputs);
		byread[i++] = *rs;
		last_re = re;

		rs = rs->next;
	}

	free(byread);
}

static void
final_pass_contig(struct re_score *scores, struct re_score *scores_revcmpl)
{

	if (Fflag)
		generate_output(scores, false);

	if (Cflag) {
		uint64_t before;

		before = rdtsc();
		reverse_complement(genome, NULL, genome_len);
		revcmpl_ticks += rdtsc() - before;
		generate_output(scores_revcmpl, true);
	}
}

static int
final_pass(char **files, u_int nfiles)
{
	struct re_score **scores, **scores_revcmpl;
	u_int i, j, hits = 0;
	fasta_t fasta;
	char *name, *sequence;

	scores = (struct re_score **)xmalloc(sizeof(*scores) * ncontigs);
	scores_revcmpl = (struct re_score**)xmalloc(sizeof(*scores) * ncontigs);
	memset(scores, 0, sizeof(*scores) * ncontigs);
	memset(scores_revcmpl, 0, sizeof(*scores_revcmpl) * ncontigs);

	/* For each contig, generate a linked list of hits from it. */
	for (i = 0; i < nreads; i++) {
		struct read_elem *re = OFFSET_TO_READ(i);

		if (re->ri->swhits == 0)
			continue;

		hits++;

		for (j = 1; j <= (u_int)re->ri->scores[0].score; j++) {
			struct re_score *rs = &re->ri->scores[j];
			u_int cn = rs->contig_num;

			assert(rs->score >= 0);
			assert(cn < ncontigs);
			assert(rs->parent == NULL);
			assert(rs->next == NULL);

			/* extra sanity */
			if (shrimp_mode != MODE_HELICOS_SPACE) {
				int thresh = get_sw_threshold(re, sw_vect_threshold, 1);
				assert(rs->score >= thresh);
			}

			/* prepend to linked list */
			if (rs->revcmpl) {
				rs->next =
				    scores_revcmpl[cn];
				scores_revcmpl[cn] = &re->ri->scores[j];
			} else {
				rs->next = scores[cn];
				scores[cn] = &re->ri->scores[j];
			}

			rs->parent = re;
		}
	}

	/* Now, do a final s-w for all reads of each contig and output them. */
	for (i = j = 0; i < nfiles; i++) {
		fprintf(stderr, "- Processing contig file [%s] (%u of %u)\n",
		    files[i], i + 1, nfiles);

		fasta = fasta_open(files[i], LETTER_SPACE);
		if (fasta == NULL) {
			fprintf(stderr, "error: failed to reload genome file "
			    "[%s] for final pass\n", files[i]);
			exit(1);
		}

		while (fasta_get_next(fasta, &name, &sequence)) {
			assert(j < ncontigs);

			if (strchr(name, '\t') != NULL || strchr(sequence, '\t') != NULL) {
				fprintf(stderr, "error: tabs are not permitted in fasta names "
				    "or sequences. Tag: [%s].\n", name);
				exit(1);
			}

			contig_name = name;
			genome_len = strlen(sequence);

			genome = fasta_sequence_to_bitfield(fasta, sequence);
			if (genome == NULL) {
				fprintf(stderr, "error: invalid contig sequence; contig: [%s]\n", name);
				exit(1);
			}
			genome_cs = NULL;

			fprintf(stderr, "  - Loaded %s letters from contig \"%s\"\n",
			    comma_integer(genome_len), name);

			free(sequence);
			final_pass_contig(scores[j], scores_revcmpl[j]);

			free(genome);
			free(contig_name);
			genome = genome_cs = NULL;
			contig_name = NULL;
			j++;
		}

		fasta_close(fasta);
	}

	free(scores);
	free(scores_revcmpl);

	return (hits);
}

static int
valid_spaced_seed()
{
	u_int seed_span, seed_weight;

	seed_span = strlen(spaced_seed);
	seed_weight = strchrcnt(spaced_seed, '1');

	if (seed_span < 1)
		return (0);

	if (strchrcnt(spaced_seed, '1') > MAX_SEED_WEIGHT)
		return (0);

	if (strchrcnt(spaced_seed, '0') != seed_span - seed_weight) 
		return (0);

	if (seed_weight < 1)
		return (0);

	return (1);
}

static void
scan_genomes_contig()
{
	uint64_t before;

	if (Fflag) {
		/*
		 * Do it forwards.
		 */
		before = rdtsc();
		scan(ncontigs, false);
		scan_ticks += (rdtsc() - before);
		reset_reads();
	}

	if (Cflag) {
		/*
		 * Do it reverse-complementwards.
		 */
		fprintf(stderr, "    - Processing reverse complement\n");

		before = rdtsc();
		if (shrimp_mode == MODE_COLOUR_SPACE)
			reverse_complement(genome, genome_cs, genome_len);
		else
			reverse_complement(genome, NULL, genome_len);
		revcmpl_ticks += (rdtsc() - before);

		before = rdtsc();
		scan(ncontigs, true);
		scan_ticks += (rdtsc() - before);
		reset_reads();
	}

	ncontigs++;
}

static int
scan_genomes(char **files, int nfiles)
{
	fasta_t fasta;
	char *name, *sequence;
	int i;

	for (i = 0; i < nfiles; i++) {
		fprintf(stderr, "- Processing contig file [%s] (%d of %d)\n",
		    files[i], i + 1, nfiles);

		fasta = fasta_open(files[i], LETTER_SPACE);
		if (fasta == NULL) {
			fprintf(stderr, "error: failed to reload genome file "
			    "[%s] for final pass\n", files[i]);
			return (1);
		}

		while (fasta_get_next(fasta, &name, &sequence)) { 
			contig_name = name;
			genome_len = strlen(sequence);
 
			if (strchr(name, '\t') != NULL || strchr(sequence, '\t') != NULL) {
				fprintf(stderr, "error: tabs are not permitted in fasta names "
				    "or sequences. Tag: [%s].\n", name);
				exit(1);
			}

			genome = fasta_sequence_to_bitfield(fasta, sequence);
			if (genome == NULL) {
				fprintf(stderr, "error: invalid contig sequence; contig: [%s]\n", name);
				exit(1);
			}
			if (shrimp_mode == MODE_COLOUR_SPACE)
				genome_cs = fasta_bitfield_to_colourspace(fasta, genome, genome_len);

			fprintf(stderr, "  - Loaded %s letters from contig \"%s\"\n",
			    comma_integer(genome_len), name);

			free(sequence);
			scan_genomes_contig();

			free(genome);
			if (shrimp_mode == MODE_COLOUR_SPACE)
				free(genome_cs);
			free(contig_name);
			genome = genome_cs = NULL;
			contig_name = NULL;
		}

		fasta_close(fasta);
	}

	return (0);
}

static u_int
set_window_lengths()
{
	uint32_t i;
	u_int max_window_len = 0;

	if (!IS_ABSOLUTE(window_len) && window_len < 100.0) {
		fprintf(stderr, "WARNING: window length is < 100%% of read "
		    "length - is this what you want?\n");
	}

	for (i = 0; i < nreads; i++) {
		struct read_elem *re = OFFSET_TO_READ(i);

		if (IS_ABSOLUTE(window_len)) {
			re->window_len = (uint16_t)-window_len;
		} else {
			int read_len = MAX(re->ri->read1_len, re->ri->read2_len);
			re->window_len = (uint16_t)ceil(((double)(read_len *
			    (window_len / 100.0))));
		}
		max_window_len = MAX(max_window_len, re->window_len);
	}

	return (max_window_len);
}

static void
print_statistics()
{
	double cellspersec, hz;
	uint64_t invocs, cells, reads_matched, total_matches, fasta_load_ticks, vticks;
	fasta_stats_t fs;
	uint32_t i;

	fs = fasta_stats();
	fasta_load_ticks = fs->total_ticks;
	free(fs);
	hz = cpuhz();

	reads_matched = total_matches = 0;
	for (i = 0; i < nreads; i++) {
		struct read_elem *re = OFFSET_TO_READ(i);

		reads_matched += (re->ri->final_matches == 0) ? 0 : 1;
		total_matches += re->ri->final_matches;
	}

	sw_vector_stats(&invocs, &cells, &vticks, &cellspersec);

	fprintf(stderr, "\nStatistics:\n");
	fprintf(stderr, "    Spaced Seed Scan:\n");
	fprintf(stderr, "        Run-time:               %.2f seconds\n",
	    (scan_ticks - vticks) / hz);
	fprintf(stderr, "        Total Kmers:            %s\n",
	    comma_integer(nkmers));
	fprintf(stderr, "        Minimal Reads/Kmer:     %s\n",
	    comma_integer(shortest_scanned_kmer_list));
	fprintf(stderr, "        Maximal Reads/Kmer:     %s\n",
	    comma_integer(longest_scanned_kmer_list));
	fprintf(stderr, "        Average Reads/Kmer:     %.2f\n",
	    (kmer_lists_scanned == 0) ? 0 :
	    (double)kmer_list_entries_scanned / (double)kmer_lists_scanned);
	fprintf(stderr, "\n");

	fprintf(stderr, "    Vector Smith-Waterman:\n");
	fprintf(stderr, "        Run-time:               %.2f seconds\n",
	    (cellspersec == 0) ? 0 : cells / cellspersec);
	fprintf(stderr, "        Invocations:            %s\n",
	    comma_integer(invocs));
	fprintf(stderr, "        Cells Computed:         %.2f million\n",
	    (double)cells / 1.0e6);
	fprintf(stderr, "        Cells per Second:       %.2f million\n",
	    cellspersec / 1.0e6);
	fprintf(stderr, "\n");

	if (shrimp_mode != MODE_HELICOS_SPACE) {
		if (shrimp_mode == MODE_COLOUR_SPACE)
			sw_full_cs_stats(&invocs, &cells, NULL, &cellspersec);
		else
			sw_full_ls_stats(&invocs, &cells, NULL, &cellspersec);

		fprintf(stderr, "    Scalar Smith-Waterman:\n");
		fprintf(stderr, "        Run-time:               %.2f seconds\n",
		    (cellspersec == 0) ? 0 : cells / cellspersec);
		fprintf(stderr, "        Invocations:            %s\n",
		    comma_integer(invocs));
		fprintf(stderr, "        Cells Computed:         %.2f million\n",
		    (double)cells / 1.0e6);
		fprintf(stderr, "        Cells per Second:       %.2f million\n",
		    cellspersec / 1.0e6);
	} else {
		struct dag_statistics *dsp = dag_get_statistics();

		fprintf(stderr, "    DAG Read-Read Alignment:\n");
		fprintf(stderr, "        Run-time:               %.2f seconds\n",
		    dsp->kmers_seconds);
		fprintf(stderr, "        Avg Kmers/Read Pair:    %.2f\n",
		    (double)dsp->kmers_total_kmers / (double)dsp->kmers_invocations);
		fasta_load_ticks = (uint64_t)MAX(0, (float)fasta_load_ticks -
		    (float)dsp->kmers_seconds * hz);
		fprintf(stderr, "\n");

		fprintf(stderr, "    DAG Pair-Genome Alignment:\n");
		fprintf(stderr, "        Run-time:               %.2f seconds\n",
		    dsp->aligner_seconds);
		fprintf(stderr, "        Invocations:            %s\n",
		    comma_integer(dsp->aligner_invocations));
		free(dsp);
	}
	fprintf(stderr, "\n");

	fprintf(stderr, "    Miscellaneous:\n");
	fprintf(stderr, "        Load Time:              %.2f seconds\n",
	    (double)fasta_load_ticks / hz);
	fprintf(stderr, "        Revcmpl. Time:          %.2f seconds\n",
	    (double)revcmpl_ticks / hz);
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
	    comma_integer(nduphits));
}

static void
usage(char *progname)
{
	char *slash, *seed;
	double vect_sw_threshold;

	switch (shrimp_mode) {
	case MODE_LETTER_SPACE:
		seed = DEF_SPACED_SEED_LS;
		vect_sw_threshold = DEF_SW_VECT_THRESHOLD;
		break;
	case MODE_COLOUR_SPACE:
		seed = DEF_SPACED_SEED_CS;
		vect_sw_threshold = DEF_SW_VECT_THRESHOLD;
		break;
	case MODE_HELICOS_SPACE:
		seed = DEF_SPACED_SEED_DAG;
		vect_sw_threshold = DEF_SW_VECT_THRESHOLD_DAG;
		break;
	default:
		assert(0);
	}

	slash = strrchr(progname, '/');
	if (slash != NULL)
		progname = slash + 1;

	fprintf(stderr, "usage: %s [parameters] [options] "
	    "reads_file genome_file1 genome_file2...\n", progname);

	fprintf(stderr, "Parameters:\n");

	fprintf(stderr,
	    "    -s    Spaced Seed                             (default: %s)\n",
	    seed);

	fprintf(stderr,
	    "    -n    Seed Matches per Window                 (default: %d)\n",
	    DEF_NUM_MATCHES);

	fprintf(stderr,
	    "    -t    Seed Hit Taboo Length                   (default: %d)\n",
	    DEF_HIT_TABOO_LEN);

	if (shrimp_mode != MODE_HELICOS_SPACE) {
		fprintf(stderr,
		    "    -9    Seed Generation Taboo Length            (default: %d)\n",
		    DEF_SEED_TABOO_LEN);
	}

	fprintf(stderr,
	    "    -w    Seed Window Length                      (default: "
	    "%.02f%%)\n", DEF_WINDOW_LEN);

	fprintf(stderr,
	    "    -o    Maximum Hits per Read                   (default: %d)\n",
	    DEF_NUM_OUTPUTS);

	fprintf(stderr,
	    "    -r    Maximum Read Length                     (default: %d)\n",
	    DEF_MAX_READ_LEN);

	fprintf(stderr,
	    "    -d    Kmer Std. Deviation Limit               (default: %d%s)"
	    "\n", DEF_KMER_STDDEV_LIMIT,
	    (DEF_KMER_STDDEV_LIMIT < 0) ? " [None]" : "");

	if (shrimp_mode == MODE_HELICOS_SPACE) {
		fprintf(stderr, "\n");
		fprintf(stderr,
		    "    -p    DAG Epsilon                             (default: %d)\n",
		    DEF_DAG_EPSILON);

		fprintf(stderr,
		    "    -1    DAG Read Match                          (default: %d)\n",
		    DEF_DAG_READ_MATCH);

		fprintf(stderr,
		    "    -y    DAG Read Mismatch                       (default: %d)\n",
		    DEF_DAG_READ_MISMATCH);

		fprintf(stderr,
		    "    -z    DAG Read Gap                            (default: %d)\n",
		    DEF_DAG_READ_GAP);

		fprintf(stderr,
		    "    -a    DAG Reference Match                     (default: %d)\n",
		    DEF_DAG_REF_MATCH);

		fprintf(stderr,
		    "    -b    DAG Reference Mismatch                  (default: %d)\n",
		    DEF_DAG_REF_MISMATCH);

		fprintf(stderr,
		    "    -c    DAG Reference Half Match                (default: %d)\n",
		    DEF_DAG_REF_HALF_MATCH);

		fprintf(stderr,
		    "    -j    DAG Reference Neither Match             (default: %d)\n",
		    DEF_DAG_REF_NEITHER_MATCH);

		fprintf(stderr,
		    "    -k    DAG Reference Match,Deletion            (default: %d)\n",
		    DEF_DAG_REF_MATCH_DELETION);

		fprintf(stderr,
		    "    -l    DAG Reference Mismatch,Deletion         (default: %d)\n",
		    DEF_DAG_REF_MISMATCH_DELETION);

		fprintf(stderr,
		    "    -u    DAG Reference Error,Insertion           (default: %d)\n",
		    DEF_DAG_REF_ERROR_INSERTION);

		fprintf(stderr,
		    "    -2    DAG Reference Weighted Threshold        (default: %.3f)\n",
		    DEF_DAG_REF_WEIGHTED_THRESHOLD);
	}

	fprintf(stderr, "\n");
	fprintf(stderr,
	    "    -m    S-W Match Value                         (default: %d)\n",
	    (shrimp_mode == MODE_HELICOS_SPACE) ? DEF_MATCH_VALUE_DAG : DEF_MATCH_VALUE);

	fprintf(stderr,
	    "    -i    S-W Mismatch Value                      (default: %d)\n",
	    (shrimp_mode == MODE_HELICOS_SPACE) ? DEF_MISMATCH_VALUE_DAG : DEF_MISMATCH_VALUE);

	fprintf(stderr,
	    "    -g    S-W Gap Open Penalty (Reference)        (default: %d)\n",
	    (shrimp_mode == MODE_HELICOS_SPACE) ? DEF_A_GAP_OPEN_DAG : DEF_A_GAP_OPEN);

	fprintf(stderr,
	    "    -q    S-W Gap Open Penalty (Query)            (default: %d)\n",
	    (shrimp_mode == MODE_HELICOS_SPACE) ? DEF_B_GAP_OPEN_DAG : DEF_B_GAP_OPEN);

	fprintf(stderr,
	    "    -e    S-W Gap Extend Penalty (Reference)      (default: %d)\n",
	    (shrimp_mode == MODE_HELICOS_SPACE) ? DEF_A_GAP_EXTEND_DAG : DEF_A_GAP_EXTEND);

	fprintf(stderr,
	    "    -f    S-W Gap Extend Penalty (Query)          (default: %d)\n",
	    (shrimp_mode == MODE_HELICOS_SPACE) ? DEF_B_GAP_EXTEND_DAG : DEF_B_GAP_EXTEND);

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

	exit(1);
}

int
main(int argc, char **argv)
{
	char *reads_file, *progname, *optstr;
	char **genome_files;
	int ch, ret, ngenome_files;
	u_int max_window_len;
	bool a_gap_open_set, b_gap_open_set;
	bool a_gap_extend_set, b_gap_extend_set;

	set_mode_from_argv(argv);

	a_gap_open_set = b_gap_open_set = a_gap_extend_set = b_gap_extend_set = false;

	/* set the appropriate defaults based on mode */
	switch (shrimp_mode) {
	case MODE_COLOUR_SPACE:
		spaced_seed = DEF_SPACED_SEED_CS;
		optstr = "s:n:t:9:w:o:r:d:m:i:g:q:e:f:x:h:v:BCFPR";
		break;
	case MODE_LETTER_SPACE:
		spaced_seed = DEF_SPACED_SEED_LS;
		optstr = "s:n:t:9:w:o:r:d:m:i:g:q:e:f:h:BCFPR";
		break;
	case MODE_HELICOS_SPACE:
		spaced_seed = DEF_SPACED_SEED_DAG;
		optstr = "s:n:t:w:o:r:d:p:1:y:z:a:b:c:j:k:l:u:2:m:i:g:q:e:f:v:BCFPR";
		match_value = DEF_MATCH_VALUE_DAG;
		mismatch_value = DEF_MISMATCH_VALUE_DAG;
		a_gap_open = DEF_A_GAP_OPEN_DAG;
		b_gap_open = DEF_B_GAP_OPEN_DAG;
		a_gap_extend = DEF_A_GAP_EXTEND_DAG;
		b_gap_extend = DEF_B_GAP_EXTEND_DAG;
		sw_vect_threshold = DEF_SW_VECT_THRESHOLD_DAG;
		break;
	default:
		assert(0);
	}

	fprintf(stderr, "--------------------------------------------------"
	    "------------------------------\n");
	fprintf(stderr, "rmapper: %s.\nSHRiMP %s [%s]\n", get_mode_string(),
	    SHRIMP_VERSION_STRING, get_compiler());
	fprintf(stderr, "--------------------------------------------------"
	    "------------------------------\n");

	progname = argv[0];

	while ((ch = getopt(argc, argv, optstr)) != -1) {
		switch (ch) {
		case 's':
			spaced_seed = xstrdup(optarg);
			break;
		case 'n':
			num_matches = atoi(optarg);
			break;
		case 't':
			hit_taboo_len = atoi(optarg);
			break;
		case '9':
			assert(shrimp_mode != MODE_HELICOS_SPACE);
			seed_taboo_len = atoi(optarg);
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
		case 'd':
			kmer_stddev_limit = atoi(optarg);
			break;
		case 'p':
			dag_epsilon = atoi(optarg);
			break;
		case '1':
			dag_read_match = atoi(optarg);
			break;
		case 'y':
			dag_read_mismatch = atoi(optarg);
			break;
		case 'z':
			dag_read_gap = atoi(optarg);
			break;
		case 'a':
			dag_ref_match = atoi(optarg);
			break;
		case 'b':
			dag_ref_mismatch = atoi(optarg);
			break;
		case 'c':
			dag_ref_half_match = atoi(optarg);
			break;
		case 'j':
			dag_ref_neither_match = atoi(optarg);
			break;
		case 'k':
			dag_ref_match_deletion = atoi(optarg);
			break;
		case 'l':
			dag_ref_mismatch_deletion = atoi(optarg);
			break;
		case 'u':
			dag_ref_error_insertion = atoi(optarg);
			break;
		case '2':
			dag_ref_weighted_thresh = atof(optarg);
			break;
		case 'm':
			match_value = atoi(optarg);
			break;
		case 'i':
			mismatch_value = atoi(optarg);
			break;
		case 'g':
			a_gap_open_set = true;
			a_gap_open = atoi(optarg);
			break;
		case 'q':
			b_gap_open_set = true;
			b_gap_open = atoi(optarg);
			break;
		case 'e':
			a_gap_extend_set = true;
			a_gap_extend = atoi(optarg);
			break;
		case 'f':
			b_gap_extend_set = true;
			b_gap_extend = atoi(optarg);
			break;
		case 'x':
			assert(shrimp_mode == MODE_COLOUR_SPACE);
			xover_penalty = atoi(optarg);
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
		case 'P':
			Pflag = true;
			break;
		case 'R':
			Rflag = true;
			break;
		default:
			usage(progname);
		}
	}
	argc -= optind;
	argv += optind;

	if (argc < 2) {
		fprintf(stderr, "error: %sgenome_file(s) not specified\n",
		    (argc == 0) ? "reads_file, " : "");
		usage(progname);
	}

	reads_file    = argv[0];
	genome_files  = &argv[1];
	ngenome_files = argc - 1;

	if (!Cflag && !Fflag)
		Cflag = Fflag = true;

	if (shrimp_mode == MODE_LETTER_SPACE)
		sw_vect_threshold = sw_full_threshold;

	if (!valid_spaced_seed()) {
		fprintf(stderr, "error: invalid spaced seed\n");
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

	if (a_gap_open > 0 || b_gap_open > 0) {
		fprintf(stderr, "error: invalid gap open penalty\n");
		exit(1);
	}

	if (a_gap_extend > 0 || b_gap_extend > 0) {
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

	if (a_gap_open_set && !b_gap_open_set ||
	    a_gap_extend_set && !b_gap_extend_set)
		fputc('\n', stderr);
	if (a_gap_open_set && !b_gap_open_set) {
		fprintf(stderr, "Notice: Gap open penalty set for reference but not query; assuming symmetry.\n");
		b_gap_open = a_gap_open;
	}
	if (a_gap_extend_set && !b_gap_extend_set) {
		fprintf(stderr, "Notice: Gap extend penalty set for reference but not query; assuming symmetry.\n");
		b_gap_extend = a_gap_extend;
	}
	if (a_gap_open_set && !b_gap_open_set ||
	    a_gap_extend_set && !b_gap_extend_set)
		fputc('\n', stderr);

	fprintf(stderr, "Settings:\n");
	fprintf(stderr, "    Spaced Seed:                          %s\n",
	    spaced_seed);
	fprintf(stderr, "    Spaced Seed Span:                     %u\n",
	    (u_int)strlen(spaced_seed));
	fprintf(stderr, "    Spaced Seed Weight:                   %u\n",
	    strchrcnt(spaced_seed, '1'));
	fprintf(stderr, "    Seed Matches per Window:              %u\n",
	    num_matches);
	fprintf(stderr, "    Seed Hit Taboo Length:                %u\n",
	    hit_taboo_len);

	if (shrimp_mode != MODE_HELICOS_SPACE) {
		fprintf(stderr, "    Seed Generation Taboo Length:         %u\n",
		    seed_taboo_len);
	}

	if (IS_ABSOLUTE(window_len)) {
		fprintf(stderr, "    Seed Window Length:                   "
		    "%u\n", (u_int)-window_len);
	} else {
		fprintf(stderr, "    Seed Window Length:                   "
		    "%.02f%%\n", window_len);
	}

	fprintf(stderr, "    Maximum Hits per Read:                %u\n",
	    num_outputs);
	fprintf(stderr, "    Kmer Std. Deviation Limit:            %d%s\n",
	    kmer_stddev_limit, (kmer_stddev_limit < 0) ? " (None)" : "");

	if (shrimp_mode == MODE_HELICOS_SPACE) {
		fprintf(stderr, "\n");
		fprintf(stderr, "    DAG Epsilon:                          "
		    "%u\n", dag_epsilon);
		fprintf(stderr, "    DAG Read Match:                       "
		    "%d\n", dag_read_match);
		fprintf(stderr, "    DAG Read Mismatch:                    "
		    "%d\n", dag_read_mismatch);
		fprintf(stderr, "    DAG Read Gap:                         "
		    "%d\n", dag_read_gap);
		fprintf(stderr, "    DAG Reference Match:                  "
		    "%d\n", dag_ref_match);
		fprintf(stderr, "    DAG Reference Mismatch:               "
		    "%d\n", dag_ref_mismatch);
		fprintf(stderr, "    DAG Reference Half Match:             "
		    "%d\n", dag_ref_half_match);
		fprintf(stderr, "    DAG Reference Neither Match:          "
		    "%d\n", dag_ref_neither_match);
		fprintf(stderr, "    DAG Reference Match,Deletion:         "
		    "%d\n", dag_ref_match_deletion);
		fprintf(stderr, "    DAG Reference Mismatch,Deletion:      "
		    "%d\n", dag_ref_mismatch_deletion);
		fprintf(stderr, "    DAG Reference Error,Insertion:        "
		    "%d\n", dag_ref_error_insertion);
		fprintf(stderr, "    DAG Reference Weighted Threshold:     "
		    "%.03f\n", dag_ref_weighted_thresh);
	}

	fprintf(stderr, "\n");
	fprintf(stderr, "    S-W Match Value:                      "
	    "%d\n", match_value);
	fprintf(stderr, "    S-W Mismatch Value:                   "
	    "%d\n", mismatch_value);
	fprintf(stderr, "    S-W Gap Open Penalty (Ref):           "
	    "%d\n", a_gap_open);
	fprintf(stderr, "    S-W Gap Open Penalty (Qry):           "
	    "%d\n", b_gap_open);
	fprintf(stderr, "    S-W Gap Extend Penalty (Ref):         "
	    "%d\n", a_gap_extend);
	fprintf(stderr, "    S-W Gap Extend Penalty (Qry):         "
	    "%d\n", b_gap_extend);

	if (shrimp_mode == MODE_COLOUR_SPACE) {
		fprintf(stderr, "    S-W Crossover Penalty:                "
		    "%d\n", xover_penalty);
	}

	if (shrimp_mode == MODE_COLOUR_SPACE || shrimp_mode == MODE_HELICOS_SPACE) {
		if (IS_ABSOLUTE(sw_vect_threshold)) {
			fprintf(stderr, "    S-W Vector Hit Threshold:     "
			    "        %u\n", (u_int)-sw_vect_threshold);
		} else {
			fprintf(stderr, "    S-W Vector Hit Threshold:     "
			    "        %.02f%%\n", sw_vect_threshold);
		}
	}
	if (shrimp_mode == MODE_COLOUR_SPACE) {
		if (IS_ABSOLUTE(sw_full_threshold)) {
			fprintf(stderr, "    S-W Full Hit Threshold:       "
			    "        %u\n", (u_int)-sw_full_threshold);
		} else {
			fprintf(stderr, "    S-W Full Hit Threshold:       "
			    "        %.02f%%\n", sw_full_threshold);
		}
	}
	if (shrimp_mode == MODE_LETTER_SPACE) {
		if (IS_ABSOLUTE(sw_full_threshold)) {
			fprintf(stderr, "    S-W Hit Threshold:            "
			    "        %u\n", (u_int)-sw_full_threshold);
		} else {
			fprintf(stderr, "    S-W Hit Threshold:            "
			    "        %.02f\n", sw_full_threshold);
		}
	}
	fprintf(stderr, "\n");

	dag_setup(dag_read_match, dag_read_mismatch, dag_read_gap, dag_ref_match,
	    dag_ref_mismatch, dag_ref_half_match, dag_ref_neither_match,
	    dag_ref_match_deletion, dag_ref_mismatch_deletion,
	    dag_ref_error_insertion);

	if (!load_reads(reads_file))
		exit(1);

	fprintf(stderr, "  - Configuring window lengths...\n");
	max_window_len = set_window_lengths();
	fprintf(stderr, "  - Maximum window length: %u\n", max_window_len);
		

	if (sw_vector_setup(max_window_len * 2, longest_read_len, a_gap_open,
	    a_gap_extend, b_gap_open, b_gap_extend, match_value, mismatch_value,
	    shrimp_mode == MODE_COLOUR_SPACE, false)) {
		fprintf(stderr, "failed to initialise vector "
		    "Smith-Waterman (%s)\n", strerror(errno));
		exit(1);
	}

	if (shrimp_mode == MODE_COLOUR_SPACE) {
/* XXX - a vs. b gap */
		ret = sw_full_cs_setup(max_window_len * 2, longest_read_len,
		    a_gap_open, a_gap_extend, match_value, mismatch_value,
		    xover_penalty, false);
	} else {
		ret = sw_full_ls_setup(max_window_len * 2, longest_read_len,
		    a_gap_open, a_gap_extend, b_gap_open, b_gap_extend,
		    match_value, mismatch_value, false);
	}
	if (ret) {
		fprintf(stderr, "failed to initialise scalar "
		    "Smith-Waterman (%s)\n", strerror(errno));
		exit(1);
	}

	readmap_prune();

	if (scan_genomes(genome_files, ngenome_files))
		exit(1);

	fprintf(stderr, "\nGenerating output...\n");
	final_pass(genome_files, ngenome_files);

	print_statistics();

	return (0);
}
