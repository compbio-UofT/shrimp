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

#include <sys/types.h>
#include <sys/stat.h>
#include <sys/time.h>

#include "../common/fasta.h"
#include "../common/sw-full-common.h"
#include "../common/sw-full-cs.h"
#include "../common/sw-full-ls.h"
#include "../common/sw-vector.h"
#include "../common/output.h"
#include "../common/util.h"
#include "../common/version.h"

#include "rmapper.h"

/* External parameters */
static char *spaced_seed	= DEF_SPACED_SEED;
static u_int window_len		= DEF_WINDOW_LEN;
static u_int num_matches	= DEF_NUM_MATCHES;
static u_int taboo_len		= DEF_TABOO_LEN;
static u_int num_outputs	= DEF_NUM_OUTPUTS;
static u_int max_read_len	= DEF_MAX_READ_LEN;
static int   kmer_stddev_limit  = DEF_KMER_STDDEV_LIMIT;

static int   match_value	= DEF_MATCH_VALUE;
static int   mismatch_value	= DEF_MISMATCH_VALUE;
static int   gap_open		= DEF_GAP_OPEN;
static int   gap_extend		= DEF_GAP_EXTEND;
static int   xover_penalty	= DEF_XOVER_PENALTY;
static u_int sw_vect_threshold	= DEF_SW_VECT_THRESHOLD;
static u_int sw_full_threshold	= DEF_SW_FULL_THRESHOLD;

/* Genomic sequence, stored in 32-bit words, first is in the LSB */
static uint32_t *genome;
static uint32_t	 genome_len;
static char     *contig_name;			/* name of contig in genome */

/*
 * Genomic sequence, in letter space, stored in 32-bit words, first is in LSB
 *
 * NB: Colour space only
 */
static uint32_t *genome_ls;

/* Reads */
static struct read_elem  *reads;		/* list of reads */
static struct read_node **readmap;		/* kmer to read map */
static uint32_t		 *readmap_len;		/* entries in each readmap[n] */

static u_int	nreads;				/* total reads loaded */
static u_int	nkmers;				/* total kmers of reads loaded*/

static u_int	longest_read_len;		/* longest read we've seen */

/* Flags */
static int Bflag = false;			/* print a progress bar */
static int Pflag = false;			/* pretty print results */
static int Rflag = false;			/* add read seq to output */

/* Scan stats */
static uint64_t shortest_scanned_kmer_list;
static uint64_t longest_scanned_kmer_list;
static uint64_t kmer_lists_scanned;
static uint64_t kmer_list_entries_scanned;

#ifdef USE_COLOURS
const bool use_colours = true;
#else
const bool use_colours = false;
#endif

#define PROGRESS_BAR(_a, _b, _c, _d)	\
    if (Bflag) progress_bar((_a), (_b), (_c), (_d))

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

/* update our score min-heap */
static void
save_score(struct read_elem *re, int score, int index, int contig_num,
    bool revcmpl)
{
	struct re_score *scores;

	scores = re->scores;

	if (score < scores[1].score)
		return;

	scores[1].score = score;
	scores[1].index = index;
	scores[1].contig_num = contig_num;
	scores[1].revcmpl = revcmpl;
	reheap(scores, 1);
}

/* compress the given kmer into an index in 'readmap' according to the seed */
static uint32_t
kmer_to_mapidx(uint32_t *kmer)
{
	static int seed_span;

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

	fprintf(stderr, "Pruning kmers... (mean: %f, stddev: %f, "
	    "stddev limit: +/- %f)\n", mean, stddev,
	    kmer_stddev_limit * stddev);
	if (mean < 1.0)
		fprintf(stderr, "WARNING: low mean - are you sure you want to "
		    "prune kmers?\n");

	for (i = 0; i < j; i++) {
		if (readmap_len[i] > mean + (kmer_stddev_limit * stddev))
			readmap[i] = NULL;
	}
}

/* scan the genome by kmers, running S-W as needed, and updating scores */
static void
scan(int contig_num, bool revcmpl)
{
	struct read_node *rn;
	struct read_elem *re;
	uint32_t *kmer;
	uint32_t i, j, idx, mapidx;
	uint32_t base, skip, score;
	uint32_t goff, glen;
	int prevhit, seed_span;

	seed_span = strlen(spaced_seed);

	PROGRESS_BAR(stderr, 0, 0, 10);

	kmer = xmalloc(sizeof(kmer[0]) * BPTO32BW(seed_span));
	memset(kmer, 0, sizeof(kmer[0]) * BPTO32BW(seed_span));
	for (i = 0; i < seed_span - 1; i++)
		bitfield_prepend(kmer, seed_span, EXTRACT(genome, i));

	skip = 0;
	for (i = seed_span - 1; i < genome_len; i++) {
		PROGRESS_BAR(stderr, i, genome_len, 10);

		base = EXTRACT(genome, i);
		bitfield_prepend(kmer, seed_span, base);

		/*
		 * Ignore kmers that contain an 'N'.
		 */
		if (base > 3)
			skip = seed_span - 1;

		if (skip > 0) {
			skip--;
			continue;
		}

		mapidx = kmer_to_mapidx(kmer);
		idx = i - (seed_span - 1);

		j = 0;
		for (rn = readmap[mapidx]; rn != NULL; rn = rn->next) {
			re = rn->read;

			prevhit = re->hits[re->prev_hit];
			if ((idx - prevhit) <= taboo_len && prevhit != -1)
				continue;

			re->hits[re->next_hit] = idx;
			re->prev_hit = re->next_hit;
			re->next_hit = (re->next_hit + 1) % num_matches;

			if ((idx - re->hits[re->next_hit]) < window_len &&
			    (i - re->last_swhit_idx) >= window_len &&
			    re->hits[re->next_hit] != -1) {
				if (i < window_len)
					goff = 0;
				else
					goff = i - window_len;

				if (goff + (window_len * 2) >= genome_len)
					glen = genome_len - goff;
				else
					glen = window_len * 2;

				if (use_colours) {
					score = sw_vector(genome, goff, glen,
					    re->read, re->read_len,
					    genome_ls, re->initbp);
				} else {
					score = sw_vector(genome, goff, glen,
					    re->read, re->read_len, NULL, -1);
				}

				if (score >= sw_vect_threshold) {
					save_score(re, score, i, contig_num,
					    revcmpl);
					re->last_swhit_idx = i;
					re->swhits++;
				}
			}

			j++;
		}

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

static void
load_genome_helper(int base, ssize_t offset, int isnewentry, char *name,
    int initbp)
{
	static bool first = true;
	static int lastbp;

	/* shut up icc */
	(void)name;
	(void)lastbp;
	(void)initbp;

	if (base == FASTA_DEALLOC)
		return;

	/* handle initial call to alloc resources */
	if (base == FASTA_ALLOC) {
		first = true;
		if (genome != NULL)
			free(genome);
		genome = xmalloc(sizeof(genome[0]) * BPTO32BW(offset));
		genome_len = 0;
		memset(genome, 0, sizeof(genome[0]) * BPTO32BW(offset));

		if (use_colours) {
			if (genome_ls != NULL)
				free(genome_ls);
			genome_ls =
			    xmalloc(sizeof(genome_ls[0]) * BPTO32BW(offset));
			memset(genome_ls, 0,
			    sizeof(genome_ls[0]) * BPTO32BW(offset));
			lastbp = BASE_T;
		}
		return;
	}

	if (first) {
		if (contig_name != NULL)
			free(contig_name);
		contig_name = xstrdup(name);
	}

	if (isnewentry && !first) {
		fprintf(stderr, "error: genome file consists of more than one "
		    "contig!\n");
		fprintf(stderr, "       rmapper expects one contig per fasta "
		    "file - use 'splittigs' to break files up.\n");
		exit(1);
	}

	assert(base >= 0 && base <= 15);
	assert(offset == genome_len);

	if (use_colours) {
		bitfield_append(genome, offset, lstocs(lastbp, base));
		genome_len++;
		bitfield_append(genome_ls, offset, base);
		lastbp = base;
	} else {
		bitfield_append(genome, offset, base);
		genome_len++;
	}

	first = false;
}

static int
load_genome(const char *file)
{
	ssize_t ret;

	ret = load_fasta(file, load_genome_helper, LETTER_SPACE);

	if (ret != -1)
		fprintf(stderr, "Loaded %u letters from contig \"%s\" [%s]\n",
		    (unsigned int)ret, contig_name, file);

	return (ret == -1);
}

static void
load_reads_helper(int base, ssize_t offset, int isnewentry, char *name,
    int initbp)
{
	static struct read_elem *re;
	static uint32_t **past_kmers;
	static uint32_t *kmer;
	static uint32_t cnt, npast_kmers, skip;

	size_t len;
	int i, seed_span, seed_weight;

	/* shut up icc */
	(void)offset;

	if (base == FASTA_DEALLOC)
		return;

	seed_span = strlen(spaced_seed);
	seed_weight = strchrcnt(spaced_seed, '1');

	if (kmer == NULL)
		kmer = xmalloc(sizeof(kmer[0]) * BPTO32BW(seed_span));

	if (past_kmers == NULL) {
		len = sizeof(past_kmers[0]) * (max_read_len - seed_span + 1);
		past_kmers = xmalloc(len);

		len = sizeof(kmer[0]) * BPTO32BW(seed_span);
		for (i = 0; i < (max_read_len - seed_span + 1); i++) {
			past_kmers[i] = xmalloc(len);
			memset(past_kmers[i], 0, len);
		}
	}

	/* handle initial call to alloc resources */
	if (base == FASTA_ALLOC) {
		len = sizeof(readmap[0]) * power(4, seed_weight);
		readmap = xmalloc(len);
		memset(readmap, 0, len);

		len = sizeof(readmap_len[0]) * power(4, seed_weight);
		readmap_len = xmalloc(len);
		memset(readmap_len, 0, len);
		return;
	}

	if (isnewentry) {
		len = sizeof(struct read_elem) +
		    (sizeof(re->hits[0]) * num_matches);
		re = xmalloc(len);
		memset(re, 0, len);

		re->name = strdup(name);
		memset(re->hits, -1, sizeof(re->hits[0]) * num_matches);

		len = sizeof(re->read[0]) * BPTO32BW(max_read_len);
		re->read = xmalloc(len);
		memset(re->read, 0, len);

		len = sizeof(struct re_score) * (num_outputs + 1);
		re->scores = xmalloc(len);
		memset(re->scores, 0, len);
		re->scores[0].score = num_outputs;

		/* init scores to negatives so we only need percolate down */
		for (i = 1; i <= re->scores[0].score; i++)
			re->scores[i].score = 0x80000000 + i;

		cnt = npast_kmers = skip = 0;
		memset(kmer, 0, sizeof(kmer[0]) * BPTO32BW(seed_span));
		re->next = reads;
		reads = re;
		nreads++;

		re->initbp = initbp;
	}

	assert(re != NULL);

	if (re->read_len == max_read_len) {
		fprintf(stderr, "error: read [%s] exceeds %u characters\n"
		    "please increase the maximum read length parameter\n",
		    re->name, max_read_len);
		exit(1);
	}

	bitfield_append(re->read, re->read_len, base);
	re->read_len++;
	longest_read_len = MAX(re->read_len, longest_read_len);

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
	if (use_colours && isnewentry)
		return;

	bitfield_prepend(kmer, seed_span, base);

	/*
	 * Ignore kmers that contain an 'N'.
	 */
	if (base > 3)
		skip = seed_span - 1;

	if (skip > 0) {
		skip--;
		return;
	}

	if (++cnt >= seed_span) {
		struct read_node *rn;
		uint32_t mapidx;
		bool unique = true;

		/*
		 * Ensure that this kmer is unique within the read. If it's
		 * not, we don't want to point from one kmer to the same
		 * read multiple times.
		 *
		 * This will only hit about 0.03% of 8mers in 24mer reads.
	 	 */
		assert(npast_kmers < (max_read_len - seed_span + 1));

		len = sizeof(kmer[0]) * BPTO32BW(seed_span);
		for (i = 0; i < npast_kmers; i++) {
			if (memcmp(kmer, past_kmers[i], len) == 0) {
				unique = false;
				break;
			}
		}

		if (unique) {
			rn = xmalloc(sizeof(struct read_node));

			mapidx = kmer_to_mapidx(kmer);
			
			rn->read = re;
			rn->next = readmap[mapidx];
			readmap[mapidx] = rn;
			readmap_len[mapidx]++;
			memcpy(past_kmers[npast_kmers++], kmer, len);
			nkmers++;
		}
	}
}

static int
load_reads(const char *file)
{
	ssize_t ret;

	if (use_colours)
		ret = load_fasta(file, load_reads_helper, COLOUR_SPACE);
	else
		ret = load_fasta(file, load_reads_helper, LETTER_SPACE);

	if (ret != -1)
		fprintf(stderr, "Loaded %u %s in %u reads (%u kmers)\n",
		    (unsigned int)ret, (use_colours) ? "colours" : "letters",
		    nreads, nkmers);

	return (ret == -1);
}

static void
generate_output(struct re_score *rs, bool revcmpl)
{
	static bool firstcall = true;

	struct sw_full_results sfr;
	char *dbalign, *qralign;
	struct read_elem *re;
	uint32_t goff, glen;

	for (; rs != NULL; rs = rs->next) {
		re = rs->parent;

		if (rs->index < window_len)
			goff = 0;
		else
			goff = rs->index - window_len;

		assert(goff < genome_len);

		if (goff + (window_len * 2) >= genome_len)
			glen = genome_len - goff;
		else
			glen = window_len * 2;

		if (use_colours) {
			sw_full_cs(genome_ls, goff, glen, re->read,
			    re->read_len, re->initbp, rs->score, &dbalign,
			    &qralign, &sfr);
		} else {
			sw_full_ls(genome, goff, glen, re->read, re->read_len,
			    rs->score, &dbalign, &qralign, &sfr);
		}

		if (sfr.score < sw_full_threshold)
			continue;

		re->final_matches++;

		if (!firstcall)
			fputc('\n', stdout);
		firstcall = false;

		output_normal(stdout, re->name, contig_name, &sfr, genome_len,
		    goff, use_colours, re->read, re->read_len, re->initbp,
		    revcmpl, Rflag);
		fputc('\n', stdout);

		if (Pflag) {
			fputc('\n', stdout);
			output_pretty(stdout, re->name, contig_name, &sfr,
			    dbalign, qralign,
			    (use_colours) ? genome_ls : genome, genome_len,
			    goff, use_colours, re->read, re->read_len,
			    re->initbp, revcmpl);
		}
	}
}

static int
final_pass(char **files, int nfiles)
{
	struct {
		struct re_score *scores;
		struct re_score *revcmpl_scores;
	} *contig_list;
	struct read_elem *re;
	int i, hits = 0;

	contig_list = xmalloc(sizeof(*contig_list) * nfiles);
	memset(contig_list, 0, sizeof(*contig_list) * nfiles);

	/* For each contig, generate a linked list of hits from it. */
	for (re = reads; re != NULL; re = re->next) {
		if (re->swhits == 0)
			continue;

		hits++;
		
		for (i = 1; i <= re->scores[0].score; i++) {
			int cn = re->scores[i].contig_num;

			if (re->scores[i].score < 0)
				continue;

			assert(cn >= 0 && cn < nfiles);
			assert(re->scores[i].parent == NULL);
			assert(re->scores[i].next == NULL);
			assert(re->scores[i].score >= sw_vect_threshold);

			if (re->scores[i].revcmpl) {
				re->scores[i].next =
				    contig_list[cn].revcmpl_scores;
				contig_list[cn].revcmpl_scores = &re->scores[i];
			} else {
				re->scores[i].next =
				    contig_list[cn].scores;
				contig_list[cn].scores = &re->scores[i];
			}

			re->scores[i].parent = re;
		}
	}

	/* Now, do a final s-w for all reads of each contig and output them. */
	for (i = 0; i < nfiles; i++) {
		if (load_genome(files[i])) {
			fprintf(stderr, "error: failed to reload genome file "
			    "[%s] for final pass\n", files[i]);
			exit(1);
		}
		
		generate_output(contig_list[i].scores, false);
		if (use_colours)
			reverse_complement(genome_ls, genome, genome_len);
		else
			reverse_complement(genome, NULL, genome_len);
		generate_output(contig_list[i].revcmpl_scores, true);
	}

	free(contig_list);

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

static int
scan_genomes(char **files, int nfiles, double *scantime, double *loadtime,
    double *revcmpltime)
{
	struct timeval tv1, tv2;
	double stime, ltime, rtime;
	int i;

	stime = ltime = rtime = 0;

	for (i = 0; i < nfiles; i++) {
		fprintf(stderr, "Processing contig file [%s] (%d of %d)\n",
		    files[i], i + 1, nfiles);

		gettimeofday(&tv1, NULL);
		if (load_genome(files[i]))
			return (1);
		gettimeofday(&tv2, NULL);

		ltime += ((double)tv2.tv_sec + (double)tv2.tv_usec / 1.0e6) -
		    ((double)tv1.tv_sec + (double)tv1.tv_usec / 1.0e6);

		/*
		 * Do it forwards.
		 */
		gettimeofday(&tv1, NULL);
		scan(i, false);
		gettimeofday(&tv2, NULL);

		stime += ((double)tv2.tv_sec + (double)tv2.tv_usec / 1.0e6) -
		    ((double)tv1.tv_sec + (double)tv1.tv_usec / 1.0e6);

		/*
		 * Do it reverse-complementwards.
		 */
		fprintf(stderr, "Processing reverse complement of contig file "
		    "[%s] (%d of %d)\n", files[i], i + 1, nfiles);

		gettimeofday(&tv1, NULL);
		if (use_colours)
			reverse_complement(genome_ls, genome, genome_len);
		else
			reverse_complement(genome, NULL, genome_len);
		gettimeofday(&tv2, NULL);

		rtime += ((double)tv2.tv_sec + (double)tv2.tv_usec / 1.0e6) -
		    ((double)tv1.tv_sec + (double)tv1.tv_usec / 1.0e6);

		gettimeofday(&tv1, NULL);
		scan(i, true);
		gettimeofday(&tv2, NULL);

		stime += ((double)tv2.tv_sec + (double)tv2.tv_usec / 1.0e6) -
		    ((double)tv1.tv_sec + (double)tv1.tv_usec / 1.0e6);
	}

	if (scantime != NULL)
		*scantime = stime;
	if (loadtime != NULL)
		*loadtime = ltime;
	if (revcmpltime != NULL)
		*revcmpltime = rtime;

	return (0);
}

static void
usage(char *progname)
{
	char *slash;

	slash = strrchr(progname, '/');
	if (slash != NULL)
		progname = slash + 1;

	fprintf(stderr, "usage: %s [parameters] [options] "
	    "reads_file genome_file1 genome_file2...\n", progname);

	fprintf(stderr, "Parameters:\n");

	fprintf(stderr,
	    "    -s    Spaced Seed                             (default: %s)\n",
	    DEF_SPACED_SEED);

	fprintf(stderr,
	    "    -n    Seed Matches per Window                 (default: %d)\n",
	    DEF_NUM_MATCHES);

	fprintf(stderr,
	    "    -t    Seed Taboo Length                       (default: %d)\n",
	    DEF_TABOO_LEN);

	fprintf(stderr,
	    "    -w    Seed Window Length                      (default: %d)\n",
	    DEF_WINDOW_LEN);

	fprintf(stderr,
	    "    -o    Maximum Hits per Read                   (default: %d)\n",
	    DEF_NUM_OUTPUTS);

	fprintf(stderr,
	    "    -r    Maximum Read Length                     (default: %d)\n",
	    DEF_MAX_READ_LEN);

	fprintf(stderr,
	    "    -d    Kmer Std. Deviation Limit               (default: %d)\n",
	    DEF_KMER_STDDEV_LIMIT);

	fprintf(stderr,
	    "    -m    S-W Match Value                         (default: %d)\n",
	    DEF_MATCH_VALUE);

	fprintf(stderr,
	    "    -i    S-W Mismatch Value                      (default: %d)\n",
	    DEF_MISMATCH_VALUE);

	fprintf(stderr,
	    "    -g    S-W Gap Open Penalty                    (default: %d)\n",
	    DEF_GAP_OPEN);

	fprintf(stderr,
	    "    -e    S-W Gap Extend Penalty                  (default: %d)\n",
	    DEF_GAP_EXTEND);

	if (use_colours) {
		fprintf(stderr,
		    "    -x    S-W Crossover Penalty                   ("
		    "default: %d)\n", DEF_XOVER_PENALTY);
	}

	if (use_colours) {
		fprintf(stderr,
		    "    -h    S-W Full Hit Threshold                  "
		    "(default: %d)\n", DEF_SW_FULL_THRESHOLD);
		fprintf(stderr,
		    "    -v    S-W Vector Hit Threshold                "
		    "(default: %d)\n", DEF_SW_VECT_THRESHOLD);
	} else {
		fprintf(stderr,
		    "    -h    S-W Hit Threshold                       "
		    "(default: %d)\n", DEF_SW_FULL_THRESHOLD);
	}

	fprintf(stderr, "\n");
	fprintf(stderr, "Options:\n");

	fprintf(stderr,
	    "    -B    Print Scan Progress Bar                 (default: "
	    "disabled)\n"); 

	fprintf(stderr,
	    "    -P    Pretty Print Alignments                 (default: "
	    "disabled)\n"); 

	fprintf(stderr,
	    "    -R    Print Reads in Output                   (default: "
	    "disabled)\n");

	exit(1);
}

int
main(int argc, char **argv)
{
	struct read_elem *re;
	char *reads_file, *progname, *optstr;
	char **genome_files;
	double cellspersec, stime, ltime, rtime;
	uint64_t invocs, cells, reads_matched, total_matches;
	int ch, ret, ngenome_files;

	fprintf(stderr, "--------------------------------------------------"
	    "------------------------------\n");
	fprintf(stderr, "rmapper: %s SPACE. (SHRiMP version %s)\n",
	    (use_colours) ? "COLOUR" : "LETTER", SHRIMP_VERSION_STRING);
	fprintf(stderr, "--------------------------------------------------"
	    "------------------------------\n");

	progname = argv[0];

	if (use_colours)
		optstr = "s:n:t:w:o:r:d:m:i:g:e:x:h:v:BPR";
	else
		optstr = "s:n:t:w:o:r:d:m:i:g:e:h:BPR";

	while ((ch = getopt(argc, argv, optstr)) != -1) {
		switch (ch) {
		case 's':
			spaced_seed = xstrdup(optarg);
			break;
		case 'n':
			num_matches = atoi(optarg);
			break;
		case 't':
			taboo_len = atoi(optarg);
			break;
		case 'w':
			window_len = atoi(optarg);
			break;
		case 'o':
			num_outputs = atoi(optarg);
			break;
		case 'r':
			max_read_len = atoi(optarg);
			break;
		case 'd':
			kmer_stddev_limit = atoi(optarg);
			break;
		case 'm':
			match_value = atoi(optarg);
			break;
		case 'i':
			mismatch_value = atoi(optarg);
			break;
		case 'g':
			gap_open = atoi(optarg);
			break;
		case 'e':
			gap_extend = atoi(optarg);
			break;
		case 'x':
			assert(use_colours);
			xover_penalty = atoi(optarg);
			break;
		case 'h':
			sw_full_threshold = atoi(optarg);
			break;
		case 'v':
			assert(use_colours);
			sw_vect_threshold = atoi(optarg);
			break;
		case 'B':
			Bflag = true;
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

	if (!use_colours)
		sw_vect_threshold = sw_full_threshold;

	if (!valid_spaced_seed()) {
		fprintf(stderr, "error: invalid spaced seed\n");
		exit(1);
	}

	if (window_len < 1 || window_len < strlen(spaced_seed)) {
		fprintf(stderr, "error: invalid window length\n");
		exit(1);
	}

	if (num_matches < 1) {
		fprintf(stderr, "error: invalid number of matches\n");
		exit(1);
	}

	if (taboo_len > window_len) {
		fprintf(stderr, "error: invalid taboo length\n");
		exit(1);
	}

	if (num_outputs < 1) {
		fprintf(stderr, "error: invalid maximum hits per read\n");
		exit(1);
	}

	if (gap_open > 0) {
		fprintf(stderr, "error: invalid gap open penalty\n");
		exit(1);
	}

	if (gap_extend > 0) {
		fprintf(stderr, "error: invalid gap extend penalty\n");
		exit(1);
	}

	fprintf(stderr, "Settings:\n");
	fprintf(stderr, "    Spaced Seed:                %s\n", spaced_seed);
	fprintf(stderr, "    Spaced Seed Span:           %u\n",
	    (u_int)strlen(spaced_seed));
	fprintf(stderr, "    Spaced Seed Weight:         %u\n",
	    strchrcnt(spaced_seed, '1'));
	fprintf(stderr, "    Seed Matches per Window:    %u\n", num_matches);
	fprintf(stderr, "    Seed Taboo Length:          %u\n", taboo_len);
	fprintf(stderr, "    Seed Window Length:         %u\n", window_len);
	fprintf(stderr, "    Maximum Hits per Read:      %u\n", num_outputs);
	fprintf(stderr, "    Maximum Read Length:        %u\n", max_read_len);
	fprintf(stderr, "    Kmer Std. Deviation Limit:  %d%s\n",
	    kmer_stddev_limit, (kmer_stddev_limit < 0) ? " (None)" : "");
	fprintf(stderr, "    S-W Match Value:            %d\n", match_value);
	fprintf(stderr, "    S-W Mismatch Value:         %d\n", mismatch_value);
	fprintf(stderr, "    S-W Gap Open Penalty:       %d\n", gap_open);
	fprintf(stderr, "    S-W Gap Extend Penalty:     %d\n", gap_extend);

	if (use_colours) {
		fprintf(stderr, "    S-W Crossover Penalty:      %d\n",
		    xover_penalty);
	}

	if (use_colours) {
		fprintf(stderr, "    S-W Vector Hit Threshold:   %u\n",
		    sw_vect_threshold);
		fprintf(stderr, "    S-W Full Hit Threshold:     %u\n",
		    sw_full_threshold);
	} else {
		fprintf(stderr, "    S-W Hit Threshold:          %u\n",
		    sw_full_threshold);
	}
	fprintf(stderr, "\n");

	if (load_reads(reads_file))
		exit(1);

	if (sw_vector_setup(window_len * 2, longest_read_len, gap_open,
	    gap_extend, match_value, mismatch_value, use_colours, false)) {
		fprintf(stderr, "failed to initialise vector "
		    "Smith-Waterman (%s)\n", strerror(errno));
		exit(1);
	}

	if (use_colours) {
		ret = sw_full_cs_setup(window_len * 2, longest_read_len,
		    gap_open, gap_extend, match_value, mismatch_value,
		    xover_penalty, false);
	} else {
		ret = sw_full_ls_setup(window_len * 2, longest_read_len,
		    gap_open, gap_extend, match_value, mismatch_value, false);
	}
	if (ret) {
		fprintf(stderr, "failed to initialise scalar "
		    "Smith-Waterman (%s)\n", strerror(errno));
		exit(1);
	}

	readmap_prune();

	if (scan_genomes(genome_files, ngenome_files, &stime, &ltime, &rtime))
		exit(1);

	fprintf(stderr, "Generating output...\n");
	final_pass(genome_files, ngenome_files);

	reads_matched = total_matches = 0;
	for (re = reads; re != NULL; re = re->next) {
		reads_matched += (re->final_matches == 0) ? 0 : 1;
		total_matches += re->final_matches;
	}

	sw_vector_stats(&invocs, &cells, NULL, &cellspersec);
	
	fprintf(stderr, "\nStatistics:\n");
	fprintf(stderr, "    Spaced Seed Scan:\n");
	fprintf(stderr, "        Run-time:               %.2f seconds\n",
	    stime - ((cellspersec == 0) ? 0 : (double)cells / cellspersec));
	fprintf(stderr, "        Total Kmers:            %u\n", nkmers);
	fprintf(stderr, "        Minimal Reads/Kmer:     %" PRIu64 "\n",
	    shortest_scanned_kmer_list);
	fprintf(stderr, "        Maximal Reads/Kmer:     %" PRIu64 "\n",
	    longest_scanned_kmer_list);
	fprintf(stderr, "        Average Reads/Kmer:     %.2f\n",
	    (kmer_lists_scanned == 0) ? 0 :
	    (double)kmer_list_entries_scanned / (double)kmer_lists_scanned);
	fprintf(stderr, "\n");

	fprintf(stderr, "    Vector Smith-Waterman:\n");
	fprintf(stderr, "        Run-time:               %.2f seconds\n",
	    (cellspersec == 0) ? 0 : cells / cellspersec);
	fprintf(stderr, "        Invocations:            %" PRIu64 "\n",invocs);
	fprintf(stderr, "        Cells Computed:         %.2f million\n",
	    (double)cells / 1.0e6);
	fprintf(stderr, "        Cells per Second:       %.2f million\n",
	    cellspersec / 1.0e6);
	fprintf(stderr, "\n");

	if (use_colours)
		sw_full_cs_stats(&invocs, &cells, NULL, &cellspersec);
	else
		sw_full_ls_stats(&invocs, &cells, NULL, &cellspersec);

	fprintf(stderr, "    Scalar Smith-Waterman:\n");
	fprintf(stderr, "        Run-time:               %.2f seconds\n",
	    (cellspersec == 0) ? 0 : cells / cellspersec);
	fprintf(stderr, "        Invocations:            %" PRIu64 "\n",invocs);
	fprintf(stderr, "        Cells Computed:         %.2f million\n",
	    (double)cells / 1.0e6);
	fprintf(stderr, "        Cells per Second:       %.2f million\n",
	    cellspersec / 1.0e6);
	fprintf(stderr, "\n");

	fprintf(stderr, "    General:\n");
	fprintf(stderr, "        Reads Matched:          %" PRIu64 "    "
	    "(%.4f%%)\n", reads_matched,
	    (nreads == 0) ? 0 : ((double)reads_matched / (double)nreads) * 100);
	fprintf(stderr, "        Total Matches:          %" PRIu64 "\n",
	    total_matches);
	fprintf(stderr, "        Avg Hits/Matched Read:  %.2f\n",
	    (total_matches == 0) ? 0 : ((double)total_matches /
	    (double)reads_matched));
	
	return (0);
}
