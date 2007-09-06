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
#include "../common/util.h"

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
static u_int sw_threshold	= DEF_SW_THRESHOLD;

/* Genomic sequence, stored in 32-bit words, first is in the LSB */
static uint32_t *genome;
static uint32_t	 genome_len;

/*
 * Genomic sequence, in letter space, stored in 32-bit words, first is in LSB
 *
 * NB: Colour space only
 */
static uint32_t *genome_ls;
static uint32_t  genome_ls_len;

/* Reads */
static struct read_elem  *reads;		/* list of reads */
static struct read_node **readmap;		/* kmer to read map */
static uint32_t		 *readmap_len;		/* entries in each readmap[n] */

static u_int	nreads;				/* total reads loaded */
static u_int	nkmers;				/* total kmers of reads loaded*/

static u_int	longest_read_len;		/* longest read we've seen */

/* Flags */
static int bflag = 0;				/* print a progress bar */
static int pflag = 0;				/* pretty print results */

/* Scan stats */
static uint64_t shortest_scanned_kmer_list;
static uint64_t longest_scanned_kmer_list;
static uint64_t kmer_lists_scanned;
static uint64_t kmer_list_entries_scanned;

#ifdef USE_COLOURS
const int use_colours = 1;
#else
const int use_colours = 0;
#endif

struct alignment_sort {
	int32_t		score;
	uint32_t	index;
	char	       *name;
	char	       *dbalign;
	char	       *qralign;
	void	       *cookie;
	int		goff;
};

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

/* qsort callback - descending score, increasing index */
static int
qsort_alignments(const void *a, const void *b)
{
	const struct alignment_sort *sa, *sb;

	sa = a;
	sb = b;

	if (sa->score < sb->score)
		return (1);

	if (sa->score > sb->score)
		return (-1);

	if (sa->index > sb->index)
		return (1);

	if (sa->index < sb->index)
		return (-1);

	return (0);
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
save_score(struct read_elem *re, int score, int index)
{
	struct re_score *scores;

	scores = re->scores;

	if (score < scores[1].score)
		return;

	scores[1].score = score;
	scores[1].index = index;
	reheap(scores, 1);
}

static void
progress_bar(uint32_t at, uint32_t of)
{
	static int lastperc, beenhere;
	static char whirly = '\\';

	char progbuf[52];
	int perc, i, j, dec;

	perc = (at * 1000) / of;

	if (beenhere && perc - lastperc != 1)
		return;

	beenhere = 1;
	lastperc = perc;

	dec = perc % 10;
	perc /= 10;

	/* any excuse to have a whirly gig */
	switch (whirly) {
	case '|':
		whirly = '/';
		break;
	case '/':
		whirly = '-';
		break;
	case '-':
		whirly = '\\';
		break;
	case '\\':
		whirly = '|';
		break;
	}
	progbuf[25] = whirly;
		
	for (i = j = 0; i <= 100; i += 2) {
		if (j != 25) {
			if (i <= perc)
				progbuf[j++] = '=';
			else
				progbuf[j++] = ' ';
		} else {
			j++;
		}
	}
	progbuf[51] = '\0';

	fprintf(stderr, "\rProgress: [%s] %3d.%1d%%", progbuf, perc, dec);
	fflush(stderr);
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
 * Prepend the low 4 bits of 'val' to the start of the bitfield in 'bf'.
 * 'entries' is the maximum number of 4-bit words to be stored in the
 * bitfield.
 */
static void
bitfield_prepend(uint32_t *bf, int entries, uint32_t val)
{
	uint32_t tmp;
	int i;

	for (i = 0; i < BPTO32BW(entries); i++) {
		tmp = bf[i] >> 28;
		bf[i] <<= 4;
		bf[i] |= val;
		val = tmp;
	}

	bf[i - 1] &= (0xffffffff >> (32 - (4 * (entries % 8))));
}

/*
 * Append the low 4 bits of 'val' to the end of the bitfield in 'bf'.
 * 'entries' is the number of 4-bit words in 'bf' prior to the append.
 */
static void
bitfield_append(uint32_t *bf, uint32_t entries, uint32_t val)
{
	uint32_t word;

	word = bf[entries / 8];
	word &= ~(0xf << (4 * (entries % 8)));
	word |= ((val & 0xf) << (4 * (entries % 8)));
	bf[entries / 8] = word;
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

	for (i = 0; i < j; i++) {
		if (readmap_len[i] > mean + (kmer_stddev_limit * stddev))
			readmap[i] = NULL;
	}
}

/* scan the genome by kmers, running S-W as needed, and updating scores */
static void
scan()
{
	struct read_node *rn;
	struct read_elem *re;
	uint32_t *kmer;
	uint32_t i, j, idx, mapidx;
	uint32_t base, skip, score;
	uint32_t goff, glen;
	int prevhit, seed_span;

	seed_span = strlen(spaced_seed);

	kmer = xmalloc(sizeof(kmer[0]) * BPTO32BW(seed_span));
	memset(kmer, 0, sizeof(kmer[0]) * BPTO32BW(seed_span));
	for (i = 0; i < seed_span - 1; i++)
		bitfield_prepend(kmer, seed_span, EXTRACT(genome, i));

	skip = 0;
	for (i = seed_span - 1; i < genome_len; i++) {
		if (bflag)
			progress_bar(i, genome_len);

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

				if (score >= sw_threshold) {
					save_score(re, score, i);
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

	if (bflag) {
		progress_bar(i, genome_len);
		fprintf(stderr, "\n\n");
	}
}

static void
load_genome_helper(int base, ssize_t offset, int isnewentry, char *name,
    int initbp)
{
	static int first = 1;
	static int last;

	/* shut up icc */
	(void)name;
	(void)last;
	(void)initbp;

	/* handle initial call to alloc resources */
	if (base == -1) {
		genome = xmalloc(sizeof(genome[0]) * BPTO32BW(offset));
		memset(genome, 0, sizeof(genome[0]) * BPTO32BW(offset));

		if (use_colours) {
			genome_ls =
			    xmalloc(sizeof(genome_ls[0]) * BPTO32BW(offset));
			memset(genome_ls, 0,
			    sizeof(genome_ls[0]) * BPTO32BW(offset));
			last = BASE_T;
		}
		return;
	}

	if (isnewentry && !first) {
		fprintf(stderr, "error: genome file consists of more than one "
		    "reftig!\n");
		exit(1);
	}

	assert(base >= 0 && base <= 5);
	assert(offset == genome_len);

	if (use_colours) {
		bitfield_append(genome, offset, lstocs(last, base));
		genome_len++;
		bitfield_append(genome_ls, offset, base);
		genome_ls_len++;
		last = base;
	} else {
		bitfield_append(genome, offset, base);
		genome_len++;
	}

	if (first)
		printf("# REFTIG: [%s]\n", name);

	first = 0;
}

static int
load_genome(const char *file)
{
	ssize_t ret;

	ret = load_fasta(file, load_genome_helper, LETTER_SPACE);

	if (ret != -1)
		fprintf(stderr, "Loaded %u letters of genome\n",
		    (unsigned int)ret);

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
	if (base == -1) {
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
		re->scores[0].score = num_outputs;
		re->scores[0].index = 0;

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
		int unique = 1;

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
				unique = 0;
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
print_pretty(const char *name, const struct sw_full_results *sfr,
    const char *dbalign, const char *qralign, uint32_t goff, int newread)
{
	static int firstcall = 1;

	int j, len;
	uint32_t *genome_ptr;
	uint32_t aoff;
	char translate[5] = { 'A', 'C', 'G', 'T', 'N' };

	if (use_colours)
		genome_ptr = genome_ls;
	else
		genome_ptr = genome;

	if (newread && !firstcall) {
		putchar('\n');
		putchar('\n');
	}
	firstcall = 0;

	if (newread) {
		printf("---------------------------------------------------\n");
		printf("READ: [%s]\n", name);
		printf("---------------------------------------------------\n");
	}

	putchar('\n');

	aoff = sfr->genome_start;

	len = strlen(dbalign);
	assert(len == strlen(qralign));

	printf("Score:   %d   (read start/end: %d/%d)\n",
	    sfr->score, sfr->read_start,
	    sfr->read_start + sfr->mapped - 1);
	printf("Index:   %u\n", goff + aoff);
	
	printf("Reftig:  ");
	for (j = 4; j > 0; j--) {
		if (j <= goff + aoff)
			putchar(translate[
			    EXTRACT(genome_ptr, goff + aoff - j)]);
	}
	printf("%s", dbalign);
	for (j = 0; j < 4; j++) {
		if ((goff + aoff + len + j) < genome_len)
			putchar(translate[EXTRACT(genome_ptr,
			    goff + aoff + len + j)]);
	}
	putchar('\n');

	printf("Match:   ");
	for (j = 4; j > 0; j--) {
		if (j <= goff + aoff)
			putchar(' ');
	}
	for (j = 0; j < len; j++) {
		if (dbalign[j] == qralign[j] && dbalign[j] != '-') {
			putchar('|');
		} else {
			assert(j != 0);
			if (dbalign[j] == toupper((int)qralign[j]))
				putchar('X');
			else if (islower((int)qralign[j]))
				putchar('x');
			else
				putchar(' ');
		}
	}
	putchar('\n');

	printf("Read:    ");
	for (j = 4; j > 0; j--) {
		if (j <= goff + aoff)
			putchar(' ');
	}
	printf("%s\n", qralign);
}

static void
print_normal(const char *name, const struct sw_full_results *sfr,
    uint32_t goff, int newread)
{
	static int firstcall = 1;

	if (!firstcall && newread)
		putchar('\n');
	firstcall = 0;

	if (use_colours) {
		printf("[%s] %d %u %d %d %d %d %d %d %d\n", name, sfr->score,
		    goff + sfr->genome_start, sfr->read_start, sfr->mapped,
		    sfr->matches, sfr->mismatches, sfr->insertions,
		    sfr->deletions, sfr->crossovers);
	} else {
		printf("[%s] %d %u %d %d %d %d %d %d\n", name, sfr->score,
		    goff + sfr->genome_start, sfr->read_start, sfr->mapped,
		    sfr->matches, sfr->mismatches, sfr->insertions,
		    sfr->deletions);
	}
}

static void
print_alignments(struct read_elem *re)
{
	struct sw_full_results *sfr;
	struct alignment_sort *qa;
	uint32_t goff, glen;
	uint32_t prev_index;
	int32_t prev_score;
	int i, newread;

	/*
	 * Allocate space and do all of the full S-W scans up front, then
	 * sort, and finally output. This is not necessary in the letter
	 * space case, but it is needed for the colour space and common
	 * code is always nicer.
	 */
	sfr = xmalloc(sizeof(*sfr) * re->scores[0].score);
	qa = xmalloc(sizeof(*qa) * re->scores[0].score);

	memset(sfr, 0, sizeof(*sfr) * re->scores[0].score);
	memset(qa, 0, sizeof(*qa) * re->scores[0].score);

	for (i = 1; i <= re->scores[0].score; i++) {
		if (re->scores[i].score < 0)
			continue;

		if (re->scores[i].index < window_len)
			goff = 0;
		else
			goff = re->scores[i].index - window_len;

		if (goff + (window_len * 2) >= genome_len)
			glen = genome_len - goff;
		else
			glen = window_len * 2;

		if (use_colours) {
			sw_full_cs(genome, goff, glen, re->read, re->read_len,
			    genome_ls, re->initbp, re->scores[i].score,
			    &qa[i-1].dbalign, &qa[i-1].qralign, &sfr[i-1]);
		} else {
			sw_full_ls(genome, goff, glen, re->read, re->read_len,
			    re->scores[i].score, &qa[i-1].dbalign,
			    &qa[i-1].qralign, &sfr[i-1]);
		}

		qa[i-1].score = sfr[i-1].score;
		qa[i-1].index = sfr[i-1].genome_start + goff;
		qa[i-1].name = re->name;
		qa[i-1].dbalign = xstrdup(qa[i-1].dbalign);
		qa[i-1].qralign = xstrdup(qa[i-1].qralign);
		qa[i-1].cookie = &sfr[i-1];
		qa[i-1].goff = goff;
	}

	qsort(qa, re->scores[0].score, sizeof(qa[0]), qsort_alignments);

	newread = 1;

	/* shut up, gcc */
	prev_index = prev_score = -1;

	for (i = 0; i < re->scores[0].score; i++) {
		if (qa[i].cookie == NULL)
			continue;
	
		/* skip duplicates */
		if (newread == 0 && qa[i].score == prev_score &&
		    qa[i].index == prev_index)
			continue;

		if (pflag) {
			print_pretty(qa[i].name, qa[i].cookie, qa[i].dbalign,
			    qa[i].qralign, qa[i].goff, newread);
		} else {
			print_normal(qa[i].name, qa[i].cookie, qa[i].goff,
			    newread);
		}

		newread = 0;
		
		free(qa[i].dbalign);
		free(qa[i].qralign);

		prev_score = qa[i].score;
		prev_index = qa[i].index;
	}

	free(sfr);
	free(qa);
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
usage(char *progname)
{
	char *slash;

	slash = strrchr(progname, '/');
	if (slash != NULL)
		progname = slash + 1;

	fprintf(stderr, "usage: %s [parameters] [options] "
	    "genome_file reads_file\n", progname);

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

	fprintf(stderr,
	    "    -h    S-W Hit Threshold                       (default: %d)\n",
	    DEF_SW_THRESHOLD);

	fprintf(stderr, "\n");
	fprintf(stderr, "Options:\n");

	fprintf(stderr,
	    "    -b    Print Scan Progress Bar                 (default: "
	    "disabled)\n"); 

	fprintf(stderr,
	    "    -p    Pretty Print Alignments                 (default: "
	    "disabled)\n"); 

	exit(1);
}

int
main(int argc, char **argv)
{
	struct timeval tv1, tv2;
	char *genome_file, *reads_file, *progname, *optstr;
	struct read_elem *re;
	double cellspersec;
	uint64_t invocs, cells;
	int ch, ret, hits = 0;

	fprintf(stderr, "------------------------------------------------\n");
	fprintf(stderr, "rmapper: %s SPACE. Version %s\n",
	    (use_colours) ? "COLOUR" : "LETTER", RMAPPER_VERSION_STR);
	fprintf(stderr, "------------------------------------------------\n\n");

	progname = argv[0];

	if (use_colours)
		optstr = "s:n:t:w:o:r:d:m:i:g:e:x:h:bp";
	else
		optstr = "s:n:t:w:o:r:d:m:i:g:e:h:bp";

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
			sw_threshold = atoi(optarg);
			break;
		case 'b':
			bflag = 1;
			break;
		case 'p':
			pflag = 1;
			break;
		default:
			usage(progname);
		}
	}
	argc -= optind;
	argv += optind;
	
	if (argc != 2) {
		fprintf(stderr, "error: %sreads_file not specified\n",
		    (argc == 0) ? "genome_file, " : "");
		usage(progname);
	}

	genome_file = argv[0];
	reads_file  = argv[1];

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

	if (load_genome(genome_file)) {
		exit(1);
	}

	if (load_reads(reads_file)) {
		exit(1);
	}

	if (sw_vector_setup(window_len * 2, longest_read_len,
	    gap_open, gap_extend, match_value, mismatch_value, use_colours)) {
		fprintf(stderr, "failed to initialise vector "
		    "Smith-Waterman (%s)\n", strerror(errno));
		exit(1);
	}

	if (use_colours) {
		ret = sw_full_cs_setup(window_len * 2, longest_read_len,
		    gap_open, gap_extend, match_value, mismatch_value,
		    xover_penalty);
	} else {
		ret = sw_full_ls_setup(window_len * 2, longest_read_len,
		    gap_open, gap_extend, match_value, mismatch_value);
	}
	if (ret) {
		fprintf(stderr, "failed to initialise scalar Smith-Waterman "
		    "(%s)\n", strerror(errno));
		exit(1);
	}

	fprintf(stderr, "\n");
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

	fprintf(stderr, "    S-W Hit Threshold:          %u\n", sw_threshold);
	fprintf(stderr, "\n");

	readmap_prune();

	gettimeofday(&tv1, NULL);
	scan();
	gettimeofday(&tv2, NULL);

	if (!pflag) {
		if (use_colours) {
			printf("#\n# [READ_NAME] score index read_start "
			    "read_mapped_length matches mismatches "
			    "insertions deletions crossovers\n");
		} else {
			printf("#\n# [READ_NAME] score index read_start "
			    "read_mapped_length matches mismatches "
			    "insertions deletions\n");
		}
	}

	for (re = reads; re != NULL; re = re->next) {
		if (re->swhits == 0)
			continue;

		print_alignments(re);
		hits++;
	}

	sw_vector_stats(&invocs, &cells, NULL, &cellspersec);
	
	fprintf(stderr, "Statistics:\n");
	fprintf(stderr, "    Spaced Seed Scan:\n");
	fprintf(stderr, "        Run-time:               %.2f seconds\n",
	    ((double)tv2.tv_sec + (double)tv2.tv_usec / 1.0e6) -
	    ((double)tv1.tv_sec + (double)tv1.tv_usec / 1.0e6) -
	    (double)cells / cellspersec);
	fprintf(stderr, "        Total Kmers:            %u\n", nkmers);
	fprintf(stderr, "        Minimal Reads/Kmer:     %" PRIu64 "\n",
	    shortest_scanned_kmer_list);
	fprintf(stderr, "        Maximal Reads/Kmer:     %" PRIu64 "\n",
	    longest_scanned_kmer_list);
	fprintf(stderr, "        Average Reads/Kmer:     %.2f\n",
	    (double)kmer_list_entries_scanned / (double)kmer_lists_scanned);
	fprintf(stderr, "\n");

	fprintf(stderr, "    Vector Smith-Waterman:\n");
	fprintf(stderr, "        Run-time:               %.2f seconds\n",
	    cells / cellspersec);
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
	    cells / cellspersec);
	fprintf(stderr, "        Invocations:            %" PRIu64 "\n",invocs);
	fprintf(stderr, "        Cells Computed:         %.2f million\n",
	    (double)cells / 1.0e6);
	fprintf(stderr, "        Cells per Second:       %.2f million\n",
	    cellspersec / 1.0e6);
	fprintf(stderr, "\n");

	fprintf(stderr, "    General:\n");
	fprintf(stderr, "        Reads Matched:          %d\n", hits);
	
	return (0);
}
