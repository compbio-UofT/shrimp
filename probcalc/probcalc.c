/*	$Id$	*/

#include <dirent.h>
#include <errno.h>
#include <fcntl.h>
#include <inttypes.h>
#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <zlib.h>

#include <sys/types.h>
#include <sys/stat.h>

#include "../common/sw-full-common.h"
#include "../common/input.h"
#include "../common/dynhash.h"
#include "../common/util.h"
#include "../common/version.h"

#include "probcalc.h"

static dynhash_t read_list;		/* list of reads and top matches */
static dynhash_t contig_cache;		/* contig name cache */
static dynhash_t read_seq_cache;	/* read sequence cache */
static dynhash_t read_edit_cache;	/* read edit string cache */

static bool Bflag = false;		/* print a progress bar */
static bool Gflag = false;		/* calculate rates and output them */
static bool Rflag = false;		/* include read sequence in output */
static bool Sflag = false;		/* do it all in one pass (save in ram)*/

static char *rates_file;		/* -g flag specifies a rates file */
static char *rates_string;		/* -r user-supplied rates */

struct readinfo {
	char	    *name;
	struct input top_matches[0];
};

struct rates {
	uint64_t samples;
	uint64_t total_len;
	uint64_t insertions;
	uint64_t deletions;
	uint64_t matches;
	uint64_t mismatches;
	uint64_t crossovers;

	double   erate;			/* error rate (crossovers) */
	double   srate;			/* substitution rate (mismatches) */
	double   irate;			/* indel rate (insertions, deletions) */
	double   mrate;			/* match rate (matches) */
};

struct readstatspval {
	struct input     *rs;
	double		  pchance;
	double		  pgenome;
	double            normodds;
};

struct pass_cb {
	struct rates	rates;
	uint64_t	nbytes;
	uint64_t	nfiles;
	uint64_t	total_files;
	int		pass;
};

enum {
	SORT_PCHANCE = 1,
	SORT_PGENOME,
	SORT_NORMODDS
};

static uint64_t total_alignments;	/* total outputs in all files */
static uint64_t total_unique_reads;	/* total unique reads in all files */

/* arguments */
static double   normodds_cutoff		= DEF_NORMODDS_CUTOFF;
static double   pchance_cutoff		= DEF_PCHANCE_CUTOFF;
static double	pgenome_cutoff	        = DEF_PGENOME_CUTOFF;
static int	top_matches		= DEF_TOP_MATCHES;
static uint64_t genome_len;
static int	sort_field		= SORT_PCHANCE;

#define PROGRESS_BAR(_a, _b, _c, _d)	\
    if (Bflag) progress_bar((_a), (_b), (_c), (_d))

#define ALMOST_ZERO	0.000000001
#define ALMOST_ONE	0.999999999

static uint32_t
keyhasher(void *k)
{
	
	return (hash_string((char *)k));
} 

static int
keycomparer(void *a, void *b)
{

	return (strcmp((char *)a, (char *)b));
}

/* percolate down in our min-heap */
static void
reheap(struct input *stats, int node)
{
	struct input tmp;
	int left, right, max;

	left  = node * 2;
	right = left + 1;
	max   = node;
	
	if (left <= stats[0].score &&
	    stats[left].score < stats[node].score)
		max = left;

	if (right <= stats[0].score &&
	    stats[right].score < stats[max].score)
		max = right;

	if (max != node) {
		tmp = stats[node];
		stats[node] = stats[max];
		stats[max] = tmp;
		reheap(stats, max);
	}
}

/* update our match min-heap */
static void
save_match(struct readinfo *ri, struct input *input)
{
	struct input *stats = ri->top_matches;

	if (input->score < stats[1].score)
		return;

	memcpy(&stats[1], input, sizeof(stats[1]));
	reheap(stats, 1);
}

/* Returns the number of bytes read */
static void
read_file(const char *fpath)
{
	struct input input;
	struct readinfo *ri;
	char *key;
	gzFile fp;
	int i;

	fp = gzopen(fpath, "r");
	if (fp == NULL) {
		fprintf(stderr, "error: could not open file [%s]: %s\n",
		    fpath, strerror(errno));
		exit(1);
	}

	while (input_parseline(fp, &input)) {
		/*
		 * If we don't want the read sequence, free it if it exists.
		 * Otherwise, cache it so we don't have duplicates hogging
		 * memory.
		 */
		if (!Rflag) {
			if (input.read_seq != NULL)
				free(input.read_seq);
			input.read_seq = NULL;
		} else if (input.read_seq != NULL) {
			if (dynhash_find(read_seq_cache, input.read_seq,
			    (void **)&key, NULL)) {
				free(input.read_seq);
				input.read_seq = key;
			} else {
				if (dynhash_add(read_seq_cache, input.read_seq,
				    NULL) == false) {
					fprintf(stderr, "error: failed to add "
					    "to read_seq cache (%d) - probably "
					    "out of memory\n", __LINE__);
					exit(1);
				}
			}
		}

		/*
		 * Cache the contig names similarly, so we don't waste memory.
		 */
		if (dynhash_find(contig_cache, input.genome, (void **)&key,
		    NULL)) {
			free(input.genome);
			input.genome = key;
		} else {
			if (dynhash_add(contig_cache, input.genome,
			    NULL) == false) {
				fprintf(stderr, "error: failed to add to "
				    "contig cache (%d) - probably out of "
				    "memory\n", __LINE__);
				exit(1);
			}
		}

		/*
		 * Finally, cache edit strings for memory happiness.
		 */
		if (dynhash_find(read_edit_cache, input.edit, (void **)&key,
		    NULL)) {
			free(input.edit);
			input.edit = key;
		} else {
			if (dynhash_add(read_edit_cache, input.edit,
			    NULL) == false) {
				fprintf(stderr, "error: failed to add to "
				    "edit cache (%d) - probably out of "
				    "memory\n", __LINE__);
				exit(1);
			}
		}

		total_alignments++;

		if (dynhash_find(read_list, input.read, (void **)&key,
		    (void **)&ri)) {
			/* use only one read name string to save space */
			free(input.read);
			input.read = key;

			/* test if score better than worst, if so, add it */
			if (input.score > ri->top_matches[1].score)
				save_match(ri, &input);
		} else {
			total_unique_reads++;

			/* alloc, init, insert in dynhash */
			ri = (struct readinfo *)xmalloc(sizeof(*ri) +
			    sizeof(ri->top_matches[0]) * (top_matches + 1));
			memset(ri, 0, sizeof(*ri) +
			    sizeof(ri->top_matches[0]) * (top_matches + 1));

			ri->name = input.read;
			ri->top_matches[0].score = top_matches;
			for (i = 1; i <= top_matches; i++)
				ri->top_matches[i].score = 0x80000000 + i;

			save_match(ri, &input);

			if (dynhash_add(read_list, ri->name, ri) == false) {
				fprintf(stderr, "error: failed to add to "
				    "dynhash - probably out of memory\n");
				exit(1);
			}
		}
	}

	gzclose(fp);
}

static double
p_chance(uint64_t l, int k, int nsubs, int nerrors, int nindels, int origlen)
{
	double r = 0;

	if (isinf(pow(4, (nsubs + nerrors))))
		r += ls_choose(k - 2, nsubs + nerrors) + ((nsubs + nerrors) * log(4));
	else
		r += ls_choose(k - 2, nsubs + nerrors) + log(pow(4, (nsubs + nerrors)) - 1);

	/*
	 * two for ins and dels, four for the four possible letters to insert.
	 * above should really take length of indel into account
	 */
	r += ls_choose(k - 1, nindels) + (nindels * log(2*4));

	/* r is now log of fraction of matched words over total words */
	r -= (k * log(4));
	r = exp(r);
	r = 1 - r;

	/* if hit not full length we need to correct */
	l *= MAX(origlen - k+1, 1);

	/* 2 for the two strands */
	r = 2.0 * l * log(r);

	return (1.0 - exp(r));
}

static void
calc_rates(void *arg, void *key, void *val)
{
	struct rates *rates = (struct rates *)arg;
	struct readinfo *ri = (struct readinfo *)val;
	struct input *rs;
	double d;
	int32_t best;
	int i, rlen;

	(void)key;

	if (Sflag) {
		static uint64_t calls;
		PROGRESS_BAR(stderr, calls, total_unique_reads, 100);
		calls++;
	}

	best = 0;
	for (i = 1; i <= ri->top_matches[0].score; i++) {
		if (best == 0 ||
		    ri->top_matches[i].score > ri->top_matches[best].score) {
			/* NB: strictly >, as we want same results for double and single passes
			       if there are two with the same score, but different alignments */
			best = i;
		}
	}

	if (!Sflag) {
		assert(ri->top_matches[0].score == 1);
		assert(ri->top_matches[1].score >= 0);
	}

	/*
	 * Only count the best score with no equivalent best scores.
	 */
	assert(best != 0);
	rs = &ri->top_matches[best];
	rlen =  rs->matches + rs->mismatches + rs->insertions + rs->deletions;

	d = p_chance(genome_len, rlen, rs->mismatches, rs->crossovers,
		rs->insertions + rs->deletions, rs->read_length);

	if (best != 0 && d < pchance_cutoff) {
		rates->samples++;
		rates->total_len  += rs->matches + rs->mismatches;
		rates->insertions += rs->insertions;
		rates->deletions  += rs->deletions;
		rates->matches    += rs->matches;
		rates->mismatches += rs->mismatches;
		rates->crossovers += rs->crossovers;
	}
}

static double
p_thissource(int k, int nerrors, double erate, int nsubs, double subrate,
    int nindels, double indelrate, int nmatches, double matchrate, int origlen)
{
	double r;

	(void)origlen;

	r  = ls_choose(k, nerrors) + (nerrors * log(erate));
	r += ls_choose(k - nerrors, nsubs) + (nsubs * log(subrate));
	r += ls_choose(k - nerrors - nsubs, nindels) +
	    (nindels * log(indelrate));
	r += nmatches * log(matchrate);
	r = exp(r);


	if (r < ALMOST_ZERO)
		r = ALMOST_ZERO;
	if (r > ALMOST_ONE)
		r = ALMOST_ONE;

	return (r);
}

static int
rspvcmp(const void *a, const void *b)
{
	const struct readstatspval *rspva = (const struct readstatspval *)a;
	const struct readstatspval *rspvb = (const struct readstatspval *)b;
	double cmp_a, cmp_b;

	switch (sort_field) {
	case SORT_PGENOME:
		/* reversed - want big first */
	        cmp_a = rspvb->pgenome;
		cmp_b = rspva->pgenome;
		break;
	case SORT_PCHANCE:
		cmp_a = rspva->pchance;
		cmp_b = rspvb->pchance;
		break;
	case SORT_NORMODDS:
		/* reversed - want big first */
	        cmp_a = rspvb->normodds;
		cmp_b = rspva->normodds;
		break;
	default:
		/* shut up, gcc */
		cmp_a = cmp_b = 0;
		assert(0);
	}

	if (cmp_a == cmp_b)
		return (0);
	else if (cmp_a < cmp_b)
		return (-1);
	else
		return (1);
}

static void
calc_probs(void *arg, void *key, void *val)
{
	static struct readstatspval *rspv;
	static bool called = false;

	struct rates *rates = (struct rates *)arg;
	struct readinfo *ri = (struct readinfo *)val;
	double s, norm;
	int i, j, rlen;

	(void)key;

	if (Sflag) {
		static uint64_t calls;
		PROGRESS_BAR(stderr, calls, total_unique_reads, 10);
		calls++;
	}

	if (rspv == NULL)
		rspv = (struct readstatspval *)xmalloc(sizeof(rspv[0]) * (top_matches + 1));

	/* 1: Calculate P(chance) and P(genome) for each hit */
	norm = 0;
	for (i = 1, j = 0; i <= top_matches; i++) {
		struct input *rs = &ri->top_matches[i];

		if (rs->score < 0)
			continue;

		rlen = rs->matches + rs->mismatches + rs->insertions +
		    rs->deletions;
		s = p_chance(genome_len, rlen, rs->mismatches, rs->crossovers,
		    rs->insertions + rs->deletions, rs->read_length);
		if (s < ALMOST_ZERO)
			s = ALMOST_ZERO;

		if (s > pchance_cutoff)
			continue;

		rspv[j].rs = rs;
		rspv[j].pchance = s;

		rlen = rs->matches + rs->mismatches;
		rspv[j].pgenome = p_thissource(rlen, rs->crossovers,
		    rates->erate, rs->mismatches, rates->srate,
		    rs->insertions + rs->deletions, rates->irate,
		    rs->matches, rates->mrate, rs->read_length);
		rspv[j].normodds = rspv[j].pgenome / rspv[j].pchance;
		norm += rspv[j].normodds;
		j++;
	}

	/* 2: Normalise our values */
	for (i = 0; i < j; i++)
		rspv[i].normodds = rspv[i].normodds / norm;

	/* 3: Sort the values ascending */
	qsort(rspv, j, sizeof(*rspv), rspvcmp);

	/* 4: Finally, print out values in ascending order until the cutoff */
	if (!called) {
		printf("#FORMAT: readname contigname strand contigstart contigend readstart readend "
		    "readlength score editstring %snormodds pgenome pchance\n",
		    (Rflag) ? "readsequence " : "");
		called = true;
	}

	for (i = 0; i < j; i++) {
		struct input *rs = rspv[i].rs;
		char *readseq = "";

		if (rspv[i].normodds < normodds_cutoff) {
			if (sort_field == SORT_NORMODDS)
				break;
			else
				continue;
		}
		if (rspv[i].pgenome < pgenome_cutoff) {
			if (sort_field == SORT_PGENOME)
				break;
			else
				continue;
		}
		if (rspv[i].pchance > pchance_cutoff) {
			if (sort_field == SORT_PCHANCE)
				break;
			else
				continue;
		}

		if (rs->read_seq != NULL && rs->read_seq[0] != '\0')
			readseq = rs->read_seq;
		else if (Rflag)
			readseq = " ";

		/* the only sane way is to reproduce the output code, sadly. */
		printf(">%s\t%s\t%c", rs->read, rs->genome,
		    (INPUT_IS_REVCMPL(rs)) ? '-' : '+');

		/* NB: internally 0 is first position, output is 1. adjust */
		printf("\t%u\t%u\t%d\t%d\t%d\t%d\t%s\t%s%s%f\t%f\t%f\n",
		    rs->genome_start + 1, rs->genome_end + 1, rs->read_start + 1,
		    rs->read_end + 1, rs->read_length, rs->score, rs->edit,
		    (Rflag) ? readseq : "", (Rflag) ? "\t" : "",
		    rspv[i].normodds, rspv[i].pgenome, rspv[i].pchance);
	}
}

static unsigned int cleanup_cb_called = 0;

static void
cleanup_cb(void *arg, void *key, void *val)
{

	/* shut up, icc */
	(void)arg;

	cleanup_cb_called++;
	free(key);
	if (val != NULL) {
		assert(arg == read_list);
		free(val);
	} else {
		assert(arg == read_seq_cache ||
		       arg == read_edit_cache);
	}
}

/* Remove everything from read_list and read_seq_cache. */
static void
cleanup()
{

	cleanup_cb_called = 0;
	dynhash_iterate(read_list, cleanup_cb, read_list);
	assert(dynhash_count(read_list) == cleanup_cb_called);
	dynhash_destroy(read_list);

	read_list = dynhash_create(keyhasher, keycomparer);
	if (read_list == NULL) {
		fprintf(stderr, "error: failed to allocate read_list\n");
		exit(1);
	}

	cleanup_cb_called = 0;
	dynhash_iterate(read_seq_cache, cleanup_cb, read_seq_cache);
	assert(dynhash_count(read_seq_cache) == cleanup_cb_called);
	dynhash_destroy(read_seq_cache);

	read_seq_cache = dynhash_create(keyhasher, keycomparer);
	if (read_seq_cache == NULL) {
		fprintf(stderr, "error: failed to allocate read_seq_cache\n");
		exit(1);
	}

	dynhash_iterate(read_edit_cache, cleanup_cb, read_edit_cache);
	dynhash_destroy(read_edit_cache);
	read_edit_cache = dynhash_create(keyhasher, keycomparer);
	if (read_edit_cache == NULL) {
		fprintf(stderr, "error: failed to allocate read_edit_cache\n");
		exit(1);
	}
}

static void
load_rates(const char *rates_file, struct rates *rates)
{
	FILE *fp;
	char buf[512];
	char *b;

	memset(rates, 0, sizeof(*rates));

	fp = fopen(rates_file, "r");
	if (fp == NULL) {
		fprintf(stderr, "error: failed to open rates file [%s]\n",
		    rates_file);
		exit(1);
	}

	while (fgets(buf, sizeof(buf), fp) != NULL) {
		int r;
		uint64_t totaligns, totunique, samples, length, insertions;
		uint64_t deletions, matches, mismatches, xovers;

		if (buf[0] != '>')
			continue;

		b = strtrim(buf + 1);

		/* tot. alignments, tot. unique reads, samples, length, insertions, deletions,
		 * matches, mismatches, xovers */
		r = sscanf(b, "%"PRIu64 "%"PRIu64 "%"PRIu64 "%"PRIu64 "%"PRIu64 "%"PRIu64 "%"PRIu64
		    "%"PRIu64 "%"PRIu64, &totaligns, &totunique, &samples, &length, &insertions,
		    &deletions, &matches, &mismatches, &xovers);

		if (r != 9) {
			fprintf(stderr, "error: failed to parse rates file line [%s]\n", b);
			exit(1);
		}

		total_alignments  += totaligns;
		total_unique_reads+= totunique;
		rates->samples    += samples;
		rates->total_len  += length;
		rates->insertions += insertions;
		rates->deletions  += deletions;
		rates->matches    += matches;
		rates->mismatches += mismatches;
		rates->crossovers += xovers;
	}

	fclose(fp);
}

static void
ratestats(struct rates *rates)
{

	fprintf(stderr, "\nCalculated rates for %" PRIu64 " reads\n",
	    total_unique_reads);

	fprintf(stderr, "\nStatistics:\n");
	fprintf(stderr, "    total matches:      %" PRIu64 "\n",
	    total_alignments);
	fprintf(stderr, "    total unique reads: %" PRIu64 "\n",
	    total_unique_reads);
	fprintf(stderr, "    total samples:      %" PRIu64 "\n",
	    rates->samples);
	fprintf(stderr, "    total length:       %" PRIu64 "\n",
	    rates->total_len); 
	fprintf(stderr, "    insertions:         %" PRIu64 "\n",
	    rates->insertions); 
	fprintf(stderr, "    deletions:          %" PRIu64 "\n",
	    rates->deletions); 
	fprintf(stderr, "    matches:            %" PRIu64 "\n",
	    rates->matches); 
	fprintf(stderr, "    mismatches:         %" PRIu64 "\n",
	    rates->mismatches); 
	fprintf(stderr, "    crossovers:         %" PRIu64 "\n",
	    rates->crossovers); 

	if (Gflag) {
		printf(">%" PRIu64 " %" PRIu64 " %"PRIu64 " %" PRIu64 " %" PRIu64
		    " %" PRIu64 " %" PRIu64 " %" PRIu64 " %" PRIu64 "\n",
		    total_alignments, total_unique_reads, rates->samples, rates->total_len,
		    rates->insertions, rates->deletions, rates->matches, rates->mismatches,
		    rates->crossovers);
		exit(0);
	}

	rates->erate = (double)rates->crossovers / (double)rates->total_len;
	rates->srate = (double)rates->mismatches / (double)rates->total_len;
	rates->irate = (double)(rates->insertions + rates->deletions) /
	    (double)rates->total_len;
	rates->mrate = (double)rates->matches / (double)rates->total_len;

	if (rates->erate == 0.0 || rates->srate == 0.0 || rates->irate == 0.0 ||
	    rates->mrate == 0.0) {
		fprintf(stderr, "NB: one or more rates changed from 0 to "
		    "%f\n", ALMOST_ZERO);
	}
	rates->erate = (rates->erate == 0.0) ? ALMOST_ZERO : rates->erate;
	rates->srate = (rates->srate == 0.0) ? ALMOST_ZERO : rates->srate;
	rates->irate = (rates->irate == 0.0) ? ALMOST_ZERO : rates->irate;
	rates->mrate = (rates->mrate == 0.0) ? ALMOST_ZERO : rates->mrate;

	fprintf(stderr, "    error rate:         %.10f\n", rates->erate);
	fprintf(stderr, "    substitution rate:  %.10f\n", rates->srate);
	fprintf(stderr, "    indel rate:         %.10f\n", rates->irate);
	fprintf(stderr, "    match rate:         %.10f\n", rates->mrate);
}

static void
parse_rates_string(char *rates_string, struct rates *rates)
{
	char *estr, *sstr, *istr, *mstr;

	estr = strtok(rates_string, ",");
	sstr = strtok(NULL, ",");
	istr = strtok(NULL, ",");
	mstr = strtok(NULL, "");

	if (estr == NULL || sstr == NULL || istr == NULL || mstr == NULL) {
		fprintf(stderr, "error: failed to parse -r argument\n");
		exit(1);
	}

	rates->erate = atof(estr);
	rates->srate = atof(sstr);
	rates->irate = atof(istr);
	rates->mrate = atof(mstr);

	/* permit 7.00 as 0.07 */
	if (rates->erate > 1.0)
		rates->erate /= 100.0; 
	if (rates->srate > 1.0)
		rates->srate /= 100.0; 
	if (rates->irate > 1.0)
		rates->irate /= 100.0; 
	if (rates->mrate > 1.0)
		rates->mrate /= 100.0; 

	if (rates->erate > 1.0 || rates->srate > 1.0 || rates->irate > 1.0 ||
	    rates->mrate > 1.0) {
		fprintf(stderr, "error: user-provided rate(s) > 1.0\n");
		exit(1);
	}

	if (rates->erate == 0.0 || rates->srate == 0.0 || rates->irate == 0.0 ||
	    rates->mrate == 0.0) {
		fprintf(stderr, "NB: one or more rates changed from 0 to "
		    "%f\n", ALMOST_ZERO);
	}
	rates->erate = (rates->erate == 0.0) ? ALMOST_ZERO : rates->erate;
	rates->srate = (rates->srate == 0.0) ? ALMOST_ZERO : rates->srate;
	rates->irate = (rates->irate == 0.0) ? ALMOST_ZERO : rates->irate;
	rates->mrate = (rates->mrate == 0.0) ? ALMOST_ZERO : rates->mrate;

	fprintf(stderr, "    error rate:         %.10f\n", rates->erate);
	fprintf(stderr, "    substitution rate:  %.10f\n", rates->srate);
	fprintf(stderr, "    indel rate:         %.10f\n", rates->irate);
	fprintf(stderr, "    match rate:         %.10f\n", rates->mrate);
}

static void
single_pass_cb(char *path, struct stat *sb, void *arg)
{
	struct pass_cb *pcb = (struct pass_cb *)arg;

	read_file(path);
	pcb->nbytes += sb->st_size;
	pcb->nfiles++;
	PROGRESS_BAR(stderr, pcb->nfiles, pcb->total_files, 100);
}

/*
 * Do everything in one pass, saving it all in memory. This is potentially
 * much faster (factor of two or so), but can eat a ton of memory.
 *
 * This mode does _not_ require that all matches for a single individual read
 * are in the same file.
 */
static void
do_single_pass(char **objs, int nobjs, uint64_t files)
{
	struct pass_cb pcb;

	fprintf(stderr, "warning: single pass mode may consume a large amount "
	    "of memory!\n\n");

	memset(&pcb, 0, sizeof(pcb));
	pcb.total_files = files;
	fprintf(stderr, "Parsing %" PRIu64 " file(s)...\n", files);
	PROGRESS_BAR(stderr, 0, 0, 100);
	file_iterator_n(objs, nobjs, single_pass_cb, &pcb);
	PROGRESS_BAR(stderr, files, files, 100);
	fprintf(stderr, "\nParsed %.2f MB in %" PRIu64 " file(s).\n",
	    (double)pcb.nbytes / (1024 * 1024), pcb.nfiles);

	/*
	 * First iterative pass: calculate rates for top matches of reads.
	 */
	if (rates_file != NULL) {
		fprintf(stderr, "\nUsing user-defined rates file...\n");
		load_rates(rates_file, &pcb.rates);
	} else if (rates_string != NULL) {
		fprintf(stderr, "\nUsing user-defined rates...\n");
	} else {
		fprintf(stderr, "\nCalculating top match rates...\n");
		PROGRESS_BAR(stderr, 0, 0, 100);
		dynhash_iterate(read_list, calc_rates, &pcb.rates);
		PROGRESS_BAR(stderr, total_unique_reads, total_unique_reads, 100);
	}

	if (total_unique_reads == 0 && rates_string == NULL) {
		fprintf(stderr, "error: no matches were found in input "
		    "file(s)\n");
		exit(1);
	}

	if (rates_string == NULL)
		ratestats(&pcb.rates);
	else
		parse_rates_string(rates_string, &pcb.rates);

	/*
	 * Second iterative pass: Determine probabilities for all reads' best
	 * alignments.
	 */
	fprintf(stderr, "\nGenerating output...\n");
	PROGRESS_BAR(stderr, 0, 0, 10);
	dynhash_iterate(read_list, calc_probs, &pcb.rates);
	PROGRESS_BAR(stderr, total_unique_reads, total_unique_reads, 10);
	if (Bflag)
		putc('\n', stderr);
}

static void
double_pass_cb(char *path, struct stat *sb, void *arg)
{
	struct pass_cb *pcb = (struct pass_cb *)arg;

	read_file(path);
	pcb->nbytes += sb->st_size;
	pcb->nfiles++;
	PROGRESS_BAR(stderr, pcb->nfiles, pcb->total_files, 100);

	/* iterate over what's been read and free stored data */
	assert(pcb->pass == 1 || pcb->pass == 2);
	if (pcb->pass == 1)
		dynhash_iterate(read_list, calc_rates, &pcb->rates);
	else
		dynhash_iterate(read_list, calc_probs, &pcb->rates);
	cleanup();
}

/*
 * Do everything in two passes, once per file to calculate rates, and another
 * time per file to do the probability calculations.
 *
 * This mode requires that all matches for a single individual read are in
 * the same file.
 */
static void
do_double_pass(char **objs, int nobjs, uint64_t files)
{
	struct pass_cb pcb;
	int tmp_top_matches;

	/*
	 * First iterative pass: Calculate rates.
	 */
	if (rates_file != NULL) {
		fprintf(stderr, "\nUsing user-defined rates file...\n");
		load_rates(rates_file, &pcb.rates);
	} else if (rates_string != NULL) {
		fprintf(stderr, "\nUsing user-defined rates...\n");
	} else {
		/* only need best match for calculating rates */
		tmp_top_matches = top_matches;
		top_matches = 1;

		memset(&pcb, 0, sizeof(pcb));
		pcb.pass = 1;
		pcb.total_files = files;

		fprintf(stderr, "PASS 1: Parsing %" PRIu64 " file(s) to calculate "
		    "rates...\n", files);
		PROGRESS_BAR(stderr, 0, 0, 100);
		file_iterator_n(objs, nobjs, double_pass_cb, &pcb);
		PROGRESS_BAR(stderr, files, files, 100);
		fprintf(stderr, "\nParsed %.2f MB in %" PRIu64 " file(s).\n",
		    (double)pcb.nbytes / (1024 * 1024), pcb.nfiles);

		/* restore desired top matches */
		top_matches = tmp_top_matches;
	}

	if (total_unique_reads == 0 && rates_string == NULL) {
		fprintf(stderr, "error: no matches were found in input "
		    "file(s)\n");
		exit(1);
	}

	if (rates_string == NULL)
		ratestats(&pcb.rates);
	else
		parse_rates_string(rates_string, &pcb.rates);

	/*
	 * Second iterative pass: Determine probabilities for all reads' best
	 * alignments.
	 */
	pcb.pass = 2;
	pcb.nbytes = pcb.nfiles = 0;

	fprintf(stderr, "\n%sParsing %" PRIu64 " file(s) to calculate "
	    "probabilities...\n", (rates_file == NULL) ? "PASS 2: " : "", files);
	PROGRESS_BAR(stderr, 0, 0, 100);
	file_iterator_n(objs, nobjs, double_pass_cb, &pcb);
	PROGRESS_BAR(stderr, files, files, 100);
	fprintf(stderr, "\nParsed %.2f MB in %" PRIu64 " file(s).\n",
	    (double)pcb.nbytes / (1024 * 1024), pcb.nfiles);
}

static void
count_files(char *path, struct stat *sb, void *arg)
{
	uint64_t *i = (uint64_t *)arg;

	/* shut up, icc */
	(void)path;
	(void)sb;

	(*i)++;
}

static void
usage(char *progname)
{
	char *slash;

	slash = strrchr(progname, '/');
	if (slash != NULL)
		progname = slash + 1;

	fprintf(stderr, "usage: %s [-g rates_file] [-n normodds_cutoff] [-o pgenome_cutoff] "
	    "[-p pchance_cutoff] [-r erate,srate,irate,mrate] [-s normodds|pgenome|pchance] "
	    "[-t top_matches] [-B] [-G] [-R] [-S] "
	    "total_genome_len results_dir1|results_file1 "
	    "results_dir2|results_file2 ...\n", progname);
	exit(1);
}

int
main(int argc, char **argv)
{
	char *progname;
	uint64_t total_files;
	int ch;

	set_mode_from_argv(argv);

	fprintf(stderr, "--------------------------------------------------"
	    "------------------------------\n");
	fprintf(stderr, "probcalc. (SHRiMP %s [%s])\n",
	    SHRIMP_VERSION_STRING, get_compiler());
	fprintf(stderr, "--------------------------------------------------"
	    "------------------------------\n");

	progname = argv[0];
	while ((ch = getopt(argc, argv, "n:o:p:g:r:s:t:BGRS")) != -1) {
		switch (ch) {
		case 'g':
			rates_file = xstrdup(optarg);
			break;
		case 'n':
			normodds_cutoff = atof(optarg);
			break;
		case 'o':
			pgenome_cutoff = atof(optarg);
			break;
		case 'p':
			pchance_cutoff = atof(optarg);
			break;
		case 'r':
			rates_string = xstrdup(optarg);
			break;
		case 's':
			if (strstr(optarg, "pchance") != NULL) {
				sort_field = SORT_PCHANCE;
			} else if (strstr(optarg, "pgenome") != NULL) {
				sort_field = SORT_PGENOME;
			} else if (strstr(optarg, "normodds") != NULL) {
				sort_field = SORT_NORMODDS;
			} else {
				fprintf(stderr, "error: invalid sort "
				    "criteria\n");
				usage(progname);
			}
			break;
		case 't':
			top_matches = atoi(optarg);
			break;
		case 'B':
			Bflag = true;
			break;
		case 'G':
			Gflag = true;
			break;
		case 'R':
			Rflag = true;
			break;
		case 'S':
			Sflag = true;
			break;
		default:
			usage(progname);
		}
	}
	argc -= optind;
	argv += optind;

	if (rates_file != NULL && rates_string != NULL) {
		fprintf(stderr, "error: -g and -r flags cannot be mixed\n");
		exit(1);
	}

	if (Gflag && rates_file != NULL) {
		fprintf(stderr, "error: -G and -g flags cannot be mixed\n");
		exit(1);
	}

	if (argc < 2)
		usage(progname);
	
	if (!is_number(argv[0])) {
		fprintf(stderr, "error: expected a number for genome length "
		    "(got: %s)\n", argv[0]);
		usage(progname);
	}

	genome_len = strtoul(argv[0], NULL, 0);

	argc--;
	argv++;

	fprintf(stderr, "\nSettings:\n");
	fprintf(stderr, "    Normodds Cutoff:  < %f\n", normodds_cutoff);
	fprintf(stderr, "    Pgenome Cutoff:   < %f\n", pgenome_cutoff);
	fprintf(stderr, "    Pchance Cutoff:   > %f\n", pchance_cutoff);
	fprintf(stderr, "    Sort Criteria:      %s\n",
	    (sort_field == SORT_PCHANCE)  ? "pchance"  :
	    (sort_field == SORT_PGENOME)  ? "pgenome"  :
	    (sort_field == SORT_NORMODDS) ? "normodds" : "<unknown>");
	fprintf(stderr, "    Top Matches:        %d\n", top_matches);
	fprintf(stderr, "    Genome Length:      %" PRIu64 "\n", genome_len);

	fprintf(stderr, "\n");

	if (Gflag)
		fprintf(stderr, "NOTICE: Calculating rates only.\n");

	read_list = dynhash_create(keyhasher, keycomparer);
	if (read_list == NULL) {
		fprintf(stderr, "error: failed to allocate read_list\n");
		exit(1);
	}

	contig_cache = dynhash_create(keyhasher, keycomparer);
	if (contig_cache == NULL) {
		fprintf(stderr, "error: failed to allocate contig_cache\n");
		exit(1);
	}

	read_seq_cache = dynhash_create(keyhasher, keycomparer);
	if (read_seq_cache == NULL) {
		fprintf(stderr, "error: failed to allocate read_seq_cache\n");
		exit(1);
	}

	read_edit_cache = dynhash_create(keyhasher, keycomparer);
	if (read_edit_cache == NULL) {
		fprintf(stderr, "error: failed to allocate read_edit_cache\n");
		exit(1);
	}

	total_files = 0;
	file_iterator_n(argv, argc, count_files, &total_files);

	if (total_files == 0) {
		fprintf(stderr, "error: no files found to parse!\n");
		exit(1);
	}

	if (Sflag)
		do_single_pass(argv, argc, total_files);
	else
		do_double_pass(argv, argc, total_files);

	return (0);
}
