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

#include <sys/types.h>
#include <sys/stat.h>

#include "../common/sw-full-common.h"
#include "../common/input.h"
#include "../common/lookup.h"
#include "../common/util.h"
#include "../common/version.h"

#include "probcalc.h"

static lookup_t read_list;		/* list of reads and top matches */
static lookup_t contig_cache;		/* contig name cache */
static lookup_t read_seq_cache;		/* read sequence cache */

static bool Bflag = false;		/* print a progress bar */
static bool Rflag = false;		/* include read sequence in output */

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

static unsigned int
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
static uint64_t
read_file(const char *dir, const char *file)
{
	char fpath[2048];
	struct input input;
	struct stat sb;
	struct readinfo *ri;
	char *key;
	FILE *fp;
	int i;

	fpath[0] = '\0';
	strcpy(fpath, dir);
	strcat(fpath, "/");
	strcat(fpath, file);

	fp = fopen(fpath, "r");
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
		} else {
			if (lookup_find(read_seq_cache, input.read_seq,
			    (void *)&key, NULL)) {
				free(input.read_seq);
				input.read_seq = key;
			} else {
				if (lookup_add(read_seq_cache, input.read_seq,
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
		if (lookup_find(contig_cache, input.genome, (void *)&key,
		    NULL)) {
			free(input.genome);
			input.genome = key;
		} else {
			if (lookup_add(contig_cache, input.genome,
			    NULL) == false) {
				fprintf(stderr, "error: failed to add to "
				    "contig cache (%d) - probably out of "
				    "memory\n", __LINE__);
				exit(1);
			}
		}

		total_alignments++;

		if (lookup_find(read_list, input.read, (void *)&key,
		    (void *)&ri)) {
			/* use only one read name string to save space */
			free(input.read);
			input.read = key;

			/* test if score better than worst, if so, add it */
			if (input.score > ri->top_matches[1].score)
				save_match(ri, &input);
		} else {
			total_unique_reads++;

			/* alloc, init, insert in lookup */
			ri = xmalloc(sizeof(*ri) +
			    sizeof(ri->top_matches[0]) * (top_matches + 1));
			memset(ri, 0, sizeof(*ri) +
			    sizeof(ri->top_matches[0]) * (top_matches + 1));

			ri->name = input.read;

			ri->top_matches[0].score = top_matches;
			for (i = 1; i <= top_matches; i++)
				ri->top_matches[i].score = 0x80000000 + i;

			save_match(ri, &input);

			if (lookup_add(read_list, ri->name, ri) == false) {
				fprintf(stderr, "error: failed to add to lookup"
				    " - probably out of memory\n");
				exit(1);
			}
		}
	}

	if (fstat(fileno(fp), &sb) == -1) {
		fprintf(stderr, "error: failed to stat file [%s]: %s\n",
		    fpath, strerror(errno));
		exit(1);
	}

	fclose(fp);

	return (sb.st_size);
}

static double
p_chance(uint64_t l, int k, int nsubs, int nerrors, int nindels, int origlen)
{
	double r = 0;

	r += ls_choose(k - 2, nsubs+nerrors) + ((nsubs + nerrors) * log(4));

	/*
	 * two for ins and dels, four for the four possible letters to insert.
	 * above should really take length of indel into account
	 */
	r += ls_choose(k - 1, nindels) + (nindels * log(2*4));

	/* r is now log of fraction of matched words over total words */
	r -= (k * log(4));
	r = pow(M_E, r);
	r = 1 - r;

	/* if hit not full length we need to correct */
	l *= MAX(origlen - k+1, 1);

	/* 2 for the two strands */
	r = 2.0 * l * log(r);

	return (1- pow(M_E, r));
}



static void
calc_rates(void *arg, void *key, void *val)
{
	static uint64_t calls;

	struct rates *rates = arg;
	struct readinfo *ri = val;
	struct input *rs;
	double d;
	int32_t best;
	int i, rlen;

	(void)key;

	PROGRESS_BAR(stderr, calls, total_unique_reads, 100);
	calls++;

	best = 0;
	for (i = 1; i <= ri->top_matches[0].score; i++) {
		if (best == 0 ||
		    ri->top_matches[i].score >= ri->top_matches[best].score) {
			best = i;
		}
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
	r = (pow(M_E, r));


	if (r < 0.00000001)
		r = 0.00000001;
	if (r > 0.99999999)
		r = 0.99999999;

	return (r);
}

static int
rspvcmp(const void *a, const void *b)
{
	const struct readstatspval *rspva = a, *rspvb = b;
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
	static uint64_t calls;

	struct rates *rates = arg;
	struct readinfo *ri = val;
	double s, norm;
	int i, j, rlen;

	(void)key;

	PROGRESS_BAR(stderr, calls, total_unique_reads, 10);
	calls++;

	if (rspv == NULL)
		rspv = xmalloc(sizeof(rspv[0]) * (top_matches + 1));

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

	/* 4: Sort the values ascending */
	qsort(rspv, j, sizeof(*rspv), rspvcmp);

	/* 5: Finally, print out values in ascending order until the cutoff */
	for (i = 0; i < j; i++) {
		struct input *rs = rspv[i].rs;

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

		/* the only sane way is to reproduce the output code, sadly. */
		printf("> r=\"%s\" g=\"%s\" g_str=%c ", rs->read, rs->genome,
		    (rs->revcmpl) ? '-' : '+');

		if (rs->read_seq != NULL && rs->read_seq[0] != '\0')
			printf("read_seq=\"%s\" ", rs->read_seq);

		/* NB: internally 0 is first position, output is 1. adjust */
		printf("score=%d g_start=%u g_end=%u r_start=%u r_end=%u "
		    "r_len=%u match=%u subs=%u ins=%u dels=%u xovers=%u "
		    "normodds=%f pgenome=%f pchance=%f\n\n", rs->score,
		    rs->genome_start + 1, rs->genome_end + 1,
		    rs->read_start + 1, rs->read_end + 1, rs->read_length,
		    rs->matches, rs->mismatches, rs->insertions, rs->deletions,
		    rs->crossovers, rspv[i].normodds, rspv[i].pgenome,
		    rspv[i].pchance);
	}
}

static void
usage(char *progname)
{
	char *slash;

	slash = strrchr(progname, '/');
	if (slash != NULL)
		progname = slash + 1;

	fprintf(stderr, "usage: %s [-n normodds_cutoff] [-o pgenome_cutoff] "
	    "[-p pchance_cutoff] [-s normodds|pgenome|pchance] "
	    "[-t top_matches] [-B] [-R] total_genome_len results_directory\n",
	    progname);
	exit(1);
}

int
main(int argc, char **argv)
{
	struct stat;
	struct rates rates;
	DIR *dp;
	char *progname, *dir;
	struct dirent *de;
	uint64_t i, bytes, files;
	int ch;

	fprintf(stderr, "--------------------------------------------------"
	    "------------------------------\n");
	fprintf(stderr, "probcalc. (SHRiMP version %s)\n",
	    SHRIMP_VERSION_STRING);
	fprintf(stderr, "--------------------------------------------------"
	    "------------------------------\n");

	progname = argv[0];
	while ((ch = getopt(argc, argv, "n:o:p:s:t:BR")) != -1) {
		switch (ch) {
		case 'n':
			normodds_cutoff = atof(optarg);
			break;
		case 'o':
			pgenome_cutoff = atof(optarg);
			break;
		case 'p':
			pchance_cutoff = atof(optarg);
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
		case 'R':
			Rflag = true;
			break;
		default:
			usage(progname);
		}
	}
	argc -= optind;
	argv += optind;

	if (argc != 2)
		usage(progname);

	genome_len = strtoul(argv[0], NULL, 0);
	dir = argv[1];

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
	
	read_list = lookup_create(keyhasher, keycomparer, false);
	if (read_list == NULL) {
		fprintf(stderr, "error: failed to allocate read_list\n");
		exit(1);
	}

	contig_cache = lookup_create(keyhasher, keycomparer, false);
	if (contig_cache == NULL) {
		fprintf(stderr, "error: failed to allocate contig_cache\n");
		exit(1);
	}

	read_seq_cache = lookup_create(keyhasher, keycomparer, false);
	if (read_seq_cache == NULL) {
		fprintf(stderr, "error: failed to allocate read_seq_cache\n");
		exit(1);
	}

	/*
	 * Our parameter is a directory, since we may have tons of files -
	 * more than can fit in argv[].
	 */
	dp = opendir(dir);
	if (dp == NULL) {
		fprintf(stderr, "error: failed to open directory [%s]: %s\n",
		    dir, strerror(errno));
		exit(1);
	}

	/* count number of entries */
	files = 0;
	while (1) {
		de = readdir(dp);
		if (de == NULL)
			break;

		if (de->d_type == DT_REG && strcmp(de->d_name, ".") &&
		    strcmp(de->d_name, ".."))
			files++;
	}

	rewinddir(dp);

	i = bytes = 0;
	fprintf(stderr, "Parsing %" PRIu64 " files...\n", files);
	PROGRESS_BAR(stderr, 0, 0, 100);
	while (1) {
		de = readdir(dp);
		if (de == NULL)
			break;

		if (de->d_type == DT_REG && strcmp(de->d_name, ".") &&
		    strcmp(de->d_name, "..")) {
			bytes += read_file(dir, de->d_name);
			i++;
		}

		PROGRESS_BAR(stderr, i, files, 100);
	}
	PROGRESS_BAR(stderr, files, files, 100);
	fprintf(stderr, "\nParsed %.2f MB in %" PRIu64 " files.\n",
	    (double)bytes / (1024 * 1024), i);

	closedir(dp);

	/*
	 * First iterative pass: calculate rates for top matches of reads.
	 */

	fprintf(stderr, "\nCalculating top match rates...\n");
	PROGRESS_BAR(stderr, 0, 0, 100);
	memset(&rates, 0, sizeof(rates));
	lookup_iterate(read_list, calc_rates, &rates);
	PROGRESS_BAR(stderr, total_unique_reads, total_unique_reads, 100);
	fprintf(stderr, "\nCalculated rates for %" PRIu64 " reads\n",
	    total_unique_reads);

	fprintf(stderr, "\nStatistics:\n");
	fprintf(stderr, "    total matches:      %" PRIu64 "\n",
	    total_alignments);
	fprintf(stderr, "    total unique reads: %" PRIu64 "\n",
	    total_unique_reads);
	fprintf(stderr, "    total samples:      %" PRIu64 "\n",
	    rates.samples);
	fprintf(stderr, "    total length:       %" PRIu64 "\n",
	    rates.total_len); 
	fprintf(stderr, "    insertions:         %" PRIu64 "\n",
	    rates.insertions); 
	fprintf(stderr, "    deletions:          %" PRIu64 "\n",
	    rates.deletions); 
	fprintf(stderr, "    matches:            %" PRIu64 "\n",
	    rates.matches); 
	fprintf(stderr, "    mismatches:         %" PRIu64 "\n",
	    rates.mismatches); 
	fprintf(stderr, "    crossovers:         %" PRIu64 "\n",
	    rates.crossovers); 

	rates.erate = (double)rates.crossovers / (double)rates.total_len;
	rates.srate = (double)rates.mismatches / (double)rates.total_len;
	rates.irate = (double)(rates.insertions + rates.deletions) /
	    (double)rates.total_len;
	rates.mrate = (double)rates.matches / (double)rates.total_len;

	if (rates.erate == 0.0 || rates.srate == 0.0 || rates.irate == 0.0 ||
	    rates.mrate == 0.0) {
		fprintf(stderr, "NB: one or more rates changed from 0 to "
		    "1.0e-9\n");
	}
	rates.erate = (rates.erate == 0.0) ? 1.0e-9 : rates.erate;
	rates.srate = (rates.srate == 0.0) ? 1.0e-9 : rates.srate;
	rates.irate = (rates.irate == 0.0) ? 1.0e-9 : rates.irate;
	rates.mrate = (rates.mrate == 0.0) ? 10.e-9 : rates.mrate;

	fprintf(stderr, "    error rate:         %.10f\n", rates.erate);
	fprintf(stderr, "    substitution rate:  %.10f\n", rates.srate);
	fprintf(stderr, "    indel rate:         %.10f\n", rates.irate);
	fprintf(stderr, "    match rate:         %.10f\n", rates.mrate);

	/*
	 * Second iterative pass: Determine probabilities for all reads' best
	 * alignments.
	 */
	fprintf(stderr, "\nGenerating output...\n");
	PROGRESS_BAR(stderr, 0, 0, 10);
	lookup_iterate(read_list, calc_probs, &rates);
	PROGRESS_BAR(stderr, total_unique_reads, total_unique_reads, 10);
	if (Bflag)
		putc('\n', stderr);

	return (0);
}
