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

#include "../common/lookup.h"
#include "../common/util.h"

#include "finalpass.h"

static lookup_t read_list;		/* list of reads and top matches */
static lookup_t reftig_cache;		/* reftig name cache */

/* fields in our output */
enum {
	FIELD_NAME = 0,
	FIELD_SCORE,
	FIELD_INDEX,
	FIELD_READ_START,
	FIELD_READ_MAPPED_LENGTH,
	FIELD_MATCHES,
	FIELD_MISMATCHES,
	FIELD_INSERTIONS,
	FIELD_DELETIONS,
	FIELD_CROSSOVERS
};
#define NFIELDS	10

struct readstats {
	char    *reftig;
	int32_t  score;
	uint32_t index;
	uint32_t read_start;
	uint32_t read_mapped_length;
	uint32_t matches;
	uint32_t mismatches;
	uint32_t insertions;
	uint32_t deletions;
	uint32_t crossovers;
};

struct readinfo {
	char		*name;
	struct readstats top_matches[0];
};

struct rates {
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
	struct readstats *rs;
	double            pval;
};

static uint64_t total_alignments;	/* total outputs in all files */
static uint64_t total_unique_reads;	/* total unique reads in all files */

/* arguments */
static double   pval_cutoff = DEF_PVAL_CUTOFF;
static int	top_matches = DEF_TOP_MATCHES;
static uint64_t genome_len;

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
reheap(struct readstats *stats, int node)
{
	struct readstats tmp;
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
save_match(struct readinfo *ri, char *reftig, int32_t score, uint32_t index,
    uint32_t read_start, uint32_t read_mapped_length, uint32_t matches,
    uint32_t mismatches, uint32_t insertions, uint32_t deletions,
    uint32_t crossovers)
{
	struct readstats *stats = ri->top_matches;

	if (score < stats[1].score)
		return;

	stats[1].reftig			= reftig;
	stats[1].score			= score;
	stats[1].index 			= index;
	stats[1].read_start		= read_start;
	stats[1].read_mapped_length	= read_mapped_length;
	stats[1].matches		= matches;
	stats[1].mismatches		= mismatches;
	stats[1].insertions		= insertions;
	stats[1].deletions		= deletions;
	stats[1].crossovers		= crossovers;

	reheap(stats, 1);
}

static void
progress_bar(uint32_t at, uint32_t of)
{
	static int lastperc, beenhere;
	static char whirly = '\\';

	char progbuf[52];
	int perc, i, j, dec;

	perc = (at * 10000) / of;

	if (beenhere && perc - lastperc != 1)
		return;

	beenhere = 1;
	lastperc = perc;

	dec = perc % 100;
	perc /= 100;

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

	fprintf(stderr, "\rProgress: [%s] %3d.%02d%%", progbuf, perc, dec);
	fflush(stderr);
}

static char *
trim_brackets(char *str)
{

	if (str[0] == '[')
		str++;
	if (str[strlen(str) - 1] == ']')
		str[strlen(str) - 1] = '\0';

	return (str);
}

static void
parse_line(char *buf, const char *fname, const int line, char *reftig)
{
	char *fields[NFIELDS];
	char *field, *name;
	struct readinfo *ri;
	int32_t score;
	int i;

	memset(fields, 0, sizeof(fields));

	field = strtok(buf, " \t\r\n");
	for (i = 0; i < (sizeof(fields) / sizeof(fields[0])); i++) {
		if (field == NULL)
			break;
		fields[i] = field;
		field = strtok(NULL, " \t\r\n");
	}

	/*
	 * The last valid field for letter space and colour space alignments
	 * must be present.
	 */
	if (fields[FIELD_DELETIONS] == NULL) {
		fprintf(stderr, "error: file format error: file [%s] line %d\n",
		    fname, line);
		exit(1);
	}

	/* Fake any colour space fields if we're in letter space */
	if (fields[FIELD_CROSSOVERS] == NULL)
		fields[FIELD_CROSSOVERS] = "0";

	total_alignments++;

	name = trim_brackets(fields[FIELD_NAME]);

	if (lookup_find(read_list, name, NULL, (void *)&ri)) {
		/* test if score better than worst, if so, add it */
		score = atoi(fields[FIELD_SCORE]);
		if (score > ri->top_matches[1].score) {
			save_match(ri, reftig, score,
			    atoi(fields[FIELD_INDEX]),
			    atoi(fields[FIELD_READ_START]),
			    atoi(fields[FIELD_READ_MAPPED_LENGTH]),
			    atoi(fields[FIELD_MATCHES]),
			    atoi(fields[FIELD_MISMATCHES]),
			    atoi(fields[FIELD_INSERTIONS]),
			    atoi(fields[FIELD_DELETIONS]),
			    atoi(fields[FIELD_CROSSOVERS]));
		}
	} else {
		total_unique_reads++;

		/* alloc, init, insert in lookup */
		ri = xmalloc(sizeof(*ri) +
		    sizeof(ri->top_matches[0]) * (top_matches + 1));
		memset(ri, 0, sizeof(*ri) +
		    sizeof(ri->top_matches[0]) * (top_matches + 1));

		ri->name = xstrdup(name);

		ri->top_matches[0].score = top_matches;
		for (i = 1; i <= top_matches; i++)
			ri->top_matches[i].score = 0x80000000 + i;

		save_match(ri, reftig,
		    atoi(fields[FIELD_SCORE]),
		    atoi(fields[FIELD_INDEX]),
		    atoi(fields[FIELD_READ_START]),
		    atoi(fields[FIELD_READ_MAPPED_LENGTH]),
		    atoi(fields[FIELD_MATCHES]),
		    atoi(fields[FIELD_MISMATCHES]),
		    atoi(fields[FIELD_INSERTIONS]),
		    atoi(fields[FIELD_DELETIONS]),
		    atoi(fields[FIELD_CROSSOVERS]));

		if (lookup_add(read_list, ri->name, ri) == false) {
			fprintf(stderr, "error: failed to add to lookup - "
			    "probably out of memory\n");
			exit(1);
		}
	}
}

/* extract the reftig name from a comment within the file */
static char *
extract_reftig(char *line, char *curreftig)
{
	char *s, *key;

	if (strncmp(line, "# REFTIG: [", 11) == 0) {
		line += 11;
		s = line;
		while (*line != ']' && *line != '\0')
			line++;
		*line = '\0';

		if (lookup_find(reftig_cache, s, (void *)&key, NULL))
			return (key);

		s = xstrdup(s);

		if (lookup_add(reftig_cache, s, NULL) == false) {
			fprintf(stderr, "error: failed to add to reftig cache "
			    "- probably out of memory\n");
			exit(1);
		}

		return (s);
	}

	return (curreftig);
}

/*
 * if no reftig name was included in the input file, derive it from the
 * filename. this is a hack for data produced before the reftig name was
 * included.
 */
static char *
derive_reftig(char *file)
{
	char *s, *dot, *key;

	s = xstrdup(file);
	dot = strstr(s, ".fa");
	if (dot != NULL)
		*dot = '\0';

	if (lookup_find(reftig_cache, s, (void *)&key, NULL)) {
		free(s);
		return (key);
	}

	if (lookup_add(reftig_cache, s, NULL) == false) {
		fprintf(stderr, "error: failed to add to reftig cache "
		    "- probably out of memory\n");
		exit(1);
	}
	
	return (s);
}

/* Returns the number of bytes read */
static uint64_t
read_file(const char *dir, const char *file)
{
	char fpath[PATH_MAX + 1];
	char buf[512];
	struct stat sb;
	char *reftig;
	FILE *fp;
	u_int line;

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

	reftig = NULL;
	for (line = 0; fgets(buf, sizeof(buf), fp) != NULL; line++) {
		if (buf[0] == '#') {
			reftig = extract_reftig(buf, reftig);
		} else if (buf[0] != '\n') {
			if (reftig == NULL)
				reftig = derive_reftig((char *)file);

			parse_line(buf, file, line, reftig);
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

static void
calc_rates(void *arg, void *key, void *val)
{
	static uint64_t calls;

	struct rates *rates = arg;
	struct readinfo *ri = val;
	struct readstats *rs;
	int32_t best, prev_best;
	int i;

	(void)key;

	progress_bar((int)calls, (int)total_unique_reads);	/* XXX */
	calls++;

	best = prev_best = 0;
	for (i = 1; i <= ri->top_matches[0].score; i++) {
		if (ri->top_matches[i].score > best) {
			prev_best = best;
			best = i;
		}
	}

	/*
	 * Only count the best score with no equivalent best scores.
	 */
	if (best != 0 && (prev_best == 0 ||
	    ri->top_matches[best].score >
	    ri->top_matches[prev_best].score)) {
		rs = &ri->top_matches[best];

		rates->total_len  += rs->matches + rs->mismatches;
		rates->insertions += rs->insertions;
		rates->deletions  += rs->deletions;
		rates->matches    += rs->matches;
		rates->mismatches += rs->mismatches;
		rates->crossovers += rs->crossovers;
	}
}

static double
p_nohits(uint64_t l, int k, int nsubs, int nerrors, int nindels)
{
	double r;

	r  = ls_choose(k - 2, nsubs) + (nsubs * log(3));
	r += ls_choose(k - 2 - nsubs, nerrors) + (nerrors * log(3));
	r += ls_choose(k - 1 - nsubs - nerrors, nindels) + (nindels * log(2));

	/* (1 - (r/4^k))^l */
	r -= (k * log(4));
	r = pow(M_E, r);
	r = 1 - r;
	r = l * log(r);

	return (pow(M_E, r));
}

static double
p_otherswrong(int k, int nerrors, double erate, int nsubs, double subrate,
    int nindels, double indelrate, int nmatches, double matchrate)
{
	double r;

	r  = ls_choose(k, nerrors) + (nerrors * log(erate));
	r += ls_choose(k - nerrors, nsubs) + (nsubs * log(subrate));
	r += ls_choose(k - nerrors - nsubs, nindels) +
	    (nindels * log(indelrate));
	r += nmatches * log(matchrate);

	return (pow(M_E, r));
}

static int
rspvcmp(const void *a, const void *b)
{
	const struct readstatspval *rspva = a, *rspvb = b;

	if (rspva->pval == rspvb->pval)
		return (0);
	else if (rspva->pval < rspva->pval)
		return (1);
	else
		return (-1);
}

static void
calc_probs(void *arg, void *key, void *val)
{
	static struct readstatspval *rspv;

	struct rates *rates = arg;
	struct readinfo *ri = val;
	double r, norm;
	int i, j, rlen;

	(void)key;

	if (rspv == NULL)
		rspv = xmalloc(sizeof(rspv[0]) * (top_matches + 1));

	/* 1: Calculate P(others wrong) for each hit */
	norm = 0;
	for (i = 1; i <= top_matches; i++) {
		struct readstats *rs = &ri->top_matches[i];

		if (rs->score < 0)
			continue;

		rspv[i - 1].rs = rs;
		rlen = rs->matches + rs->mismatches;
		rspv[i - 1].pval = p_otherswrong(rlen, rs->crossovers,
		    rates->erate, rs->mismatches, rates->srate,
		    rs->insertions + rs->deletions, rates->irate,
		    rs->matches, rates->mrate);
		norm += rspv[i - 1].pval;
	}

	/* 2: Normalise our values */
	for (i = 1; i <= top_matches; i++) {
		struct readstats *rs = &ri->top_matches[i];

		if (rs->score < 0)
			continue;

		rspv[i - 1].pval /= norm;
	}

	/* 3: For each match calculate the probably it's bad */
	j = 0;
	for (i = 1; i <= top_matches; i++) {
		struct readstats *rs = &ri->top_matches[i];

		if (rs->score < 0)
			continue;

		rlen = rs->matches + rs->mismatches;
		r = p_nohits(genome_len, rlen, rs->mismatches, rs->crossovers,
		    rs->insertions + rs->deletions);

		r *= rspv[i - 1].pval;
		r = 1 - r;

		rspv[j].rs = rs;
		rspv[j++].pval = r;
	}

	/* 4: Sort the values ascending */
	if (j > 0)
		qsort(rspv, j, sizeof(rspv), rspvcmp);

	/* 5: Finally, print out values in ascending order until the cutoff */
	for (i = 0; i < j; i++) {
		struct readstats *rs = rspv[i].rs;

		if (rspv[i].pval >= pval_cutoff)
			break;

		/* handle letter and colour spaces separately */
		if (rates->crossovers == 0) {
			printf("[%s] [%s] %.8f %d %u %u %u %u %u %u %u\n",
			    ri->name, rs->reftig, rspv[i].pval, rs->score,
			    rs->index, rs->read_start, rs->read_mapped_length,
			    rs->matches, rs->mismatches, rs->insertions,
			    rs->deletions);
		} else {
			printf("[%s] [%s] %.8f %d %u %u %u %u %u %u %u %u\n",
			    ri->name, rs->reftig, rspv[i].pval, rs->score,
			    rs->index, rs->read_start, rs->read_mapped_length,
			    rs->matches, rs->mismatches, rs->insertions,
			    rs->deletions, rs->crossovers);
		}
	}
}

static void
usage(char *progname)
{
	char *slash;

	slash = strrchr(progname, '/');
	if (slash != NULL)
		progname = slash + 1;

	fprintf(stderr, "usage: %s [-p pval] [-t top_matches] genome_len "
	    "results_directory\n", progname);
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

	progname = argv[0];
	while ((ch = getopt(argc, argv, "p:t:")) != -1) {	
		switch (ch) {
		case 'p':
			pval_cutoff = atof(optarg);
			break;
		case 't':
			top_matches = atoi(optarg);
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

	read_list = lookup_create(keyhasher, keycomparer, false);
	if (read_list == NULL) {
		fprintf(stderr, "error: failed to allocate read_list\n");
		exit(1);
	}

	reftig_cache = lookup_create(keyhasher, keycomparer, false);
	if (reftig_cache == NULL) {
		fprintf(stderr, "error: failed to allocate reftig_cache\n");
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
	while (1) {
		de = readdir(dp);
		if (de == NULL)
			break;

		if (de->d_type == DT_REG && strcmp(de->d_name, ".") &&
		    strcmp(de->d_name, ".."))
			bytes += read_file(dir, de->d_name);

		progress_bar(i, files);

		i++;
	}
	fprintf(stderr, "\nParsed %.2f MB in %" PRIu64 " files.\n",
	    (double)bytes / (1024 * 1024), i);

	closedir(dp);

	/*
	 * First iterative pass: calculate rates for top matches of reads.
	 */

	fprintf(stderr, "\nCalculating top match rates...\n");
	memset(&rates, 0, sizeof(rates));
	lookup_iterate(read_list, calc_rates, &rates);
	fprintf(stderr, "\nCalculated rates for %" PRIu64 " reads\n",
	    total_unique_reads);

	fprintf(stderr, "total matches:      %" PRIu64 "\n", total_alignments);
	fprintf(stderr, "total unique reads: %" PRIu64 "\n",total_unique_reads);
	fprintf(stderr, "total length:       %" PRIu64 "\n", rates.total_len); 
	fprintf(stderr, "insertions:         %" PRIu64 "\n", rates.insertions); 
	fprintf(stderr, "deletions:          %" PRIu64 "\n", rates.deletions); 
	fprintf(stderr, "matches:            %" PRIu64 "\n", rates.matches); 
	fprintf(stderr, "mismatches:         %" PRIu64 "\n", rates.mismatches); 
	fprintf(stderr, "crossovers:         %" PRIu64 "\n", rates.crossovers); 

	rates.erate = (double)rates.crossovers / (double)rates.total_len;
	rates.srate = (double)rates.mismatches / (double)rates.total_len;
	rates.irate = (double)(rates.insertions + rates.deletions) /
	    (double)rates.total_len;
	rates.mrate = (double)rates.matches / (double)rates.total_len;

	fprintf(stderr, "error rate:         %.10f\n", rates.erate);
	fprintf(stderr, "substitution rate:  %.10f\n", rates.srate);
	fprintf(stderr, "indel rate:         %.10f\n", rates.irate);
	fprintf(stderr, "match rate:         %.10f\n", rates.mrate);

	/*
	 * Second iterative pass: Determine probabilities for all reads' best
	 * alignments.
	 */

	/* handle letter and colour spaces separately */
	if (rates.crossovers == 0) {
		printf("# [READ_NAME] [REFTIG_NAME] pval score index "
		    "read_start read_mapped_length matches mismatches "
		    "insertions deletions\n");
	} else {
		printf("# [READ_NAME] [REFTIG_NAME] pval score index "
		    "read_start read_mapped_length matches mismatches "
		    "insertions deletions crossovers\n");
	}
	lookup_iterate(read_list, calc_probs, &rates);

	return (0);
}
