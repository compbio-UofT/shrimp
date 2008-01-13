/*	$Id$	*/

#include <assert.h>
#include <errno.h>
#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <unistd.h>

#include <sys/time.h>
#include <sys/types.h>

#include "fasta.h"
#include "util.h"

uint64_t
rdtsc() {
	uint32_t lo, hi;

	__asm__ __volatile__ ("rdtsc" : "=a" (lo), "=d" (hi));

	return (((uint64_t)hi << 32) | lo);
}

double
cpuhz()
{
	uint64_t before, after;
	struct timeval tv1, tv2;
	int diff;

	/* XXX - abusive, poor man's calc; needs good (2ms) clock granularity */
	gettimeofday(&tv1, NULL);
	before = rdtsc();
	do {
		gettimeofday(&tv2, NULL);

		diff = tv2.tv_usec - tv1.tv_usec;
		if (diff < 0)
			diff = 1000000 - tv1.tv_usec + tv2.tv_usec;
	} while (diff < 2000);
	after = rdtsc();

	return (((double)(after - before) / diff) * 1.0e6);
}

u_int
strchrcnt(const char *str, const char c)
{
	int i;

	i = 0;
	while (*str != '\0') {
		if (*str++ == c)
			i++;
	}

	return (i);
}

bool
is_number(const char *str)
{

	while (*str != '\0')
		if (!isdigit((int)*str++))
			return (false);

	return (true);
}

void
xstat(const char *path, struct stat *sbp)
{
	
	if (stat(path, sbp) != 0) {
		fprintf(stderr, "error: failed to stat [%s]: %s\n", path,
		    strerror(errno));
		exit(1);
	}
}

void *
xmalloc(size_t size)
{
	void *ptr;

	ptr = malloc(size);
	if (ptr == NULL) {
		fprintf(stderr, "error: malloc failed: %s\n", strerror(errno));
		exit(1);
	}

	return (ptr);
} 

void *
xrealloc(void *ptr, size_t size)
{

	ptr = realloc(ptr, size);
	if (ptr == NULL) {
		fprintf(stderr, "error: realloc failed: %s\n", strerror(errno));
		exit(1);
	}

	return (ptr);
}

char *
xstrdup(const char *str)
{
	char *dup;

	assert(str != NULL);

	dup = strdup(str);
	if (dup == NULL) {
		fprintf(stderr, "error: strdup failed: %s\n", strerror(errno));
		exit(1);
	}

	return (dup);
}

uint32_t
hash_string(const char *x)
{
	uint32_t hash = 0;

	while (*x != '\0')
		hash = 31 * hash + *x++;

	return (hash);
}

/* factorial using stirling's approximation after 20 */
double
ls_factorial(u_int n)
{
	const double fact[21] = {
		1.0,
		1.0,
		2.0,
		6.0,
		24.0,
		120.0,
		720.0,
		5040.0,
		40320.0,
		362880.0,
		3628800.0,
		39916800.0,
		479001600.0,
		6227020800.0,
		87178291200.0,
		1307674368000.0,
		20922789888000.0,
		355687428096000.0,
		6402373705728000.0,
		121645100408832000.0,
		2432902008176640000.0
	};
	double a, b;

	if (n <= 20)
		return log (fact[n]);

	a = log( sqrt(2 * M_PI * n));
	b = n * log(n / M_E);

	return (a + b);
}

/* choose in log space */
double
ls_choose(int64_t n, int64_t k)
{
	double a, b, c;

	if (k < 0 || k > n)
		return (0);

	a = ls_factorial(n);
	b = ls_factorial(k);
	c = ls_factorial(n - k);

	return (a - (b + c));
}

char *
trim_brackets(char *str)
{

	if (str[0] == '[')
		str++;
	if (str[strlen(str) - 1] == ']')
		str[strlen(str) - 1] = '\0';

	return (str);
}

/*
 * Prepend the low 4 bits of 'val' to the start of the bitfield in 'bf'.
 * 'entries' is the maximum number of 4-bit words to be stored in the
 * bitfield.
 */
void
bitfield_prepend(uint32_t *bf, uint32_t entries, uint32_t val)
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
void
bitfield_append(uint32_t *bf, uint32_t entries, uint32_t val)
{
	uint32_t word;

	word = bf[entries / 8];
	word &= ~(0xf << (4 * (entries % 8)));
	word |= ((val & 0xf) << (4 * (entries % 8)));
	bf[entries / 8] = word;
}

void
progress_bar(FILE *output, uint64_t at, uint64_t of, uint incr)
{
	static int lastperc, beenhere;
	static char whirly = '\\';

	char progbuf[52];
	int perc, i, j, dec;

	if (at == 0 && of == 0) {
		beenhere = lastperc = 0;
		whirly = '\\';
		return;
	}

	perc = (at * 100 * incr) / of;

	if (beenhere && perc == lastperc)
		return;

	beenhere = 1;
	lastperc = perc;

	dec = perc % incr;
	perc /= incr;

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
	if (at >= of)
		whirly = '|';

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

	fprintf(output, "\rProgress: [%s] %3d.%02d%%", progbuf, perc, dec);
	fflush(output);
}

static inline uint32_t
swap_nibbles(uint32_t i)
{

	return (((i & 0xf0000000) >> 28) |
		((i & 0x0f000000) >> 20) |
		((i & 0x00f00000) >> 12) |
		((i & 0x000f0000) >>  4) |
		((i & 0x0000f000) <<  4) |
		((i & 0x00000f00) << 12) |
		((i & 0x000000f0) << 20) |
		((i & 0x0000000f) << 28));
}

/*
 * genome <- reverse_complement(genome)
 *
 * This is straightforward on purpose; it's suboptimal, but we won't be getting
 * any performance improvements by speeding it up.
 */
void
reverse_complement(uint32_t *g_ls, uint32_t *g_cs, uint32_t g_len)
{
	uint32_t i, j, tmp, fudge, up, down;

	assert(g_len != 0);
	assert(g_ls != NULL);

	/* First, swap all words and the nibbles within and complement them. */
	for (i = 0, j = BPTO32BW(g_len) - 1; i <= j; i++, j--) {
		tmp = g_ls[i];
		g_ls[i] = swap_nibbles(g_ls[j]);
		g_ls[j] = swap_nibbles(tmp);

		g_ls[i]=(complement_base((g_ls[i] & 0x0000000f) >>  0) <<  0) |
			(complement_base((g_ls[i] & 0x000000f0) >>  4) <<  4) |
			(complement_base((g_ls[i] & 0x00000f00) >>  8) <<  8) |
			(complement_base((g_ls[i] & 0x0000f000) >> 12) << 12) |
			(complement_base((g_ls[i] & 0x000f0000) >> 16) << 16) |
			(complement_base((g_ls[i] & 0x00f00000) >> 20) << 20) |
			(complement_base((g_ls[i] & 0x0f000000) >> 24) << 24) |
			(complement_base((g_ls[i] & 0xf0000000) >> 28) << 28);

		g_ls[j]=(complement_base((g_ls[j] & 0x0000000f) >>  0) <<  0) |
			(complement_base((g_ls[j] & 0x000000f0) >>  4) <<  4) |
			(complement_base((g_ls[j] & 0x00000f00) >>  8) <<  8) |
			(complement_base((g_ls[j] & 0x0000f000) >> 12) << 12) |
			(complement_base((g_ls[j] & 0x000f0000) >> 16) << 16) |
			(complement_base((g_ls[j] & 0x00f00000) >> 20) << 20) |
			(complement_base((g_ls[j] & 0x0f000000) >> 24) << 24) |
			(complement_base((g_ls[j] & 0xf0000000) >> 28) << 28);
	}

	/*
	 * If (g_len % 8) != 0, we need to shift all words down since
	 * the last word had some zeroed fields.
	 */
	fudge = 8 - (g_len % 8);
	if (fudge != 8) {
		down = fudge * 4;
		up = 32 - down;

		for (i = 0; i < BPTO32BW(g_len); i++) {
			if (i > 0)
				g_ls[i - 1] |= (g_ls[i] << up);
			g_ls[i] >>= down;
		}
	}

	/* If necessary, regenerate the colour-space genome */
	if (g_cs != NULL) {
		int base, lastbp;

		lastbp = BASE_T;
		for (i = 0; i < g_len; i++) {
			base = EXTRACT(g_ls, i);
			bitfield_append(g_cs, i, lstocs(lastbp, base));
			lastbp = base;
		}
	}
}

/*
 * Given a path, if it's a regular file (or symlink to one), call fh on it.
 * If it's a directory, call fh on all regular files within it (or symlinks to
 * regular files).
 *
 * Returns the number of files fh was called on.
 */
uint64_t
file_iterator(char *path, void (*fh)(char *, struct stat *, void *),
    void *arg)
{
	char fpath[2048];
	struct stat sb;
	DIR *dp;
	struct dirent *de;
	uint64_t files;

	/* is a regular file... */
	xstat(path, &sb);
	if (S_ISREG(sb.st_mode)) {
		fh(path, &sb, arg);
		return (1);
	}

	/* is (hopefully) a directory... */
	dp = opendir(path);
	if (dp == NULL) {
		fprintf(stderr, "error: failed to open directory [%s]: %s\n",
		    path, strerror(errno));
		exit(1);
	}

	files = 0;
	while (1) {
		de = readdir(dp);
		if (de == NULL)
			break;

		if (de->d_type != DT_REG && de->d_type != DT_LNK)
			continue;

		strcpy(fpath, path);
		if (fpath[strlen(path) - 1] != '/')
			strcat(fpath, "/");
		strcat(fpath, de->d_name);

		/* ensure it's a regular file or link to one */
		xstat(fpath, &sb);
		if (S_ISREG(sb.st_mode)) {
			fh(fpath, &sb, arg);
			files++;
		} else {
			fprintf(stderr, "warning: [%s] is neither a regular "
			    "file, nor a link to one; skipping...", fpath);
			continue;
		}
	}

	closedir(dp);

	return (files);
}

uint64_t
file_iterator_n(char **paths, int npaths,
    void (*fh)(char *, struct stat *, void *), void *arg)
{
	uint64_t files;
	int i;

	for (i = files = 0; i < npaths; i++)
		files += file_iterator(paths[i], fh, arg);

	return (files);
}
