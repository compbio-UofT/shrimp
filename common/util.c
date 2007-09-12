/*	$Id$	*/

#include <assert.h>
#include <errno.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <unistd.h>

#include <sys/time.h>
#include <sys/types.h>

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

unsigned int
hash_string(const char *x)
{
	unsigned int hash = 0;

	while (*x != '\0')
		hash = 31 * hash + *x++;

	return (hash);
}

/* factorial using stirling's approximation after 20 */
double
factorial(u_int n)
{
	const double fact[21] = {
		1,
		1,
		2,
		6,
		24,
		120,
		720,
		5040,
		40320,
		362880,
		3628800,
		39916800,
		479001600,
		6227020800,
		87178291200,
		1307674368000,
		20922789888000,
		355687428096000,
		6402373705728000,
		121645100408832000,
		2432902008176640000
	};
	double a, b;

	if (n <= 20)
		return (fact[n]);

	a = sqrt(2 * M_PI * n);
	b = pow(n / M_E, n);

	return (a * b);
}

/* choose in log space */
double
ls_choose(int64_t n, int64_t k)
{
	double a, b, c;

	if (k < 0 || k > n)
		return (0);

	a = log(factorial(n));
	b = log(factorial(k));
	c = log(factorial(n - k));

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

/*
 * Extract the reftig name from a comment within the file.
 *
 * NB: This function will modify the string in place if it's a 'REFTIG: '
 *     line.
 */
char *
extract_reftig(char *line)
{
	char *s;

	if (strncmp(line, "# REFTIG: [", 11) == 0) {
		line += 11;
		s = line;

		while (*line != ']' && *line != '\0')
			line++;
		*line = '\0';

		return (s);
	}

	return (NULL);
}

