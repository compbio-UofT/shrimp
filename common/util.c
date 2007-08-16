/*	$Id$	*/

#include <stdint.h>
#include <string.h>

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
