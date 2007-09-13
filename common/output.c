/*	$Id$	*/

#include <assert.h>
#include <ctype.h>
#include <inttypes.h>
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include "../common/sw-full-common.h"
#include "../common/output.h"
#include "../common/util.h"

void
output_pretty(FILE *fp, const char *readname, const char *reftigname,
    const struct sw_full_results *sfr, const char *dbalign, const char *qralign,
    uint32_t *genome, uint32_t genome_len, uint32_t goff, int newread)
{
	static int firstcall = 1;

	int j, len;
	uint32_t aoff;
	char translate[5] = { 'A', 'C', 'G', 'T', 'N' };

	if (newread && !firstcall) {
		putc('\n', fp);
		putc('\n', fp);
	}
	firstcall = 0;

	if (newread) {
		fprintf(fp,
		    "------------------------------------------------------\n");
		fprintf(fp, "READ: [%s]\n", readname);
		fprintf(fp,
		    "------------------------------------------------------\n");
	}

	putc('\n', fp);

	aoff = sfr->genome_start;

	len = strlen(dbalign);
	assert(len == strlen(qralign));

	fprintf(fp, "Score:   %d   (read start/end: %d/%d)\n",
	    sfr->score, sfr->read_start,
	    sfr->read_start + sfr->rmapped - 1);
	fprintf(fp, "Index:   %u   (end: %u%s%s%s)\n", goff + aoff,
	    goff + sfr->genome_start + sfr->gmapped - 1,
	    (reftigname != NULL) ? ", reftig: [" : "",
	    (reftigname != NULL) ? reftigname : "",
	    (reftigname != NULL) ? "]" : "");
	
	fprintf(fp, "Reftig:  ");
	for (j = 4; j > 0; j--) {
		if (j <= goff + aoff)
			putc(translate[
			    EXTRACT(genome, goff + aoff - j)], fp);
	}
	fprintf(fp, "%s", dbalign);
	for (j = 0; j < 4; j++) {
		if ((goff + aoff + len + j) < genome_len)
			putc(translate[EXTRACT(genome,
			    goff + aoff + len + j)], fp);
	}
	putc('\n', fp);

	fprintf(fp, "Match:   ");
	for (j = 4; j > 0; j--) {
		if (j <= goff + aoff)
			putc(' ', fp);
	}
	for (j = 0; j < len; j++) {
		if (dbalign[j] == qralign[j] && dbalign[j] != '-') {
			putc('|', fp);
		} else {
			assert(j > 0 || dbalign[j] == toupper((int)qralign[j]));
			if (dbalign[j] == toupper((int)qralign[j]))
				putc('X', fp);
			else if (islower((int)qralign[j]))
				putc('x', fp);
			else
				putc(' ', fp);
		}
	}
	putc('\n', fp);

	fprintf(fp, "Read:    ");
	for (j = 4; j > 0; j--) {
		if (j <= goff + aoff)
			putc(' ', fp);
	}
	fprintf(fp, "%s\n", qralign);
}

void
output_normal(FILE *fp, const char *readname, const char *reftigname,
    const struct sw_full_results *sfr, uint32_t goff, int use_colours,
    int newread)
{
	static int firstcall = 1;

	/* shut up, icc */
	(void)reftigname;

	if (!firstcall && newread)
		putc('\n', fp);
	firstcall = 0;

	if (use_colours) {
		fprintf(fp, "[%s] %d %u %u %d %d %d %d %d %d %d\n", readname,
		    sfr->score, goff + sfr->genome_start,
		    goff + sfr->genome_start + sfr->gmapped - 1,
		    sfr->read_start, sfr->rmapped, sfr->matches,
		    sfr->mismatches, sfr->insertions, sfr->deletions,
		    sfr->crossovers);
	} else {
		fprintf(fp, "[%s] %d %u %u %d %d %d %d %d %d\n", readname,
		    sfr->score, goff + sfr->genome_start,
		    goff + sfr->genome_start + sfr->gmapped - 1,
		    sfr->read_start, sfr->rmapped, sfr->matches,
		    sfr->mismatches, sfr->insertions, sfr->deletions);
	}
}
