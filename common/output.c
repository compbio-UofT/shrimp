/*	$Id$	*/

#include <assert.h>
#include <ctype.h>
#include <inttypes.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include "../common/fasta.h"
#include "../common/sw-full-common.h"
#include "../common/output.h"
#include "../common/util.h"

char *
readtostr(uint32_t *read, u_int len, bool use_colours, int initbp)
{
	static char *buf;
	static u_int buflen, i, j;

	if (buflen < len) {
		if (buf != NULL)
			free(buf);
		buf = (char *)xmalloc(len + 1);
		buflen = len + 1;
	}

	i = 0;

	if (use_colours)
		buf[i++] = base_translate(initbp, false);

	for (j = 0; j < len; j++) {
		if (use_colours)
			buf[i + j] = base_translate(EXTRACT(read, j), true);
		else
			buf[i + j] = base_translate(EXTRACT(read, j), false); 
	}

	buf[i + j] = '\0';

	return (buf);
}

/* This is so ugly it hurts. */
void
output_pretty(FILE *fp, const char *readname, const char *contigname,
    const struct sw_full_results *sfr, const char *dbalign, const char *qralign,
    uint32_t *genome, uint32_t genome_len, uint32_t goff, bool use_colours,
    uint32_t *read, u_int readlen, int initbp, bool revcmpl)
{
	char *gpre, *gpost, *lspre, *lspost, *mpre, *nospace = "";
	int j, len;
	uint32_t genome_start, genome_end, read_start, read_end;

	/* shut up, icc */
	(void)readname;
	(void)genome;
	(void)contigname;

	genome_start = goff + sfr->genome_start;
	genome_end   = goff + sfr->genome_start + sfr->gmapped - 1;

	assert(genome_len > genome_start && genome_len > genome_end);

	if (revcmpl) {
		uint32_t tmp = genome_start;
		genome_start = genome_len - genome_end - 1;
		genome_end = genome_len - tmp - 1;
	}

	read_start = sfr->read_start;
	read_end   = sfr->read_start + sfr->rmapped - 1;

	len = strlen(dbalign);
	assert(len == strlen(qralign));

	gpre = gpost = lspre = lspost = mpre = nospace;
	if (read_start > 0) {
		gpre  = (char *)xmalloc(read_start + 1);
		lspre = (char *)xmalloc(read_start + 1);
		mpre  = (char *)xmalloc(read_start + 1);
		for (j = 0; j < read_start; j++) {
			if (genome_start + j > read_start)
				gpre[j] = base_translate(EXTRACT(genome,
				    genome_start - read_start + j), false);
			else
				gpre[j] = '-';
			lspre[j] = '-';
			mpre[j] = ' ';
		}
		gpre[j] = lspre[j] = mpre[j] = '\0';
	}
	if (read_end < (readlen - 1)) {
		gpost  = (char *)xmalloc(readlen - read_end);
		lspost = (char *)xmalloc(readlen - read_end);
		for (j = 0; j < (readlen - read_end - 1); j++) {
			if (genome_end + 1 + j < genome_len)
				gpost[j] = base_translate(EXTRACT(genome,
				    genome_end + 1 + j), false);
			else
				gpost[j] = '-';
			lspost[j] = '-';
		}
		gpost[j] = lspost[j] = '\0';
	}

	/* NB: internally 0 is first position, output uses 1. adjust. */

	fprintf(fp, "G: %10" PRId64 "    %s%s%s    %-10" PRId64 "\n",
	    (revcmpl) ? (int64_t)genome_end + 1 : (int64_t)genome_start + 1,
	    gpre, dbalign, gpost,
	    (revcmpl) ? (int64_t)genome_start + 1 : (int64_t)genome_end + 1);

	fprintf(fp, "%16s %s", "", mpre);
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

	if (use_colours) {
		fprintf(fp, "T: %10s    %s%s%s\n", "", lspre, qralign, lspost);
	} else {
		fprintf(fp, "R: %10u    %s%s%s    %-10u\n", read_start + 1,
		    lspre, qralign, lspost, read_end + 1);
	}

	if (use_colours) {
		char *rstr;

		fprintf(fp, "R: %10u   ", read_start + 1);
		rstr = readtostr(read, readlen, use_colours, initbp);
		assert(strlen(rstr) > 1);
		fputc(*rstr++, fp);
		for (j = 0; j < read_start; j++) {
			assert(*rstr != '\0');
			fputc(*rstr++, fp);
		}
		for (j = 0; *rstr != '\0';) {
			if (qralign[j] == '-')
				fputc('-', fp);
			else
				fputc(*rstr++, fp);
			if (qralign[j] != '\0')
				j++;
		}
		fprintf(fp, "    %-10u\n", read_end + 1);
	}

	if (gpre != nospace)
		free(gpre);
	if (gpost != nospace)
		free(gpost);
	if (lspre != nospace)
		free(lspre);
	if (lspost != nospace)
		free(lspost);
	if (mpre != nospace)
		free(mpre);
}

/* Escape '=' and '"' using '\', otherwise our parser will have trouble. */
static char *
escapestr(const char *str)
{
	char *buf;
	int i, j, len;
	bool escape = false;

	for (i = len = 0; str[i] != '\0'; i++) {
		if (str[i] == '=' || str[i] == '"') {
			len++;
			escape = true;
		}
		len++;
	}
	len++;

	if (!escape)
		return (NULL);

	buf = (char *)xmalloc(len);

	for (i = j = 0; str[i] != '\0'; i++, j++) {
		if (str[i] == '=' || str[i] == '"')
			buf[j++] = '\\';
		buf[j] = str[i];
	}
	buf[j] = '\0';

	assert(j >= 0 && j < len);

	return (buf);
}

void
output_normal(FILE *fp, const char *readname, const char *contigname,
    const struct sw_full_results *sfr, uint32_t genome_len, uint32_t goff,
    bool use_colours, uint32_t *read, u_int readlen, int initbp, bool revcmpl,
    bool print_readseq)
{
	uint32_t genome_start, genome_end;
	const char *cname, *rname;

	assert(fp != NULL && readname != NULL);
	assert(contigname != NULL && read != NULL);

	genome_start = goff + sfr->genome_start;
	genome_end = goff + sfr->genome_start + sfr->gmapped - 1;

	assert(genome_len > genome_start && genome_len > genome_end);

	if (revcmpl) {
		uint32_t tmp = genome_start;
		genome_start = genome_len - genome_end - 1;
		genome_end = genome_len - tmp - 1;
	}

	rname = escapestr(readname);
	if (rname == NULL)
		rname = readname;
	
	cname = escapestr(contigname);
	if (cname == NULL)
		cname = contigname;

	fprintf(fp, "> r=\"%s\" g=\"%s\" g_str=%c ", rname, cname,
	    (revcmpl) ? '-' : '+');

	if (print_readseq) {
		fprintf(fp, "r_seq=\"%s\" ",
		    readtostr(read, readlen, use_colours, initbp));
	}

	/* NB: internally 0 is first position, output uses 1. adjust. */
	fprintf(fp, "score=%d g_start=%u g_end=%u r_start=%d r_end=%d r_len=%u "
	    "match=%d subs=%d ins=%d dels=%d", sfr->score, genome_start + 1,
	    genome_end + 1, sfr->read_start + 1,
	    sfr->read_start + sfr->rmapped - 1 + 1, readlen, sfr->matches,
	    sfr->mismatches, sfr->insertions, sfr->deletions);
	  
	if (use_colours)
		fprintf(fp, " xovers=%d", sfr->crossovers);

	if (rname != readname)
		free((void *)rname);
	if (cname != contigname)
		free((void *)cname);
}
