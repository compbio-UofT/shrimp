/*	$Id$	*/

#include <assert.h>
#include <ctype.h>
#include <errno.h>
#include <fcntl.h>
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
#include "../common/util.h"

fasta_t
fasta_open(const char *file, int space)
{
	fasta_t fasta;
	struct stat sb;
	FILE *fp;

	assert(space == COLOUR_SPACE || space == LETTER_SPACE);

	if (stat(file, &sb))
		return (NULL);

	if (!S_ISREG(sb.st_mode))
		return (NULL);

	fp = fopen(file, "r");
	if (fp == NULL)
		return (NULL);

	fasta = (fasta_t)xmalloc(sizeof(*fasta));
	memset(fasta, 0, sizeof(*fasta));

	fasta->fp = fp;
	fasta->file = xstrdup(file);
	fasta->space = space;
	memset(fasta->translate, -1, sizeof(fasta->translate));

	if (space == COLOUR_SPACE) {
		fasta->translate['0'] = BASE_0;
		fasta->translate['1'] = BASE_1;
		fasta->translate['2'] = BASE_2;
		fasta->translate['3'] = BASE_3;
		fasta->translate['4'] = BASE_N;
		fasta->translate['N'] = BASE_N;
		fasta->translate['.'] = BASE_N;
		fasta->translate['X'] = BASE_X;
	} else {	
		fasta->translate['A'] = BASE_A;
		fasta->translate['C'] = BASE_C;
		fasta->translate['G'] = BASE_G;
		fasta->translate['T'] = BASE_T;
		fasta->translate['U'] = BASE_U;
		fasta->translate['M'] = BASE_M;
		fasta->translate['R'] = BASE_R;
		fasta->translate['W'] = BASE_W;
		fasta->translate['S'] = BASE_S;
		fasta->translate['Y'] = BASE_Y;
		fasta->translate['K'] = BASE_K;
		fasta->translate['V'] = BASE_V;
		fasta->translate['H'] = BASE_H;
		fasta->translate['D'] = BASE_D;
		fasta->translate['B'] = BASE_B;
		fasta->translate['N'] = BASE_N;
		fasta->translate['.'] = BASE_N;
		fasta->translate['X'] = BASE_X;
	}

	return (fasta);
}

void
fasta_close(fasta_t fasta)
{

	fclose(fasta->fp);
	free(fasta->file);
	free(fasta);
}

static char *
extract_name(char *buffer)
{

	assert(buffer[0] == '>');

	return (xstrdup(strtrim(&buffer[1])));
}

bool
fasta_get_next(fasta_t fasta, char **name, char **sequence)
{
	int i;
	bool gotname = false;
	uint32_t sequence_length = 0;
	uint32_t max_sequence_length = 0;

	*name = *sequence = NULL;

	if (fasta->leftover) {
		*name = extract_name(fasta->buffer);
		fasta->leftover = false;
		gotname = true;
	}

	while (fgets(fasta->buffer, sizeof(fasta->buffer), fasta->fp) != NULL) {
		if (fasta->buffer[0] == '#')
			continue;

		if (fasta->buffer[0] == '>') {
			if (gotname) {
				fasta->leftover = true;
				break;
			}

			*name = extract_name(fasta->buffer);
			gotname = true;
			continue;
		}

		for (i = 0; fasta->buffer[i] != '\0'; i++) {
			if (fasta->buffer[i] == '\n') {
				fasta->buffer[i] = '\0';
				break;
			} else if (isspace((int)fasta->buffer[i])) {
				continue;
			}

			if (sequence_length == max_sequence_length) {
				if (max_sequence_length == 0)
					max_sequence_length = 1024*1024;
				else if (max_sequence_length >= 128*1024*1024)
					max_sequence_length += 128*1024*1024;
				else
					max_sequence_length *= 2;
				*sequence = (char *)xrealloc(*sequence, max_sequence_length);
			}
			(*sequence)[sequence_length++] = fasta->buffer[i];
		}
	}

	if (!gotname) {
		if (*sequence != NULL)
			free(*sequence);
		return (false);
	}

	/* ensure nul-termination */
	if (sequence_length == max_sequence_length)
		*sequence = (char *)xrealloc(*sequence, max_sequence_length + 1);
	(*sequence)[sequence_length] = '\0';

	assert(*name != NULL);
	assert(*sequence != NULL);
	return (true);
}

int
fasta_get_initial_base(fasta_t fasta, char *sequence)
{
	char c;

	assert(fasta->space == COLOUR_SPACE);

	c = (char)toupper((int)*sequence);
	switch (c) {
	case 'A': return (BASE_A); break;
	case 'C': return (BASE_C); break;
	case 'G': return (BASE_G); break;
	case 'T': return (BASE_T); break;
	}

	return (-1);
}

uint32_t *
fasta_bitfield_to_colourspace(fasta_t fasta, uint32_t *source, uint32_t length)
{
	int a, lastbp = BASE_T;
	uint32_t *dst;
	uint32_t i;

	assert(fasta->space == LETTER_SPACE);

	dst = (uint32_t *)xmalloc(BPTO32BW(length) * sizeof(uint32_t));

	for (i = 0; i < length; i++) {
		a = EXTRACT(source, i);
		bitfield_insert(dst, i, lstocs(lastbp, a));
		lastbp = a;
	}

	return (dst);
}

uint32_t *
fasta_sequence_to_bitfield(fasta_t fasta, char *sequence)
{
	uint32_t i, length, idx;
	uint32_t *bitfield;
	int a;
	char c;

	length = strlen(sequence);
	bitfield = (uint32_t *)xmalloc(BPTO32BW(length) * sizeof(uint32_t));

	for (i = idx = 0; i < length; i++) {
		c = (char)toupper((int)sequence[i]);

		if (i == 0 && fasta->space == COLOUR_SPACE) {
			if (c != 'A' && c == 'C' && 
			    c != 'G' && c != 'T') {
				free(bitfield);
				return (NULL);
			}

			continue;
		}

		a = fasta->translate[(int)c];
		if (a == -1) {
			fprintf(stderr, "error: invalid character ");
			if (isprint((int)c))
				fprintf(stderr, "(%c) ", c);
			else
				fprintf(stderr, "(0x%x) ", c);
			fprintf(stderr, "in input file [%s]\n", fasta->file);
			fprintf(stderr, "       (Did you mix up letter "
			    "space and colour space programs?)\n");
			exit(1);
		}

		if (fasta->space == COLOUR_SPACE) {
			assert((a >= BASE_CS_MIN && a <= BASE_CS_MAX) ||
			    (a == BASE_N || a == BASE_X));
		} else {
			assert((a >= BASE_LS_MIN && a <= BASE_LS_MAX) ||
			    (a == BASE_N || a == BASE_X));
		}

		bitfield_append(bitfield, idx++, a);
	}

	if (fasta->space == COLOUR_SPACE)
		assert(idx == length - 1);
	else
		assert(idx == length);

	return (bitfield);
}

/*
 * Give BASE_x, return the appropriate character.
 *
 * NB: Since we're limited to 4-bits, BASE_X returns 'N'.
 */
char
base_translate(int base, bool use_colours)
{
	/*
	 * NB: colour-space only valid for 0-3 and BASE_N/BASE_X
	 *     BASE_N is reported as a skipped cycle: '.' in CS.
	 */
	char cstrans[] = { '0', '1', '2', '3', '!', '@', '#', '$',
			   '%', '^', '&', '*', '?', '~', ';', '.' };
	char lstrans[] = { 'A', 'C', 'G', 'T', 'U', 'M', 'R', 'W',
			   'S', 'Y', 'K', 'V', 'H', 'D', 'B', 'N' };

	if (use_colours) {
		assert((base >= BASE_CS_MIN && base <= BASE_CS_MAX) ||
		    (base == BASE_N || base == BASE_X));
		return (cstrans[base]);
	} else {
		assert((base >= BASE_LS_MIN && base <= BASE_LS_MAX) ||
		    (base == BASE_N || base == BASE_X));
		return (lstrans[base]);
	}
}
