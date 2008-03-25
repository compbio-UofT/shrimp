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

ssize_t
load_fasta(const char *file, void (*bf)(int, ssize_t, int, char *, int), int s)
{
	char buf[4096], name[512];
	char translate[256];

	struct stat sb;
	FILE *fp;
	ssize_t len;
	int i, isnewentry, initbp;

	initbp = 0;

	assert(s == COLOUR_SPACE || s == LETTER_SPACE);

	if (stat(file, &sb)) {
		fprintf(stderr, "error: failed to stat file [%s]: %s\n", file,
		    strerror(errno));
		return (-1);
	}

	if (!S_ISREG(sb.st_mode)) {
		fprintf(stderr, "error: [%s] is not a regular file\n", file);
		return (-1);
	}

	/* tell consumer how many bytes worth of data we can maximally have */
	bf(FASTA_ALLOC, sb.st_size, -1, NULL, -1);

	fp = fopen(file, "r");
	if (fp == NULL) {
		fprintf(stderr, "error: failed to open file [%s]: %s\n", file,
		    strerror(errno));
		return (-1);
	}

	memset(translate, -1, sizeof(translate));

	if (s == COLOUR_SPACE) {
		translate['0'] = BASE_0;
		translate['1'] = BASE_1;
		translate['2'] = BASE_2;
		translate['3'] = BASE_3;
		translate['4'] = BASE_N;
		translate['N'] = BASE_N;
		translate['.'] = BASE_N;
		translate['X'] = BASE_X;
	} else {	
		translate['A'] = BASE_A;
		translate['C'] = BASE_C;
		translate['G'] = BASE_G;
		translate['T'] = BASE_T;
		translate['U'] = BASE_U;
		translate['M'] = BASE_M;
		translate['R'] = BASE_R;
		translate['W'] = BASE_W;
		translate['S'] = BASE_S;
		translate['Y'] = BASE_Y;
		translate['K'] = BASE_K;
		translate['V'] = BASE_V;
		translate['H'] = BASE_H;
		translate['D'] = BASE_D;
		translate['B'] = BASE_B;
		translate['N'] = BASE_N;
		translate['.'] = BASE_N;
		translate['X'] = BASE_X;
	}

	len = 0;
	isnewentry = 0;
	name[0] = '\0';
	while (fgets(buf, sizeof(buf), fp) != NULL) {
		bool have_initbp = false;

		if (buf[0] == '#')
			continue;

		if (buf[0] == '>') {
			char *nl;
			strncpy(name, buf + 1, sizeof(name) - 1);
			name[sizeof(name) - 1] = '\0';
			nl = strchr(name, '\n');
			if (nl != NULL)
				*nl = '\0';

			/* strip extra whitespace (e.g. dos newlines) */
			i = strlen(name) - 1;
			while (i >= 0) {
				if (!isspace((int)name[i]))
					break;
				name[i] = '\0';
			}

			isnewentry = 1;
			continue;
		}

		for (i = 0; buf[i] != '\n' && buf[i] != '\0'; i++) {
			int a;

			if (isspace((int)buf[i]))
				continue;

			buf[i] = (char)toupper((int)buf[i]);

			if (!have_initbp && s == COLOUR_SPACE) {
				if (buf[i] == 'A' || buf[i] == 'C' ||
				    buf[i] == 'G' || buf[i] == 'T') {
					switch (buf[i]) {
					case 'A': initbp = BASE_A; break;
					case 'C': initbp = BASE_C; break;
					case 'G': initbp = BASE_G; break;
					case 'T': initbp = BASE_T; break;
					}
					have_initbp = true;
					continue;
				}
			}

			a = translate[(int)buf[i]];
			if (a == -1) {
				fprintf(stderr, "error: invalid character ");
				if (isprint((int)buf[i]))
					fprintf(stderr, "(%c) ", buf[i]);
				else
					fprintf(stderr, "(0x%x) ", buf[i]);
				fprintf(stderr, "in input file [%s]\n", file);
				fprintf(stderr, "       (Did you mix up letter "
				    "space and colour space binaries?)\n");
				exit(1);
			}

			if (s == COLOUR_SPACE)
				assert(a >= BASE_CS_MIN && a <= BASE_CS_MAX);
			else
				assert(a >= BASE_LS_MIN && a <= BASE_LS_MAX);

			bf(a, len, isnewentry, name, initbp);
			isnewentry = 0;
			len++;
		}
	}

	/* tell consumer that we're done so they can cleanup, etc */
	bf(FASTA_DEALLOC, -1, -1, NULL, -1);

	fclose(fp);

	return (len);
}

/*
 * Give BASE_x, return the appropriate character.
 *
 * NB: Since we're limited to 4-bits, BASE_X returns 'N'.
 */
char
base_translate(int base, bool use_colours)
{
	char cstrans[] = { '0', '1', '2', '3', 'N' };
	char lstrans[] = { 'A', 'C', 'G', 'T', 'U', 'M', 'R', 'W',
			   'S', 'Y', 'K', 'V', 'H', 'D', 'B', 'N' };

	if (use_colours) {
		assert(base >= BASE_CS_MIN && base <= BASE_CS_MAX);
		return (cstrans[base]);
	} else {
		assert(base >= BASE_LS_MIN && base <= BASE_LS_MAX);
		return (lstrans[base]);
	}
}
