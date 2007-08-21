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

#include "fasta.h"

ssize_t
load_fasta(const char *file, void (*bf)(int, ssize_t, int, char *, char), int s)
{
	char buf[512], name[512];
	char translate[256];

	struct stat sb;
	FILE *fp;
	ssize_t len;
	int i, isnewentry;
	char initbp;

	initbp = 0;

	assert(s == COLOUR_SPACE || s == LETTER_SPACE);

	if (stat(file, &sb)) {
		fprintf(stderr, "error: failed to stat file [%s]: %s\n", file,
		    strerror(errno));
		return (-1);
	}

	/* tell consumer how many bytes worth of data we can maximally have */
	bf(-1, sb.st_size, -1, NULL, -1);

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
	} else {	
		translate['A'] = BASE_A;
		translate['C'] = BASE_C;
		translate['G'] = BASE_G;
		translate['T'] = BASE_T;
		translate['N'] = BASE_N;
	}

	len = 0;
	isnewentry = 0;
	name[0] = '\0';
	while (fgets(buf, sizeof(buf), fp) != NULL) {
		if (buf[0] == '#')
			continue;

		if (buf[0] == '>') {
			char *nl;
			strncpy(name, buf + 1, sizeof(name) - 1);
			name[sizeof(name) - 1] = '\0';
			nl = strchr(name, '\n');
			if (nl != NULL)
				*nl = '\0';
			isnewentry = 1;
			continue;
		}
			
		for (i = 0; buf[i] != '\n' && buf[i] != '\0'; i++) {
			int a;

			buf[i] = (char)toupper((int)buf[i]);

			if (s == COLOUR_SPACE) {
				if (buf[i] == 'A' || buf[i] == 'C' ||
				    buf[i] == 'G' || buf[i] == 'T') {
					initbp = buf[i];
					continue;
				}
			}

			a = translate[(int)buf[i]];
			if (a == -1) {
				fprintf(stderr, "error: invalid character (%c) "
				    "in input file [%s]\n", buf[i], file);
				exit(1);
			}

			assert(a >= 0 && a <= 7);

			bf(a, len, isnewentry, name, initbp);
			isnewentry = 0;
			len++;
		}
	}

	fclose(fp);

	return (len);
}
