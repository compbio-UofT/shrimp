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
load_fasta(const char *file, void (*base_func)(int, ssize_t, int, char *)) {
	char buf[512], name[512];
	struct stat sb;
	FILE *fp;
	ssize_t len;
	int i, isnewentry;

	if (stat(file, &sb)) {
		fprintf(stderr, "error: failed to stat file [%s]: %s\n", file,
		    strerror(errno));
		return (-1);
	}

	/* tell consumer how many bytes worth of data we can maximally have */
	base_func(-1, sb.st_size, -1, NULL);

	fp = fopen(file, "r");
	if (fp == NULL) {
		fprintf(stderr, "error: failed to open file [%s]: %s\n", file,
		    strerror(errno));
		return (-1);
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

			if (buf[i] == 'T')
				continue;

			a = buf[i] - '0';
			assert(a >= 0 && a <= 4);

			base_func(a, len, isnewentry, name);
			isnewentry = 0;
			len++;
		}
	}

	fclose(fp);

	return (len);
}
