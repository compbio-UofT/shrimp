/*	$Id$	*/

#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include <sys/types.h>
#include <sys/stat.h>

#define WIDTH 50

static char cmptbl[255];

/*
 * print the comment and fasta identifier lines out as they appear, except
 * append a '_revcmpl' string to the fasta identifier.
 *
 * XXX - remove existing '_revcmpl' is we're going back from an already
 * reversed file?
 */
static void
print_comments(char *str)
{
	int fastaline = 0, printed = 0;
	char last = '\0';

	while (*str != '\0') {
		if (last == '\n' || last == '\0')
			fastaline = (*str == '>') ? 1 : 0;
		
		if (fastaline && (!isprint((int)*str) || isspace((int)*str))) {
			puts("_revcmpl");
			printed = 1;
		}

		putchar(*str);
		last = *str++;
	}

	if (fastaline && !printed)
		fputs("_revcmpl", stdout);

	putchar('\n');
}

int
main(int argc, char **argv)
{
	FILE *fp;
	char *buf;
	struct stat sb;
	int i, j, skipline, skipidx, cpl;

	if (argc != 2) {
		fprintf(stderr, "need a filename, bud\n");
		exit(1);
	}

	if (stat(argv[1], &sb) == -1) {
		perror("stat");
		exit(1);
	}

	buf = (char *)malloc(sb.st_size);
	if (buf == NULL) {
		perror("malloc");
		exit(1);
	}

	fp = fopen(argv[1], "r");
	if (fp == NULL) {
		perror("fopen");
		exit(1);
	}

	if (fread(buf, sb.st_size, 1, fp) != 1) {
		perror("fread");
		exit(1);
	}

	memset(cmptbl, -1, sizeof(cmptbl));
	cmptbl['A'] = 'T';	cmptbl['a'] = 't';
	cmptbl['C'] = 'G';	cmptbl['c'] = 'g';
	cmptbl['G'] = 'C';	cmptbl['g'] = 'c';
	cmptbl['T'] = 'A';	cmptbl['t'] = 'a';
	cmptbl['N'] = 'N';	cmptbl['n'] = 'n';

	cmptbl['\0'] = cmptbl['\n'] = cmptbl['\r'] = 0;

	/* find start */
	skipidx = skipline = 0;
	for (i = 0; i < sb.st_size; i++) {
		if (buf[i] == '>' || buf[i] == '#') {
			skipline = 1;
			skipidx = i;
		}

		if (buf[i] == '\n') {
			/* preserve comments and identifiers */
			if (skipline) {
				buf[i] = '\0';
				print_comments(&buf[skipidx]);
			}
			skipline = 0;
		}

		if (!skipline)
			break;
	}

	/* work backwards */
	cpl = 0;
	for (j = sb.st_size - 1; j >= i; j--) {
		char cmpl = cmptbl[(int)buf[j]];

		if (cmpl == -1) {
			fprintf(stderr, "error translating to complement: %d\n",
			    buf[j]);
			exit(1);
		}

		if (cmpl != 0) {
			putchar(cmpl);
			if (++cpl == WIDTH) { 
				putchar('\n');
				cpl = 0;
			}
		}
	}
	if (cpl != 0)
		putchar('\n');

	fclose(fp);

	return (0);
}
