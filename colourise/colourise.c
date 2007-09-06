/*	$Id$	*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include <sys/types.h>

#include "fasta.h"

#define MAXCOLS		70
#define BASE_N_COLOUR	4

static int colourmat[5][5] = {
  {	0, 1, 2, 3,	BASE_N_COLOUR	},
  {	1, 0, 3, 2,	BASE_N_COLOUR	},
  {	2, 3, 0, 1,	BASE_N_COLOUR	},
  {	3, 2, 1, 0,	BASE_N_COLOUR	},
  { BASE_N_COLOUR, BASE_N_COLOUR, BASE_N_COLOUR, BASE_N_COLOUR, BASE_N_COLOUR }
};

static void
usage(char *progname)
{
	char *slash;

	slash = strrchr(progname, '/');
	if (slash != NULL)
		progname = slash + 1;

	fprintf(stderr, "usage: %s input_file\n", progname);
	exit(1);
}

static void
helper(int letter, ssize_t offset, int isnewentry, char *name)
{
	static int col = 0;
	static int last = 0;
	static int first = 1;

	/* shut up, icc */
	(void)offset;

	/* ignore resource allocation call */
	if (letter == -1)
		return;

	if (isnewentry) {
		col = 0;
		if (!first)
			printf("\n");
		printf(">%s\n", name);
		printf("T");
		last = BASE_T;
	}

	if (col++ == MAXCOLS) {
		col = 0;
		printf("\n");
	}

	printf("%d", colourmat[last][letter]);

	first = 0;
	last = letter;
}

int
main(int argc, char **argv)
{

	if (argc != 2)
		usage(argv[0]);

	return (load_fasta(argv[1], helper) == -1);
}
