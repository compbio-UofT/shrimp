/*	$Id$	*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

int
main(int argc, char **argv)
{
	char buf[512];
	FILE *in, *out;
	unsigned int readsdone, readsper;

	if (argc != 3) {
		fprintf(stderr, "error: need 2 parameters\n");
		exit(1);
	}

	readsper = atoi(argv[1]);

	in = fopen(argv[2], "r");
	if (in == NULL) {
		perror("fopen(input)");
		exit(1);
	}

	out = NULL;
	readsdone = 0;
	while (fgets(buf, sizeof(buf), in) != NULL) {
		if (buf[0] == '#')
			continue;

		if (buf[0] == '>' && (readsdone % readsper) == 0) {
			char fname[512];

			if (out != NULL)	
				fclose(out);

			snprintf(fname, sizeof(fname), "%u_to_%u.csfasta",
			    readsdone, readsdone + readsper - 1);

			out = fopen(fname, "w");
			if (out == NULL) {
				perror("fopen(out)");
				exit(1);
			}
		}

		if (fwrite(buf, strlen(buf), 1, out) != 1) {
			perror("fwrite(out)");
			exit(1);
		}

		if (buf[0] == '>')
			readsdone++;
	}

	return (0);
}
