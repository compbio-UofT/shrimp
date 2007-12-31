/*	$Id$	*/

#define LETTER_SPACE	1
#define COLOUR_SPACE	2

/*
 * We're presently using 4 bits, so we max out at 16 different bases. This just
 * works out if we treat N and X as the same.
 *
 * from: http://genome.ucsc.edu/FAQ/FAQdownloads#download5
 */

#define BASE_A		0		/* Adenine */
#define BASE_C		1		/* Cytosine */
#define BASE_G		2		/* Guanine */
#define BASE_T		3		/* Thymine */
#define BASE_U		4		/* Uracil */
#define BASE_M		5		/* A or C */
#define BASE_R		6		/* A or G */
#define BASE_W		7		/* A or T */
#define BASE_S		8		/* C or G */
#define BASE_Y		9		/* C or T */
#define BASE_K		10		/* G or T */
#define BASE_V		11		/* A or C or G */
#define BASE_H		12		/* A or C or T */
#define BASE_D		13		/* A or G or T */
#define BASE_B		14		/* C or G or T */
#define BASE_X		15		/* G or A or T or C */
#define BASE_N		15		/* G or A or T or C */

#define BASE_LS_MIN	BASE_A
#define BASE_LS_MAX	BASE_N

#define BASE_0		0
#define BASE_1		1
#define BASE_2		2
#define BASE_3		3

#define BASE_CS_MIN	BASE_0
#define BASE_CS_MAX	BASE_3

/* overrides for the base argument to indicate beginning and end of file */
#define FASTA_ALLOC	-1
#define FASTA_DEALLOC	-2

ssize_t load_fasta(const char *, void (*)(int, ssize_t, int, char *, int), int);
char    base_translate(int, bool);
