/*	$Id$	*/

#define LETTER_SPACE	1
#define COLOUR_SPACE	2

#define BASE_A		0
#define BASE_C		1
#define BASE_G		2
#define BASE_T		3

#define BASE_0		0
#define BASE_1		1
#define BASE_2		2
#define BASE_3		3

#define BASE_N		4

ssize_t load_fasta(const char *, void (*)(int, ssize_t, int, char *, int), int);
