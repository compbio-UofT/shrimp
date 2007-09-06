/*	$Id$	*/

#define BASE_A	0
#define BASE_C	1
#define BASE_G	2
#define BASE_T	3
#define BASE_N	4

ssize_t load_fasta(const char *, void (*)(int, ssize_t, int, char *));
