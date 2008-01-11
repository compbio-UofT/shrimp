/*	$Id$	*/

#include <assert.h>
#include <ctype.h>
#include <dirent.h>
#include <unistd.h>

#include <sys/types.h>
#include <sys/stat.h>

#define MAX(_a, _b) ((_a) > (_b) ? (_a) : (_b))
#define MIN(_a, _b) ((_a) < (_b) ? (_a) : (_b))

/*
 * NB: This entire software collection assumes 2 bases packed into one byte.
 */
#define EXTRACT(_genome, _i) (((_genome)[(_i) / 8] >> (4 * ((_i) % 8))) & 0xf)
#define BPTO32BW(_x) (((_x) + 7) / 8)

uint64_t	rdtsc(void);
double		cpuhz(void);
u_int		strchrcnt(const char *, const char);
bool		is_number(const char *);
void		xstat(const char *, struct stat *);
void	       *xmalloc(size_t);
void	       *xrealloc(void *, size_t);
char	       *xstrdup(const char *);
uint32_t	hash_string(const char *);
double		ls_factorial(u_int);
double		ls_choose(int64_t, int64_t);
char	       *trim_brackets(char *);
void		bitfield_prepend(uint32_t *, uint32_t, uint32_t);
void		bitfield_append(uint32_t *, uint32_t, uint32_t);
void		progress_bar(FILE *, uint64_t, uint64_t, uint);
void		reverse_complement(uint32_t *, uint32_t *, uint32_t);
uint64_t	file_iterator(char *, void (*)(char *, struct stat *, void *),
		    void *);
uint64_t	file_iterator_n(char **, int,
		    void (*)(char *, struct stat *, void *), void *);

/* for optarg (and to shut up icc) */
extern char *optarg;
extern int   optind;

static inline int
complement_base(int base)
{
	int cmpl[5] = { 3, 2, 1, 0, 4 }; /* A->T, G->C, G->C, T->A, N->N */

	/* XXX - convert anything non-{A,C,G,T} to N */
	if (base > 4)
		base = 4;

	assert(base >= 0 && base <= 5);

	return (cmpl[base]);
}

/*
 * Given the first letter of a pair corresponding to a colour, and the colour
 * itself, obtain the next letter.
 */
static inline int
cstols(int first_letter, int colour)
{
	assert(first_letter >= 0 && first_letter <= 3);
	assert(colour >= 0 && colour <= 3);

	if ((first_letter % 2) == 0)
		return ((4 + first_letter + colour) % 4);
	else
		return ((4 + first_letter - colour) % 4);
}

static inline int
lstocs(int first_letter, int second_letter)
{
	const int colourmat[5][5] = {
		{ 0, 1, 2, 3, 4 },
		{ 1, 0, 3, 2, 4 },
		{ 2, 3, 0, 1, 4 },
		{ 3, 2, 1, 0, 4 },
		{ 4, 4, 4, 4, 4 }
	};

	assert(first_letter  >= 0 && first_letter  <= 15);
	assert(second_letter >= 0 && second_letter <= 15);

	/* XXX - convert anything non-{A,C,G,T} to N */
	if (first_letter > 4)
		first_letter = 4;
	if (second_letter > 4)
		second_letter = 4;

	return (colourmat[first_letter][second_letter]);
}
