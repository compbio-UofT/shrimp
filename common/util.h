/*	$Id$	*/

#include <assert.h>

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
void	       *xmalloc(size_t);
char	       *xstrdup(const char *);
unsigned int	hash_string(const char *);
double		factorial(u_int);
double		ls_choose(int64_t, int64_t);
char	       *trim_brackets(char *);
void		bitfield_prepend(uint32_t *, uint32_t, uint32_t);
void		bitfield_append(uint32_t *, uint32_t, uint32_t);
char	       *extract_reftig(char *);

/* for optarg (and to shut up icc) */
extern char *optarg;
extern int   optind;

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

	assert(first_letter  >= 0 && first_letter  <= 4);
	assert(second_letter >= 0 && second_letter <= 4);

	return (colourmat[first_letter][second_letter]);
}
