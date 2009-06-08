/*	$Id$	*/

#include <assert.h>
#include <ctype.h>
#include <dirent.h>
#include <unistd.h>
#include <zlib.h>

#include <sys/types.h>
#include <sys/stat.h>

#include "../common/fasta.h"

typedef enum {
	MODE_LETTER_SPACE = 1,
	MODE_COLOUR_SPACE = 2,
	MODE_HELICOS_SPACE= 3
} shrimp_mode_t;

extern shrimp_mode_t shrimp_mode;

#ifdef __GNUC__
#define __predict_false(_x)	__builtin_expect((_x), 0)
#define __predict_true(_x)	__builtin_expect((_x), 1)
#else
#define __predict_false(_x)	(_x)
#define __predict_true(_x)	(_x)
#endif

#define MAX(_a, _b) ((_a) > (_b) ? (_a) : (_b))
#define MIN(_a, _b) ((_a) < (_b) ? (_a) : (_b))

/*
 * NB: This entire software collection assumes 2 bases packed into one byte.
 */
#define EXTRACT(_genome, _i) (((_genome)[(_i) / 8] >> (4 * ((_i) % 8))) & 0xf)
#define BPTO32BW(_x) (((_x) + 7) / 8)

struct _strbuf_t {
	char   *string;
	u_int	string_length;
	u_int	string_alloced;
};
typedef struct _strbuf_t * strbuf_t;

void		set_mode_from_argv(char **);
const char     *get_mode_string(void);
uint64_t	gettimeinusecs(void);
uint64_t	rdtsc(void);
double		cpuhz(void);
u_int		strchrcnt(const char *, const char);
bool		is_number(const char *);
bool		is_whitespace(const char *);
void		xstat(const char *, struct stat *);
void	       *xmalloc(size_t);
void	       *xrealloc(void *, size_t);
char	       *xstrdup(const char *);
uint32_t	hash_string(const char *);
double		ls_factorial(u_int);
double		ls_choose(int64_t, int64_t);
char	       *trim_brackets(char *);
void		bitfield_prepend(uint32_t *, uint32_t, uint32_t);
void		bitfield_insert(uint32_t *, uint32_t, uint32_t);
void		bitfield_append(uint32_t *, uint32_t, uint32_t);
void		progress_bar(FILE *, uint64_t, uint64_t, uint);
void		reverse_complement(uint32_t *, uint32_t *, uint32_t);
uint64_t	file_iterator(char *, void (*)(char *, struct stat *, void *),
		    void *);
uint64_t	file_iterator_n(char **, int,
		    void (*)(char *, struct stat *, void *), void *);
char           *get_compiler(void);
char	       *strrev(char *);
char	       *strtrim(char *);
strbuf_t	strbuf_create(void);
char	       *strbuf_string(strbuf_t, int *);
void	        strbuf_append(strbuf_t, char *, ...);
void		strbuf_destroy(strbuf_t);
char	       *fast_gzgets(gzFile, char *, int);
char	       *comma_integer(uint64_t);

/* for optarg (and to shut up icc) */
extern char *optarg;
extern int   optind;

static inline int
complement_base(int base)
{
	static const int cmpl[16] = {
		BASE_T,		/* A -> T */
		BASE_G,		/* C -> G */
		BASE_C,		/* G -> C */
		BASE_A,		/* T -> A */
		BASE_A,		/* U -> A */
		BASE_K,		/* M -> K */
		BASE_Y,		/* R -> Y */
		BASE_W,		/* W -> W */
		BASE_S,		/* S -> S */
		BASE_R,		/* Y -> R */
		BASE_M,		/* K -> M */
		BASE_B,		/* V -> B */
		BASE_D,		/* H -> D */
		BASE_H,		/* D -> H */
		BASE_V,		/* B -> V */
		BASE_N,		/* X,N,- -> N */
	};

	assert((base >= BASE_LS_MIN && base <= BASE_LS_MAX) ||
	    (base == BASE_X || base == BASE_N));

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
	const int colourmat[4][4] = {
		{ 0, 1, 2, 3 },
		{ 1, 0, 3, 2 },
		{ 2, 3, 0, 1 },
		{ 3, 2, 1, 0 },
	};

	assert(first_letter  >= 0 && first_letter  <= 15);
	assert(second_letter >= 0 && second_letter <= 15);

	/* XXX - convert anything non-{A,C,G,T} to N */
	if (first_letter > BASE_T || second_letter > BASE_T)
		return (BASE_N);

	return (colourmat[first_letter][second_letter]);
}
