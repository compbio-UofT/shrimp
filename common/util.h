/*	$Id: util.h,v 1.22 2009/06/16 23:26:21 rumble Exp $	*/
#ifndef _UTIL_H
#define _UTIL_H

/*
 * Force use of C linking for util.c, even if using g++.
 */
#ifdef __cplusplus
//extern "C" {
#endif


#include <assert.h>
#include <ctype.h>
#include <dirent.h>
#include <unistd.h>
#include <zlib.h>

#include <sys/types.h>
#include <sys/stat.h>
#include "../gmapper/gmapper-definitions.h"
#include "../common/fasta.h"
#include "../common/stats.h"
#include "../common/hash.h"
#include "../common/my-alloc.h"

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

/*
 * If window_len, sw_vect_threshold, sw_full_threshold are absolute values,
 * we'll set them negative to distinguish.
 */
#define IS_ABSOLUTE(x)  ((x) < 0)

  //static inline double abs_or_pct(double x, double base) {
  //return IS_ABSOLUTE(x) ? -x : base * (x / 100.0);
  //}
#define abs_or_pct(x, base) (IS_ABSOLUTE(x) ? -(x) : (base) * ((x) / 100.0))

#define PASS_THRESHOLD_1000(abs_val, pct_val_1000, thres) \
  (IS_ABSOLUTE((thres))? (abs_val) >= (int)(-(thres)) : (pct_val_1000) >= (int)((thres) * 1000))

#define KMER_TO_MAPIDX(kmer, sn) (Hflag? kmer_to_mapidx_hash((kmer), (sn)) : kmer_to_mapidx_orig((kmer), (sn)))


struct _strbuf_t {
	char   *string;
	u_int	string_length;
	u_int	string_alloced;
};
typedef struct _strbuf_t * strbuf_t;

void		set_mode_from_argv(char **, shrimp_mode_t *);
const char     *get_mode_string(shrimp_mode_t);
uint64_t	gettimeinusecs(void);
uint64_t	rdtsc(void);
double		cpuhz(void);
u_int		strchrcnt(const char *, const char);
bool		is_number(const char *);
bool		is_whitespace(const char *);
void		xstat(const char *, struct stat *);
void *		xmalloc(size_t);
void *		xmalloc_m(size_t, char const *);
void *		xmalloc_c(size_t, count_t *);
void *		xcalloc(size_t);
void *		xcalloc_m(size_t, char const *);
void *		xcalloc_c(size_t, count_t *);
void	       *xrealloc(void *, size_t);
void	       *xrealloc_c(void *, size_t, size_t, count_t *);
char	       *xstrdup(const char *);
uint32_t	hash_string(const char *);
double		ls_factorial(u_int);
double		ls_choose(int64_t, int64_t);
char	       *trim_brackets(char *);
void		bitfield_prepend(uint32_t *, uint32_t, uint32_t);
void		bitfield_insert(uint32_t *, uint32_t, uint32_t);
void		bitfield_append(uint32_t *, uint32_t, uint32_t);
void		progress_bar(FILE *, uint64_t, uint64_t, uint);
void		reverse_complement(uint32_t *, uint32_t *, uint32_t, bool);
uint32_t *	reverse_complement_read_cs(uint32_t *, int8_t, int8_t, uint32_t, bool);
uint32_t *	reverse_complement_read_ls(uint32_t *, uint32_t, bool);
void		reverse_complement_read_ls_text(char *, char *);
uint64_t	file_iterator(char *, void (*)(char *, struct stat *, void *), void *);
uint64_t	file_iterator_n(char **, int, void (*)(char *, struct stat *, void *), void *);
char const *	get_compiler(void);
char *		strrev(char *);
char	*	strtrim(char *);
strbuf_t	strbuf_create(void);
char *		strbuf_string(strbuf_t, int *);
void		strbuf_append(strbuf_t, char const *, ...);
void		strbuf_destroy(strbuf_t);
char *		fast_gzgets(gzFile, char*, int);
char *		fast_gzgets_safe(fasta_t);
char *		comma_integer(uint64_t);
void		xgzwrite(gzFile, voidp, unsigned);
void		xgzread(gzFile, voidp, size_t);
void		edit2cigar(char *, uint16_t, uint16_t, uint16_t, char *);
size_t		removedups(void *, size_t, size_t, int (*)(void const *, void const *));
void		crash(int, int, char const *, ...);
void		logit(int, char const *, ...);
long long	nchoosek(int, int);
double		log_nchoosek(int, int);
void		cat(FILE *, FILE *);


/* for optarg (and to shut up icc) */
//extern char *optarg;
//extern int   optind;

static inline int
complement_base(int base, bool is_rna)
{
	static const int cmpl[16] = {
		BASE_T,		/* A -> T (or U) */
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

	return ((is_rna && cmpl[base] == BASE_T) ? BASE_U : cmpl[base]);
}

/*
 * Given the first letter of a pair corresponding to a colour, and the colour
 * itself, obtain the next letter.
 */
static inline int
cstols(int first_letter, int colour, bool is_rna)
{
	//TODO NOT SURE IF THIS IS CORRECT WAY TO HANDLE
	if (first_letter==BASE_N || !(colour>=0 && colour<=3)) {
		return BASE_N;
	}
	assert((first_letter >= 0 && first_letter <= 3 && !is_rna) ||
	       (first_letter >= 0 && first_letter <= 4 &&  is_rna));
	assert(colour >= 0 && colour <= 3);
	int ret;

	if (is_rna && first_letter == BASE_U)
		first_letter = BASE_T;

	if ((first_letter % 2) == 0)
		ret = ((4 + first_letter + colour) % 4);
	else
		ret = ((4 + first_letter - colour) % 4);

	if (is_rna && ret == BASE_T)
		return (BASE_U);
	return (ret);
}

static inline int
lstocs(int first_letter, int second_letter, bool is_rna)
{
	const int colourmat[4][4] = {
		{ 0, 1, 2, 3 },
		{ 1, 0, 3, 2 },
		{ 2, 3, 0, 1 },
		{ 3, 2, 1, 0 },
	};

	assert(first_letter  >= 0 && first_letter  <= 15);
	assert(second_letter >= 0 && second_letter <= 15);

	if (is_rna) {
		if (first_letter == BASE_U)
			first_letter = BASE_T;
		if (second_letter == BASE_U)
			second_letter = BASE_T;
	}

	/* XXX - convert anything non-{A,C,G,T} to N */
	if (first_letter > BASE_T || second_letter > BASE_T)
		return (BASE_N);

	return (colourmat[first_letter][second_letter]);
}

static inline uint
ceil_div(uint a, uint b) {
  assert(b > 0);

  if (a == 0) return 0;
  return ((a - 1) / b) + 1;
}

/*
 * Compute a hash value for a genomic region.
 */
static inline uint32_t
hash_genome_window(uint32_t * genome, uint goff, uint glen) {
  static uint const bases_per_buffer = 16;

  uint i, j;
  uint base;
  uint32_t key = 0;
  uint32_t buffer;

  for (i = 0; i < ceil_div(glen, bases_per_buffer); i++) {
    buffer = 0;
    for (j = 0; j < bases_per_buffer && i * bases_per_buffer + j < glen; j++) {
      base = EXTRACT(genome, goff + i * bases_per_buffer + j);
      buffer <<= 2;
      buffer |= (base & 0x03);
    }
    hash_accumulate(&key, buffer);
  }
  hash_finalize(&key);

  return key;
}

/* compute base^power */
static inline size_t
power(size_t base, size_t exp)
{
  size_t result = 1;

  while (exp > 0) {
    if ((exp % 2) == 1)
      result *= base;
    base *= base;
    exp /= 2;
  }

  return result;
}

/* compute 4^power */
static inline llint
power4(int exp) {
  return (llint)1 << (2 * exp);
}


/* standard quality value from probability of error */
static inline int
qv_from_pr_err(double pr_err)
{
  if (pr_err > .99999999)
    return 0;
  else if (pr_err < 1E-25)
    return 250;
  else
    return (int)(-10.0 * log(pr_err) / log(10.0));
}

static inline int
qv_from_pr_corr(double pr_corr)
{
  return qv_from_pr_err(1 - pr_corr);
}

static inline double
pr_err_from_qv(int qv)
{
  if (qv <= 0)
    return .99999999;
  else if (qv >= 250)
    return 1E-25;
  else
    return pow(10.0, -(double)qv/10.0);
}


static inline int
double_to_neglog(double x, int shift = 1000)
{
  return (int)((double)shift * -log(x));
}


static inline double
neglog_to_double(int y, int shift = 1000)
{
  return exp(-(double)y / (double)shift);
}


static inline double
normal_cdf(double x, double mean, double stddev)
{
  double y = (x - mean) / stddev;
  if (y < 0) y = -y;
  double b0 = 0.2316419;
  double b1 = 0.319381530;
  double b2 = -0.356563782;
  double b3 = 1.781477937;
  double b4 = -1.821255978;
  double b5 = 1.330274429;
  double pi = 3.141592653589;
  double t = 1.0 / (1.0 + b0 * y);
  double res = (exp(- y * y / 2) / sqrt(2.0 * pi)) * ((((b5 * t + b4) * t + b3) * t + b2) * t + b1) * t;
  if (x > mean) res = 1 - res;
  return res;
}



#ifdef __cplusplus
//} /* extern "C" */
#endif

#endif
