/*	$Id: util.c,v 1.25 2009/06/16 23:26:21 rumble Exp $	*/

#include <assert.h>
#include <inttypes.h>
#include <errno.h>
#include <math.h>
#include <stdarg.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <unistd.h>
#include <zlib.h>
#include <limits.h>
//#include <pthread.h>

#include <sys/time.h>
#include <sys/types.h>

#include "../common/fasta.h"
#include "../common/util.h"

//#define IMPORT_INLINES
//#include "../common/stats.c"


void
set_mode_from_argv(char **argv, shrimp_mode_t * shrimp_mode)
{

	if (strstr(argv[0], "-hs") != NULL || strstr(argv[0], "-HS") != NULL)
		*shrimp_mode = MODE_HELICOS_SPACE;
	else if (strstr(argv[0], "-cs") != NULL || strstr(argv[0], "-CS") != NULL)
		*shrimp_mode = MODE_COLOUR_SPACE;
	else
		*shrimp_mode = MODE_LETTER_SPACE;
}

const char *
get_mode_string(shrimp_mode_t shrimp_mode)
{

	switch (shrimp_mode) {
	case MODE_COLOUR_SPACE:
		return ("COLOUR SPACE (AB SOLiD)");
	case MODE_LETTER_SPACE:
		return ("LETTER SPACE (454,Illumina/Solexa,etc.)");
	case MODE_HELICOS_SPACE:
		return ("HELICOS SPACE (Helicos Single Molecule)");
	}

	return ("--bloody no idea!--");
}

uint64_t
gettimeinusecs()
{
	struct timeval tv;
	gettimeofday(&tv, NULL);
	return ((uint64_t)tv.tv_sec * 1000000 + tv.tv_usec);	
}

uint64_t
rdtsc()
{
	uint32_t lo, hi;

#ifdef __GNUC__
	__asm__ __volatile__ ("rdtsc" : "=a" (lo), "=d" (hi));
#else
	asm("rdtsc" : "=a" (lo), "=d" (hi));
#endif

	return (((uint64_t)hi << 32) | lo);
}

double
cpuhz()
{
	uint64_t before;
	struct timeval tv1, tv2;
	int diff;

	/* XXX - abusive, poor man's calc; needs good (2ms) clock granularity */
	gettimeofday(&tv1, NULL);
	before = rdtsc();
	do {
		gettimeofday(&tv2, NULL);

		diff = tv2.tv_usec - tv1.tv_usec;
		if (diff < 0)
			diff = 1000000 - tv1.tv_usec + tv2.tv_usec;
	} while (diff < 2000);

	return (((double)(rdtsc() - before) / diff) * 1.0e6);
}

u_int
strchrcnt(const char *str, const char c)
{
	int i;

	i = 0;
	while (*str != '\0') {
		if (*str++ == c)
			i++;
	}

	return (i);
}

bool
is_number(const char *str)
{

	while (*str != '\0')
		if (!isdigit((int)*str++))
			return (false);

	return (true);
}

bool
is_whitespace(const char *str)
{

	while (*str != '\0')
		if (!isspace((int)*str++))
			return (false);

	return (true);
}

void
xstat(const char *path, struct stat *sbp)
{
	
	if (stat(path, sbp) != 0) {
		fprintf(stderr, "error: failed to stat [%s]: %s\n", path,
		    strerror(errno));
		exit(1);
	}
}

void *
xmalloc(size_t size)
{
	void *ptr;

	ptr = malloc(size);
	if (ptr == NULL) {
		fprintf(stderr, "error: malloc failed: %s\n", strerror(errno));
		exit(1);
	}

	return (ptr);
}

void * xmalloc_m(size_t size, char const * error_msg)
{
  void * p;

  p = malloc(size);
  if (p == NULL) {
    fprintf(stderr, "error: malloc failed: %s [%s]\n", strerror(errno), error_msg);
    exit(1);
  }
  return p;
}


void *
xmalloc_c(size_t size, count_t * c)
{
  if (c != NULL)
    count_add(c, size);
  return xmalloc(size);
}

void *
xcalloc(size_t size)
{
  void *ptr;

  ptr = calloc(size, 1);
  if (ptr == NULL) {
    fprintf(stderr, "error: calloc failed: %s\n", strerror(errno));
    exit(1);
  }

  return ptr;
}

void * xcalloc_m(size_t size, char const * error_msg)
{
  void * p;

  p = calloc(size, 1);
  if (p == NULL) {
    fprintf(stderr, "error: calloc failed: %s [%s]\n", strerror(errno), error_msg);
    exit(1);
  }
  return p;
}

void *
xcalloc_c(size_t size, count_t * c)
{
  if (c != NULL)
    count_add(c, size);
  return xcalloc(size);
}

void *
xrealloc(void *ptr, size_t size)
{
  ptr = realloc(ptr, size);
  if (ptr == NULL) {
    fprintf(stderr, "error: realloc failed: %s\n", strerror(errno));
    exit(1);
  }

  return (ptr);
}

void *
xrealloc_c(void *ptr, size_t size, size_t old_size, count_t * c)
{
  if (c != NULL)
    count_add(c, (int64_t)size - (int64_t)old_size);
  return xrealloc(ptr, size);
}


char *
xstrdup(const char *str)
{
	char *dup;

	assert(str != NULL);

	dup = strdup(str);
	if (dup == NULL) {
		fprintf(stderr, "error: strdup failed: %s\n", strerror(errno));
		exit(1);
	}

	return (dup);
}

uint32_t
hash_string(const char *x)
{
	uint32_t hash = 0;

	while (*x != '\0')
		hash = 31 * hash + *x++;

	return (hash);
}

/* factorial using stirling's approximation after 20 */
double
ls_factorial(u_int n)
{
	const double fact[21] = {
		1.0,
		1.0,
		2.0,
		6.0,
		24.0,
		120.0,
		720.0,
		5040.0,
		40320.0,
		362880.0,
		3628800.0,
		39916800.0,
		479001600.0,
		6227020800.0,
		87178291200.0,
		1307674368000.0,
		20922789888000.0,
		355687428096000.0,
		6402373705728000.0,
		121645100408832000.0,
		2432902008176640000.0
	};
	double a, b;

	if (n <= 20)
		return log (fact[n]);

	a = log(sqrt(2 * M_PI * n));
	b = n * log(n / M_E);

	return (a + b);
}

/* choose in log space */
double
ls_choose(int64_t n, int64_t k)
{
	double a, b, c;

	if (k < 0 || k > n)
		return (0);

	a = ls_factorial(n);
	b = ls_factorial(k);
	c = ls_factorial(n - k);

	return (a - (b + c));
}

char *
trim_brackets(char *str)
{

	if (str[0] == '[')
		str++;
	if (str[strlen(str) - 1] == ']')
		str[strlen(str) - 1] = '\0';

	return (str);
}

/*
 * Prepend the low 4 bits of 'val' to the start of the bitfield in 'bf'.
 * 'entries' is the maximum number of 4-bit words to be stored in the
 * bitfield.
 */
void
bitfield_prepend(uint32_t *bf, uint32_t entries, uint32_t val)
{
	uint32_t tmp;
	u_int i;

	for (i = 0; i < BPTO32BW(entries); i++) {
		tmp = bf[i] >> 28;
		bf[i] <<= 4;
		bf[i] |= val;
		val = tmp;
	}

	bf[i - 1] &= (0xffffffff >> (32 - (4 * (entries % 8))));
}

/*
 * Insert the low 4 bits of 'val' into the bitfield in 'bf' at
 * 'index', where 'index' is count at 4-bit fields.
 */
void
bitfield_insert(uint32_t *bf, uint32_t index, uint32_t val)
{

	bitfield_append(bf, index, val);
}

/*
 * Append the low 4 bits of 'val' to the end of the bitfield in 'bf'.
 * 'entries' is the number of 4-bit words in 'bf' prior to the append.
 */
void
bitfield_append(uint32_t *bf, uint32_t entries, uint32_t val)
{
	uint32_t word;

	word = bf[entries / 8];
	word &= ~(0xf << (4 * (entries % 8)));
	word |= ((val & 0xf) << (4 * (entries % 8)));
	bf[entries / 8] = word;
}

void
progress_bar(FILE *out, uint64_t at, uint64_t of, uint incr)
{
	static int lastperc, beenhere;
	static char whirly = '\\';

	char progbuf[52 + 16];
	int perc, i, j, dec;

	if (at == 0 && of == 0) {
		beenhere = lastperc = 0;
		whirly = '\\';
		return;
	}

	perc = (at * 100 * incr) / of;

	if (beenhere && perc == lastperc)
		return;

	beenhere = 1;
	lastperc = perc;

	dec = perc % incr;
	perc /= incr;

	/* any excuse to have a whirly gig */
	switch (whirly) {
	case '|':
		whirly = '/';
		break;
	case '/':
		whirly = '-';
		break;
	case '-':
		whirly = '\\';
		break;
	case '\\':
		whirly = '|';
		break;
	}
	if (at >= of)
		whirly = '|';

	progbuf[25] = whirly;
		
	for (i = j = 0; i <= 100; i += 2) {
		if (j != 25) {
			if (i <= perc)
				progbuf[j++] = '=';
			else
				progbuf[j++] = ' ';
		} else {
			j++;
		}
	}
	memset(&progbuf[51], 0, 16);	/* XXX - valgrind */

	if (incr == 100)
		fprintf(out, "\rProgress: [%s] %3d.%02d%%", progbuf, perc, dec);
	else if (incr == 10)
		fprintf(out, "\rProgress: [%s] %3d.%d%%", progbuf, perc, dec);
	else
		fprintf(out, "\rProgress: [%s] %3d%%", progbuf, perc);
	
	fflush(out);
}

static inline uint32_t
swap_nibbles(uint32_t i)
{

	return (((i & 0xf0000000) >> 28) |
		((i & 0x0f000000) >> 20) |
		((i & 0x00f00000) >> 12) |
		((i & 0x000f0000) >>  4) |
		((i & 0x0000f000) <<  4) |
		((i & 0x00000f00) << 12) |
		((i & 0x000000f0) << 20) |
		((i & 0x0000000f) << 28));
}

/*
 * genome <- reverse_complement(genome)
 *
 * This is straightforward on purpose; it's suboptimal, but we won't be getting
 * any performance improvements by speeding it up.
 */
void
reverse_complement(uint32_t *g_ls, uint32_t *g_cs, uint32_t g_len, bool is_rna)
{
	uint32_t i, j, tmp, fudge, up, down;

	assert(g_len != 0);
	assert(g_ls != NULL);

	/* First, swap all words and the nibbles within and complement them. */
	for (i = 0, j = BPTO32BW(g_len) - 1; i <= j; i++, j--) {
		tmp = g_ls[i];
		g_ls[i] = swap_nibbles(g_ls[j]);
		if (i != j)
			g_ls[j] = swap_nibbles(tmp);

		g_ls[i]=(complement_base((g_ls[i] & 0x0000000f) >>  0, is_rna) <<  0) |
			(complement_base((g_ls[i] & 0x000000f0) >>  4, is_rna) <<  4) |
			(complement_base((g_ls[i] & 0x00000f00) >>  8, is_rna) <<  8) |
			(complement_base((g_ls[i] & 0x0000f000) >> 12, is_rna) << 12) |
			(complement_base((g_ls[i] & 0x000f0000) >> 16, is_rna) << 16) |
			(complement_base((g_ls[i] & 0x00f00000) >> 20, is_rna) << 20) |
			(complement_base((g_ls[i] & 0x0f000000) >> 24, is_rna) << 24) |
			(complement_base((g_ls[i] & 0xf0000000) >> 28, is_rna) << 28);

		/* Don't swap twice if we're at the same index */
		if (i == j) {
			// if i == j and i == 0, then we're done. Else j would
			// underflow in the for loop above with unsigned badness
			if (i == 0)
				break;
			else
				continue;
		}

		g_ls[j]=(complement_base((g_ls[j] & 0x0000000f) >>  0, is_rna) <<  0) |
			(complement_base((g_ls[j] & 0x000000f0) >>  4, is_rna) <<  4) |
			(complement_base((g_ls[j] & 0x00000f00) >>  8, is_rna) <<  8) |
			(complement_base((g_ls[j] & 0x0000f000) >> 12, is_rna) << 12) |
			(complement_base((g_ls[j] & 0x000f0000) >> 16, is_rna) << 16) |
			(complement_base((g_ls[j] & 0x00f00000) >> 20, is_rna) << 20) |
			(complement_base((g_ls[j] & 0x0f000000) >> 24, is_rna) << 24) |
			(complement_base((g_ls[j] & 0xf0000000) >> 28, is_rna) << 28);
	}

	/*
	 * If (g_len % 8) != 0, we need to shift all words down since
	 * the last word had some zeroed fields.
	 */
	fudge = 8 - (g_len % 8);
	if (fudge != 8) {
		down = fudge * 4;
		up = 32 - down;

		for (i = 0; i < BPTO32BW(g_len); i++) {
			if (i > 0)
				g_ls[i - 1] |= (g_ls[i] << up);
			g_ls[i] >>= down;
		}
	}

	/* If necessary, regenerate the colour-space genome */
	if (g_cs != NULL) {
		int base, lastbp;

		lastbp = BASE_T;
		for (i = 0; i < g_len; i++) {
			base = EXTRACT(g_ls, i);
			bitfield_append(g_cs, i, lstocs(lastbp, base, is_rna));
			lastbp = base;
		}
	}

}
/*
 * compute and return the reverse complement of the read in letterspace
 */
uint32_t *
reverse_complement_read_ls(uint32_t * read,uint32_t len, bool is_rna){
	uint32_t * rev_cmp;
	uint32_t i, j, fudge, up, down;

	if (len == 0) return NULL;

	rev_cmp = (uint32_t *)xmalloc(sizeof(uint32_t)*BPTO32BW(len));
	for (i = 0, j = BPTO32BW(len) - 1; i <= j; i++, j--) {
		rev_cmp[i] = swap_nibbles(read[j]);
		if (i != j)
			rev_cmp[j] = swap_nibbles(read[i]);

		rev_cmp[i]=(complement_base((rev_cmp[i] & 0x0000000f) >>  0, is_rna) <<  0) |
				(complement_base((rev_cmp[i] & 0x000000f0) >>  4, is_rna) <<  4) |
				(complement_base((rev_cmp[i] & 0x00000f00) >>  8, is_rna) <<  8) |
				(complement_base((rev_cmp[i] & 0x0000f000) >> 12, is_rna) << 12) |
				(complement_base((rev_cmp[i] & 0x000f0000) >> 16, is_rna) << 16) |
				(complement_base((rev_cmp[i] & 0x00f00000) >> 20, is_rna) << 20) |
				(complement_base((rev_cmp[i] & 0x0f000000) >> 24, is_rna) << 24) |
				(complement_base((rev_cmp[i] & 0xf0000000) >> 28, is_rna) << 28);

		/* Don't swap twice if we're at the same index */
		if (i == j) {
			// if i == j and i == 0, then we're done. Else j would
			// underflow in the for loop above with unsigned badness
			if (i == 0)
				break;
			else
				continue;
		}

		rev_cmp[j]=(complement_base((rev_cmp[j] & 0x0000000f) >>  0, is_rna) <<  0) |
				(complement_base((rev_cmp[j] & 0x000000f0) >>  4, is_rna) <<  4) |
				(complement_base((rev_cmp[j] & 0x00000f00) >>  8, is_rna) <<  8) |
				(complement_base((rev_cmp[j] & 0x0000f000) >> 12, is_rna) << 12) |
				(complement_base((rev_cmp[j] & 0x000f0000) >> 16, is_rna) << 16) |
				(complement_base((rev_cmp[j] & 0x00f00000) >> 20, is_rna) << 20) |
				(complement_base((rev_cmp[j] & 0x0f000000) >> 24, is_rna) << 24) |
				(complement_base((rev_cmp[j] & 0xf0000000) >> 28, is_rna) << 28);
	}
	/*
	 * If (g_len % 8) != 0, we need to shift all words down since
	 * the last word had some zeroed fields.
	 */
	fudge = 8 - (len % 8);
	if (fudge != 8) {
		down = fudge * 4;
		up = 32 - down;

		for (i = 0; i < BPTO32BW(len); i++) {
			if (i > 0)
				rev_cmp[i - 1] |= (rev_cmp[i] << up);
			rev_cmp[i] >>= down;
		}
	}
	return rev_cmp;

}

uint32_t *
reverse_complement_read_cs(uint32_t * read, int8_t initbp, int8_t initbp_rc, uint32_t len, bool is_rna){
	uint32_t * read_rc;
	uint i;
	int8_t base;

	assert(len > 0);
	read_rc = (uint32_t *)xmalloc(sizeof(uint32_t)*BPTO32BW(len));

	base = (int8_t)cstols(initbp, EXTRACT(read, 0), is_rna);
	for (i = 1; i < len; i++){
	  base = (int8_t)cstols(base, EXTRACT(read,i), is_rna);
		bitfield_insert(read_rc, len - i, EXTRACT(read,i));
	}
	bitfield_insert(read_rc, 0, lstocs(base, complement_base(initbp_rc, is_rna), is_rna));

	return read_rc;
}

void
reverse_complement_read_ls_text(char * read, char *ret){
	int len = strlen(read);
	int c;
	for (c = 0; c < len; c ++){
		char base = read[len - 1 - c];
		if (base == 'g' || base == 'G'){
			base = 'c';
		} else if (base == 'a' || base == 'A'){
			base = 't';
		} else if (base == 't' || base == 'T'){
			base = 'a';
		} else if (base == 'c' || base == 'C'){
			base = 'g';
		}
		ret[c] = base;
	}
	ret[len] = '\0';
}

/*
 * Given a path, if it's a regular file (or symlink to one), call fh on it.
 * If it's a directory, call fh on all regular files within it (or symlinks to
 * regular files).
 *
 * Returns the number of files fh was called on.
 */
uint64_t
file_iterator(char *path, void (*fh)(char *, struct stat *, void *),
    void *arg)
{
	char fpath[2048];
	struct stat sb;
	DIR *dp;
	struct dirent *de;
	uint64_t files;

	/* is a regular file... */
	xstat(path, &sb);
	if (S_ISREG(sb.st_mode)) {
		fh(path, &sb, arg);
		return (1);
	}

	/* is (hopefully) a directory... */
	dp = opendir(path);
	if (dp == NULL) {
		fprintf(stderr, "error: failed to open directory [%s]: %s\n",
		    path, strerror(errno));
		exit(1);
	}

	files = 0;
	while (1) {
		de = readdir(dp);
		if (de == NULL)
			break;

		strcpy(fpath, path);
		if (fpath[strlen(path) - 1] != '/')
			strcat(fpath, "/");
		strcat(fpath, de->d_name);
		xstat(fpath, &sb);

#if defined(DT_REG) && defined(DT_LNK)
		if (de->d_type != DT_REG && de->d_type != DT_LNK)
			continue;
#else
		if (!S_ISREG(sb.st_mode) && !S_ISLNK(sb.st_mode))
			continue;
#endif

		/* ensure it's a regular file or link to one */
		if (S_ISREG(sb.st_mode)) {
			fh(fpath, &sb, arg);
			files++;
		} else {
			fprintf(stderr, "warning: [%s] is neither a regular "
			    "file, nor a link to one; skipping...", fpath);
			continue;
		}
	}

	closedir(dp);

	return (files);
}

uint64_t
file_iterator_n(char **paths, int npaths,
    void (*fh)(char *, struct stat *, void *), void *arg)
{
	uint64_t files;
	int i;

	for (i = files = 0; i < npaths; i++)
		files += file_iterator(paths[i], fh, arg);

	return (files);
}

char const *
get_compiler()
{

#if defined(__GNUC__)
	if (strstr(__VERSION__, "Intel(R)"))
		return ("ICC " __VERSION__);
	else
		return ("GCC " __VERSION__);
#elif defined(__SUNPRO_C)
	return ("Sun Pro C");
#elif defined(__SUNPRO_CC)
	return ("Sun Pro C++");
#elif defined(__cplusplus)
	return ("unknown C++");
#else
	return ("unknown C");
#endif
}

/* reverse the string `str' in place */
char *
strrev(char *str)
{
	char c;
	int i, j;

	j = strlen(str) - 1;
	for (i = 0; i < j; i++, j--) {
		c = str[j];
		str[j] = str[i];
		str[i] = c;
	}

	return (str);
}

/* trim whitespace in `str' in place at beginning and end */
char *
strtrim(char *str)
{
	char *ret;

        assert(str != NULL);

        while (isspace((int)*str) && *str != '\0')
                str++;

        ret = str;

	if (*str != '\0') {
		while (*str != '\0')
			str++;

		str--;
		
		while (isspace((int)*str))
			str--;
		str++;
	}

        *str = '\0';

        return (ret);
}

strbuf_t
strbuf_create()
{
	strbuf_t sbp;

	sbp = (strbuf_t)xmalloc(sizeof(*sbp));
	memset(sbp, 0, sizeof(*sbp));
	sbp->string_alloced = 4096;
	sbp->string = (char *)xmalloc(4096);

	return (sbp);
}

char *
strbuf_string(strbuf_t sbp, int *length)
{

	assert(sbp->string_length == strlen(sbp->string));
	if (length != NULL)
		*length = sbp->string_length;
	return (xstrdup(sbp->string));
}

void
strbuf_append(strbuf_t sbp, char const *fmt, ...)
{
	va_list ap;
	int bytes;

	assert(sbp->string_length < sbp->string_alloced);
	if (sbp->string_alloced == sbp->string_length) {
		sbp->string_alloced += 4096;
		sbp->string = (char *)xrealloc(sbp->string, sbp->string_alloced);
	}

	va_start(ap, fmt);
	do {
		bytes = vsnprintf(&sbp->string[sbp->string_length],
		    sbp->string_alloced - sbp->string_length, fmt, ap);

		/* wasn't enough space. resize and try again */
		if (sbp->string_length + bytes >= sbp->string_alloced) {
			sbp->string_alloced += 4096;
			sbp->string = (char *)xrealloc(sbp->string, sbp->string_alloced);
		}
	} while (sbp->string_length + bytes >= sbp->string_alloced);
	va_end(ap);

	sbp->string_length += bytes;
	assert(sbp->string_length < sbp->string_alloced);
}

void
strbuf_destroy(strbuf_t sbp)
{

	free(sbp->string);
	free(sbp);
}

/*
 * Fast and loose replacement for gzgets, for use when one is just
 * sequentially calling gzgets on the whole file. Note that this is
 * _not_ reentrant!
 */
char *
fast_gzgets(gzFile gz, char *buf, int len)
{
	static char *save_buf;
	static int   save_len;
	static int   save_bytes;
	static int   save_skip;

	int ret, i, total;

	assert(len >= 0);

	if (len == 0)
		goto out;

	if (save_len < len) {
		save_len = len;
		save_buf = (char *)xrealloc(save_buf, save_len);
	}

	assert(save_bytes >= 0);
	assert(save_bytes < save_len);

	ret = -1;
	total = 0;
	do {
		bool nl = false;

		assert(save_skip < save_len);

		if (save_bytes == 0) {
			ret = gzread(gz, save_buf, save_len - 1);
			if (ret < 0)
				break;
			save_skip = 0;
			save_bytes = ret;
		}

		for (i = 0; i < save_bytes && (total + i) < (save_len - 1); i++) {
			if ((buf[total + i] = save_buf[save_skip + i]) == '\n') {
				nl  = true;
				i++;
				break;
			}
		}

		total += i;
		buf[total] = '\0';

		assert(save_bytes >= 0);
		assert(total < save_len);

		if (nl) {
			assert(i > 0);
			assert(i <= save_bytes);
			assert(save_buf[save_skip + i - 1] == '\n');
			save_bytes -= i;
			save_skip += i;
			return (buf);
		} else if (total == (save_len - 1)) {
			if (i == total) { 
				save_bytes = save_skip = 0;
			} else {
				assert(i < total);
				save_bytes -= i;
				save_skip += i;
			}
			return (buf);
		} else if (save_bytes == 0) {
			if (total == 0)
				break;
			assert(ret == 0);
			save_bytes = save_skip = 0;
			return (buf);
		} else if (i == save_bytes) {
			assert(i > 0);
			save_bytes -= i;
			save_skip += i;
		} else {
			assert(0);
		}
	} while (true);

 out:
	if (save_buf != NULL)
		free(save_buf);
	save_buf = NULL;
	save_len = save_bytes = save_skip = 0;

	return (NULL);
}
/*
 * Fast and loose replacement for gzgets, for use when one is just
 * sequentially calling gzgets on the whole file. Note that this is
 * _not_ reentrant!
 */

char *
fast_gzgets_safe(fasta_t f)
{
	
	int ret, i, total;

	assert(sizeof(f->buffer) >= 0);

	if (sizeof(f->buffer) == 0) {
		fprintf(stderr,"Buffer has size 0!!!\n");
		exit(1);	
	}
	if (f->save_len < (int)sizeof(f->buffer)) {
		f->save_len = sizeof(f->buffer);
		f->save_buf = (char *)xrealloc(f->save_buf, f->save_len);
	}

	assert(f->save_bytes >= 0);
	assert(f->save_bytes < f->save_len);

	ret = -1;
	total = 0;
	do {
		bool nl = false;

		assert(f->save_skip < f->save_len);

		if (f->save_bytes == 0) {
			ret = gzread(f->fp, f->save_buf, f->save_len - 1);
			if (ret < 0)
				break;
			f->save_skip = 0;
			f->save_bytes = ret;
		}
		int offset=0;
		for (i = 0; i < f->save_bytes && (total + i) < (f->save_len - 1); i++) {
			if ((f->buffer[total + i - offset] = f->save_buf[f->save_skip + i]) == '\n' ) {
				if (f->header && f->buffer[0]=='#') {
					f->buffer[total+i-offset]='\0';
					offset=i+1;
					total=0;
					continue;
				}
				f->header=false;
				nl  = true;
				i++;
				break;
			}
		}

		total += i;
		f->buffer[total] = '\0';

		assert(f->save_bytes >= 0);
		assert(total < f->save_len);

		if (nl) {
			assert(i > 0);
			assert(i <= f->save_bytes);
			assert(f->save_buf[f->save_skip + i - 1] == '\n');
			f->save_bytes -= i;
			f->save_skip += i;
			return (f->buffer);
		} else if (total == (f->save_len - 1)) {
			if (i == total) { 
				f->save_bytes = f->save_skip = 0;
			} else {
				assert(i < total);
				f->save_bytes -= i;
				f->save_skip += i;
			}
			return (f->buffer);
		} else if (f->save_bytes == 0) {
			if (total == 0)
				break;
			assert(ret == 0);
			f->save_bytes = f->save_skip = 0;
			return (f->buffer);
		} else if (i == f->save_bytes) {
			assert(i > 0);
			f->save_bytes -= i;
			f->save_skip += i;
		} else {
			assert(0);
		}
	} while (true);

	if (f->save_buf != NULL)
		free(f->save_buf);
	f->save_buf = NULL;
	f->save_len = f->save_bytes = f->save_skip = 0;

	return (NULL);
}

#if 0
/*
 * Try to assure our fast_gzgets works as fgets and gzgets.
 * I think gzgets doesn't fully conform to fgets when TESTBUFLEN = 0.
 */
void
fast_gztest(int nfiles, char **files)
{
	int i;
#define TESTBUFLEN 9731
	for (i = 0; i < nfiles; i++) {
		FILE *fp;
		gzFile gz1, gz2;
		char *buf  = (char *)xmalloc(TESTBUFLEN);
		char *buf2 = (char *)xmalloc(TESTBUFLEN);
		char *buf3 = (char *)xmalloc(TESTBUFLEN);

		printf("=== %s ===\n", files[i]);

		fp = fopen(files[i], "r");
		gz1 = gzopen(files[i], "r");
		gz2 = gzopen(files[i], "r");
		if (fp == NULL && gz1 == NULL && gz2 == NULL) {
			printf("skipping [%s]\n", files[i]);
			continue;
		}
		assert(fp != NULL && gz1 != NULL && gz2 != NULL);

		while (gzgets(gz1, buf, TESTBUFLEN) != NULL) {
			char *ret = fast_gzgets(gz2, buf2, TESTBUFLEN);
			char *ret2 = fgets(buf3, TESTBUFLEN, fp);
			assert(ret != NULL);
			assert(ret2 != NULL);
			if (strcmp(buf, buf2) != 0) {
				printf("buf [%s]\n !=\n buf2 [%s]\n", buf, buf2);
				assert(0);
			}
			if (strcmp(buf, buf3) != 0) {
				printf("buf [%s]\n !=\n buf3 [%s]\n", buf, buf2);
				assert(0);
			}
		}
		assert(fgets(buf3, sizeof(buf3), fp) == NULL);
		assert(fast_gzgets(gz2, buf2, sizeof(buf2)) == NULL);

		fclose(fp);
		gzclose(gz1);
		gzclose(gz2);
	}
}
#endif

void xgzread(gzFile fp, voidp buf, size_t len) {
        size_t total = 0;
        while (total < len) {
                int res;
                int to_read = (len - total <= (size_t)INT_MAX? (int)(len - total) : INT_MAX);
                res = gzread(fp, (char *)buf+total, to_read);
                if (res <= 0) {
                        fprintf(stderr, "error: gzread returned %u\n", res);
                        exit(0);
                }
                total += (size_t)res;
        }

}

void xgzwrite(gzFile fp, voidp buf, unsigned len){
	uint total = 0;
	while (total < len){
		uint res;
		res = gzwrite(fp,(char *)buf+total,len-total);
		if (res <= 0){
			fprintf(stderr,"error: gzwrite returned %u\n",res);
			exit(0);
		}
		total += res;
	}
}

/*
 * Return a string on the stack corresponding to an unsigned integer that also
 * features commas. E.g.: int 1000 yields "1,000".
 *
 * There is a pool of sequentially allocated buffers returned, so this should
 * be safe to use multiple times in function arguments.
 */
char *
comma_integer(uint64_t val)
{
	static char rets[50][32];	// no malloc, allow uses in fn args, etc
	static int col = 0;

	char *ret = rets[(col++ % (sizeof(rets) / sizeof(rets[0])))];
	char str[sizeof(rets[0])];
	int skip, i, j;

	memset(str, 0, sizeof(str));	// XXX - shut up, valgrind
	snprintf(str, sizeof(str), "%" PRIu64, val);

	skip = 3 - (strlen(str) % 3);
	for (i = j = 0; str[i] != '\0'; i++) {
		if ((i + skip) % 3 == 0 && i != 0)
			ret[j++] = ',';
		ret[j++] = str[i];
	}
	ret[j] = '\0';

	return (ret);
}
enum mode {start,match,mismatch,ref_gap,read_gap};

static int
finish_mode(enum mode mode,int count,char *res){
	if (mode == match || mode == mismatch) return sprintf(res,"%iM",count);
	else if (mode == ref_gap) return sprintf(res,"%iI",count);
	else if (mode == read_gap) return sprintf(res,"%iD",count);
	else return 0;
}

void
edit2cigar(char * edit,uint16_t read_start,uint16_t read_end,uint16_t read_length,char *res){
	char * current=edit;
	enum mode mode = start;
	int count = 0;
	int last_count = 0;

	if(read_start != 0){
	  res += sprintf(res,"%uS",(unsigned int)read_start);
	}

	while(*current != '\0'){
		if (*current <= '9' && *current >= '0'){
			if(mode != match){
				if (mode == mismatch){
					last_count += count;
				} else {
					res += finish_mode(mode,count + last_count, res);
					last_count = 0;
				}
				count = 0;
			}
			int t;
			t = *current - '0';
			assert(t >= 0 && t <= 9);
			count = count*10 + t;
			mode = match;
		} else if(*current == '('){
			res += finish_mode(mode, count + last_count, res);
			count = 0;
			last_count = 0;
			mode = ref_gap;
		} else if(*current == 'G' || *current == 'A' || *current == 'T' || *current == 'C'){
			if (mode == ref_gap){
				count++;
			} else{
				if (mode != mismatch){
					if (mode == match){
						last_count += count;
					} else{
						res += finish_mode(mode,count + last_count,res);
						last_count = 0;
					}
					count = 0;
				}
				count++;
				mode = mismatch;

			}

		} else if (*current == ')'){
			res += finish_mode(mode, count + last_count, res);
			count = 0;
			last_count = 0;
			mode = start;
		} else if (*current == '-'){
			if (mode != read_gap){
				res += finish_mode(mode, count + last_count, res);
				count = 0;
				last_count = 0;
			}
			count ++;
			mode = read_gap;
		} else if (*current == 'x'){
			last_count += count;
			count = 0;
		}
		current++;
	}
	res += finish_mode(mode,count+ last_count,res);
	if(read_end + 1 != read_length){
	  sprintf(res,"%uS",(unsigned int)(read_length - read_end - 1));
	}
}


// pre: array is (q)sorted
size_t removedups(void * a, size_t n, size_t sz, int (*cmp)(void const *, void const *))
{
  char * p = (char *)a;
  size_t m = 0;
  size_t i, j;
  i = 0;
  while (i < n) {
    j = i+1;
    while (j < n && !cmp(p+i*sz, p+j*sz)) j++;
    if (m < i)
      memcpy(p+m*sz, p+i*sz, sz);
    m++;
    i = j;
  }
  return m;
}


void
crash(int exit_code, int display_errno, char const * msg, ...)
{
  va_list fmtargs;
  char new_msg[strlen(msg) + 1000];

  sprintf(new_msg, "error: %s%s%s\n", msg,
	  display_errno? ": " : "",
	  display_errno? strerror(errno) : "");

  va_start(fmtargs, msg);
  vfprintf(stderr, new_msg, fmtargs);
  va_end(fmtargs);

  exit(exit_code);
}


void
logit(int display_errno, char const * msg, ...)
{
  va_list fmtargs;
  char new_msg[strlen(msg) + 1000];

  sprintf(new_msg, "note: %s%s%s\n", msg,
	  display_errno? ": " : "",
	  display_errno? strerror(errno) : "");

  va_start(fmtargs, msg);
  vfprintf(stderr, new_msg, fmtargs);
  va_end(fmtargs);
}


long long
nchoosek(int n, int k)
{
  long long res = 1;
  int i;
  for (i = 0; i < k; i++)
    res *= (n - i);
  for (i = 0; i < k; i++)
    res /= (i + 1);
  return res;
}


double
log_nchoosek(int n, int k)
{
  double res = 0.0;
  int i;
  for (i = 0; i < k; i++)
    res += log(n - i) - log(i + 1);
  return res;
}

void
cat(FILE * src, FILE * dest)
{
  size_t buffer_size = 2046;
  char buffer[buffer_size];
  size_t read;
  bool ends_in_newline = true;
  while ((read = fread(buffer, 1, buffer_size - 1, src))) {
    buffer[read] = '\0';
    fprintf(dest, "%s", buffer);
    if (buffer[read - 1] == '\n') {
      ends_in_newline = true;
    } else {
      ends_in_newline = false;
    }
  }
  if (!ends_in_newline) {
    fprintf(dest,"\n");
  }
}
