/*	$Id$	*/

#include <assert.h>
#include <ctype.h>
#include <errno.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <unistd.h>

#include <mmintrin.h>	/* MMX */
#include <xmmintrin.h>	/* SSE */
#include <emmintrin.h>	/* SSE2 */
//#include <pmmintrin.h>/* SSE3 */

#include <sys/time.h>

#include "../common/util.h"

#include "sw-vector.h"

static int	initialised;
static int8_t  *db, *db_ls, *qr;
static int	dblen, qrlen;
static int16_t *nogap, *b_gap;
static int	gap_open, gap_ext;
static int	match, mismatch;
static uint64_t swticks, swcells, swinvocs;
static int	use_colours;

/*
 * Calculate the Smith-Waterman score.
 *
 * This is basically an SSE2 version of Wozniak's vectored implementation, but
 * without a score table. Further, we assume a fixed database and query size,
 * so *nogap and *b_gap must be pre-allocated (the malloc overhead for very
 * small scans is _huge_).
 *
 * NOTE THE FOLLOWING:
 *
 *	1) seqA must be padded with 7 bytes at the beginning and end. The first
 *	   element of seqA should be the first pad byte.
 *
 *	2) seqB must be padded with bytes on the end up to mod 8 characters.
 *	   The first element of seqB should be (of course) the first character.
 *
 *	3) seqA and seqB's padding _must_ be different, otherwise our logic will
 *	   consider the padding as matches!
 */
static int
vect_sw(int8_t *seqA, int lena, int8_t *seqB, int lenb,
    int8_t *ls_seqA, int initbp)
{
	int i, j, score = 0;
	__m128i v_score, v_zero, v_gap_ext, v_gap_open_ext;
	__m128i v_match, v_mismatch;
	__m128i v_a_gap, v_b_gap, v_nogap;
	__m128i v_last_nogap, v_prev_nogap, v_seq_a, v_seq_b;
	__m128i v_tmp;
	int16_t w[8];

	/* shut up icc */
	(void)ls_seqA;
	(void)initbp;

#define SET16(a, e7, e6, e5, e4, e3, e2, e1, e0)      \
	_mm_set_epi16((int16_t)a[e7], (int16_t)a[e6], \
		      (int16_t)a[e5], (int16_t)a[e4], \
		      (int16_t)a[e3], (int16_t)a[e2], \
		      (int16_t)a[e1], (int16_t)a[e0])

	v_score		= _mm_setzero_si128();
	v_zero		= _mm_setzero_si128();
	v_match		= SET16((&match), 0, 0, 0, 0, 0, 0, 0, 0);
        v_mismatch	= SET16((&mismatch), 0, 0, 0, 0, 0, 0, 0, 0);
	v_gap_ext	= SET16((&gap_ext), 0, 0, 0, 0, 0, 0, 0, 0);
	v_gap_open_ext	= SET16((&gap_open), 0, 0, 0, 0, 0, 0, 0, 0);
	v_gap_open_ext	= _mm_add_epi16(v_gap_open_ext, v_gap_ext);

        for (i = 0; i < lena + 14; i++) {
                nogap[i] = 0;
                b_gap[i] = (int16_t)-gap_open;
        }

	/*
	 * When using colour space reads, we must handle the first row
	 * specially. This is because the read will begin with some marker
	 * base, which will affect matching against the genome.
	 *
	 * For 25mer reads, this actually makes things faster, because our
	 * vectorised portion becomes evenly divisible by 8 again. Yey.
	 */
	if (use_colours) {
		int a_gap, prev_nogap, last_nogap;

		a_gap = -gap_open;
		last_nogap = prev_nogap = 0;
		for (i = 7; i < (lena + 7); i++) {
			int a, ms;

			a_gap = MAX((last_nogap - gap_open - gap_ext),
			    (a_gap - gap_ext));
			b_gap[i] =(uint16_t)MAX((nogap[i] - gap_open - gap_ext),
			    (b_gap[i] - gap_ext));

			a = lstocs(ls_seqA[i], initbp);
			ms = (a == seqB[0]) ? match : mismatch;

			last_nogap = MAX((prev_nogap + ms), 0);
			last_nogap = MAX(last_nogap, a_gap);
			last_nogap = MAX(last_nogap, b_gap[i]);
			prev_nogap = nogap[i];
			nogap[i] = (uint16_t)last_nogap;
			score = MAX(score, last_nogap);
		}

		v_score = SET16((&score), 0, 0, 0, 0, 0, 0, 0, 0);
		score = 0;
		seqB++;
		lenb--;

		assert(lenb != 0);
	}

	for (i = 0; i < (lenb + 7)/8; i++) {
		int k = i * 8;

		v_b_gap = SET16(b_gap, 6, 6, 5, 4, 3, 2, 1, 0);
		v_nogap = SET16(nogap, 6, 6, 5, 4, 3, 2, 1, 0);
		v_seq_a = SET16(seqA, 0, 0, 1, 2, 3, 4, 5, 6);
		v_seq_b = SET16(seqB, k+7, k+6, k+5, k+4, k+3, k+2, k+1, k+0);

		v_a_gap = v_gap_ext;
		v_a_gap = _mm_sub_epi16(v_a_gap, v_gap_open_ext);

		v_last_nogap = _mm_setzero_si128();
		v_prev_nogap = _mm_setzero_si128();

		for (j = 0; j < (lena + 7); j++) {
			v_b_gap = _mm_slli_si128(v_b_gap, 2);
			v_b_gap = _mm_insert_epi16(v_b_gap, b_gap[j+7], 0);

			v_nogap = _mm_slli_si128(v_nogap, 2);
			v_nogap = _mm_insert_epi16(v_nogap, nogap[j+7], 0);

			v_seq_a = _mm_slli_si128(v_seq_a, 2);
			v_seq_a = _mm_insert_epi16(v_seq_a, seqA[j+7], 0);

			v_tmp = _mm_sub_epi16(v_last_nogap, v_gap_open_ext);
			v_a_gap = _mm_sub_epi16(v_a_gap, v_gap_ext);
			v_a_gap = _mm_max_epi16(v_a_gap, v_tmp);

			v_tmp = _mm_sub_epi16(v_nogap, v_gap_open_ext);
			v_b_gap = _mm_sub_epi16(v_b_gap, v_gap_ext);
			v_b_gap = _mm_max_epi16(v_b_gap, v_tmp);

			/* compute the score (v_last_nogap is a tmp variable) */
			v_last_nogap = _mm_cmpeq_epi16(v_seq_a, v_seq_b);
			v_tmp = _mm_and_si128(v_last_nogap, v_match);
			v_last_nogap = _mm_cmpeq_epi16(v_last_nogap, v_zero);
			v_last_nogap = _mm_and_si128(v_last_nogap, v_mismatch);
			v_tmp = _mm_or_si128(v_tmp, v_last_nogap);

			v_last_nogap = _mm_add_epi16(v_prev_nogap, v_tmp);
			v_last_nogap = _mm_max_epi16(v_last_nogap, v_zero);
			v_last_nogap = _mm_max_epi16(v_last_nogap, v_a_gap);
			v_last_nogap = _mm_max_epi16(v_last_nogap, v_b_gap);
			
			v_prev_nogap = v_nogap;
			v_nogap = v_last_nogap;

			b_gap[j] = (int16_t)_mm_extract_epi16(v_b_gap, 7);
			nogap[j] = (int16_t)_mm_extract_epi16(v_nogap, 7);

			v_score = _mm_max_epi16(v_score, v_last_nogap);
		}
	}

	_mm_store_si128((__m128i *)w, v_score);
	for (i = 0; i < 8; i++)
		score = MAX(score, w[i]);

	return (score);
}

int
sw_vector_setup(int _dblen, int _qrlen, int _gap_open, int _gap_ext,
    int _match, int _mismatch, int _use_colours)
{

	dblen = _dblen;
	db = malloc((dblen + 14) * sizeof(db[0]));
	if (db == NULL)
		return (1);

	db_ls = malloc((dblen + 14) * sizeof(db_ls[0]));
	if (db_ls == NULL)
		return (1);

	qrlen = _qrlen;
	qr = malloc((qrlen + 14) * sizeof(qr[0]));
	if (qr == NULL)
		return (1);

	nogap = malloc((dblen + 14) * sizeof(nogap[0]));
	if (nogap == NULL)
		return (1);

	b_gap = malloc((dblen + 14) * sizeof(b_gap[0]));
	if (b_gap == NULL)
		return (1);

	gap_open = -(_gap_open);
	gap_ext = -(_gap_ext);
	match = _match;
	mismatch = _mismatch;
	use_colours = _use_colours;

	initialised = 1;

	return (0);
}

void
sw_vector_stats(uint64_t *invoc, uint64_t *cells, uint64_t *ticks,
    double *cellspersec)
{
	
	if (invoc != NULL)
		*invoc = swinvocs;
	if (cells != NULL)
		*cells = swcells;
	if (ticks != NULL)
		*ticks = swticks;
	if (cellspersec != NULL)
		*cellspersec = (double)swcells / ((double)swticks / cpuhz());
}

int
sw_vector(uint32_t *genome, int goff, int glen, uint32_t *read, int rlen,
    uint32_t *genome_ls, int initbp)
{
	uint64_t before, after;
	int i, score;

	before = rdtsc();

	if (!initialised)
		abort();

	swinvocs++;

	assert(glen > 0 && glen <= dblen);
	assert(rlen > 0 && rlen <= qrlen);

	memset(db, -1, (dblen + 14) * sizeof(db[0]));
	memset(qr, -2, (qrlen + 14) * sizeof(qr[0]));

	for (i = 0; i < glen; i++)
		db[i+7] = (int8_t)EXTRACT(genome, goff + i);

	if (genome_ls != NULL) {
		for (i = 0; i < glen; i++)
			db_ls[i+7] = (int8_t)EXTRACT(genome_ls, goff + i);
	}

	for (i = 0; i < rlen; i++)
		qr[i+7] = (int8_t)EXTRACT(read, i);

	score = vect_sw(&db[0], glen, &qr[7], rlen, &db_ls[0], initbp);

	swcells += (glen * rlen);
	after = rdtsc();
	swticks += (after - before);

	return (score);
}
