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

#include <sys/types.h>
#include <sys/stat.h>
#include <sys/time.h>

#include "../common/util.h"

#include "rmapper.h"
#include "sw-full-common.h"
#include "sw-full-cs.h"

struct swcell {
	int score_n;
	int back_n;

	int score_w;
	int back_w;

	int score_nw;
	int back_nw;
};

enum {
	FROM_NORTH_NORTH = 0,
	FROM_NORTH_NORTHWEST,

	FROM_WEST_NORTHWEST,
	FROM_WEST_WEST,

	FROM_NORTHWEST_NORTH,
	FROM_NORTHWEST_NORTHWEST,
	FROM_NORTHWEST_WEST
};

enum {
	FROM_A_NORTH_NORTH = 1,
	FROM_A_NORTH_NORTHWEST,

	FROM_A_WEST_NORTHWEST,
	FROM_A_WEST_WEST,

	FROM_A_NORTHWEST_NORTH,
	FROM_A_NORTHWEST_NORTHWEST,
	FROM_A_NORTHWEST_WEST,

	FROM_B_NORTH_NORTH,
	FROM_B_NORTH_NORTHWEST,

	FROM_B_WEST_NORTHWEST,
	FROM_B_WEST_WEST,

	FROM_B_NORTHWEST_NORTH,
	FROM_B_NORTHWEST_NORTHWEST,
	FROM_B_NORTHWEST_WEST,

	FROM_C_NORTH_NORTH,
	FROM_C_NORTH_NORTHWEST,

	FROM_C_WEST_NORTHWEST,
	FROM_C_WEST_WEST,

	FROM_C_NORTHWEST_NORTH,
	FROM_C_NORTHWEST_NORTHWEST,
	FROM_C_NORTHWEST_WEST,

	FROM_D_NORTH_NORTH,
	FROM_D_NORTH_NORTHWEST,

	FROM_D_WEST_NORTHWEST,
	FROM_D_WEST_WEST,

	FROM_D_NORTHWEST_NORTH,
	FROM_D_NORTHWEST_NORTHWEST,
	FROM_D_NORTHWEST_WEST
};

enum {
	BACK_INSERTION = 1,
	BACK_A_DELETION,
	BACK_B_DELETION,
	BACK_C_DELETION,
	BACK_D_DELETION,
	BACK_A_MATCH_MISMATCH,
	BACK_B_MATCH_MISMATCH,
	BACK_C_MATCH_MISMATCH,
	BACK_D_MATCH_MISMATCH
};

static int		initialised;
static int8_t	       *db, *qr[4];
static int		dblen, qrlen;
static int		gap_open, gap_ext;
static int		match, mismatch;
static int		xover_penalty;
static uint64_t		swticks, swcells, swinvocs;
static struct swcell   *swmatrix[4];
static uint8_t	       *backtrace;
static char	       *dbalign, *qralign;

#define SWM_x(_x, _i, _j) swmatrix[(_x)][((_i) + 1) * (lena + 1) + (_j) + 1]

#define BT_CROSSOVER		0x80
#define BT_ISCROSSOVER(_x)	((_x) & BT_CROSSOVER)
#define BT_TYPE(_x)		((_x) & 0x0f)

/*
 * Perform a full Smith-Waterman alignment. For the colour case, this means
 * computing each possible letter space read string and doing a four layer
 * scan.
 */
static int
full_sw(int lena, int lenb, int maxscore, int *iret, int *jret, int *kret)
{
	int i, j, k, l, n, max_i, max_j, max_k;
	int score, ms, go, ge, tmp, tmp2, resetval;
	const int neighbours[4][3] = {
		{ 1, 2, 3 },
		{ 0, 2, 3 },
		{ 0, 1, 3 },
		{ 0, 1, 2 }
	};
	const int from[4][7] = {
		{ FROM_A_NORTH_NORTH,		FROM_A_NORTH_NORTHWEST,
		  FROM_A_WEST_NORTHWEST,	FROM_A_WEST_WEST,
		  FROM_A_NORTHWEST_NORTH,	FROM_A_NORTHWEST_NORTHWEST,
		  FROM_A_NORTHWEST_WEST },

		{ FROM_B_NORTH_NORTH,		FROM_B_NORTH_NORTHWEST,
		  FROM_B_WEST_NORTHWEST,	FROM_B_WEST_WEST,
		  FROM_B_NORTHWEST_NORTH,	FROM_B_NORTHWEST_NORTHWEST,
		  FROM_B_NORTHWEST_WEST },

		{ FROM_C_NORTH_NORTH,		FROM_C_NORTH_NORTHWEST,
		  FROM_C_WEST_NORTHWEST,	FROM_C_WEST_WEST,
		  FROM_C_NORTHWEST_NORTH,	FROM_C_NORTHWEST_NORTHWEST,
		  FROM_C_NORTHWEST_WEST },

		{ FROM_D_NORTH_NORTH,		FROM_D_NORTH_NORTHWEST,
		  FROM_D_WEST_NORTHWEST,	FROM_D_WEST_WEST,
		  FROM_D_NORTHWEST_NORTH,	FROM_D_NORTHWEST_NORTHWEST,
		  FROM_D_NORTHWEST_WEST }
	};

	/* shut up gcc */
	max_i = max_j = max_k = j = 0;

	score = 0;
	go = gap_open;
	ge = gap_ext;

	for (i = 0; i < 4; i++) {
		for (j = -1; j < lenb; j++) {
			SWM_x(i, j, -1).score_nw = 0;
			SWM_x(i, j, -1).score_n = -go;
			SWM_x(i, j, -1).score_w = -go;
		}

		for (j = -1; j < lena; j++) {
			SWM_x(i, -1, j).score_nw = 0;
			SWM_x(i, -1, j).score_n  = -go;
			SWM_x(i, -1, j).score_w  = -go;
		}
	}

	for (i = 0; i < lenb; i++) {
		for (j = 0; j < lena; j++) {
			for (k = 0; k < 4; k++) {
				/*
				 * Do not permit a free crossover on the first
				 * letter of the read.
				 */
				if (i != 0)
					resetval = (2 * xover_penalty) - 1;
				else
					resetval = 0;

				/*
				 * northwest
				 */
				ms = (db[j] == qr[k][i]) ? match : mismatch;

				tmp  = SWM_x(k, i-1, j-1).score_nw + ms;
				tmp2 = from[k][FROM_NORTHWEST_NORTHWEST];

				if (SWM_x(k, i-1, j-1).score_n + ms > tmp) {
					tmp  = SWM_x(k, i-1, j-1).score_n + ms;
					tmp2 = from[k][FROM_NORTHWEST_NORTH];
				}

				if (SWM_x(k, i-1, j-1).score_w + ms > tmp) {
					tmp  = SWM_x(k, i-1, j-1).score_w + ms;
					tmp2 = from[k][FROM_NORTHWEST_WEST];
				}

				/* check neighbours' NW */
				for (l = 0; l < 3; l++) {
					n = neighbours[k][l];
					ms = (db[j] == qr[n][i]) ?
					    match : mismatch;

					if (SWM_x(n, i-1, j-1).score_nw + ms +
					    xover_penalty > tmp) {
						tmp  =
						    SWM_x(n, i-1, j-1).score_nw
						    + ms + xover_penalty;
						tmp2 = from[n][
						    FROM_NORTHWEST_NORTHWEST];
					}
				}

				/* check neighbours' N */
				for (l = 0; l < 3; l++) {
					n = neighbours[k][l];

					if (SWM_x(n, i-1, j-1).score_n + ms +
					    xover_penalty > tmp) {
						tmp  =
						    SWM_x(n, i-1, j-1).score_n
						    + ms + xover_penalty;
						tmp2 = from[n][
						    FROM_NORTHWEST_NORTH];
					}
				}

				/* check neighbour's W */
				for (l = 0; l < 3; l++) {
					n = neighbours[k][l];

					if (SWM_x(n, i-1, j-1).score_w + ms +
					    xover_penalty > tmp) {
						tmp  =
						    SWM_x(n, i-1, j-1).score_w
						    + ms + xover_penalty;
						tmp2 = from[n][
						    FROM_NORTHWEST_WEST];
					}
				}

				if (tmp <= 0)
					tmp = tmp2 = resetval;

				SWM_x(k, i, j).score_nw = tmp;
				SWM_x(k, i, j).back_nw  = tmp2;


				/*
				 * north
				 */
				tmp  = SWM_x(k, i-1, j).score_nw - go - ge;
				tmp2 = from[k][FROM_NORTH_NORTHWEST];

				if (SWM_x(k, i-1, j).score_n - ge > tmp) {
					tmp  = SWM_x(k, i-1, j).score_n - ge;
					tmp2 = from[k][FROM_NORTH_NORTH];
				}

				/* check neighbours' NW */
				for (l = 0; l < 3; l++) {
					n = neighbours[k][l];

					if (SWM_x(n, i-1, j).score_nw - go - ge
					    + xover_penalty > tmp) {
						tmp  =
						    SWM_x(n, i-1, j).score_nw
						    - go - ge + xover_penalty;
						tmp2 = from[n][
						    FROM_NORTH_NORTHWEST];
					}
				}

				/* check neighbours' N */
				for (l = 0; l < 3; l++) {
					n = neighbours[k][l];

					if (SWM_x(n, i-1, j).score_n - ge
					    + xover_penalty > tmp) {
						tmp  =
						    SWM_x(n, i-1, j).score_n
						    - ge + xover_penalty;
						tmp2 = from[n][
						    FROM_NORTH_NORTH];
					}
				}


				if (tmp <= 0)
					tmp = tmp2 = resetval;
					
				SWM_x(k, i, j).score_n = tmp;
				SWM_x(k, i, j).back_n  = tmp2;

				
				/*
				 * west
				 */
				tmp  = SWM_x(k, i, j-1).score_nw - go - ge;
				tmp2 = from[k][FROM_WEST_NORTHWEST];

				if (SWM_x(k, i, j-1).score_w - ge > tmp) {
					tmp  = SWM_x(k, i, j-1).score_w - ge;
					tmp2 = from[k][FROM_WEST_WEST];
				}

/* it doesn't make sense to cross over on a genomic gap */
#if 0
				/* check neighbours' NW */
				for (l = 0; l < 3; l++) {
					n = neighbours[k][l];

					if (SWM_x(n, i, j-1).score_nw - go - ge
					    + xover_penalty > tmp) {
						tmp  =
						    SWM_x(n, i, j-1).score_nw
						    - go - ge + xover_penalty;
						tmp2 = from[n][
						    FROM_WEST_NORTHWEST];
					}
				}

				/* check neighbours' W */
				for (l = 0; l < 3; l++) {
					n = neighbours[k][l];

					if (SWM_x(n, i, j-1).score_w - ge
					    + xover_penalty > tmp) {
						tmp  =
						    SWM_x(n, i-1, j-1).score_w
						    - ge + xover_penalty;
						tmp2 = from[n][
						    FROM_WEST_WEST];
					}
				}
#endif
				if (tmp <= 0)
					tmp = tmp2 = resetval;

				SWM_x(k, i, j).score_w = tmp;
				SWM_x(k, i, j).back_w  = tmp2;


				/*
				 * max score
				 */
				if (SWM_x(k, i, j).score_nw > score) {
					score = SWM_x(k, i, j).score_nw;
					max_i = i, max_j = j, max_k = k;
				}
				if (SWM_x(k, i, j).score_n > score) {
					score = SWM_x(k, i, j).score_n;
					max_i = i, max_j = j, max_k = k;
				}
				if (SWM_x(k, i, j).score_w > score) {
					score = SWM_x(k, i, j).score_w;
					max_i = i, max_j = j, max_k = k;
				}
			}
		}
	}

	// not necessarily true, but typically so, and should be close
	// if not. can we define a correct assertion here?
	//assert(score >= maxscore);
	(void)maxscore;

	*iret = max_i;
	*jret = max_j;
	*kret = max_k;

	return (score);
}

/*
 * Fill in the backtrace in order to do a pretty printout.
 *
 * Returns the beginning matrix cell (i, j) in 'sfr->read_start' and
 * 'sfr->genome_start'.
 *
 * The return value is the first valid offset in the backtrace buffer.
 */
static int
do_backtrace(int lena, int i, int j, int k, struct sw_full_results *sfr)
{
	int off, from, fromscore;

	from = SWM_x(k, i, j).back_nw;
	fromscore = SWM_x(k, i, j).score_nw;
	if (SWM_x(k, i, j).score_w > fromscore) {
		from = SWM_x(k, i, j).back_w;
		fromscore = SWM_x(k, i, j).score_w;
	}
	if (SWM_x(k, i, j).score_n > fromscore)
		from = SWM_x(k, i, j).back_n;


	/* fill out the backtrace */
	off = (dblen + qrlen) - 1;
	while (i >= 0 && j >= 0) {
		assert(off >= 0);

		if (SWM_x(k, i, j).score_n <= 0 &&
		    SWM_x(k, i, j).score_nw <= 0 &&
		    SWM_x(k, i, j).score_w <= 0)
			break;

		/* common operations first */
		switch (from) {
		case FROM_A_NORTH_NORTH:
		case FROM_A_NORTH_NORTHWEST:
		case FROM_B_NORTH_NORTH:
		case FROM_B_NORTH_NORTHWEST:
		case FROM_C_NORTH_NORTH:
		case FROM_C_NORTH_NORTHWEST:
		case FROM_D_NORTH_NORTH:
		case FROM_D_NORTH_NORTHWEST:
			sfr->deletions++;
			sfr->read_start = i--;
			break;

		case FROM_A_WEST_WEST:
		case FROM_A_WEST_NORTHWEST:
		case FROM_B_WEST_WEST:
		case FROM_B_WEST_NORTHWEST:
		case FROM_C_WEST_WEST:
		case FROM_C_WEST_NORTHWEST:
		case FROM_D_WEST_WEST:
		case FROM_D_WEST_NORTHWEST:
			backtrace[off] = BACK_INSERTION;
			sfr->insertions++;
			sfr->genome_start = j--;
			break;

		case FROM_A_NORTHWEST_NORTH:
		case FROM_A_NORTHWEST_NORTHWEST:
		case FROM_A_NORTHWEST_WEST:
		case FROM_B_NORTHWEST_NORTH:
		case FROM_B_NORTHWEST_NORTHWEST:
		case FROM_B_NORTHWEST_WEST:
		case FROM_C_NORTHWEST_NORTH:
		case FROM_C_NORTHWEST_NORTHWEST:
		case FROM_C_NORTHWEST_WEST:
		case FROM_D_NORTHWEST_NORTH:
		case FROM_D_NORTHWEST_NORTHWEST:
		case FROM_D_NORTHWEST_WEST:
			if (db[j] == qr[k][i])
				sfr->matches++;
			else
				sfr->mismatches++;
			sfr->read_start = i--;
			sfr->genome_start = j--;
			break;

		default:
			assert(0);
		}

		/* handle match/mismatch and north */
		switch (from) {
		case FROM_A_NORTH_NORTH:
		case FROM_A_NORTH_NORTHWEST:
			backtrace[off] = BACK_A_DELETION;
			break;

		case FROM_B_NORTH_NORTH:
		case FROM_B_NORTH_NORTHWEST:
			backtrace[off] = BACK_B_DELETION;
			break;

		case FROM_C_NORTH_NORTH:
		case FROM_C_NORTH_NORTHWEST:
			backtrace[off] = BACK_C_DELETION;
			break;

		case FROM_D_NORTH_NORTH:
		case FROM_D_NORTH_NORTHWEST:
			backtrace[off] = BACK_D_DELETION;
			break;

		case FROM_A_NORTHWEST_NORTH:
		case FROM_A_NORTHWEST_NORTHWEST:
		case FROM_A_NORTHWEST_WEST:
			backtrace[off] = BACK_A_MATCH_MISMATCH;
			break;

		case FROM_B_NORTHWEST_NORTH:
		case FROM_B_NORTHWEST_NORTHWEST:
		case FROM_B_NORTHWEST_WEST:
			backtrace[off] = BACK_B_MATCH_MISMATCH;
			break;

		case FROM_C_NORTHWEST_NORTH:
		case FROM_C_NORTHWEST_NORTHWEST:
		case FROM_C_NORTHWEST_WEST:
			backtrace[off] = BACK_C_MATCH_MISMATCH;
			break;

		case FROM_D_NORTHWEST_NORTH:
		case FROM_D_NORTHWEST_NORTHWEST:
		case FROM_D_NORTHWEST_WEST:
			backtrace[off] = BACK_D_MATCH_MISMATCH;
			break;
		}

		/* continue backtrace (nb: i and j have already been changed) */
		switch (from) {
		case FROM_A_NORTH_NORTH:
		case FROM_B_NORTH_NORTH:
		case FROM_C_NORTH_NORTH:
		case FROM_D_NORTH_NORTH:
			from = SWM_x(k, i, j).back_n;
			break;

		case FROM_A_NORTH_NORTHWEST:
		case FROM_B_NORTH_NORTHWEST:
		case FROM_C_NORTH_NORTHWEST:
		case FROM_D_NORTH_NORTHWEST:
			from = SWM_x(k, i, j).back_nw;
			break;

		case FROM_A_WEST_WEST:
		case FROM_B_WEST_WEST:
		case FROM_C_WEST_WEST:
		case FROM_D_WEST_WEST:
			from = SWM_x(k, i, j).back_w;
			break;

		case FROM_A_WEST_NORTHWEST:
		case FROM_B_WEST_NORTHWEST:
		case FROM_C_WEST_NORTHWEST:
		case FROM_D_WEST_NORTHWEST:
			from = SWM_x(k, i, j).back_nw;
			break;

		case FROM_A_NORTHWEST_NORTH:
		case FROM_B_NORTHWEST_NORTH:
		case FROM_C_NORTHWEST_NORTH:
		case FROM_D_NORTHWEST_NORTH:
			from = SWM_x(k, i, j).back_n;
			break;

		case FROM_A_NORTHWEST_NORTHWEST:
		case FROM_B_NORTHWEST_NORTHWEST:
		case FROM_C_NORTHWEST_NORTHWEST:
		case FROM_D_NORTHWEST_NORTHWEST:
			from = SWM_x(k, i, j).back_nw;
			break;

		case FROM_A_NORTHWEST_WEST:
		case FROM_B_NORTHWEST_WEST:
		case FROM_C_NORTHWEST_WEST:
		case FROM_D_NORTHWEST_WEST:
			from = SWM_x(k, i, j).back_w;
			break;

		default:
			assert(0);
		}

		/* set k */
		switch (from) {
		case FROM_A_NORTH_NORTH:
		case FROM_A_NORTH_NORTHWEST:
		case FROM_A_WEST_WEST:
		case FROM_A_WEST_NORTHWEST:
		case FROM_A_NORTHWEST_NORTH:
		case FROM_A_NORTHWEST_NORTHWEST:
		case FROM_A_NORTHWEST_WEST:
			if (k != 0) {
				backtrace[off] |= BT_CROSSOVER;
				sfr->crossovers++;
				k = 0;
			}
			break;

		case FROM_B_NORTH_NORTH:
		case FROM_B_NORTH_NORTHWEST:
		case FROM_B_WEST_WEST:
		case FROM_B_WEST_NORTHWEST:
		case FROM_B_NORTHWEST_NORTH:
		case FROM_B_NORTHWEST_NORTHWEST:
		case FROM_B_NORTHWEST_WEST:
			if (k != 1) {
				backtrace[off] |= BT_CROSSOVER;
				sfr->crossovers++;
				k = 1;
			}
			break;

		case FROM_C_NORTH_NORTH:
		case FROM_C_NORTH_NORTHWEST:
		case FROM_C_WEST_WEST:
		case FROM_C_WEST_NORTHWEST:
		case FROM_C_NORTHWEST_NORTH:
		case FROM_C_NORTHWEST_NORTHWEST:
		case FROM_C_NORTHWEST_WEST:
			if (k != 2) {
				backtrace[off] |= BT_CROSSOVER;
				sfr->crossovers++;
				k = 2;
			}
			break;

		case FROM_D_NORTH_NORTH:
		case FROM_D_NORTH_NORTHWEST:
		case FROM_D_WEST_WEST:
		case FROM_D_WEST_NORTHWEST:
		case FROM_D_NORTHWEST_NORTH:
		case FROM_D_NORTHWEST_NORTHWEST:
		case FROM_D_NORTHWEST_WEST:
			if (k != 3) {
				backtrace[off] |= BT_CROSSOVER;
				sfr->crossovers++;
				k = 3;
			}
			break;
		}

		off--;
	}

	return (off + 1);
}

/*
 * Pretty print our alignment of 'db' and 'qr' in 'dbalign' and 'qralign'.
 *
 * i, j represent the beginning cell in the matrix.
 * k is the first valid offset in the backtrace buffer.
 */
static void
pretty_print(int i, int j, int k)
{
	char *d, *q;
	int l;
	char utranslate[5] = { 'A', 'C', 'G', 'T', 'N' };
	char ltranslate[5] = { 'a', 'c', 'g', 't', 'n' };

	d = dbalign;
	q = qralign;

	for (l = k; l < (dblen + qrlen); l++) {
		switch (BT_TYPE(backtrace[l])) {
		case BACK_A_DELETION:
			*d++ = '-';
			*q++ = utranslate[(int)qr[0][i++]];
			break;

		case BACK_B_DELETION:
			*d++ = '-';
			*q++ = utranslate[(int)qr[1][i++]];
			break;

		case BACK_C_DELETION:
			*d++ = '-';
			*q++ = utranslate[(int)qr[2][i++]];
			break;

		case BACK_D_DELETION:
			*d++ = '-';
			*q++ = utranslate[(int)qr[3][i++]];
			break;

		case BACK_INSERTION:
			*d++ = utranslate[(int)db[j++]];
			*q++ = '-';
			break;

		case BACK_A_MATCH_MISMATCH:
			*d++ = utranslate[(int)db[j++]];
			if (BT_ISCROSSOVER(backtrace[l]))
				*q++ = ltranslate[(int)qr[0][i++]];
			else
				*q++ = utranslate[(int)qr[0][i++]];
			break;

		case BACK_B_MATCH_MISMATCH:
			*d++ = utranslate[(int)db[j++]];
			if (BT_ISCROSSOVER(backtrace[l]))
				*q++ = ltranslate[(int)qr[1][i++]];
			else
				*q++ = utranslate[(int)qr[1][i++]];
			break;

		case BACK_C_MATCH_MISMATCH:
			*d++ = utranslate[(int)db[j++]];
			if (BT_ISCROSSOVER(backtrace[l]))
				*q++ = ltranslate[(int)qr[2][i++]];
			else
				*q++ = utranslate[(int)qr[2][i++]];
			break;

		case BACK_D_MATCH_MISMATCH:
			*d++ = utranslate[(int)db[j++]];
			if (BT_ISCROSSOVER(backtrace[l]))
				*q++ = ltranslate[(int)qr[3][i++]];
			else
				*q++ = utranslate[(int)qr[3][i++]];
			break;
	
		default:
			printf("%x\n", backtrace[l]);
			assert(0);
		}
	}

	*d = *q = '\0';
}

int
sw_full_cs_setup(int _dblen, int _qrlen, int _gap_open, int _gap_ext,
    int _match, int _mismatch, int _xover_penalty)
{
	int i;

	dblen = _dblen;
	db = malloc(dblen * sizeof(db[0]));
	if (db == NULL)
		return (1);

	qrlen = _qrlen;
	for (i = 0; i < 4; i++) {
		qr[i] = malloc(qrlen * sizeof(qr[0]));
		if (qr[i] == NULL)
			return (1);
	}

	for (i = 0; i < 4; i++) {
		swmatrix[i] = malloc((dblen + 1) * (qrlen + 1) *
		    sizeof(swmatrix[0][0]));
		if (swmatrix[i] == NULL)
			return (1);
	}

	backtrace = malloc((dblen + qrlen) * sizeof(backtrace[0]));
	if (backtrace == NULL)
		return (1);

	dbalign = malloc((dblen + qrlen + 1) * sizeof(dbalign[0]));
	if (dbalign == NULL)
		return (1);

	qralign = malloc((dblen + qrlen + 1) * sizeof(dbalign[0]));
	if (qralign == NULL)
		return (1);

	gap_open = -(_gap_open);
	gap_ext = -(_gap_ext);
	match = _match;
	mismatch = _mismatch;
	xover_penalty = _xover_penalty;

	initialised = 1;

	return (0);
}

void
sw_full_cs_stats(uint64_t *invoc, uint64_t *cells, uint64_t *ticks,
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

void
sw_full_cs(uint32_t *genome, int goff, int glen, uint32_t *read, int rlen,
    uint32_t *genome_ls, int initbp, int maxscore, char **dbalignp,
    char **qralignp, struct sw_full_results *sfr)
{
	struct sw_full_results scratch;
	uint64_t before, after;
	int i, j, k;

	before = rdtsc();

	if (!initialised)
		abort();

	swinvocs++;

	assert(glen > 0 && glen <= dblen);
	assert(rlen > 0 && rlen <= qrlen);

	if (sfr == NULL)
		sfr = &scratch;

	memset(sfr, 0, sizeof(*sfr));
	memset(backtrace, 0, (dblen + qrlen) * sizeof(backtrace[0]));

	dbalign[0] = qralign[0] = '\0';

	if (dbalignp != NULL)
		*dbalignp = dbalign;
	if (qralignp != NULL)
		*qralignp = qralign;
	

	if (genome_ls == NULL) {
		for (i = 0; i < glen; i++)
			db[i] = (int8_t)EXTRACT(genome, goff + i);
	} else {
		for (i = 0; i < glen; i++)
			db[i] = (int8_t)EXTRACT(genome_ls, goff + i);
	}

	/*
	 * Generate each possible letter space sequence from the colour space
	 * read. qr[0] corresponds to initbp, which is given initial preference.
	 */
	assert(initbp >= 0 && initbp <= 3);
	for (i = 0; i < 4; i++) {
		int letter = (i + initbp) % 4;

		for (j = 0; j < rlen; j++) {
			qr[i][j] = (int8_t)cstols(letter, EXTRACT(read, j));
			letter = qr[i][j];
		}
	}

	sfr->score = full_sw(glen, rlen, maxscore, &i, &j, &k);
	k = do_backtrace(glen, i, j, k, sfr);
	pretty_print(sfr->read_start, sfr->genome_start, k);
	sfr->mapped = i - sfr->read_start + 1;

	swcells += (glen * rlen);
	after = rdtsc();
	swticks += (after - before);
}
