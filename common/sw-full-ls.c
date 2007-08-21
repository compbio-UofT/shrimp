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

#include "rmapper.h"
#include "sw-full.h"
#include "util.h"

struct swcell {
	int score_north;
	int back_north;

	int score_west;
	int back_west;

	int score_northwest;
	int back_northwest;
};

enum {
	FROM_NORTH_NORTH = 1,
	FROM_NORTH_NORTHWEST,

	FROM_WEST_NORTHWEST,
	FROM_WEST_WEST,

	FROM_NORTHWEST_NORTH,
	FROM_NORTHWEST_NORTHWEST,
	FROM_NORTHWEST_WEST
};

enum {
	BACK_INSERTION = 1,
	BACK_DELETION,
	BACK_MATCH_MISMATCH
};

static int		initialised;
static int8_t	       *db, *qr;
static int		dblen, qrlen;
static int		gap_open, gap_ext;
static int		match, mismatch;
static uint64_t		swticks, swcells, swinvocs;
static struct swcell   *swmatrix;
static int8_t	       *backtrace;
static char	       *dbalign, *qralign;

#define SWM(_i, _j) swmatrix[((_i) + 1) * (lena + 1) + (_j) + 1]

static void
full_sw(int lena, int lenb, int maxscore, int *iret, int *jret)
{
	int i, j;
	int score, ms, go, ge, tmp, tmp2;

	/* shut up gcc */
	j = 0;

	score = 0;
	go = gap_open;
	ge = gap_ext;

	for (i = -1; i < lenb; i++) {
		SWM(i, -1).score_northwest = 0;
		SWM(i, -1).score_north = -go;
		SWM(i, -1).score_west = -go;
	}

	for (i = -1; i < lena; i++) {
		SWM(-1, i).score_northwest = 0;
		SWM(-1, i).score_north = -go;
		SWM(-1, i).score_west = -go;
	}

	for (i = 0; i < lenb; i++) {
		for (j = 0; j < lena; j++) {

			/*
			 * northwest
			 */
			ms = (db[j] == qr[i]) ? match : mismatch;

			tmp  = SWM(i-1, j-1).score_northwest + ms;
			tmp2 = FROM_NORTHWEST_NORTHWEST;

			if (SWM(i-1, j-1).score_north + ms > tmp) {
				tmp  = SWM(i-1, j-1).score_north + ms;
				tmp2 = FROM_NORTHWEST_NORTH;
			}

			if (SWM(i-1, j-1).score_west + ms > tmp) {
				tmp  = SWM(i-1, j-1).score_west + ms;
				tmp2 = FROM_NORTHWEST_WEST;
			}

			if (tmp <= 0)
				tmp = tmp2 = 0;

			SWM(i, j).score_northwest = tmp;
			SWM(i, j).back_northwest  = tmp2;


			/*
			 * north
			 */
			tmp  = SWM(i-1, j).score_northwest - go - ge;
			tmp2 = FROM_NORTH_NORTHWEST;

			if (SWM(i-1, j).score_north - ge > tmp) {
				tmp  = SWM(i-1, j).score_north - ge;
				tmp2 = FROM_NORTH_NORTH;
			}

			if (tmp <= 0)
				tmp = tmp2 = 0;
				
			SWM(i, j).score_north = tmp;
			SWM(i, j).back_north  = tmp2;

			
			/*
			 * west
			 */
			tmp  = SWM(i, j-1).score_northwest - go - ge;
			tmp2 = FROM_WEST_NORTHWEST;

			if (SWM(i, j-1).score_west - ge > tmp) {
				tmp  = SWM(i, j-1).score_west - ge;
				tmp2 = FROM_WEST_WEST;
			}

			if (tmp <= 0)
				tmp = tmp2 = 0;

			SWM(i, j).score_west = tmp;
			SWM(i, j).back_west  = tmp2;


			/*
			 * max score
			 */
			score = MAX(score, SWM(i, j).score_northwest);
			score = MAX(score, SWM(i, j).score_north);
			score = MAX(score, SWM(i, j).score_west);

			if (score == maxscore)
				break;
		}

		if (score == maxscore)
			break;
	}
	if (score != maxscore) {
		fprintf(stderr, "error: full_sw failed to find maxscore!\n");
		exit(1);
	}

	*iret = i;
	*jret = j;
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
do_backtrace(int lena, int i, int j, struct sw_full_results *sfr)
{
	int k, from, fromscore;

	from = SWM(i, j).back_northwest;
	fromscore = SWM(i, j).score_northwest;
	if (SWM(i, j).score_west > fromscore) {
		from = SWM(i, j).back_west;
		fromscore = SWM(i, j).score_west;
	}
	if (SWM(i, j).score_north > fromscore)
		from = SWM(i, j).back_north;


	/* fill out the backtrace */
	k = (dblen + qrlen) - 1;
	while (i >= 0 && j >= 0) {
		assert(k >= 0);

		if (SWM(i, j).score_north <= 0 &&
		    SWM(i, j).score_northwest <= 0 &&
		    SWM(i, j).score_west <= 0)
			break;

		/* common operations first */
		switch (from) {
		case FROM_NORTH_NORTH:
		case FROM_NORTH_NORTHWEST:
			backtrace[k] = BACK_DELETION;
			sfr->deletions++;
			sfr->read_start = i--;
			break;

		case FROM_WEST_WEST:
		case FROM_WEST_NORTHWEST:
			backtrace[k] = BACK_INSERTION;
			sfr->insertions++;
			sfr->genome_start = j--;
			break;

		case FROM_NORTHWEST_NORTH:
		case FROM_NORTHWEST_NORTHWEST:
		case FROM_NORTHWEST_WEST:
			backtrace[k] = BACK_MATCH_MISMATCH;
			if (db[j] == qr[i])
				sfr->matches++;
			else
				sfr->mismatches++;
			sfr->read_start = i--;
			sfr->genome_start = j--;
			break;

		default:
			assert(0);
		}

		/* continue backtrace (nb: i and j have already been changed) */
		switch (from) {
		case FROM_NORTH_NORTH:
			from = SWM(i, j).back_north;
			break;

		case FROM_NORTH_NORTHWEST:
			from = SWM(i, j).back_northwest;
			break;

		case FROM_WEST_WEST:
			from = SWM(i, j).back_west;
			break;

		case FROM_WEST_NORTHWEST:
			from = SWM(i, j).back_northwest;
			break;

		case FROM_NORTHWEST_NORTH:
			from = SWM(i, j).back_north;
			break;

		case FROM_NORTHWEST_NORTHWEST:
			from = SWM(i, j).back_northwest;
			break;

		case FROM_NORTHWEST_WEST:
			from = SWM(i, j).back_west;
			break;

		default:
			assert(0);
		}

		k--;
	}

	return (k + 1);
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
	int l, done;
#ifdef USE_COLORS
	char translate[5] = { '0', '1', '2', '3', 'N' };
#else
	char translate[5] = { 'A', 'C', 'G', 'T', 'N' };
#endif


	d = dbalign;
	q = qralign;

	done = 0;
	for (l = k; l < (dblen + qrlen); l++) {
		switch (backtrace[l]) {
		case BACK_DELETION:
			*d++ = '-';
			*q++ = (char)(translate[(int)qr[i++]]);
			break;

		case BACK_INSERTION:
			*d++ = (char)(translate[(int)db[j++]]);
			*q++ = '-';
			break;

		case BACK_MATCH_MISMATCH:
			*d++ = (char)(translate[(int)db[j++]]);
			*q++ = (char)(translate[(int)qr[i++]]);
			break;

		default:
			done = 1;
		}
		
		if (done)
			break;
	}

	*d = *q = '\0';
}

int
sw_full_setup(int _dblen, int _qrlen, int _gap_open, int _gap_ext,
    int _match, int _mismatch)
{

	dblen = _dblen;
	db = malloc(dblen * sizeof(db[0]));
	if (db == NULL)
		return (1);

	qrlen = _qrlen;
	qr = malloc(qrlen * sizeof(qr[0]));
	if (qr == NULL)
		return (1);

	swmatrix = malloc((dblen + 1) * (qrlen + 1) * sizeof(swmatrix[0]));
	if (swmatrix == NULL)
		return (1);

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

	initialised = 1;

	return (0);
}

void
sw_full_stats(uint64_t *invoc, uint64_t *cells, uint64_t *ticks,
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
sw_full(uint32_t *genome, int goff, int glen, uint32_t *read, int rlen,
    int maxscore, char **dbalignp, char **qralignp, struct sw_full_results *sfr)
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
	
	for (i = 0; i < glen; i++)
		db[i] = (int8_t)EXTRACT(genome, goff + i);

	for (i = 0; i < rlen; i++)
		qr[i] = (int8_t)EXTRACT(read, i);

	full_sw(glen, rlen, maxscore, &i, &j);
	k = do_backtrace(glen, i, j, sfr);
	pretty_print(sfr->read_start, sfr->genome_start, k);
	sfr->mapped = i - sfr->read_start + 1;
	sfr->score = maxscore;

	swcells += (glen * rlen);
	after = rdtsc();
	swticks += (after - before);
}
