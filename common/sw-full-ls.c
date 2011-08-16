/*	$Id: sw-full-ls.c,v 1.16 2009/06/12 21:27:35 rumble Exp $	*/

#include <assert.h>
#include <ctype.h>
#include <errno.h>
#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <unistd.h>
#include <zlib.h>

#include <sys/types.h>
#include <sys/stat.h>
#include <sys/time.h>
#include <limits.h>

#include "../common/fasta.h"
#include "../common/sw-full-common.h"
#include "../common/sw-full-ls.h"
#include "../common/util.h"
#include "../common/time_counter.h"

typedef struct swcell {
	int	score_north;
	int	score_west;
	int	score_northwest;

	int8_t	back_north;
	int8_t	back_west;
	int8_t	back_northwest;
} swcell;

#define FROM_NORTH_NORTH		0x1
#define FROM_NORTH_NORTHWEST		0x2
#define	FROM_WEST_NORTHWEST		0x3
#define	FROM_WEST_WEST			0x4
#define FROM_NORTHWEST_NORTH		0x5
#define FROM_NORTHWEST_NORTHWEST	0x6
#define FROM_NORTHWEST_WEST		0x7

#define BACK_INSERTION			0x1
#define BACK_DELETION			0x2
#define BACK_MATCH_MISMATCH		0x3

static int		initialised;
static int8_t	       *db, *qr;
static int		dblen, qrlen;
static int		a_gap_open, a_gap_ext;
static int		b_gap_open, b_gap_ext;
static int		match, mismatch;
static struct swcell   *swmatrix;
static int8_t	       *backtrace;
static char	       *dbalign, *qralign;
static int		anchor_width;

/* statistics */
static uint64_t		swcells, swinvocs;
static time_counter	sw_tc;

#pragma omp threadprivate(initialised,db,qr,dblen,qrlen,a_gap_open,a_gap_ext,b_gap_open,b_gap_ext,\
		match,mismatch,swmatrix,backtrace,dbalign,qralign,anchor_width,sw_tc,swcells,swinvocs)


inline static void
init_cell(int idx, int local_alignment) {
  if (local_alignment) {
	  swmatrix[idx].score_northwest = 0;
	  swmatrix[idx].score_north = -b_gap_open;
	  swmatrix[idx].score_west = -a_gap_open;
  } else {
	  swmatrix[idx].score_northwest = -INT_MAX/2;
	  swmatrix[idx].score_north = -INT_MAX/2;
	  swmatrix[idx].score_west = -INT_MAX/2;
  }
  swmatrix[idx].back_northwest = 0;
  swmatrix[idx].back_north = 0;
  swmatrix[idx].back_west = 0;
}


#ifdef DEBUG_SW
static void print_sw(int lena, int lenb) {
	int i,j;
	printf("      %5s ","-");
	for (j=1; j< lenb+1; j++) {
		printf("%5c ",base_translate(qr[j-1],false)); 
	}
	printf("\n");
	//rows
	for (i=0; i<lena+1; i++) {
		//cols
		if (i==0) {
			printf("    - ");
		} else {
			printf("%5c ",base_translate(db[i-1],false));
		}
		for (j=0; j<lenb+1; j++) {
			swcell curr=swmatrix[j*(lena+1)+i];
			int tmp=0;
			tmp=MAX(curr.score_north,curr.score_west);
			tmp=MAX(tmp,curr.score_northwest);
			if (tmp<-1000) {
				printf("%5d ",-99);
			} else {
			  //printf("%5d ",tmp);
			  printf("%5d/%5d/%5d ", curr.score_north,curr.score_west,curr.score_northwest);
			}
		}
		printf("\n");
	}
}

static void print_sw_backtrace(int lena, int lenb) {
	int i,j;
	printf("      %5s ","-");
	for (j=1; j< lenb+1; j++) {
		printf("%5c ",base_translate(qr[j-1],false)); 
	}
	printf("\n");
	//rows
	for (i=0; i<lena+1; i++) {
		//cols
		if (i==0) {
			printf("    - ");
		} else {
			printf("%5c ",base_translate(db[i-1],false));
		}
		for (j=0; j<lenb+1; j++) {
			swcell curr=swmatrix[j*(lena+1)+i];
			int btrace[3]={0,0,0};
			int maxscore=0;
			maxscore=MAX(curr.score_north,curr.score_west);
			maxscore=MAX(maxscore,curr.score_northwest);
			if (curr.score_west==maxscore) {
				btrace[0]=curr.back_west;
			}	
			if (curr.score_northwest==maxscore) {
				btrace[1]=curr.back_northwest;
			}
			if (curr.score_north==maxscore) {
				btrace[2]=curr.back_north;
			}
			printf("%d/%d/%d ",btrace[0],btrace[1],btrace[2]);
		}
		printf("\n");
	}
}
#endif


static int
full_sw(int lena, int lenb, int threshscore, int maxscore, int *iret, int *jret, bool revcmpl,
	struct anchor * anchors, int anchors_cnt, int local_alignment)
{
  //fprintf(stderr,"Executing full_sw\n");
  int max_i=0; int max_j=0;
  int i, j;
  //int sw_band, ne_band;
  int score, ms, a_go, a_ge, b_go, b_ge, tmp;
  int8_t tmp2;
  struct anchor rectangle;

  /* shut up gcc */
  j = 0;

  score = 0;
  a_go = a_gap_open;
  a_ge = a_gap_ext;
  b_go = b_gap_open;
  b_ge = b_gap_ext;

  if (anchors != NULL && anchor_width >= 0) {
    anchor_join(anchors, anchors_cnt, &rectangle);
    anchor_widen(&rectangle, anchor_width);
  } else {
    struct anchor tmp_anchors[2];

    tmp_anchors[0].x = 0;
    tmp_anchors[0].y = (lenb * match - threshscore) / match;
    tmp_anchors[0].length = 1;
    tmp_anchors[0].width = 1;

    tmp_anchors[1].x = lena-1;
    tmp_anchors[1].y = lenb-1-tmp_anchors[0].y;
    tmp_anchors[1].length = 1;
    tmp_anchors[1].width = 1;

    anchor_join(tmp_anchors, 2, &rectangle);
  }

  for (j = 0; j < lena + 1; j++) {
    init_cell(j,1);
  }

  /*
  for (i = 0; i < lenb + 1; i++) {
    int idx = i * (lena + 1);

    swmatrix[idx].score_northwest = 0;
    swmatrix[idx].score_north = 0;
    swmatrix[idx].score_west = 0;

    swmatrix[idx].back_northwest = 0;
    swmatrix[idx].back_north = 0;
    swmatrix[idx].back_west = 0;
  }
  */

  /*
   * Figure out our band.
   *   We can actually skip computation of a significant number of
   *   cells, which could never be part of an alignment corresponding
   *   to our threshhold score.
   */
  //sw_band = ((lenb * match - threshscore + match - 1) / match) + 1;
  //ne_band = lena - (lenb - sw_band);

  for (i = 0; i < lenb; i++) {
    /*
     * computing row i of virtual matrix, stored in row i+1
     */
    int x_min, x_max;

    anchor_get_x_range(&rectangle, lena, lenb, i, &x_min, &x_max);
    if (!local_alignment) {
	//init_cell((i + 1) * (lena + 1) + (x_min - 1) + 1, x_min == 0 ?  1 : 0);
	init_cell((i + 1) * (lena + 1) + (x_min - 1) + 1, 0);
    } else {
    	init_cell((i + 1) * (lena + 1) + (x_min - 1) + 1,1);
    }
    //if (x_min > 0) {
    //}

    swcells += x_max - x_min + 1;

    for (j = x_min; j <= x_max; j++) {
      /*
       * computing column j of virtual matrix, stored in column j+1
       */
      struct swcell *cell_nw, *cell_n, *cell_w, *cell_cur;

      cell_nw  = &swmatrix[i * (lena + 1) + j];
      cell_n   = cell_nw + 1; 
      cell_w   = cell_nw + (lena + 1);
      cell_cur = cell_w + 1;

      /* banding */
      //if (i >= sw_band + j) {
      //memset(cell_cur, 0, sizeof(*cell_cur));
      //continue;
      //}
      //if (j >= ne_band + i) {
      //memset(cell_cur, 0, sizeof(*cell_cur));
      //break;
      //}

      /*
       * northwest
       */
      ms = (db[j] == qr[i]) ? match : mismatch;

      if (!revcmpl) {
	tmp  = cell_nw->score_northwest + ms;
	tmp2 = FROM_NORTHWEST_NORTHWEST;

	if (cell_nw->score_north + ms > tmp) {
	  tmp  = cell_nw->score_north + ms;
	  tmp2 = FROM_NORTHWEST_NORTH;
	}

	if (cell_nw->score_west + ms > tmp) {
	  tmp  = cell_nw->score_west + ms;
	  tmp2 = FROM_NORTHWEST_WEST;
	}
      } else {
	tmp  = cell_nw->score_west + ms;
	tmp2 = FROM_NORTHWEST_WEST;

	if (cell_nw->score_north + ms > tmp) {
	  tmp  = cell_nw->score_north + ms;
	  tmp2 = FROM_NORTHWEST_NORTH;
	}

	if (cell_nw->score_northwest + ms > tmp) {
	  tmp  = cell_nw->score_northwest + ms;
	  tmp2 = FROM_NORTHWEST_NORTHWEST;
	}
      }

      if (tmp <= 0 && local_alignment)
	tmp = tmp2 = 0;

      cell_cur->score_northwest = tmp;
      cell_cur->back_northwest  = tmp2;


      /*
       * north
       */
      if (!revcmpl) {
	tmp  = cell_n->score_northwest - b_go - b_ge;
	tmp2 = FROM_NORTH_NORTHWEST;

	if (cell_n->score_north - b_ge > tmp) {
	  tmp  = cell_n->score_north - b_ge;
	  tmp2 = FROM_NORTH_NORTH;
	}
      } else {
	tmp  = cell_n->score_north - b_ge;
	tmp2 = FROM_NORTH_NORTH;

	if (cell_n->score_northwest - b_go - b_ge > tmp) {
	  tmp  = cell_n->score_northwest - b_go - b_ge;
	  tmp2 = FROM_NORTH_NORTHWEST;
	}
      }

      if (tmp <= 0 && local_alignment)
	tmp = tmp2 = 0;
				
      cell_cur->score_north = tmp;
      cell_cur->back_north  = tmp2;

			
      /*
       * west
       */
      if (!revcmpl) {
	tmp  = cell_w->score_northwest - a_go - a_ge;
	tmp2 = FROM_WEST_NORTHWEST;

	if (cell_w->score_west - a_ge > tmp) {
	  tmp  = cell_w->score_west - a_ge;
	  tmp2 = FROM_WEST_WEST;
	}
      } else {
	tmp  = cell_w->score_west - a_ge;
	tmp2 = FROM_WEST_WEST;

	if (cell_w->score_northwest - a_go - a_ge > tmp) {
	  tmp  = cell_w->score_northwest - a_go - a_ge;
	  tmp2 = FROM_WEST_NORTHWEST;
	}
      }

      if (tmp <= 0 && local_alignment)
	 tmp = tmp2 = 0;

      cell_cur->score_west = tmp;
      cell_cur->back_west  = tmp2;


      /*
       * max score
       */
      if (local_alignment || i==lenb-1) {
	      int tmp;
	      tmp = MAX(cell_cur->score_north, cell_cur->score_northwest);
	      tmp = MAX(tmp, cell_cur->score_west);
	      if (tmp>score) {
		score=tmp;
	      	max_i = i;
	      	max_j = j;
	      }
      }

      if (score == maxscore && local_alignment)
	break;
    }

    if (score == maxscore && local_alignment)
      break;
 
   if (i+1 < lenb) {
      int next_x_min, next_x_max;

      anchor_get_x_range(&rectangle, lena, lenb, i+1, &next_x_min, &next_x_max);
      for (j = x_max + 1; j <= next_x_max; j++) {
	init_cell((i + 1) * (lena + 1) + (j + 1),local_alignment);
      }
    }
  }

  *iret = max_i;
  *jret = max_j;
#ifdef DEBUG_SW
  fprintf(stderr,"Returning i = %d, j= %d, score= %d , maxscore=%d\n",i,j,score,maxscore);
  print_sw(lena,lenb);
  print_sw_backtrace(lena,lenb);
  fprintf(stderr,"Final score is %d\n",score);
#endif
  if (score == maxscore || !local_alignment)
    return score;
  else if (anchors != NULL)
    return full_sw(lena, lenb, threshscore, maxscore, iret, jret, revcmpl, NULL, 0,local_alignment);
  else {
    assert(0);
    return 0;
  }
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
	struct swcell *cell;
	int k, from, fromscore;

	cell = &swmatrix[(i + 1) * (lena + 1) + j + 1];
	from = cell->back_northwest;
	fromscore = cell->score_northwest;
	if (cell->score_west > fromscore) {
		from = cell->back_west;
		fromscore = cell->score_west;
	}
	if (cell->score_north > fromscore)
		from = cell->back_north;

	assert(from != 0);

	/* fill out the backtrace */
	k = (dblen + qrlen) - 1;
	while (i >= 0 && j >= 0) {
		//printf("Got cell %d , %d for backtrace\n",i+1,j+1);
		assert(k >= 0);

		cell = NULL;

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
			fprintf(stderr, "INTERNAL ERROR: from = %d\n", from);
			assert(0);
		}

		/* continue backtrace (nb: i and j have already been changed) */
		cell = &swmatrix[(i + 1) * (lena + 1) + j + 1];

		switch (from) {
		case FROM_NORTH_NORTH:
			from = cell->back_north;
			break;

		case FROM_NORTH_NORTHWEST:
			from = cell->back_northwest;
			break;

		case FROM_WEST_WEST:
			from = cell->back_west;
			break;

		case FROM_WEST_NORTHWEST:
			from = cell->back_northwest;
			break;

		case FROM_NORTHWEST_NORTH:
			from = cell->back_north;
			break;

		case FROM_NORTHWEST_NORTHWEST:
			from = cell->back_northwest;
			break;

		case FROM_NORTHWEST_WEST:
			from = cell->back_west;
			break;

		default:
			fprintf(stderr, "INTERNAL ERROR: from = %d\n", from);
			assert(0);
		}

		k--;

		if (from == 0)
			break;
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

	d = dbalign;
	q = qralign;

	done = 0;
	for (l = k; l < (dblen + qrlen); l++) {
		switch (backtrace[l]) {
		case BACK_DELETION:
			*d++ = '-';
			*q++ = base_translate(qr[i++], false);
			break;

		case BACK_INSERTION:
			*d++ = base_translate(db[j++], false);
			*q++ = '-';
			break;

		case BACK_MATCH_MISMATCH:
			*d++ = base_translate(db[j++], false);
			*q++ = base_translate(qr[i++], false);
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
sw_full_ls_cleanup(void)
{
	free(db);
	free(qr);
	free(swmatrix);
	free(backtrace);
	free(dbalign);
	free(qralign);
	return (0);
}
int
sw_full_ls_setup(int _dblen, int _qrlen, int _a_gap_open, int _a_gap_ext,
    int _b_gap_open, int _b_gap_ext, int _match, int _mismatch,
		 bool reset_stats, int _anchor_width)
{

	dblen = _dblen;
	db = (int8_t *)malloc(dblen * sizeof(db[0]));
	if (db == NULL)
		return (1);

	qrlen = _qrlen;
	qr = (int8_t *)malloc(qrlen * sizeof(qr[0]));
	if (qr == NULL)
		return (1);

	swmatrix = (struct swcell *)malloc((dblen + 1) * (qrlen + 1) * sizeof(swmatrix[0]));
	if (swmatrix == NULL)
		return (1);

	backtrace = (int8_t *)malloc((dblen + qrlen) * sizeof(backtrace[0]));
	if (backtrace == NULL)
		return (1);

	dbalign = (char *)malloc((dblen + qrlen + 1) * sizeof(dbalign[0]));
	if (dbalign == NULL)
		return (1);

	qralign = (char *)malloc((dblen + qrlen + 1) * sizeof(dbalign[0]));
	if (qralign == NULL)
		return (1);

	a_gap_open = -(_a_gap_open);
	a_gap_ext = -(_a_gap_ext);
	b_gap_open = -(_b_gap_open);
	b_gap_ext = -(_b_gap_ext);
	match = _match;
	mismatch = _mismatch;

	if (reset_stats) {
	  swcells = swinvocs = 0;
	  sw_tc.type = DEF_FAST_TIME_COUNTER;
	  sw_tc.counter = 0;
	}

	anchor_width = _anchor_width;

	initialised = 1;

	return (0);
}

void
sw_full_ls_stats(uint64_t *invocs, uint64_t *cells, double *secs)
{
	
	if (invocs != NULL)
		*invocs = swinvocs;
	if (cells != NULL)
		*cells = swcells;
	if (secs != NULL)
	  *secs = time_counter_get_secs(&sw_tc);
}

void
sw_full_ls(uint32_t *genome, int goff, int glen, uint32_t *read, int rlen,
    int threshscore, int maxscore, struct sw_full_results *sfr, bool revcmpl,
    struct anchor * anchors, int anchors_cnt, int local_alignment)
{
	struct sw_full_results scratch;
	int i, j, k;

	//llint before = rdtsc(), after;
	TIME_COUNTER_START(sw_tc);

	if (!initialised)
		abort();

	swinvocs++;

	assert(glen > 0 && glen <= dblen);
	assert(rlen > 0 && rlen <= qrlen);

	if (sfr == NULL) {
	  sfr = &scratch;
	  memset(sfr, 0, sizeof(*sfr));
	}
	memset(backtrace, 0, (dblen + qrlen) * sizeof(backtrace[0]));

	dbalign[0] = qralign[0] = '\0';

	for (i = 0; i < glen; i++)
		db[i] = (int8_t)EXTRACT(genome, goff + i);

	for (i = 0; i < rlen; i++)
		qr[i] = (int8_t)EXTRACT(read, i);

	sfr->score = full_sw(glen, rlen, threshscore, maxscore, &i, &j,revcmpl, anchors, anchors_cnt,local_alignment);
	k = do_backtrace(glen, i, j, sfr);
	pretty_print(sfr->read_start, sfr->genome_start, k);
	sfr->gmapped = j - sfr->genome_start + 1;
	sfr->genome_start += goff;
	sfr->rmapped = i - sfr->read_start + 1;
	sfr->dbalign = xstrdup(dbalign);
	sfr->qralign = xstrdup(qralign);

	//swcells += (glen * rlen);
	//after = rdtsc();
	//swticks += MAX(after - before, 0);
	TIME_COUNTER_STOP(sw_tc);
}
