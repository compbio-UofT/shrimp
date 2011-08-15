/*	$Id: sw-full-cs.c,v 1.15 2009/06/16 23:26:21 rumble Exp $	*/

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
#include <limits.h>

#include <sys/types.h>
#include <sys/stat.h>
#include <sys/time.h>

#include "../common/util.h"
#include "../common/fasta.h"
#include "../common/sw-full-common.h"
#include "../common/sw-full-cs.h"
#include "../common/time_counter.h"

typedef struct swcell {
  struct {
    int	score_n;
    int	score_w;
    int	score_nw;

    int8_t	back_n;
    int8_t	back_w;
    int8_t	back_nw;
  } from[4];
} swcell;

#define FROM_A	0x00
#define FROM_B	0x01 
#define FROM_C	0x02
#define FROM_D	0x03

#define FROM_NORTH_NORTH		0x01
#define FROM_NORTH_NORTHWEST		0x02
#define FROM_WEST_NORTHWEST		0x03
#define FROM_WEST_WEST			0x04
#define FROM_NORTHWEST_NORTH		0x05
#define FROM_NORTHWEST_NORTHWEST	0x06
#define FROM_NORTHWEST_WEST		0x07

#define FROM_x(_mat, _dir)		(int8_t)(((_dir) << 2) | (_mat))

#define FROM_A_NORTH_NORTH		FROM_x(FROM_A, FROM_NORTH_NORTH)
#define FROM_A_NORTH_NORTHWEST		FROM_x(FROM_A, FROM_NORTH_NORTHWEST)
#define	FROM_A_WEST_NORTHWEST		FROM_x(FROM_A, FROM_WEST_NORTHWEST)
#define FROM_A_WEST_WEST		FROM_x(FROM_A, FROM_WEST_WEST)
#define FROM_A_NORTHWEST_NORTH		FROM_x(FROM_A, FROM_NORTHWEST_NORTH)
#define FROM_A_NORTHWEST_NORTHWEST	FROM_x(FROM_A, FROM_NORTHWEST_NORTHWEST)
#define FROM_A_NORTHWEST_WEST		FROM_x(FROM_A, FROM_NORTHWEST_WEST)

#define FROM_B_NORTH_NORTH		FROM_x(FROM_B, FROM_NORTH_NORTH)
#define FROM_B_NORTH_NORTHWEST		FROM_x(FROM_B, FROM_NORTH_NORTHWEST)
#define	FROM_B_WEST_NORTHWEST		FROM_x(FROM_B, FROM_WEST_NORTHWEST)
#define FROM_B_WEST_WEST		FROM_x(FROM_B, FROM_WEST_WEST)
#define FROM_B_NORTHWEST_NORTH		FROM_x(FROM_B, FROM_NORTHWEST_NORTH)
#define FROM_B_NORTHWEST_NORTHWEST	FROM_x(FROM_B, FROM_NORTHWEST_NORTHWEST)
#define FROM_B_NORTHWEST_WEST		FROM_x(FROM_B, FROM_NORTHWEST_WEST)

#define FROM_C_NORTH_NORTH		FROM_x(FROM_C, FROM_NORTH_NORTH)
#define FROM_C_NORTH_NORTHWEST		FROM_x(FROM_C, FROM_NORTH_NORTHWEST)
#define	FROM_C_WEST_NORTHWEST		FROM_x(FROM_C, FROM_WEST_NORTHWEST)
#define FROM_C_WEST_WEST		FROM_x(FROM_C, FROM_WEST_WEST)
#define FROM_C_NORTHWEST_NORTH		FROM_x(FROM_C, FROM_NORTHWEST_NORTH)
#define FROM_C_NORTHWEST_NORTHWEST	FROM_x(FROM_C, FROM_NORTHWEST_NORTHWEST)
#define FROM_C_NORTHWEST_WEST		FROM_x(FROM_C, FROM_NORTHWEST_WEST)

#define FROM_D_NORTH_NORTH		FROM_x(FROM_D, FROM_NORTH_NORTH)
#define FROM_D_NORTH_NORTHWEST		FROM_x(FROM_D, FROM_NORTH_NORTHWEST)
#define	FROM_D_WEST_NORTHWEST		FROM_x(FROM_D, FROM_WEST_NORTHWEST)
#define FROM_D_WEST_WEST		FROM_x(FROM_D, FROM_WEST_WEST)
#define FROM_D_NORTHWEST_NORTH		FROM_x(FROM_D, FROM_NORTHWEST_NORTH)
#define FROM_D_NORTHWEST_NORTHWEST	FROM_x(FROM_D, FROM_NORTHWEST_NORTHWEST)
#define FROM_D_NORTHWEST_WEST		FROM_x(FROM_D, FROM_NORTHWEST_WEST)

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
static int		a_gap_open, a_gap_ext, b_gap_open, b_gap_ext;
static int		match, mismatch;
static int		global_xover_penalty;
static struct swcell   *swmatrix;
static uint8_t	       *backtrace;
static char	       *dbalign, *qralign;
static int		anchor_width;
static int		indel_taboo_len;

/* statistics */
static uint64_t		swcells, swinvocs;
static time_counter	sw_tc;

#pragma omp threadprivate(initialised,db,qr,dblen,qrlen,a_gap_open,a_gap_ext,b_gap_open,b_gap_ext,match,mismatch,global_xover_penalty,\
			  swmatrix,backtrace,dbalign,qralign,sw_tc,swcells,swinvocs,indel_taboo_len)

#define BT_CROSSOVER		0x80
#define BT_CLIPPED		0xf0
#define BT_ISCROSSOVER(_x)	((_x) & BT_CROSSOVER)
#define BT_TYPE(_x)		((_x) & 0x0f)

#ifdef DEBUG_CROSSOVERS
static int _glen;
static int _rlen;
#endif


#ifdef DEBUG_SW
static void print_sw(int lena, int lenb) {
	printf("len a %d, len b %d\n",lena,lenb);
	int k;
	for (k=0; k<4; k++) {
		int i,j;
		printf("      %5s ","-");
		for (j=1; j< lena+1; j++) {
			printf("%5c ",base_translate(db[j-1],false)); 
		}
		printf("\n");
		//rows
		for (i=0; i<lenb+1; i++) {
			//cols
			if (i==0) {
				printf("    - ");
			} else {
				printf("%5c ",base_translate(qr[k][i-1],false));
			}
			for (j=0; j<lena+1; j++) {
				swcell curr=swmatrix[i*(lena+1)+j];
				int tmp=0;
				tmp=MAX(curr.from[k].score_n,curr.from[k].score_w);
				tmp=MAX(tmp,curr.from[k].score_nw);
				if (tmp<-99) {	
					tmp=-99;	
				} 
				printf("%5d ",tmp);
			}
			printf("\n");
		}
	}
}
#endif
/*
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
	1}
		printf("\n");
	}
}*/

inline static void
init_cell(int idx, int local_alignment, int xover_penalty) {
  if (local_alignment) {
	  swmatrix[idx].from[0].score_nw = 0;
	  swmatrix[idx].from[0].score_n  = -b_gap_open;
	  swmatrix[idx].from[0].score_w  = -a_gap_open;
	  swmatrix[idx].from[1].score_nw = xover_penalty;
	  swmatrix[idx].from[1].score_n  = -b_gap_open + xover_penalty;
	  swmatrix[idx].from[1].score_w  = -a_gap_open + xover_penalty;
	  swmatrix[idx].from[2].score_nw = xover_penalty;
	  swmatrix[idx].from[2].score_n  = -b_gap_open + xover_penalty;
	  swmatrix[idx].from[2].score_w  = -a_gap_open + xover_penalty;
	  swmatrix[idx].from[3].score_nw = xover_penalty;
	  swmatrix[idx].from[3].score_n  = -b_gap_open + xover_penalty;
	  swmatrix[idx].from[3].score_w  = -a_gap_open + xover_penalty;
  } else {
	  swmatrix[idx].from[0].score_nw = -INT_MAX/2;
	  swmatrix[idx].from[0].score_n  = -INT_MAX/2;
	  swmatrix[idx].from[0].score_w  = -INT_MAX/2;
	  swmatrix[idx].from[1].score_nw = -INT_MAX/2;
	  swmatrix[idx].from[1].score_n  = -INT_MAX/2;
	  swmatrix[idx].from[1].score_w  = -INT_MAX/2;
	  swmatrix[idx].from[2].score_nw = -INT_MAX/2;
	  swmatrix[idx].from[2].score_n  = -INT_MAX/2;
	  swmatrix[idx].from[2].score_w  = -INT_MAX/2;
	  swmatrix[idx].from[3].score_nw = -INT_MAX/2;
	  swmatrix[idx].from[3].score_n  = -INT_MAX/2;
	  swmatrix[idx].from[3].score_w  = -INT_MAX/2;

  }

  swmatrix[idx].from[0].back_nw = 0;
  swmatrix[idx].from[0].back_n  = 0;
  swmatrix[idx].from[0].back_w  = 0;
  swmatrix[idx].from[1].back_nw = 0;
  swmatrix[idx].from[1].back_n  = 0;
  swmatrix[idx].from[1].back_w  = 0;
  swmatrix[idx].from[2].back_nw = 0;
  swmatrix[idx].from[2].back_n  = 0;
  swmatrix[idx].from[2].back_w  = 0;
  swmatrix[idx].from[3].back_nw = 0;
  swmatrix[idx].from[3].back_n  = 0;
  swmatrix[idx].from[3].back_w  = 0;
}

/*
 * Perform a full Smith-Waterman alignment. For the colour case, this means
 * computing each possible letter space read string and doing a four layer
 * scan.
 */
//lena - genome length
static int
full_sw(int lena, int lenb, int threshscore, int *iret, int *jret,
	int *kret, bool revcmpl,
	struct anchor * anchors, int anchors_cnt, int local_alignment, int * crossover_score)
{
  int i, j, k, l, max_i, max_j, max_k;
  int score, ms, tmp, resetval, xover_penalty;
  //int go, ge;
  //int sw_band, ne_band;
  int8_t tmp2;
  struct anchor rectangle;


  /* shut up gcc */
  max_i = max_j = max_k = j = 0;

  score = 0;
  //go = gap_open;
  //ge = gap_ext;

  for (j = 0; j < lena + 1; j++) {
    init_cell(j, 1, global_xover_penalty);
  }
  //for (j = 0; j < lenb + 1; j++) {
  //init_cell(j * (lena + 1));
  //}

  /*
   * Figure out our band.
   *   We can actually skip computation of a significant number of
   *   cells, which could never be part of an alignment corresponding
   *   to our threshhold score.
   */
  //sw_band = ((lenb * match - threshscore + match - 1) / match) + 1;
  //ne_band = lena - (lenb - sw_band);

  if (anchors != NULL && anchor_width >= 0) {
    anchor_join(anchors, anchors_cnt, &rectangle);
    anchor_widen(&rectangle, anchor_width);
  } else {
    struct anchor tmp_anchors[2];

    tmp_anchors[0].x = 0;
    tmp_anchors[0].y = (lenb * match - threshscore) / match;
    tmp_anchors[0].length = 1;
    tmp_anchors[0].width = 1;

    tmp_anchors[1].x = lena - 1;
    tmp_anchors[1].y = lenb - 1 - tmp_anchors[0].y;
    tmp_anchors[1].length = 1;
    tmp_anchors[1].width = 1;

    anchor_join(tmp_anchors, 2, &rectangle);
  }



  for (i = 0; i < lenb; i++) {
    /*
     * computing row i of virtual matrix, stored in row i+1
     */
    int x_min, x_max;

    xover_penalty = (crossover_score == NULL? global_xover_penalty : crossover_score[i]);

    anchor_get_x_range(&rectangle, lena, lenb, i, &x_min, &x_max);
    if (!local_alignment) {
    	//x_max=MIN(lena,x_max); x_min=MAX(0,x_min-lenb/40); 
    	//init_cell(i * (lena + 1) + x_max  + 1,  0);
    	//init_cell((i + 1) * (lena + 1) + (x_min - 1) + 1, x_min == 0  ?  1 : 0);
      init_cell((i + 1) * (lena + 1) + (x_min - 1) + 1, 0, xover_penalty);
    } else {
      init_cell((i + 1) * (lena + 1) + (x_min - 1) + 1, 1, xover_penalty);
    }
    //if (x_min > 0) {
    //fprintf(stderr,"INIT cell %d , %d, %d\n",i+1, (x_min-1)+1,x_max);
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

      for (k = 0; k < 4; k++) {
	if (k != 0)
	  resetval = xover_penalty;
	else
	  resetval = 0;

	/*
	 * northwest
	 */
	if (db[j] == BASE_N || qr[k][i] == BASE_N)
	  ms = 0;
	else
	  ms = (db[j] == qr[k][i]) ? match : mismatch;

	if (!revcmpl) {
	  tmp  = cell_nw->from[k].score_nw + ms;
	  tmp2 = FROM_x(k, FROM_NORTHWEST_NORTHWEST);

	  // end of an insertion: not in taboo zone
	  if (i < lenb - indel_taboo_len && cell_nw->from[k].score_n + ms > tmp) {
	    tmp  = cell_nw->from[k].score_n + ms;
	    tmp2 = FROM_x(k, FROM_NORTHWEST_NORTH);
	  }

	  // end of a deletion
	  if (cell_nw->from[k].score_w + ms > tmp) {
	    tmp  = cell_nw->from[k].score_w + ms;
	    tmp2 = FROM_x(k, FROM_NORTHWEST_WEST);
	  }
	} else {
	  //end of a deletion
	  tmp  = cell_nw->from[k].score_w + ms;
	  tmp2 = FROM_x(k, FROM_NORTHWEST_WEST);

	  //end of an insertion: not in taboo zone
	  if (i < lenb - indel_taboo_len && cell_nw->from[k].score_n + ms > tmp) {
	    tmp  = cell_nw->from[k].score_n + ms;
	    tmp2 = FROM_x(k, FROM_NORTHWEST_NORTH);
	  }

	  if (cell_nw->from[k].score_nw + ms > tmp) {
	    tmp  = cell_nw->from[k].score_nw + ms;
	    tmp2 = FROM_x(k, FROM_NORTHWEST_NORTHWEST);
	  }
	}

	/* check neighbours */
	for (l = 0; l < 4; l++) {
	  if (l == k)
	    continue;

	  if (!revcmpl) {
	    /* northwest */
	    if (cell_nw->from[l].score_nw + ms + xover_penalty > tmp) {
	      tmp  = cell_nw->from[l].score_nw + ms + xover_penalty;
	      tmp2 = FROM_x(l, FROM_NORTHWEST_NORTHWEST);
	    }

	    /* north */ // end of insertion, not in taboo zone
	    if (i < lenb - indel_taboo_len && cell_nw->from[l].score_n + ms + xover_penalty > tmp) {
	      tmp  = cell_nw->from[l].score_n + ms + xover_penalty;
	      tmp2 = FROM_x(l, FROM_NORTHWEST_NORTH);
	    }

	    /* west */
	    if (cell_nw->from[l].score_w + ms + xover_penalty > tmp) {
	      tmp  = cell_nw->from[l].score_w + ms + xover_penalty;
	      tmp2 = FROM_x(l, FROM_NORTHWEST_WEST);
	    }
	  } else {
	    /* west */
	    if (cell_nw->from[l].score_w + ms + xover_penalty > tmp) {
	      tmp  = cell_nw->from[l].score_w + ms + xover_penalty;
	      tmp2 = FROM_x(l, FROM_NORTHWEST_WEST);
	    }

	    /* north */ // end of insertion, not in taboo zone
	    if (i < lenb - indel_taboo_len && cell_nw->from[l].score_n + ms + xover_penalty > tmp) {
	      tmp  = cell_nw->from[l].score_n + ms + xover_penalty;
	      tmp2 = FROM_x(l, FROM_NORTHWEST_NORTH);
	    }

	    /* northwest */
	    if (cell_nw->from[l].score_nw + ms + xover_penalty > tmp) {
	      tmp  = cell_nw->from[l].score_nw + ms + xover_penalty;
	      tmp2 = FROM_x(l, FROM_NORTHWEST_NORTHWEST);
	    }
	  }
	}

	if (tmp <= resetval && local_alignment) {
	  tmp = resetval;
	  tmp2 = 0;
	}

	cell_cur->from[k].score_nw = tmp;
	cell_cur->from[k].back_nw  = tmp2;


	/*
	 * north
	 */
	if (!revcmpl) {
	  // insertion start
	  tmp  = cell_n->from[k].score_nw - b_gap_open - b_gap_ext;
	  tmp2 = FROM_x(k, FROM_NORTH_NORTHWEST);

	  if (!(i < lenb - indel_taboo_len) || cell_n->from[k].score_n - b_gap_ext > tmp) {
	    tmp  = cell_n->from[k].score_n - b_gap_ext;
	    tmp2 = FROM_x(k, FROM_NORTH_NORTH);
	  }
	} else {
	  tmp  = cell_n->from[k].score_n - b_gap_ext;
	  tmp2 = FROM_x(k, FROM_NORTH_NORTH);

	  // insertion start
	  if (i < lenb - indel_taboo_len && cell_n->from[k].score_nw - b_gap_open - b_gap_ext > tmp) {
	    tmp  = cell_n->from[k].score_nw - b_gap_open - b_gap_ext;
	    tmp2 = FROM_x(k, FROM_NORTH_NORTHWEST);
	  }
	}

	/* check neighbours */
	for (l = 0; l < 4; l++) {
	  if (l == k)
	    continue;

	  if (!revcmpl) {
	    /* northwest */ // insertion start
	    if (i < lenb - indel_taboo_len && cell_n->from[l].score_nw - b_gap_open - b_gap_ext + xover_penalty > tmp) {
	      tmp  = cell_n->from[l].score_nw - b_gap_open - b_gap_ext + xover_penalty;
	      tmp2 = FROM_x(l, FROM_NORTH_NORTHWEST);
	    }

	    /* north */
	    if (cell_n->from[l].score_n - b_gap_ext + xover_penalty > tmp) {
	      tmp  = cell_n->from[l].score_n - b_gap_ext + xover_penalty;
	      tmp2 = FROM_x(l, FROM_NORTH_NORTH);
	    }
	  } else {
	    /* north */
	    if (cell_n->from[l].score_n - b_gap_ext + xover_penalty > tmp) {
	      tmp  = cell_n->from[l].score_n - b_gap_ext + xover_penalty;
	      tmp2 = FROM_x(l, FROM_NORTH_NORTH);
	    }

	    /* northwest */ // insertion start
	    if (i < lenb - indel_taboo_len && cell_n->from[l].score_nw - b_gap_open - b_gap_ext + xover_penalty > tmp) {
	      tmp  = cell_n->from[l].score_nw - b_gap_open - b_gap_ext + xover_penalty;
	      tmp2 = FROM_x(l, FROM_NORTH_NORTHWEST);
	    }
	  }
	}

	if (tmp <= resetval && local_alignment) {
	  tmp = resetval;
	  tmp2 = 0;
	}
					
	cell_cur->from[k].score_n = tmp;
	cell_cur->from[k].back_n  = tmp2;

				
	/*
	 * west
	 */
	if (!revcmpl) {
	  // deletion start
	  tmp  = cell_w->from[k].score_nw - a_gap_open - a_gap_ext;
	  tmp2 = FROM_x(k, FROM_WEST_NORTHWEST);

	  if (!(i < lenb - indel_taboo_len) || cell_w->from[k].score_w - a_gap_ext > tmp) {
	    tmp  = cell_w->from[k].score_w - a_gap_ext;
	    tmp2 = FROM_x(k, FROM_WEST_WEST);
	  }
	} else {
	  tmp  = cell_w->from[k].score_w - a_gap_ext;
	  tmp2 = FROM_x(k, FROM_WEST_WEST);

	  // deletion start
	  if (i < lenb - indel_taboo_len && cell_w->from[k].score_nw - a_gap_open - a_gap_ext > tmp) {
	    tmp  = cell_w->from[k].score_nw - a_gap_open - a_gap_ext;
	    tmp2 = FROM_x(k, FROM_WEST_NORTHWEST);
	  }
	}

	/*
	 * NB: It doesn't make sense to cross over on a
	 *     genomic gap, so we won't.
	 */

	if (tmp <= resetval && local_alignment) {
	  tmp = resetval;
	  tmp2 = 0;
	}

	cell_cur->from[k].score_w = tmp;
	cell_cur->from[k].back_w  = tmp2;


	/*
	 * max score
	 */
	if (local_alignment || i==lenb-1) {
	if (!revcmpl) {
	  if (cell_cur->from[k].score_nw > score) {
	    score = cell_cur->from[k].score_nw;
	    max_i = i, max_j = j, max_k = k;
	  }
	  if (cell_cur->from[k].score_n > score) {
	    score = cell_cur->from[k].score_n;
	    max_i = i, max_j = j, max_k = k;
	  }
	  if (cell_cur->from[k].score_w > score) {
	    score = cell_cur->from[k].score_w;
	    max_i = i, max_j = j, max_k = k;
	  }
	} else {
	  if (cell_cur->from[k].score_w > score) {
	    score = cell_cur->from[k].score_w;
	    max_i = i, max_j = j, max_k = k;
	  }
	  if (cell_cur->from[k].score_n > score) {
	    score = cell_cur->from[k].score_n;
	    max_i = i, max_j = j, max_k = k;
	  }
	  if (cell_cur->from[k].score_nw > score) {
	    score = cell_cur->from[k].score_nw;
	    max_i = i, max_j = j, max_k = k;
	  }
	}
	}

#ifdef DEBUG_SW
	fprintf(stderr, "i:%d j:%d k:%d score_nw:%d [%u,%s] score_n:%d [%u,%s] score_w:%d [%u,%s] xover_penalty:%d\n", i+1, j+1, k,
		cell_cur->from[k].score_nw, cell_cur->from[k].back_nw & 0x3,
		(cell_cur->from[k].back_nw >> 2 == 0 ? "!" :
		 (cell_cur->from[k].back_nw >> 2 == FROM_NORTHWEST_NORTH ? "n" :
		  (cell_cur->from[k].back_nw >> 2 == FROM_NORTHWEST_NORTHWEST ? "nw" : "w"))),

		cell_cur->from[k].score_n, cell_cur->from[k].back_n & 0x3,
		(cell_cur->from[k].back_n >> 2 == 0 ? "!" :
		 (cell_cur->from[k].back_n >> 2 == FROM_NORTH_NORTH ? "n" : "nw")),

		cell_cur->from[k].score_w, cell_cur->from[k].back_w & 0x3,
		(cell_cur->from[k].back_w >> 2 == 0 ? "!" :
		 (cell_cur->from[k].back_w >> 2 == FROM_WEST_NORTHWEST ? "nw" : "w")),

		xover_penalty);
#endif

      }
    }

    if (i+1 < lenb) {
      int next_x_min, next_x_max;

      anchor_get_x_range(&rectangle, lena, lenb, i+1, &next_x_min, &next_x_max);
      for (j = x_max + 1; j <= next_x_max; j++) {
	//fprintf(stderr,"Init cell %d , , %d\n",i+1,j+1);
	init_cell((i + 1) * (lena + 1) + (j + 1), local_alignment, xover_penalty); // still xover on i-th color
      }
    }
  }

#ifdef DEBUG_SW
  fprintf(stderr, "max_i:%d max_j:%d max_k:%d\n", max_i+1, max_j+1, max_k);
#endif

  *iret = max_i;
  *jret = max_j;
  *kret = max_k;
  //print_sw(lena,lenb);
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
  struct swcell *cell;
  int off, from, fromscore;

  off = (dblen + qrlen) - 1;

  cell = &swmatrix[(i + 1) * (lena + 1) + j + 1];

  from = cell->from[k].back_nw;
  fromscore = cell->from[k].score_nw;

  if (cell->from[k].score_w > fromscore) {
    from = cell->from[k].back_w;
    fromscore = cell->from[k].score_w;
  }
  if (cell->from[k].score_n > fromscore)
    from = cell->from[k].back_n;

  if (from == 0) {
    int l, base;
    fprintf(stderr, "Assertion failed.\nQr:");
    for (l = 1, base = qr[0][0]; l < qrlen; l++) {
      fprintf(stderr, "%d", lstocs(base, qr[0][l], false));
      base = qr[0][l];
    }
    fprintf(stderr, "\n");
  }
  assert(from != 0);

  /* fill out the backtrace */
  while (i >= 0 && j >= 0) {
	//printf("i %d, j %d\n",i,j);
    assert(off >= 0);

    cell = NULL;

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
      if (db[j] == qr[k][i] || db[j] == BASE_N || qr[k][i] == BASE_N)
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

    /* handle match/mismatch and north */
    switch (from) {
    case FROM_A_NORTH_NORTH:
    case FROM_A_NORTH_NORTHWEST:
    case FROM_B_NORTH_NORTH:
    case FROM_B_NORTH_NORTHWEST:
    case FROM_C_NORTH_NORTH:
    case FROM_C_NORTH_NORTHWEST:
    case FROM_D_NORTH_NORTH:
    case FROM_D_NORTH_NORTHWEST:
      switch(k) {
      case 0:
	backtrace[off]= BACK_A_DELETION;
	break;
      case 1:
	backtrace[off]= BACK_B_DELETION;
	break;
      case 2:
	backtrace[off]= BACK_C_DELETION;
	break;
      case 3:
	backtrace[off]= BACK_D_DELETION;
	break;
      default:
	fprintf(stderr, "INTERNAL ERROR: k = %d\n", k);
	assert(0);
      }
      break;

    case FROM_A_WEST_WEST:
    case FROM_A_WEST_NORTHWEST:
    case FROM_B_WEST_WEST:
    case FROM_B_WEST_NORTHWEST:
    case FROM_C_WEST_WEST:
    case FROM_C_WEST_NORTHWEST:
    case FROM_D_WEST_WEST:
    case FROM_D_WEST_NORTHWEST:
      /* doesn't make sense to cross over on a genomic gap */
      backtrace[off] = BACK_INSERTION;
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
      switch(k) {
      case 0:
	backtrace[off] = BACK_A_MATCH_MISMATCH;
	break;
      case 1:
	backtrace[off] = BACK_B_MATCH_MISMATCH;
	break;
      case 2:
	backtrace[off] = BACK_C_MATCH_MISMATCH;
	break;
      case 3:
	backtrace[off] = BACK_D_MATCH_MISMATCH;
	break;
      default:
	fprintf(stderr, "INTERNAL ERROR: k = %d\n", k);
	assert(0);
      }
      break;

    default:
      fprintf(stderr, "INTERNAL ERROR: from = %d\n", from);
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

    default:
      fprintf(stderr, "INTERNAL ERROR: from = %d\n", from);
      assert(0);
    }


    /*
     * Continue backtrace (nb: i,j and k have already been changed).
     */
    cell = &swmatrix[(i + 1) * (lena + 1) + j + 1];

    switch (from) {
    case FROM_A_NORTH_NORTH:
    case FROM_B_NORTH_NORTH:
    case FROM_C_NORTH_NORTH:
    case FROM_D_NORTH_NORTH:
      from = cell->from[k].back_n;
      break;

    case FROM_A_NORTH_NORTHWEST:
    case FROM_B_NORTH_NORTHWEST:
    case FROM_C_NORTH_NORTHWEST:
    case FROM_D_NORTH_NORTHWEST:
      from = cell->from[k].back_nw;
      break;

    case FROM_A_WEST_WEST:
    case FROM_B_WEST_WEST:
    case FROM_C_WEST_WEST:
    case FROM_D_WEST_WEST:
      from = cell->from[k].back_w;
      break;

    case FROM_A_WEST_NORTHWEST:
    case FROM_B_WEST_NORTHWEST:
    case FROM_C_WEST_NORTHWEST:
    case FROM_D_WEST_NORTHWEST:
      from = cell->from[k].back_nw;
      break;

    case FROM_A_NORTHWEST_NORTH:
    case FROM_B_NORTHWEST_NORTH:
    case FROM_C_NORTHWEST_NORTH:
    case FROM_D_NORTHWEST_NORTH:
      from = cell->from[k].back_n;
      break;

    case FROM_A_NORTHWEST_NORTHWEST:
    case FROM_B_NORTHWEST_NORTHWEST:
    case FROM_C_NORTHWEST_NORTHWEST:
    case FROM_D_NORTHWEST_NORTHWEST:
      from = cell->from[k].back_nw;
      break;

    case FROM_A_NORTHWEST_WEST:
    case FROM_B_NORTHWEST_WEST:
    case FROM_C_NORTHWEST_WEST:
    case FROM_D_NORTHWEST_WEST:
      from = cell->from[k].back_w;
      break;

    default:
      fprintf(stderr, "INTERNAL ERROR: from = %d\n", from);
      assert(0);
    }

    off--;

    if (from == 0)
      break;		  
  }

  off++;

  if (k != 0) {
    backtrace[off] |= BT_CROSSOVER;
    sfr->crossovers++;
  }

  return (off);
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

  d = dbalign;
  q = qralign;

  for (l = k; l < (dblen + qrlen); l++) {

#ifdef DEBUG_CROSSOVERS
    int a;
    if (BT_ISCROSSOVER(backtrace[l])
	&& (BT_TYPE(backtrace[l]) == BACK_A_DELETION
	    || BT_TYPE(backtrace[l]) == BACK_B_DELETION
	    || BT_TYPE(backtrace[l]) == BACK_C_DELETION
	    || BT_TYPE(backtrace[l]) == BACK_D_DELETION)) {
      fprintf(stderr, "sw-full-cs: crossover in \"deletion\" (really, insertion):\n");
      fprintf(stderr, "db:");
      for (a = 0; a < _glen; a++)
	fprintf(stderr, "%c", base_translate(db[a], false));
      fprintf(stderr, "\n");
      fprintf(stderr, "qr[0]:");
      for (a = 0; a < _rlen; a++)
	fprintf(stderr, "%c", base_translate(qr[0][a], false));
      fprintf(stderr, "\n");
    }
#endif

    switch (BT_TYPE(backtrace[l])) {
    case BACK_A_DELETION:
      *d++ = '-';
      if (BT_ISCROSSOVER(backtrace[l]))
	*q++ = (char)tolower((int)base_translate(qr[0][i++], false));
      else	
	*q++ = base_translate(qr[0][i++], false);
      break;

    case BACK_B_DELETION:
      *d++ = '-';
      if (BT_ISCROSSOVER(backtrace[l]))
	*q++ = (char)tolower((int)base_translate(qr[1][i++], false));
      else	
	*q++ = base_translate(qr[1][i++], false);
      break;

    case BACK_C_DELETION:
      *d++ = '-';
      if (BT_ISCROSSOVER(backtrace[l]))
	*q++ = (char)tolower((int)base_translate(qr[2][i++], false));
      else	
	*q++ = base_translate(qr[2][i++], false);
      break;

    case BACK_D_DELETION:
      *d++ = '-';
      if (BT_ISCROSSOVER(backtrace[l]))
	*q++ = (char)tolower((int)base_translate(qr[3][i++], false));
      else	
	*q++ = base_translate(qr[3][i++], false);
      break;

    case BACK_INSERTION:
      *d++ = base_translate(db[j++], false);
      *q++ = '-';
      break;

    case BACK_A_MATCH_MISMATCH:
      *d++ = base_translate(db[j++], false);
      if (BT_ISCROSSOVER(backtrace[l]))
	*q++ = (char)tolower((int)base_translate(qr[0][i++], false));
      else	
	*q++ = base_translate(qr[0][i++], false);
      break;

    case BACK_B_MATCH_MISMATCH:
      *d++ = base_translate(db[j++], false);
      if (BT_ISCROSSOVER(backtrace[l]))
	*q++ = (char)tolower((int)base_translate(qr[1][i++], false));
      else	
	*q++ = base_translate(qr[1][i++], false);
      break;

    case BACK_C_MATCH_MISMATCH:
      *d++ = base_translate(db[j++], false);
      if (BT_ISCROSSOVER(backtrace[l]))
	*q++ = (char)tolower((int)base_translate(qr[2][i++], false));
      else	
	*q++ = base_translate(qr[2][i++], false);
      break;

    case BACK_D_MATCH_MISMATCH:
      *d++ = base_translate(db[j++], false);
      if (BT_ISCROSSOVER(backtrace[l]))
	*q++ = (char)tolower((int)base_translate(qr[3][i++], false));
      else	
	*q++ = base_translate(qr[3][i++], false);
      break;
	
    default:
      fprintf(stderr, "INTERNAL ERROR: backtrace[l] = 0x%x\n", backtrace[l]);
      assert(0);
    }
    if ((BT_TYPE(backtrace[l]) == BACK_A_MATCH_MISMATCH || BT_TYPE(backtrace[l]) == BACK_B_MATCH_MISMATCH
	 || BT_TYPE(backtrace[l]) == BACK_C_MATCH_MISMATCH || BT_TYPE(backtrace[l]) == BACK_D_MATCH_MISMATCH)
	&& (*(q-1) == 'n' || *(q-1) == 'N')) {
      if (BT_ISCROSSOVER(backtrace[l]))
	*(q-1) = (char)tolower(*(d-1));
      else
	*(q-1) = *(d-1);
    }
  }

  *d = *q = '\0';
}

int
sw_full_cs_cleanup(void) {
	free(db);
	int i;
	for (i=0; i<4; i++) {
		free(qr[i]);
	}
	free(swmatrix);
	free(backtrace);
	free(dbalign);
	free(qralign);
	return 0;
}

int
sw_full_cs_setup(int _dblen, int _qrlen, int _a_gap_open, int _a_gap_ext, int _b_gap_open, int _b_gap_ext,
		 int _match, int _mismatch, int _global_xover_penalty, bool reset_stats,
		 int _anchor_width, int _indel_taboo_len)
{
  int i;

  dblen = _dblen;
  db = (int8_t *)malloc(dblen * sizeof(db[0]));
  if (db == NULL)
    return (1);

  qrlen = _qrlen;
  for (i = 0; i < 4; i++) {
    qr[i] = (int8_t *)malloc(qrlen * sizeof(qr[0]));
    if (qr[i] == NULL)
      return (1);
  }

  swmatrix = (struct swcell *)malloc((dblen + 1) * (qrlen + 1) *
				     sizeof(swmatrix[0]));
  if (swmatrix == NULL)
    return (1);

  backtrace = (uint8_t *)malloc((dblen + qrlen) * sizeof(backtrace[0]));
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
  global_xover_penalty = _global_xover_penalty;

  if (reset_stats) {
    swcells = swinvocs = 0;
    sw_tc.type = DEF_FAST_TIME_COUNTER;
    sw_tc.counter = 0;
  }

  anchor_width = _anchor_width;
  indel_taboo_len = _indel_taboo_len;

  initialised = 1;

  return (0);
}

void
sw_full_cs_stats(uint64_t *invocs, uint64_t *cells, double *secs)
{
	
  if (invocs != NULL)
    *invocs = swinvocs;
  if (cells != NULL)
    *cells = swcells;
  if (secs != NULL)
    *secs = time_counter_get_secs(&sw_tc);
}

void
sw_full_cs(uint32_t *genome_ls, int goff, int glen, uint32_t *read, int rlen,
	   int initbp, int threshscore, struct sw_full_results *sfr, bool revcmpl, bool is_rna,
	   struct anchor * anchors, int anchors_cnt, int local_alignment, int * crossover_score)
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
    db[i] = (int8_t)EXTRACT(genome_ls, goff + i);

  /*
   * Generate each possible letter space sequence from the colour space
   * read. qr[0] corresponds to initbp, which is given initial preference.
   */
  assert(initbp >= 0 && initbp <= 3);
  for (i = 0; i < 4; i++) {
    int letter = (i + initbp) % 4;

    for (j = 0; j < rlen; j++) {
      int base = EXTRACT(read, j);

      if (base == BASE_N || base == BASE_X) {
	qr[i][j] = BASE_N;
	letter = (i + initbp) % 4;
      } else {
	qr[i][j] = (int8_t)cstols(letter, base, is_rna);
	letter = qr[i][j];
      }
    }
  }

#ifdef DEBUG_SW
  fprintf(stderr, "db: ");
  for (j = 0; j < glen; j++)
    fprintf(stderr, "%c", base_translate(db[j], false));
  fprintf(stderr, "\n");
  for (i = 0; i < 4; i++) {
    fprintf(stderr, "qr[%u]: ", i);
    for (j = 0; j < rlen; j++)
      fprintf(stderr, "%c", base_translate(qr[i][j], false));
    fprintf(stderr, "\n");
  }
#endif

#ifdef DEBUG_CROSSOVERS
  _glen = glen;
  _rlen = rlen;
#endif

  sfr->score = full_sw(glen, rlen, threshscore, &i, &j, &k, revcmpl, anchors, anchors_cnt, local_alignment, crossover_score);
  if (sfr->score >= 0 && sfr->score >= threshscore) {
    k = do_backtrace(glen, i, j, k, sfr);
    pretty_print(sfr->read_start, sfr->genome_start, k);
    sfr->gmapped = j - sfr->genome_start + 1;
    sfr->genome_start += goff;
    sfr->rmapped = i - sfr->read_start + 1;
    sfr->dbalign = xstrdup(dbalign);
    sfr->qralign = xstrdup(qralign);
  } else {
    sfr->score = 0;
  }

#ifdef DEBUG_SW
  fprintf(stderr, "reported alignment:\n\t%s\n\t%s\n", sfr->dbalign, sfr->qralign);
#endif

  //swcells += (glen * rlen);
  //after = rdtsc();
  //swticks += MAX(after - before, 0);
  TIME_COUNTER_STOP(sw_tc);
}
