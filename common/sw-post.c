/*
 * In this module, we are given the results of a full SW run,
 * and we compute two things:
 *
 * 1. The probability the location produced the read (over all possible alignments).
 *    Currently, we only sum over all alignments respecting the current gaps.
 *    As a result, this is useless to do in letter space, where the mapper
 *    cannot distinguish between errors and SNPs.
 *
 * 2. For each output letter, the probability that it is correct.
 *    Only for color space. In letter space, this is given by the base quality value.
 *
 * Both are computed as scores.
 */

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

#include "../common/fasta.h"
#include "../common/util.h"
#include "../common/sw-post.h"
#include "../common/sw-full-common.h"
#include "../common/time_counter.h"


static int	initialized;

static double	pr_snp;
static double	pr_xover;
static double	pr_del_open;
static double	pr_del_extend;
static double	pr_ins_open;
static double	pr_ins_extend;

static bool	use_read_qvs;
static bool	use_sanger_qvs;
static int	default_qual;		// if no qvs, use this instead
static int	qual_vector_offset;	// i.e. is there a useless qv for the initial base in cs?
static int	qual_delta;		// how much to subtract from chars to get the int

static int	init_bp;
static int	len;

//static double	neglogsixteenth;
//static double	neglogfourth;

typedef struct column{
  double forwards[16]; //we'll misuse this for the viterbi
  double backwards[16];
  double forwscale; //adding numerical stability
  double backscale; //adding numerical stability
  int ncols;
  int nlets;
  int letssize;
  int colssize;
  int* lets;
  int* cols;
  double* letserrrate;
  double* colserrrate;
  char backpointer[16]; //previous state for viterbi
  double posterior[4];
  int max_posterior;
  int base_call;
} states;

static struct column *	columns;
static int max_len;

static uint64_t		cells, invocs;
static time_counter	tc;

static int check;

#pragma omp threadprivate(initialized,\
			  pr_snp,pr_xover,pr_del_open,pr_del_extend,pr_ins_open,pr_ins_extend,\
			  use_read_qvs,use_sanger_qvs,default_qual,qual_vector_offset,qual_delta,\
			  init_bp,len,columns,max_len,\
			  tc,cells,invocs,check)


/*********************************************************************************
 *
 * BEGIN Forward-backward code
 *
 *********************************************************************************/

#define left(i) ( ((i) >> 2) & 3)
#define right(i) ( (i) & 3)
#define MIN2(i,j) ( ((i)< (j))? (i):(j))

/* In order to understand any of the code below, you need to understand color-space;
   Specifically that LETTER ^ LETTER = COLOR : T (00) ^ C (10) = 2 (10), etc. 
   And that LETTER ^ COLOR = NEXTLETTER: T (00) ^ 3 (11) = A (11). */

/* compute prior probability of the node given the emissions. letters are thought to be at the
   left "side" of the pair emitted by the node */

double nodePrior(states* allstates, int i, int j) { //i  is state, j is node in the state
  double val = 0;
  double errrate;
  int let, col, k;
  for (k = 0; k < allstates[i].nlets; k++) {
    let = allstates[i].lets[k];
    errrate = allstates[i].letserrrate[k];
    if (right(j) == let) {
      val = val - log(1-errrate);
    }
    else {
      val = val - log(errrate/3.0);
    }
  }
  //fprintf(stderr, "nodeprior: %g", val);
  for (k = 0; k < allstates[i].ncols; k++) {
    col = allstates[i].cols[k];
    errrate = allstates[i].colserrrate[k];
    if ((left(j) ^ right(j)) == col) {
      val = val - log(1-errrate);
    }
    else {
      val = val - log(errrate/3.0);
    }
    //fprintf(stderr, " %g\n", val);
  }
  return val;
}

/* Little helper for debugging */

void printStates(states* allstates, int stateslen, FILE* stream) {
  int i,j,k;
  fprintf(stream, "\nCONTIG %d", stateslen);

  for (i=0; i< stateslen; i++) {
    fprintf(stream, "\nCOLORS[%d] ",i);
    for (k = 0; k < allstates[i].ncols; k++) {
      fprintf(stream, "%d (%g)",allstates[i].cols[k], allstates[i].colserrrate[k]);
    }
  }
  for (i=0; i< stateslen; i++) {
    fprintf(stream, "\nFORWARDSS[%d] ",i);
    for (j=0; j< 16; j++) {
      fprintf(stream, "%.5g ",allstates[i].forwards[j] + allstates[i].forwscale);
    }    
  }
  for (i=0; i< stateslen; i++) {
    fprintf(stream, "\nBACKWARDSS[%d] ",i);
    for (j=0; j< 16; j++) {
      fprintf(stream, "%.5g ",allstates[i].backwards[j] + allstates[i].backscale);
    }    
  }

  for (i=0; i< stateslen; i++) {
    fprintf(stream, "\nLETS[%d] ",i);
    for (k = 0; k < allstates[i].nlets; k++) {
      fprintf(stream, "%d ",allstates[i].lets[k]);
    }
  
    fprintf(stream, "%c",base_to_char(allstates[i].max_posterior, LETTER_SPACE));
    fprintf(stream, " %.5g %.5g %.5g %.5g",
	    allstates[i].posterior[0],allstates[i].posterior[1],allstates[i].posterior[2],allstates[i].posterior[3]);
  }
  
  fprintf(stream, "\n");
}

/*maximum posterior traceback */

void post_traceback (states* allstates, int stateslen, double norm_px) {
  int i = 0, j, maxval;

  for (i = 0; i < stateslen; i++) {
    for (j=0; j< 4; j++) allstates[i].posterior[j] = 0;
    for (j = 0; j < 16; j++) {
      //     fprintf(stderr, "%g %g %g\n",allstates[i].forwards[j], allstates[i].backwards[j], norm_px); 
      allstates[i].posterior[right(j)] += exp(-1 * (allstates[i].forwards[j] + allstates[i].backwards[j] + allstates[i].forwscale + allstates[i].backscale - norm_px)); 
      //      fprintf(stderr, "distrib[%d,%d] = %g\n", i,j, exp(-1 * (allstates[i].forwards[j] + allstates[i].backwards[j] + allstates[i].forwscale + allstates[i].backscale - norm_px)));
    }
    maxval = 0;
    for (j=1; j< 4; j++)  {
      //      fprintf(stderr, "let_distrib[%d,%d] = %g\n", i,j, distrib[j]);

      if (allstates[i].posterior[j] >allstates[i].posterior[maxval]) 
	maxval = j;
    }
    //    fprintf (stderr, "\n");
    //if (allstates[i].posterior[maxval] > confrate) {
    allstates[i].max_posterior = maxval;
    //}
    //else {
    //  allstates[i].max_posterior = BASE_N;
    //}
  }

}


/*viterbi traceback */
/*
char* vit_traceback (states* allstates, int stateslen) {
  char* result = (char*) calloc (stateslen + 1, 1);
  int i,j;
  int minval, prev;

  assert(0); // not changed to right letter emission

  for (i = stateslen -1; i >= 0; i--) {
    minval = 0;
    for (j = 0; j< 16; j++) {
      if (allstates[i].forwards[j] < allstates[i].forwards[minval]) {
	minval = j;
      }
    }
    prev = allstates[i].backpointer[minval];
    if (i && (left(minval) != right (prev))) {
      fprintf (stderr, "BACKTRACE error %d %d %d\n", i, minval, prev);
      exit(2);
    }
    result[i] = letmap[left(minval)];
  }
  return result;
}


void viterbi (states* allstates, int stateslen) {
  int i,j,k,let,col;
  int minback;
  double valback;
  double val;

  assert(0); // not changed to right letter emission

  i = 0;
  for (j = 0; j < 16; j++) {
    allstates[i].forwards[j] = nodePrior(allstates,i,j);
  }

  for (i=1; i < stateslen; i++) {    
    for (j = 0; j < 16; j++) {
      allstates[i].forwards[j] = nodePrior(allstates,i,j);

      minback = left(j);
      for (k = 1; k < 16; k++) {
	if (left(j) == right(k)) {
	  if (allstates[i-1].forwards[k] < allstates[i-1].forwards[minback]) {
	    minback = k;
	  }
	}
      }
      valback = allstates[i-1].forwards[minback];
      allstates[i].forwards[j] += valback;
      allstates[i].backpointer[j] = minback;
    }
  }
}
*/

double do_backwards (states* allstates, int stateslen) {
  int i,j,k; //,let,col;
  double val;
  
  i = stateslen-1;
  allstates[i].backscale = 999999999;
  for (j = 0; j < 16; j++) {
    allstates[i].backwards[j] = 0; // matei change: bug fix
    allstates[i].backscale = MIN2 (allstates[i].backscale, allstates[i].backwards[j]);
  }
  for (j = 0; j < 16; j++) {
    allstates[i].backwards[j] -= allstates[i].backscale;
  }

  for (i = stateslen-2; i >=0; i--) {    
    allstates[i].backscale = 999999999;
    memset(allstates[i].backwards, 0, 16 * sizeof(allstates[i].backwards[0])); // matei: bug fix
    for (j = 0; j < 16; j++) {
      for (k = 0; k < 16; k++) {
	if (right(j) == left(k)) {
	  val = nodePrior(allstates,i+1,k);
	  allstates[i].backwards[j] += exp(-1*(val + allstates[i+1].backwards[k]));
	}
      }
      //      fprintf(stdout, "bw was [%d, %d] = %g\n", i, j, allstates[i].backwards[j]); 

      allstates[i].backwards[j] = -log(allstates[i].backwards[j]); // + neglogfourth;
      allstates[i].backscale = MIN2 (allstates[i].backscale, allstates[i].backwards[j]);
    }
    for (j = 0; j < 16; j++) {

      allstates[i].backwards[j] -= allstates[i].backscale;
      //      fprintf(stdout, "bw is [%d, %d] = %g\n", i, j, allstates[i].backwards[j]); 
    }
    allstates[i].backscale += allstates[i+1].backscale;

  }
  val = 0;
  i = 0;
  for (j = 0; j < 16; j++) {
    if (left(j) == init_bp) { // matei change: second letter emission
      val += exp(-1*(allstates[i].backwards[j] + nodePrior(allstates,i,j))); // + neglogfourth));
    }
  }
  return -log(val) + allstates[0].backscale;
}

double do_forwards (states* allstates, int stateslen) {
  int i,j,k; //,let,col;
  double val;
  
  i = 0; j = 0;
  allstates[i].forwscale = 999999999;
  for (j = 0; j < 16; j++) {
    if (left(j) == init_bp) { // matei change: second letter emission
      allstates[i].forwards[j] = nodePrior(allstates,i,j); // + neglogfourth;
      allstates[i].forwscale = MIN2 (allstates[i].forwscale, allstates[i].forwards[j]);
    } else {
      allstates[i].forwards[j] = HUGE_VAL;
    }
  }
  for (j = 0; j < 16; j++) {
    allstates[i].forwards[j] -= allstates[i].forwscale;
  }

  for (i=1; i < stateslen; i++) {    
    allstates[i].forwscale = 999999999;
    memset(allstates[i].forwards, 0, 16 * sizeof(allstates[i].forwards[0])); // matei: bug fix
    for (j = 0; j < 16; j++) {
      val = nodePrior(allstates,i,j);
      for (k = 0; k < 16; k++) {
	if (left(j) == right(k)) {
	  allstates[i].forwards[j] += exp(-1*(allstates[i-1].forwards[k]));
	}
      }
      allstates[i].forwards[j] = val - log(allstates[i].forwards[j]); //+ neglogfourth;
      allstates[i].forwscale = MIN2 (allstates[i].forwscale, allstates[i].forwards[j]);
    }
    for (j = 0; j < 16; j++) {
      allstates[i].forwards[j] -= allstates[i].forwscale;
    }
    allstates[i].forwscale += allstates[i-1].forwscale;
  }

  val = 0;
  i = stateslen-1;
  for (j = 0; j < 16; j++) {
    val += exp(-1*(allstates[i].forwards[j])); // matei change: bug fix
  }
  return -log(val)+ allstates[i].forwscale;
}

double forward_backward (states* allstates, int stateslen) {
  double no1, no2;
  no1 = do_forwards(allstates, stateslen);
  no2 = do_backwards(allstates, stateslen);

#ifdef DEBUG_POST_SW
  fprintf (stderr, "SANITY CHECK: no1 == no2 %g %g\n", no1, no2);
#endif

  // don't really want a hard assert due to precision issues
  return no1;
}


/*********************************************************************************
 *
 * END Forward-backward code
 *
 *********************************************************************************/


int
post_sw_setup(int _max_len,
	      double _pr_snp, double _pr_xover, double _pr_del_open, double _pr_del_extend, double _pr_ins_open, double _pr_ins_extend,
	      bool _use_read_qvs, bool _use_sanger_qvs, int _qual_vector_offset, int _qual_delta,
	      bool reset_stats)
{
  assert(0 == BASE_0);
  assert((BASE_A ^ BASE_C) == BASE_1);
  assert((BASE_A ^ BASE_G) == BASE_2);
  assert((BASE_A ^ BASE_T) == BASE_3);
  assert((BASE_C ^ BASE_G) == BASE_3);
  assert((BASE_C ^ BASE_T) == BASE_2);
  assert((BASE_G ^ BASE_T) == BASE_1);

  pr_snp = _pr_snp;
  pr_xover = _pr_xover;
  pr_del_open = _pr_del_open;
  pr_del_extend = _pr_del_extend;
  pr_ins_open = _pr_ins_open;
  pr_ins_extend = _pr_ins_extend;

  qual_delta = _qual_delta;
  use_read_qvs = _use_read_qvs;
  use_sanger_qvs = _use_sanger_qvs;
  if (!use_read_qvs) {
    default_qual = qv_from_pr_err(pr_xover);
    //pr_xover = pr_err_from_qv(default_qual);
  } else {
    qual_vector_offset = _qual_vector_offset;
  }

  //neglogsixteenth = -log(1.0/16.0);
  //neglogfourth = -log(1.0/4.0);

  max_len = _max_len;
  columns = (struct column *)xmalloc(max_len * sizeof(columns[0]));
  for (int i = 0; i < max_len; i++) {
    columns[i].lets = (int *)xmalloc(1 * sizeof(columns[i].lets[0]));
    columns[i].cols = (int *)xmalloc(1 * sizeof(columns[i].cols[0]));
    columns[i].letserrrate = (double *)xmalloc(1 * sizeof(columns[i].letserrrate[0]));
    columns[i].colserrrate = (double *)xmalloc(1 * sizeof(columns[i].colserrrate[0]));
  }

  if (reset_stats) {
    cells = invocs = 0;
    tc.type = DEF_FAST_TIME_COUNTER;
    tc.counter = 0;
  }

  initialized = 1;

  check = 0;

  return 1;
}


int
post_sw_cleanup()
{
  for (int i = 0; i < max_len; i++) {
    free(columns[i].lets);
    free(columns[i].cols);
    free(columns[i].letserrrate);
    free(columns[i].colserrrate);
  }
  free(columns);
  return 1;
}


int
post_sw_stats(uint64_t * _invocs, uint64_t * _cells, double * _secs)
{
  if (_invocs != NULL)
    *_invocs = invocs;
  if (_cells != NULL)
    *_cells = cells;
  if (_secs != NULL)
    *_secs = time_counter_get_secs(&tc);

  return 1;
}


/*
 * Extract genome sequence, read, and qvs of interest.
 */
static void
load_local_vectors(uint32_t * read, int _init_bp, char * qual, struct sw_full_results * sfrp)
{
  int start_run, col;
  int min_qv;
  int i, j;

  start_run = 0;
  min_qv = 10000;
  for (j = 0; j < sfrp->read_start; j++) {
    col = EXTRACT(read, j);
    if (col == BASE_N) {
      start_run = BASE_N;
      min_qv = 0;
      j = sfrp->read_start;
      break;
    }
    start_run ^= col;
    if (use_read_qvs)
      min_qv = MIN(min_qv, (int)qual[qual_vector_offset+j]);
  }

  len = 0;
  for (i = 0; sfrp->dbalign[i] != 0; i++) {
    if (sfrp->qralign[i] != '-') { // ow, it's a deletion; nothing to do
      if (sfrp->dbalign[i] != '-') { // MATCH
	columns[len].nlets = 1;
	columns[len].lets[0] = fasta_get_initial_base(COLOUR_SPACE, &sfrp->dbalign[i]); // => BASE_A/C/G/T
	columns[len].letserrrate[0] = pr_snp;
      } else {
	columns[len].nlets = 0;
      }

      // MATCH or INSERTION
      columns[len].ncols = 1;
      col = EXTRACT(read, j);
      if ((len == 0 && start_run == BASE_N) || col == BASE_N) {
	//columns[len].ncols = 0; // no emission
	columns[len].cols[0] = BASE_0;
	columns[len].colserrrate[0] = .75;
      } else {
	columns[len].cols[0] = EXTRACT(read, j) ^ (len == 0? start_run : 0);
	if (use_read_qvs) {
	  columns[len].colserrrate[0] = pr_err_from_qv((len == 0? MIN(min_qv, (int)qual[qual_vector_offset + j]) : (int)qual[qual_vector_offset + j]) - qual_delta);
	  if (!use_sanger_qvs) {
	    columns[len].colserrrate[0] /= (1 + columns[len].colserrrate[0]);
	  }
	  if (columns[len].colserrrate[0] > .75) columns[len].colserrrate[0] = .75;
	} else {
	  columns[len].colserrrate[0] = pr_xover;
	} 
      }
      columns[len].base_call = char_to_base(sfrp->qralign[i]);
      assert(base_to_char(columns[len].base_call, LETTER_SPACE) == toupper(sfrp->qralign[i]));

      len++;
      j++;
    }
  }
  init_bp = _init_bp;

#ifdef DEBUG_POST_SW
  int _i;
  fprintf(stderr, "db:  ");
  for (_i = 0; _i < len; _i++) {
    fprintf(stderr, "    %c", columns[_i].nlets > 0 ? base_to_char(columns[_i].lets[0], LETTER_SPACE) : '-');
  }
  fprintf(stderr, "\n");
  fprintf(stderr, "qr: %c", base_to_char(init_bp, LETTER_SPACE));
  for (_i = 0; _i < len; _i++) {
    fprintf(stderr, "  %c  ", (columns[_i].ncols > 0 ? base_to_char(columns[_i].cols[0], COLOUR_SPACE) : '-'));
  }
  fprintf(stderr, "\n");
  fprintf(stderr, "qv:  ");
  for (_i = 0; _i < len; _i++) {
    fprintf(stderr, "%3d  ", qv_from_pr_err(columns[_i].colserrrate[0]));
  }
  fprintf(stderr, "\n");
#endif
}


static void
fix_base_calls(uint32_t * read, struct sw_full_results * sfrp,
	       states* allstates, int stateslen)
{
  int i = 0;
  int j = 0;
  int prev_base = init_bp;
  sfrp->matches = 0;
  sfrp->mismatches = 0;
  sfrp->crossovers = 0;
  while (sfrp->qralign[i] != 0) {
    if (sfrp->qralign[i] != '-') {
      // position i in qralign corresponds to pos j in allstates
      int crt_base = allstates[j].max_posterior;
      sfrp->qralign[i] = base_to_char(crt_base, LETTER_SPACE);
      if ((prev_base ^ crt_base) == allstates[j].cols[0]) {
	sfrp->qralign[i] = toupper(base_to_char(crt_base, LETTER_SPACE));
      } else {
	sfrp->qralign[i] = tolower(base_to_char(crt_base, LETTER_SPACE));
	++(sfrp->crossovers);
      }
      if (sfrp->dbalign[i] != '-') {
	if (toupper(sfrp->dbalign[i]) == toupper(sfrp->qralign[i])) {
	  ++(sfrp->matches);
	} else {
	  ++(sfrp->mismatches);
	}
      }
      prev_base = crt_base;
      ++j;
    }
    ++i;
  }
  assert(j == stateslen);
}


static void
get_base_qualities(struct sw_full_results * sfrp)
{
  int i, k;

  sfrp->qual = (char *)xmalloc((strlen(sfrp->qralign) + 1) * sizeof(sfrp->qual[0]));
  for (i = 0, k = 0; sfrp->qralign[i] != 0; i++) {
    if (sfrp->qralign[i] != '-') {
      int tmp = columns[k].base_call != BASE_N ? qv_from_pr_corr(columns[k].posterior[columns[k].base_call]) : 0;
      if (tmp > 40)
	tmp = 40;
      sfrp->qual[k] = 33 + tmp; // always 33+ in SAM
      k++;
    }
  }
  assert(k == len);
  sfrp->qual[k] = 0;
}


static double
get_posterior(struct sw_full_results * sfrp, double total_score)
{
  int i;
  double res;

  res = exp(-total_score); // - len * neglogfourth));
  for (i = 0; sfrp->dbalign[i] != 0; i++) {
    if (sfrp->dbalign[i] == '-') {
      res *= pr_ins_extend;
      if (i == 0 || sfrp->dbalign[i-1] != '-') {
	res *= pr_ins_open;
      }
    } else if (sfrp->qralign[i] == '-') {
      res *= pr_del_extend;
      if (i == 0 || sfrp->qralign[i-1] != '-') {
	res *= pr_del_open;
      }
    }
  }

  return res;
}


/*
 * Main method, called after full SW.
 */
void
post_sw(uint32_t * read, int _init_bp, char * qual,
	struct sw_full_results * sfrp)
{
  double total_score;

  //llint before = rdtsc(), after;
  TIME_COUNTER_START(tc);

  invocs++;

  assert(sfrp != NULL);
  assert(sfrp->dbalign != NULL);

  if (!initialized)
    abort();

#ifdef DEBUG_POST_SW
  int _i, _j, _last_base, _new_base;
  char const * spaces = "                        ";
  fprintf(stderr, "Post SW\n");
  fprintf(stderr, "dbalign: %s%s\n", spaces + strlen(spaces) - sfrp->read_start - 1, sfrp->dbalign);
  fprintf(stderr, "qralign: %s%s (offset: %d)\n", spaces + strlen(spaces) - sfrp->read_start - 1, sfrp->qralign, sfrp->read_start);
  fprintf(stderr, "read cs: %c", base_to_char(_init_bp, LETTER_SPACE));
  for (_i = 0, _j = 0; _i < (int)sfrp->read_start + (int)strlen(sfrp->qralign); _i++) {
    if (_j < sfrp->read_start) {
      fprintf(stderr, "%c", base_to_char(EXTRACT(read, _j), COLOUR_SPACE));
      _j++;
    } else {
      if (sfrp->qralign[_i - sfrp->read_start] == '-') {
	fprintf(stderr, "-");
      } else {
	fprintf(stderr, "%c", base_to_char(EXTRACT(read, _j), COLOUR_SPACE));
	_j++;
      }
    }
  }
  fprintf(stderr, "\n");
  fprintf(stderr, "read ls:  ");
  _last_base = _init_bp;
  for (_i = 0, _j = 0; _i < (int)sfrp->read_start + (int)strlen(sfrp->qralign); _i++) {
    if (_j < sfrp->read_start) {
      _new_base = cstols(_last_base, EXTRACT(read, _j), false);
      fprintf(stderr, "%c", base_to_char(_new_base, LETTER_SPACE));
      _last_base = _new_base;
      _j++;
    } else {
      if (sfrp->qralign[_i - sfrp->read_start] == '-') {
	fprintf(stderr, "-");
      } else {
	_new_base = cstols(_last_base, EXTRACT(read, _j), false);
	fprintf(stderr, "%c", base_to_char(_new_base, LETTER_SPACE));
	_last_base = _new_base;
	_j++;
      }
    }
  }
  fprintf(stderr, "\n");
  fprintf(stderr, "read qv:  ");
  for (_i = 0, _j = 0; _i < (int)sfrp->read_start + (int)strlen(sfrp->qralign); _i++) {
    if (_j < sfrp->read_start) {
      fprintf(stderr, "%c", use_read_qvs? qual[_j] : qual_delta + default_qual);
      _j++;
    } else {
      if (sfrp->qralign[_i - sfrp->read_start] == '-') {
	fprintf(stderr, " ");
      } else {
	fprintf(stderr, "%c", use_read_qvs? qual[_j] : qual_delta + default_qual);
	_j++;
      }
    }
  }
  fprintf(stderr, "\n");
#endif

  load_local_vectors(read, _init_bp, qual, sfrp);
  total_score = forward_backward(columns, len);
  post_traceback(columns, len, total_score);
  fix_base_calls(read, sfrp, columns, len);
  get_base_qualities(sfrp);
  sfrp->posterior = get_posterior(sfrp, total_score);

#ifdef DEBUG_POST_SW
  fprintf(stderr, "don: ");
  for (_i = 0; _i < len; _i++) {
    fprintf(stderr, "    %c", base_to_char(columns[_i].max_posterior, LETTER_SPACE));
  }
  fprintf(stderr, "\n");
  fprintf(stderr, "bqv: ");
  for (_i = 0; _i < len; _i++) {
    int res = columns[_i].posterior[columns[_i].max_posterior] > 1 - .00000001? 80 :
      (int)(-10.0*(log(1 - columns[_i].posterior[columns[_i].max_posterior])/log(10.0)));
    fprintf(stderr, "  %3d", res);
  }
  fprintf(stderr, "\n");

  fprintf(stderr, "qralign: ");
  for (_i = 0, _j = 0; sfrp->qralign[_i] != 0; _i++) {
    if (sfrp->qralign[_i] != '-') {
      fprintf(stderr, "  %c", sfrp->qralign[_i]);
    }
  }
  fprintf(stderr, "\n");
  fprintf(stderr, "bqv:     ");
  for (_i = 0, _j = 0; sfrp->qralign[_i] != 0; _i++) {
    if (sfrp->qralign[_i] != '-') {
      fprintf(stderr, "%3d", (int)(sfrp->qual[_j] - qual_delta));
      _j++;
    }
  }
  fprintf(stderr, "\n");

  printStates(columns, len, stderr);
#endif

  cells += 16*len;
  //after = rdtsc();
  //ticks += MAX(after - before, 0);
  TIME_COUNTER_STOP(tc);
}
