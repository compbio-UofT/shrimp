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


static int	initialized;

static double	pr_snp;
static double	pr_error;
static bool	use_read_qvs;
static int	default_qual;		// if no qvs, use this instead

static int *	db;
static int *	qr;
static int *	qv;
static int	init_bp;

static int	qual_vector_offset;	// i.e. is there a useless qv for the initial base in cs?
static int	qual_delta;		// how much to subtract from chars to get the int

#pragma omp threadprivate(initialized,pr_snp,pr_error,use_read_qvs,default_qual,db,qr,qv,qual_vector_offset,qual_delta)


void
post_sw_setup(double _pr_snp, bool _use_read_qvs, int max_db_len, int max_qr_len,
	      int * _qual_vector_offset, int * _qual_delta)
{
  assert(0 == BASE_0);
  assert(BASE_A ^ BASE_C == BASE_1);
  assert(BASE_A ^ BASE_G == BASE_2);
  assert(BASE_A ^ BASE_T == BASE_3);
  assert(BASE_C ^ BASE_G == BASE_3);
  assert(BASE_C ^ BASE_T == BASE_2);
  assert(BASE_G ^ BASE_T == BASE_1);

  pr_snp = _pr_snp;
  use_read_qvs = _use_read_qvs;
  if (!use_read_qvs) {
    pr_error = 0.01;
    default_qual = 20;
  } else {
    assert(_qual_vector_offset != NULL && _qual_delta != NULL);

    qual_vector_offset = *_qual_vector_offset;
    qual_delta = *_qual_delta;
  }

  db = (int *)xmalloc(max_db_len * sizeof(db[0]));
  qr = (int *)xmalloc(max_qr_len * sizeof(qr[0]));
  qv = (int *)xmalloc(max_qr_len * sizeof(qv[0])); // qvs of qr

  initialized = 1;
}


/*
 * Extract genome sequence, read, and qvs of interest.
 */
void
load_local_vectors(uint32_t * read, int _init_bp, char * qual,
		   struct sw_full_results * sfrp)
{
  int prev_run, min_qv;
  int i, j, k;

  prev_run = 0;
  min_qv = 10000;
  for (j = 0; j < sfrp->read_start; j++) {
    prev_run ^= EXTRACT(read, j);
    if (use_read_qvs)
      min_qv = MIN(min_qv, (int)qual[qual_vector_offset+j] - qual_delta);
  }

  for (i = 0; sfrp->dbalign[i] != 0; i++) {
    if (sfrp->qralign[i] != '-') { // ow, it's a deletion; nothing to do
      if (sfrp->dbalign[i] != '-') { // match
	db[k] = (int8_t)fasta_get_initial_base(COLOUR_SPACE, &sfrp->dbalign[i]); // => BASE_A/C/G/T
	qr[k] = EXTRACT(read, j) ^ prev_run;
	if (use_read_qvs) {
	  qv[k] = MIN(min_qv, (int)qual[qual_vector_offset+j] - qual_delta);
	  min_qv = 10000;
	} else {
	  qv[k] = default_qual;
	}
	prev_run = 0;
	j++;
	k++;
      } else { // insertion: accumulate bases
	prev_run ^= EXTRACT(read, j);
	if (use_read_qvs)
	  min_qv = MIN(min_qv, (int)qual[qual_vector_offset+j] - qual_delta);
	j++;
      }
    }
  }
  init_bp = _init_bp;
}  


/*
 * Main method, called after full SW.
 */
void
post_sw(uint32_t * read, int _init_bp, char * qual,
	struct sw_full_results * sfrp)
{
  assert(sfrp != NULL);
  assert(sfrp->dbalign != NULL);

  if (!initialized)
    abort();

  load_local_vectors(read, init_bp, qual, sfrp);

}




static 
