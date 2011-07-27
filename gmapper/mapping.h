#ifndef _MAPPING_H
#define _MAPPING_H

#ifdef __cplusplus
//extern "C" {
#endif

#include "gmapper.h"
#include "../common/util.h"
#include "../common/fasta.h"

#undef EXTERN
#undef STATIC
#ifdef _MODULE_MAPPING
#define EXTERN(_type, _id, _init_val) _type _id = _init_val
#define STATIC(_type, _id, _init_val) static _type _id = _init_val
#else
#define EXTERN(_type, _id, _init_val) extern _type _id
#define STATIC(_type, _id, _init_val)
#endif


void		handle_read(read_entry *, struct read_mapping_options_t *, int);
void		handle_readpair(pair_entry *, struct readpair_mapping_options_t *, int);
int		get_insert_size(read_hit *, read_hit *);


static inline double
get_pr_missed(read_entry * re_p)
{
  if (re_p->read_len < 40)
    return 1e-10;
  else if (re_p->read_len < 60)
    return 1e-14;
  else
    return 1e-16;
}

static inline double
pr_random_mapping_given_score(read_entry * re_p, int score)
{
  int read_len = re_p->read_len;
  if (score > read_len * match_score)
    return 1e-200;
  if (shrimp_mode == MODE_COLOUR_SPACE)
  {
    // how many crossovers would give this score?
    int n_xovers = ceil_div(read_len * match_score - score, abs(crossover_score));
    // pr of observing those by chance
    double tmp = -log_nchoosek(read_len, n_xovers) - n_xovers * log(3) + read_len * log(4);
    return exp(-tmp);
  }
  else // LETTER_SPACE
  {
    // how many errors/mismatches would give this score?
    int n_mismatches = ceil_div(read_len * match_score - score, abs(mismatch_score - match_score));
    // pr of observing those by chance
    double tmp = -log_nchoosek(read_len, n_mismatches) - n_mismatches * log(3) + read_len * log(4);
    return exp(-tmp);
  }
}


#ifdef __cplusplus
//} /* extern "C" */
#endif

#endif
