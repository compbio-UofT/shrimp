#ifndef _SEEDS_H
#define _SEEDS_H

#ifdef __cplusplus
//extern "C" {
#endif

#include "gmapper.h"
#include "gmapper-defaults.h"

#undef EXTERN
#undef STATIC
#ifdef _MODULE_SEEDS
#define EXTERN(_type, _id, _init_val) _type _id = _init_val
#define STATIC(_type, _id, _init_val) static _type _id = _init_val
#else
#define EXTERN(_type, _id, _init_val) extern _type _id
#define STATIC(_type, _id, _init_val)
#endif


STATIC(int,			default_seeds_cs_min_weight,	DEF_DEF_SEEDS_CS_MIN_WEIGHT);
STATIC(int,			default_seeds_cs_max_weight,	DEF_DEF_SEEDS_CS_MAX_WEIGHT);
STATIC(int,			default_seeds_cs_weight,	DEF_DEF_SEEDS_CS_WEIGHT);
STATIC(int,			default_seeds_cs_cnt[9],	DEF_DEF_SEEDS_CS_CNT);
STATIC(char const *,		default_seeds_cs[9][5],		DEF_DEF_SEEDS_CS);
STATIC(int,			default_seeds_ls_min_weight,	DEF_DEF_SEEDS_LS_MIN_WEIGHT);
STATIC(int,			default_seeds_ls_max_weight,	DEF_DEF_SEEDS_LS_MAX_WEIGHT);
STATIC(int,			default_seeds_ls_weight,	DEF_DEF_SEEDS_LS_WEIGHT);
STATIC(int,			default_seeds_ls_cnt[9],	DEF_DEF_SEEDS_LS_CNT);
STATIC(char const *,		default_seeds_ls[9][5],		DEF_DEF_SEEDS_LS);
STATIC(int,			default_seeds_mirna_cnt,	DEF_DEF_SEEDS_MIRNA_CNT);
STATIC(char const *,		default_seeds_mirna[5],		DEF_DEF_SEEDS_MIRNA);


bool		add_spaced_seed(char const *);
void		load_default_mirna_seeds();
bool		load_default_seeds(int);
void		init_seed_hash_mask();
char *		seed_to_string(int);
bool		valid_spaced_seeds();


#ifdef __cplusplus
//} /* extern "C" */
#endif

#endif
