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


#ifdef __cplusplus
//} /* extern "C" */
#endif

#endif
