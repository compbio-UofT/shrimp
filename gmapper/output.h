#ifndef _OUTPUT_H
#define _OUTPUT_H

#ifdef __cplusplus
//extern "C" {
#endif

#include "gmapper.h"
#include "../common/util.h"
#include "../common/fasta.h"

#undef EXTERN
#undef STATIC
#ifdef _MODULE_OUTPUT
#define EXTERN(_type, _id, _init_val) _type _id = _init_val
#define STATIC(_type, _id, _init_val) static _type _id = _init_val
#else
#define EXTERN(_type, _id, _init_val) extern _type _id
#define STATIC(_type, _id, _init_val)
#endif


void	hit_output(struct read_entry *, struct read_hit *, struct read_hit *,
		   char * *, char * *, bool, int *, int);
void	read_output(struct read_entry *, struct read_hit * *, int *);
void	readpair_output(struct read_entry *, struct read_entry *, struct read_hit_pair *, int *);


#ifdef __cplusplus
//} /* extern "C" */
#endif

#endif
