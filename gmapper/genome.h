#ifndef _GENOME_H
#define _GENOME_H

#ifdef __cplusplus
//extern "C" {
#endif

#include "gmapper.h"
#include "../common/util.h"
#include "../common/fasta.h"

#undef EXTERN
#undef STATIC
#ifdef _MODULE_GENOME
#define EXTERN(_type, _id, _init_val) _type _id = _init_val
#define STATIC(_type, _id, _init_val) static _type _id = _init_val
#else
#define EXTERN(_type, _id, _init_val) extern _type _id
#define STATIC(_type, _id, _init_val)
#endif

  
bool		save_genome_map_seed(const char *, int);
bool		load_genome_map_seed(const char *);
bool		save_genome_map(const char *);
bool		load_genome_map(const char *);
void		print_genomemap_stats();
void		free_genome();
bool		load_genome(char **, int);
void		trim_genome();
bool		genome_load_map_save_mmap(char *, char const *);
bool		genome_load_mmap(char const *);


#ifdef __cplusplus
//} /* extern "C" */
#endif

#endif
