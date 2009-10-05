/*
 * mapper.h
 *
 *  Created on: 2009-10-02
 *      Author: dlister
 */

#ifndef MAPPER_H_
#define MAPPER_H_

#include "../common/bitmap.h"
#include <stdlib.h>

#define HASH_TABLE_POWER	12	/* 4^HASH_POWER entries in table */

extern struct seed_type *seed;
extern uint32_t **seed_hash_mask;
extern uint max_seed_span;
extern uint n_seeds;
extern u_int	nkmers;				/* total kmers of reads loaded*/


struct seed_type {
  bitmap_type	mask[1];	/* a bitmask, least significant bit = rightmost match */
  uint32_t	span;		/* max 64 (could be uint8_t) */
  uint32_t	weight;		/* max 64 */
  // aligned to 8B
};

extern size_t
power(size_t base, size_t exp);

extern uint32_t
kmer_to_mapidx_hash(uint32_t *kmerWindow, u_int sn);


#endif /* MAPPER_H_ */
