/*
 * mapper.c
 *
 *  Created on: 2009-10-02
 *      Author: dlister
 *
 *
 *  Used to refactor common items from rmapper and gmapper
 */

#include <stdio.h>

#include "../mapper/mapper.h"
#include "../common/util.h"


struct seed_type *seed = NULL;
uint32_t **seed_hash_mask = NULL;
uint max_seed_span = 0;
uint n_seeds = 0;
u_int	nkmers = 0;				/* total kmers of reads loaded*/

size_t
power(size_t base, size_t exp)
{
	size_t result = 1;

	while (exp > 0) {
		if ((exp % 2) == 1)
			result *= base;
		base *= base;
		exp /= 2;
	}

	return (result);
}


/* pulled off the web; this may or may not be any good */
static uint32_t
hash(uint32_t a)
{
	a = (a+0x7ed55d16) + (a<<12);
	a = (a^0xc761c23c) ^ (a>>19);
	a = (a+0x165667b1) + (a<<5);
	a = (a+0xd3a2646c) ^ (a<<9);
	a = (a+0xfd7046c5) + (a<<3);
	a = (a^0xb55a4f09) ^ (a>>16);
	return (a);
}

/* hash-based version or kmer -> map index function for larger seeds */
uint32_t
kmer_to_mapidx_hash(uint32_t *kmerWindow, u_int sn)
{
	static uint32_t maxidx = ((uint32_t)1 << 2*HASH_TABLE_POWER) - 1;

	uint32_t mapidx = 0;
	uint i;

	assert(seed_hash_mask != NULL);

	for (i = 0; i < BPTO32BW(max_seed_span); i++)
		mapidx = hash((kmerWindow[i] & seed_hash_mask[sn][i]) ^ mapidx);

	return mapidx & maxidx;
}
