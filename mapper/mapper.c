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
#include <string.h>

#include "../mapper/mapper.h"
#include "../common/util.h"


struct seed_type *seed = NULL;
uint32_t **seed_hash_mask = NULL;
uint max_seed_span = 0;
uint min_seed_span = MAX_SEED_SPAN;
uint32_t n_seeds = 0;
u_int	nkmers = 0;				/* total kmers of reads loaded*/
uint avg_seed_span = 0;

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

/*
 * Compress the given kmer into an index in 'readmap' according to the seed.
 * While not optimal, this is only about 20% of the spaced seed scan time.
 *
 * This is the original version for smaller seeds.
 *
 * XXX- This algorithm only considers bases 0-3, which implies overlap
 *      when we have other bases (mainly uracil, but also wobble codes).
 *      This won't affect sensitivity, but may cause extra S-W calls.
 */
uint32_t
kmer_to_mapidx_orig(uint32_t *kmerWindow, u_int sn)
{
	bitmap_type a = seed[sn].mask[0];
	uint32_t mapidx = 0;
	int i = 0;

	do {
		if ((a & 0x1) == 0x1) {
			mapidx <<= 2;
			mapidx |= ((kmerWindow[i/8] >> (i%8)*4) & 0x3);
		}
		a >>= 1;
		i++;

	} while (a != 0x0);

	assert(mapidx < power(4, seed[sn].weight));

	return mapidx;
}

bool
add_spaced_seed(const char *seedStr)
{
	uint i;

	seed = (struct seed_type *)xrealloc(seed, sizeof(struct seed_type) * (n_seeds + 1));
	seed[n_seeds].mask[0] = 0x0;
	seed[n_seeds].span = strlen(seedStr);
	seed[n_seeds].weight = strchrcnt(seedStr, '1');

	if (seed[n_seeds].span < 1
			|| seed[n_seeds].span > MAX_SEED_SPAN
			|| seed[n_seeds].weight < 1
			|| strchrcnt(seedStr, '0') != seed[n_seeds].span - seed[n_seeds].weight)
		return false;

	for (i = 0; i < seed[n_seeds].span; i++)
		bitmap_prepend(seed[n_seeds].mask, 1, (seedStr[i] == '1' ? 1 : 0));

	if (seed[n_seeds].span > max_seed_span)
		max_seed_span = seed[n_seeds].span;

	if (seed[n_seeds].span < min_seed_span)
		min_seed_span = seed[n_seeds].span;

	n_seeds++;

	avg_seed_span = 0;
	for(i =0; i < n_seeds;i++){
		avg_seed_span += seed[i].span;
	}
	avg_seed_span = avg_seed_span/n_seeds;

	return true;
}


void
load_default_seeds() {
	int i;

	n_seeds = 0;
	switch(shrimp_mode) {
	case MODE_COLOUR_SPACE:
		for (i = 0; i < default_spaced_seeds_cs_cnt; i++)
			add_spaced_seed(default_spaced_seeds_cs[i]);
		break;
	case MODE_LETTER_SPACE:
		for (i = 0; i < default_spaced_seeds_ls_cnt; i++)
			add_spaced_seed(default_spaced_seeds_ls[i]);
		break;
	case MODE_HELICOS_SPACE:
		for (i = 0; i < default_spaced_seeds_hs_cnt; i++)
			add_spaced_seed(default_spaced_seeds_hs[i]);
		break;
	}
}

void
init_seed_hash_mask(void)
{
	uint sn;
	int i;

	seed_hash_mask = (uint32_t **)xmalloc(sizeof(seed_hash_mask[0])*n_seeds);
	for (sn = 0; sn < n_seeds; sn++) {
		seed_hash_mask[sn] = (uint32_t *)xcalloc(sizeof(seed_hash_mask[sn][0])*BPTO32BW(max_seed_span));

		for (i = seed[sn].span - 1; i >= 0; i--)
			bitfield_prepend(seed_hash_mask[sn], max_seed_span,
					bitmap_extract(seed[sn].mask, 1, i) == 1? 0xf : 0x0);
	}
}

char *
seed_to_string(uint sn)
{
	static char buffer[100];
	bitmap_type tmp;
	int i;

	assert(sn < n_seeds);

	buffer[seed[sn].span] = 0;
	for (i = seed[sn].span - 1, tmp = seed[sn].mask[0];
			i >= 0;
			i--, tmp >>= 1) {
		if (bitmap_extract(&tmp, 1, 0) == 1)
			buffer[i] = '1';
		else
			buffer[i] = '0';
	}

	return buffer;
}

