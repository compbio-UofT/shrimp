//TODO refactor all of this file
#include "../common/bitmap.h"

#define HASH_TABLE_POWER	12	/* 4^HASH_POWER entries in table */

struct seed_type {
  bitmap_type	mask[1];	/* a bitmask, least significant bit = rightmost match */
  uint32_t	span;		/* max 64 (could be uint8_t) */
  uint32_t	weight;		/* max 64 */
  // aligned to 8B
};
