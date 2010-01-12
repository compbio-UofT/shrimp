#ifndef _SW_GAPLESS_H
#define _SW_GAPLESS_H


#include <stdbool.h>
#include <stdlib.h>
#include <stdint.h>


int	sw_gapless_setup(int, int, bool);
void	sw_gapless_stats(uint64_t *, uint64_t *, uint64_t *);
int	sw_gapless(uint32_t *, int, uint32_t *, int, int, int,
		   uint32_t *, int, bool);


#endif
