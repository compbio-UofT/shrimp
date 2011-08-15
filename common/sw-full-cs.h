/*	$Id: sw-full-cs.h,v 1.6 2009/06/16 23:26:21 rumble Exp $	*/

#ifndef _SW_FULL_CS_H
#define _SW_FULL_CS_H

#include "anchors.h"
int	sw_full_cs_cleanup(void);
int	sw_full_cs_setup(int, int, int, int, int, int, int, int, int, bool, int, int = 0);
void	sw_full_cs_stats(uint64_t *, uint64_t *, double *);
void	sw_full_cs(uint32_t *, int, int, uint32_t *, int, int, int,
		   struct sw_full_results *, bool, bool, struct anchor *, int, int, int * = NULL);

#endif
