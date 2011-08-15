/*	$Id: sw-full-ls.h,v 1.8 2009/06/12 21:27:35 rumble Exp $	*/

#ifndef _SW_FULL_LS_H
#define _SW_FULL_LS_H


#include "anchors.h"

int	sw_full_ls_setup(int, int, int, int, int, int, int, int, bool, int);
int	sw_full_ls_cleanup(void);
void	sw_full_ls_stats(uint64_t *, uint64_t *, double *);
void	sw_full_ls(uint32_t *, int, int, uint32_t *, int, int, int,
		   struct sw_full_results *, bool, struct anchor *, int, int);


#endif
