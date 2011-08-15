#ifndef _SW_POST_H
#define _SW_POST_H

#include <stdint.h>
#include "../common/sw-full-common.h"


int	post_sw_setup(int, double, double, double, double, double, double, bool, bool, int, int, bool);
int	post_sw_cleanup();
int	post_sw_stats(uint64_t *, uint64_t *, double *);

void	post_sw(uint32_t *, int, char *, struct sw_full_results *);


#endif
