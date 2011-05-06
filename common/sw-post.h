#ifndef _SW_POST_H
#define _SW_POST_H

#include <stdint.h>
#include "../common/sw-full-common.h"


void	post_sw_setup(double, bool, int, int, int);
void	post_sw(uint32_t *, int, char *, struct sw_full_results *);


#endif
