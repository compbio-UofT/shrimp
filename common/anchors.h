#ifndef _ANCHORS_H
#define _ANCHORS_H

#include <stdint.h>
#include <sys/types.h>
#include "../gmapper/gmapper-definitions.h"
#include "util.h"


void	anchor_join(struct anchor const *, int, struct anchor *);
void	anchor_widen(struct anchor *, int);
void	anchor_get_x_range(struct anchor const *, int, int, int, int *, int *);
void	anchor_uw_join(struct anchor *, struct anchor const *);
int	anchor_uw_cmp(void const *, void const *);


static inline bool
anchor_uw_colinear(struct anchor const * a1, struct anchor const * a2) {
  return a1->x - a1->y == a2->x - a2->y;
}

static inline bool
anchor_uw_intersect(struct anchor const * a1, struct anchor const * a2) {
  return anchor_uw_colinear(a1, a2)
    && ((a1->x == a2->x)
	|| (a1->x < a2->x && a2->x <= a1->x + a1->length)
	|| (a2->x < a1->x && a1->x <= a2->x + a2->length));
}

static inline void
anchor_reverse(struct anchor * a, int x_len, int y_len) {
  a->x = -a->x + (x_len - 1) - (a->length - 1) - (a->width - 1);
  a->y = -a->y + (y_len - 1) - (a->length - 1) + (a->width - 1);
}

static inline void
anchor_to_relative(struct anchor * a, llint goff) {
  a->x -= goff;
}

#endif
