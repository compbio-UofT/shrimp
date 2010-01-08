#ifndef _ANCHORS_H
#define _ANCHORS_H

#include <stdint.h>
#include <sys/types.h>

struct anchor {
  int16_t	x;
  int16_t	y;
  uint8_t	y_alt;
  uint8_t	length;
  uint8_t	width;
  uint8_t	more_than_once;
};

struct uw_anchor {
  uint32_t	x;
  uint8_t	y;
  uint8_t	length;
  uint16_t	weight;
  uint32_t	cn;
};


void join_anchors(struct anchor *, uint, struct anchor *);
void widen_anchor(struct anchor *, uint);
void get_x_range(struct anchor *, uint, uint, int, int *, int *);

void uw_anchors_join(struct uw_anchor *, struct uw_anchor const *);


static inline bool
uw_anchors_colinear(struct uw_anchor const * a1, struct uw_anchor const * a2) {
  return (int64_t)(a1->x) - (int64_t)(a1->y) == (int64_t)(a2->x) - (int64_t)(a2->y);
}


static inline bool
uw_anchors_intersect(struct uw_anchor const * a1, struct uw_anchor const * a2) {
  return uw_anchors_colinear(a1, a2)
    && ((a1->x == a2->x)
	|| (a1->x < a2->x && a2->x <= a1->x + (uint32_t)a1->length)
	|| (a2->x < a1->x && a1->x <= a2->x + (uint32_t)a2->length));
}



    


#endif
