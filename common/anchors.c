#include <stdlib.h>
#include <stdio.h>
#include <limits.h>
#include <assert.h>
#include <stdbool.h>
#include "anchors.h"


void anchor_join(struct anchor const * anchors, int anchors_cnt,
		 struct anchor * dest) {
  llint border_nw_min, border_sw_min, border_ne_max, border_se_max;
  int i;

  assert(anchors != NULL && dest != NULL);

  border_nw_min = INT_MAX;
  border_sw_min = INT_MAX;
  border_ne_max = INT_MIN;
  border_se_max = INT_MIN;

  dest->weight = 0;
  dest->cn = anchors[0].cn;

  for (i = 0; i < anchors_cnt; i++) {
    llint border_nw, border_sw, border_ne, border_se;

    border_nw = anchors[i].x + anchors[i].y;
    border_sw = anchors[i].x - anchors[i].y;
    border_ne = border_sw + 2*(anchors[i].width - 1);
    border_se = border_nw + 2*(anchors[i].length - 1);

    border_nw_min = MIN(border_nw_min, border_nw);
    border_sw_min = MIN(border_sw_min, border_sw);
    border_ne_max = MAX(border_ne_max, border_ne);
    border_se_max = MAX(border_se_max, border_se);

    dest->weight += anchors[i].weight;

#ifdef DEBUG_ANCHORS
    fprintf(stderr, "i:%d border_nw_min=%d border_sw_min=%d border_ne_max=%d border_se_max=%d\n",
	    i, border_nw_min, border_sw_min, border_ne_max, border_se_max);
#endif
  }

  if ((border_nw_min + border_sw_min) % 2 != 0) border_nw_min--;
  dest->x = (border_nw_min + border_sw_min)/2;
  dest->y = border_nw_min - dest->x;

  if ((border_ne_max - border_sw_min) % 2 != 0) border_ne_max++;
  dest->width = (border_ne_max - border_sw_min)/2 + 1;

  if ((border_se_max - border_nw_min) % 2 != 0) border_se_max++;
  dest->length = (border_se_max - border_nw_min)/2 + 1;
}


void anchor_widen(struct anchor * anchor, int width) {
  assert(anchor != NULL);

  anchor->x -= width/2;
  anchor->y += width/2;
  anchor->width += width;
}


void anchor_get_x_range(struct anchor const * anchor, int x_len, int y_len, int y,
			int * x_min, int * x_max) {
  assert(anchor != NULL && x_min != NULL && x_max != NULL);

  if (y < anchor->y) {
    *x_min = 0;
  } else if (y <= anchor->y + (anchor->length - 1)) {
    *x_min = (int)(anchor->x + (y - anchor->y));
  } else {
    *x_min = (int)(anchor->x + anchor->length);
  }

  if (*x_min < 0)
    *x_min = 0;
  if (*x_min >= x_len)
    *x_min = x_len - 1;

  if (y < anchor->y - (anchor->width - 1)) {
    *x_max = (int)(anchor->x + (anchor->width - 1) - 1);
  } else if (y <= anchor->y - (anchor->width - 1) + (anchor->length - 1)) {
    *x_max = (int)(anchor->x + (anchor->width - 1) + (y - (anchor->y - (anchor->width - 1))));
  } else {
    *x_max = x_len - 1;
  }
  
  if (*x_max < 0)
    *x_max = 0;
  if (*x_max >= x_len)
    *x_max = x_len - 1;
}


void
anchor_uw_join(struct anchor * dest, struct anchor const * src)
{
  assert(src->width == 1 && dest->width == 1);
  assert(anchor_uw_colinear(dest, src));

  if (src->x < dest->x) {
    llint tmp = dest->x;
    dest->x = src->x;
    dest->y = src->y;
    if (src->x + src->length > tmp + dest->length) {
      dest->length = src->length;
    } else {
      dest->length += (int)(tmp - dest->x);
    }
  } else {
    if (src->x + src->length > dest->x + dest->length) {
      dest->length = (int)(src->x - dest->x + src->length);
    }
  }
  dest->weight += src->weight;
}


int
anchor_uw_cmp(void const * p1, void const * p2)
{
  if (((struct anchor *)p1)->x < ((struct anchor *)p2)->x)
    return -1;
  else if (((struct anchor *)p1)->x > ((struct anchor *)p2)->x)
    return 1;
  else
    return 0;
}
