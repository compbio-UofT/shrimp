#ifndef __READ_HIT_HEAP__
#define __READ_HIT_HEAP__

#include <stdint.h>
#include "../gmapper/gmapper.h"
#include "util.h"
#include "anchors.h"

typedef struct read_hit_heap_e read_hit_heap_e;
struct read_hit_heap_e {
        int32_t  score;
	int32_t  max_score;
	int32_t  isize_score;
	int32_t  isize;
        read_hit * rh[2];
};

typedef struct read_hit_heap read_hit_heap;
struct read_hit_heap {
        struct read_hit_heap_e *array;
        int32_t capacity;
        int32_t load;	
	bool alignment_compare;
	bool (*cmp)(read_hit_heap_e * a, read_hit_heap_e *b);
};

void read_hit_heap_init(read_hit_heap * h, uint32_t capacity, bool alignment_compare);
void read_hit_heap_insert(read_hit_heap * h, read_hit_heap_e * e);
void read_hit_heap_get_min(read_hit_heap * h, read_hit_heap_e * dest);
void read_hit_heap_extract_min(read_hit_heap * h, read_hit_heap_e * dest);
void read_hit_heap_destroy(struct read_hit_heap * h);
void read_hit_heap_insert_bounded_strata(struct read_hit_heap * h, struct read_hit_heap_e *e);
void read_hit_heap_insert_bounded(struct read_hit_heap * h, struct read_hit_heap_e *e);
void read_hit_heap_heapify(struct read_hit_heap * h);
void read_hit_heap_sort(struct read_hit_heap * h);
void read_hit_heap_remove_dups_and_sort(struct read_hit_heap *h);
void read_hit_heap_remove_dups_and_sort_strata(struct read_hit_heap * rhh);
#endif
