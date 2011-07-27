#ifndef __MERGE_SAM_HEAP__
#define __MERGE_SAM_HEAP__

#include "sam2pretty_lib.h"

typedef struct heap_pa_elem heap_pa_elem;
struct heap_pa_elem {
        uint32_t  score;
	uint32_t  isize_score;
        pretty * rest;
};

typedef struct heap_pa heap_pa;
struct heap_pa {
        struct heap_pa_elem *array;
        int32_t capacity;
        int32_t load;
};

void heap_pa_init(heap_pa * h, uint32_t capacity);
void heap_pa_insert(heap_pa * h, heap_pa_elem * e);
void heap_pa_get_min(heap_pa * h, heap_pa_elem * dest);
void heap_pa_extract_min(heap_pa * h, heap_pa_elem * dest);
void heap_pa_destroy(struct heap_pa * h);
void heap_pa_insert_bounded_strata(struct heap_pa * h, struct heap_pa_elem *e);
void heap_pa_insert_bounded(struct heap_pa * h, struct heap_pa_elem *e);
#endif
