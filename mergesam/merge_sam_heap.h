#ifndef __MERGE_SAM_HEAP__
#define __MERGE_SAM_HEAP__

#include "sam2pretty_lib.h"

typedef struct heap_pa_elem heap_pa_elem;
typedef struct heap_pa heap_pa;
struct heap_pa_elem {
        long long score;
	int idist;
        pretty * rest;
};

struct heap_pa {
        struct heap_pa_elem *array;
        uint32_t capacity;
        uint32_t load;
};

void heap_pa_init(heap_pa * h, uint32_t capacity);
void heap_pa_insert(heap_pa * h, heap_pa_elem * e);
void heap_pa_get_max(heap_pa * h, heap_pa_elem * dest);
void heap_pa_extract_max(heap_pa * h, heap_pa_elem * dest);
void heap_pa_destroy(struct heap_pa * h);
#endif
