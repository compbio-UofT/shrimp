
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include "sam2pretty_lib.h"
#include "merge_sam_heap.h"


static bool e_compare(heap_pa_elem * a, heap_pa_elem * b) {
	if (a->score>b->score) {
		return true;
	} else if (a->score==b->score) {
		if (a->isize<b->isize) {
			return true;
		} else {
			return false;
		}
	}
	return false;
}


void heap_pa_init(struct heap_pa * h, uint32_t capacity) {
	assert(h != NULL);
	h->array = (struct heap_pa_elem *)malloc(capacity * sizeof(struct heap_pa_elem)); 
	if (h->array==NULL) { 
		fprintf(stderr,"Failed to allocate heap!\n"); 
		exit(1); 
	} 
	h->capacity = capacity;
	h->load = 0;
}

void heap_pa_destroy(struct heap_pa * h) {
	assert(h != NULL);

	free(h->array);
}

void heap_pa_percolate_up(struct heap_pa * h, uint32_t node) {
	struct heap_pa_elem tmp;
	uint32_t parent;

	parent = node / 2;
	while (node > 1 && e_compare(h->array+node-1,h->array+parent-1)) {
		tmp = h->array[parent-1];
		h->array[parent-1] = h->array[node-1];
		h->array[node-1] = tmp;

		node = parent;
		parent = node / 2;
	}
}

void heap_pa_percolate_down(struct heap_pa * h, uint32_t node) {
	struct heap_pa_elem tmp;
	uint32_t left, right, max;

	do {
		left = node * 2;
		right = left + 1;
		max = node;

		if (left <= h->load && e_compare(h->array+left-1,h->array+node-1)) 
		max = left;

		if (right <= h->load && e_compare(h->array+right-1,h->array+max-1)) 
		max = right;

		if (max == node)
		break;

		tmp = h->array[max-1];
		h->array[max-1] = h->array[node-1];
		h->array[node-1] = tmp;

		node = max;
	} while (1);
}

void heap_pa_extract_max(struct heap_pa * h, struct heap_pa_elem * dest)  {
	assert(h != NULL && h->load > 0);

	*dest = h->array[0];
	h->load--;
	if (h->load > 0) {
		h->array[0] = h->array[h->load];

		heap_pa_percolate_down(h, 1);
	}
}

void heap_pa_replace_max(struct heap_pa * h, struct heap_pa_elem * e) {
	assert(h != NULL && h->load > 0);

	h->array[0] = *e;
	heap_pa_percolate_down(h, 1);
}


void heap_pa_get_max(struct heap_pa * h, struct heap_pa_elem * dest) {
	assert(h != NULL && h->load > 0);

	*dest = h->array[0];
}

void heap_pa_insert(struct heap_pa * h, struct heap_pa_elem * e) {
	assert(h != NULL && h->load < h->capacity);

	h->array[h->load] = *e;
	h->load++;

	heap_pa_percolate_up(h, h->load);
}

void heap_pa_heapify(struct heap_pa * h){
	assert(h != NULL);

	uint32_t node;
	for (node = h->load / 2; node >= 1; node--) {
		heap_pa_percolate_down(h, node);
	}
}

void heap_pa_heapsort(struct heap_pa * h){ 
	assert(h != NULL);

	uint32_t orig_load = h->load;
	heap_pa_elem tmp;

	while (h->load > 1) {
		tmp = h->array[0];
		h->array[0] = h->array[h->load-1];
		h->array[h->load-1] = tmp;
		h->load--;

		heap_pa_percolate_down(h, 1);
	}
	h->load = orig_load;
}


