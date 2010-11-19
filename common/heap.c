#include <stdlib.h>
#include <string.h>

#include "heap.h"


struct heap *
heap_init(struct heap * h, uint capacity, uint elem_size, uint key_offset) {
  if (h == NULL) {
    h = (struct heap *)xmalloc(sizeof(struct heap));
  }

  h->array = (char *)xmalloc(capacity * elem_size);
  h->result = (char *)xmalloc(elem_size);
  h->capacity = capacity;
  h->elem_size = elem_size;
  h->key_offset = key_offset;
  h->load = 0;

  return h;
}


void
heap_destroy(struct heap * h) {
  free(h->array);
  free(h->result);
}


void
heap_percolate_up(struct heap * h, uint node) {
  uint parent;

  parent = node / 2;
  while (node > 1
	 //&& compare(&h->array[(parent - 1)*h->elem_size], &h->array[(node - 1)*h->elem_size]) > 0
	 *(uint32_t *)(&h[(node - 1)*h->elem_size + h->key_offset]) < *(uint32_t *)(&h[(parent - 1)*h->elem_size + h->key_offset])
	 ) {
    // swap heap elements
    memcpy(h->result, &h->array[(parent - 1)*h->elem_size], h->elem_size);
    memcpy(&h->array[(parent - 1)*h->elem_size], &h->array[(node - 1)*h->elem_size], h->elem_size);
    memcpy(&h->array[(node - 1)*h->elem_size], h->result, h->elem_size);

    node = parent;
    parent = node / 2;
  }
} 

void
heap_percolate_down(struct heap * h, uint node) {
  uint left, right, min;

  do {
    left = node * 2;
    right = left + 1;
    min = node;

    if (left <= h->load
	//&& compare(&h->array[(left - 1)*h->elem_size], &h->array[(node - 1)*h->elem_size]) < 0
	*(uint32_t *)(&h[(left - 1)*h->elem_size + h->key_offset]) < *(uint32_t *)(&h[(node - 1)*h->elem_size + h->key_offset])
	)
      min = left;

    if (right <= h->load
	//&& compare(&h->array[(right - 1)*h->elem_size], &h->array[(min - 1)*h->elem_size]) < 0
	*(uint32_t *)(&h[(right - 1)*h->elem_size + h->key_offset]) < *(uint32_t *)(&h[(min - 1)*h->elem_size + h->key_offset])
	)
      min = right;

    if (min == node)
      break;

    memcpy(h->result, &h->array[(min - 1)*h->elem_size], h->elem_size);
    memcpy(&h->array[(min - 1)*h->elem_size], &h->array[(node - 1)*h->elem_size], h->elem_size);
    memcpy(&h->array[(node - 1)*h->elem_size], h->result, h->elem_size);
    node = min;
  } while (1);
}


void *
heap_extract_min(struct heap * h) {
  assert(h->load > 0);

  if (h->load == 1) {
    h->load = 0;
    return &h->array[0];
  }

  // swap element at the top with last
  memcpy(h->result, &h->array[(h->load - 1)*h->elem_size], h->elem_size);
  memcpy(&h->array[(h->load - 1)*h->elem_size], &h->array[0], h->elem_size);
  memcpy(&h->array[0], h->result, h->elem_size);

  // remove it conceptually
  h->load--;

  // percolate new top elem down
  heap_percolate_down(h, 1);

  return &h->array[(h->load)*h->elem_size];
}


void
heap_insert(struct heap *h, void * elem) {
  assert(h->load < h->capacity);

  // add element at the bottom of the heap
  memcpy(&h->array[(h->load)*h->elem_size], elem, h->elem_size);
  h->load++;

  // percolate it up
  heap_percolate_up(h, h->load);
}
