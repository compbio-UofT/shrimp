#ifndef _HEAP_H
#define _HEAP_H

#include <stdlib.h>
#include <string.h>
#include <assert.h>


#define DEF_STRUCT_HEAP_ELEM(_key_t,_rest_t,_id)	\
  struct heap_##_id##_elem {				\
    _key_t key;						\
    _rest_t rest;					\
  };

#define DEF_STRUCT_HEAP(_id)						\
  struct heap_##_id {							\
    struct heap_##_id##_elem *	array;					\
    uint		capacity;					\
    uint		load;						\
  };

#define DEF_HEAP_INIT(_id)						\
  static void								\
  heap_##_id##_init(struct heap_##_id * h, uint capacity)		\
  {									\
    assert(h != NULL);							\
									\
    h->array = (struct heap_##_id##_elem *)xmalloc(capacity * sizeof(struct heap_##_id##_elem)); \
    h->capacity = capacity;						\
    h->load = 0;							\
  }

#define DEF_HEAP_DESTROY(_id)			\
  static void					\
  heap_##_id##_destroy(struct heap_##_id * h)	\
  {						\
    assert(h != NULL);				\
						\
    free(h->array);				\
  }

#define DEF_HEAP_PERCOLATE_UP(_id)					\
  static void								\
  heap_##_id##_percolate_up(struct heap_##_id * h, uint node)	\
  {									\
    struct heap_##_id##_elem tmp;					\
    uint parent;							\
									\
    parent = node / 2;							\
    while (node > 1 && h->array[node-1].key < h->array[parent-1].key) {	\
      tmp = h->array[parent-1];						\
      h->array[parent-1] = h->array[node-1];				\
      h->array[node-1] = tmp;						\
									\
      node = parent;							\
      parent = node / 2;						\
    }									\
  }

#define DEF_HEAP_PERCOLATE_DOWN(_id)					\
  static void								\
  heap_##_id##_percolate_down(struct heap_##_id * h, uint node)	\
  {									\
    struct heap_##_id##_elem tmp;					\
    uint left, right, min;						\
									\
    do {								\
      left = node * 2;							\
      right = left + 1;							\
      min = node;							\
									\
      if (left <= h->load && h->array[left-1].key < h->array[node-1].key) \
	min = left;							\
									\
      if (right <= h->load && h->array[right-1].key < h->array[min-1].key) \
	min = right;							\
									\
      if (min == node)							\
	break;								\
									\
      tmp = h->array[min-1];						\
      h->array[min-1] = h->array[node-1];				\
      h->array[node-1] = tmp;						\
									\
      node = min;							\
    } while (1);							\
  }

#define DEF_HEAP_EXTRACT_MIN(_id)					\
  static void								\
  heap_##_id##_extract_min(struct heap_##_id * h, struct heap_##_id##_elem * dest) \
  {									\
    assert(h != NULL && h->load > 0);					\
									\
    *dest = h->array[0];						\
    h->load--;								\
    if (h->load > 0) {							\
      h->array[0] = h->array[h->load];					\
      									\
      heap_##_id##_percolate_down(h, 1);				\
    }									\
  }

#define DEF_HEAP_INSERT(_id)						\
  static void								\
  heap_##_id##_insert(struct heap_##_id * h, struct heap_##_id##_elem * e) \
  {									\
    assert(h != NULL && h->load < h->capacity);				\
									\
    h->array[h->load] = *e;						\
    h->load++;								\
									\
    heap_##_id##_percolate_up(h, h->load);				\
  }

#define DEF_HEAP(_key_t,_rest_t,_id)		\
  DEF_STRUCT_HEAP_ELEM(_key_t,_rest_t,_id)	\
  DEF_STRUCT_HEAP(_id)				\
  DEF_HEAP_INIT(_id)				\
  DEF_HEAP_DESTROY(_id)				\
  DEF_HEAP_PERCOLATE_UP(_id)			\
  DEF_HEAP_PERCOLATE_DOWN(_id)			\
  DEF_HEAP_EXTRACT_MIN(_id)			\
  DEF_HEAP_INSERT(_id)

DEF_HEAP(uint32_t,uint,uu)


#endif
