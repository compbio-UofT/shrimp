#ifndef _HEAP_H
#define _HEAP_H

#include <stdlib.h>
#include <string.h>
#include <assert.h>


#define DEF_STRUCT_HEAP_ELEM(_key_t,_rest_t,_id)	\
  typedef struct heap_##_id##_elem {			\
    _key_t key;						\
    _rest_t rest;					\
  } heap_##_id##_elem;

#define DEF_STRUCT_HEAP(_id)						\
  typedef struct heap_##_id {						\
    struct heap_##_id##_elem *	array;					\
    uint		capacity;					\
    uint		load;						\
  } heap_##_id;

#define DEF_HEAP_INIT(_id)						\
  static inline void							\
  heap_##_id##_init(struct heap_##_id * h, uint capacity)		\
  {									\
    assert(h != NULL);							\
									\
    h->array = (struct heap_##_id##_elem *)xmalloc(capacity * sizeof(struct heap_##_id##_elem)); \
    h->capacity = capacity;						\
    h->load = 0;							\
  }

#define DEF_HEAP_DESTROY(_id)			\
  static inline void					\
  heap_##_id##_destroy(struct heap_##_id * h)	\
  {						\
    assert(h != NULL);				\
						\
    free(h->array);				\
    h->capacity = 0;				\
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
  static inline void								\
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

#define DEF_HEAP_REPLACE_MIN(_id)					\
  static inline void								\
  heap_##_id##_replace_min(struct heap_##_id * h, struct heap_##_id##_elem * e) \
  {									\
    assert(h != NULL && h->load > 0);					\
									\
    h->array[0] = *e;							\
    heap_##_id##_percolate_down(h, 1);					\
  }


#define DEF_HEAP_GET_MIN(_id)						\
  static inline void								\
  heap_##_id##_get_min(struct heap_##_id * h, struct heap_##_id##_elem * dest) \
  {									\
    assert(h != NULL && h->load > 0);					\
									\
    *dest = h->array[0];						\
  }

#define DEF_HEAP_INSERT(_id)						\
  static inline void								\
  heap_##_id##_insert(struct heap_##_id * h, struct heap_##_id##_elem * e) \
  {									\
    assert(h != NULL && h->load < h->capacity);				\
									\
    h->array[h->load] = *e;						\
    h->load++;								\
									\
    heap_##_id##_percolate_up(h, h->load);				\
  }

#define DEF_HEAP_HEAPIFY(_id)				\
  static void						\
  heap_##_id##_heapify(struct heap_##_id * h)		\
  {							\
    assert(h != NULL);					\
							\
    uint node;						\
    for (node = h->load / 2; node >= 1; node--) {	\
      heap_##_id##_percolate_down(h, node);		\
    }							\
  }

#define DEF_HEAP_HEAPSORT(_id)			\
  static void					\
  heap_##_id##_heapsort(struct heap_##_id * h)	\
  {						\
    assert(h != NULL);				\
						\
    uint orig_load = h->load;			\
    heap_##_id##_elem tmp;			\
						\
    while (h->load > 1) {			\
      tmp = h->array[0];			\
      h->array[0] = h->array[h->load-1];	\
      h->array[h->load-1] = tmp;		\
      h->load--;				\
						\
      heap_##_id##_percolate_down(h, 1);	\
    }						\
    h->load = orig_load;			\
  }

#define DEF_HEAP_ELEM_CMP(_id)						\
  static int								\
  heap_##_id##_elem_cmp(void const * p1, void const * p2)		\
  {									\
    int64_t res = (int64_t)((heap_##_id##_elem *)p1)->key -		\
      (int64_t)((heap_##_id##_elem *)p2)->key;				\
    return res > (int64_t)0 ? -1 :					\
      res < (int64_t)0? 1 : 0;						\
  }

#define DEF_HEAP_QSORT(_id)					\
  static inline void						\
  heap_##_id##_qsort(struct heap_##_id * h)			\
  {								\
    qsort(h->array, h->load, sizeof(struct heap_##_id##_elem),	\
	  heap_##_id##_elem_cmp);				\
  }

#define DEF_HEAP_ELEM_CMP_DOUBLE(_id)						\
  static int								\
  heap_##_id##_elem_cmp_double(void const * p1, void const * p2)		\
  {									\
    double res = (double)((heap_##_id##_elem *)p1)->key -		\
      (double)((heap_##_id##_elem *)p2)->key;				\
    return res > 0.0 ? -1 :					\
      res < 0.0 ? 1 : 0;						\
  }

#define DEF_HEAP_QSORT_DOUBLE(_id)					\
  static inline void						\
  heap_##_id##_qsort_double(struct heap_##_id * h)			\
  {								\
    qsort(h->array, h->load, sizeof(struct heap_##_id##_elem),	\
	  heap_##_id##_elem_cmp_double);				\
  }

#define DEF_HEAP(_key_t,_rest_t,_id)		\
  DEF_STRUCT_HEAP_ELEM(_key_t,_rest_t,_id)	\
  DEF_STRUCT_HEAP(_id)				\
  DEF_HEAP_INIT(_id)				\
  DEF_HEAP_DESTROY(_id)				\
  DEF_HEAP_PERCOLATE_UP(_id)			\
  DEF_HEAP_PERCOLATE_DOWN(_id)			\
  DEF_HEAP_EXTRACT_MIN(_id)			\
  DEF_HEAP_GET_MIN(_id)				\
  DEF_HEAP_REPLACE_MIN(_id)			\
  DEF_HEAP_INSERT(_id)				\
  DEF_HEAP_ELEM_CMP(_id)			\
  DEF_HEAP_ELEM_CMP_DOUBLE(_id)			\
  DEF_HEAP_QSORT(_id)				\
  DEF_HEAP_QSORT_DOUBLE(_id)
  //DEF_HEAP_HEAPIFY(_id)
  //DEF_HEAP_HEAPSORT(_id)

//DEF_HEAP(uint32_t,uint,uu)


#define DEF_EXTHEAP_PERCOLATE_UP(_data_t,_id)				\
  static inline void							\
  extheap_##_id##_percolate_up(_data_t * a, int * load, int node)	\
  {									\
    _data_t tmp;							\
    int parent;								\
									\
    parent = node / 2;							\
    while (node > 1 && EXTHEAP_##_id##_CMP(a[node-1], a[parent-1])) {	\
      tmp = a[parent-1];						\
      a[parent-1] = a[node-1];						\
      a[node-1] = tmp;							\
									\
      node = parent;							\
      parent = node / 2;						\
    }									\
  }


#define DEF_EXTHEAP_PERCOLATE_DOWN(_data_t,_id)				\
  static inline void							\
  extheap_##_id##_percolate_down(_data_t * a, int * load, int node)	\
  {									\
    _data_t tmp;							\
    int left, right, min;						\
									\
    do {								\
      left = node * 2;							\
      right = left + 1;							\
      min = node;							\
									\
      if (left <= *load && EXTHEAP_##_id##_CMP(a[left-1], a[node-1]))	\
	min = left;							\
									\
      if (right <= *load && EXTHEAP_##_id##_CMP(a[right-1], a[min-1]))	\
	min = right;							\
									\
      if (min == node)							\
	break;								\
									\
      tmp = a[min-1];							\
      a[min-1] = a[node-1];						\
      a[node-1] = tmp;							\
									\
      node = min;							\
    } while (1);							\
  }

#define DEF_EXTHEAP_DELETE_MIN(_data_t,_id)				\
  static inline void							\
  extheap_##_id##_extract_min(_data_t * a, int * load)			\
  {									\
    assert(a != NULL && load != NULL && *load > 0);			\
									\
    (*load)--;								\
    if (*load > 0) {							\
      a[0] = a[*load];							\
      extheap_##_id##_percolate_down(a, load, 1);			\
    }									\
  }

#define DEF_EXTHEAP_REPLACE_MIN(_data_t,_id)				\
  static inline void							\
  extheap_##_id##_replace_min(_data_t * a, int * load, _data_t e)	\
  {									\
    assert(a != NULL && load != NULL && *load > 0);			\
									\
    a[0] = e;								\
    extheap_##_id##_percolate_down(a, load, 1);				\
  }


#define DEF_EXTHEAP_INSERT(_data_t,_id)					\
  static inline void							\
  extheap_##_id##_insert(_data_t * a, int * load, _data_t e)		\
  {									\
    assert(a != NULL);							\
									\
    a[*load] = e;							\
    (*load)++;								\
    extheap_##_id##_percolate_up(a, load, *load);			\
  }

#define DEF_EXTHEAP_HEAPIFY(_data_t,_id)		\
  static inline void					\
  extheap_##_id##_heapify(_data_t * a, int * load)	\
  {							\
    assert(a != NULL);					\
							\
    int node;						\
    for (node = *load / 2; node >= 1; node--) {		\
      extheap_##_id##_percolate_down(a, load, node);	\
    }							\
  }

#define DEF_EXTHEAP(_data_t,_id)		\
  DEF_EXTHEAP_PERCOLATE_UP(_data_t,_id)		\
  DEF_EXTHEAP_PERCOLATE_DOWN(_data_t,_id)	\
  DEF_EXTHEAP_DELETE_MIN(_data_t,_id)		\
  DEF_EXTHEAP_REPLACE_MIN(_data_t,_id)		\
  DEF_EXTHEAP_INSERT(_data_t,_id)		\
  DEF_EXTHEAP_HEAPIFY(_data_t,_id)



#endif
