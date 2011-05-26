#ifndef _MY_ALLOC_H
#define _MY_ALLOC_H

#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <assert.h>
#include <stdarg.h>
#include "../common/stats.h"


#define MYALLOC_COUNT_IT	0x1
#define MYALLOC_WARN_MAX	0x2
#define MYALLOC_ERR_MAX		0x4
#define MYALLOC_WARN_FAIL	0x8
#define MYALLOC_ERR_FAIL	0xF

extern size_t max_mem;
extern size_t crt_mem;

extern bool my_alloc_initialized;
extern bool warned_max;
extern bool warned_fail;


static inline void
my_alloc_init(size_t _max_mem)
{
  assert(sizeof(size_t) == sizeof(void *));
  max_mem = _max_mem;
  crt_mem = 0;
  warned_max = false;
  warned_fail = false;
  my_alloc_initialized = true;
}

static inline void *
my_malloc(size_t size, count_t * counter, int options, char const * msg, ...)
{
#ifndef NDEBUG
  va_list fmtargs;
#endif
  void * res;

#pragma omp critical (my_alloc)
  {
    if (options & MYALLOC_COUNT_IT) {
      if (crt_mem + size > max_mem) {
	if ((options & MYALLOC_WARN_MAX) && !warned_max) {
	  warned_max = true;
	  fprintf(stderr, "my-alloc warning: exceeding maximum memory"
#ifdef NDEBUG
		  "\n");
#else
	  ": ");
	  va_start(fmtargs, msg);
	  vfprintf(stderr, msg, fmtargs);
	  va_end(fmtargs);
#endif
	}
	if (options & MYALLOC_ERR_MAX) {
	  fprintf(stderr, "my-alloc error: exceeding maximum memory"
#ifdef NDEBUG
		  "\n");
#else
	  ": ");
	  va_start(fmtargs, msg);
	  vfprintf(stderr, msg, fmtargs);
	  va_end(fmtargs);
#endif
	  exit(1);
	}
      }
    }

    res = malloc(size);

    if (res == NULL) {
      if ((options & MYALLOC_WARN_FAIL) && !warned_fail) {
	warned_fail = true;
	fprintf(stderr, "my-alloc warning: malloc failed"
#ifdef NDEBUG
		"\n");
#else
	": ");
	va_start(fmtargs, msg);
	vfprintf(stderr, msg, fmtargs);
	va_end(fmtargs);
#endif
      }
      if (options & MYALLOC_ERR_FAIL) {
	fprintf(stderr, "my-alloc error: malloc failed"
#ifdef NDEBUG
		"\n");
#else
	": ");
	va_start(fmtargs, msg);
	vfprintf(stderr, msg, fmtargs);
	va_end(fmtargs);
#endif
	exit(1);
      }
    } else {
      if (options & MYALLOC_COUNT_IT)
	crt_mem += size;
      if (counter != NULL)
	count_add(counter, size);
    }

  }
  return res;
}


static inline void *
my_malloc(size_t size, count_t * counter, char const * msg, ...)
{
#ifndef NDEBUG
  va_list fmtargs;
#endif
  void * res;

  assert(my_alloc_initialized);

#pragma omp critical (my_alloc)
  {
    if (crt_mem + size > max_mem) {
      if (!warned_max) {
	warned_max = true;
	fprintf(stderr, "my-alloc warning: exceeding maximum memory"
#ifdef NDEBUG
		"\n");
#else
	": ");
	va_start(fmtargs, msg);
	vfprintf(stderr, msg, fmtargs);
	va_end(fmtargs);
#endif
      }
    }

    res = malloc(size);

    if (res == NULL) {
      fprintf(stderr, "my-alloc error: malloc failed"
#ifdef NDEBUG
		"\n");
#else
	": ");
      va_start(fmtargs, msg);
      vfprintf(stderr, msg, fmtargs);
      va_end(fmtargs);
#endif
      exit(1);
    } else {
      crt_mem += size;
      if (counter != NULL)
	count_add(counter, (int64_t)size);
    }
  }

  return res;
}


static inline void *
my_calloc(size_t size, count_t * counter, char const * msg, ...)
{
#ifndef NDEBUG
  va_list fmtargs;
#endif
  void * res;

  assert(my_alloc_initialized);

#pragma omp critical (my_alloc)
  {
    if (crt_mem + size > max_mem) {
      if (!warned_max) {
	warned_max = true;
	fprintf(stderr, "my-alloc warning: exceeding maximum memory"
#ifdef NDEBUG
		"\n");
#else
	": ");
	va_start(fmtargs, msg);
	vfprintf(stderr, msg, fmtargs);
	va_end(fmtargs);
#endif
      }
    }

  res = calloc(size, 1);

    if (res == NULL) {
      fprintf(stderr, "my-alloc error: malloc failed"
#ifdef NDEBUG
		"\n");
#else
	": ");
      va_start(fmtargs, msg);
      vfprintf(stderr, msg, fmtargs);
      va_end(fmtargs);
#endif
      exit(1);
    } else {
      crt_mem += size;
      if (counter != NULL)
	count_add(counter, (int64_t)size);
    }
  }

  return res;
}


static inline void *
my_realloc(void * p, size_t size, size_t old_size, count_t * counter, char const * msg, ...)
{
#ifndef NDEBUG
  va_list fmtargs;
#endif
  void * res;

  assert(my_alloc_initialized);

#pragma omp critical (my_alloc)
  {
    if ((crt_mem - old_size) + size > max_mem) {
      if (!warned_max) {
	warned_max = true;
	fprintf(stderr, "my-alloc warning: exceeding maximum memory"
#ifdef NDEBUG
		"\n");
#else
	": ");
	va_start(fmtargs, msg);
	vfprintf(stderr, msg, fmtargs);
	va_end(fmtargs);
#endif
      }
    }

  res = realloc(p, size);

    if (res == NULL) {
      fprintf(stderr, "my-alloc error: malloc failed"
#ifdef NDEBUG
		"\n");
#else
	": ");
      va_start(fmtargs, msg);
      vfprintf(stderr, msg, fmtargs);
      va_end(fmtargs);
#endif
      exit(1);
    } else {
      crt_mem -= old_size;
      crt_mem += size;
      if (counter != NULL)
	count_add(counter, (int64_t)size - (int64_t)old_size);
    }
  }

  return res;
}


static inline void
my_free(void * p, size_t size, count_t * counter)
{
#pragma omp critical (myalloc)
  {
    assert(size <= crt_mem);
    free(p);
    crt_mem -= size;
    if (counter != NULL)
      count_add(counter, -((int64_t)size));
  }
}





#endif
