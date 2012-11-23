#ifndef _MY_ALLOC_H
#define _MY_ALLOC_H

#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <assert.h>
#include <stdarg.h>
#include <string.h>
#include "../common/stats.h"

#ifdef NDEBUG
#define MYALLOC_DISABLE_CRT
#define MYALLOC_DISABLE_ALERT
#endif

#ifndef MYALLOC_DISABLE_CRT
#define MYALLOC_ENABLE_CRT
#endif

#ifndef MYALLOC_DISABLE_ALERT
#define MYALLOC_ENABLE_ALERT
#endif


#define MYALLOC_COUNT_IT	0x1
#define MYALLOC_WARN_MAX	0x2
#define MYALLOC_ERR_MAX		0x4
#define MYALLOC_WARN_FAIL	0x8
#define MYALLOC_ERR_FAIL	0xF


extern bool my_alloc_initialized;
#ifdef MYALLOC_ENABLE_CRT
extern size_t max_mem;
extern size_t crt_mem;
extern bool warned_max;
#endif
#ifdef MYALLOC_ENABLE_ALERT
extern size_t alert_mem;
#endif
extern bool warned_fail;


static inline void
my_alloc_init(size_t _max_mem, size_t _alert_mem)
{
  assert(sizeof(size_t) == sizeof(void *));
#ifdef MYALLOC_ENABLE_CRT
  max_mem = _max_mem;
  crt_mem = 0;
  warned_max = false;
#endif
#ifdef MYALLOC_ENABLE_ALERT
  alert_mem = _alert_mem;
#endif
  warned_fail = false;
  my_alloc_initialized = true;
}


static inline void *
my_malloc(size_t size, count_t * counter, char const * msg, ...)
{
#ifndef NDEBUG
  va_list fmtargs;
#endif
  void * res;

  assert(my_alloc_initialized);
  //assert(size > 0);

#ifdef DEBUG_MY_ALLOC
  {
    fprintf(stderr, "my_malloc: +%lld: ", (long long)size);
    char new_msg[strlen(msg) + 1 + 200];
    strcpy(new_msg, msg);
    strcat(new_msg, "\n");
    va_start(fmtargs, msg);
    vfprintf(stderr, new_msg, fmtargs);
    va_end(fmtargs);
  }
#endif

#ifdef MYALLOC_ENABLE_CRT
#pragma omp critical (cs_my_alloc)
  {
    if (crt_mem + size > max_mem) {
      if (!warned_max) {
	warned_max = true;
#ifdef NDEBUG
	fprintf(stderr, "my_malloc warning: exceeding maximum memory\n");
#else
	char new_msg[strlen(msg) + 1 + 200];
	strcpy(new_msg, "my_malloc warning: exceeding maximum memory: ");
	strcat(new_msg, msg);
	strcat(new_msg, "\n");
	va_start(fmtargs, msg);
	vfprintf(stderr, new_msg, fmtargs);
	va_end(fmtargs);
#endif
      }
    }
#endif

#ifdef MYALLOC_ENABLE_ALERT
    if (size > alert_mem) {
#ifdef NDEBUG
      fprintf(stderr, "my_malloc alert: size=%lld\n", (long long)size);
#else
      char new_msg[strlen(msg) + 1 + 200];
      sprintf(new_msg, "my_malloc alert: size=%lld: %s\n", (long long)size, msg);
      va_start(fmtargs, msg);
      vfprintf(stderr, new_msg, fmtargs);
      va_end(fmtargs);
#endif
    }
#endif

#ifdef MYALLOC_ENABLE_CS
#pragma omp critical(alloc_cs)
    {
#endif
      res = malloc(size);
#ifdef MYALLOC_ENABLE_CS
    }
#endif

    if (size > 0 && res == NULL) {
      // spit error message and crash
#ifdef NDEBUG
	fprintf(stderr, "my_malloc error: malloc failed\n");
#else
	char new_msg[strlen(msg) + 1 + 200];
	strcpy(new_msg, "my_malloc warning: malloc failed: ");
	strcat(new_msg, msg);
	strcat(new_msg, "\n");
	va_start(fmtargs, msg);
	vfprintf(stderr, new_msg, fmtargs);
	va_end(fmtargs);
#endif
      exit(1);
    }
#ifdef MYALLOC_ENABLE_CRT
    else {
      crt_mem += size;
      if (counter != NULL)
	count_add(counter, (int64_t)size);
    }
  } // end of critical section
#endif

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

#ifdef DEBUG_MY_ALLOC
  {
    fprintf(stderr, "my_calloc: +%lld: ", (long long)size);
    char new_msg[strlen(msg) + 1 + 200];
    strcpy(new_msg, msg);
    strcat(new_msg, "\n");
    va_start(fmtargs, msg);
    vfprintf(stderr, new_msg, fmtargs);
    va_end(fmtargs);
  }
#endif

#ifdef MYALLOC_ENABLE_CRT
#pragma omp critical (cs_my_alloc)
  {
    if (crt_mem + size > max_mem) {
      if (!warned_max) {
	warned_max = true;
#ifdef NDEBUG
	fprintf(stderr, "my_calloc warning: exceeding maximum memory\n");
#else
	char new_msg[strlen(msg) + 1 + 200];
	strcpy(new_msg, "my_calloc warning: exceeding maximum memory: ");
	strcat(new_msg, msg);
	strcat(new_msg, "\n");
	va_start(fmtargs, msg);
	vfprintf(stderr, new_msg, fmtargs);
	va_end(fmtargs);
#endif
      }
    }
#endif

#ifdef MYALLOC_ENABLE_ALERT
    if (size > alert_mem) {
#ifdef NDEBUG
      fprintf(stderr, "my_calloc alert: size=%lld\n", (long long)size);
#else
      char new_msg[strlen(msg) + 1 + 200];
      sprintf(new_msg, "my_calloc alert: size=%lld: %s\n", (long long)size, msg);
      va_start(fmtargs, msg);
      vfprintf(stderr, new_msg, fmtargs);
      va_end(fmtargs);
#endif
    }
#endif

#ifdef MYALLOC_ENABLE_CS
#pragma omp critical(alloc_cs)
    {
#endif
      res = calloc(size, 1);
#ifdef MYALLOC_ENABLE_CS
    }
#endif

    if (size > 0 && res == NULL) {
      // spit error message and crash
#ifdef NDEBUG
      fprintf(stderr, "my_calloc error: calloc failed\n");
#else
      char new_msg[strlen(msg) + 1 + 200];
      strcpy(new_msg, "my_calloc warning: calloc failed: ");
      strcat(new_msg, msg);
      strcat(new_msg, "\n");
      va_start(fmtargs, msg);
      vfprintf(stderr, new_msg, fmtargs);
      va_end(fmtargs);
#endif
      exit(1);
    }
#ifdef MYALLOC_ENABLE_CRT
    else {
      crt_mem += size;
      if (counter != NULL)
	count_add(counter, (int64_t)size);
    }
  } // end of critical section
#endif

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

#ifdef DEBUG_MY_ALLOC
  {
    fprintf(stderr, "my_realloc: -%lld +%lld: ", (long long)old_size, (long long)size);
    char new_msg[strlen(msg) + 1 + 200];
    strcpy(new_msg, msg);
    strcat(new_msg, "\n");
    va_start(fmtargs, msg);
    vfprintf(stderr, new_msg, fmtargs);
    va_end(fmtargs);
  }
#endif

#ifdef MYALLOC_ENABLE_CRT
#pragma omp critical (cs_my_alloc)
  {
    if ((crt_mem - old_size) + size > max_mem) {
      if (!warned_max) {
	warned_max = true;
#ifdef NDEBUG
	fprintf(stderr, "my_realloc warning: exceeding maximum memory\n");
#else
	char new_msg[strlen(msg) + 1 + 200];
	strcpy(new_msg, "my_realloc warning: exceeding maximum memory: ");
	strcat(new_msg, msg);
	strcat(new_msg, "\n");
	va_start(fmtargs, msg);
	vfprintf(stderr, new_msg, fmtargs);
	va_end(fmtargs);
#endif
      }
    }
#endif

#ifdef MYALLOC_ENABLE_ALERT
    if ((long long)size - (long long)old_size > (long long)alert_mem) {
#ifdef NDEBUG
      fprintf(stderr, "my_realloc alert: size=%lld\n", (long long)size - (long long)old_size);
#else
      char new_msg[strlen(msg) + 1 + 200];
      sprintf(new_msg, "my_realloc alert: size=%lld: %s\n", (long long)size - (long long)old_size, msg);
      va_start(fmtargs, msg);
      vfprintf(stderr, new_msg, fmtargs);
      va_end(fmtargs);
#endif
    }
#endif

#ifdef MYALLOC_ENABLE_CS
#pragma omp critical(alloc_cs)
    {
#endif
      res = realloc(p, size);
#ifdef MYALLOC_ENABLE_CS
    }
#endif

    if (size > 0 && res == NULL) {
      // spit error message and crash
#ifdef NDEBUG
      fprintf(stderr, "my_realloc error: realloc failed\n");
#else
      char new_msg[strlen(msg) + 1 + 200];
      strcpy(new_msg, "my_realloc warning: realloc failed: ");
      strcat(new_msg, msg);
      strcat(new_msg, "\n");
      va_start(fmtargs, msg);
      vfprintf(stderr, new_msg, fmtargs);
      va_end(fmtargs);
#endif
      exit(1);
    }
#ifdef MYALLOC_ENABLE_CRT
    else {
      crt_mem -= old_size;
      crt_mem += size;
      if (counter != NULL)
	count_add(counter, (int64_t)size - (int64_t)old_size);
    }
  }
#endif

  return res;
}


static inline void
my_free(void * p, size_t size, count_t * counter, char const * msg = NULL, ...)
{
#ifndef NDEBUG
  va_list fmtargs;
#endif

#ifdef DEBUG_MY_ALLOC
  {
    fprintf(stderr, "my_free: -%lld: ", (long long)size);
    char new_msg[strlen(msg) + 1 + 200];
    strcpy(new_msg, msg);
    strcat(new_msg, "\n");
    va_start(fmtargs, msg);
    vfprintf(stderr, new_msg, fmtargs);
    va_end(fmtargs);
  }
#endif

#ifdef MYALLOC_ENABLE_CRT
#pragma omp critical (cs_my_alloc)
  {
#ifndef NDEBUG
    if (size > crt_mem) {
      fprintf(stderr, "my_free: crashing: p=%p size=%lld crt_mem=%lld counter.crt=%lld counter.max=%lld: ",
	      p, (long long)size, (long long)crt_mem, (long long)counter->crt, (long long) counter->max);
      char new_msg[strlen(msg) + 1 + 200];
      strcpy(new_msg, msg);
      strcat(new_msg, "\n");
      va_start(fmtargs, msg);
      vfprintf(stderr, new_msg, fmtargs);
      va_end(fmtargs);
    }
#endif
    assert(size <= crt_mem);
#endif

#ifdef MYALLOC_ENABLE_CS
#pragma omp critical(alloc_cs)
    {
#endif
      free(p);
#ifdef MYALLOC_ENABLE_CS
    }
#endif

#ifdef MYALLOC_ENABLE_CRT
    crt_mem -= size;
    if (counter != NULL)
      count_add(counter, -((int64_t)size));
  }
#endif
}





#endif
