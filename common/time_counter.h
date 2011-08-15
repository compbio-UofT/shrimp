#ifndef _TIME_COUNTER_H
#define _TIME_COUNTER_H

#include <stdlib.h>
#include "../common/util.h"

/*
 * We use two types of time counters:
 * rdtsc() is faster but less reliable;
 * gettimeinusecs() is slower but accurate
 */

typedef struct time_counter {
  long long int counter;
  int type;
} time_counter;

#ifdef NDEBUG
#define DEF_FAST_TIME_COUNTER	0
#else
#define DEF_FAST_TIME_COUNTER	1
#endif


static inline long long int
time_counter_check(time_counter const * tc)
{
  if (tc->type == 0)
    return rdtsc();
  else
    return gettimeinusecs();
}

static inline void
time_counter_add(time_counter * tc, long long int before, long long int after = -1)
{
  if (after == -1)
    after = time_counter_check(tc);
  if ((tc->type == 0 && after >= before) || tc->type == 1)
    tc->counter += after - before;
}

static inline double
time_counter_get_secs(time_counter const * tc)
{
  if (tc->type == 0)
    return (double)tc->counter / cpuhz();
  else
    return (double)tc->counter / 1.0e6;
}


#define TIME_COUNTER_START(tc)			\
  long long int before = time_counter_check(&(tc));

#define TIME_COUNTER_STOP(tc)			\
  time_counter_add(&(tc), before);


#endif
