/*	$Id$	*/

#ifndef _LOOKUP_H_
#define _LOOKUP_H_

#include <pthread.h>

#include "util.h"

#ifdef USE_PTHREADS
#ifdef USE_SPINLOCKS
typedef pthread_spinlock_t	lookup_lock_t;
#define lookup_lock_create(_x)	pthread_spin_init(_x, PTHREAD_PROCESS_SHARED)
#define lookup_lock_destroy(_x)	pthread_spin_destroy(_x)
#define lookup_lock_acquire(_x)	pthread_spin_lock(_x)
#define lookup_lock_release(_x)	pthread_spin_unlock(_x)
#else
typedef pthread_mutex_t		lookup_lock_t;
#define lookup_lock_create(_x)	pthread_mutex_init(_x, NULL)
#define lookup_lock_destroy(_x)	pthread_mutex_destroy(_x)
#define lookup_lock_acquire(_x)	pthread_mutex_lock(_x)
#define lookup_lock_release(_x)	pthread_mutex_unlock(_x)
#endif
#else
typedef int			lookup_lock_t;
#define lookup_lock_create(_x)	(0)
#define lookup_lock_destroy(_x) (0)
#define lookup_lock_acquire(_x) (0)
#define lookup_lock_release(_x) (0)
#endif

/*
 * A slight optimisation allows for a hash bucket with one element to directly
 * contain the key/value pair.
 *
 * Otherwise, another lookup structure (probably a red/black tree) exists and
 * is used for the lookup.
 */
struct hash_bucket {
	bool is_tree;
	lookup_lock_t lock;

	union {
		void *lookup_backend;
		struct {
			void *key;
			void *value;
		} kv_pair;
	} u_hb;
};

struct lookup {
	unsigned int count;
	lookup_lock_t count_lock;

	int (*keycompare_func)(void *, void *);
	unsigned int (*hash_func)(void *);
	struct hash_bucket *hash_table;
	size_t n_buckets;
	bool synchronised;
};

typedef struct lookup * lookup_t;

lookup_t lookup_create(unsigned int (*)(void *), int (*)(void *, void *), bool);
void lookup_destroy(lookup_t);
bool lookup_find(lookup_t, void *, void **, void **);
bool lookup_add(lookup_t, void *, void *);
bool lookup_remove(lookup_t, void *, void **, void **);
unsigned int lookup_count(lookup_t);
void lookup_iterate(lookup_t, void (*)(void *, void *, void *), void *);
void lookup_lock_key(lookup_t, void *);
void lookup_unlock_key(lookup_t, void *);

#endif	/* !_LOOKUP_H_ */
