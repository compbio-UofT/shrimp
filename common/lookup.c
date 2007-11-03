/*	$Id$	*/

#include <assert.h>
#include <pthread.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

//#include "atomic.h"
#include "lookup.h"
#include "red_black_tree.h"
#include "util.h"

#define USE_RBT
//#define USE_ATOMIC_OPS

/*
 * Abstract away our choice of storage, since we may well change to something
 * new and don't want to have to make api changes in mapreduce.c all the time.
 *
 * We use two layers. The first is a constant-sized hash table. This allows for
 * very fast lookups, but we can't know how many different keys will exist in
 * advance. To deal with collisions, we use red/black trees in each hash bucket
 * rather than a linked list. This should keep things from degenerating to
 * O(n), though sizing the hash table appropriately is still very important.
 * O(1) is still way faster than O(logn) for us when n is big.
 *
 * XXX - We probably would like a runtime switch.
 */

#ifdef USE_ATOMIC_OPS

static inline void
count_incr(lookup_t lu)
{

	/* XXX - test more expensive than lock? */
	if (lu->synchronised)
		x86_atomic_incr(&lu->count, 1);
	else
		lu->count++;
}

static inline void
count_decr(lookup_t lu)
{

	/* XXX - test more expensive than lock? */
	if (lu->synchronised)
		x86_atomic_decr(&lu->count, 1);
	else
		lu->count--;
}

#else

static inline void
count_incr(lookup_t lu)
{
	if (lu->synchronised) {
		(void)lookup_lock_acquire(&lu->count_lock);
		lu->count++;
		(void)lookup_lock_release(&lu->count_lock);
	} else
		lu->count++;
}

static inline void
count_decr(lookup_t lu)
{
	if (lu->synchronised) {
		(void)lookup_lock_acquire(&lu->count_lock);
		lu->count--;
		(void)lookup_lock_release(&lu->count_lock);
	} else
		lu->count--;
}

#endif

/*
 * Create a new lookup structure, with whatever backend we're compiled to
 * support.
 */
lookup_t
lookup_create(unsigned int (*hash_func)(void *),
    int (*keycompare_func)(void *, void *), bool synchronised)
{
	size_t i;
	struct lookup *lu;

	if (hash_func == NULL || keycompare_func == NULL)
		return (NULL);

	lu = malloc(sizeof(*lu));
	if (lu == NULL)
		return (NULL);

	memset(lu, 0, sizeof(*lu));

	lu->n_buckets = (64 * 1024 * 1024) / sizeof(struct hash_bucket);
	lu->hash_table = malloc(lu->n_buckets * sizeof(struct hash_bucket));
	if (lu->hash_table == NULL) {
		free(lu);
		return (NULL);
	}

	memset(lu->hash_table, 0, lu->n_buckets * sizeof(struct hash_bucket)); 

	lu->hash_func = hash_func;
	lu->keycompare_func = keycompare_func;
	lu->synchronised = synchronised;

	if (synchronised && lookup_lock_create(&lu->count_lock) != 0) {
		free(lu->hash_table);
		free(lu);
		return (NULL);
	}

	for (i = 0; synchronised && i < lu->n_buckets; i++) {
		int ret;

		ret = lookup_lock_create(&lu->hash_table[i].lock);
		if (ret != 0) {
			int j;

			for (j = 0; j < i; j++)
				(void)lookup_lock_destroy(
				    &lu->hash_table[i].lock);

			(void)lookup_lock_destroy(&lu->count_lock);
			free(lu->hash_table);
			free(lu);
			return (NULL);
		}
	}

	return (lu);
}

/*
 * We permit this to be called on a non-empty structure, so the caller had
 * better be sure that they won't leave any dangling pointers.
 *
 * The caller must also be sure that nobody could be holding locks on this
 * structure anymore.
 */
void
lookup_destroy(lookup_t lu)
{
	int i;

#ifdef USE_RBT
	for (i = 0; i < lu->n_buckets; i++) {
		if (lu->hash_table[i].is_tree) {
			assert(lu->hash_table[i].u_hb.lookup_backend != NULL);
			RBTreeDestroy(lu->hash_table[i].u_hb.lookup_backend);
		}
	}
#endif

	if (lu->synchronised)
		(void)lookup_lock_destroy(&lu->count_lock);

	for (i = 0; lu->synchronised && i < lu->n_buckets; i++)
		(void)lookup_lock_destroy(&lu->hash_table[i].lock);

	free(lu->hash_table);
	free(lu);
}

/*
 * Return the value associated with 'key', if it exists.
 */
bool
lookup_find(lookup_t lu, void *key, void **rkey, void **rvalue)
{
	bool ret;
	unsigned int hash;

	assert(key != NULL);

	hash = lu->hash_func(key) % lu->n_buckets;
	
	if (lu->hash_table[hash].is_tree == false) {
		if (lu->hash_table[hash].u_hb.kv_pair.key != NULL) {
			if (lu->keycompare_func(key,
			    lu->hash_table[hash].u_hb.kv_pair.key) == 0) {
				assert(lu->count > 0);

				if (rkey != NULL)
					*rkey =
					  lu->hash_table[hash].u_hb.kv_pair.key;
				
				if (rvalue != NULL)
					*rvalue =
					lu->hash_table[hash].u_hb.kv_pair.value;

				return (true);
			}
		}

		return (false);
	}
	
#ifdef USE_RBT
	assert(lu->hash_table[hash].u_hb.lookup_backend != NULL);
	ret = RBTreeFind(lu->hash_table[hash].u_hb.lookup_backend,
	    key, rkey, rvalue);
#endif

	return (ret);
}

/*
 * Add a key/value pair to the lookup structure.
 */
bool
lookup_add(lookup_t lu, void *key, void *value)
{
	unsigned int hash;

	assert(key != NULL);

	hash = lu->hash_func(key) % lu->n_buckets;

	/*
	 * Handle the optimised case first; if our hash table doesn't have
	 * collisions, don't even bother with secondary lookup structures.
	 */
	if (lu->hash_table[hash].is_tree == false) {
		if (lu->hash_table[hash].u_hb.kv_pair.key == NULL) {
			/*
			 * Bucket is empty. Add key/value pair.
			 */
			lu->hash_table[hash].u_hb.kv_pair.key = key;
			lu->hash_table[hash].u_hb.kv_pair.value = value;

			count_incr(lu);

			return (true);
		} else {
			/*
			 * Bucket has no tree, only one key/value pair.
			 * Create a new tree, add the old value, and fall
			 * through.
			 */
			void *old_key, *old_value, *lookup_backend;

			old_key = lu->hash_table[hash].u_hb.kv_pair.key;
			old_value = lu->hash_table[hash].u_hb.kv_pair.value;

#ifdef USE_RBT
			lookup_backend = RBTreeCreate(lu->keycompare_func);
			if (lookup_backend == NULL)
				return (false);

			if (RBTreeInsert(lookup_backend, old_key, old_value)
			    == NULL) {
				RBTreeDestroy(lookup_backend);
				return (false);
			}
#endif

			lu->hash_table[hash].u_hb.kv_pair.key = NULL;
			lu->hash_table[hash].u_hb.kv_pair.value = NULL;
			lu->hash_table[hash].is_tree = true;
			lu->hash_table[hash].u_hb.lookup_backend=lookup_backend;

			/* Fall through to add new key below... */
		}
	}

	/*
	 * Bucket contains a lookup structure, use it.
	 */

	/*
	 * Only one of each key allowed.
	 *
	 * XXX - Add this functionality to RBTreeInsert and avoid traversing
	 *	 twice.
	 */
	if (lookup_find(lu, key, NULL, NULL))
		return (false);

#ifdef USE_RBT
	assert(lu->hash_table[hash].u_hb.lookup_backend != NULL);
	if (RBTreeInsert(lu->hash_table[hash].u_hb.lookup_backend,
	    key, value) == NULL)
		return (false);
#endif

	count_incr(lu);	

	return (true);
}

/*
 * Remove a key from the lookup structure, returning the original key and value,
 * which may be pointers to comparable, but physically distinct keys/values.
 */
bool
lookup_remove(lookup_t lu, void *key, void **rkey, void **rvalue)
{
	bool ret;
	unsigned int hash;

	/* shut up, icc */
	ret = false;
	(void)ret;

	assert(key != NULL);

	hash = lu->hash_func(key) % lu->n_buckets;

	if (lu->hash_table[hash].is_tree == false) {
		if (lu->hash_table[hash].u_hb.kv_pair.key != NULL) {
			if (lu->keycompare_func(key,
			    lu->hash_table[hash].u_hb.kv_pair.key) == 0) {
				assert(lu->count > 0);

				if (rkey != NULL)
					*rkey =
					  lu->hash_table[hash].u_hb.kv_pair.key;
				
				if (rvalue != NULL)
					*rvalue =
					lu->hash_table[hash].u_hb.kv_pair.value;

				lu->hash_table[hash].u_hb.kv_pair.key = NULL;
				lu->hash_table[hash].u_hb.kv_pair.value = NULL;

				count_decr(lu);	

				return (true);
			}
		}

		return (false);
	}

	if (lookup_find(lu, key, rkey, rvalue) == false)
		return (false);

	assert(lu->count > 0);

#ifdef USE_RBT
	assert(lu->hash_table[hash].u_hb.lookup_backend != NULL);
	ret = RBDelete(lu->hash_table[hash].u_hb.lookup_backend, key);
	assert(ret);
#endif

	count_decr(lu);	

	return (true);
}

/*
 * Return the number of elements in our structure.
 */
unsigned int
lookup_count(lookup_t lu)
{

	return (lu->count);
}

/*
 * Process all key,value pairs in our structure with the provided
 * callback function.
 *
 * NB: This function assumes that it will be called when no synchonrisation is
 * needed.
 */
void
lookup_iterate(lookup_t lu, void (*iter_func)(void *, void *, void *),
    void *arg)
{
	size_t i, count;

	/*
	 * Since mr5 tries to merge into a single structure and keys are
	 * naturally related to their position in the hash table, we want to
	 * avoid as much thread contention as possible.
	 *
	 * Starting from 0 for all threads would be foolish, as we'd be
	 * potentially contending the entire way through. Instead, start at some
	 * random offset and hope that this keeps us out of trouble. 
	 */
	i = random() % lu->n_buckets;

	for (count = 0; count < lu->n_buckets; count++) {
		if (lu->hash_table[i].is_tree == false) {
			if (lu->hash_table[i].u_hb.kv_pair.key != NULL)
				iter_func(arg,
				    lu->hash_table[i].u_hb.kv_pair.key,
				    lu->hash_table[i].u_hb.kv_pair.value);
		} else {
#ifdef USE_RBT

			assert(lu->hash_table[i].u_hb.lookup_backend != NULL);
			RBTreeIterate(lu->hash_table[i].u_hb.lookup_backend,
			    iter_func, arg);
#endif
		}

		i = (i + 1) % lu->n_buckets;
	}
}

void
lookup_lock_key(lookup_t lu, void *key)
{
	unsigned int hash = 0;

	/* shut up, icc */
	(void)hash;

	assert(lu->synchronised);
	hash = lu->hash_func(key) % lu->n_buckets;
	(void)lookup_lock_acquire(&lu->hash_table[hash].lock);
}

void
lookup_unlock_key(lookup_t lu, void *key)
{
	unsigned int hash = 0;

	/* shut up, icc */
	(void)hash;

	assert(lu->synchronised);
	hash = lu->hash_func(key) % lu->n_buckets;
	(void)lookup_lock_release(&lu->hash_table[hash].lock);
}
