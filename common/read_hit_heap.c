
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include "read_hit_heap.h"


static inline int compare_scores(read_hit_heap_e * a, read_hit_heap_e * b) {
	if (a->score<b->score || (a->score==b->score && a->isize_score>b->isize_score)) {
		return 1;
	} else if (a->score==b->score && a->isize_score==b->isize_score) {
		return 0;
	} 
	return -1;
}

static bool rhh_compare_scores(read_hit_heap_e * a, read_hit_heap_e * b) {
	return compare_scores(a,b)==1 ? true  : false;
}


/* return 1 if a<b, 0 if a==b and -1 otherwise */

static inline int compare_read_hit(read_hit *a, read_hit * b) {
	if (a->g_off<b->g_off) {
		return 1;
	} else if (a->g_off==b->g_off) {
		if (a->cn<b->cn) {
			return 1;	
		} else if (a->cn==b->cn) {
			if (a->st<b->st) {
				return 1;
			} else if (a->st==b->st) {
				if (a->gen_st<b->gen_st) {
					return 1;
				} else if (a->gen_st==b->gen_st) {
					return 0;
				}
			}
		}
	}
	return -1;
}

static inline bool equal_alignments(read_hit_heap_e * a, read_hit_heap_e * b) {
	assert(a->rh[0]!=NULL);
	assert(b->rh[0]!=NULL);
	int result = compare_read_hit(a->rh[0],b->rh[0]);
	if (result==0 && a->rh[1]!=NULL && b->rh[1]!=NULL) {
		result = compare_read_hit(a->rh[1],b->rh[1]);
	}
	if (result==0) {
		return true;
	}
	return false;
}

static bool rhh_compare_alignments(read_hit_heap_e * a, read_hit_heap_e * b) {
	assert(a->rh[0]!=NULL);
	assert(b->rh[0]!=NULL);
	int result = compare_read_hit(a->rh[0],b->rh[0]);
	if (result==1) { 
		return true;
	} else if (result==-1) {
		return false;
	} else if (a->rh[1]!=NULL && b->rh[1]!=NULL) {
		result = compare_read_hit(a->rh[1],b->rh[1]);
		if (result==1) {
			return true;
		} else if (result==-1) {
			return false;
		} else {
			return compare_scores(a,b)==1 ? true : false;
		}
	}
	return false;
}

void read_hit_heap_init(read_hit_heap * rhh, uint32_t capacity,bool alignment_compare) {
	assert(capacity>0);
	assert(rhh != NULL);
	rhh->array = (struct read_hit_heap_e *)malloc(capacity * sizeof(struct read_hit_heap_e)); 
	if (rhh->array==NULL) { 
		fprintf(stderr,"Failed to allocate heap!\n"); 
		exit(1); 
	} 
	rhh->capacity = capacity;
	rhh->load = 0;
	rhh->alignment_compare=alignment_compare;
	if (alignment_compare) {
		rhh->cmp=rhh_compare_scores;
	} else {
		rhh->cmp=rhh_compare_alignments;
	}
}

void read_hit_heap_destroy(struct read_hit_heap * rhh) {
	assert(rhh != NULL);
	free(rhh->array);
}

static inline void read_hit_heap_percolate_up(struct read_hit_heap * rhh, uint32_t node) {
	read_hit_heap_e tmp;
	uint32_t parent;
	parent = node / 2;
	while (node > 1 && rhh->cmp(rhh->array+node-1,rhh->array+parent-1)) {
		tmp = rhh->array[parent-1];
		rhh->array[parent-1] = rhh->array[node-1];
		rhh->array[node-1] = tmp;
		node = parent;
		parent = node / 2;
	}
}

static inline void read_hit_heap_percolate_down(struct read_hit_heap * rhh, int32_t node) {
	read_hit_heap_e tmp;
	int32_t left, right, max;
	do {
		left = node * 2;
		right = left + 1;
		max = node;
		if (left <= rhh->load && rhh->cmp(rhh->array+left-1,rhh->array+node-1)) 
		max = left;
		if (right <= rhh->load && rhh->cmp(rhh->array+right-1,rhh->array+max-1)) 
		max = right;
		if (max == node)
		break;
		tmp = rhh->array[max-1];
		rhh->array[max-1] = rhh->array[node-1];
		rhh->array[node-1] = tmp;
		node = max;
	} while (1);
}

void read_hit_heap_extract_min(struct read_hit_heap * rhh, struct read_hit_heap_e * dest)  {
	assert(rhh != NULL && rhh->load > 0);
	*dest = rhh->array[0];
	rhh->load--;
	if (rhh->load > 0) {
		rhh->array[0] = rhh->array[rhh->load];
		read_hit_heap_percolate_down(rhh, 1);
	}
}

void read_hit_heap_replace_min(struct read_hit_heap * rhh, struct read_hit_heap_e * e) {
	assert(rhh != NULL && rhh->load > 0);
	rhh->array[0] = *e;
	read_hit_heap_percolate_down(rhh, 1);
}

void read_hit_heap_get_min(struct read_hit_heap * rhh, struct read_hit_heap_e * dest) {
	assert(rhh != NULL && rhh->load > 0);
	*dest = rhh->array[0];
}

void read_hit_heap_insert(struct read_hit_heap * rhh, struct read_hit_heap_e * e) {
	assert(rhh != NULL && rhh->load < rhh->capacity);
	rhh->array[rhh->load] = *e;
	rhh->load++;
	read_hit_heap_percolate_up(rhh, rhh->load);
}

void read_hit_heap_insert_bounded(struct read_hit_heap * rhh, struct read_hit_heap_e *e) {
	if (rhh->load<rhh->capacity) {
		if (rhh->load==0) {
			rhh->array[0]=*e;
			rhh->load=1;
		} else {
			read_hit_heap_insert(rhh,e);
		}
	} else if (rhh->cmp(rhh->array+0,e)) {
		assert(rhh->load==rhh->capacity);
		read_hit_heap_replace_min(rhh,e);
	}
}

void read_hit_heap_insert_bounded_strata(struct read_hit_heap * rhh, struct read_hit_heap_e *e) {
	if (rhh->load==0) {
		rhh->array[0]=*e;
		rhh->load=1;
	} else if (rhh->cmp(e,rhh->array+0)) {
		return;
	} else if (rhh->cmp(rhh->array+0,e)) {
		//means that we should wipe heap and take the new one
		rhh->load=1;
		rhh->array[0]=*e;
	} else if (rhh->load<rhh->capacity) {
		//equal
		rhh->array[rhh->load++]=*e;
	}
}


void read_hit_heap_heapify(struct read_hit_heap * rhh){
	assert(rhh != NULL);
	uint32_t node;
	for (node = rhh->load / 2; node >= 1; node--) {
		read_hit_heap_percolate_down(rhh, node);
	}
}

void read_hit_heap_alignment_compare(struct read_hit_heap *rhh) {
	if (rhh->load!=0) {
		assert(1==0);
	}
	rhh->alignment_compare=true;
}

void read_hit_heap_score_compare(struct read_hit_heap *rhh) {
	if (rhh->load!=0) {
		assert(1==0);
	}
	rhh->alignment_compare=false;
}

void read_hit_heap_sort(struct read_hit_heap * rhh) {
	assert(rhh!=NULL);
	uint32_t orig_load=rhh->load;
	read_hit_heap_e tmp;

	while(rhh->load>1) {
		tmp=rhh->array[0];
		rhh->array[0]=rhh->array[rhh->load-1];
		rhh->array[rhh->load-1]=tmp;
		rhh->load--;
		read_hit_heap_percolate_down(rhh,1);
	}
	rhh->load=orig_load;
}

//precondition must be sorted by alignment
static inline void read_hit_heap_remove_dups(struct read_hit_heap * rhh) {
	assert(rhh->alignment_compare==true);
	read_hit_heap_sort(rhh);
	int32_t i=0; int32_t index=0;
	while (i<rhh->load) {
		int32_t end_of_match;
		for (end_of_match=i; (end_of_match+1)<rhh->load && equal_alignments(rhh->array+i,rhh->array+end_of_match+1); end_of_match++);
		if (index!=i) {
			rhh->array[index]=rhh->array[i];
		}
		index++;
		i=end_of_match+1;
	}
	rhh->load=index;
	rhh->alignment_compare=false;
}
void read_hit_heap_remove_dups_and_sort(struct read_hit_heap * rhh) {
	read_hit_heap_remove_dups(rhh);
	read_hit_heap_heapify(rhh);
	read_hit_heap_sort(rhh);	
}
void read_hit_heap_remove_dups_and_sort_strata(struct read_hit_heap * rhh) {
	read_hit_heap_remove_dups(rhh);
	int32_t orig_load=rhh->load;
	rhh->load=0;
	read_hit_heap_e tmp;
	int32_t i;
	for (i=0; i<orig_load; i++) {
		tmp=rhh->array[i];
		read_hit_heap_insert_bounded_strata(rhh,&tmp);
	}
	read_hit_heap_sort(rhh);	
}
