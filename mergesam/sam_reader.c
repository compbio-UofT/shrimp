#include <stdio.h>
#include <assert.h>
#include <string.h>
#include <stdlib.h>
#include "mergesam_heap.h"
#include "sam_reader.h"

bool found_sam_headers;

static inline void pp_ll_zero(pp_ll * ll) {
	ll->head=NULL;
	ll->tail=NULL;
	ll->length=0;
}

static inline void pp_ll_append(pp_ll* ll,pretty * pa) {
	if ( ll->length==0 ) {
		ll->head=pa;
		ll->tail=pa;
		pa->next=NULL;
		ll->length=1;
	} else {
		ll->tail->next=pa;
		ll->tail=pa;
		pa->next=NULL;
		ll->length++;
	}
}

static inline void pp_ll_set(pp_ll * ll,pretty *  pa) {
	ll->length=1;
	ll->head=pa;
	ll->tail=pa;
	pa->next=NULL;
}

static inline void pp_ll_combine_and_check_single(pp_ll * m_ll,pp_ll ** ll,int offset,pretty ** unaligned_pa,heap_pa * h) {
	//heap_pa * h=thread_heaps+omp_get_thread_num();
	h->load=0;
	heap_pa_elem e;
	pretty * pa=NULL;
	int i;
	for (i=0; i<options.number_of_sam_files; i++) {
		pp_ll * local_ll = ll[i]+offset;
		if (local_ll->length>0) {
			pa = local_ll->head;
			while (pa!=NULL) {
				e.score=pa->score+(pa->mate_pair!=NULL ? pa->mate_pair->score : 0);
				e.isize_score=MAX_INT32;
				e.rest=pa;
				if (options.strata) {
					heap_pa_insert_bounded_strata(h,&e);
				} else {
					heap_pa_insert_bounded(h,&e);
				}
				pa=pa->next;
			}
		}
	}
	if (h->load==0 || (options.max_alignments>0 && h->load>options.max_alignments)) {
		//m_ll->length=0;
	} else {
		assert(options.max_outputs>0);
		i=h->load>options.max_outputs ? 1 : 0;
		if (m_ll->length==0) {
			m_ll->head=h->array[i].rest;
			m_ll->tail=h->array[h->load-1].rest;
		} else {
			m_ll->tail->next=h->array[i].rest;
			m_ll->tail=h->array[h->load-1].rest;
		}
		for (; i<h->load-1; i++) {
			h->array[i].rest->next=h->array[i+1].rest;
		}
		h->array[i].rest->next=NULL;
		m_ll->length+=h->load- (h->load>options.max_outputs ? 1 : 0);
	}
	if (options.sam_unaligned && h->load>0) {
		*unaligned_pa=h->array[0].rest;
	}
}

static inline void remove_offending_fields(pretty * pa) {
	pa->has_ih=false;
	pa->has_hi=false;
	pa->has_h0=false;
	pa->has_h1=false;
	pa->has_h2=false;
}

void pp_ll_combine_and_check(pp_ll * m_ll,pp_ll ** ll,heap_pa *h) {
	//want to make a heap
	//check first if paired of not
	assert(m_ll->length==0);
	//heap_pa * h=thread_heaps+omp_get_thread_num();
	h->load=0;
	heap_pa_elem e;
	pretty * pa = NULL;
	pretty * unaligned_pa=NULL;
	if (options.paired) {
		assert(!options.unpaired);
		//check paired first, and then go lower
		int i;
		for (i=0; i<options.number_of_sam_files; i++) {
			pp_ll * paired_ll = ll[i]+PAIRED;
			if (paired_ll->length>0) {
				pa = paired_ll->head;
				while (pa!=NULL) {
					e.score=pa->score+pa->mate_pair->score;
					assert(pa->mapped && pa->mate_pair!=NULL);
					assert(pa->mate_pair->mapped &&  pa->mp_mapped);
					e.isize_score=options.expected_insert_size >= 0 ? abs(options.expected_insert_size-pa->isize) : 0;
					e.rest=pa;
					if (options.strata) {
						heap_pa_insert_bounded_strata(h,&e);
					} else {
						heap_pa_insert_bounded(h,&e);
					}
					pa=pa->next;
				}
			}
		}
		if (h->load==0 || (options.max_alignments>0 && h->load>options.max_alignments)) {
			//m_ll->length=0;
		} else {
			assert(options.max_outputs>0);
			i=h->load>options.max_outputs ? 1 : 0;
			m_ll->head=h->array[i].rest;
			m_ll->tail=h->array[h->load-1].rest;
			for (; i<h->load-1; i++) {
				h->array[i].rest->next=h->array[i+1].rest;
			}
			h->array[i].rest->next=NULL;
			m_ll->length=h->load- (h->load>options.max_outputs ? 1 : 0);
		}
		if (options.sam_unaligned && h->load>0) {
			unaligned_pa=h->array[0].rest;
		}
	}
	if (h->load==0) {
		if (options.paired && options.half_paired) {
			pp_ll_combine_and_check_single(m_ll,ll,FIRST_LEG,&unaligned_pa,h);
			pp_ll_combine_and_check_single(m_ll,ll,SECOND_LEG,&unaligned_pa,h);
		} else if (options.unpaired) {
			pp_ll_combine_and_check_single(m_ll,ll,UNPAIRED,&unaligned_pa,h);
		}
	}
	if (m_ll->length==0 && options.sam_unaligned) {
		if (unaligned_pa==NULL) {
			pp_ll_combine_and_check_single(m_ll,ll,UNMAPPED,&unaligned_pa,h);
		}
		if (unaligned_pa!=NULL) {
			m_ll->length=1;
			m_ll->head=unaligned_pa;
			m_ll->tail=unaligned_pa;
			unaligned_pa->next=NULL;
			pretty_from_aux_inplace(unaligned_pa);
			pretty_print_sam_unaligned(unaligned_pa,true);
			if (unaligned_pa->paired_sequencing) {
				pretty_from_aux_inplace(unaligned_pa->mate_pair);
				pretty_print_sam_unaligned(unaligned_pa->mate_pair,true);
			}
		}
	} else {
		pa=m_ll->head;
		while(pa!=NULL) {
			pretty_from_aux_inplace(pa);
			remove_offending_fields(pa);
			pretty_print_sam_update(pa, true);
			//revert_sam_string(pa);
			if (pa->mate_pair!=NULL) {
				pretty_from_aux_inplace(pa->mate_pair);
				remove_offending_fields(pa->mate_pair);
				pretty_print_sam_update(pa->mate_pair, true);
				//revert_sam_string(pa->mate_pair);
			}
			pa=pa->next;
		}
	}	
	return;	
	
		
}

static inline void pp_ll_append_and_check(pp_ll* ll,pretty * pa) {
	if (pa->paired_sequencing) {
		if (pa->proper_pair) {
			assert(pa->mapped && pa->mp_mapped);
			//both sides mapped
			pp_ll_append(ll+PAIRED,pa);
		} else if (options.half_paired && (pa->mapped || pa->mp_mapped)) {
			//lets figure out which leg
			if (pa->mapped) {
				if (pa->first_in_pair) {
					pp_ll_append(ll+FIRST_LEG,pa);
				} else {
					pp_ll_append(ll+SECOND_LEG,pa);
				}
			} else {
				if (pa->first_in_pair) {
					pp_ll_append(ll+SECOND_LEG,pa);
				} else {
					pp_ll_append(ll+FIRST_LEG,pa);
				}
			}
		} else if (options.sam_unaligned && !pa->mapped && !pa->mp_mapped) {
			//unmapped paired
			pp_ll_append(ll+UNMAPPED,pa);
		}
	} else {
		if (pa->mapped) {
			//unpaired and mapped
			pp_ll_append(ll+UNPAIRED,pa);
		} else if (options.sam_unaligned) {
			//unpaired and unmapped
			pp_ll_append(ll+UNMAPPED,pa);
		}
	}
}
void grow_sam_pretty(sam_reader * sr) {
	size_t new_pretty_stack_size=(size_t)(sr->pretty_stack_size*GROWTH_FACTOR+1);
	assert(new_pretty_stack_size>sr->pretty_stack_size);
	fprintf(stderr,"Growing %lu to %lu entries\n",sr->pretty_stack_size,new_pretty_stack_size);
	size_t old_size = sr->pretty_stack_size;
	char * old = (char*)sr->pretty_stack;
	sr->pretty_stack = (pretty*)realloc(sr->pretty_stack,new_pretty_stack_size*sizeof(pretty));
	if (sr->pretty_stack==NULL) {
		fprintf(stderr,"Failed realloc of pretty stack!\n");
		exit(1);
	}
	memset(sr->pretty_stack + old_size , 0, (new_pretty_stack_size-old_size)*sizeof(pretty));
	sr->pretty_stack_size=new_pretty_stack_size;	
	char * new_p = (char*)sr->pretty_stack;
	int i;
	for (i=0; i<LL_ALL*options.read_rate; i++) {
		pp_ll * ll= sr->pp_lls+i;
		if (ll->head!=NULL) {
			ll->head=(pretty*)(((char*)ll->head)+(new_p-old));
		}
		if (ll->tail!=NULL) {
			ll->tail=(pretty*)(((char*)ll->tail)+(new_p-old));
		}
	}
	size_t ui;
	for (ui=0; ui<old_size; ui++) {
		pretty * pa=sr->pretty_stack+ui;
		if (pa->next!=NULL) {
			pa->next=(pretty*)(((char*)pa->next)+(new_p-old));
		}
		if (pa->mate_pair!=NULL) {
			pa->mate_pair=(pretty*)(((char*)pa->mate_pair)+(new_p-old));
		}
	}
}

static inline void revert_sam_string(pretty * pa) {
	char * nill = (char*)memchr(pa->sam_string,'\0',pa->sam_string_length);
	while (nill!=NULL) {
		*nill='\t';
		nill=(char*)memchr(pa->sam_string,'\0',pa->sam_string_length-(nill-pa->sam_string));
	}
}


static inline int sam_header_field_sort(char * a , char * b, char * check_for) {
	bool a_has=(a[1]==check_for[0] && a[2]==check_for[1]);
	bool b_has=(b[1]==check_for[0] && b[2]==check_for[1]);
	if (a_has && b_has) {
		return strcmp(a,b);
	} else if (a_has) {
		return -1;
	} else if (b_has) {
		return 1;
	} else {
		return 0;
	}
}

int sam_header_sort(const void * a,const void * b) {
	char * s1=((char**)a)[0];
	char * s2=((char**)b)[0];
	if (strlen(s1)<4 || strlen(s2)<4) {
		fprintf(stderr,"Failed to sort sam header line!\n");
		fprintf(stderr,"%s and %s\n",s1,s2);
		exit(1);
	}
	if (s1[0]!='@' || s2[0]!='@') {
		fprintf(stderr,"These two lines are not sam-header lines! %s and %s\n",s1,s2);
		exit(1);
	}
	int ret=0;
	//check for HD
	ret=sam_header_field_sort(s1,s2,"HD");
	if (ret!=0) { return ret; };
	//check for SQ
	ret=sam_header_field_sort(s1,s2,"SQ");
	if (ret!=0) { return ret; };
	//check for RG
	ret=sam_header_field_sort(s1,s2,"RG");
	if (ret!=0) { return ret; };
	//check for PG
	ret=sam_header_field_sort(s1,s2,"PG");
	if (ret!=0) { return ret; };
	//check for PG
	ret=sam_header_field_sort(s1,s2,"CO");
	return ret;
	
	
	
}


void sam_close(sam_reader * sr) {
	fb_close(sr->fb);
	free(sr->pretty_stack);
	free(sr->pp_lls);
	free(sr);
}

sam_reader * sam_open(char * sam_filename,fastx_readnames * fxrn) {
	sam_reader * sr = (sam_reader*)malloc(sizeof(sam_reader));
	if (sr==NULL) {
		fprintf(stderr,"thread_open_sam : failed to allocate memory for thread info structure\n");
		exit(1);
	}
	sr->fb=fb_open(sam_filename,options.buffer_size,options.read_size);
	sr->pretty_stack_filled=0;
	sr->pretty_stack_size=fxrn->reads_inmem*10;
	sr->pretty_stack=(pretty*)malloc(sizeof(pretty)*sr->pretty_stack_size);
	if (sr->pretty_stack==NULL) {
		fprintf(stderr,"Failed to allocate memory for pretty_stack\n");
		exit(1);
	}
	memset(sr->pretty_stack,0,sizeof(pretty)*sr->pretty_stack_size);
	sr->pp_lls=(pp_ll*)malloc(sizeof(pp_ll)*LL_ALL*options.read_rate);
	if (sr->pp_lls==NULL) {
		fprintf(stderr,"Failed to allocate memory for pp_lls.\n");
	}
	memset(sr->pp_lls,0,sizeof(pp_ll)*LL_ALL*options.read_rate);
	return sr;
}
void parse_sam(sam_reader * sr,size_t read_amount,fastx_readnames * fxrn) {
	char * current_newline=NULL;
	size_t last_tested=0;
	while (sr->fb->unseen_end!=sr->fb->unseen_start && last_tested<read_amount) {
        	size_t unseen_end_mod=sr->fb->unseen_end%sr->fb->size;
        	size_t unseen_start_mod=sr->fb->unseen_start%sr->fb->size;
		size_t space=(unseen_start_mod>=unseen_end_mod ? sr->fb->size : unseen_end_mod) - unseen_start_mod;
		current_newline=(char*)memchr(sr->fb->base+unseen_start_mod,'\n',space);
		if (current_newline==NULL) {
			sr->fb->unseen_start+=space;
		} else {
			size_t current_newline_index=current_newline-sr->fb->base;
			if (current_newline_index==0 || sr->fb->base[current_newline_index-1]=='\0') {
				*current_newline='\0';
				current_newline_index++;
				for(; current_newline_index<sr->fb->size && sr->fb->base[current_newline_index]=='\0'; current_newline_index++) {sr->fb->base[current_newline_index]='\0';};
				sr->fb->unseen_start+=current_newline_index-unseen_start_mod;
				unseen_start_mod=sr->fb->unseen_start%sr->fb->size;
			} else {
				//fprintf(stderr,"Using %lu\n",unseen_start_mod);
				char * line=sr->fb->base+unseen_start_mod;
				assert(sr->pretty_stack_size>sr->pretty_stack_filled);
				if (line[0]!='@') {
					size_t length_of_string = current_newline-line;
					if (sr->pretty_stack_size==sr->pretty_stack_filled) {	
						grow_sam_pretty(sr);
					}
					for (; last_tested<read_amount; last_tested++) {
						char * first_tab = (char*)memchr(line,'\t',length_of_string);
						if (first_tab==NULL) {
							fprintf(stderr,"CANNOT FIND FIRST TAB!\n");
							exit(1);
						}
						size_t compare_length = first_tab-line;
						char * hit_list_read_name=fxrn->read_names+(last_tested+fxrn->reads_seen)*SIZE_READ_NAME*sizeof(char);
						
						if (strncmp(line,hit_list_read_name,compare_length)==0 ) {//&& !isdigit(hit_list_read_name[compare_length])) {	
							pretty * pa=sr->pretty_stack+sr->pretty_stack_filled++;
							pa->read_id=last_tested;
							pa->sam_header=false;
							pretty_from_string_inplace(line,length_of_string,pa);
							//fprintf(stderr,"Found a hit for read |%s| vs |%s|\n",sr->pretty_stack[sr->pretty_stack_filled-1].read_name,read_names+(last_tested+reads_offset)*SIZE_READ_NAME*sizeof(char));
							//check whether paired, unpaired or first or second leg
							//fprintf(stdout,"FOUND score is %d\n",pa->score);
							if (pa->paired_sequencing) {
								options.paired=true;
								if (sr->pretty_stack_filled!=1 && !(pa-1)->sam_header && (pa-1)->mate_pair==NULL) {
									(pa-1)->mate_pair=pa; pa->mate_pair=pa-1;
									//fprintf(stderr,"Added as pair!\n");
									pp_ll_append_and_check(sr->pp_lls+LL_ALL*pa->read_id,pa);
								} else {
									pa->mate_pair=NULL;
								}
							} else {
								options.unpaired=true;
								pp_ll_append_and_check(sr->pp_lls+LL_ALL*pa->read_id,pa);
							}
							*current_newline='\0';
							break;
						} 	
					}
					
				} else {
					found_sam_headers=true;
					pretty * pa=sr->pretty_stack+sr->pretty_stack_filled++;
					pa->sam_string=line;
					*current_newline='\0';
					pp_ll_append(sr->sam_headers,pa);
					pa->sam_header=true;
				}
				if (last_tested!=read_amount) {
					sr->fb->unseen_start+=current_newline_index-unseen_start_mod+1;
					unseen_start_mod=sr->fb->unseen_start%sr->fb->size;
				}
			}
			sr->fb->exhausted=false;
		}
	}
	if (last_tested<read_amount && (sr->fb->frb.eof==0 || sr->fb->unseen_end!=sr->fb->unseen_start)) {
		fprintf(stderr,"failed to read into buffer, please increase buffer size!\n");
		exit(1);
	}
}
