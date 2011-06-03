#include <stdio.h>
#include <assert.h>
#include <string.h>
#include <stdlib.h>
#include "mergesam_heap.h"
#include "sam_reader.h"
#include "../common/util.h"

bool found_sam_headers;

static inline int pa_to_mapq(pretty * pa) {
	if (!pa->has_z1 || !pa->has_z0) {
		//return 255;
		return 0;
	}
	if (pa->z1==0.0) {
		return 255;
	}
	return qv_from_pr_corr(pa->z0/pa->z1);
}

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
	//fprintf(stderr,"%p length of %lu\n",ll,ll->length);
}

static inline void pp_ll_set(pp_ll * ll,pretty *  pa) {
	ll->length=1;
	ll->head=pa;
	ll->tail=pa;
	pa->next=NULL;
}

static inline void pp_ll_combine_and_check_heap(pp_ll * m_ll,pp_ll ** ll,int offset,pretty ** unaligned_pa,heap_pa * h,double * z1s,bool paired) {
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
				if (paired) {
					e.score=pa->score+pa->mate_pair->score;
					assert(pa->mapped && pa->mate_pair!=NULL);
					assert(pa->mate_pair->mapped &&  pa->mp_mapped);
					e.isize_score=options.expected_insert_size >= 0 ? abs(options.expected_insert_size-pa->isize) : 0;
				} else {
					//Add in MAPQ code to grab stats here?
					e.score=pa->score+(pa->mate_pair!=NULL ? pa->mate_pair->score : 0);
					e.isize_score=MAX_INT32;
				}
				e.rest=pa;
				if (pa->has_z1 && z1s!=NULL) {
					z1s[pa->fileno]=pa->z1;
				}
				if (options.strata) {
					heap_pa_insert_bounded_strata(h,&e);
				} else {
					heap_pa_insert_bounded(h,&e);
				}
				pa=pa->next;
			}
		}
	}
	if (h->load>0 && (options.max_alignments==-1 || h->load<=options.max_alignments)) {
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
		m_ll->length+= h->load - (h->load>options.max_outputs ? 1 : 0);
	}
	if ((options.sam_unaligned || options.unaligned_fastx) && h->load>0) {
		*unaligned_pa=h->array[0].rest;
	}
}
static inline void pp_ll_combine_and_check_heap_single(pp_ll * m_ll,pp_ll ** ll,int offset,pretty ** unaligned_pa,heap_pa * h,double * z1s) {
	pp_ll_combine_and_check_heap(m_ll,ll,offset,unaligned_pa,h,z1s,false);
}
static inline void pp_ll_combine_and_check_heap_paired(pp_ll * m_ll,pp_ll ** ll,int offset,pretty ** unaligned_pa,heap_pa * h,double * z1s) {
	pp_ll_combine_and_check_heap(m_ll,ll,offset,unaligned_pa,h,z1s,true);
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
	double z1s[LL_ALL*options.number_of_sam_files];
	memset(z1s,0,sizeof(double)*LL_ALL*options.number_of_sam_files);
	double z1_sums[LL_ALL];
	memset(z1_sums,0,sizeof(double)*LL_ALL);	
	if (options.paired) {
		pp_ll_combine_and_check_heap_paired(m_ll,ll,PAIRED,&unaligned_pa,h,z1s+PAIRED*options.number_of_sam_files);
		if (options.half_paired) { 
			h->load=0;
			pp_ll_combine_and_check_heap_single(m_ll,ll,FIRST_LEG,&unaligned_pa,h,z1s+FIRST_LEG*options.number_of_sam_files);
			h->load=0;
			pp_ll_combine_and_check_heap_single(m_ll,ll,SECOND_LEG,&unaligned_pa,h,z1s+SECOND_LEG*options.number_of_sam_files);
		}
	} else if (options.unpaired) {
		pp_ll_combine_and_check_heap_single(m_ll,ll,UNPAIRED,&unaligned_pa,h,z1s+UNPAIRED*options.number_of_sam_files);
	}
	int i;
	for (i=0; i<options.number_of_sam_files; i++) {
		z1_sums[PAIRED]+=z1s[PAIRED*options.number_of_sam_files+i];
		z1_sums[FIRST_LEG]+=z1s[FIRST_LEG*options.number_of_sam_files+i];
		z1_sums[SECOND_LEG]+=z1s[SECOND_LEG*options.number_of_sam_files+i];
		z1_sums[UNPAIRED]+=z1s[UNPAIRED*options.number_of_sam_files+i];
	}
	if (m_ll->length==0 && (options.sam_unaligned || options.unaligned_fastx)) {
		if (unaligned_pa==NULL) {
			if (!options.half_paired) { 
				h->load=0;
				pp_ll_combine_and_check_heap_single(m_ll,ll,FIRST_LEG,&unaligned_pa,h,z1s+FIRST_LEG*options.number_of_sam_files);
				h->load=0;
				pp_ll_combine_and_check_heap_single(m_ll,ll,SECOND_LEG,&unaligned_pa,h,z1s+SECOND_LEG*options.number_of_sam_files);
			}
			h->load=0;
			pp_ll_combine_and_check_heap_single(m_ll,ll,UNMAPPED,&unaligned_pa,h,NULL);
		}
		if (unaligned_pa!=NULL) {
			m_ll->length=1;
			m_ll->head=unaligned_pa;
			m_ll->tail=unaligned_pa;
			unaligned_pa->next=NULL;
			
			pretty_from_aux_inplace(unaligned_pa);
			if (options.unaligned_fastx) {
				pretty_print_sam_fastx(unaligned_pa, true);
			} else {
				pretty_print_sam_unaligned(unaligned_pa,true);
			}
			if (unaligned_pa->paired_sequencing) {
				pretty_from_aux_inplace(unaligned_pa->mate_pair);
				if (options.unaligned_fastx) {
					pretty_print_sam_fastx(unaligned_pa->mate_pair, true);
				} else {
					pretty_print_sam_unaligned(unaligned_pa->mate_pair,true);
				}
			}
		}
	} else if (!options.unaligned_fastx) {
		pa=m_ll->head;
		while(pa!=NULL) {
			if (pa->has_z1) {
				if (pa->paired_sequencing) {
					if (pa->proper_pair) {
						pa->z1=z1_sums[PAIRED];
						pa->mate_pair->z1=z1_sums[PAIRED];
					} else if (pa->first_in_pair) {
						pa->z1=z1_sums[FIRST_LEG];
						pa->mate_pair->z1=z1_sums[SECOND_LEG];
					} else if (pa->second_in_pair) {
						pa->z1=z1_sums[SECOND_LEG];
						pa->mate_pair->z1=z1_sums[FIRST_LEG];
					} else {
						assert(1==0);
					}
					if (pa->mate_pair->mapped) {
						pa->mate_pair->mapq=pa_to_mapq(pa->mate_pair);
					}
				} else {
					pa->z1=z1_sums[UNPAIRED];
				}
				pa->mapq=pa_to_mapq(pa);
			}
			//pretty_from_aux_inplace(pa);
			if (options.aligned_fastx) {
				pa->next=NULL;
				pretty_from_aux_inplace(pa);
				pretty_print_sam_fastx(pa, true);
			} else {
				remove_offending_fields(pa);
				pretty_print_sam_update(pa, true);
			}
			//revert_sam_string(pa);
			if (pa->mate_pair!=NULL) {
				assert(pa->mate_pair->mate_pair!=NULL);
				if (options.aligned_fastx) {
					pretty_from_aux_inplace(pa->mate_pair);
					pretty_print_sam_fastx(pa->mate_pair, true);
				} else {
					remove_offending_fields(pa->mate_pair);
					pretty_print_sam_update(pa->mate_pair, true);
				}
				//revert_sam_string(pa->mate_pair);
			}
			pa=pa->next;
		}
	} else {
		m_ll->length=0; m_ll->head=NULL;
	}	
	return;	
	
		
}

static inline void pp_ll_append_and_check(pp_ll* ll,pretty * pa) {
	if (pa->paired_sequencing) {
		if (pa->proper_pair) {
			assert(pa->mapped && pa->mp_mapped);
			//both sides mapped
			pp_ll_append(ll+PAIRED,pa);
		} else if ((options.half_paired || options.sam_unaligned || options.unaligned_fastx)  && (pa->mapped || pa->mp_mapped)) {
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
		} else if ((options.sam_unaligned || options.unaligned_fastx) && !pa->mapped && !pa->mp_mapped) {
			//unmapped paired
			pp_ll_append(ll+UNMAPPED,pa);
		}
	} else {
		if (pa->mapped) {
			//unpaired and mapped
			pp_ll_append(ll+UNPAIRED,pa);
		} else if (options.sam_unaligned || options.unaligned_fastx) {
			//unpaired and unmapped
			pp_ll_append(ll+UNMAPPED,pa);
		}
	}
}
void grow_sam_pretty(sam_reader * sr) {
	fprintf(stderr,"Cannot double size!\n");
	exit(1);
	/*size_t new_pretty_stack_size=(size_t)(sr->pretty_stack_size*2);
	assert(new_pretty_stack_size>sr->pretty_stack_size);
	fprintf(stderr,"Growing %lu to %lu entries\n",sr->pretty_stack_size,new_pretty_stack_size);
	size_t old_size = sr->pretty_stack_size;
	char * old = (char*)sr->pretty_stack;
	sr->pretty_stack = (pretty*)realloc(sr->pretty_stack,new_pretty_stack_size*sizeof(pretty));
	if (sr->pretty_stack==NULL) {
		fprintf(stderr,"Failed realloc of pretty stack!\n");
		exit(1);
	}
	//memset(sr->pretty_stack + old_size , 0, (new_pretty_stack_size-old_size)*sizeof(pretty));
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
	}*/
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
	sr->pretty_stack_start=0;
	sr->pretty_stack_end=0;
	fprintf(stderr,"Starting a alignments stack with size %lu\n",options.alignments_stack_size);
	sr->pretty_stack_size=options.alignments_stack_size;
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
	sr->inter_offsets=(size_t*)malloc(sizeof(size_t)*options.read_rate);
	sr->pretty_stack_ends=(size_t*)malloc(sizeof(size_t)*options.read_rate);
	fprintf(stderr,"Size of inter_offsets is %d\n",options.read_rate);
	if (sr->inter_offsets==NULL) {
		fprintf(stderr,"Failed to allocate memory for inter_offsets\n");
		exit(1);
	}
	sr->last_tested=0;
	return sr;
}
void parse_sam(sam_reader * sr,fastx_readnames * fxrn) {
	char * current_newline=NULL;
	while (sr->fb->unseen_end!=sr->fb->unseen_inter && sr->last_tested<options.read_rate+fxrn->reads_seen) {
		const size_t pretty_stack_start_mod=sr->pretty_stack_start%sr->pretty_stack_size;
		const size_t pretty_stack_end_mod=sr->pretty_stack_end%sr->pretty_stack_size;
		if (sr->pretty_stack_start>sr->pretty_stack_end && pretty_stack_end_mod==pretty_stack_start_mod) {
			fprintf(stderr,"Stack out of space, need to flush!\n");
			return;
		}
        	size_t unseen_end_mod=sr->fb->unseen_end%sr->fb->size;
        	size_t unseen_inter_mod=sr->fb->unseen_inter%sr->fb->size;
		const size_t space=(unseen_inter_mod>=unseen_end_mod ? sr->fb->size : unseen_end_mod) - unseen_inter_mod;
		current_newline=(char*)memchr(sr->fb->base+unseen_inter_mod,'\n',space);
		if (current_newline==NULL) {
			sr->fb->unseen_inter+=space;
		} else {
			size_t current_newline_index=current_newline-sr->fb->base;
			if (current_newline_index==0 || sr->fb->base[current_newline_index-1]=='\0') {
				*current_newline='\0';
				current_newline_index++;
				for(; current_newline_index<sr->fb->size && sr->fb->base[current_newline_index]=='\0'; current_newline_index++) {sr->fb->base[current_newline_index]='\0';};
				sr->fb->unseen_inter+=current_newline_index-unseen_inter_mod;
				//sr->fb->unseen_start+=current_newline_index-unseen_start_mod;
				unseen_inter_mod=sr->fb->unseen_inter%sr->fb->size;
			} else {
				//fprintf(stderr,"Using %lu\n",unseen_start_mod);
				char * const line=sr->fb->base+unseen_inter_mod;
				//assert(sr->pretty_stack_size>sr->pretty_stack_filled);
				if (line[0]!='@') {
					const size_t length_of_string = current_newline-line;
					//if (sr->pretty_stack_size==sr->pretty_stack_filled) {	
					//	grow_sam_pretty(sr);
					//}
					for (; sr->last_tested<options.read_rate+fxrn->reads_seen; sr->last_tested++) {
						if (sr->last_tested==fxrn->reads_filled) {
							return;
						}
						//clear the row before using it!!!!!!!!!!
						const size_t read_id = sr->last_tested%options.read_rate;
						//ERRRROROROROROOROROROROROOROROROROORORORRRRR@!!!!
						char * const first_tab = (char*)memchr(line,'\t',length_of_string);
						if (first_tab==NULL) {
							fprintf(stderr,"CANNOT FIND FIRST TAB!\n");
							exit(1);
						}
						const size_t compare_length = first_tab-line;
						assert(compare_length!=0);
						char * const hit_list_read_name=fxrn->read_names+(sr->last_tested%fxrn->reads_inmem)*SIZE_READ_NAME*sizeof(char);
						//char buffer[SIZE_READ_NAME];
						//strncpy(buffer,hit_list_read_name,compare_length);
						//buffer[compare_length]='\0';
						//fprintf(stderr,"%p HL: %s,",sr,buffer);
						//strncpy(buffer,line,compare_length+15);
						//buffer[compare_length+15]='\0';
						//fprintf(stderr,"L: %s, %lu\n",buffer,read_id);
						//fprintf(stderr,"id: %d name: ||%s|| %lu seen %lu last\n",read_id,buffer,fxrn->reads_seen,sr->last_tested);
						//fprintf(stderr,"%s vs %s\n",buffer,hit_list_read_name);	
						if (strncmp(line,hit_list_read_name,compare_length)==0 ) {//&& !isdigit(hit_list_read_name[compare_length])) {	
							//fprintf(stderr,"XXXX %lu END, %lu read id\n",sr->pretty_stack_end,read_id);
							if (sr->pretty_stack_end-sr->pretty_stack_start>=sr->pretty_stack_size) {
								if (sr->last_tested==fxrn->reads_seen) {
									fprintf(stderr,"Alignments stack size is too small! Please use a larger size then '%lu' using '--alignments-stack-size'\n",options.alignments_stack_size);
									exit(1);
								}
								return;	
							}
							const size_t pa_index=(sr->pretty_stack_end++%sr->pretty_stack_size);
							//fprintf(stderr,"Put into index %lu, but read_id %lu\n",pa_index,read_id);
							pretty * const pa=sr->pretty_stack+pa_index;
							//pretty from string inplace memsets entry to 0!!! , assign after otherwise gets wiped!
							pretty_from_string_inplace(line,length_of_string,pa);
							pa->read_id=read_id;
							pa->sam_header=false;
							pa->fileno=sr->fileno;	
							assert(pa->read_name!=NULL);
							assert(pa->mate_pair==NULL);
							if (pa->paired_sequencing) {
								options.paired=true;
								pretty * const previous=sr->pretty_stack + (pa_index+sr->pretty_stack_size-1)%sr->pretty_stack_size;
								//fprintf(stderr,"previous %p\n",previous);
								if (sr->pretty_stack_start+1!=sr->pretty_stack_end && !previous->sam_header && previous->mate_pair==NULL) {
									previous->mate_pair=pa; pa->mate_pair=previous;
									assert(pa->read_name!=NULL);
									assert(pa->mate_pair->read_name!=NULL);
									pp_ll_append_and_check(sr->pp_lls+LL_ALL*pa->read_id,pa);
									assert(pa->read_name!=NULL);
									assert(pa->mate_pair->read_name!=NULL);
								} else {
									pa->mate_pair=NULL;
								}
							} else {
								options.unpaired=true;
								pp_ll_append_and_check(sr->pp_lls+LL_ALL*pa->read_id,pa);
							}
							*current_newline='\0';
							break;
						} else {	
							sr->inter_offsets[read_id]=sr->fb->unseen_start;
							//fprintf(stderr,"Setting interoffset %lu\n",sr->fb->unseen_start);
							sr->pretty_stack_ends[read_id]=sr->pretty_stack_end;
							sr->fb->unseen_start=sr->fb->unseen_inter;
							//fprintf(stderr,"Asign start %lu, %lu\n",read_id,sr->fb->unseen_start);
							//if ((sr->last_tested+1)<(options.read_rate+fxrn->reads_seen)) {
							//	const size_t next_read_id=(read_id+1)%options.read_rate;
							//	memset(sr->pp_lls+LL_ALL*next_read_id,0,sizeof(pp_ll)*LL_ALL);
							//	fprintf(stderr,"Clearing spot %lu\n",next_read_id);
							//}
						}
					}
					
				} else {
					found_sam_headers=true;
					pretty * pa=sr->pretty_stack+sr->pretty_stack_end++%sr->pretty_stack_size;
					pa->sam_string=line;
					*current_newline='\0';
					pp_ll_append(sr->sam_headers,pa);
					pa->sam_header=true;
				}
				if (sr->last_tested!=options.read_rate+fxrn->reads_seen) {
					sr->fb->unseen_inter+=current_newline_index-unseen_inter_mod+1;
					unseen_inter_mod=sr->fb->unseen_inter%sr->fb->size;
				}
			}
			sr->fb->exhausted=false;
		}
	}
	//fprintf(stderr,"RETURN GRACE %lu!=%lu,%lu<%d+%lu \n",sr->fb->unseen_end,sr->fb->unseen_inter, sr->last_tested,options.read_rate,fxrn->reads_seen);
	if (sr->fb->frb.eof==1 && sr->fb->unseen_end==sr->fb->unseen_inter && sr->last_tested==fxrn->reads_seen) {
		//sr->last_tested++;
		fprintf(stderr,"EOF\n");
	}
		//fprintf(stderr,"RETURN %llu - %llu < %llu \n",sr->last_tested, fxrn->reads_seen, read_amount);
	/*if (last_tested<read_amount && (sr->fb->frb.eof==0 || sr->fb->unseen_end!=sr->fb->unseen_start)) {
		fprintf(stderr,"failed to read into buffer, please increase buffer size!\n");
		exit(1);
	}*/
}
