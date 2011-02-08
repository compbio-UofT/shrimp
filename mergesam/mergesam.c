#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>
#include <assert.h>
#include "sam2pretty_lib.h"
#include "file_buffer.h"
#include "mergesam_heap.h"
#include <ctype.h>
#include <omp.h>

#define MIN(a,b) (((a) < (b)) ? (a) : (b))
#define MAX_INT32	2147483647

#define PAIRED		0
#define UNPAIRED	1
#define FIRST_LEG	2
#define SECOND_LEG	3
#define UNMAPPED	4
#define LL_ALL		5

//IO SETTINGS
#define DEF_READ_SIZE	1024*1024*5
#define DEF_BUFFER_SIZE	1024*1024*15
#define DEF_READ_RATE 10
size_t read_rate=DEF_READ_RATE;
size_t buffer_size=DEF_BUFFER_SIZE; 
size_t read_size=DEF_READ_SIZE;

//Variables for IO of read names
#define SIZE_READ_NAME 255
char * reads_filename=NULL; 
FILE * reads_file=NULL;
size_t reads_offset=0;
size_t reads_seen=0;
size_t reads_unseen=0;
size_t reads_inmem=0;
size_t reads_filled=0;
size_t reads_line_number=0;
bool reads_exhausted=false;
char * read_names;
bool fastq_seen_name=false;
bool fastq_seen_plus=false;
int64_t fastq_seen_seq=0;
int64_t fastq_seen_qual=0;

//Pretty linked list
typedef struct pp_ll pp_ll;
struct pp_ll {
	pretty * head;
	pretty * tail;
	size_t length;
};
pp_ll * sam_headers;
pp_ll ** pp_ll_index;
pp_ll * master_ll;

//SAM file ios
#define GROWTH_FACTOR 1.3
heap_pa * thread_heaps;
int number_of_sam_files=0;
typedef struct sam_info sam_info;
struct sam_info {
	file_buffer * fb;
	pretty * pretty_stack;
	size_t pretty_stack_size;
	size_t pretty_stack_filled;	
	pp_ll * pp_lls;
	pp_ll * sam_headers;
};
sam_info ** sam_files;
char * sam_header_filename=NULL;

//Filtering parameters
#define DEF_MAX_ALIGNMENTS	-1
#define DEF_MAX_OUTPUTS	20
#define DEF_INSERT_SIZE -1
int32_t max_alignments=DEF_MAX_ALIGNMENTS;
int32_t max_outputs=DEF_MAX_OUTPUTS;
int expected_insert_size=DEF_INSERT_SIZE;

bool fastq=false;
bool strata=false;
bool half_paired=false;
bool sam_unaligned=false;
bool found_sam_headers=false;
bool sam_format=false;

bool paired=false;
bool unpaired=false; 
bool colour_space=false;
bool letter_space=false;


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

static inline void revert_sam_string(pretty * pa) {
	char * nill = (char*)memchr(pa->sam_string,'\0',pa->sam_string_length);
	while (nill!=NULL) {
		*nill='\t';
		nill=(char*)memchr(pa->sam_string,'\0',pa->sam_string_length-(nill-pa->sam_string));
	}
}

static inline void pp_ll_combine_and_check_single(pp_ll * m_ll,pp_ll ** ll,int offset,pretty ** unaligned_pa) {
	heap_pa * h=thread_heaps+omp_get_thread_num();
	h->load=0;
	heap_pa_elem e;
	pretty * pa=NULL;
	int i;
	for (i=0; i<number_of_sam_files; i++) {
		pp_ll * local_ll = ll[i]+offset;
		if (local_ll->length>0) {
			pa = local_ll->head;
			while (pa!=NULL) {
				e.score=pa->score+(pa->mate_pair!=NULL ? pa->mate_pair->score : 0);
				e.isize_score=MAX_INT32;
				e.rest=pa;
				if (strata) {
					heap_pa_insert_bounded_strata(h,&e);
				} else {
					heap_pa_insert_bounded(h,&e);
				}
				pa=pa->next;
			}
		}
	}
	if (h->load==0 || (max_alignments>0 && h->load>max_alignments)) {
		//m_ll->length=0;
	} else {
		assert(max_outputs>0);
		i=h->load>max_outputs ? 1 : 0;
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
		m_ll->length+=h->load- (h->load>max_outputs ? 1 : 0);
	}
	if (sam_unaligned && h->load>0) {
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

static inline void pp_ll_combine_and_check(pp_ll * m_ll,pp_ll ** ll) {
	//want to make a heap
	//check first if paired of not
	assert(m_ll->length==0);
	heap_pa * h=thread_heaps+omp_get_thread_num();
	h->load=0;
	heap_pa_elem e;
	pretty * pa = NULL;
	pretty * unaligned_pa=NULL;
	if (paired) {
		assert(!unpaired);
		//check paired first, and then go lower
		int i;
		for (i=0; i<number_of_sam_files; i++) {
			pp_ll * paired_ll = ll[i]+PAIRED;
			if (paired_ll->length>0) {
				pa = paired_ll->head;
				while (pa!=NULL) {
					e.score=pa->score+pa->mate_pair->score;
					assert(pa->mapped && pa->mate_pair!=NULL);
					assert(pa->mate_pair->mapped &&  pa->mp_mapped);
					e.isize_score=expected_insert_size >= 0 ? abs(expected_insert_size-pa->isize) : 0;
					e.rest=pa;
					if (strata) {
						heap_pa_insert_bounded_strata(h,&e);
					} else {
						heap_pa_insert_bounded(h,&e);
					}
					pa=pa->next;
				}
			}
		}
		if (h->load==0 || (max_alignments>0 && h->load>max_alignments)) {
			//m_ll->length=0;
		} else {
			assert(max_outputs>0);
			i=h->load>max_outputs ? 1 : 0;
			m_ll->head=h->array[i].rest;
			m_ll->tail=h->array[h->load-1].rest;
			for (; i<h->load-1; i++) {
				h->array[i].rest->next=h->array[i+1].rest;
			}
			h->array[i].rest->next=NULL;
			m_ll->length=h->load- (h->load>max_outputs ? 1 : 0);
		}
		if (sam_unaligned && h->load>0) {
			unaligned_pa=h->array[0].rest;
		}
	}
	if (h->load==0) {
		if (paired && half_paired) {
			pp_ll_combine_and_check_single(m_ll,ll,FIRST_LEG,&unaligned_pa);
			pp_ll_combine_and_check_single(m_ll,ll,SECOND_LEG,&unaligned_pa);
		} else if (unpaired) {
			pp_ll_combine_and_check_single(m_ll,ll,UNPAIRED,&unaligned_pa);
		}
	}
	if (m_ll->length==0 && sam_unaligned) {
		if (unaligned_pa==NULL) {
			pp_ll_combine_and_check_single(m_ll,ll,UNMAPPED,&unaligned_pa);
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
		} else if (half_paired && (pa->mapped || pa->mp_mapped)) {
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
		} else if (sam_unaligned && !pa->mapped && !pa->mp_mapped) {
			//unmapped paired
			pp_ll_append(ll+UNMAPPED,pa);
		}
	} else {
		if (pa->mapped) {
			//unpaired and mapped
			pp_ll_append(ll+UNPAIRED,pa);
		} else if (sam_unaligned) {
			//unpaired and unmapped
			pp_ll_append(ll+UNMAPPED,pa);
		}
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


void sam_close(sam_info * si) {
	fb_close(si->fb);
	free(si->pretty_stack);
	free(si->pp_lls);
	free(si);
}

sam_info * sam_open(char * sam_filename) {
	sam_info * si = (sam_info*)malloc(sizeof(sam_info));
	if (si==NULL) {
		fprintf(stderr,"thread_open_sam : failed to allocate memory for thread info structure\n");
		exit(1);
	}
	si->fb=fb_open(sam_filename,buffer_size,read_size);
	si->pretty_stack_filled=0;
	si->pretty_stack_size=reads_inmem*10;
	si->pretty_stack=(pretty*)malloc(sizeof(pretty)*si->pretty_stack_size);
	if (si->pretty_stack==NULL) {
		fprintf(stderr,"Failed to allocate memory for pretty_stack\n");
		exit(1);
	}
	memset(si->pretty_stack,0,sizeof(pretty)*si->pretty_stack_size);
	si->pp_lls=(pp_ll*)malloc(sizeof(pp_ll)*LL_ALL*read_rate);
	if (si->pp_lls==NULL) {
		fprintf(stderr,"Failed to allocate memory for pp_lls.\n");
	}
	memset(si->pp_lls,0,sizeof(pp_ll)*LL_ALL*read_rate);
	return si;
}

void usage(char * s) {
	fprintf(stderr, 
	"usage: %s [options/parameters] [--colour-space/--letter-space] <r> <s1> <s2> ...\n", s);
	fprintf(stderr,
	"   <r>     Reads filename, if paired then one of the two paired files\n");
	fprintf(stderr,
	"   <s?>    A SAM file input for mergigng\n");
	fprintf(stderr,
	"Parameters:      (all sizes are in bytes unless specified)\n");
	fprintf(stderr,
	"      --buffer-size    File buffer size in memory per file   (Default: %d)\n",DEF_BUFFER_SIZE);
	fprintf(stderr,
	"      --read-size      Read size, read into buffer with this (Default: %d)\n",DEF_READ_SIZE);
	fprintf(stderr,
	"      --read-rate      How many reads to process at once     (Default: %d)\n",DEF_READ_RATE);
	fprintf(stderr,
	"      --expected-isize Expected insert size, for tie-break   (Default: disabled)\n");
	fprintf(stderr,
	"   -N/--threads        The number of threads to use          (Default: 1)\n");
	fprintf(stderr,
	"   -o/--report         The maximum alignments to report      (Default: %d)\n",DEF_MAX_OUTPUTS);
	fprintf(stderr,
	"      --max-alignments Max. align. per read  (-1=all)        (Default: %d)\n",DEF_MAX_ALIGNMENTS); 
	fprintf(stderr,
	"      --sam-header     Use file as SAM header\n");
	fprintf(stderr,"\nOptions:\n");
	fprintf(stderr,
	"      --help           This usage screen\n");
	fprintf(stderr,
	"   -Q/--fastq          Reads are in fastq format             (Default: disabled)\n");
	fprintf(stderr,
	"   -E/--sam            Output in SAM format                  (Default: disabled)\n");
	fprintf(stderr,
	"      --sam-unaligned  Unaligned reads in SAM output         (Default: disabled)\n");
	fprintf(stderr,
	"      --half-paired    Output half mapped read pairs         (Default: disabled)\n");
	fprintf(stderr,
	"      --strata         Print only the best scoring hits\n");
	fprintf(stderr,
	"      --colour-space   Reads file contains ABSOLiD data\n");
	fprintf(stderr,
	"      --letter-space   Reads file contains non-ABSOLiD data\n");	
	exit(1);
}


struct option long_op[] =
        {
		{"fastq", 0, 0, 'Q'},
		{"sam",1,0,'E'},
		{"threads",1,0,'N'},
		{"report",1,0,'o'},
		{"max-alignments",1,0,200},
		{"sam-header",1,0,201},
		{"colour-space",0,0,202},
		{"letter-space",0,0,203},
                {"help", 0, 0, 204},
		{"buffer-size", 1, 0, 205},
		{"read-size", 1, 0, 206},
		{"read-rate",1,0,207},
		{"expected-isize",1,0,208},
		{"half-paired",0,0,209},
		{"sam-unaligned",0,0,210},
		{"strata",0,0,211},
                {0,0,0,0}
        };

static inline void fill_fb(file_buffer * fb) {
	while (!fb->exhausted) {
		fill_read_buffer(&fb->frb);
		add_read_buffer_to_main(fb);
		if (!fb->exhausted && !fb->changed && fb->frb.eof==0) {
			fprintf(stderr,"too small buffer!\n");
			exit(1);
		}
	}
	//fprintf(stdout,"Filled %lu to %lu of %lu |%s|\n",fb->unseen_start, fb->unseen_end, fb->size,fb->base);
}


static inline void reconsolidate_reads(void) {
	assert(reads_seen+reads_unseen==reads_filled);
	if (reads_unseen!=reads_filled) {
		//fprintf(stderr,"reconsolidate\n");
		memmove(read_names,read_names+reads_seen*SIZE_READ_NAME,(reads_filled-reads_seen)*SIZE_READ_NAME);
		//reads_offset+=reads_seen;
		reads_filled=reads_unseen;
		reads_seen=0;
		reads_exhausted=false;
	}
}

void parse_reads(file_buffer * fb) {
	assert(reads_seen+reads_unseen==reads_filled);
	reconsolidate_reads();
	char * current_newline=NULL;
	while (fb->unseen_end!=fb->unseen_start && reads_filled<reads_inmem) {
        	size_t unseen_end_mod=fb->unseen_end%fb->size;
        	size_t unseen_start_mod=fb->unseen_start%fb->size;
		size_t space=(unseen_start_mod>=unseen_end_mod ? fb->size : unseen_end_mod) - unseen_start_mod;
		current_newline=(char*)memchr(fb->base+unseen_start_mod,'\n',space);
		if (current_newline==NULL) {
			fb->unseen_start+=space;
		} else {
			*current_newline='\0';
			size_t current_newline_index=current_newline-fb->base;
			if (current_newline_index==0 || fb->base[current_newline_index-1]=='\0') {
				current_newline_index++;
				for(; current_newline_index<fb->size && fb->base[current_newline_index]=='\0'; current_newline_index++) {fb->base[current_newline_index]='\0';};
				fb->unseen_start+=current_newline_index-unseen_start_mod;
				unseen_start_mod=fb->unseen_start%fb->size;
			} else {
				char * read_name=fb->base+unseen_start_mod;
				if (fastq) {
					if (!fastq_seen_name) {
						if (read_name[0]=='@') {
							size_t read_name_length=strlen(read_name);
							memmove(read_names+reads_filled*SIZE_READ_NAME,read_name+1,read_name_length);
							read_names[reads_filled*SIZE_READ_NAME+read_name_length]='\0';
							char * space=strchr(read_names+reads_filled*SIZE_READ_NAME,' ');
							if (space!=NULL) {*space='\0';};
							char * tab=strchr(read_names+reads_filled*SIZE_READ_NAME,'\t');
							if (tab!=NULL) {*tab='\0';};
							reads_unseen++;
							reads_filled++;
							fastq_seen_name=true;
							fastq_seen_seq=colour_space ? -1 : 0;
							fastq_seen_qual=0;
							fastq_seen_plus=false;
						}
					} else {
						if (!fastq_seen_plus) {
							if (read_name[0]=='+') {
								fastq_seen_plus=true;
							} else {
								fastq_seen_seq+=strlen(read_name)-1;
							}
						} else {
							fastq_seen_qual+=strlen(read_name)-1;
							assert(fastq_seen_qual<=fastq_seen_seq);
							if (fastq_seen_qual==fastq_seen_seq) {
								fastq_seen_name=false;
							}	
						}
					}
				} else {
					if (read_name[0]=='>') {
						size_t read_name_length=strlen(read_name);
						memmove(read_names+reads_filled*SIZE_READ_NAME,read_name+1,read_name_length);
						read_names[reads_filled*SIZE_READ_NAME+read_name_length]='\0';
						char * space=strchr(read_names+reads_filled*SIZE_READ_NAME,' ');
						if (space!=NULL) {*space='\0';};
						char * tab=strchr(read_names+reads_filled*SIZE_READ_NAME,'\t');
						if (tab!=NULL) {*tab='\0';};
						reads_unseen++;
						reads_filled++;
					}
				}
				fb->unseen_start+=current_newline_index-unseen_start_mod+1;
				unseen_start_mod=fb->unseen_start%fb->size;
			}
			fb->exhausted=false;
		}
	}
	if (fb->exhausted || reads_filled==reads_inmem) {	
		reads_exhausted=true;
	} else {
		reads_exhausted=false;
	}
}


static inline void grow_sam_pretty(sam_info * si) {
	size_t new_pretty_stack_size=(size_t)(si->pretty_stack_size*GROWTH_FACTOR+1);
	assert(new_pretty_stack_size>si->pretty_stack_size);
	fprintf(stderr,"Growing %lu to %lu entries\n",si->pretty_stack_size,new_pretty_stack_size);
	size_t old_size = si->pretty_stack_size;
	char * old = (char*)si->pretty_stack;
	si->pretty_stack = (pretty*)realloc(si->pretty_stack,new_pretty_stack_size*sizeof(pretty));
	if (si->pretty_stack==NULL) {
		fprintf(stderr,"Failed realloc of pretty stack!\n");
		exit(1);
	}
	memset(si->pretty_stack + old_size , 0, (new_pretty_stack_size-old_size)*sizeof(pretty));
	si->pretty_stack_size=new_pretty_stack_size;	
	char * new_p = (char*)si->pretty_stack;
	size_t i;
	for (i=0; i<LL_ALL*read_rate; i++) {
		pp_ll * ll= si->pp_lls+i;
		if (ll->head!=NULL) {
			ll->head=(pretty*)(((char*)ll->head)+(new_p-old));
		}
		if (ll->tail!=NULL) {
			ll->tail=(pretty*)(((char*)ll->tail)+(new_p-old));
		}
	}
	for (i=0; i<old_size; i++) {
		pretty * pa=si->pretty_stack+i;
		if (pa->next!=NULL) {
			pa->next=(pretty*)(((char*)pa->next)+(new_p-old));
		}
		if (pa->mate_pair!=NULL) {
			pa->mate_pair=(pretty*)(((char*)pa->mate_pair)+(new_p-old));
		}
	}
}

static void print_buffer(file_buffer * fb) {
	size_t start = fb->unseen_start;
	size_t end = fb->unseen_end;
	size_t index;
	fprintf(stderr,"|FB|");
	fprintf(stderr,"%lu %lu\n",fb->unseen_start, fb->unseen_end);
	for (index=start; index<end; index++) {
		fprintf(stderr,"%c",fb->base[index%fb->size]);
	}
	fprintf(stderr,"|END FB|\n");
}

static void print_frb_buffer(file_buffer * fb) {
	size_t start=fb->frb.seen;
	size_t index;
	fprintf(stderr,"|FRB|\n");
	for (index=start; index<fb->frb.size; index++) {
		fprintf(stderr,"%c",fb->frb.base[index]);
	}
	fprintf(stderr,"|END FRB|\n");
}

static void parse_sam(sam_info * si,size_t read_amount) {
	char * current_newline=NULL;
	size_t last_tested=0;
	//fprintf(stderr,"|FRB||%s|||\n",si->fb->frb.base);
	//fprintf(stderr,"|FB||%s|||\n",si->fb->base);
	//print_buffer(si->fb);
	//print_frb_buffer(si->fb);
	while (si->fb->unseen_end!=si->fb->unseen_start && last_tested<read_amount) {
        	size_t unseen_end_mod=si->fb->unseen_end%si->fb->size;
        	size_t unseen_start_mod=si->fb->unseen_start%si->fb->size;
		size_t space=(unseen_start_mod>=unseen_end_mod ? si->fb->size : unseen_end_mod) - unseen_start_mod;
		current_newline=(char*)memchr(si->fb->base+unseen_start_mod,'\n',space);
		if (current_newline==NULL) {
			si->fb->unseen_start+=space;
		} else {
			size_t current_newline_index=current_newline-si->fb->base;
			if (current_newline_index==0 || si->fb->base[current_newline_index-1]=='\0') {
				*current_newline='\0';
				current_newline_index++;
				for(; current_newline_index<si->fb->size && si->fb->base[current_newline_index]=='\0'; current_newline_index++) {si->fb->base[current_newline_index]='\0';};
				si->fb->unseen_start+=current_newline_index-unseen_start_mod;
				unseen_start_mod=si->fb->unseen_start%si->fb->size;
			} else {
				//fprintf(stderr,"Using %lu\n",unseen_start_mod);
				char * line=si->fb->base+unseen_start_mod;
				assert(si->pretty_stack_size>si->pretty_stack_filled);
				if (line[0]!='@') {
					size_t length_of_string = current_newline-line;
					if (si->pretty_stack_size==si->pretty_stack_filled) {	
						grow_sam_pretty(si);
					}
					for (; last_tested<read_amount; last_tested++) {
						char * first_tab = (char*)memchr(line,'\t',length_of_string);
						if (first_tab==NULL) {
							fprintf(stderr,"CANNOT FIND FIRST TAB!\n");
							exit(1);
						}
						size_t compare_length = first_tab-line;
						char * hit_list_read_name=read_names+(last_tested+reads_offset)*SIZE_READ_NAME*sizeof(char);
						//char buffer[compare_length+1];
						//strncpy(buffer,line,compare_length);
						//buffer[compare_length]='\0';
						//	*current_newline='\0';
						//fprintf(stdout,"B: |%s| vs |%s|\n",hit_list_read_name,line);
							//*current_newline='\n';
						
						if (strncmp(line,hit_list_read_name,compare_length)==0 ) {//&& !isdigit(hit_list_read_name[compare_length])) {	
							pretty * pa=si->pretty_stack+si->pretty_stack_filled++;
							pa->read_id=last_tested;
							pa->sam_header=false;
							pretty_from_string_inplace(line,length_of_string,pa);
							//fprintf(stderr,"Found a hit for read |%s| vs |%s|\n",si->pretty_stack[si->pretty_stack_filled-1].read_name,read_names+(last_tested+reads_offset)*SIZE_READ_NAME*sizeof(char));
							//check whether paired, unpaired or first or second leg
							//fprintf(stdout,"FOUND score is %d\n",pa->score);
							if (pa->paired_sequencing) {
								paired=true;
								if (si->pretty_stack_filled!=1 && !(pa-1)->sam_header && (pa-1)->mate_pair==NULL) {
									(pa-1)->mate_pair=pa; pa->mate_pair=pa-1;
									//fprintf(stderr,"Added as pair!\n");
									pp_ll_append_and_check(si->pp_lls+LL_ALL*pa->read_id,pa);
								} else {
									pa->mate_pair=NULL;
								}
							} else {
								unpaired=true;
								pp_ll_append_and_check(si->pp_lls+LL_ALL*pa->read_id,pa);
							}
							*current_newline='\0';
							break;
						} 	
					}
					
				} else {
					found_sam_headers=true;
					pretty * pa=si->pretty_stack+si->pretty_stack_filled++;
					pa->sam_string=line;
					*current_newline='\0';
					pp_ll_append(si->sam_headers,pa);
					pa->sam_header=true;
				}
				if (last_tested!=read_amount) {
					si->fb->unseen_start+=current_newline_index-unseen_start_mod+1;
					unseen_start_mod=si->fb->unseen_start%si->fb->size;
				}
			}
			si->fb->exhausted=false;
		}
	}
	if (last_tested<read_amount && (si->fb->frb.eof==0 || si->fb->unseen_end!=si->fb->unseen_start)) {
		fprintf(stderr,"failed to read into buffer, please increase buffer size!\n");
		exit(1);
	}
	//fprintf(stdout,"wanted to read in %d but last tested is %d\n",read_amount,last_tested);
	/*if (si->fb->exhausted || reads_filled==reads_inmem) {	
		reads_exhausted=true;
	} else {
		reads_exhausted=false;
	}*/
}


static size_t inline string_to_byte_size(char * s) {
	char * x=s;
	while (isdigit(x[0])) {x++;};
	size_t multiplier=1;
	if (*x=='K') {
		multiplier=1024;
	} else if (*x=='M') {
		multiplier=1024*1024;
	} else if (*x=='G') {
		multiplier=1024*1024*1024;
	}
	char old_x=*x;
	*x='\0';
	int ret=atoi(s);
	*x=old_x;
	if (ret<=0) {
		return 0;
	}
	return ret*multiplier;
}

int main (int argc, char ** argv) {
        int op_id;
        char short_op[] = "o:QN:E";
	buffer_size=DEF_BUFFER_SIZE;
	read_size=DEF_READ_SIZE;
	read_rate=DEF_READ_RATE;
        char c = getopt_long(argc, argv, short_op, long_op, &op_id);
	int threads=1;
        while (c != EOF) {
		switch (c) {
		case 'Q':
			fastq=true;
			break;
		case 205:
			buffer_size=string_to_byte_size(optarg);
			break;
		case 206:
			read_size=string_to_byte_size(optarg);
			break;
		case 207:
			read_rate=atol(optarg);
			break;
		case 208:
			expected_insert_size=atoi(optarg);
			if (expected_insert_size<0) {
				fprintf(stderr,"Please specify a insert size >= 0!\n");
				usage(argv[0]);
			}
			break;
		case 'N':
			threads=atoi(optarg);
			break;
		case 209:
			half_paired=true;
			break;
		case 210:
			sam_unaligned=true;
			break;
		case 211:
			strata=true;
			break;
		case 'o':
			max_outputs=atoi(optarg);
			if (max_outputs<=0) {
				fprintf(stderr,"Please specify a max_output that is positive!\n");
				usage(argv[0]);
			}
			break;
		case 200:
			max_alignments=atoi(optarg);
			if (max_alignments<=0) {
				fprintf(stderr,"Please specify a max_alignments that is positive!\n");
				usage(argv[0]);
			}
			break;
		case 204:
			usage(argv[0]);
			break;
		case 201:
			{
			sam_header_filename=optarg;
			FILE * sam_header_file = fopen(sam_header_filename,"r");
			if (sam_header_file==NULL) {
				perror("Failed to open sam header file ");
				usage(argv[0]);
			}
			size_t buffer_size=2046;
			char buffer[buffer_size];
			size_t read; bool ends_in_newline=true;
			while ((read=fread(buffer,1,buffer_size-1,sam_header_file))) {
				buffer[read]='\0';
				fprintf(stdout,"%s",buffer);
				if (buffer[read-1]=='\n') {
					ends_in_newline=true;
				} else {
					ends_in_newline=false;
				}
			}
			if (!ends_in_newline) {
				fprintf(stdout,"\n");
			}
			}
		case 'E':
			sam_format=true;
			break;
		case 202:
			colour_space=true;
			break;
		case 203:
			letter_space=true;
			break;
		default:
			usage(argv[0]);
			break;
		}
        	c = getopt_long(argc, argv, short_op, long_op, &op_id);
	}


	if ((colour_space && letter_space) || (!colour_space && !letter_space)) {
		fprintf(stderr,"can only enter either --colour-space or --letter-space not both!\n");
		usage(argv[0]);
	}
	if (!sam_format) {
		fprintf(stderr,"%s currently only supports sam format, please use '--sam' or '-E'\n",argv[0]);
		usage(argv[0]);
	}

	omp_set_num_threads(threads); 
	fprintf(stderr,"Set to %d threads!\n",threads);
		
	
	if (argc<=optind+1) {
		fprintf(stderr,"Please specify reads file and at least one sam file!\n");
		usage(argv[0]);
	}

	fprintf(stderr,"Size of %lu\n",sizeof(pretty));
	reads_inmem=20*read_rate;
	read_names=(char*)malloc(sizeof(char)*reads_inmem*SIZE_READ_NAME);
	if (read_names==NULL) {
		fprintf(stderr,"failed to allocate memory for read_names\n");
		exit(1);
	}


	argc-=optind;
	argv+=optind;
	
	reads_filename=argv[0];
	fprintf(stderr,"Using %s as reads filename\n",reads_filename);
	argc--;
	argv++;

	number_of_sam_files=argc;
	//Open each sam input file	
	sam_files=(sam_info**)malloc(sizeof(sam_info*)*number_of_sam_files);
	if (sam_files==NULL) {
		fprintf(stderr,"failed to allocate memory for sam_files!\n");
		exit(1);
	}

	master_ll = (pp_ll*)malloc(sizeof(pp_ll)*read_rate);
	if (master_ll==NULL) {
		fprintf(stderr,"Failed to allocate memory for master_ll\n");
		exit(1);
	}


	//allocate memory for sam_headers
	sam_headers=(pp_ll*)malloc(sizeof(pp_ll)*number_of_sam_files);
	if (sam_headers==NULL) {
		fprintf(stderr,"Failed to allocate memory for sam_headers\n");
		exit(1);
	}
	

	//index first the read then the file number	
	pp_ll_index = (pp_ll**)malloc(sizeof(pp_ll*)*read_rate*number_of_sam_files);
	if (pp_ll_index==NULL) {
		fprintf(stderr,"Failed to allocate memory for pp_ll_index!\n");
		exit(1);
	}	
	int i;
	for (i=0; i<number_of_sam_files; i++) {
		sam_files[i]=sam_open(argv[i]);
		sam_files[i]->sam_headers=sam_headers+i;
		int j; 
		for (j=0; j<(int)read_rate; j++) {
			pp_ll_index[j*number_of_sam_files+i]=sam_files[i]->pp_lls+j*LL_ALL;
		}
	}

	//calculate alignments cutoff
	int32_t alignments_cutoff=max_alignments==-1 ? max_outputs : MIN(max_alignments,max_outputs);


	//max the heaps for each thread
	fprintf(stderr,"There are %d threads\n",threads);
	thread_heaps=(heap_pa* )malloc(sizeof(heap_pa)*threads);
	if (thread_heaps==NULL) {
		fprintf(stderr,"Failed to allocate memory for thread_heaps!\n");
		exit(1);
	}
	for (i=0; i<threads; i++ ) {
		heap_pa_init(thread_heaps+i,alignments_cutoff+1);
		fprintf(stderr,"INIT THREAD_HEAP %p\n",thread_heaps+i);
	}	

	size_t reads_processed=0;
	//get the hit list, process it, do it again!
	fprintf(stderr,"Setting up buffer with size %lu and read_size %lu\n",buffer_size,read_size);
	file_buffer * fb = fb_open(reads_filename,buffer_size,read_size);
	while (!reads_exhausted) {
		//Populate the hitlist to as large as possible
		while (!reads_exhausted) {
			fill_fb(fb);
			parse_reads(fb);	
		}
		//clear the sam_headers memory
		memset(sam_headers,0,sizeof(pp_ll)*number_of_sam_files);	
	
		//Hit list in memory start processing!
		int i;
		int reads_to_process=0;
		/*for (i=0; i<(int)reads_filled; i++) {
			fprintf(stdout,"FRN: %s\n",read_names+i*SIZE_READ_NAME);	
		}
		reads_unseen=0;
		reads_seen=reads_filled;
		reads_offset=0;
		reconsolidate_reads();*/
		
		//while (reads_offset+read_rate<reads_filled || (fb->frb.eof==1 && reads_offset<reads_filled)) {
		while (reads_to_process>=0 && reads_offset<reads_filled) {
			assert(reads_to_process<=(int)read_rate);
			if (paired && unpaired) {
				fprintf(stderr,"FAIL! can't have both paired and unpaired data in input file!\n");
				exit(1);
			}
			#pragma omp parallel for //schedule(static,10)
			for (i=0; i<reads_to_process; i++) {	
				pp_ll_combine_and_check(master_ll+i,pp_ll_index+i*number_of_sam_files);
			}
			//one thread - deal with sam headers before the reads themselves
			if (found_sam_headers && sam_header_filename==NULL) {
				int header_entries=0;
				int i;
				for (i=0; i<number_of_sam_files; i++) {
					header_entries+=sam_files[i]->sam_headers->length;
				}
				//get pointers to each header line
				char ** sam_lines=(char**)malloc(sizeof(char*)*header_entries);
				if (sam_lines==NULL) {
					fprintf(stderr,"Failed to allocate memory for sam_header entries!\n");
					exit(1);
				}	
				int index=0;
				for (i=0; i<number_of_sam_files; i++) {
					pretty * pa = sam_files[i]->sam_headers->head;
					while(pa!=NULL) {
						sam_lines[index++]=pa->sam_string;
						pa=pa->next;
					}
				}
				//want to sort the headers here
				qsort(sam_lines, header_entries, sizeof(char*),sam_header_sort);	
				//want to print the headers here
				assert(index>0);
				fprintf(stdout,"%s\n",sam_lines[0]);
				for (i=1; i<index; i++) {
					int ret=strcmp(sam_lines[i],sam_lines[i-1]);
					if (ret!=0) {
						fprintf(stdout,"%s\n",sam_lines[i]);
					}	
				}	
				free(sam_lines);
				memset(sam_headers,0,sizeof(pp_ll)*number_of_sam_files);	
				found_sam_headers=false;
			}
			
			//one thread
			for (i=0; i<reads_to_process; i++) {	
				pretty * pa=master_ll[i].head;
				while (pa!=NULL) {
					if (pa->mate_pair!=NULL) {
						fprintf(stdout,"%s\n",pa->mate_pair->sam_string);		
					}
					fprintf(stdout,"%s\n",pa->sam_string);		
					pa=pa->next;
				}
			}
			reads_processed+=reads_to_process;
			reads_offset+=reads_to_process;	
			memset(master_ll,0,sizeof(pp_ll)*read_rate);
			if (reads_offset+read_rate<reads_filled) {
				reads_to_process=read_rate;
			} else if (fb->frb.eof==1 && reads_offset<reads_filled) {
				reads_to_process=reads_filled-reads_offset;
			} else {
				break;
			}
			//parallel
			#pragma omp parallel for 
			for (i=0; i<number_of_sam_files; i++) {
				memset(sam_files[i]->pretty_stack,0,sizeof(pretty)*sam_files[i]->pretty_stack_filled);
				sam_files[i]->pretty_stack_filled=0;
				memset(sam_files[i]->pp_lls,0,sizeof(pp_ll)*LL_ALL*read_rate);
				//fprintf(stderr,"THREAD %d\n",omp_get_thread_num());
				fill_fb(sam_files[i]->fb);	
				parse_sam(sam_files[i],reads_to_process);
			}
		}
		//Should be done processing
		//keep unseen reads, and do it again
		reads_unseen=reads_filled-reads_offset;
		reads_seen=reads_filled-reads_unseen;
		reads_offset=0;
		reconsolidate_reads();
	}
	fprintf(stderr,"Processed %lu reads\n",reads_processed);
	free(master_ll);
	free(sam_headers);
	free(pp_ll_index);
	for (i=0; i<number_of_sam_files; i++) {
		sam_close(sam_files[i]);
	}
	free(sam_files);
	for (i=0; i<threads; i++) {
		heap_pa_destroy(thread_heaps+i);
	}
	free(thread_heaps);
	fb_close(fb);
	free(read_names);
	return 0;
}




