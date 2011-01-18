#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <omp.h>
#include <math.h>
#include <ctype.h>
#include <zlib.h>
#include <limits.h>
#include <assert.h>
#include "merge_sam.h"
#include "merge_sam_heap.h"
#include "sam2pretty_lib.h"

#define ABS(a)	   (((a) < 0) ? -(a) : (a))

typedef struct pp_ll pp_ll;
struct pp_ll {
	pretty * head;
	pretty * tail;
	size_t length;
};

int id=0;

static int qsort_sam_header_strcmp(const void *a, const void *b)  { 
    const char **ia = (const char **)a;
    const char **ib = (const char **)b;
    const char *s = *ia;
    const char *t = *ib;
    if (s[1]=='H' && s[2]=='D' && (t[1]!='H' && t[2]!='D')) {
	return -1;
    }
    if (t[1]=='H' && t[2]=='D' && (s[1]!='H' && s[2]!='D')) {
	return 1;
    }
    if (s[1]=='S' && s[2]=='Q' && (t[1]!='S' && t[2]!='Q')) {
	return -1;
    }
    if (t[1]=='S' && t[2]=='Q' && (s[1]!='S' && s[2]!='Q')) {
	return 1;
    }
    if (s[1]=='R' && s[2]=='G' && (t[1]!='R' && t[2]!='G')) {
	return -1;
    }
    if (t[1]=='R' && t[2]=='G' && (s[1]!='R' && s[2]!='G')) {
	return 1;
    }
    if (s[1]=='P' && s[2]=='G' && (t[1]!='P' && t[2]!='G')) {
	return -1;
    }
    if (t[1]=='P' && t[2]=='G' && (s[1]!='P' && s[2]!='G')) {
	return 1;
    }
    return -strcmp(s, t);
} 


size_t file_buffer_size = DEF_FILE_BUFFER_SIZE;
uint64_t current_read = 0; //1 based, 0 is invalid

file_buffer reads_file_buffer;
file_buffer * file_buffers;

pp_ll ** file_ll;
pp_ll * master_ll;
char ** read_names;
uint64_t read_names_size;
int64_t read_names_filled;


static int phred_scale(double p_value) {
        double p_error = 1 - p_value;
        double phred_score = -10*log(p_error)/log(10);
	if (isinf(phred_score)) {
		if (phred_score>0) {
			return 254;
		} else {
			return 0;
		}
	}
	int ret = (int)(round(phred_score));
	if (ret<0) {
		ret=0;
	} else if (ret>254) {
		ret=254;
	}
        return ret;
}

static double boltzman(int score, mapq_info * mqi) {
	return exp(((double)score)*mqi->score_matrix_lambda/mqi->sw_scale_constant);
}

/*
static void compute_mapq(heap_pa * hpa, mapq_info * mqi, output_filter * of, int length, int top_score, int best_isize) { 
	if (mqi->top_hits<of->number_outputs) {
		fprintf(stderr,"mapq tophits must be greater then or qual to number of outputs!\n");
		exit(1);
	}
	double mapq_norm_first=0.0;
	double mapq_norm_second=0.0;
	double mapq_norm_unpaired=0.0;
	heap_pa mapq_heap;
	heap_pa_init(&mapq_heap, length);
	//first round of normalization
	memcpy(mapq_heap.array,hpa->array,sizeof(heap_pa_elem)*length);
	mapq_heap.load=hpa->load;
	int looked_at=0;
	while (mapq_heap.load>0 && (looked_at<mqi->top_hits || mqi->top_hits==0)) {
		heap_pa_elem tmp;
		heap_pa_extract_max(&mapq_heap,&tmp);
		pretty * pa=tmp.rest;
		if (!mqi->strata || (tmp.score==top_score && tmp.isize==best_isize)) {
			if (pa->paired_sequencing && pa->mapped && pa->mate_pair->mapped) {
				pretty * pa_first = pa; 
				pretty * pa_second = pa->mate_pair;
				if (pa->second_in_pair) {
					pa_first = pa->mate_pair;
					pa_second = pa;
				}
				pa_first->temp_mapq=boltzman(pa_first->score, mqi);
				pa_second->temp_mapq=boltzman(pa_second->score, mqi);
				mapq_norm_first+=pa_first->temp_mapq;
				mapq_norm_second+=pa_second->temp_mapq;
				looked_at++;
			} else if (pa->mapped) {
				pa->temp_mapq=boltzman(pa->score, mqi);
				mapq_norm_unpaired+=pa->temp_mapq;
				looked_at++;
			}	
		} else {
			break;
		}
	}
	//paired normalization
	double mapq_norm_paired=0.0;
	memcpy(mapq_heap.array,hpa->array,sizeof(heap_pa_elem)*length);
	mapq_heap.load=hpa->load;
	looked_at=0;
	while (mapq_heap.load>0 && (looked_at<mqi->top_hits || mqi->top_hits==0)) {
		heap_pa_elem tmp;
		heap_pa_extract_max(&mapq_heap,&tmp);
		pretty * pa=tmp.rest;
		if (!mqi->strata || (tmp.score==top_score && tmp.isize==best_isize)) {
			if (pa->paired_sequencing && pa->mapped && pa->mate_pair->mapped) {
				pretty * pa_first = pa; 
				pretty * pa_second = pa->mate_pair;
				if (pa->second_in_pair) {
					pa_first = pa->mate_pair;
					pa_second = pa;
				}
				mapq_norm_paired+=
					pa_first->temp_mapq*pa_second->temp_mapq/(mapq_norm_first*mapq_norm_second);
				looked_at++;
			} else if (pa->mapped) {
				looked_at++;
			}	
		} else {
			break;
		}
	}
	//assignment
	memcpy(mapq_heap.array,hpa->array,sizeof(heap_pa_elem)*length);
	mapq_heap.load=hpa->load;
	looked_at=0;
	while (mapq_heap.load>0 && (looked_at<mqi->top_hits || mqi->top_hits==0)) {
		heap_pa_elem tmp;
		heap_pa_extract_max(&mapq_heap,&tmp);
		pretty * pa=tmp.rest;
		if (!mqi->strata || (tmp.score==top_score && tmp.isize==best_isize)) {
			if (pa->paired_sequencing && pa->mapped && pa->mate_pair->mapped) {
				pretty * pa_first = pa; 
				pretty * pa_second = pa->mate_pair;
				if (pa->second_in_pair) {
					pa_first = pa->mate_pair;
					pa_second = pa;
				}
				int mapq=phred_scale((pa_first->temp_mapq*pa_second->temp_mapq/(mapq_norm_first*mapq_norm_second))/mapq_norm_paired);
				pa_first->mapq=mapq;	
				pa_second->mapq=mapq;		
				looked_at++;
			} else if (pa->mapped) {
				pa->mapq=phred_scale(pa->temp_mapq/mapq_norm_unpaired);
				looked_at++;
			}
		} else {
			break;
		}
	}
	heap_pa_destroy(&mapq_heap);
}
*/

static void pp_ll_zero(pp_ll * ll) {
	ll->head=NULL;
	ll->tail=NULL;
	ll->length=0;
}

static void pp_ll_append(pp_ll * ll,pretty * pa) {
	if (ll->head==NULL) {
		ll->head=pa;
		ll->length=0;
	} else {
		ll->tail->next=pa;
	}
	ll->tail=pa;
	pa->next=NULL;
	ll->length++;	
} 

static void pp_ll_combine(pp_ll * dest, pp_ll * source) {
	if (dest->length==0) {
		*dest=*source;
	} else if (source->length!=0) {
		dest->length+=source->length;
		dest->tail->next=source->head;
		dest->tail=source->tail;
	}
	source->length=0;
	source->head=NULL;
	source->tail=NULL;
}

static void reads_to_heap(heap_pa * hpa, pretty * pa, int * best_score, int * best_alignments, int * best_idist, output_filter * of) {
	*best_score=-INT_MAX;
	*best_alignments=0;
	*best_idist=INT_MAX;
	heap_pa_elem heap_e;
	while (pa!=NULL) {
		heap_e.rest=pa;
		heap_e.score=pa->score + ((pa->mate_pair!=NULL) ? pa->mate_pair->score : 0);
		heap_e.idist=of->use_isize ? ABS(ABS(pa->isize)-of->isize) : 0;
		if (heap_e.score>*best_score || (heap_e.score==*best_score && heap_e.idist<*best_idist)) {
			*best_score=heap_e.score;	
			*best_idist=heap_e.idist;
			*best_alignments=1;
		} else if (heap_e.score==*best_score && heap_e.idist==*best_idist) {
			(*best_alignments)++;
		}
		heap_pa_insert(hpa, &heap_e);
		pa=pa->next;
	}
}

static void add_marked_free_rest(pp_ll * dest, pp_ll * src) {
	pretty * pa=src->head;
	while (pa!=NULL) {
		pretty * next=pa->next;
		if (pa->mark) {
			pp_ll_append(dest,pa);
		} else {
			//free it and its mate if possible
			if (pa->mate_pair!=NULL) {
				pretty_free(pa->mate_pair);
				pa->mate_pair=NULL;
			}
			pretty_free(pa);
		}
		pa=next;
	}
}


static void process_read(pp_ll * ll, output_filter * of, mapq_info * mqi) {
	//Do one pass over reads to pair up and get a count
	//For how many entrys there are
	pp_ll paired;
	pp_ll_zero(&paired);
	pp_ll only_first_mapped;
	pp_ll_zero(&only_first_mapped);
	pp_ll only_second_mapped;
	pp_ll_zero(&only_second_mapped);
	pp_ll unpaired;
	pp_ll_zero(&unpaired);
	pp_ll unmapped;
	pp_ll_zero(&unmapped);
	int paired_read=-1;
	pretty * pa =ll->head;
	assert(pa!=NULL);
	while (pa!=NULL) {
		pa->mark=false;
		pretty * next = pa->next;
		if (pa->paired_sequencing) {
			if (paired_read==0) {
				fprintf(stderr,"merge_sam : process_read, read_name %s, cannot be both paired in sequencing and not!\n",pa->read_name);
				pretty_print_sam(stderr,pa);
				exit(1);
			}	
			paired_read=1;
			//set up mate_pair links
			pa->mate_pair=pa->next;
			if (pa->mate_pair==NULL) {
				fprintf(stderr,"merge_sam : process_read, mate_pair is null!\n");
				pretty_print_sam(stderr,pa);
				exit(1);
			}
			pa->mate_pair->mate_pair=pa;
			next=pa->mate_pair->next;
			if (pa->second_in_pair) {
				pa=pa->mate_pair;
			}
			assert(pa->first_in_pair);
			if (pa->isize!=-pa->mate_pair->isize) {
				fprintf(stderr,"merge_sam : process_read, pairs insert sizes are not right!\n");
				pretty_print_sam(stderr,pa);
				pretty_print_sam(stderr,pa->mate_pair);
				exit(1);	
			}	
			//check reference name too? too complicated for check at this point
			//leave to user	
			if (pa->mate_genome_start_unpadded!=pa->mate_pair->genome_start_unpadded ||
			    pa->genome_start_unpadded!=pa->mate_pair->mate_genome_start_unpadded) { 
				fprintf(stderr,"merge_sam : process_read, mate_pos and pos disagree!\n");
				pretty_print_sam(stderr,pa);
				pretty_print_sam(stderr,pa->mate_pair);
				exit(1);
			}
			//categorize it properly
			assert(pa->first_in_pair);
			assert(pa->mate_pair->second_in_pair);
			if (!pa->mapped && !pa->mate_pair->mapped) {
				//goes into unmapped
				if (unmapped.length==0) {
					pp_ll_append(&unmapped,pa);
				} else {
					if (pa->mate_pair!=NULL) {
						pretty_free(pa->mate_pair);
					}				
					pretty_free(pa);
				}
			} else if (!pa->mapped) { 
				//this means that the mate pair is mapped
				pp_ll_append(&only_second_mapped,pa);
			} else if (!pa->mate_pair->mapped) {
				//goes into half paired
				pp_ll_append(&only_first_mapped,pa);
			} else {
				//goes into fully mapped
				pp_ll_append(&paired,pa);
			}
		} else {
			//goes into unpaired
			if (paired_read==1) {
				fprintf(stderr,"merge_sam : process_read, read_name %s, cannot be both paired in sequencing and not!\n",pa->read_name);
				pretty_print_sam(stderr,pa);
				exit(1);
			}
			paired_read=0;
			if (pa->mapped) {
				pp_ll_append(&unpaired,pa);
			} else {
				//goes into unmapped
				pp_ll_append(&unmapped,pa);
			}
		}
		pa=next;
	}
	//assert(length==ll->length);
	/*pa=paired.head;
	while(pa!=NULL) {
		pretty_print_sam(stdout,pa);
		pretty_print_sam(stdout,pa->mate_pair);
		pa=pa->next;
	}
	fprintf(stdout,"Only first mapped entires:\n");
	pa=only_first_mapped.head;
	while(pa!=NULL) {
		pretty_print_sam(stdout,pa);
		pretty_print_sam(stdout,pa->mate_pair);
		pa=pa->next;
	}
	fprintf(stdout,"Only second mapped entires:\n");
	pa=only_second_mapped.head;
	while(pa!=NULL) {
		pretty_print_sam(stdout,pa);
		pretty_print_sam(stdout,pa->mate_pair);
		pa=pa->next;
	}
	fprintf(stdout,"UnPaired entires:\n");
	pa=unpaired.head;
	while(pa!=NULL) {
		pretty_print_sam(stdout,pa);
		pa=pa->next;
	}
	fprintf(stdout,"UnMappped entires:\n");
	pa=unmapped.head;
	while(pa!=NULL) {
		pretty_print_sam(stdout,pa);
		if (pa->mate_pair!=NULL) {
			pretty_print_sam(stdout,pa->mate_pair);
		}
		pa=pa->next;
	}*/	
	assert(only_first_mapped.length==only_second_mapped.length);
	//If there is at least one unpaired int sequencing and mapped read
	int best_score[2]={-INT_MAX,-INT_MAX};
	int best_alignments[2]={0,0};
	int best_idist[2]={INT_MAX,INT_MAX};
	int alignments[2]= {0,0};
	heap_pa hpa[2];
	hpa[0].array=NULL; hpa[1].array=NULL;
	hpa[0].load=0; hpa[1].load=0;
	assert(only_first_mapped.length==only_second_mapped.length);
	if (unpaired.length>0) {
		heap_pa_init(hpa+0,unpaired.length);
		assert(paired_read==0);
		reads_to_heap(hpa+0, unpaired.head, best_score+0, best_alignments+0, best_idist+0, of);
		alignments[0]=unpaired.length;
	} else if (paired.length>0) {
		//paired mode
		heap_pa_init(hpa+0,paired.length);
		assert(paired_read==1);
		reads_to_heap(hpa+0, paired.head, best_score+0, best_alignments+0, best_idist+0, of);
		alignments[0]=paired.length;
	} else if (only_first_mapped.length>0) {
		//first in pair
		heap_pa_init(hpa+0,only_first_mapped.length);
		//second in pair
		heap_pa_init(hpa+1,only_second_mapped.length);
		if (of->half_paired) {
			reads_to_heap(hpa+0, only_first_mapped.head, best_score+0, best_alignments+0, best_idist+0, of);
			alignments[0]=only_first_mapped.length;
			reads_to_heap(hpa+1, only_second_mapped.head, best_score+1, best_alignments+1, best_idist+1, of);
			alignments[1]=only_second_mapped.length;
		}
	} else if (unmapped.length>0) {
		//unmapped
		assert(unmapped.length==1);
		heap_pa_init(hpa+0,unmapped.length);
		reads_to_heap(hpa+0, unmapped.head, best_score+0, best_alignments+0, best_idist+0, of);
		if (of->unaligned) {
			hpa[0].load=1;
			alignments[0]=1;
		} else {
			hpa[0].load=0;
			alignments[0]=0;
		}
	} else {
		fprintf(stderr,"merge_sam: process_read, a fatal error has occured, cannot find read!\n");
		exit(1);
	}
	//If we should analyze this read
	//If max alignments set to all OR (not strata and found less then or equal to MAX)
	//OR (strata and best_alignments less then or equal to MAX)
	int i;
	for (i=0; i<2; i++) {
		if ( of->max_alignments==0 || (!of->strata && alignments[i]<=of->max_alignments) || (of->strata && best_alignments[i]<=of->max_alignments)) {
			//If we need to compute MAPQ, lets give it a shot
			//if (mqi->calculate && has_mapping) {
			//	compute_mapq(&hpa,mqi,of,length,top_score,best_isize);
			//}
			//Lets find out what to print and set the SAM strings accordingly
			int printed_alignments=0;
			while (hpa[i].load>0 && ( printed_alignments<of->number_outputs || of->number_outputs==0) && (!of->strata || printed_alignments<best_alignments[i]) ) {
				heap_pa_elem heap_e;
				heap_pa_extract_max(hpa+i,&heap_e);
				pretty * pa=heap_e.rest;
				pa->mark=true; 
				pretty_print_sam_update(pa);
				if (pa->mate_pair!=NULL) {
					pa->mate_pair->mark=true; 
					pretty_print_sam_update(pa->mate_pair);
				}
				printed_alignments++;
			}
		}
	}
	heap_pa_destroy(hpa+0);
	heap_pa_destroy(hpa+1);
	//clean up all alignments that will not be printed
	pp_ll new_ll;
	pp_ll_zero(&new_ll);
	if (unpaired.length>0) {
		add_marked_free_rest(&new_ll, &unpaired);
	}
	if (paired.length>0) {
		add_marked_free_rest(&new_ll, &paired);
	}
	if (only_first_mapped.length>0) {
		add_marked_free_rest(&new_ll, &only_first_mapped);
		add_marked_free_rest(&new_ll, &only_second_mapped);
	}
	if (unmapped.length>0) {
		add_marked_free_rest(&new_ll, &unmapped);
	}
	*ll=new_ll;	
}

static int read_more(file_buffer * fb) {
	//shift over what we have no not used
	size_t i;
	for (i=0; i<fb->left_over; i++) {
		if (fb->buffer[fb->filled-fb->left_over+i]=='\0') {
			fb->buffer[fb->filled-fb->left_over+i]='\n';
		}
		fb->buffer[i]=fb->buffer[fb->filled-fb->left_over+i];
	}
	fb->buffer[i]='\0'; //just nice for debugging
	//fprintf(stderr,"BUFFER: |%s| %d\n",fb->buffer,fb->left_over);
	//read in more
	int ret = fb->left_over + gzread(fb->file,fb->buffer+i,sizeof(char)*fb->size-i);
	fb->left_over=0;
	/*if (ret<fb->size-i && !gzeof(fb->file)) {
		fprintf(stderr,"A fatal error has occured in reading a file, %s\n");
		perror("");
		exit(1);
	}*/
	fb->filled=ret;
	fb->buffer[ret]='\0';
	//fprintf(stderr,"BUFFERAFT: |%s|\n",fb->buffer);
	if (gzeof(fb->file)) {
		fb->eof=1;
	}
	return ret;
}


static int get_sam_header_read(file_buffer * fb, pretty ** pa,int * last_used) {
	//Read in more from the alginments file
	read_more(fb);
	int last=0; size_t current=0; char * line=fb->buffer;
	//While we have not processed to end of filled buffer, keep processing
	pretty * tail;
	int headers_found=0;
	while (current<fb->filled) {
		//find the next newline
		while (current<fb->filled && fb->buffer[current]!='\n') { current++; };
		//if we have a full SAM record in memory
		if (current<fb->filled) {
			fb->buffer[current]='\0';
			//line starts after where last line finished
			line = fb->buffer + last;
			if ( line[0]!='@' && line[0]!='\n' ) {
				fb->left_over=fb->filled - last;
				return headers_found;
			} else {
				if (*pa==NULL) {
					*pa=pretty_new();
					tail=*pa;
					tail->sam_string=strdup(line);
				} else {
					tail->next=pretty_new();
					tail=tail->next;
					tail->sam_string=strdup(line);
				}
				headers_found++;
			}
			current++;
			last=current;
		}
	}
	fb->left_over=fb->filled - last;
	if (current>0 && headers_found==0) {
		fprintf(stderr,"Failed to fit a SAM alignment into memory! Please increase buffer size!\n");	
		exit(1);
	}
	return -1;
}
static int get_sam_header(file_buffer * fb, pretty ** pas) {
	int last_used=0;
	int ret;
	while (fb->eof==0) {
		if ((ret=get_sam_header_read(fb,pas,&last_used))>-1) {
			return ret;
		}
	}
	if (fb->eof==1) {
		if ((ret=get_sam_header_read(fb,pas,&last_used))>-1) {
			return ret;
		}
	}
	fb->left_over=0;
	return 0;
}
static int get_sam_entries_read(file_buffer * fb, pp_ll * ll,int * last_used) {
	//Read in more from the alginments file
	read_more(fb);
	int last=0; size_t current=0; char * line=fb->buffer;
	int read_alignments=0;
	//While we have not processed to end of filled buffer, keep processing
	while (current<fb->filled) {
		//find the next newline
		while (current<fb->filled && fb->buffer[current]!='\n') { current++; };
		//if we have a full SAM record in memory
		if (current<fb->filled) {
			read_alignments++;
			fb->buffer[current]='\0';
			//line starts after where last line finished
			line = fb->buffer + last;
			//if its a SAM header skip it
			if ( line[0]!='@' ) {
				//make a copy of the line for sam2pretty library
				char * line_copy=strdup(line);
				pretty * pa = pretty_from_string(line_copy);
				if (pa==NULL) {
					fprintf(stderr,"A fatal error has occured in parsing the following SAM line, \n%s\n",line);
					exit(1);
				}	
				int j; int found=0;
				//since alignments are in same order as reads, no need to compare to previous,
				//just resume where we left of last time
				for (j=*last_used; found==0 && j<read_names_filled; j++) {
					//TODO store length!
					size_t length=strlen(pa->read_name);
					//printf("%s vs %s\n",read_names[j],pa->read_name);
					if (strncmp(read_names[j],pa->read_name,length)==0 && (!isdigit(read_names[j][length]) || !isdigit(read_names[j][length-1]))) {//read_names[j][length-1]=='F' || read_names[j][length-1]=='R')) {
						pp_ll_append(ll+j,pa);
						found=1;
						*last_used=j;
					}
				}
				//free the line copy
				free(line_copy);
				if (found==0) {
					pretty_free(pa);
					fb->left_over=fb->filled - (line-fb->buffer);
					return 0;
				}
			} else {
				fprintf(stderr,"Expecting alignments, found header!\n");
				exit(1);
			}
			current++;
			last=current;
		}
	}
	fb->left_over=fb->filled - last;
	if (current>0 && read_alignments==0) {
		fprintf(stderr,"Failed to fit a SAM alignment into memory! Please increase buffer size!\n");	
		exit(1);
	}
	return 1;
}

static int get_sam_entries(file_buffer * fb, pp_ll * ll) {
	int last_used=0;
	while (fb->eof==0) {
		if (get_sam_entries_read(fb,ll,&last_used)==0) {
			return 0;
		}
	}
	if (fb->eof==1) {
		if (get_sam_entries_read(fb,ll,&last_used)==0) {
			return 0;
		}
	}
	fb->left_over=0;
	return 0;
}


static int index_new_read_names(bool fastq_reads,bool colour_space) {
	int last=0; size_t current=0; int index=0;
	while(reads_file_buffer.buffer[current]=='\n') { current ++; last++;};
	while (current<reads_file_buffer.filled) {
		while (current<reads_file_buffer.filled && reads_file_buffer.buffer[current]!='\n') { current++; };
		if (current<reads_file_buffer.filled) {
			reads_file_buffer.buffer[current]='\0';
		//	fprintf(stderr,"R:%s\n",reads_file_buffer.buffer + last);
			read_names[index++]=reads_file_buffer.buffer + last;
			current++;
			last=current;
		}
	}
	read_names[index++]=reads_file_buffer.buffer + last;
	//figure out which ones are valid
	int line=0; read_names_filled=0; int in_header=1;
	if (fastq_reads) {
		while (line<index) {
			if (in_header==1) {
				if (read_names[line][0]=='#' || read_names[line][0]=='\n') {
					//skip this
					line++;
				} else {
					in_header=0;
				}
			}
			if (in_header==0) {
				read_names[read_names_filled]=read_names[line++]+1;
				//fprintf(stderr,"Found read %s\n",read_names[read_names_filled-1]);	
				//get the sequence
				int seq_length=0;
				while (line<index && read_names[line][0]!='+') {
					size_t length = strlen(read_names[line]);
					if (isdigit(read_names[line][length-1])) {
						colour_space=true;
					}
					seq_length+=length;
					line++;
				}
				//skip this line
				line++;
				int qual_length=0;
				while (line<index && !(qual_length==seq_length || (colour_space && qual_length+1==seq_length))) {
					qual_length+=strlen(read_names[line++]);
				}
				if (qual_length==0 || seq_length==0 || (line>=index && !(qual_length==seq_length || (colour_space && qual_length+1==seq_length)) )) {
					//didn't get all of last read
					reads_file_buffer.left_over=reads_file_buffer.filled - (read_names[read_names_filled]-1-reads_file_buffer.buffer);
					return 0;
				} else {
					//got this read keep going
					while((*read_names[read_names_filled]==' ' || *read_names[read_names_filled]=='\t') && *read_names[read_names_filled]!='\n') { read_names[read_names_filled]++; };
					read_names_filled++;
					in_header=1;	
				}
			}
		}
	} else {
		while (line<index) {
			if (in_header==1) {
				if (read_names[line][0]=='#' || read_names[line][0]=='\n') {
					//skip this
					line++;
				} else {
					in_header=0;
				}
			}
			if (in_header==0) {
				read_names[read_names_filled]=read_names[line++]+1;
				while (line<index && read_names[line][0]!='>' && read_names[line][0]!='\n') { line++; };
				if (line>=index && reads_file_buffer.eof!=1) {
					//didn't get all of last read
					reads_file_buffer.left_over=reads_file_buffer.filled - (read_names[read_names_filled]-1-reads_file_buffer.buffer);
					return 0;
				} else {
					//got this read keep going
					while((*read_names[read_names_filled]==' ' || *read_names[read_names_filled]=='\t') && *read_names[read_names_filled]!='\n') { read_names[read_names_filled]++; };
					read_names_filled++;
					in_header=1;	
				}
			}
		}
	}	
	return 0;
}

static int find_last_read(bool fastq_reads,bool colour_space)  {
	index_new_read_names(fastq_reads,colour_space);
	int64_t i;
	//fprintf(stderr,"have %d names filled \n",read_names_filled);
	if (read_names_filled==0) {
		return 0;
	}
	int new_index=0;
	for (i=1; i<read_names_filled; i++) {
		if (strcmp(read_names[new_index],read_names[i])!=0) {
			read_names[++new_index]=read_names[i];
		}
	}
	read_names_filled=new_index+1;
	for (i=0; i<read_names_filled; i++) {
	}
	return 0;
}

void merge_sam(char * reads_filename, bool fastq_reads, bool colour_space, int read_threads, int compute_threads, int number_of_files, char ** filenames, FILE * output_file, output_filter of, mapq_info mqi,int buffer_size) {
	file_buffer_size=buffer_size*1024*1024;
	int i;
	//allocate file buffers
	file_buffers=(file_buffer*)malloc(sizeof(file_buffer)*number_of_files);
	if (file_buffers==NULL) {
		fprintf(stderr,"Failed to malloc space for file_buffer structures\n");
		exit(1);
	}
	memset(file_buffers,0,sizeof(file_buffer)*number_of_files);
	//open the files
	pretty * pretty_sam_headers[number_of_files];
	memset(pretty_sam_headers,0,sizeof(pretty*)*number_of_files);
	int headers_found=0;
	for (i=0; i<number_of_files; i++) {
		file_buffers[i].buffer=(char*)malloc(file_buffer_size+1);
		if (file_buffers[i].buffer==NULL) {
			fprintf(stderr,"Failed to allocate buffer %d, in merge_sam.c\n",i);
			exit(1);
		}
		file_buffers[i].size=file_buffer_size;
		file_buffers[i].filled=0;
		file_buffers[i].left_over=0;
		file_buffers[i].last_read=0;
		file_buffers[i].eof=0;
		file_buffers[i].file=gzopen(filenames[i],"r");	
		if (file_buffers[i].file==NULL) {
			fprintf(stderr,"Failed to open file %s\n",filenames[i]);
			perror("");
			exit(1);
		}
		headers_found+=get_sam_header(&file_buffers[i],pretty_sam_headers+i);
	}
	if (headers_found>0){
		char * sam_headers[headers_found];
		int index=0;
		for (i=0; i<number_of_files; i++) {
			pretty * pa = pretty_sam_headers[i];
			while(pa!=NULL) {
				sam_headers[index++]=strdup(pa->sam_string);
				pretty * tmp = pa->next;
				pretty_free(pa);
				pa=tmp;
			}	
		}
		if (!of.detect_isize && !of.header_provided) {
			qsort(sam_headers,headers_found,sizeof(char*),qsort_sam_header_strcmp);
			fprintf(output_file,"%s\n",sam_headers[0]);
			for (i=1; i<headers_found; i++) {
				if (strcmp(sam_headers[i-1],sam_headers[i])!=0) {
					fprintf(output_file,"%s\n",sam_headers[i]);
				}
			}
		}
		for (i=0; i<headers_found; i++) {
			free(sam_headers[i]);
		}
	}
	//allocate memory for reads_file_buffer	
	reads_file_buffer.buffer=(char*)malloc(file_buffer_size+1);
	if (reads_file_buffer.buffer==NULL) {
		fprintf(stderr,"Failed to malloc memory for reads file buffer\n");
		exit(1);
	}
	reads_file_buffer.size=file_buffer_size;
	reads_file_buffer.left_over=0;
	reads_file_buffer.filled=0;
	reads_file_buffer.last_read=0;
	reads_file_buffer.eof=0;
	//open the reads file
        reads_file_buffer.file=gzopen(reads_filename,"r");
	if (reads_file_buffer.file==NULL) {
		fprintf(stderr,"Failed to open reads file, %s\n",reads_filename);
		perror("");
		exit(1);
	}
	//allocate space for read names
	read_names_size=file_buffer_size;
	read_names_filled=0;
	read_names=(char**)malloc(sizeof(char*)*read_names_size);
	if (read_names==NULL) {
		fprintf(stderr,"Failed to allocate memory for read_names\n");
		exit(1);
	}	
	memset(read_names,0,read_names_size*sizeof(char*));	
	file_ll = (pp_ll**)malloc(sizeof(pp_ll*)*number_of_files);
	if (file_ll==NULL) {
		fprintf(stderr,"Failed to allocate memory for file_ll\n");
		exit(1);
	}
	memset(file_ll,0,sizeof(pp_ll*)*number_of_files);
	//allocate space for histogram if needed
	uint64_t * isizes = (uint64_t*)malloc(sizeof(uint64_t)*of.max_isize);
	if (isizes==NULL) {
		fprintf(stderr,"Failed to allocate memory for isizes buffer\n");
		exit(1);
	}
	memset(isizes,0,sizeof(uint64_t)*of.max_isize);
	long long reads_processed=0;
	//Keep processing until no more reads avaliable to process!
	while (read_more(&reads_file_buffer)) {
		find_last_read(fastq_reads,colour_space);
		reads_processed+=read_names_filled;
		if (read_names_filled==0) {
			break;
		}
		int i;
		//for (i=0; i<read_names_filled; i++) {
		//	fprintf(stderr,"%d R:%s\n",read_names_filled,read_names[i]);
		//}
		//exit(1);
		fprintf(stderr,"Processed %lld reads so far\n",reads_processed);
		#pragma omp parallel for schedule(static,1) num_threads(read_threads)
		for (i=0; i<number_of_files; i++) {
			//fprintf(stderr,"%d: with file %d\n",omp_get_thread_num(),i);
			file_ll[i]=(pp_ll*)malloc(sizeof(pp_ll)*read_names_filled);
			if (file_ll[i]==NULL) {
				fprintf(stderr,"Failed to allocate memory for file_ll[i]\n");
				exit(1);
			}
			memset(file_ll[i],0,sizeof(pp_ll)*read_names_filled);
			get_sam_entries(&file_buffers[i],file_ll[i]);	
			int j;
			for (j=0; j<read_names_filled; j++) {
			}
		}
		master_ll=(pp_ll*)malloc(sizeof(pp_ll)*read_names_filled);
		if (master_ll==NULL) {
			fprintf(stderr,"Failed to allocate memory for master_ll\n");
			exit(1);
		}
		memset(master_ll,0,sizeof(pp_ll)*read_names_filled);
		for (i=0; i<read_names_filled; i++) {
			int j;
			for (j=0; j<number_of_files; j++) {
				pp_ll_combine(master_ll+i,file_ll[j]+i);
			}
		}
		#pragma omp parallel for num_threads(compute_threads) schedule(static,50000)
		for (i=0; i<read_names_filled; i++) {
			pp_ll * ll = master_ll+i;
			if (ll->length>0) {
				process_read(ll,&of,&mqi);
			}
		}
		for (i=0; i<read_names_filled; i++) {
			pretty * pa = master_ll[i].head;
			while (pa!=NULL) {
				if (of.detect_isize) {
					if (pa->mate_pair!=NULL && pa->mate_pair->mapped && pa->mapped) {
						isizes[abs(pa->isize)]++;
					}
					if (pa->mate_pair!=NULL) {
						pretty_free(pa->mate_pair);
					}
				} else { 
					fprintf(output_file,"%s",pa->sam_string);
					if (pa->mate_pair!=NULL) {
						fprintf(output_file,"%s",pa->mate_pair->sam_string);
						pretty_free(pa->mate_pair);
					}	
				}
				pretty * tmp = pa->next; 
				pretty_free(pa);
				pa=tmp;
			}
		}
		for (i=0; i<number_of_files; i++) {
			free(file_ll[i]);
			file_ll[i]=NULL;
		}
		free(master_ll);
	}	
	for (i=0; i<number_of_files; i++) {
		free(file_buffers[i].buffer);
		gzclose(file_buffers[i].file);	
	}
	fprintf(stderr,"Processed %lld reads, reached-eof: %s\n",reads_processed,reads_file_buffer.eof==1 ? "Yes" : "No");
	free(file_ll);
	free(file_buffers);	
	gzclose(reads_file_buffer.file);
	free(reads_file_buffer.buffer);
	free(read_names);
	if (of.detect_isize) {
		uint64_t sum=0;
		uint64_t entries=0;
		uint64_t i;
		for (i=0; i<of.max_isize; i++) {
			sum+=i*isizes[i];
			entries+=isizes[i];
		}
		double mean = ((double)sum)/((double)entries);
		double std_sum=0.0;
		for (i=0; i<of.max_isize; i++) {
			std_sum+=isizes[i]*(i-mean)*(i-mean);
		}
		double stddev=sqrt(std_sum/entries);
		printf("Mean: %.5e,\tStandard-Deviation: %.5e\n",mean,stddev);
	}
	free(isizes);
	return;
}
