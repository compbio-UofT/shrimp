#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <omp.h>
#include <math.h>
#include <ctype.h>
#include <zlib.h>
#include <limits.h>
#include "merge_sam.h"
#include "merge_sam_heap.h"
#include "sam2pretty_lib.h"


static int qsort_sam_header_strcmp(const void *a, const void *b)  { 
    const char **ia = (const char **)a;
    const char **ib = (const char **)b;
    const char *s = *ia;
    const char *t = *ib;
    if (s[1]=='H' && s[2]=='D' && (t[1]!='H' && t[2]!='D')) {
	return 1;
    }
    if (t[1]=='H' && t[2]=='D' && (s[1]!='H' && s[2]!='D')) {
	return -1;
    }
    if (s[1]=='S' && s[2]=='Q' && (t[1]!='S' && t[2]!='Q')) {
	return 1;
    }
    if (t[1]=='S' && t[2]=='Q' && (s[1]!='S' && s[2]!='Q')) {
	return -1;
    }
    if (s[1]=='R' && s[2]=='G' && (t[1]!='R' && t[2]!='G')) {
	return 1;
    }
    if (t[1]=='R' && t[2]=='G' && (s[1]!='R' && s[2]!='G')) {
	return -1;
    }
    if (s[1]=='P' && s[2]=='G' && (t[1]!='P' && t[2]!='G')) {
	return 1;
    }
    if (t[1]=='P' && t[2]=='G' && (s[1]!='P' && s[2]!='G')) {
	return -1;
    }
    return strcmp(s, t);
} 


size_t file_buffer_size = DEF_FILE_BUFFER_SIZE;
uint64_t current_read = 0; //1 based, 0 is invalid

file_buffer reads_file_buffer;
file_buffer * file_buffers;
pretty *** file_pa;
pretty ** master_pa;
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



static void process_read(pretty ** ppa, output_filter * of, mapq_info * mqi) {
	pretty * pa = *ppa;
	int length=0;
	pretty * head = pa;
	//Do one pass over reads to pair up and get a count
	//For how many entrys there are
	while (pa!=NULL) {
		length++;
		pa->mark=false;
		if (pa->paired_sequencing) {
			pa->mate_pair=pa->next;
			if (pa->mate_pair==NULL) {
				fprintf(stderr,"merge_sam : process_read, mate_pair is null!\n");
				pretty_print_sam(stdout,pa);
				exit(1);
			}
			if (pa->isize!=-pa->mate_pair->isize) {
				fprintf(stderr,"merge_sam : process_read, pairs insert sizes are not right!\n");
				pretty_print_sam(stdout,pa);
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
			pa->next=pa->next->next;
		}
		pa=pa->next;
	}
	//can actually make heap smaller and drop off elements? TODO
	//Put everything onto a heap now that we know the size
	heap_pa hpa;
	heap_pa_init(&hpa, length);
	pa=head;
	int top_score=-100; //TODO
	bool has_mapping=false;
	int best_alignments=0;
	int best_isize=INT_MAX;
	while (pa!=NULL) {
		bool look_at_hit=true;
		if ( !of->half_paired && pa->mapped && pa->mate_pair!=NULL && !pa->mate_pair->mapped) {
			look_at_hit=false;
		}
		if ( !of->half_paired && !pa->mapped && pa->mate_pair!=NULL && pa->mate_pair->mapped) {
			look_at_hit=false;
		}
		if ( !of->unaligned && !pa->mapped && pa->mate_pair==NULL) {
			look_at_hit=false;
		}
		if ( !of->unaligned && !pa->mapped && pa->mate_pair!=NULL && !pa->mate_pair->mapped) {
			look_at_hit=false;
		}
		if (look_at_hit) {
			heap_pa_elem tmp;
			tmp.rest=pa;
			tmp.score=pa->score;
			tmp.isize=of->use_isize ? abs(abs(pa->isize)-of->isize) : 0;
			has_mapping=has_mapping || pa->mapped;	
			if (pa->mate_pair!=NULL && pa->mate_pair->has_score) {
				tmp.score+=pa->mate_pair->score;
				has_mapping=has_mapping || pa->mate_pair->mapped;
			}
			if (tmp.score>top_score || (tmp.score==top_score && tmp.isize<best_isize)) {
				top_score=tmp.score;
				best_isize=tmp.isize;
				best_alignments=1;
			} else if (tmp.score==top_score && tmp.isize==best_isize) {
				best_alignments++;
			}
			heap_pa_insert(&hpa, &tmp);
		}
		pa=pa->next;
	}
	//If we should analyze this read
	//If max alignments set to all OR (not strata and found less then or equal to MAX)
	//OR (strata and best_alignments less then or equal to MAX)
	if ( of->max_alignments==0 || (!of->strata && length<=of->max_alignments) || (of->strata && best_alignments<=of->max_alignments)) {
		//If we need to compute MAPQ, lets give it a shot
		if (mqi->calculate && has_mapping) {
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
			memcpy(mapq_heap.array,hpa.array,sizeof(heap_pa_elem)*length);
			mapq_heap.load=hpa.load;
			int looked_at=0;
			while (mapq_heap.load>0 && (looked_at<mqi->top_hits || mqi->top_hits==0)) {
				heap_pa_elem tmp;
				heap_pa_extract_max(&mapq_heap,&tmp);
				pa=tmp.rest;
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
			memcpy(mapq_heap.array,hpa.array,sizeof(heap_pa_elem)*length);
			mapq_heap.load=hpa.load;
			looked_at=0;
			while (mapq_heap.load>0 && (looked_at<mqi->top_hits || mqi->top_hits==0)) {
				heap_pa_elem tmp;
				heap_pa_extract_max(&mapq_heap,&tmp);
				pa=tmp.rest;
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
			memcpy(mapq_heap.array,hpa.array,sizeof(heap_pa_elem)*length);
			mapq_heap.load=hpa.load;
			looked_at=0;
			while (mapq_heap.load>0 && (looked_at<mqi->top_hits || mqi->top_hits==0)) {
				heap_pa_elem tmp;
				heap_pa_extract_max(&mapq_heap,&tmp);
				pa=tmp.rest;
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
		//Lets find out what to print and set the SAM strings accordingly
		int printed_alignments=0;
		while (hpa.load>0 && ( printed_alignments<of->number_outputs || of->number_outputs==0) ) {
			heap_pa_elem tmp;
			heap_pa_extract_max(&hpa,&tmp);
			pa=tmp.rest;
			if (!of->strata || (tmp.score==top_score && tmp.isize==best_isize)) {
				if (pa->paired_sequencing && !of->half_paired && (!pa->mapped || !pa->mate_pair->mapped)) {
					continue;
				}
				pa->mark=true; 
				pretty_print_sam_update(pa);
				if (pa->mate_pair!=NULL) {
					pa->mate_pair->mark=true; 
					pretty_print_sam_update(pa->mate_pair);
				}
				printed_alignments++;
				if (!has_mapping) {
					break;
				}
			} else {
				break;
			}
		}
	}
	heap_pa_destroy(&hpa);
	//clean up all alignments that will not be printed
	pretty * new_head=NULL;
	pretty * new_tail=NULL;
	pa=head;
	while (pa!=NULL) {
		pretty * tmp=pa->next;
		if (pa->mark) {
			if (new_head==NULL) {
				new_head=pa;
			} else {
				new_tail->next=pa;
			}
			new_tail=pa;
			pa->next=NULL;
		} else {
			if (pa->mate_pair!=NULL) {
				pretty_free(pa->mate_pair);
			}
			pretty_free(pa);
		}
		pa=tmp;
	}
	*ppa=new_head;	
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
static int get_sam_entries_read(file_buffer * fb, pretty ** pas,int * last_used) {
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
					if (strncmp(read_names[j],pa->read_name,length)==0 && (!isdigit(read_names[j][length]) || read_names[j][length-1]=='F' || read_names[j][length-1]=='R')) {
						if (pas[j]==NULL) {
							pas[j]=pa;
						} else {
							pa->next=pas[j];
							pas[j]=pa;
						}
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

static int get_sam_entries(file_buffer * fb, pretty ** pas) {
	int last_used=0;
	while (fb->eof==0) {
		if (get_sam_entries_read(fb,pas,&last_used)==0) {
			return 0;
		}
	}
	if (fb->eof==1) {
		if (get_sam_entries_read(fb,pas,&last_used)==0) {
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
				read_names[read_names_filled++]=read_names[line++]+1;
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
					read_names_filled--;
					reads_file_buffer.left_over=reads_file_buffer.filled - (read_names[read_names_filled]-1-reads_file_buffer.buffer);
					return 0;
				} else {
					//got this read keep going
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
				read_names[read_names_filled++]=read_names[line++]+1;
				while (line<index && read_names[line][0]!='>' && read_names[line][0]!='\n') { line++; };
				if (line>=index && reads_file_buffer.eof!=1) {
					//didn't get all of last read
					read_names_filled--;
					reads_file_buffer.left_over=reads_file_buffer.filled - (read_names[read_names_filled]-1-reads_file_buffer.buffer);
					return 0;
				} else {
					//got this read keep going
					in_header=1;	
				}
			}
		}
	}	
	return 0;
}

static int find_last_read(bool fastq_reads,bool colour_space)  {
	index_new_read_names(fastq_reads,colour_space);
	int64_t i; int new_index=0;
	//fprintf(stderr,"have %d names filled \n",read_names_filled);
	if (read_names_filled==0) {
		return 0;
	}
	for (i=1; i<read_names_filled; i++) {
		if (strcmp(read_names[new_index],read_names[i])!=0) {
			read_names[++new_index]=read_names[i];
		}
	}
	read_names_filled=new_index+1;
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
		qsort(sam_headers,headers_found,sizeof(char*),qsort_sam_header_strcmp);
		fprintf(output_file,"%s\n",sam_headers[0]);
		for (i=1; i<headers_found; i++) {
			if (strcmp(sam_headers[i-1],sam_headers[i])!=0) {
				fprintf(output_file,"%s\n",sam_headers[i]);
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
	file_pa = (pretty***)malloc(sizeof(pretty**)*number_of_files);
	if (file_pa==NULL) {
		fprintf(stderr,"Failed to allocate memory for file_pa\n");
		exit(1);
	}
	memset(file_pa,0,sizeof(pretty*)*number_of_files);
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
		#pragma omp parallel for schedule(static,1) num_threads(read_threads)
		for (i=0; i<number_of_files; i++) {
			//fprintf(stderr,"%d: with file %d\n",omp_get_thread_num(),i);
			file_pa[i]=(pretty**)malloc(sizeof(pretty*)*read_names_filled);
			if (file_pa[i]==NULL) {
				fprintf(stderr,"Failed to allocate memory for file_pa[i]\n");
				exit(1);
			}
			memset(file_pa[i],0,sizeof(pretty*)*read_names_filled);
			get_sam_entries(&file_buffers[i],file_pa[i]);	
		}
		master_pa=(pretty**)malloc(sizeof(pretty*)*read_names_filled);
		if (master_pa==NULL) {
			fprintf(stderr,"Failed to allocate memory for master_pa\n");
			exit(1);
		}
		memset(master_pa,0,sizeof(pretty*)*read_names_filled);
		for (i=0; i<read_names_filled; i++) {
			int j;
			for (j=0; j<number_of_files; j++) {
				if (file_pa[j][i]!=NULL) {
					pretty * pa = file_pa[j][i];
					if (master_pa[i]==NULL) {
						master_pa[i]=pa;
					} else {
						pretty * tail = pa;
						while(tail->next!=NULL) { tail=tail->next; };
						tail->next=master_pa[i];
						master_pa[i]=pa;
					}
				}	
			}
		}
		#pragma omp parallel for num_threads(compute_threads) schedule(static,50000)
		for (i=0; i<read_names_filled; i++) {
			pretty ** pa = master_pa+i;
			process_read(pa,&of,&mqi);
		}
		for (i=0; i<read_names_filled; i++) {
			pretty * pa = master_pa[i];
			while (pa!=NULL) {
				fprintf(output_file,"%s",pa->sam_string);
				if (pa->mate_pair!=NULL) {
					fprintf(output_file,"%s",pa->mate_pair->sam_string);
					pretty_free(pa->mate_pair);
				}
				pretty * tmp = pa->next; 
				pretty_free(pa);
				pa=tmp;
			}
		}
		for (i=0; i<number_of_files; i++) {
			free(file_pa[i]);
			file_pa[i]=NULL;
		}
		free(master_pa);
	}	
	for (i=0; i<number_of_files; i++) {
		free(file_buffers[i].buffer);
		gzclose(file_buffers[i].file);	
	}
	fprintf(stderr,"Processed %lld reads, %d\n",reads_processed,reads_file_buffer.eof);
	free(file_pa);
	free(file_buffers);	
	gzclose(reads_file_buffer.file);
	free(reads_file_buffer.buffer);
	free(read_names);
	return;
}
