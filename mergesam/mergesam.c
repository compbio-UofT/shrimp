#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>
#include <assert.h>
#include "sam2pretty_lib.h"
#include "file_buffer.h"
#include "mergesam.h"
#include "mergesam_heap.h"
#include "fastx_readnames.h"
#include "sam_reader.h"
#include "../common/util.h"
#include <ctype.h>
#include <omp.h>

runtime_options options;


fastx_readnames fxrn;

//Pretty linked list
pp_ll * sam_headers;
pp_ll ** pp_ll_index;
pp_ll * master_ll;

//SAM file ios
heap_pa * thread_heaps;
typedef struct sam_reader sam_reader;
sam_reader ** sam_files;
char * sam_header_filename=NULL;


void process_sam_headers() {
	if (found_sam_headers && sam_header_filename==NULL) {
		int header_entries=0;
		int i;
		for (i=0; i<options.number_of_sam_files; i++) {
			header_entries+=sam_files[i]->sam_headers->length;
		}
		//get pointers to each header line
		char ** sam_lines=(char**)malloc(sizeof(char*)*header_entries);
		if (sam_lines==NULL) {
			fprintf(stderr,"Failed to allocate memory for sam_header entries!\n");
			exit(1);
		}	
		int index=0;
		for (i=0; i<options.number_of_sam_files; i++) {
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
		memset(sam_headers,0,sizeof(pp_ll)*options.number_of_sam_files);	
		found_sam_headers=false;
	}
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
//	fprintf(stderr,
//	"      --read-size      Read size, read into buffer with this (Default: %d)\n",DEF_READ_SIZE);
	fprintf(stderr, 
	"   -s/--stack-size     Input alignment stack size            (Default: %d)\n",DEF_ALIGNMENTS_STACK_SIZE);
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
	fprintf(stderr,
	"      --unaligned-fastx\n");
	fprintf(stderr,
	"      --aligned-fastx\n");
	exit(1);
}


struct option long_op[] =
        {
		{"fastq", 0, 0, 'Q'},
		{"sam",0,0,'E'},
		{"threads",1,0,'N'},
		{"report",1,0,'o'},
		{"max-alignments",1,0,1},
		{"sam-header",1,0,2},
		{"colour-space",0,0,3},
		{"letter-space",0,0,4},
                {"help", 0, 0, 5},
		{"buffer-size", 1, 0, 6},
		{"read-size", 1, 0, 7},
		{"read-rate",1,0,8},
		{"expected-isize",1,0,9},
		{"half-paired",0,0,10},
		{"sam-unaligned",0,0,11},
		{"strata",0,0,12},
		{"stack-size",1,0,'s'},
		{"unaligned-fastx",0,0,'u'},
		{"aligned-fastx",0,0,'a'},
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
	options.max_alignments=DEF_MAX_ALIGNMENTS;
	options.max_outputs=DEF_MAX_OUTPUTS;
	options.expected_insert_size=DEF_INSERT_SIZE;	
	options.fastq=false;
	options.strata=false;
	options.half_paired=false;
	options.sam_unaligned=false;
	options.sam_format=false;
	options.paired=false;
	options.unpaired=false;
	options.colour_space=false;
	options.letter_space=false;
	options.buffer_size=DEF_BUFFER_SIZE;
	options.read_size=DEF_READ_SIZE;
	options.read_rate=DEF_READ_RATE;
	options.alignments_stack_size=DEF_ALIGNMENTS_STACK_SIZE;
	options.threads=1;
	found_sam_headers=false;
        int op_id;
        char short_op[] = "o:QN:Es:au";
        char c = getopt_long(argc, argv, short_op, long_op, &op_id);
        while (c != EOF) {
		switch (c) {
		case 'Q':
			options.fastq=true;
			break;
		case 6:
			options.buffer_size=string_to_byte_size(optarg);
			break;
		case 7:
			options.read_size=string_to_byte_size(optarg);
			break;
		case 8:
			options.read_rate=atol(optarg);
			break;
		case 9:
			options.expected_insert_size=atoi(optarg);
			if (options.expected_insert_size<0) {
				fprintf(stderr,"Please specify a insert size >= 0!\n");
				usage(argv[0]);
			}
			break;
		case 'N':
			options.threads=atoi(optarg);
			break;
		case 10:
			options.half_paired=true;
			break;
		case 11:
			options.sam_unaligned=true;
			break;
		case 12:
			options.strata=true;
			break;
		case 'o':
			options.max_outputs=atoi(optarg);
			if (options.max_outputs<=0) {
				fprintf(stderr,"Please specify a max_output that is positive!\n");
				usage(argv[0]);
			}
			break;
		case 1:
			options.max_alignments=atoi(optarg);
			if (options.max_alignments<=0) {
				fprintf(stderr,"Please specify a max_alignments that is positive!\n");
				usage(argv[0]);
			}
			break;
		case 5:
			usage(argv[0]);
			break;
		case 2:
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
			options.sam_format=true;
			break;
		case 3:
			options.colour_space=true;
			break;
		case 4:
			options.letter_space=true;
			break;
		case 's':
			options.alignments_stack_size=atoi(optarg);
			assert(options.alignments_stack_size>0);
			break;
		case 'a':
			options.aligned_fastx=true;
			break;
		case 'u':
			options.unaligned_fastx=true;
			break;
		default:
			fprintf(stderr,"%d : %c , %d is not an option!\n",c,(char)c,op_id);
			usage(argv[0]);
			break;
		}
        	c = getopt_long(argc, argv, short_op, long_op, &op_id);
	}


	if ((options.colour_space && options.letter_space) || (!options.colour_space && !options.letter_space)) {
		fprintf(stderr,"can only enter either --colour-space or --letter-space not both!\n");
		usage(argv[0]);
	}
	if (!options.sam_format) {
		fprintf(stderr,"%s currently only supports sam format, please use '--sam' or '-E'\n",argv[0]);
		usage(argv[0]);
	}

	omp_set_num_threads(options.threads); 
	fprintf(stderr,"Set to %d threads!\n",options.threads);
		
	
	if (argc<=optind+1) {
		fprintf(stderr,"Please specify reads file and at least one sam file!\n");
		usage(argv[0]);
	}

	fprintf(stderr,"Size of %lu\n",sizeof(pretty));
	memset(&fxrn,0,sizeof(fastx_readnames));
	fxrn.reads_inmem=20*options.read_rate;
	fxrn.read_names=(char*)malloc(sizeof(char)*fxrn.reads_inmem*SIZE_READ_NAME);
	if (fxrn.read_names==NULL) {
		fprintf(stderr,"failed to allocate memory for read_names\n");
		exit(1);
	}
	//memset(fxrn.read_names,'Z',sizeof(char)*fxrn.reads_inmem*SIZE_READ_NAME);


	argc-=optind;
	argv+=optind;
	
	//Variables for IO of read names
	char * reads_filename=argv[0];
	fprintf(stderr,"Using %s as reads filename\n",reads_filename);
	argc--;
	argv++;

	options.number_of_sam_files=argc;
	//Open each sam input file	
	sam_files=(sam_reader**)malloc(sizeof(sam_reader*)*options.number_of_sam_files);
	if (sam_files==NULL) {
		fprintf(stderr,"failed to allocate memory for sam_files!\n");
		exit(1);
	}

	master_ll = (pp_ll*)malloc(sizeof(pp_ll)*options.read_rate);
	if (master_ll==NULL) {
		fprintf(stderr,"Failed to allocate memory for master_ll\n");
		exit(1);
	}
	memset(master_ll,0,sizeof(pp_ll)*options.read_rate);


	//allocate memory for sam_headers
	sam_headers=(pp_ll*)malloc(sizeof(pp_ll)*options.number_of_sam_files);
	if (sam_headers==NULL) {
		fprintf(stderr,"Failed to allocate memory for sam_headers\n");
		exit(1);
	}
	

	//index first the read then the file number	
	pp_ll_index = (pp_ll**)malloc(sizeof(pp_ll*)*options.read_rate*options.number_of_sam_files);
	if (pp_ll_index==NULL) {
		fprintf(stderr,"Failed to allocate memory for pp_ll_index!\n");
		exit(1);
	}	
	int i;
	for (i=0; i<options.number_of_sam_files; i++) {
		sam_files[i]=sam_open(argv[i],&fxrn);
		sam_files[i]->sam_headers=sam_headers+i;
		sam_files[i]->fileno=i;
		int j; 
		for (j=0; j<options.read_rate; j++) {
			pp_ll_index[j*options.number_of_sam_files+i]=sam_files[i]->pp_lls+j*LL_ALL;
		}
	}

	//calculate alignments cutoff
	int32_t alignments_cutoff=options.max_alignments==-1 ? options.max_outputs : MIN(options.max_alignments,options.max_outputs);


	//max the heaps for each thread
	fprintf(stderr,"There are %d threads\n",options.threads);
	thread_heaps=(heap_pa* )malloc(sizeof(heap_pa)*options.threads);
	if (thread_heaps==NULL) {
		fprintf(stderr,"Failed to allocate memory for thread_heaps!\n");
		exit(1);
	}
	for (i=0; i<options.threads; i++ ) {
		heap_pa_init(thread_heaps+i,alignments_cutoff+1);
		fprintf(stderr,"INIT THREAD_HEAP %p\n",thread_heaps+i);
	}	

	size_t reads_processed=0;
	clock_t start_time=clock();
	long iterations=0;
	//get the hit list, process it, do it again!
	fprintf(stderr,"Setting up buffer with size %lu and read_size %lu\n",options.buffer_size,options.read_size);
	bool have_non_eof_file=true;
	fxrn.fb=fb_open(reads_filename,options.buffer_size,options.read_size);
	while (!fxrn.reads_exhausted && have_non_eof_file) {
		//fprintf(stderr,"LOOP START\n");
		//fprintf(stderr,"Reads seen %lu, reads unseen %lu, reads filled %lu\n",fxrn.reads_seen,fxrn.reads_unseen,fxrn.reads_filled);
		//Populate the hitlist to as large as possible
		while (!fxrn.reads_exhausted) {
			fill_fb(fxrn.fb);
			parse_reads(&fxrn);	
		}
		//fprintf(stderr,"YReads seen %lu, reads unseen %lu, reads filled %lu\n",fxrn.reads_seen,fxrn.reads_unseen,fxrn.reads_filled);
		//clear the sam_headers memory
		memset(sam_headers,0,sizeof(pp_ll)*options.number_of_sam_files);	
	
		//Hit list in memory start processing!
		
		int i;	
		clock_t last_time = clock ();
		while (fxrn.reads_seen<fxrn.reads_filled) {
			iterations++;
			int reads_to_process=options.read_rate;

	
			//READ IN DATA
			#pragma omp parallel for schedule(guided)
			for (i=0; i<options.number_of_sam_files; i++) {
				if (sam_files[i]->fb->frb.eof!=1 || sam_files[i]->fb->unseen_end!=sam_files[i]->fb->unseen_inter) {
					fill_fb(sam_files[i]->fb);	
					parse_sam(sam_files[i],&fxrn);
				}
			}
			
			//find the minimum number of complete read alignments in memory
			have_non_eof_file=false;
			for (i=0; i<options.number_of_sam_files; i++) {
				if (sam_files[i]->fb->frb.eof!=1 || sam_files[i]->fb->unseen_end!=sam_files[i]->fb->unseen_inter) {
					reads_to_process=MIN(reads_to_process,sam_files[i]->last_tested-fxrn.reads_seen);
					have_non_eof_file=true;
				}
			}
			//fprintf(stderr,"Processing %d reads entries on this interation..\n",reads_to_process);
			if (reads_to_process>0) {
				for (i=0; i<options.number_of_sam_files; i++) {
					const size_t read_id=(fxrn.reads_seen+reads_to_process-1)%options.read_rate;
					//fprintf(stderr,"read rate %lu, using read_id %lu %lu\n",options.read_rate,read_id,sam_files[i]->inter_offsets[read_id]);
					sam_files[i]->fb->unseen_start=sam_files[i]->inter_offsets[read_id];
					sam_files[i]->pretty_stack_start=sam_files[i]->pretty_stack_ends[read_id];
				}
			} else {
				fprintf(stderr,"AN ERROR HAS OCCURED!\n");
				exit(1);	
			}
	
			assert(reads_to_process<=options.read_rate);

			if (options.paired && options.unpaired) {
				fprintf(stderr,"FAIL! can't have both paired and unpaired data in input file!\n");
				exit(1);
			}
			
			//prepare output
			#pragma omp parallel for schedule(guided) //schedule(static,10)
			for (i=0; i<reads_to_process; i++) {	
				const size_t read_id=(fxrn.reads_seen+i)%options.read_rate;
				pp_ll_combine_and_check(master_ll+i,pp_ll_index+read_id*options.number_of_sam_files,thread_heaps+omp_get_thread_num());
				int j;
				for (j=0; j<options.number_of_sam_files; j++) {
					memset(pp_ll_index[read_id*options.number_of_sam_files+j],0,sizeof(pp_ll)*LL_ALL);
				}
			}

			//output
			process_sam_headers();	
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
			memset(master_ll,0,sizeof(pp_ll)*reads_to_process);
			//update counters
			if (reads_to_process>0) {
				fxrn.reads_exhausted=false;
				fxrn.reads_seen+=reads_to_process;
				fxrn.reads_unseen-=reads_to_process;
			}
			reads_processed+=reads_to_process;
			//fprintf(stderr,"XReads seen %lu, reads unseen %lu, reads filled %lu\n",fxrn.reads_seen,fxrn.reads_unseen,fxrn.reads_filled);
			if ( (clock()-last_time)/options.threads > CLOCKS_PER_SEC/4) {
				double reads_per_second=reads_processed/( (double)(clock()-start_time)/(CLOCKS_PER_SEC*options.threads));
				double reads_per_iteration=reads_processed/(double)iterations;
				fprintf(stderr,"Processing overall at %lf reads / second, %lf reads / iteration, processed %lu\n",reads_per_second,reads_per_iteration,reads_processed);
				last_time=clock();
			} 
		}
	}
	fprintf(stderr,"Processed %lu reads\n",reads_processed);
	free(master_ll);
	free(sam_headers);
	free(pp_ll_index);
	for (i=0; i<options.number_of_sam_files; i++) {
		sam_close(sam_files[i]);
	}
	free(sam_files);
	for (i=0; i<options.threads; i++) {
		heap_pa_destroy(thread_heaps+i);
	}
	free(thread_heaps);
	fb_close(fxrn.fb);
	free(fxrn.read_names);
	return 0;
}




