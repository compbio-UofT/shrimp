#define _MODULE_GMAPPER
#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>
#include <assert.h>
#include <ctype.h>
#include <omp.h>
#include "mergesam.h"
#include "sam2pretty_lib.h"
#include "file_buffer.h"
#include "fastx_readnames.h"
#include "sam_reader.h"
#include "../common/util.h"
#include "../gmapper/gmapper-defaults.h"
#include "../gmapper/gmapper.h"

#include "mergesam_heap.h"

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
//char * sam_header_filename=NULL;
//sam_header_filename=NULL;

int64_t genome_length=0;


char * command_line=NULL;


//char * reads_filename;
//reads_filename=NULL;


int64_t genome_length_from_headers(char ** sam_lines, int header_entries) {
	int64_t g_length=0;
	int i; 
	for (i=0; i<header_entries; i++) {
		char * line = sam_lines[i];
		if (line[0]=='@' && line[1]=='S' && line[2]=='Q') {
			int x;
			for (x=4; line[x]!='\0'; x++) {
				if (line[x]=='\t' && line[x+1]=='L' && line[x+2]=='N') {
					x++;
					break;
				}
			}
			if (line[x]=='\0' || line[x]!='L') {
				fprintf(stderr,"Invalid sam header format\n");
				exit(1);
			}	
			assert(line[x]=='L');
			x++;
			assert(line[x]!='\0'); //N
			x++;
			assert(line[x]!='\0'); //:
			x++;
			assert(line[x]!='\0');
			int y=x;
			while (line[y]!='\t' && line[y]!='\0' && line[y]!='\n') {
				y++;
			}
			char old_char = line[y];
			line[y]='\0';
			g_length+=atoi(line+x);
			line[y]=old_char;
		}
	}
	return g_length;
} 


void process_sam_headers() {
	if (found_sam_headers) {
		int header_entries=0;
		int i;
		for (i=0; i<options.number_of_sam_files; i++) {
			header_entries+=sam_files[i]->sam_headers->length;
		}
		//get pointers to each header line
		char ** sam_lines=(char**)malloc(sizeof(char*)*(header_entries+1)); //add in the command line
		if (sam_lines==NULL) {
			fprintf(stderr,"Failed to allocate memory for sam_header entries!\n");
			exit(1);
		}	
		int index=0;
		int pg_id=0;
		for (i=0; i<options.number_of_sam_files; i++) {
			pretty * pa = sam_files[i]->sam_headers->head;
			while(pa!=NULL) {
				char * s = pa->sam_string; 
				if (strncmp(s,"@PG	ID:",strlen("@PG	ID:"))==0) {
					assert(pg_id<100000000);
					char * x = (char*)malloc(sizeof(char)*(strlen(s)+13));
					sprintf(x,"%s%d-%s","@PG	ID:",pg_id++,s+strlen("@PG	ID:"));
					pa->sam_string=x;	
				}
				sam_lines[index++]=pa->sam_string;
				pa=pa->next;
			}
		}
		assert(command_line!=NULL);
		genome_length=genome_length_from_headers(sam_lines, header_entries);	
		fprintf(stderr,"Calculated genome length to be , %ld\n",genome_length);
		//sam_lines[index++]=command_line;
		//want to sort the headers here
		qsort(sam_lines, header_entries, sizeof(char*),sam_header_sort);
		//want to print the headers here
		assert(index>0);
		fprintf(stdout,"%s\n",sam_lines[0]);
		bool printed_pg_self=false;
		if (sam_header_filename==NULL) { 
		for (i=1; i<index; i++) {
			int ret=sam_lines[i-1]!=NULL ? strcmp(sam_lines[i],sam_lines[i-1]) : 1;
			if (!printed_pg_self && strncmp(sam_lines[i],"@PG",strlen("@PG"))==0) {
				fprintf(stdout,"%s\n",command_line);
				printed_pg_self=true;
			}
			if (ret!=0) {
				fprintf(stdout,"%s\n",sam_lines[i]);
			}	
			if (strncmp(sam_lines[i],"@PG	ID:",strlen("@PG	ID:"))==0) {
				free(sam_lines[i]);
				sam_lines[i]=NULL;
			}
		}	
		}
		if (!printed_pg_self) {
			fprintf(stdout,"%s\n",command_line);
			printed_pg_self=true;
		}
		//get the genome length!!!
		free(sam_lines);
		memset(sam_headers,0,sizeof(pp_ll)*options.number_of_sam_files);	
		found_sam_headers=false;
	}
}

void usage(char * s) {
	fprintf(stderr, 
	"usage: %s [options/parameters] <r> <s1> <s2> ...\n", s);
	fprintf(stderr,
	"   <r>     Reads filename, if paired then one of the two paired files\n");
	fprintf(stderr,
	"   <s?>    A SAM file input for mergigng\n");
	fprintf(stderr,
	"Runtime:      (all sizes are in bytes unless specified)\n");
	fprintf(stderr,
	"      --buffer-size    File buffer size in memory per file   (Default: %d)\n",DEF_BUFFER_SIZE);
	fprintf(stderr,
	"      --read-size      Read size, read into buffer with this (Default: %d)\n",DEF_READ_SIZE);
	fprintf(stderr, 
	"   -s/--stack-size     Input alignment stack size            (Default: %d)\n",DEF_ALIGNMENTS_STACK_SIZE);
	fprintf(stderr,
	"      --read-rate      How many reads to process at once     (Default: %d)\n",DEF_READ_RATE);
	fprintf(stderr,
	"Output options:\n");
	fprintf(stderr,
	"      --un                    Output unaligned FAST(A/Q) file       (Default: disabled)\n");
	fprintf(stderr,
	"      --al                    Output aligned FAST(A/Q) file         (Default: disabled)\n");
	fprintf(stderr,
	"      --sam-unaligned         Unaligned reads in SAM output         (Default: disabled)\n");
	fprintf(stderr,
	"   -o/--report                The maximum alignments to report      (Default: %d)\n",DEF_MAX_OUTPUTS);
	fprintf(stderr,
	"   -N/--threads               The number of threads to use          (Default: 1)\n");
	fprintf(stderr,
	"   -E/--sam                   Output in SAM format                  (Default: disabled)\n");
	fprintf(stderr,
	"   -Q/--fastq                 Reads are in fastq format             (Default: auto-detect)\n");
	fprintf(stderr,
	"      --strata                Print only the best scoring hits\n");
	fprintf(stderr,
	"      --max-alignments        Max. align. per read  (0=all)        (Default: %d)\n",DEF_MAX_ALIGNMENTS); 
	fprintf(stderr,
	"      --no-half-paired        Do not output half-paired alignments (Default: disabled)\n");
	fprintf(stderr,
	"      --insert-size-dist      Use mean,stdev for insert size dist. (Default: %d,%d)\n",DEF_INSERT_SIZE_MEAN,DEF_INSERT_SIZE_STDDEV); 
	fprintf(stderr,
	"      --single-best-mapping   See documentation, output onlt *best (Default: disabled)\n");
	fprintf(stderr,
	"      --min-mapq              Minimum mapping quality              (Default: 0)\n");
	fprintf(stderr,
	"      --all-contigs           Output given no more merging after   (Default: disabled)\n");
	fprintf(stderr,
	"      --half-paired           Output half-paired mappings          (Default: enabled)\n");
	fprintf(stderr,
	"      --no-mapping-qualities  Do not compute mapping qualities     (Default: disabled)\n");
	fprintf(stderr,
	"      --leave-mapq            Leave the mapq field as it was       (Default: disabled)\n");
	fprintf(stderr,
	"      --sam-header            Use file as SAM header\n");
	fprintf(stderr,
	"      --no-improper-mappings  Do not pair up half-paired mappings  (Default: disabled)\n");
	fprintf(stderr,
	"      --no-autodetect-input   Do not try to auto-detect FAST(A/Q)  (Default: disabled)\n");
	fprintf(stderr,
	"      --help                  This usage screen\n");
	exit(1);
}


struct option long_op[] =
        {
		{"un",1,0,10},
		{"al",1,0,11},
		{"sam-unaligned",0,0,12},
		{"report",1,0,'o'},
		{"threads",1,0,'N'},
		{"sam",0,0,'E'},
		{"fastq", 0, 0, 'Q'},
		{"strata",0,0,9},
		{"max-alignments",1,0,14},
		{"no-half-paired",0,0,19},
		{"insert-size-dist",1,0,25},
		{"single-best-mapping",0,0,37},
		{"all-contigs",0,0,38},
		{"half-paired",0,0,41},
		{"no-mapping-qualities",0,0,39},
		{"sam-header",1,0,2},
		{"no-improper-mappings",0,0,42},
		{"no-autodetect-input",0,0,48},
		
		{"leave-mapq-untouched",0,0,3},
                {"help", 0, 0, 5},
		{"buffer-size", 1, 0, 6},
		{"read-size", 1, 0, 7},
		{"read-rate",1,0,8},
		{"stack-size",1,0,'s'},
		{"min-mapq",1,0,4},
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
	int i;
	char pg_line_prefix[]="@PG	ID:mergesam	VN:2.2.0	CL:";
	size_t command_line_length=strlen(pg_line_prefix);
	for (i=0; i<argc; i++) {
		command_line_length+=strlen(argv[i])+2;
	}
	command_line=(char*)malloc(sizeof(char)*command_line_length);
	if (command_line==NULL) {
		fprintf(stderr,"Failed to allocate memory to store command line arguments\n");
		exit(1);
	}
	command_line[0]='\0';
	strcpy(command_line,pg_line_prefix);
	for (i=0; i<argc; i++) {
		strcat(command_line,argv[i]);
		strcat(command_line," ");
	}	
	options.unaligned_reads_file=NULL;
	options.aligned_reads_file=NULL;
	options.sam_unaligned=sam_unaligned;
	options.max_outputs=DEF_MAX_OUTPUTS;
	options.threads=DEF_NUM_THREADS;
	options.sam_format=Eflag;
	options.fastq=Qflag;
	options.fastq_set=autodetect_input ? false : true; //has the fastq flag been set, i.e. has it been set to fastq=false or fastq=true 
	options.strata=strata_flag;
	options.max_alignments=DEF_MAX_ALIGNMENTS;
	options.half_paired=half_paired; //for no half paired default setting
	options.single_best=single_best_mapping;
	options.insert_size_mean=DEF_INSERT_SIZE_MEAN;
	options.insert_size_stddev=DEF_INSERT_SIZE_STDDEV;
	options.all_contigs=all_contigs;
	/* see above */ //options.half_paired=half_paired; //for no half paired default setting
	options.no_mapping_qualities=compute_mapping_qualities ? false : true;
	sam_header_filename=NULL;
	options.no_improper_mappings=improper_mappings ? false : true;
	options.no_autodetect_input=autodetect_input ? false : true;

	options.leave_mapq=false;

	options.paired=false;
	options.unpaired=false;
	options.colour_space=false;
	options.letter_space=false;
	options.mode_set=false;
	options.buffer_size=DEF_BUFFER_SIZE;
	options.read_size=DEF_READ_SIZE;
	options.read_rate=DEF_READ_RATE;
	options.alignments_stack_size=DEF_ALIGNMENTS_STACK_SIZE;
	options.min_mapq=0;
	found_sam_headers=false;
        int op_id;
        char short_op[] = "o:QN:Es:au";
        char c = getopt_long(argc, argv, short_op, long_op, &op_id);
        while (c != EOF) {
		switch (c) {
		//fastq
		case 'Q':
			options.fastq=true;
			break;
		//no-auto-detect
		case 48:
			options.no_autodetect_input=true;
			break;
		//single-best
		case 37:
			options.single_best=true;
			break;
		//buffer-size
		case 6:
			options.buffer_size=string_to_byte_size(optarg);
			break;
		//read-size
		case 7:
			options.read_size=string_to_byte_size(optarg);
			break;
		//read-rate
		case 8:
			options.read_rate=atol(optarg);
			break;
		//threads
		case 'N':
			options.threads=atoi(optarg);
			break;
		//half-paired
		case 41:
			options.half_paired=true;
			break;
		//no-half-paired
		case 19:
			options.half_paired=false;	
			break;
		//sam-unaligned	
		case 12:
			options.sam_unaligned=true;
			break;
		//strata
		case 9:
			options.strata=true;
			break;
		//report
		case 'o':
			options.max_outputs=atoi(optarg);
			if (options.max_outputs<=0) {
				fprintf(stderr,"Please specify a max_output that is positive!\n");
				usage(argv[0]);
			}
			break;
		//max-alignments
		case 14:
			options.max_alignments=atoi(optarg);
			if (options.max_alignments<=0) {
				fprintf(stderr,"Please specify a max_alignments that is positive!\n");
				usage(argv[0]);
			}
			break;
		//help
		case 5:
			usage(argv[0]);
			break;
		//sam-header
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
		//sam format
		case 'E':
			options.sam_format=true;
			break;
		//stack-size
		case 's':
			options.alignments_stack_size=atoi(optarg);
			assert(options.alignments_stack_size>0);
			break;
		//al
		case 11:
			options.aligned_reads_file=fopen(optarg,"w");
			if (options.aligned_reads_file==NULL) {
				fprintf(stderr,"Failed to open file for writting %s\n",optarg);
				exit(1);
			}
			break;
		//un
		case 10:
			options.unaligned_reads_file=fopen(optarg,"w");
			if (options.unaligned_reads_file==NULL) {
				fprintf(stderr,"Failed to open file for writting %s\n",optarg);
				exit(1);
			}
			break;
		//insert-size-dist
		case 25:
			{
                        char * c = strtok(optarg, ",");
                        if (c == NULL) {
                          fprintf(stderr, "argmuent for insert-size-dist should be \"mean,stddev\" [%s]", optarg);
			  exit(1);
			}
                        options.insert_size_mean = atof(c);
                        c = strtok(NULL, ",");
                        if (c == NULL) {
                          fprintf(stderr, "argmuent for insert-size-dist should be \"mean,stddev\" [%s]", optarg);
			  exit(1);
			}
                        options.insert_size_stddev = atof(c);
			}
			break;
		//all-contigs
		case 38:
			options.all_contigs=true;
			break;
		//no-mapping-qualities
		case 39:
			options.no_mapping_qualities=true;
			break;
		//no-improper-mappings
		case 42:
			options.no_improper_mappings=true;
			break;	
		//leave-mapq
		case 3:
			options.leave_mapq=true;
			break;
		//min mapq
		case 4:
			options.min_mapq=atoi(optarg);
			break;
		default:
			fprintf(stderr,"%d : %c , %d is not an option!\n",c,(char)c,op_id);
			usage(argv[0]);
			break;
		}
        	c = getopt_long(argc, argv, short_op, long_op, &op_id);
	}


	/* SANITY CHECK ARGUMENTS */
	/*	{"un",1,0,10},
		{"al",1,0,11},
		{"sam-unaligned",0,0,12},
		{"report",1,0,'o'},
		{"threads",1,0,'N'},
		{"sam",0,0,'E'},
		{"fastq", 0, 0, 'Q'},
		{"strata",0,0,9},
		{"max-alignments",1,0,14},
		{"no-half-paired",0,0,19},
		{"insert-size-dist",1,0,25},
		{"single-best-mapping",0,0,37},
		{"all-contigs",0,0,38},
		{"half-paired",0,0,41},
		{"no-mapping-qualities",0,0,39},
		{"sam-header",1,0,2},
		{"no-improper-mappings",0,0,42},
		{"no-autodetect-input",0,0,48}, */
	if (options.unaligned_reads_file!=NULL && options.aligned_reads_file!=NULL) {
		fprintf(stderr," ! Please, '--un' xor '--al' == 1!\n");
		exit(1);
	}
	
	if (!options.sam_format && options.unaligned_reads_file==NULL && options.aligned_reads_file==NULL) {
		fprintf(stderr," ! Mergesam currently only supports output in SAM or FAST(A/Q) format, please use one of '--un','--al',or '--sam'\n");
		exit(1);
	}

	if (options.single_best && options.no_mapping_qualities) {
		fprintf(stderr," ! '--single-best' cannot be used in combination with '--no-mapping-qualities'\n");
		exit(1);
	}

	if (options.single_best) {
		options.max_outputs=1;
		fprintf(stderr," + Setting max outputs per class to 1, because of single_best.\n");	
	}

	/* END OF SANITY CHECKING ... EVERYTHING IS SANE .. MAYBE ... */

	omp_set_num_threads(options.threads); 
	fprintf(stderr," + Running with %d threads!\n",options.threads);
	
	if (argc<=optind+1) {
		fprintf(stderr," ! Please specify reads file and at least one sam file!\n");
		usage(argv[0]);
	}

	memset(&fxrn,0,sizeof(fastx_readnames));
	fxrn.reads_inmem=20*options.read_rate;
	fxrn.read_names=(char*)malloc(sizeof(char)*fxrn.reads_inmem*SIZE_READ_NAME);
	if (fxrn.read_names==NULL) {
		fprintf(stderr," ! Failed to allocate memory for read_names\n");
		exit(1);
	}
	//memset(fxrn.read_names,'Z',sizeof(char)*fxrn.reads_inmem*SIZE_READ_NAME);


	argc-=optind;
	argv+=optind;
	
	//Variables for IO of read names
	reads_filename=argv[0];
	fprintf(stderr," + Using %s as reads filename\n",reads_filename);
	argc--;
	argv++;

	options.number_of_sam_files=argc;
	//Open each sam input file	
	sam_files=(sam_reader**)malloc(sizeof(sam_reader*)*options.number_of_sam_files);
	if (sam_files==NULL) {
		fprintf(stderr," ! Failed to allocate memory for sam_files!\n");
		exit(1);
	}

	master_ll = (pp_ll*)malloc(sizeof(pp_ll)*options.read_rate);
	if (master_ll==NULL) {
		fprintf(stderr," ! Failed to allocate memory for master_ll\n");
		exit(1);
	}
	memset(master_ll,0,sizeof(pp_ll)*options.read_rate);


	//allocate memory for sam_headers
	sam_headers=(pp_ll*)malloc(sizeof(pp_ll)*options.number_of_sam_files);
	if (sam_headers==NULL) {
		fprintf(stderr," ! Failed to allocate memory for sam_headers\n");
		exit(1);
	}
	

	//index first the read then the file number	
	pp_ll_index = (pp_ll**)malloc(sizeof(pp_ll*)*options.read_rate*options.number_of_sam_files);
	if (pp_ll_index==NULL) {
		fprintf(stderr," ! Failed to allocate memory for pp_ll_index!\n");
		exit(1);
	}	
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
	int32_t alignments_cutoff=options.max_alignments==0 ? options.max_outputs : MIN(options.max_alignments,options.max_outputs);


	//max the heaps for each thread
	thread_heaps=(heap_pa* )malloc(sizeof(heap_pa)*options.threads);
	if (thread_heaps==NULL) {
		fprintf(stderr," ! Failed to allocate memory for thread_heaps!\n");
		exit(1);
	}
	for (i=0; i<options.threads; i++ ) {
		heap_pa_init(thread_heaps+i,alignments_cutoff+(options.single_best ? 0 : 1));
		fprintf(stderr," + Initializing thread_heap for thread %d at address %p\n",i,thread_heaps+i);
	}	

	//initialize the thread buffers for each thread
	output_buffer obs[options.threads];
	for (i=0; i<options.threads; i++) {
		obs[i].size=((size_t)(options.buffer_size*GROWTH_FACTOR))+1;
		obs[i].base=(char*)malloc(sizeof(char)*obs[i].size);
		if (obs[i].base==NULL) {
			fprintf(stderr," ! Failed to allocate memory for the output buffers!\n");
			exit(1);
		}
		obs[i].used=0;
	}

	size_t reads_processed=0;
	clock_t start_time=clock();
	long iterations=0;
	//get the hit list, process it, do it again!
	fprintf(stderr," + Setting up buffer with size %lu and read_size %lu\n",options.buffer_size,options.read_size);
	bool have_non_eof_file=true;
	fxrn.fb=fb_open(reads_filename,options.buffer_size,options.read_size);
	if (!options.no_autodetect_input && !options.fastq_set) {
		options.fastq=auto_detect_fastq(fxrn.fb->frb.file,reads_filename);
	}

	assert(options.unaligned_reads_file==NULL || options.aligned_reads_file==NULL);
	FILE * output_file=(options.unaligned_reads_file!=NULL ? options.unaligned_reads_file : (options.aligned_reads_file!=NULL ? options.aligned_reads_file : stdout));
	while (!fxrn.reads_exhausted && have_non_eof_file) {
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
				obs[omp_get_thread_num()].used=0;
				if (sam_files[i]->fb->frb.eof!=1 || sam_files[i]->fb->unseen_end!=sam_files[i]->fb->unseen_inter) {
					fill_fb(sam_files[i]->fb);	
					parse_sam(sam_files[i],&fxrn);
				}
			}
			process_sam_headers();	
			
			//find the minimum number of complete read alignments in memory
			have_non_eof_file=false;
			for (i=0; i<options.number_of_sam_files; i++) {
				if (sam_files[i]->fb->frb.eof!=1 || sam_files[i]->fb->unseen_end!=sam_files[i]->fb->unseen_inter) {
					//fprintf(stderr,"MIN(%d,%d)\n",reads_to_process,sam_files[i]->last_tested-fxrn.reads_seen);
					reads_to_process=MIN(reads_to_process,sam_files[i]->last_tested-fxrn.reads_seen);
					have_non_eof_file=true;
				}
			}
			if (!have_non_eof_file) {
				reads_to_process=fxrn.reads_filled-fxrn.reads_seen;
				for (i=0; i<options.number_of_sam_files; i++) {
					if (sam_files[i]->last_tested>fxrn.reads_seen) {
						reads_to_process=MAX(reads_to_process,sam_files[i]->last_tested-fxrn.reads_seen);
					}
				}
				if (reads_to_process==0) {
					break;
				}
				reads_to_process=MIN(reads_to_process,options.read_rate);
			}
			assert(reads_to_process<=options.read_rate);
			//fprintf(stderr,"Processing %d reads entries on this iteration..\n",reads_to_process);
			if (reads_to_process>0) {
				for (i=0; i<options.number_of_sam_files; i++) {
					const size_t read_id=(fxrn.reads_seen+reads_to_process-1)%options.read_rate;
					//fprintf(stderr,"read rate %lu, using read_id %lu %lu\n",options.read_rate,read_id,sam_files[i]->inter_offsets[read_id]);
					sam_files[i]->fb->unseen_start=sam_files[i]->inter_offsets[read_id];
					sam_files[i]->pretty_stack_start=sam_files[i]->pretty_stack_ends[read_id];
				}
			} else {
				fprintf(stderr,"AN ERROR HAS OCCURED! - try increasing buffer size?\n");
				exit(1);	
			}
	

			if (options.paired && options.unpaired) {
				fprintf(stderr,"FAIL! can't have both paired and unpaired data in input file!\n");
				exit(1);
			}
			//fprintf(stderr,"reads to process %lu\n",reads_to_process);
			//prepare output
			#pragma omp parallel for schedule(guided) //schedule(static,10)
			for (i=0; i<reads_to_process; i++) {	
				const size_t read_id=(fxrn.reads_seen+i)%options.read_rate;
				int thread_num = omp_get_thread_num();
				pp_ll_combine_and_check(master_ll+i,pp_ll_index+read_id*options.number_of_sam_files,thread_heaps+thread_num,obs+thread_num);
				int j;
				for (j=0; j<options.number_of_sam_files; j++) {
					memset(pp_ll_index[read_id*options.number_of_sam_files+j],0,sizeof(pp_ll)*LL_ALL);
				}
			}

			//output
			for (i=0; i<reads_to_process; i++) {
				pretty * pa=master_ll[i].head;
				while (pa!=NULL) {
					if (pa->paired_sequencing) {
						if (pa->first_in_pair) {
							fprintf(output_file,"%s\n",pa->sam_string);		
							if (pa->mate_pair!=NULL) {
								fprintf(output_file,"%s\n",pa->mate_pair->sam_string);		
							}
						} else {
							if (pa->mate_pair!=NULL) {
								fprintf(output_file,"%s\n",pa->mate_pair->sam_string);		
							}
							fprintf(output_file,"%s\n",pa->sam_string);		
						}
					} else {
						fprintf(output_file,"%s\n",pa->sam_string);		
					}
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
		free(obs[i].base);
	}
	free(thread_heaps);
	fb_close(fxrn.fb);
	free(fxrn.read_names);
	return 0;
}




