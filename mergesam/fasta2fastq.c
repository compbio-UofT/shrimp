#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>
#include <assert.h>
#include "file_buffer.h"
#include "fasta_reader.h"
#include "fasta2fastq.h"
#include <ctype.h>
#include <omp.h>

runtime_options options;

fasta_reader fard;
char * threads_fasta;
char * threads_qual;

void usage(char * s) {
	fprintf(stderr, 
	"usage: %s [options/parameters] <fasta> <qual>\n", s);
	fprintf(stderr,
	"   <fasta> The fasta filename to read in\n");
	fprintf(stderr,
	"   <qual>  The qual file corresponding to the fasta file\n");
	fprintf(stderr,
	"Parameters:      (all sizes are in bytes unless specified)\n");
	fprintf(stderr,
	"      --buffer-size    File buffer size in memory per file   (Default: %d)\n",DEF_BUFFER_SIZE);
	fprintf(stderr,
	"      --read-size      Read size, read into buffer with this (Default: %d)\n",DEF_READ_SIZE);
	fprintf(stderr,
	"   -N/--threads        The number of threads to use          (Default: 1)\n");
	fprintf(stderr,"\nOptions:\n");
	fprintf(stderr,
	"      --help           This usage screen\n");
	fprintf(stderr,
	"      --colour-space   Reads file contains ABSOLiD data\n");
	fprintf(stderr,
	"      --letter-space   Reads file contains non-ABSOLiD data\n");	
	exit(1);
}


struct option long_op[] =
        {
		{"threads",1,0,'N'},
		{"colour-space",0,0,3},
		{"letter-space",0,0,4},
                {"help", 0, 0, 5},
		{"buffer-size", 1, 0, 6},
		{"read-size", 1, 0, 7},
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
	options.threads=1;
        int op_id;
        char short_op[] = "N:";
        char c = getopt_long(argc, argv, short_op, long_op, &op_id);
        while (c != EOF) {
		switch (c) {
		case 6:
			options.buffer_size=string_to_byte_size(optarg);
			break;
		case 7:
			options.read_size=string_to_byte_size(optarg);
			break;
		case 'N':
			options.threads=atoi(optarg);
			break;
		case 5:
			usage(argv[0]);
			break;
		case 3:
			options.colour_space=true;
			break;
		case 4:
			options.letter_space=true;
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

	omp_set_num_threads(options.threads); 
	fprintf(stderr,"Set to %d threads!\n",options.threads);
	
	if (argc<=optind+1) {
		fprintf(stderr,"Please specify reads file and at least one sam file!\n");
		usage(argv[0]);
	}

	argc-=optind;
	argv+=optind;

	if (argc!=2) {
		fprintf(stderr,"Please specify both a fasta and qual file!\n");
		usage(argv[0]);
	}

	
	//Variables for IO of read names
	char * fasta_filename=argv[0];
	char * qual_filename=argv[1];
	fprintf(stderr,"Using %s as fasta reads filename and %s as qual filename\n",fasta_filename,qual_filename);
	argc-=2;
	argv+=2;


	//max the heaps for each thread
	fprintf(stderr,"There are %d threads\n",options.threads);
	threads_fasta=(char* )malloc(sizeof(char)*options.buffer_size*options.threads);
	if (threads_fasta==NULL) {
		fprintf(stderr,"Failed to allocate memory for threads_fasta!\n");
		exit(1);
	}
	threads_qual=(char* )malloc(sizeof(char)*options.buffer_size*options.threads);
	if (threads_qual==NULL) {
		fprintf(stderr,"Failed to allocate memory for threads_qual!\n");
		exit(1);
	}
	

	size_t reads_processed=0;
	//get the hit list, process it, do it again!
	fprintf(stderr,"Setting up buffer with size %lu and read_size %lu\n",options.buffer_size,options.read_size);
	fard.fb=fb_open(fasta_filename,options.buffer_size,options.read_size);
	while (!fard.exhausted) {
		//Populate the hitlist to as large as possible
		fill_fb(fard.fb);
		fprintf(stderr,"|%s|\n",fard.fb->base);
		exit(1);
			
		//Should be done processing
		//keep unseen reads, and do it again
	}
	fb_close(fard.fb);
	return 0;
}




