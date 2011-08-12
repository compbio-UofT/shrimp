#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>
#include <assert.h>
#include "file_buffer.h"
#include "fasta2fastq.h"
#include "lineindex_lib.h"
#include <ctype.h>
#include <omp.h>

runtime_options options;


void usage(char * s) {
	fprintf(stderr, 
	"usage: %s [options/parameters] <filename>\n", s);
	fprintf(stderr,
	"   <filename>\n");
	fprintf(stderr,
	"Parameters:      (all sizes are in bytes unless specified)\n");
	fprintf(stderr,
	"      --buffer-size    File buffer size in memory per file   (Default: %d)\n",DEF_BUFFER_SIZE);
	fprintf(stderr,
	"      --read-size      Read size, read into buffer with this (Default: %d)\n",DEF_READ_SIZE);
	exit(1);
}


struct option long_op[] =
        {
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
		case 5:
			usage(argv[0]);
			break;
		case 'N':
			options.threads=atoi(optarg);
			break;
		default:
			fprintf(stderr,"%d : %c , %d is not an option!\n",c,(char)c,op_id);
			usage(argv[0]);
			break;
		}
        	c = getopt_long(argc, argv, short_op, long_op, &op_id);
	}


	omp_set_num_threads(options.threads);
	fprintf(stderr,"Set to %d threads!\n",options.threads);
	
	argc-=optind;
	argv+=optind;

	if (argc!=1) {
		fprintf(stderr,"Please specify a filename!\n");
		usage(argv[0]);
	}

	
	//Variables for IO of read names
	char * filename=argv[0];
	fprintf(stderr,"Using %s as filename\n",filename);
	argc-=1;
	argv+=1;

	lineindex_table * thread_lineindexes[options.threads];
	int i;
	for (i=0; i<options.threads; i++) {
		thread_lineindexes[i]=lineindex_init(1);
	}
	//master table
	lineindex_table * li = lineindex_init(1);

	//get the hit list, process it, do it again!
	fprintf(stderr,"Setting up buffer with size %lu and read_size %lu\n",options.buffer_size,options.read_size);
	file_buffer * fb = fb_open(filename,options.buffer_size,options.read_size); 
	while (!fb->exhausted) {
		//fill buffer
		fill_fb(fb);
		//index lines
		size_t newly_added=add_lineindex_from_memory_threaded(li, thread_lineindexes,fb->base, fb->unseen_end-fb->unseen_start,options.threads, '#');
		if (newly_added>0) {
			fb->exhausted=false;
		}
		fb->unseen_start=0;
		fb->unseen_end=0;
		//print lines
		for (i=0; i<li->end; i++) {
			fprintf(stderr,"%d |%s|\n",i,li->table[i]);
		}
		li->start=0;
		li->end=0;
	}
	fb_close(fb);
	return 0;
}




