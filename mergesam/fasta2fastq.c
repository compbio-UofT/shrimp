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

void to_ascii_33(char * s) {
	char * prev=s;
	char * current=s;
	int index=0;
	while (isspace(*current) && *current!='\0') { current++; };
	while (*current!='\0') {
		if (*current==' ') {
			*current='\0';
			int temp=atoi(prev);
			s[index++]=(char)(34+(char)temp);
			current++;
			prev=current;
			while (isspace(*current) && *current!='\0') { current++; };
		} else {
			current++;
		}
		
	}
	//last one
	if (*prev!='\0') {
		int temp=atoi(prev);
		s[index++]=(char)(34+(char)temp);
	}
	s[index]='\0';
	return;
}

static void  inline update_last_line(lineindex_table * li, file_buffer * fb, int lines_processed ) {
	//size_t lines_processed=li->end-li->start;
	char * qual_last_line=li->table[(li->start+lines_processed-1)%li->size];
	//fprintf(stderr,"last line |%s|\n",qual_last_line);
	//fprintf(stderr,"last %lu\n",qual_last_line+strlen(qual_last_line)-fb->base);
	size_t qual_last_line_mod=qual_last_line-fb->base;
	size_t qual_start_mod=fb->unseen_start%fb->size;
	if (qual_last_line_mod >= qual_start_mod) {
		fb->unseen_start+=qual_last_line_mod-qual_start_mod;
	} else {
		fb->unseen_start+=fb->size-(qual_start_mod-qual_last_line_mod);
	}
	while(fb->base[fb->unseen_start%fb->size]!='\0' && fb->unseen_start<fb->unseen_end) {
		fb->unseen_start++;
	}
	while (fb->base[fb->unseen_start%fb->size]=='\0' && fb->unseen_start<fb->unseen_end) {
		fb->unseen_start++;
	}
	return;
}

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
		{"colour-space",0,0,3},
		{"letter-space",0,0,4},
                {"help", 0, 0, 5},
		{"buffer-size", 1, 0, 6},
		{"read-size", 1, 0, 7},
                {0,0,0,0}
        };

static inline bool fill_fb(file_buffer * fb) {
	bool has_changed=false;
	while (!fb->exhausted) {
		fill_read_buffer(&fb->frb);
		add_read_buffer_to_main(fb);
		if (!fb->exhausted && !fb->changed && fb->frb.eof==0) {
			fprintf(stderr,"too small buffer!\n");
			exit(1);
		}
		has_changed=has_changed || fb->changed;
	}
	return has_changed;
	//fprintf(stdout,"Filled %lu to %lu of %lu |%s|\n",fb->unseen_start, fb->unseen_end, fb->size,fb->base);
}

static void inline fill_fb_and_index(lineindex_table * li, lineindex_table ** thread_lineindexes, file_buffer * fb) {
	size_t old_em=fb->unseen_end%fb->size;	
	if (fill_fb(fb)) {
		size_t newly_added;
		size_t em=fb->unseen_end%fb->size;
		if (em > old_em) {
			newly_added=add_lineindex_from_memory_threaded(li, thread_lineindexes,fb->base+old_em, em-old_em,options.threads, '#');
		} else {
			newly_added=add_lineindex_from_memory_threaded(li, thread_lineindexes,fb->base+old_em, fb->size-old_em,options.threads, '#');
			newly_added+=add_lineindex_from_memory_threaded(li, thread_lineindexes,fb->base, em,options.threads, '#');
		}	
		if (newly_added>0) {
			fb->exhausted=false;
		}
	}
	return;
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
		case 5:
			usage(argv[0]);
			break;
		case 'N':
			options.threads=atoi(optarg);
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

	lineindex_table * qual_thread_lineindexes[options.threads];
	int i;
	for (i=0; i<options.threads; i++) {
		qual_thread_lineindexes[i]=lineindex_init(1);
	}
	//master table
	lineindex_table * qual_li = lineindex_init(1);

	lineindex_table * fasta_thread_lineindexes[options.threads];
	for (i=0; i<options.threads; i++) {
		fasta_thread_lineindexes[i]=lineindex_init(1);
	}
	//master table
	lineindex_table * fasta_li = lineindex_init(1);

	size_t reads_processed=0;
	//get the hit list, process it, do it again!
	fprintf(stderr,"Setting up buffer with size %lu and read_size %lu\n",options.buffer_size,options.read_size);
	file_buffer * qual_fb = fb_open(qual_filename,options.buffer_size,options.read_size); 
	file_buffer * fasta_fb = fb_open(fasta_filename,options.buffer_size,options.read_size); 
	int lines_processed=-1;
	while (lines_processed!=0) {
		//index lines
		fill_fb_and_index(qual_li, qual_thread_lineindexes,qual_fb);
		fill_fb_and_index(fasta_li, fasta_thread_lineindexes,fasta_fb);
		lines_processed=qual_li->end-qual_li->start;
		//print lines
		for (i=0; i<lines_processed; i++) {
			fprintf(stderr,"%s\n",fasta_li->table[(fasta_li->start+i)%fasta_li->size]);
			char * qual_string=qual_li->table[(qual_li->start+i)%qual_li->size];
			if (qual_string[0]!='>') {
				to_ascii_33(qual_li->table[(qual_li->start+i)%qual_li->size]);
				fprintf(stderr,"%s\n",qual_li->table[(qual_li->start+i)%qual_li->size]);
			}
		}
		update_last_line(qual_li,qual_fb,lines_processed);
		update_last_line(fasta_li,fasta_fb,lines_processed);
		//fprintf(stderr,"%lu %lu\n",qual_fb->unseen_start,qual_fb->unseen_end);
		//fprintf(stderr,"END OF IT |%s|\n",qual_fb->base+qual_fb->unseen_start%qual_fb->size);
		qual_li->start+=lines_processed;
		fasta_li->start+=lines_processed;
	}
	return 0;
}



