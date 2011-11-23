#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>
#include <assert.h>
#include <time.h>
#include "file_buffer.h"
#include "fasta2fastq.h"
#include "lineindex_lib.h"
#include <ctype.h>
#include <omp.h>

runtime_options options;


typedef struct {
	size_t size;
	size_t used;
	char * base;
} thread_buffer;

thread_buffer * thread_buffers;

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
	//while(fb->base[fb->unseen_start%fb->size]=='\0' && fb->unseen_start<fb->unseen_end) {
	//	fb->unseen_start++;
	//}
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
	"Required:\n");
	fprintf(stderr,
	"      --qv-offset      The ASCII offset for the integer values in the qual file\n");
	fprintf(stderr,
	"Parameters:      (all sizes are in bytes unless specified)\n");
	fprintf(stderr,
	"      --buffer-size    File buffer size in memory per file   (Default: %d)\n",DEF_BUFFER_SIZE);
	fprintf(stderr,
	"      --read-size      Read size, read into buffer with this (Default: %d)\n",DEF_READ_SIZE);
	fprintf(stderr,"\nOptions:\n");
	fprintf(stderr,
	"      --help           This usage screen\n");
	exit(1);
}


struct option long_op[] =
        {
                {"help", 0, 0, 5},
		{"buffer-size", 1, 0, 6},
		{"read-size", 1, 0, 7},
		{"qv-offset",1,0,8},
                {0,0,0,0}
        };

static inline bool fill_fb(file_buffer * fb) {
	time_t io_start_time=time(NULL);
	fprintf(stderr,"IO start ... ");
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
	fprintf(stderr,"IO end ... %lu seconds\n",(time(NULL)-io_start_time));
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
	options.buffer_size=DEF_BUFFER_SIZE;
	options.read_size=DEF_READ_SIZE;
	options.threads=1;
	options.qv_offset=DEF_QV_OFFSET;
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
		case 8:
			options.qv_offset=atoi(optarg);
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

	if (options.qv_offset<=0) {
		fprintf(stderr,"Please specify a qv_offset. This is used when converting qual files into fastq format.\nFor SOLiD data this value will be most likely 34.\nFor Illumina data this value will be most likely 64, except for Illumina 1.8+ when it is 33.\n");		
		usage(argv[0]);
		exit(1);
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


	//set up the thread_buffers
	thread_buffers=(thread_buffer*)malloc(sizeof(thread_buffer)*options.threads);
	if (thread_buffers==NULL) {
		fprintf(stderr,"Failed to malloc memory for thread buffers!\n");
		exit(1);
	}
	for (i=0; i<options.threads; i++ ) {
		thread_buffers[i].size=(options.buffer_size/options.threads+1000)*1.3;
		thread_buffers[i].base=(char*)malloc(sizeof(char)*thread_buffers[i].size);
		if (thread_buffers[i].base==NULL) {
			fprintf(stderr,"Failed to allocate memory for thread buffers!\n");
			exit(1);
		}
	}

	//get the hit list, process it, do it again!
	fprintf(stderr,"Setting up buffer with size %lu and read_size %lu\n",options.buffer_size,options.read_size);
	file_buffer * qual_fb = fb_open(qual_filename,options.buffer_size,options.read_size); 
	file_buffer * fasta_fb = fb_open(fasta_filename,options.buffer_size,options.read_size); 
	size_t lines_processed=0;
	clock_t start_time=clock();
	clock_t last_time=clock();
	size_t iterations=0;
	bool first_loop=true;
	while (lines_processed!=0 || first_loop) {
		first_loop=false;
		//index lines
		fill_fb_and_index(qual_li, qual_thread_lineindexes,qual_fb);
		fill_fb_and_index(fasta_li, fasta_thread_lineindexes,fasta_fb);
		lines_processed=qual_li->end-qual_li->start;
		//print lines
		//figure out which thread handles which
		int lines_to_print[options.threads];
		int start[options.threads];
		for (i=0; i<options.threads; i++) {
			start[i]=(i==0 ?  0 : start[i-1]+lines_to_print[i-1]);
			lines_to_print[i]=lines_processed/options.threads+(lines_processed%options.threads > i ? 1 : 0);
			//fprintf(stderr,"Thread %d , start %d , lines %d\n",i,start[i],lines_to_print[i]);
		}
		#pragma omp parallel  
		{
		int thread_id = omp_get_thread_num();
		thread_buffer * ob = thread_buffers+thread_id;
		ob->used=0;
		ob->base[0]='\0';
		int i;
		for (i=start[thread_id]; i<start[thread_id]+lines_to_print[thread_id]; i++) {
			//fprintf(stderr,"Running %d on %d, %d\n",thread_id,i,(fasta_li->start+i)%fasta_li->size);
			//fprintf(stdout,"@%s\n",fasta_li->table[(fasta_li->start+i)%fasta_li->size]+1);
			char * to_print=fasta_li->table[(fasta_li->start+i)%fasta_li->size];
			char * qual_string=qual_li->table[(qual_li->start+i)%qual_li->size];
			//fprintf(stderr,"F |%s| vs |%s| \n",to_print,qual_string);
			while (strlen(to_print)+strlen(qual_string)+ob->used>ob->size) {
				ob->size*=1.3;
				ob->base=(char*)realloc(ob->base,sizeof(char)*ob->size);
				if (ob->base==NULL) {
					fprintf(stderr,"Failed to allocate memory for thread_buffer expand\n");
					exit(1);
				}
			}
			//qual string moves between read names and quals, check if we are printing a read name or qual
			if (qual_string[0]=='>') {
				ob->used+=sprintf(ob->base+ob->used,"@%s\n",fasta_li->table[(fasta_li->start+i)%fasta_li->size]+1);
			} else {
				to_ascii_33(qual_li->table[(qual_li->start+i)%qual_li->size]);
				ob->used+=sprintf(ob->base+ob->used,"%s\n",fasta_li->table[(fasta_li->start+i)%fasta_li->size]);
				ob->used+=sprintf(ob->base+ob->used,"+\n%s\n",qual_li->table[(qual_li->start+i)%qual_li->size]);
			}
			//fprintf(stderr,"%s and %s\n",to_print,qual_string);
			//assert(qual_string[0]!='>');
		}
		}
		
		//print the blocks
		for (i=0; i<options.threads; i++) {
			fprintf(stdout,"%s",thread_buffers[i].base);
		}
		
		update_last_line(qual_li,qual_fb,lines_processed);
		update_last_line(fasta_li,fasta_fb,lines_processed);
		//fprintf(stderr,"%lu %lu\n",qual_fb->unseen_start,qual_fb->unseen_end);
		//fprintf(stderr,"END OF IT |%s|\n",qual_fb->base+qual_fb->unseen_start%qual_fb->size);
		qual_li->start+=lines_processed;
		fasta_li->start+=lines_processed;
		if (lines_processed>0) {
			qual_fb->exhausted=false;
			fasta_fb->exhausted=false;
		}
		iterations++;
		if ( (clock()-last_time)/options.threads > CLOCKS_PER_SEC/4) {
			double lines_per_second=qual_li->start/( (double)(clock()-start_time)/(CLOCKS_PER_SEC*options.threads));
			double lines_per_iteration=qual_li->start/(double)iterations;
			fprintf(stderr,"Processing overall at %lf reads / second, %lf reads / iteration, processed %lu, lines on this iteration %lu\n",lines_per_second,lines_per_iteration,qual_li->start,lines_processed);
			last_time=clock();
		}

	}
	//free the line-indexes
	for (i=0; i<options.threads; i++) {
		lineindex_destroy(qual_thread_lineindexes[i]);
		lineindex_destroy(fasta_thread_lineindexes[i]);
	}
	//free the master index
	lineindex_destroy(qual_li);
	lineindex_destroy(fasta_li);
	//close the file_buffers
	fb_close(qual_fb);
	fb_close(fasta_fb);
	return 0;
}




