#include <stdio.h>
#include <ctype.h>
#include <assert.h>
#include <string.h>
#include <stdlib.h>
#include <getopt.h>
#include <omp.h>
#include "sam2pretty_lib.h"
#include "merge_sam.h"


#define DEF_MAX_ISIZE 50000

void usage(char * s , int e) {
	FILE * f = stderr;
	fprintf(f,"mergesam [%s] - usage\n",s);
	fprintf(f,"-- Required -----------------------------\n");
	fprintf(f,"%s -r reads_filename [--colour-space --letter-space] file1.sam [file2.sam file3.sam ...]\n",s);
	fprintf(f,"  input reads file and SAM files can be gzipped or not\n");
	fprintf(f,"-----------------------------------------\n");
	fprintf(f,"-r/--reads           Filename of input reads file *\n");
	fprintf(f,"-Q/--fastq           Reads are in fastq format\n");
	fprintf(f,"   --colour-space\n");
	fprintf(f,"   --letter-space\n");
	fprintf(f,"-- Optional -----------------------------\n");
	fprintf(f,"-o/--output          Output filename                   Default: stdout\n");
	fprintf(f,"-O/--report          Report this many hits per read    Default: All\n");
	fprintf(f,"-X/--max-alignments  Skip read if more than this hits  Default: All\n");
	fprintf(f,"   --strata          Only consider highest score hits  Default: No\n");
	fprintf(f,"   --sam-unaligned   Output unaligned entires          Default: No\n");
	fprintf(f,"   --half-paired     Output halfpaired entries         Default: No\n");
	fprintf(f,"   --mapq            Compute MAPQ values               Default: No\n");
	fprintf(f,"   --mapq-hits       MAPQ based on top this many hits  Default: (report)+1\n");
	fprintf(f,"   --mapq-strata     MAPQ based on strata hits         Default: No\n");
	fprintf(f,"-i/--insert-size     Insert size for score tie-break   Default: Off\n");
        fprintf(f,"-d/--detect-isize    Detect the insert size and exit   Default: Off\n");
	fprintf(f,"   --read-isize      Read isize from file              Default: Off\n");
	fprintf(f,"   --max-isize       When detecting use this as max    Default: %d\n",DEF_MAX_ISIZE);
	fprintf(f,"   --sam-header      Use this file for SAM header      Default: Off\n");
	fprintf(f,"-- Runtime Settings----------------------\n");
	fprintf(f,"-R/--read-threads    Should be less then #SAM files    Default: 1\n");
	fprintf(f,"-C/--compute-threads Should be low, try 1,2,3          Default: 1\n");
	fprintf(f,"-M/--buffer-size     Buffer size in MB **              Default: 10\n");
	fprintf(f,"-- Program Information ------------------\n");
	fprintf(f,"-h/--help            This help screen\n");
	fprintf(f,"-v/--version         Print the version and exit\n");
	fprintf(f,"-----------------------------------------\n");
	fprintf(f,"*  This can be gzipped or not. It is assumed the format is FASTA unless\n");
	fprintf(f,"   -Q is used in the command line. When in paired mode, only one of the\n");
	fprintf(f,"   input files is required to be passed.\n");
	fprintf(f,"** This sets how big the input file buffers should be, most memory is\n");
	fprintf(f,"   is used in storing the internal SAM representation of alignments.\n");
	fprintf(f,"   Therefore changing this parameter by a bit may dramatically effect\n");
	fprintf(f,"   memory usage.\n");
	if (e==1) {
		exit(1);
	}
}


struct option long_op[] =
	{   
		{"help", 0, 0, 'h'},
		{"version", 0, 0, 'v'},
		{"reads",1,0,'r'},
		{"threads",1,0,'N'},
		{"report",1,0,'O'},
		{"max-alignments",1,0,'X'},
		{"strata",0,0,3},
		{"sam-unaligned",0,0,4},
		{"half-paired",0,0,5},
		{"fastq",0,0,'Q'},
		{"mapq",0,0,6},
		{"mapq-hits",1,0,7},
		{"mapq-strata",0,0,8},
		{"colour-space",0,0,9},
		{"letter-space",0,0,10},
		{"read-threads",1,0,'R'},
		{"compute-threads",1,0,'C'},
		{"buffer-size",1,0,'M'},
		{"output",1,0,'o'},
		{"insert-size",1,0,'i'},
		{"detect-isize",0,0,'d'},
		{"max-isize",1,0,11},
		{"read-isize",1,0,12},
		{"sam-header",1,0,13},
		{0,0,0,0}
	};
		


int main(int argc, char** argv) {
	char * reads_filename=NULL;
	bool reads_fastq=false;
	bool colour_space=false;
	bool letter_space=false;
	char * output_filename="stdout";
	FILE * output_file=stdout;
	output_filter of;
	memset(&of,0,sizeof(output_filter));
	of.max_isize=DEF_MAX_ISIZE; 
	mapq_info mqi;
	mqi.score_matrix_lambda = 0.1824502;
	mqi.sw_scale_constant = 0.5;
	mqi.calculate=false;
	mqi.top_hits=-2; //-2 means unset, -1 means same as of, 0 means all
	mqi.strata=false;
	int read_threads=1;
	int compute_threads=1;
	int buffer_size=10;
	int op_id;
	char short_op[] = "r:hvN:O:X:Qi:C:R:M:o:d";
	char c = getopt_long(argc, argv, short_op, long_op, &op_id);
	while (c != EOF) {
		switch (c) {
		case 'h':
			usage(argv[0],1);
			break;
		case 'v':
			fprintf(stderr,"mergesam part of SHRiMP 2.1\n");
			break;
		case 'r':
			reads_filename=optarg;
			break;
		case 'N':
			{
				int threads=atoi(optarg);
				omp_set_num_threads(threads);
			}
			break;
		case 'i':
			of.isize=atoi(optarg);
			of.use_isize=true;
			break;
		case 'd':
			of.detect_isize=true;
			break;
		case 3:
			of.strata=true;
			break;
		case 4:
			of.unaligned=true;
			break;
		case 5:
			of.half_paired=true;
			break;
		case 'X':
			of.max_alignments=atoi(optarg);
			break;
		case 'O':
			of.number_outputs=atoi(optarg);
			break;
		case 'Q':
			reads_fastq=true;
			break;
		case 6:
			mqi.calculate=true;
			break;
		case 7:
			mqi.top_hits=atoi(optarg);
			break;
		case 8:
			mqi.strata=true;
			break;			
		case 9:
			colour_space=true;
			break;
		case 10:
			letter_space=true;
			break;
		case 11:
			of.max_isize=atoi(optarg);
			break;
		case 'C':
			compute_threads=atoi(optarg);
			break;
		case 'R':
			read_threads=atoi(optarg);
			break;	
		case 'M':
			buffer_size=atoi(optarg);
			break;
		case 'o':
			output_filename=optarg;
			output_file=fopen(optarg,"w");
			if (output_file==NULL) {
				fprintf(stderr,"There has been an error opening output file %s\n",output_filename);
				perror("");
				exit(1);
			}	
			break;
		case 12:
			{
			FILE * fptr = fopen(optarg,"r");
			if (fptr==NULL) {
				fprintf(stderr,"Failed to open file %s for reading!\n",optarg);
				perror("");
				exit(1);
			}
			double mean;
			int ret=fscanf(fptr,"Mean: %lf,",&mean);
			if (ret!=1) {
				fprintf(stderr,"Failed to read mean from file %s!\n",optarg);
				exit(1);
			}
			of.isize=(int)mean;
			fclose(fptr);
			}	
			break;
		case 13:
			{
			char * sam_header_filename=optarg;
                        FILE * sam_header_file = fopen(sam_header_filename,"r");
                        if (sam_header_file==NULL) {
                                perror("Failed to open sam header file ");
                                exit(1);
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
			of.header_provided=true;
			break;
		default:
			fprintf(stderr,"invalid argument!\n");
			usage(argv[0],1);
			break;
		}
		c = getopt_long(argc, argv, short_op, long_op, &op_id);
	}
	if (reads_filename==NULL) {
		fprintf(stderr,"Please specify an input reads filename\n");
		usage(argv[0],1);
	}
	if (argc==optind) {
		fprintf(stderr,"Please specify at least one sam input file!\n");
		usage(argv[0],1);
	}
	if (colour_space && letter_space) {
		fprintf(stderr,"Please specify whether the input reads are colour space or letter space\n");
		usage(argv[0],1);
	}
	if (!colour_space && !letter_space) {
		fprintf(stderr,"Please specify whether the input reads are colour space or letter space\n");
		usage(argv[0],1);
	}
	//if (compute_threads>1) {
	//	fprintf(stderr,"NEED TO FIX CODE, currently dosen't work, TODO :  fix the print order \n");
	//	exit(1);
	//}
	//print stats to user
	if (mqi.top_hits==-2) {
		if (of.number_outputs==0) {
			mqi.top_hits=0;
		} else {
			mqi.top_hits=of.number_outputs+1;
		}
	}	
	fprintf(stderr,"Input reads file: %s, in %s format\n",reads_filename,reads_fastq ? "FASTQ" : "FASTA");
	fprintf(stderr,"Output file: %s\n",output_filename);
	//print out output filter
	fprintf(stderr,"Output Filter:\n\tReport:\t%d\tMax-alignments:\t%d\n\tStrata:\t%s\tSam-unaligned:\t%s\tHalf_paired:\t%s\n\tUse isize: %s\tisize: %d\n",
		of.number_outputs, of.max_alignments, of.strata ? "Yes" : "No", of.unaligned ? "Yes" : "No" , of.half_paired ? "Yes" : "No",of.use_isize ?  "Yes" : "No",of.isize );
	//print out mapq info
	fprintf(stderr,"MAPQ Info:\n\tCompute MAPQ:\t%s\n\tScoreMatrixLambda:\t%e\tSW-ScaleConstant:\t%e\n\t#TopHits:\t%d\tStrata:\t%s\n",
		mqi.calculate ? "Yes" : "No" , mqi.score_matrix_lambda, mqi.sw_scale_constant, mqi.top_hits, mqi.strata ? "Yes" : "No");
	argc-=optind;
	argv+=optind;
	merge_sam(reads_filename,reads_fastq,colour_space,read_threads,compute_threads,argc,argv,output_file,of,mqi,buffer_size);
	return 0;
}


