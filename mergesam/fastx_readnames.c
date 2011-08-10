#include <stdio.h>
#include <assert.h>
#include <string.h>
#include <stdlib.h>
//#include "mergesam.h"
//#include "mergesam_heap.h"
#include "../common/util.h"
#include "fastx_readnames.h"
void reconsolidate_reads(fastx_readnames * fxrn) {
	assert(fxrn->reads_seen+fxrn->reads_unseen==fxrn->reads_filled);
	if (fxrn->reads_unseen!=fxrn->reads_filled) {
		memmove(fxrn->read_names,fxrn->read_names+fxrn->reads_seen*SIZE_READ_NAME,(fxrn->reads_filled-fxrn->reads_seen)*SIZE_READ_NAME);
		fxrn->reads_filled=fxrn->reads_unseen;
		fxrn->reads_seen=0;
		fxrn->reads_exhausted=false;
	}
}


void parse_reads(fastx_readnames * fxrn) {
	assert(fxrn->reads_seen+fxrn->reads_unseen==fxrn->reads_filled);
	//reconsolidate_reads(fxrn);
	char * current_newline=NULL;
	while (fxrn->fb->unseen_end!=fxrn->fb->unseen_start && (fxrn->reads_filled-fxrn->reads_seen)<fxrn->reads_inmem) {
        	size_t unseen_end_mod=fxrn->fb->unseen_end%fxrn->fb->size;
        	size_t unseen_start_mod=fxrn->fb->unseen_start%fxrn->fb->size;
		size_t space=(unseen_start_mod>=unseen_end_mod ? fxrn->fb->size : unseen_end_mod) - unseen_start_mod;
		current_newline=(char*)memchr(fxrn->fb->base+unseen_start_mod,'\n',space);
		if (current_newline==NULL) {
			fxrn->fb->unseen_start+=space;
		} else {
			*current_newline='\0';
			size_t current_newline_index=current_newline-fxrn->fb->base;
			if (current_newline_index==0 || fxrn->fb->base[current_newline_index-1]=='\0') {
				current_newline_index++;
				for(; current_newline_index<fxrn->fb->size && fxrn->fb->base[current_newline_index]=='\0'; current_newline_index++) {fxrn->fb->base[current_newline_index]='\0';};
				fxrn->fb->unseen_start+=current_newline_index-unseen_start_mod;
				unseen_start_mod=fxrn->fb->unseen_start%fxrn->fb->size;
			} else {
				char * read_name=fxrn->fb->base+unseen_start_mod;
				if (options.fastq) {
					if (!fxrn->fastq_seen_name) {
						if (read_name[0]=='@') {
							const size_t slot=fxrn->reads_filled%fxrn->reads_inmem;
							const size_t read_name_length=MIN(SIZE_READ_NAME-1,strlen(read_name));
							memmove(fxrn->read_names+slot*SIZE_READ_NAME,read_name+1,read_name_length);
							fxrn->read_names[slot*SIZE_READ_NAME+read_name_length]='\0';
							char * const  space=strchr(fxrn->read_names+slot*SIZE_READ_NAME,' ');
							if (space!=NULL) {*space='\0';};
							char * const tab=strchr(fxrn->read_names+slot*SIZE_READ_NAME,'\t');
							if (tab!=NULL) {*tab='\0';};
							fxrn->reads_unseen++;
							fxrn->reads_filled++;
							fxrn->fastq_seen_name=true;
							//fxrn->fastq_seen_seq=options.colour_space ? -1 : 0;
							fxrn->fastq_seen_seq=0;
							fxrn->fastq_seen_qual=0;
							fxrn->fastq_seen_plus=false;
						}
					} else {
						if (!fxrn->fastq_seen_plus) {
							if (read_name[0]=='+') {
								fxrn->fastq_seen_plus=true;
							} else {
								if (options.mode_set==false) {
									if (read_name[0]!='\0' && read_name[1]!='\0') {
										if (read_name[1]<63) {
											fprintf(stderr,"Detected colour-space!\n");
											options.colour_space=true;
										} else {
											fprintf(stderr,"Detected letter-space!\n");
											options.colour_space=false;
										}
										options.mode_set=true;

									}
								}
								fxrn->fastq_seen_seq+=strlen(read_name)-1;
							}
						} else {
							fxrn->fastq_seen_qual+=strlen(read_name)-1;
							assert(fxrn->fastq_seen_qual<=fxrn->fastq_seen_seq+(options.colour_space ? -1 : 0));
							assert(options.mode_set);
							if (fxrn->fastq_seen_qual==fxrn->fastq_seen_seq+(options.colour_space ? -1 : 0)) {
								fxrn->fastq_seen_name=false;
							}	
						}
					}
				} else {
					if (read_name[0]=='>') {
						const size_t slot=fxrn->reads_filled%fxrn->reads_inmem;
						const size_t read_name_length=MIN(SIZE_READ_NAME-1,strlen(read_name));
						memmove(fxrn->read_names+slot*SIZE_READ_NAME,read_name+1,read_name_length);
						fxrn->read_names[slot*SIZE_READ_NAME+read_name_length]='\0';
						char * space=strchr(fxrn->read_names+slot*SIZE_READ_NAME,' ');
						if (space!=NULL) {*space='\0';};
						char * tab=strchr(fxrn->read_names+slot*SIZE_READ_NAME,'\t');
						if (tab!=NULL) {*tab='\0';};
						fxrn->reads_unseen++;
						fxrn->reads_filled++;
					}
				}
				fxrn->fb->unseen_start+=current_newline_index-unseen_start_mod+1;
				unseen_start_mod=fxrn->fb->unseen_start%fxrn->fb->size;
			}
			fxrn->fb->exhausted=false;
		}
	}
	if (fxrn->fb->exhausted || fxrn->reads_filled==fxrn->reads_inmem) {	
		fxrn->reads_exhausted=true;
	} else {
		fxrn->reads_exhausted=false;
	}
}
