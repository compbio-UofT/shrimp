#include <stdio.h>
#include <assert.h>
#include <string.h>
#include <stdlib.h>
#include "fasta_reader.h"


void fasta_move(file_buffer * fb, char * dest) {
        size_t unseen_end_mod=fb->unseen_end%fb->size;
        size_t unseen_start_mod=fb->unseen_start%fb->size;
	size_t space=(unseen_start_mod>unseen_end_mod ? unseen_end_mod : fb->size ) - unseen_end_mod;
}

void parse_fasta(fasta_reader * fard) {
	/*assert(fard->reads_seen+fard->reads_unseen==fard->reads_filled);
	char * current_newline=NULL;
	while (fard->fb->unseen_end!=fard->fb->unseen_start && fard->reads_filled<fard->reads_inmem) {
        	size_t unseen_end_mod=fard->fb->unseen_end%fard->fb->size;
        	size_t unseen_start_mod=fard->fb->unseen_start%fard->fb->size;
		size_t space=(unseen_start_mod>=unseen_end_mod ? fard->fb->size : unseen_end_mod) - unseen_start_mod;
		current_newline=(char*)memchr(fard->fb->base+unseen_start_mod,'\n',space);
		if (current_newline==NULL) {
			fard->fb->unseen_start+=space;
		} else {
			*current_newline='\0';
			size_t current_newline_index=current_newline-fard->fb->base;
			if (current_newline_index==0 || fard->fb->base[current_newline_index-1]=='\0') {
				current_newline_index++;
				for(; current_newline_index<fard->fb->size && fard->fb->base[current_newline_index]=='\0'; current_newline_index++) {fard->fb->base[current_newline_index]='\0';};
				fard->fb->unseen_start+=current_newline_index-unseen_start_mod;
				unseen_start_mod=fard->fb->unseen_start%fard->fb->size;
			} else {
				char * read_name=fard->fb->base+unseen_start_mod;
					if (read_name[0]=='>') {
						size_t read_name_length=strlen(read_name);
						memmove(fard->read_names+fard->reads_filled*SIZE_READ_NAME,read_name+1,read_name_length);
						fard->read_names[fard->reads_filled*SIZE_READ_NAME+read_name_length]='\0';
						char * space=strchr(fard->read_names+fard->reads_filled*SIZE_READ_NAME,' ');
						if (space!=NULL) {*space='\0';};
						char * tab=strchr(fard->read_names+fard->reads_filled*SIZE_READ_NAME,'\t');
						if (tab!=NULL) {*tab='\0';};
						fard->reads_unseen++;
						fard->reads_filled++;
					}
				fard->fb->unseen_start+=current_newline_index-unseen_start_mod+1;
				unseen_start_mod=fard->fb->unseen_start%fard->fb->size;
			}
			fard->fb->exhausted=false;
		}
	}
	if (fard->fb->exhausted || fard->reads_filled==fard->reads_inmem) {	
		fard->reads_exhausted=true;
	} else {
		fard->reads_exhausted=false;
	}*/
}
