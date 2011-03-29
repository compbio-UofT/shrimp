#include <stdio.h>
#include <assert.h>
#include <string.h>
#include <stdlib.h>
#include "file_buffer.h"
#include "fasta_reader.h"

void fasta_last_entry(char * base,size_t used, char ** start,char ** end) {
	char * ptr=(char*)memrchr(base,'>',used);
	if (ptr==NULL) {
		fprintf(stderr,"please make buffer larger!\n");
		exit(1);
	}
	*start=ptr;
	ptr=(char*)memchr(ptr,'\n',used-(ptr-base+1));
	if (ptr==NULL) {
		fprintf(stderr,"Are you sure this is fasta and qual format?\n");	
		exit(1);
	}
	*end=ptr;
	**end='\0';
	fprintf(stderr,"last read |%s|\n",*start);
	return;
}

size_t fasta_move(file_buffer * fb, char * dest) {
        size_t unseen_end_mod=fb->unseen_end%fb->size;
        size_t unseen_start_mod=fb->unseen_start%fb->size;
	size_t used=0;
	if (unseen_start_mod >= unseen_end_mod) {
		size_t to_move=fb->size-unseen_start_mod;
		void * ptr=memmove(dest+used,fb->base+unseen_start_mod,to_move);
		if (ptr==NULL) {
			fprintf(stderr,"failed to move fasta_move!\n");
			exit(1);
		}
		used+=to_move;
		to_move=unseen_end_mod;
		ptr=memmove(dest+used,fb->base,to_move);
		if (ptr==NULL) {
			fprintf(stderr,"failed to move fasta_move!\n");
			exit(1);
		}
		used+=to_move;
	} else {
		size_t to_move=unseen_end_mod - unseen_start_mod;
		void * ptr=memmove(dest+used,fb->base+unseen_start_mod,to_move);
		if (ptr==NULL) {
			fprintf(stderr,"failed to move fasta_move!\n");
			exit(1);
		}
		used+=to_move;
	}
	if (fb->frb.eof!=0 && fb->frb.exhausted) {
		fprintf(stderr,"returning from fasta_move because hit EOF\n");
		return used;	
	}
	//find out where to put the end
	char * ptr=(char*)memrchr(dest,'>',used);
	if (ptr==NULL) {
		fprintf(stderr,"fasta_move failed because not enough memory to fit one fasta entry in the buffer!\n");
		exit(1);
	}
	used=ptr-dest; //how many bytes are part of the full fasta entry , also subtract the '>'
	fb->unseen_start+=used;
	return used;
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
