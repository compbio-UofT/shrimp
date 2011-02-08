
#include <zlib.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include "file_buffer.h"

#define MIN(a,b) (((a) < (b)) ? (a) : (b))


void fb_close(file_buffer * fb) {
	free(fb->base);
	free(fb->frb.base);
	gzclose(fb->frb.file);
	free(fb);
}

file_buffer * fb_open(char * path,size_t buffer_size,size_t read_size) {
	//allocate some memory
	file_buffer * fb = (file_buffer*)malloc(sizeof(file_buffer));
	if (fb==NULL) {
		fprintf(stderr,"file_buffer : failed to allocate base structure\n");
		exit(1);
	}
	memset(fb,0,sizeof(file_buffer));
	//allocate the main buffer
	fb->base=(char*)malloc(buffer_size);
	if (fb->base==NULL) {
		fprintf(stderr,"file_buffer : failed to allocate memory for base buffer, tried to allocate %lu bytes\n",buffer_size);
		exit(1);
	}
	fb->size=buffer_size;
	fb->unseen_start=0;
	fb->unseen_end=0;
	//allocate the read buffer
	fb->frb.base=(char*)malloc(read_size);
	if (fb->frb.base==NULL) {
		fprintf(stderr,"file_buffer : failed to allocate memory for base buffer, tried to allocate %lu bytes\n",read_size);
		exit(1);
	}
	fb->frb.size=read_size;
	fb->frb.filled=0;
	fb->frb.unseen=0;
	fb->frb.eof=0;
	//open the file
	fb->frb.file=(strcmp(path,"-")==0) ? gzdopen(fileno(stdin),"r") : gzopen(path,"r");
	if (fb->frb.file==NULL) {
		fprintf(stderr,"file_buffer : failed to open file %s\n",path);
		exit(1);
	}	
	return fb;
}

void fill_read_buffer(file_read_buffer * frb) {
	if (frb->unseen==frb->size) {
		return;
	}
	//move whatever was left over in buffer to front
	//fprintf(stderr,"unseen is %lu, seen is %lu\n",frb->unseen,frb->seen);
	//assert(frb->unseen==0 || !frb->pad);
	memmove(frb->base,frb->base+frb->filled-frb->unseen,frb->unseen);
	//read into the rest of the buffer
	int ret = gzread(frb->file,frb->base+frb->unseen,frb->size-frb->unseen);
	assert(frb->size-frb->unseen!=0);
	//fprintf(stderr,"trying to read %lu\n",frb->size-frb->unseen);
	if (ret<0) {
		fprintf(stderr,"A gzread error has occured\n");
		exit(1);
	}	
	frb->eof=gzeof(frb->file);
	//fprintf(stderr,"EOF %d ret %d\n",frb->eof,ret);
	if (ret==0 && frb->eof==0) {
		fprintf(stderr,"A error has occured in reading\n");
		exit(1);
	}
	
	frb->filled=ret+frb->unseen;
	frb->unseen=frb->filled;
	if (frb->eof==0 || frb->unseen!=0 ) {
		frb->exhausted=false;	
	} else {
		frb->exhausted=true;
	}
	frb->seen=0;
	//fprintf(stderr,"fill_read_buffer : %lu seen, %lu unseen, %lu filled, %lu size, EOF %d\n",frb->seen,frb->unseen,frb->filled, frb->size, frb->eof);
}


void mark_partial_line(file_read_buffer * frb,size_t offset) {
	assert(frb->seen+frb->unseen==frb->filled);
	//fprintf(stderr,"mark_partial_line : %lu seen, %lu unseen, %lu filled, %d eof\n",frb->seen,frb->unseen,frb->filled, frb->eof);
	//char buffer[129];
	//memcpy(buffer,frb->base,frb->filled);
	//buffer[frb->filled]='\0';
	//fprintf(stderr,"|%s|\n",buffer);
	frb->exhausted=offset>=frb->unseen;
	//if (frb->eof==1) {
	//	return;
	//}
	size_t unseen=frb->exhausted ? frb->unseen : offset;
	//fprintf(stderr,"mark_partial_line : unseen %lu, %lu offset, %lu unseen\n",unseen,offset,frb->unseen);
	char * ptr = (char*)memrchr(frb->base+frb->seen,'\n',unseen);
	if (ptr==NULL) {
		//fprintf(stderr,"mark_partial_line : HIT NULL\n");
		//frb->unseen=frb->filled;
	} else {
		size_t offset=(ptr-frb->base)+1;
		frb->unseen=frb->filled-offset;
		//fprintf(stderr,"mark_partial_line : HIT NOT NULL\n");
	}
}


//unseen_start is always 0 based actual location of start
//unseen_end is always 0 based location of next to be written data
void add_read_buffer_to_main(file_buffer * fb) {
	//fprintf(stderr,"ADDING\n");
	size_t unseen_end_mod=fb->unseen_end%fb->size;
	size_t unseen_start_mod=fb->unseen_start%fb->size;
	fb->changed=false;
	assert(fb->frb.seen+fb->frb.unseen==fb->frb.filled);
	assert(fb->unseen_start<=fb->unseen_end);
	assert(fb->unseen_start+fb->size>=fb->unseen_end);
	while (fb->unseen_end!=fb->unseen_start+fb->size && !fb->frb.exhausted) {
		//fprintf(stderr,"%lu unseen_end_mod, %lu unseen_startmod, %s exhausted\n",unseen_end_mod,unseen_start_mod,fb->frb.exhausted ? "YES" : "NO");
		//fprintf(stderr,"%lu unseen_end, %lu unseen_start, %s exhausted\n",fb->unseen_end,fb->unseen_start,fb->frb.exhausted ? "YES" : "NO");
		assert(fb->unseen_start+fb->size>=fb->unseen_end);
		assert(fb->frb.seen+fb->frb.unseen==fb->frb.filled);
		assert(fb->unseen_start<=fb->unseen_end);
		//fprintf(stderr,"%lu - %lu , %lu - %lu\n",unseen_start_mod,unseen_end_mod,fb->size,unseen_end_mod);
		size_t space=(unseen_start_mod>unseen_end_mod ? unseen_start_mod : fb->size ) - unseen_end_mod;
		//fprintf(stderr,"space %lu\n",space);
		assert(space>0);
		mark_partial_line(&fb->frb,space);
		size_t just_seen = fb->frb.filled-fb->frb.seen-fb->frb.unseen;
		//fprintf(stderr,"add_read_buffer_to_main : %lu just_seen, %lu frb.filled, %lu frb.seen, %lu frb.unseen\n",just_seen,fb->frb.filled,fb->frb.seen,fb->frb.unseen);
		if (just_seen!=0) {
			//fprintf(stderr,"COPY IN\n");
			//means we can copy some in
			memmove(fb->base+unseen_end_mod,fb->frb.base+fb->frb.seen,just_seen);
			fb->changed=true;
			fb->frb.seen+=just_seen;
			assert(fb->frb.seen+fb->frb.unseen==fb->frb.filled);
			fb->unseen_end+=just_seen;
			unseen_end_mod=fb->unseen_end%fb->size;
		} else if (!fb->frb.exhausted) {
			//fprintf(stderr,"add_read : 3\n");
			fb->exhausted=true;
		}
		//fprintf(stderr,"%lu unseen_end_mod, %lu unseen_startmod, %s exhausted\n",unseen_end_mod,unseen_start_mod,fb->frb.exhausted ? "YES" : "NO");
		//fprintf(stderr,"%lu unseen_end, %lu unseen_start, %s exhausted\n",fb->unseen_end,fb->unseen_start,fb->frb.exhausted ? "YES" : "NO");
		if (!fb->frb.exhausted) {
			if (unseen_start_mod<unseen_end_mod || fb->unseen_start==fb->unseen_end) { 
				//fprintf(stderr,"PADDING\n");
				assert(fb->size>=unseen_end_mod);
				memset(fb->base+unseen_end_mod,'\0',fb->size-unseen_end_mod);
				fb->changed=true;
				fb->unseen_end+=fb->size-unseen_end_mod;
				unseen_end_mod=fb->unseen_end%fb->size;
			} else {
				//fprintf(stderr,"add_read : 4\n");
				fb->exhausted=true;	
				return;
			}
		}
		unseen_end_mod=fb->unseen_end%fb->size;
	}
	if (fb->unseen_end==fb->unseen_start+fb->size || (fb->frb.exhausted && fb->frb.eof==1)) {
		//fprintf(stderr,"add_read : 5\n");
		fb->exhausted=true;
	}
	assert(fb->frb.seen+fb->frb.unseen==fb->frb.filled);
	assert(fb->unseen_start<=fb->unseen_end);
		assert(fb->unseen_start+fb->size>=fb->unseen_end);
	//fprintf(stderr,"condition %lu + %lu != %lu\n",fb->unseen_start,fb->size,fb->unseen_end);
}
	/*if (unseen_end_mod>unseen_start_mod) {
		if (unseen_end_mod==0) {
			assert(fb->frb.seen+fb->frb.unseen==fb->frb.filled);
			if (fb->frb.pad) {
				//fill top
				size_t space_top=unseen_start_mod;
				if (space_top>0) {
					mark_partial_line(&fb->frb,space_top);
					size_t just_seen = fb->frb.filled-fb->frb.seen-fb->frb.unseen;
					if (just_seen!=0) {
						//means we can copy some in
						memmove(fb->base,fb->frb.base+fb->frb.seen,just_seen);
						fb->frb.seen+=just_seen;
						assert(fb->frb.seen+fb->frb.unseen==fb->frb.filled);
						fb->unseen_end+=just_seen;
						unseen_end_mod=fb->unseen_end%fb->size;
					}
				}
			}
		}

	} else {
		//fill only middle
		size_t space_middle = (fb->unseen_end==fb->unseen_start) ? fb->size : unseen_start_mod-unseen_end_mod;
		if (space_middle!=0) {
			mark_partial_line(&fb->frb,space_middle);	
			size_t just_seen = fb->frb.filled-fb->frb.seen-fb->frb.unseen;
			if (just_seen!=0) {
				//means we can copy some in
				memmove(fb->base+unseen_end_mod,fb->frb.base+fb->frb.seen,just_seen);
				fb->frb.seen+=just_seen;
				assert(fb->frb.seen+fb->frb.unseen==fb->frb.filled);
				fb->unseen_end+=just_seen;
				unseen_end_mod+=just_seen;
				assert(unseen_end_mod<fb->size);
			}
			if (fb->frb.pad) { 
				size_t pad_size = 
				assert(unseen_start_mod>=unseen_end_mod);
				fprintf(stderr,"padding starting at %lu , size of %lu + %lu = %lu\n",unseen_end_mod,unseen_start_mod,unseen_end_mod, unseen_start_mod-unseen_end_mod);
				memset(fb->base+unseen_end_mod,'\n',unseen_start_mod-unseen_end_mod);
				fb->unseen_end+=unseen_start_mod-unseen_end_mod;
			}
		}
	}*/	

