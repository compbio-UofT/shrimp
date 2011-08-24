
#include <zlib.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include "file_buffer.h"

//#define MIN(a,b) (((a) < (b)) ? (a) : (b))
#ifdef __APPLE__
void * memrchr(const void *s, int c, size_t n) {
    const unsigned char *cp;
    if (n != 0) {
	cp = (unsigned char *)s + n;
	do {
	    if (*(--cp) == (unsigned char)c)
		return((void *)cp);
 	} while (--n != 0);
    }
    return((void *)0);
}
#endif


bool auto_detect_fastq(gzFile fp, char * reads_filename) {
	// autodetect fasta/fastq format
	bool fastq=false;
	char c;
	c = gzgetc(fp);
	while (c == '#' || c == ';') {
		// discard this line
		while (c != -1 && c != '\n')
			c = gzgetc(fp);

		if (gzeof(fp))
			break;
		if (c == -1) {
			fprintf(stderr, "did not find the end of a comment line in the input file [%s]. try disabling input autodetection\n", reads_filename);
			exit(1);
		}

		c = gzgetc(fp);
	}
	if (!gzeof(fp)) {
		if (c == -1) {
			fprintf(stderr, "did not find a non-comment line in the input file [%s]. try disabling input autodetection\n", reads_filename);
			exit(1);
		}
		if (c == '@') {
			fprintf(stderr,"detected fastq format in input file [%s]\n", reads_filename);
			fastq = true;
		} else if (c == '>') {
			fprintf(stderr, "detected fasta format in input file [%s]\n", reads_filename);
			fastq = false;
		} else {
			fprintf(stderr, "unrecognized character [%c] in input file [%s]. try disabling input autodetection\n", (char)c, reads_filename);
			exit(1);
		}
		gzungetc(c, fp);
	}
	return fastq;
}

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
	fb->unseen_inter=0;
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

//find the last new line in the area and 
void mark_partial_line(file_read_buffer * frb,size_t offset) {
	assert(frb->seen+frb->unseen==frb->filled);
	//find out if the offset is bigger then what has not been seen yet
	// if it is use this instead of the offset
	//for example if we have 20 unseen characters and user is requesting offset of 40 characters, only look at the last 20
	frb->exhausted=offset>=frb->unseen;
	size_t unseen=frb->exhausted ? frb->unseen : offset;
	//find the last newline, set the unseen count accordingly
	char * ptr = (char*)memrchr(frb->base+frb->seen,'\n',unseen);
	if (ptr!=NULL) {
		size_t offset=(ptr-frb->base)+1;
		frb->unseen=frb->filled-offset;
	} else if (frb->eof==1) {	
		fprintf(stderr,"EOF\n");
		frb->unseen=0;
	}
}


//unseen_start is always 0 based actual location of start
//unseen_end is always 0 based location of next to be written data
void add_read_buffer_to_main(file_buffer * fb) {
	size_t unseen_end_mod=fb->unseen_end%fb->size;
	size_t unseen_start_mod=fb->unseen_start%fb->size;
	fb->changed=false;
	assert(fb->frb.seen+fb->frb.unseen==fb->frb.filled);
	assert(fb->unseen_start<=fb->unseen_end);
	if (fb->unseen_start+fb->size<fb->unseen_end) {
		fprintf(stderr,"%lu+%lu<%lu\n",fb->unseen_start,fb->size,fb->unseen_end);
	}
	assert(fb->unseen_start+fb->size>=fb->unseen_end);
	//while the buffer is not filled to the top and we have stuff left to read in the frb
	while (fb->unseen_end!=fb->unseen_start+fb->size && !fb->frb.exhausted) {
		assert(fb->unseen_start+fb->size>=fb->unseen_end);
		assert(fb->frb.seen+fb->frb.unseen==fb->frb.filled);
		assert(fb->unseen_start<=fb->unseen_end);
		//find out how much continous space we have to the end of the buffer or start
		size_t space=(unseen_start_mod>unseen_end_mod ? unseen_start_mod : fb->size ) - unseen_end_mod;
		assert(space>0);
		//find out how far to look in buffer for last newline
		mark_partial_line(&fb->frb,space);
		size_t just_seen = fb->frb.filled-fb->frb.seen-fb->frb.unseen;
		//unseen was just updated in mark_partial_line if some new lines were found
		//just_seen being non-zero means that we found some new lines to move over, so do it
		if (just_seen!=0) {
			memmove(fb->base+unseen_end_mod,fb->frb.base+fb->frb.seen,just_seen);
			fb->changed=true;
			fb->frb.seen+=just_seen;
			assert(fb->frb.seen+fb->frb.unseen==fb->frb.filled);
			fb->unseen_end+=just_seen;
			unseen_end_mod=fb->unseen_end%fb->size;
		//in this case we did not find any new lines to move over and there is more
		//data waiting on the frb, so fb needs to go into exhausted state
		} else if (!fb->frb.exhausted) {
			fb->exhausted=true;
		}
		if (!fb->frb.exhausted) {
			//if the frb is not exhausted that means fb space was limiting, so fill in the rest with 0's
			//check if start_mod < end_mod or of the case if buffer has 0 bytes unseen
			if (unseen_start_mod<unseen_end_mod || fb->unseen_start==fb->unseen_end) { 
				assert(fb->size>=unseen_end_mod);
				//fill to end with zeros, since can't fit any more in, because !fb->frb.exhausted
				memset(fb->base+unseen_end_mod,'\0',fb->size-unseen_end_mod);
				fb->changed=true;
				fb->unseen_end+=fb->size-unseen_end_mod;
				unseen_end_mod=fb->unseen_end%fb->size;
			} else {
				fb->exhausted=true;	
				assert(fb->unseen_start+fb->size>=fb->unseen_end);
				return;
			}
		}
		unseen_end_mod=fb->unseen_end%fb->size;
	}
	if (fb->unseen_end==fb->unseen_start+fb->size || (fb->frb.exhausted && fb->frb.eof==1)) {
		fb->exhausted=true;
	}
	assert(fb->frb.seen+fb->frb.unseen==fb->frb.filled);
	assert(fb->unseen_start<=fb->unseen_end);
	assert(fb->unseen_start+fb->size>=fb->unseen_end);
}


