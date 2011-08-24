#ifndef __FILE_BUFFER__
#define __FILE_BUFFER__
#include <stdbool.h>
#include <zlib.h>

typedef struct file_read_buffer {
        gzFile file;
        char * base;
        size_t size;
        size_t filled;
        size_t unseen;
	size_t seen;
        int eof;
	bool exhausted;
} file_read_buffer;

typedef struct file_buffer {
        char * base;
        size_t size;
        size_t unseen_start;
        size_t unseen_end;
	size_t unseen_inter;
	bool changed;
	bool exhausted;
        file_read_buffer frb;
} file_buffer;

bool auto_detect_fastq(gzFile fp, char * reads_filename);
void fb_close(file_buffer * fb);
file_buffer * fb_open(char * path,size_t buffer_size,size_t read_size);
void fill_read_buffer(file_read_buffer * frb);
void mark_partial_line(file_read_buffer * frb, size_t);
void add_read_buffer_to_main(file_buffer * fb);


#ifdef __APPLE__
void * memrchr(const void *s, int c, size_t n);
#endif
#endif
