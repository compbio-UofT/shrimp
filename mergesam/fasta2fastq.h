#ifndef __MERGESAM_H__
#define __MERGESAM_H__
#define GROWTH_FACTOR 1.3
#define SIZE_READ_NAME 255
#define MIN(a,b) (((a) < (b)) ? (a) : (b))
#define MAX_INT32	2147483647

#define DEF_QV_OFFSET	-1


//IO SETTINGS
#define DEF_READ_SIZE	1024*10
#define DEF_BUFFER_SIZE	1024*1024*50
#define DEF_READ_RATE 6000
typedef struct runtime_options {
	//input options
	size_t buffer_size;
	size_t read_size;
	int threads;
	int qv_offset;
};
extern runtime_options options;
#endif
