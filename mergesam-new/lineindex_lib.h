#ifndef __LINEINDEX_H__
#define __LINEINDEX_H__

typedef struct lineindex_table lineindex_table;
struct lineindex_table {
	size_t size;
	size_t start;
	size_t end;
	char ** table;
};
lineindex_table * lineindex_init(size_t size);
size_t add_lineindex_from_memory_threaded(lineindex_table * lineindex, lineindex_table ** thread_lineindexes,char *data, size_t data_size, int threads,char comment_start);
void lineindex_destroy(lineindex_table * lt);
#endif
