#include <assert.h>
#include <ctype.h>
#include <errno.h>
#include <fcntl.h>
#include <inttypes.h>
#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <unistd.h>
#include <zlib.h>
#include <omp.h>

#include <sys/types.h>
#include <sys/stat.h>
#include <sys/time.h>

#include "../common/fasta.h"
#include "../common/util.h"

#define MAX_CONTIGS 1000000

char * contig_name[MAX_CONTIGS];
int contig_size[MAX_CONTIGS];
long long unsigned chunk_size[MAX_CONTIGS];

int n_seeds = 3;
int weight[20] = { 12, 12, 12 };

double overhead_mem = 1.5;

struct contig {
  char * name;
  unsigned int size;
  int chunk;
} contig[MAX_CONTIGS];

int n_contigs;


int cmp(const void *p1, const void *p2) {
  return ((struct contig *)p2)->size - ((struct contig *)p1)->size;
}


void read_contigs_from_fasta_file(char const * file_name) {
  fasta_t fasta_file;
  char * seq;

  fasta_file = fasta_open(file_name, MODE_LETTER_SPACE, false);
  if (fasta_file == NULL) {
    fprintf(stderr, "error: could not open genome file [%s]\n", file_name);
    exit(1);
  }

  fprintf(stderr, "scanning contigs...\n");
  while (n_contigs < MAX_CONTIGS && fasta_get_next_contig(fasta_file, &contig[n_contigs].name, &seq, NULL)) {
    contig[n_contigs].size = strlen(seq);
    fprintf(stderr, "%s %d\n", contig[n_contigs].name, contig[n_contigs].size);
    free(seq);
    n_contigs++;
  }
  if (n_contigs >= MAX_CONTIGS) {
    fprintf(stderr, "error: too many contigs\n");
    exit(1);
  }

  fasta_close(fasta_file);
}

void read_contig_list_from_stdin() {
  char buff[1000];
  while (n_contigs < MAX_CONTIGS && fscanf(stdin, "%s", buff) > 0) {
    contig[n_contigs].name = strdup(buff);
    int x = fscanf(stdin, "%d", &contig[n_contigs].size);
    if (x == 0){
    	//TODO check fscanf actually read a contig
    }
    n_contigs++;
  }
  if (n_contigs >= MAX_CONTIGS) {
    fprintf(stderr, "error: too many contigs\n");
    exit(1);
  }
}


int greedy_fit(long long unsigned target_len) {
  bool more = true;
  int n_chunks = 0;
  long long unsigned tmp;
  int i;

  for (i = 0; i < n_contigs; i++)
    contig[i].chunk = -1;

  while (more) {
    more = false;

    for (i = 0; i < n_contigs && contig[i].chunk >= 0; i++);
    assert(i < n_contigs);

    n_chunks++;
    contig[i].chunk = n_chunks-1;
    tmp = contig[i].size;

    for ( ; i < n_contigs; i++) {
      if (contig[i].chunk < 0) {
	if (tmp + contig[i].size < target_len) {
	  contig[i].chunk = n_chunks-1;
	  tmp += contig[i].size;
	} else {
	  more = true;
	}
      }
    }
  }

  return n_chunks;
}

int
main(int argc, char *argv[]) {
  double target_size;
  char * c;
  int i, j;
  double index_size = 0.0;

  if (argc < 3) {
    fprintf(stderr, "Usage: %s <genome_file> <target_RAM_size_in_GB> [<seed_weights>]\n", argv[0]);
    exit(1);
  }

  target_size = atof(argv[2]);
  if (target_size < 0.5 || target_size > 256.0) {
    fprintf(stderr, "error: the target memory size doesn't seem right in GB [%s]\n", argv[2]);
    exit(1);
  }

  if (argc >= 4) {
    n_seeds = 0;
    for (c = strtok(argv[3], ","); c != NULL && n_seeds < 16; c = strtok(NULL, ","), n_seeds++) {
      weight[n_seeds] = atoi(c);
      if (weight[n_seeds] < 5 || weight[n_seeds] > 16) {
	fprintf(stderr, "error: the seed weight [%s] is not in [5,16]\n", c);
	exit(1);
      }
    }
    if (n_seeds < 1 || c != NULL) {
      fprintf(stderr, "error: the number of seeds is not in [1,16]\n");
      exit(1);
    }
  }

  fprintf(stderr, "number of seeds: %d\n", n_seeds);
  fprintf(stderr, "seed weights: ");
  for (i = 0; i < n_seeds; i++) {
    fprintf(stderr, "%d%s", weight[i], i < n_seeds - 1? "," : "\n");
  }

  for (i = 0; i < n_seeds; i++) {
    double crt;
    long long int entries;
    entries = 1ll << (2 * weight[i]);
    crt = ((double)(entries * (sizeof(void *) + sizeof(uint32_t))))/(1024.0 * 1024.0 * 1024.0);
    fprintf(stderr, "seed %d: number of entries: %lld, index size: %.3f GB\n", i, entries, crt);
    index_size += crt;
  }

  // save some memory (1.5GB) for keeping reads and other data, including O/S
  if (target_size < overhead_mem + index_size) {
    fprintf(stderr, "error: not enough memory for current settings\n");
    exit(1);
  }

  long long unsigned target_len = (long long unsigned)(((target_size - overhead_mem - index_size)/(double)n_seeds) * 1024.0 * 1024.0 * 1024.0)/sizeof(uint32_t);

  fprintf(stderr, "target genome length per chunk: %llu\n", target_len);

  if (strcmp(argv[1], "-"))
    read_contigs_from_fasta_file(argv[1]);
  else
    read_contig_list_from_stdin();

  qsort(contig, n_contigs, sizeof(contig[0]), cmp);

  /*
  for (i = 0; i < n_contigs; i++) {
    fprintf(stdout, "%s\t%d\n", contig[i].name, contig[i].size);
  }
  */

  if (contig[0].size > target_len) {
    fprintf(stderr, "error: the largest contig [%s,%d] does not fit in target memory\n", contig[0].name, contig[0].size);
    exit(1);
  }

  int target_chunks;
  target_chunks = greedy_fit(target_len);

  fprintf(stderr, "target chunks: %d\n", target_chunks);

  if (target_chunks > 1) {
    fprintf(stderr, "rebalancing...\n");
    long long unsigned try_len;
    int try_chunks = 0;
    do {
      try_len = target_len;
      for (j = 0; j < 10; j++) {
	try_len = try_len - (target_len / 1000);
	fprintf(stderr, "trying chunk size: %llu.. ", try_len);
        if (try_len < contig[0].size) {
          fprintf(stderr, "no; too small\n");
          break;
        }
	try_chunks = greedy_fit(try_len);
	if (try_chunks <= target_chunks) {
	  fprintf(stderr, "yes\n");
	  break;
	}
	fprintf(stderr, "no\n");
      }
      if (j < 10 && !(try_len < contig[0].size)) {
	if (try_chunks < target_chunks) {
	  target_chunks = try_chunks;
	  fprintf(stderr, "target chunks: %d\n", target_chunks);
	}
	target_len = try_len;
      } else {
	break;
      }
    } while (1);
  }

  fprintf(stderr, "target genome length per chunk: %llu\n", target_len);
  fprintf(stderr, "estimated memory usage per chunk: %.3f GB\n",
	  (double)(target_len * sizeof(uint32_t) * n_seeds)/(1024.0 * 1024.0 * 1024.0) + index_size + overhead_mem);
	  

  target_chunks = greedy_fit(target_len);
  for (i = 0; i < target_chunks; i++) {
    chunk_size[i] = 0;
    fprintf(stdout, "chunk %d:\n", i+1);
    for (j = 0; j < n_contigs; j++) {
      if (contig[j].chunk == i) {
	chunk_size[i] += contig[j].size;
	fprintf(stdout, "%s\t%d\n", contig[j].name, contig[j].size);
      }
    }
  }

  for (i = 0; i < target_chunks; i++) {
    fprintf(stderr, "chunk %d: %llu\n", i+1, chunk_size[i]);
  }

  return 0;
}
