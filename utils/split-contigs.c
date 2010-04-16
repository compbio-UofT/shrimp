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

struct contig {
  char * name;
  unsigned int size;
  int used;
} contig[MAX_CONTIGS];


int cmp(const void *p1, const void *p2) {
  return ((struct contig *)p2)->size - ((struct contig *)p1)->size;
}


int
main(int argc, char *argv[]) {
  fasta_t fasta_file;
  double target_size;
  int n_contigs = 0;
  int n_seeds = 3;
  int seed_weight = 12;
  char * seq;
  int i;

  if (argc < 3) {
    fprintf(stderr, "Usage: %s <genome_file> <target_RAM_size_in_GB> [<n_seeds> = 3] [<seed_weight> = 12]\n", argv[0]);
    exit(1);
  }

  fasta_file = fasta_open(argv[1], LETTER_SPACE);
  if (fasta_file == NULL) {
    fprintf(stderr, "error: could not open genome file [%s]\n", argv[1]);
    exit(1);
  }

  target_size = atof(argv[2]);
  if (target_size < 0.5 || target_size > 256.0) {
    fprintf(stderr, "error: the target memory size doesn't seem right in GB [%s]\n", argv[2]);
    exit(1);
  }

  if (argc >= 4) {
    n_seeds = atoi(argv[3]);
    if (n_seeds < 1 || n_seeds > 16) {
      fprintf(stderr, "error: the number of seeds [%s] is not in [1,16]\n", argv[3]);
      exit(1);
    }
  }

  if (argc >= 5) {
    seed_weight = atoi(argv[4]);
    if (seed_weight < 7 || seed_weight > 16) {
      fprintf(stderr, "error: the seed weight [%s] is not in [7,16]\n", argv[4]);
      exit(1);
    }
  }

  fprintf(stderr, "number of seeds: %d\n", n_seeds);
  fprintf(stderr, "seed weight: %d\n", seed_weight);

  long long unsigned entries = 1llu << (2 * seed_weight);
  fprintf(stderr, "number of entries: %llu\n", entries);

  double index_size = ((double)(entries * (sizeof(void *) + sizeof(uint32_t))))/(1024.0 * 1024.0 * 1024.0);

  fprintf(stderr, "per-seed index size: %.3f GB\n", index_size);

  // save 0.5GB for keeping reads and other data
  long long unsigned target_len = (long long unsigned)((((target_size - 0.5)/(double)n_seeds) - index_size) * 1024.0 * 1024.0 * 1024.0)/sizeof(uint32_t);

  fprintf(stderr, "target genome length per chunk: %llu\n", target_len);

  fprintf(stderr, "scanning contigs...\n");
  while (n_contigs < MAX_CONTIGS && fasta_get_next(fasta_file, &contig[n_contigs].name, &seq, NULL)) {
    contig[n_contigs].size = strlen(seq);
    contig[n_contigs].used = 0;
    free(seq);
    n_contigs++;
  }
  if (n_contigs >= MAX_CONTIGS) {
    fprintf(stderr, "error: too many contigs\n");
    exit(1);
  }

  fasta_close(fasta_file);

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

  bool more = true;
  int chunk = 0;
  long long unsigned tmp;

  while (more) {
    more = false;

    for (i = 0; i < n_contigs && contig[i].used; i++);
    if (i == n_contigs)
      break;

    chunk++;
    fprintf(stdout, "chunk %d:\n", chunk);
    fprintf(stdout, "%s\t%d\n", contig[i].name, contig[i].size);
    contig[i].used = 1;
    tmp = contig[i].size;

    for ( ; i < n_contigs; i++) {
      if (!contig[i].used) {
	if (tmp + contig[i].size < target_len) {
	  fprintf(stdout, "%s\t%d\n", contig[i].name, contig[i].size);
	  contig[i].used = 1;
	  tmp += contig[i].size;
	} else {
	  more = true;
	}
      }
    }
  }


  return 0;
}
