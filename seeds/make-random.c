#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>

#include "../common/util.h"


int a[64];

int main(int argc, char * argv[]) {
  int n, m, i;

  if (argc < 3) {
    fprintf(stderr, "use: %s <length> <weight>\n", argv[0]);
    exit(1);
  }

  n = atoi(argv[1]);
  if (n <= 0 || n > 63) {
    fprintf(stderr, "error: invalid length [%s]\n", argv[1]);
    exit(1);
  }

  m = atoi(argv[2]);
  if (m <= 0 || m > n) {
    fprintf(stderr, "error: invalid weight [%s]\n", argv[2]);
    exit(1);
  }

  srand48((int)gettimeinusecs());
  for (i = 0; i < m; i++) {
    int j;
    do {
      j = (int)(drand48()*n);
    } while (a[j] == 1);
    a[j] = 1;
  }

  for (i = 0; i < n; i++) {
    fprintf(stdout, "%s", a[i] == 1? "1" : "0");
  }
  fprintf(stdout, "\n");

  return 0;
}


      
