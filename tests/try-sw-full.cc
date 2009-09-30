#include <iostream>
#include <cstring>

using namespace std;

extern "C" {

#include "../common/util.h"
#include "../common/sw-full-common.h"
#include "../common/sw-full-cs.h"

}


int costMatch = 10;
int costMismatch = -15;
int costGapOpen = -40;
int costGapExtend = -7;
int costCrossover = -14;

uint32_t *db, *qr;
int dbLen, qrLen;
int initbp;

fasta_t fasta;
char *name, *sequence;

struct sw_full_results sfr;

int main(int argc, char *argv[]) {
  if (argc < 2) {
    cerr << "Use: " << argv[0] << " <db-seq> <qr-seq>" << endl;
    exit(1);
  }

  fasta = fasta_open(argv[1], LETTER_SPACE);
  if (fasta == NULL) {
    cerr << "Could not open file " << argv[1] << endl;
    exit(1);
  }
  if (!fasta_get_next(fasta, &name, &sequence, NULL)) {
    cerr << "Could not get sequence from file " << argv[1] << endl;
  }
  dbLen = strlen(sequence);
  db = fasta_sequence_to_bitfield(fasta, sequence);
  fasta_close(fasta);
  free(name);
  free(sequence);

  fasta = fasta_open(argv[2], COLOUR_SPACE);
  if (fasta == NULL) {
    cerr << "Could not open file " << argv[2] << endl;
    exit(1);
  }
  if (!fasta_get_next(fasta, &name, &sequence, NULL)) {
    cerr << "Could not get sequence from file " << argv[2] << endl;
  }
  initbp = fasta_get_initial_base(fasta, sequence);
  qrLen = strlen(sequence) - 1;
  qr = fasta_sequence_to_bitfield(fasta, sequence);
  fasta_close(fasta);
  free(name);
  free(sequence);

  sw_full_cs_setup(dbLen, qrLen, costGapOpen, costGapExtend, costMatch, costMismatch, costCrossover, false);

  sw_full_cs(db, 0, dbLen, qr, qrLen, initbp, 0, &sfr, false, false);

  cout << sfr.dbalign << endl;
  cout << sfr.qralign << endl;
  cout << "score: " << sfr.score << endl;

  return 0;
}
