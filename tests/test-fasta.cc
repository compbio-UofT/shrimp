#include <iostream>
#include <zlib.h>

#include "../common/fasta.h"
#include "../common/util.h"

using namespace std;


int file_format, out_format, len=10;
fasta_t fasta_file;
char * name, * seq;
bool is_rna;
uint32_t * genome_ls, * genome_cs;
int genome_len;
uint8_t base;


int main(int argc, char * argv[]) {
  if (argc < 4) {
    cerr << "Use: " << argv[0] << " <file> <file format: 0=LS/1=CS> <out format: 0=LS/1=CS> [<length>=10]" << endl;
    exit(1);
  }

  file_format = (argv[2][0] == '1' ? COLOUR_SPACE : LETTER_SPACE);
  out_format = (argv[3][0] == '1' ? COLOUR_SPACE : LETTER_SPACE);
  if (argc >= 5) {
    len = atoi(argv[4]);
    if (len <= 0) {
      cerr << "error: length too small: \"" << argv[4] << "\"" << endl;
      exit(1);
    }
  }

  fasta_file = fasta_open(argv[1], file_format);
  if (fasta_file == NULL) {
    cerr << "error: could not open file \"" << argv[1] << "\"" << endl;
    exit(1);
  }

  if (!fasta_get_next(fasta_file, &name, &seq, &is_rna)) {
    cerr << "error: could not read fasta sequence" << endl;
    exit(1);
  }

  genome_len = strlen(seq);

  if (file_format == LETTER_SPACE) {
    genome_ls = fasta_sequence_to_bitfield(fasta_file, seq);
    genome_cs = fasta_bitfield_to_colourspace(fasta_file, genome_ls, genome_len, is_rna);
  } else {
    assert(0);
  }

  cout << "LS: ";
  for (int i = 0; i < len; i++) {
    base = EXTRACT(genome_ls, i);
    cout << base_translate(base, false);
  }
  cout << endl;

  cout << "CS: ";
  for (int i = 0; i < len; i++) {
    base = EXTRACT(genome_cs, i);
    cout << base_translate(base, true);
  }
  cout << endl;

  free(genome_cs);
  free(genome_ls);
  free(name);
  free(seq);
  fasta_close(fasta_file);

  return 0;
}
