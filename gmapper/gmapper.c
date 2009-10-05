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

#include <xmmintrin.h>	// for _mm_prefetch

#include <sys/types.h>
#include <sys/stat.h>
#include <sys/time.h>

#include "../common/hash.h"
#include "../common/fasta.h"
#include "../common/util.h"
#include "../gmapper/gmapper.h"
#include "../mapper/mapper.h"

/* Kmer to genome index */
static uint32_t ***genomemap;
static uint32_t **genomemap_len;

/* offset info for genome contigs */
static uint32_t *genome_offsets;
static uint32_t num_contigs;

static count_t mem_genomemap;

extern size_t
power(size_t base, size_t exp);

extern uint32_t
kmer_to_mapidx_hash(uint32_t *kmerWindow, u_int sn);

/*
 * index the kmers in the genome contained in.
 * This can then be used to align reads against.
 */
static bool
load_genome_lscs(const char *file)
{
	fasta_t fasta;
	int space; //colour space or letter space
	size_t seqlen,  bytes;
	uint32_t *read;
	char *seq, *name;
	int16_t initbp;
	uint32_t *kmerWindow;
	uint sn;

	//allocate memory for the genome map
	genomemap = (uint32_t ***) xmalloc_c(n_seeds * sizeof(genomemap[0]),
			&mem_genomemap);
	genomemap_len = (uint32_t **)xmalloc_c(n_seeds * sizeof(genomemap_len[0]),
			&mem_genomemap);


	for (sn = 0; sn < n_seeds; sn++) {

		bytes = sizeof(uint32_t) * power(4, HASH_TABLE_POWER);

		genomemap[sn] = (uint32_t **)xcalloc_c(bytes, &mem_genomemap);
		genomemap_len[sn] = (uint32_t *)xcalloc_c(bytes,&mem_genomemap);
	}

	if (shrimp_mode == MODE_LETTER_SPACE)
		space = LETTER_SPACE;
	else
		space = COLOUR_SPACE;

	//open the fasta file and check for errors
	fasta = fasta_open(file,space);
	if (fasta == NULL) {
		fprintf(stderr,"error: failded to open genome file [%s]\n",file);
		return (false);
	} else {
		fprintf(stderr,"- Processing genome file [%s]\n",file);
	}

	//Read the contigs and record their sizes
	num_contigs = 0;
	u_int i = 0;
	while(fasta_get_next(fasta, &name, &seq, NULL)){

		num_contigs++;
		genome_offsets = (uint32_t *)xrealloc(genome_offsets,sizeof(uint32_t)*num_contigs);
		genome_offsets[num_contigs] = i;

		if (strchr(name, '\t') != NULL || strchr(seq, '\t') != NULL) {
			fprintf(stderr, "error: tabs are not permitted in fasta names "
					"or sequences. Tag: [%s].\n", name);
			exit(1);
		}

		seqlen = strlen(seq);
		if (seqlen == 0) {
			fprintf(stderr, "error: genome [%s] had no sequence!\n",
					name);
			exit(1);
		}
//		if (seqlen > MAX_READ_LENGTH) {
//			fprintf(stderr, "error: genome [%s] had unreasonable "
//					"length!\n", name);
//		}
		if (shrimp_mode == MODE_COLOUR_SPACE) {
			/* the sequence begins with the initial letter base */
			if (seqlen < 1) {
				fprintf(stderr, "error: genome [%s] had sequence "
						"with no colours!\n", name);
				exit(1);
			}
			seqlen--;
		}

		read = fasta_sequence_to_bitfield(fasta,seq);
		if (read == NULL) {
			fprintf(stderr, "error: invalid sequence; tag: [%s]\n", name);
			exit(1);
		}

		if (shrimp_mode == MODE_COLOUR_SPACE){
			initbp = fasta_get_initial_base(fasta,seq);
		}

		kmerWindow = (uint32_t *)xcalloc(sizeof(kmerWindow[0])*BPTO32BW(max_seed_span));

		u_int load;

		for (load = 0; i < seqlen + genome_offsets[num_contigs]; i++) {
			uint base, sn;

			base = EXTRACT(read, i);
			bitfield_prepend(kmerWindow, max_seed_span, base);

			//skip past any Ns or Xs
			if (base == BASE_N || base == BASE_X)
				load = 0;
			else if (load < max_seed_span)
				load++;

			for (sn = 0; sn < n_seeds; sn++) {
				if (load < seed[sn].span)
					continue;

				/*
				 * For simplicity we throw out the first kmer when in colour space. If
				 * we did not do so, we'd run into a ton of colour-letter space
				 * headaches. For instance, how should the first read kmer match
				 * against a kmer from the genome? The first colour of the genome kmer
				 * depends on the previous letter in the genome, so we may have a
				 * matching read, but the colour representation doesn't agree due to
				 * different initialising bases.
				 *
				 * If we wanted to be complete, we could compute the four permutations
				 * and add them, but I'm not so sure that'd be a good idea. Perhaps
				 * this should be investigated in the future.
				 */
				if (shrimp_mode == MODE_COLOUR_SPACE && i == seed[sn].span - 1)
					continue;

				uint32_t mapidx = kmer_to_mapidx_hash(kmerWindow, sn);

				genomemap_len[sn][mapidx]++;
				genomemap[sn][mapidx] = (uint32_t *)xrealloc_c(genomemap[sn][mapidx],
						sizeof(uint32_t) * (genomemap_len[sn][mapidx] + 1),
						sizeof(uint32_t) * genomemap_len[sn][mapidx],
						&mem_genomemap);
				genomemap[sn][mapidx][genomemap_len[sn][mapidx]] = i;
				nkmers++;

			}
		}

		free(seq);
		free(kmerWindow);
	}
	fasta_close(fasta);

	fprintf(stderr,"Loaded Genome\n");

	return (true);

}

//testing main
//TODO ensure variables are initialized properly.
int main(int argc, char **argv){
	char *genome_file;
	if (argc > 1){
		genome_file = argv[1];
		set_mode_from_argv(argv);
		load_genome_lscs(genome_file);
	}
}
