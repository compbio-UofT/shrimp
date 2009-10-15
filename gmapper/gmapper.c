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

#include <xmmintrin.h>	// for _mm_prefetch

#include <sys/types.h>
#include <sys/stat.h>
#include <sys/time.h>

#include "../common/hash.h"
#include "../common/fasta.h"
#include "../common/util.h"
#include "../gmapper/gmapper.h"
#include "../mapper/mapper.h"
#include "../common/version.h"

/* Kmer to genome index */
static uint32_t ***genomemap;
static uint32_t **genomemap_len;

/* offset info for genome contigs */
static uint32_t *contig_offsets;
static char **contig_names = NULL;
static uint32_t num_contigs;

static count_t mem_genomemap;

static uint numThreads = omp_get_num_procs();
static uint chunkSize = 1000;

extern size_t
power(size_t base, size_t exp);

extern uint32_t
kmer_to_mapidx_hash(uint32_t *kmerWindow, u_int sn);

void print_info(){
	fprintf(stderr,"num_contigs = %u\n",num_contigs);
	uint i;
	for (i = 0; i < num_contigs; i++){
		fprintf(stderr,"contig %u: name=%s offset = %u\n",i,contig_names[i],contig_offsets[i]);
	}
	fprintf(stderr,"n_seeds = %u\n",n_seeds);
	for (i=0; i < n_seeds;i++){
		fprintf(stderr,"seed %u: span=%u\n",i,seed[i].span);
	}

}

static bool save_genome_map(const char *file) {
	gzFile fp = gzopen(file, "wb");
	if (fp == NULL){
		return false;
	}

	//write the header
	//fprintf(stderr,"saving num_contgs: %u\n",num_contigs);
	uint32_t i;
	gzwrite(fp,&num_contigs,sizeof(uint32_t));
	//fprintf(stderr,"saved num_contigs\n");
	for (i = 0; i < num_contigs; i++) {
		gzwrite(fp,&contig_offsets[i],sizeof(uint32_t));
		uint32_t len = strlen(contig_names[i]);
		gzwrite(fp,&len,sizeof(uint32_t));
		gzwrite(fp,contig_names[i],len +1);
	}
	//fprintf(stderr,"saving seeds\n");
	//write the seeds and genome_maps
	gzwrite(fp,&n_seeds,sizeof(uint32_t));
	for (i = 0; i < n_seeds; i++) {
		gzwrite(fp,&seed[i], sizeof(seed_type));

		//write the genome_map for this seed
		uint32_t j;
		uint32_t p = power(4, HASH_TABLE_POWER);
		//fprintf(stderr,"saving index\n");
		gzwrite(fp,&p, sizeof(uint32_t));
		for (j = 0; j < p; j++) {
			uint32_t len = genomemap_len[i][j];
			gzwrite(fp, &len, sizeof(uint32_t));
			gzwrite(fp, genomemap[i][j], sizeof(uint32_t) * len);
		}
	}
	gzclose(fp);
	return true;
}

static bool load_genome_map(const char *file){
	gzFile fp = gzopen(file,"rb");
	if (fp == NULL){
		return false;
	}

	uint32_t i;
	gzread(fp,&num_contigs,sizeof(uint32_t));
	//fprintf(stderr,"num_contigs = %u\n",num_contigs);

	contig_names = (char **)xrealloc(contig_names,sizeof(char *)*num_contigs);
	contig_offsets = (uint32_t *)xrealloc(contig_offsets,sizeof(uint32_t)*num_contigs);

	for (i = 0; i < num_contigs; i++){
		gzread(fp,&contig_offsets[i],sizeof(uint32_t));

		uint32_t len;
		gzread(fp,&len,sizeof(uint32_t));
		contig_names[i] = (char *)xrealloc(contig_names[i],sizeof(char)*len);
		gzread(fp,contig_names[i],len+1);
		//fprintf(stderr,"contig %u: name=%s offset = %u\n",i,contig_names[i],contig_offsets[i]);
	}
	gzread(fp,&n_seeds,sizeof(uint32_t));
	//fprintf(stderr,"n_seeds = %u\n",n_seeds);
	seed =(seed_type *)xrealloc(seed,sizeof(seed_type)*n_seeds);
	genomemap_len = (uint32_t **)xrealloc(genomemap_len,sizeof(uint32_t *)*n_seeds);
	genomemap = (uint32_t ***)xrealloc(genomemap,sizeof(uint32_t **) * n_seeds);
	for (i = 0; i < n_seeds; i++){
		gzread(fp,&seed[i], sizeof(seed_type));
		//fprintf(stderr,"seed %u: span=%u\n",i,seed[i].span);
		uint32_t j;
		uint32_t p;
		gzread(fp,&p,sizeof(uint32_t));
		genomemap_len[i] = (uint32_t *)xrealloc(genomemap_len[i],sizeof(uint32_t)*p);
		genomemap[i] = (uint32_t **)xrealloc(genomemap[i],sizeof(uint32_t *)*p);
		for (j = 0; j < p; j++){
			gzread(fp,&genomemap_len[i][j],sizeof(uint32_t));
			genomemap[i][j] = (uint32_t *)xrealloc(genomemap[i][j],
					sizeof(uint32_t)*genomemap_len[i][j]);
			gzread(fp,genomemap[i][j],sizeof(uint32_t) * genomemap_len[i][j]);
		}
	}
	gzclose(fp);
	return true;
}

/*
 * Launch the threads that will scan the reads
 */

static bool
launch_scan_threads(const char *file){
	fasta_t fasta;
	char **seq, **name;
	int space;
	int s = chunkSize*numThreads;

	seq = (char **)xmalloc(sizeof(char *)*s);
	name = (char **)xmalloc(sizeof(char *)*s);

	//open the fasta file and check for errors
	if (shrimp_mode == MODE_LETTER_SPACE)
			space = LETTER_SPACE;
		else
			space = COLOUR_SPACE;
	fasta = fasta_open(file,space);
	if (fasta == NULL) {
		fprintf(stderr,"error: failded to open read file [%s]\n",file);
		return (false);
	} else {
		fprintf(stderr,"- Processing read file [%s]\n",file);
	}

	//read the fasta file, s sequences at a time, and process in threads.
	bool more = true;
	while (more){
		int i;
		for(i = 0; i < s; i++){
			if(!fasta_get_next(fasta, name + i, seq + i, NULL)){
				more = false;
				break;
			}
		}
#pragma omp parallel shared(seq,name,i) num_threads(numThreads)
		{
			int j;
#pragma omp for
			for (j = 0; j < i; j++){
				//TODO work here
			}
		}

	}
	return true;
}

/*
 * index the kmers in the genome contained in the file.
 * This can then be used to align reads against.
 */
static bool
load_genome_lscs(char **files, int nfiles)
{
	fasta_t fasta;
	size_t seqlen,  bytes;
	uint32_t *read;
	char *seq, *name;
	int16_t initbp;
	uint32_t *kmerWindow;
	uint sn;
	char *file;

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

	int cfile;
	for(cfile = 0; cfile < nfiles; cfile++){
		file = files[cfile];
		//open the fasta file and check for errors
		fasta = fasta_open(file,LETTER_SPACE);
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
			contig_offsets = (uint32_t *)xrealloc(contig_offsets,sizeof(uint32_t)*num_contigs);
			contig_offsets[num_contigs - 1] = i;
			contig_names = (char **)xrealloc(contig_names,sizeof(char *)*num_contigs);
			contig_names[num_contigs - 1] = name;

			fprintf(stderr,"- Processing contig %s\n",name);

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
			DEBUG("indexing sequnce");
			u_int load;
			DEBUG("looping seq");
			for (load = 0; i < seqlen + contig_offsets[num_contigs]; i++) {
				uint base, sn;

				base = EXTRACT(read, i);
				bitfield_prepend(kmerWindow, max_seed_span, base);

				//skip past any Ns or Xs
				if (base == BASE_N || base == BASE_X)
					load = 0;
				else if (load < max_seed_span)
					load++;
				DEBUG("looping seeds");
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
					DEBUG("hashing");
					uint32_t mapidx = kmer_to_mapidx_hash(kmerWindow, sn);
					//increase the match count and store the location of the match
					DEBUG("updating len");
					genomemap_len[sn][mapidx]++;
					DEBUG("reallocing map");
					genomemap[sn][mapidx] = (uint32_t *)xrealloc_c(genomemap[sn][mapidx],
							sizeof(uint32_t) * (genomemap_len[sn][mapidx]),
							sizeof(uint32_t) * (genomemap_len[sn][mapidx] - 1),
							&mem_genomemap);
					if (genomemap[sn][mapidx] == NULL){
						DEBUG("realloc returned null");
					}
					DEBUG("updateing map");
					genomemap[sn][mapidx][genomemap_len[sn][mapidx] - 1] = i;
					nkmers++;

				}
			}
			DEBUG("done indexing sequence");
			free(seq);
			free(name);
			seq = name = NULL;

			free(kmerWindow);
		}
		fasta_close(fasta);
	}
	fprintf(stderr,"Loaded Genome\n");

	return (true);

}

void usage(char *progname,bool full_usage){
	char *slash;
	uint sn;

	load_default_seeds();
	slash = strrchr(progname, '/');
	if (slash != NULL)
		progname = slash + 1;
	fprintf(stderr, "usage: %s [parameters] [options] "
			"genome_file1\n", progname);

	fprintf(stderr, "Parameters:\n");

	fprintf(stderr,
			"    -s    Spaced Seed(s)                          (default: ");
	for (sn = 0; sn < n_seeds; sn++) {
		if (sn > 0)
			fprintf(stderr, "                                                            ");
		fprintf(stderr, "%s%s\n", seed_to_string(sn), (sn == n_seeds - 1? ")" : ","));
	}
	fprintf(stderr,
			"    -?    Full list of parameters and options\n");

}

//testing main
//int main(int argc, char **argv){
//	char *genome_file;
//	if (argc > 1){
//		genome_file = argv[1];
//		set_mode_from_argv(argv);
//		if (n_seeds == 0)
//				load_default_seeds();
//		init_seed_hash_mask();
//		if (0){
//			fprintf(stderr,"loading gneomfile\n");
//			load_genome_lscs(genome_file);
//			fprintf(stderr,"saving compressed index\n");
//			save_genome_map("testfile.gz");
//			print_info();
//		}
//		if (1){
//			fprintf(stderr,"loading compressed index\n");
//			load_genome_map("testfile.gz");
//		}
//		fprintf(stderr,"done\n");
//	}
//}

int main(int argc, char **argv){
	char *genome_file;
	char *progname = argv[0];
	char const * optstr;
	char *c;
	int ch;

	set_mode_from_argv(argv);

	switch(shrimp_mode){
	case MODE_COLOUR_SPACE:
		optstr = "?s:";
		break;
	case MODE_LETTER_SPACE:
		optstr = "?s:";
		break;
	case MODE_HELICOS_SPACE:
		fprintf(stderr,"Helicose currently unsuported\n");
		exit(1);
		break;
	default:
		assert(0);
	}

	fprintf(stderr, "--------------------------------------------------"
			"------------------------------\n");
	fprintf(stderr, "gmapper: %s.\nSHRiMP %s\n[%s]\n", get_mode_string(),
			SHRIMP_VERSION_STRING, get_compiler());
	fprintf(stderr, "--------------------------------------------------"
			"------------------------------\n");

	while ((ch = getopt(argc,argv,optstr)) != -1){
		switch (ch) {
		case 's':
			if (strchr(optarg, ',') == NULL) { // allow comma-separated seeds
				if (!add_spaced_seed(optarg)) {
					fprintf(stderr, "error: invalid spaced seed \"%s\"\n", optarg);
					exit (1);
				}
			} else {
				c = strtok(optarg, ",");
				do {
					if (!add_spaced_seed(c)) {
						fprintf(stderr, "error: invalid spaced seed \"%s\"\n", c);
						exit (1);
					}
					c = strtok(NULL, ",");
				} while (c != NULL);
			}
			break;
		case '?':
			usage(progname, true);
			break;
		default:
			usage(progname, false);
		}
	}

	argc -= optind;
	argv += optind;

	if (n_seeds == 0)
			load_default_seeds();

	init_seed_hash_mask();

	if (argc < 1){
		fprintf(stderr, "Genome File not specified\n");
		exit(1);
	}

	genome_file = argv[0];

	fprintf(stderr,"loading gneomfile\n");
	load_genome_lscs(&genome_file,1);
	fprintf(stderr, "        Genomemap:                %s\n",
					comma_integer(count_get_count(&mem_genomemap)));
	fprintf(stderr,"saving compressed index\n");
	save_genome_map("testfile.gz");
	print_info();
	fprintf(stderr,"loading compressed index\n");
	load_genome_map("testfile.gz");

	launch_scan_threads(genome_file);

}
