#define _MODULE_GENOME

#include "genome.h"


/*
 * Loading and saving the genome projection.
 */
bool save_genome_map_seed(const char *file, int sn)
{
  /*
   * This function writes the .seed.X file appropriate for loading
   * file is the file name to write to and sn is the seed number to write
   *
   * The file format is a gziped binary format as follows
   *
   * uint32_t				: shrimp_mode
   * uint32_t				: Hflag
   * seed_type			: Seed
   * uint32_t				: capacity
   * uint32_t * capacity	: genomemap_len
   * uint32_t				: total (= sum from 0 to capacity - 1 of genomemap_len)
   * uint32_t * total		: genomemap (each entry of length genomemap_len)
   *
   */
  gzFile fp = gzopen(file, "wb");
  if (fp == NULL){
    return false;
  }

  // shrimp_mode
  uint32_t m;
  m = (uint32_t)shrimp_mode;
  xgzwrite(fp, &m, sizeof(uint32_t));

  // Hflag
  uint32_t h = (uint32_t)Hflag;
  xgzwrite(fp, &h, sizeof(uint32_t));

  // Seed
  xgzwrite(fp, &seed[sn], sizeof(seed_type));

  // genomemap_len
  uint32_t capacity;
  if (Hflag) {
    capacity = (uint32_t)power4(HASH_TABLE_POWER);
  } else {
    capacity = (uint32_t)power4(seed[sn].weight);
  }
  xgzwrite(fp, genomemap_len[sn], sizeof(genomemap_len[0][0]) * capacity);

  // total
  uint32_t total = 0;
  uint32_t j;
  for (j = 0; j < capacity; j++){
    total += genomemap_len[sn][j];
  }
  xgzwrite(fp, &total, sizeof(uint32_t));

  // genome_map
  for (j = 0; j < capacity; j++) {
    xgzwrite(fp, (void *)genomemap[sn][j], sizeof(genomemap[0][0][0]) * genomemap_len[sn][j]);
  }

  gzclose(fp);
  return true;
}


bool load_genome_map_seed(const char *file)
{
  /*
   * This function reads the .seed.X file
   * file is the file name to read
   *
   * The file format is a gziped binary format as follows
   *
   * uint32_t				: shrimp_mode
   * uint32_t				: Hflag
   * seed_type			: Seed
   * uint32_t				: capacity
   * uint32_t * capacity	: genomemap_len
   * uint32_t				: total (= sum from 0 to capacity - 1 of genomemap_len)
   * uint32_t * total		: genomemap (each entry of length genomemap_len)
   *
   */
  int i;
  uint32_t j;
  uint32_t total;

  gzFile fp = gzopen(file, "rb");
  if (fp == NULL){
    fprintf(stderr,"Could not open file [%s]\n",file);
    return false;
  }

  // shrimp_mode
  uint32_t m;
  xgzread(fp, &m, sizeof(uint32_t));
  if(m != (uint32_t)shrimp_mode) {
    fprintf(stderr,"Shrimp mode in file %s does not match\n",file);
  }

  // Hflag
  uint32_t h;
  xgzread(fp, &h, sizeof(uint32_t));
  if (h != (uint32_t)Hflag){
    fprintf(stderr,"Hash settings do not match in file %s\n",file);
  }

  // Seed
  int sn = n_seeds;
  n_seeds++;
  seed = (seed_type *)xrealloc(seed, sizeof(seed_type) * n_seeds);
  genomemap_len = (uint32_t **)xrealloc_c(genomemap_len,
					  sizeof(genomemap_len[0]) * n_seeds,
					  sizeof(genomemap_len[0]) * (n_seeds - 1),
					  &mem_genomemap);
  genomemap = (uint32_t ***)xrealloc_c(genomemap,
				       sizeof(genomemap[0]) * n_seeds,
				       sizeof(genomemap[0]) * (n_seeds - 1),
				       &mem_genomemap);
  xgzread(fp,seed + sn,sizeof(seed_type));

  max_seed_span = MAX(max_seed_span, seed[sn].span);
  min_seed_span = MIN(min_seed_span, seed[sn].span);

  avg_seed_span = 0;
  for(i = 0; i < n_seeds; i++) {
    avg_seed_span += seed[i].span;
  }
  avg_seed_span = avg_seed_span/n_seeds;

  // genomemap_len
  uint32_t capacity;
  if(Hflag) {
    capacity = (uint32_t)power4(HASH_TABLE_POWER);
  } else {
    capacity = (uint32_t)power4(seed[sn].weight);
  }
  genomemap_len[sn] = (uint32_t *)xmalloc_c(sizeof(genomemap_len[0][0]) * capacity, &mem_genomemap);
  genomemap[sn] = (uint32_t **)xmalloc_c(sizeof(genomemap[0][0]) * capacity, &mem_genomemap);
  xgzread(fp, genomemap_len[sn], sizeof(uint32_t) * capacity);

  // total
  xgzread(fp, &total, sizeof(uint32_t));

  // genome_map
  uint32_t * map;
  map = (uint32_t *)xmalloc_c(sizeof(uint32_t) * total, &mem_genomemap);
  xgzread(fp,map,sizeof(uint32_t) * total);
  uint32_t * ptr;
  ptr = map;

  for (j = 0; j < capacity; j++) {
    genomemap[sn][j] = ptr;
    ptr += genomemap_len[sn][j];
  }

  gzclose(fp);
  return true;
}


bool save_genome_map(const char *prefix)
{
  /*
   * This function writes the .genome and .seed.X files appropriate for loading
   * prefix will be used to name the files eg. prefix.genome, prefix.seed.0, ...
   *
   * The file format for .seed.X files is descriped in save_genome_map_seed
   *
   * The file format for the .genome file is a gziped binary format as follows
   *
   * uint32_t					: shrimp_mode
   * uint32_t					: Hflag
   * uint32_t 				: num_contigs
   * uint32_t * num_contigs	: genome_len (the length of each contig)
   * uint32_t * num_contigs	: contig_offsets
   * per contig
   * 		uint32_t					: name_length
   * 		char * (name_length + 1)	: name including null termination
   * uint32_t					: total (= sum of BPTO32BW(genome_len)
   * per contig
   * 		uint32_t * BPTO32BW(contig_len)	: genome_contigs
   * per contig
   * 		uint32_t * BPTO32BW(contig_len)	: genome_contigs_rc
   * if colour space
   * 		per contig
   * 			uint32_t * BPTO32BW(contig_len) : genome_cs_cntigs
   * 		uint32_t * num_contigs	: genome_initbp
   */
  char * name;
  name = (char *)xmalloc(strlen(prefix) + n_seeds + 10);

  int sn;
  for(sn = 0; sn < n_seeds; sn++) {
    sprintf(name,"%s.seed.%d", prefix, sn);
    save_genome_map_seed(name, sn);
  }

  sprintf(name, "%s.genome", prefix);
  gzFile fp = gzopen(name, "wb");
  if (fp == NULL){
    return false;
  }

  //shrimp mode
  uint32_t m;
  m = (uint32_t)(shrimp_mode);
  xgzwrite(fp,&m,sizeof(uint32_t));

  //Hflag
  uint32_t h = (uint32_t)Hflag;
  xgzwrite(fp,&h,sizeof(uint32_t));

  // num contigs
  xgzwrite(fp,&num_contigs,sizeof(uint32_t));

  // genome_len
  xgzwrite(fp,genome_len,sizeof(uint32_t)*num_contigs);

  // contig_offsets
  xgzwrite(fp,contig_offsets,sizeof(uint32_t)*num_contigs);

  //names / total
  int i;
  uint32_t total = 0;
  for(i = 0; i < num_contigs; i++){
    uint32_t len = (uint32_t)strlen(contig_names[i]);
    xgzwrite(fp, &len, sizeof(uint32_t));
    xgzwrite(fp, contig_names[i], len + 1);
    total += BPTO32BW(genome_len[i]);
  }
  xgzwrite(fp,&total,sizeof(uint32_t));

  for (i = 0; i < num_contigs; i++) {
    xgzwrite(fp, (void *)genome_contigs[i], BPTO32BW(genome_len[i]) * sizeof(uint32_t));
  }
  for (i = 0; i < num_contigs; i++) {
    xgzwrite(fp, (void *)genome_contigs_rc[i], BPTO32BW(genome_len[i]) * sizeof(uint32_t));
  }
  if (shrimp_mode == MODE_COLOUR_SPACE) {
    for (i = 0; i < num_contigs; i++) {
      xgzwrite(fp, (void *)genome_cs_contigs[i], BPTO32BW(genome_len[i]) * sizeof(uint32_t));
    }
    xgzwrite(fp, (void *)genome_initbp, num_contigs * sizeof(uint32_t));
  }

  gzclose(fp);
  return true;
}


bool load_genome_map(const char *file)
{
  /*
   * This function loads the .genome
   * file is the name of the file to load from
   *
   * The file format for the .genome file is a gziped binary format as follows
   *
   * uint32_t					: shrimp_mode
   * uint32_t					: Hflag
   * uint32_t 				: num_contigs
   * uint32_t * num_contigs	: genome_len (the length of each contig)
   * uint32_t * num_contigs	: contig_offsets
   * per contig
   * 		uint32_t					: name_length
   * 		char * (name_length + 1)	: name including null termination
   * uint32_t					: total (= sum of BPTO32BW(genome_len)
   * per contig
   * 		uint32_t * BPTO32BW(contig_len)	: genome_contigs
   * per contig
   * 		uint32_t * BPTO32BW(contig_len)	: genome_contigs_rc
   * if colour space
   * 		per contig
   * 			uint32_t * BPTO32BW(contig_len) : genome_cs_cntigs
   * 		uint32_t * num_contigs	: genome_initbp
   */
  int i;

  gzFile fp = gzopen(file, "rb");
  if (fp == NULL){
    return false;
  }

  //shrimp mode
  uint32_t m;
  xgzread(fp, &m, sizeof(uint32_t));
  if (shrimp_mode != (shrimp_mode_t)m) {
    fprintf(stderr, "error: shrimp mode does not match genome file (%s)\n", file);
    exit(1);
  }

  //Hflag
  uint32_t h;
  xgzread(fp, &h, sizeof(uint32_t));
  Hflag = h;

  // num_contigs
  xgzread(fp, &num_contigs, sizeof(uint32_t));

  genome_len = (uint32_t *)xmalloc(sizeof(uint32_t) * num_contigs);
  contig_offsets = (uint32_t *)xmalloc(sizeof(uint32_t) * num_contigs);
  contig_names = (char **)xmalloc(sizeof(char *) * num_contigs);

  genome_contigs = (uint32_t **)xmalloc(sizeof(uint32_t *) * num_contigs);
  genome_contigs_rc = (uint32_t **)xmalloc(sizeof(uint32_t *) * num_contigs);
  if (shrimp_mode == MODE_COLOUR_SPACE) {
    genome_cs_contigs = (uint32_t **)xmalloc(sizeof(genome_cs_contigs[0]) * num_contigs);
    genome_initbp = (int *)xmalloc(sizeof(genome_initbp[0]) * num_contigs);
  }

  //genome_len
  xgzread(fp, genome_len, sizeof(uint32_t) * num_contigs);

  // contig_offfsets
  xgzread(fp, contig_offsets, sizeof(uint32_t) * num_contigs);

  // names / total
  for (i = 0; i < num_contigs; i++) {
    uint32_t len;
    xgzread(fp, &len,sizeof(uint32_t));
    contig_names[i] = (char *)xmalloc(sizeof(char) * (len + 1));
    xgzread(fp, contig_names[i], len + 1);
  }

  uint32_t total;
  xgzread(fp, &total, sizeof(uint32_t));

  //genome_contigs / genome_contigs_rc / genome_cs_contigs / genome_initbp
  uint32_t *gen, *gen_rc, *gen_cs, *ptr1, *ptr2, *ptr3 = NULL;
  gen = (uint32_t *)xmalloc(sizeof(uint32_t) * total);
  xgzread(fp, gen, sizeof(uint32_t) * total);
  ptr1 = gen;
  gen_rc = (uint32_t *)xmalloc(sizeof(uint32_t) * total);
  xgzread(fp, gen_rc, sizeof(uint32_t) * total);
  ptr2 = gen_rc;
  if (shrimp_mode == MODE_COLOUR_SPACE) {
    gen_cs = (uint32_t *)xmalloc(sizeof(uint32_t) * total);
    xgzread(fp, gen_cs, sizeof(uint32_t) * total);
    ptr3 = gen_cs;
    xgzread(fp, genome_initbp, sizeof(uint32_t) * num_contigs);
  }
  for (i = 0; i < num_contigs; i++) {
    genome_contigs[i] = ptr1;
    ptr1 += BPTO32BW(genome_len[i]);
    genome_contigs_rc[i] = ptr2;
    ptr2 += BPTO32BW(genome_len[i]);
    if (shrimp_mode == MODE_COLOUR_SPACE) {
      genome_cs_contigs[i] = ptr3;
      ptr3 += BPTO32BW(genome_len[i]);
    }
  }

  if (shrimp_mode == MODE_COLOUR_SPACE) {
    genome_cs_contigs_rc = (uint32_t **)xmalloc(sizeof(genome_cs_contigs_rc[0]) * num_contigs);
    for (i = 0; i < num_contigs; i++) {
      genome_cs_contigs_rc[i] = bitfield_to_colourspace(genome_contigs_rc[i], genome_len[i], false);
    }
  }

  gzclose(fp);
  return true;
}


void print_genomemap_stats()
{
  stat_t list_size, list_size_non0;
  int sn;
  uint64_t capacity, mapidx;
  uint max;

  uint64_t histogram[100];
  uint64_t cummulative_histogram[100];
  int bucket_size;
  int i, bucket;


  fprintf(stderr, "Genome Map stats:\n");

  for (sn = 0; sn < n_seeds; sn++) {
    if (Hflag)
      capacity = power(4, HASH_TABLE_POWER);
    else
      capacity = power(4, seed[sn].weight);

    stat_init(&list_size);
    stat_init(&list_size_non0);
    max = 0;
    for (mapidx = 0; mapidx < capacity; mapidx++) {
      if (genomemap_len[sn][mapidx] > list_cutoff) {
	stat_add(&list_size, 0);
	continue;
      }

      stat_add(&list_size, genomemap_len[sn][mapidx]);
      if (genomemap_len[sn][mapidx] > 0)
	stat_add(&list_size_non0, genomemap_len[sn][mapidx]);

      if (genomemap_len[sn][mapidx] > max)
	max = genomemap_len[sn][mapidx];
    }

    fprintf(stderr, "sn:%d weight:%d total_kmers:%llu lists:%llu (non-zero:%llu) list_sz_avg:%.2f (%.2f) list_sz_stddev:%.2f (%.2f) max:%u\n",
	    sn, seed[sn].weight, (long long unsigned int)stat_get_sum(&list_size),
	    (long long unsigned int)capacity, (long long unsigned int)stat_get_count(&list_size_non0),
	    stat_get_mean(&list_size), stat_get_mean(&list_size_non0),
	    stat_get_sample_stddev(&list_size), stat_get_sample_stddev(&list_size_non0), max);

    for (i = 0; i < 100; i++) {
      histogram[i] = 0;
    }

    bucket_size = ceil_div((max+1), 100); // values in [0..max]
    for (mapidx = 0; mapidx < capacity; mapidx++) {
      if (genomemap_len[sn][mapidx] > list_cutoff) {
	bucket = 0;
      } else {
	bucket = genomemap_len[sn][mapidx] / bucket_size;
	if (bucket >= 100)
	  bucket = 99;
      }
      histogram[bucket]++;
    }

    cummulative_histogram[0] = histogram[0];
    for (i = 1; i < 100; i++) {
      cummulative_histogram[i] = histogram[i] + cummulative_histogram[i-1];
    }

    for (i = 0; i < 100; i++) {
      fprintf(stderr, "[%d-%d]: %llu (cummulative: %.4f%%)\n", i*bucket_size, (i+1)*bucket_size,
	      (long long unsigned int)histogram[i], (((double)cummulative_histogram[i])/capacity)*100.0);
    }

  }
}


void free_genome(void)
{
  int sn, capacity;
  for (sn = 0; sn < n_seeds; sn++){
    if (Hflag) {
      capacity = power(4, HASH_TABLE_POWER);
    } else {
      capacity = power(4, seed[sn].weight);
    }
    //uint32_t mapidx = kmer_to_mapidx(kmerWindow, sn);
    int i;
    for (i = 0; i < capacity; i++) {
      if (genomemap[sn][i] != NULL && (i == 0 || load_file == NULL)) {
	free(genomemap[sn][i]);
      }
    }	
    free(genomemap[sn]);
    free(genomemap_len[sn]);
  }
  free(genomemap);
  free(genomemap_len);
}


/*
 * index the kmers in the genome contained in the file.
 * This can then be used to align reads against.
 */
bool load_genome(char **files, int nfiles)
{
  fasta_t fasta;
  size_t seqlen, capacity;
  uint32_t *read;
  char *seq, *name;
  uint32_t *kmerWindow;
  int sn;
  char *file;
  bool is_rna;

  //allocate memory for the genome map
  genomemap = (uint32_t ***) xmalloc_c(n_seeds * sizeof(genomemap[0]), &mem_genomemap);
  genomemap_len = (uint32_t **)xmalloc_c(n_seeds * sizeof(genomemap_len[0]), &mem_genomemap);

  for (sn = 0; sn < n_seeds; sn++) {
    capacity = power(4, Hflag? HASH_TABLE_POWER : seed[sn].weight);

    genomemap[sn] = (uint32_t **)xcalloc_c(sizeof(uint32_t *) * capacity, &mem_genomemap);
    memset(genomemap[sn],0,sizeof(uint32_t *) * capacity);
    genomemap_len[sn] = (uint32_t *)xcalloc_c(sizeof(uint32_t) * capacity, &mem_genomemap);

  }
  num_contigs = 0;
  u_int i = 0;
  int cfile;
  for(cfile = 0; cfile < nfiles; cfile++){
    file = files[cfile];
    //open the fasta file and check for errors
    fasta = fasta_open(file, MODE_LETTER_SPACE, false);
    if (fasta == NULL) {
      fprintf(stderr,"error: failded to open genome file [%s]\n",file);
      return (false);
    } else {
      fprintf(stderr,"- Processing genome file [%s]\n",file);
    }

    //Read the contigs and record their sizes

    while(fasta_get_next_contig(fasta, &name, &seq, &is_rna)){
      genome_is_rna = is_rna;
      num_contigs++;
      contig_offsets = (uint32_t *)xrealloc(contig_offsets,sizeof(uint32_t)*num_contigs);
      contig_offsets[num_contigs - 1] = i;
      contig_names = (char **)xrealloc(contig_names,sizeof(char *)*num_contigs);
      contig_names[num_contigs - 1] = name;

      fprintf(stderr,"- Processing contig %s\n",name);

      if (strchr(name, '\t') != NULL || strchr(seq, '\t') != NULL) {
	fprintf(stderr, "error: tabs are not permitted in fasta names "
		"or sequences. Tag: [%s].\n", name);
	return false;
      }

      seqlen = strlen(seq);
      if (seqlen == 0) {
	fprintf(stderr, "error: genome [%s] had no sequence!\n",
		name);
	return false;
      }

      read = fasta_sequence_to_bitfield(fasta,seq);

      if (read == NULL) {
	fprintf(stderr, "error: invalid sequence; tag: [%s]\n", name);
	return false;
      }
      genome_contigs = (uint32_t **)xrealloc(genome_contigs,sizeof(uint32_t *)*num_contigs);
      genome_contigs[num_contigs-1] = read;
      genome_contigs_rc = (uint32_t **)xrealloc(genome_contigs_rc,sizeof(uint32_t *)*num_contigs);
      genome_contigs_rc[num_contigs-1] = reverse_complement_read_ls(read,seqlen,is_rna);
      if (shrimp_mode == MODE_COLOUR_SPACE){
	genome_cs_contigs = (uint32_t **)xrealloc(genome_cs_contigs,sizeof(uint32_t *)*num_contigs);
	genome_cs_contigs[num_contigs-1] = fasta_bitfield_to_colourspace(fasta,genome_contigs[num_contigs -1],seqlen,is_rna);
	genome_cs_contigs_rc = (uint32_t **)xrealloc(genome_cs_contigs_rc,sizeof(uint32_t *)*num_contigs);
	genome_cs_contigs_rc[num_contigs - 1] = fasta_bitfield_to_colourspace(fasta, genome_contigs_rc[num_contigs -1], seqlen, is_rna);
      }
      genome_len = (uint32_t *)xrealloc(genome_len,sizeof(uint32_t)*num_contigs);
      genome_len[num_contigs -1] = seqlen;

      if (shrimp_mode == MODE_COLOUR_SPACE){
	genome_initbp = (int *)xrealloc(genome_initbp,sizeof(genome_initbp[0])*num_contigs);
	genome_initbp[num_contigs - 1] = BASE_T; // EXTRACT(read,0); NOT REALLY NEEDED
	read = fasta_bitfield_to_colourspace(fasta,read,seqlen,is_rna);
      }
      kmerWindow = (uint32_t *)xcalloc(sizeof(kmerWindow[0])*BPTO32BW(max_seed_span));
      int load = 0;
      for ( ; i < seqlen + contig_offsets[num_contigs-1]; i++) {
	int base;

	base = EXTRACT(read, i - contig_offsets[num_contigs-1]);
	bitfield_prepend(kmerWindow, max_seed_span, base);

	//skip past any Ns or Xs
	if (base == BASE_N || base == BASE_X)
	  load = 0;
	else if (load < max_seed_span)
	  load++;;
	for (sn = 0; sn < n_seeds; sn++) {
	  if (load < seed[sn].span)
	    continue;

	  uint32_t mapidx = KMER_TO_MAPIDX(kmerWindow, sn);
	  //increase the match count and store the location of the match
	  genomemap_len[sn][mapidx]++;
	  genomemap[sn][mapidx] = (uint32_t *)xrealloc_c(genomemap[sn][mapidx],
							 sizeof(uint32_t) * (genomemap_len[sn][mapidx]),
							 sizeof(uint32_t) * (genomemap_len[sn][mapidx] - 1),
							 &mem_genomemap);
	  genomemap[sn][mapidx][genomemap_len[sn][mapidx] - 1] = i - seed[sn].span + 1;

	}
      }
      if (shrimp_mode==MODE_COLOUR_SPACE) {
	free(read);
      }
			
      free(seq);
      seq = NULL;
      name = NULL;

      free(kmerWindow);
    }
    fasta_close(fasta);
  }

  fprintf(stderr,"Loaded Genome\n");
  return (true);
}


/*
 * Trim long genome lists
 */
void trim_genome()
{
  uint capacity;
  int sn;
  uint32_t mapidx;

  for (sn = 0; sn < n_seeds; sn++) {
    capacity = power(4, Hflag? HASH_TABLE_POWER : seed[sn].weight);

    for (mapidx = 0; mapidx < capacity; mapidx++) {
      if (genomemap_len[sn][mapidx] > list_cutoff) {
	genomemap_len[sn][mapidx] = 0;
	if (load_file == NULL) {
	  free(genomemap[sn][mapidx]);
	} // otherwise, this memory is block-allocated
	genomemap[sn][mapidx] = NULL;
      }
    }
  }
}
