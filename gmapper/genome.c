#define _MODULE_GENOME

#include <sys/types.h> //shm_open
#include <sys/mman.h>
#include <fcntl.h>
#include "genome.h"
#include "seeds.h"

#define MMAP_ALIGN 8


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
  uint32_t capacity = (uint32_t)power4(Hflag? HASH_TABLE_POWER : seed[sn].weight);
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
  //uint32_t total;

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
  seed = (seed_type *)
    //xrealloc(seed, sizeof(seed_type) * n_seeds);
    my_realloc(seed, sizeof(seed_type) * n_seeds, (n_seeds - 1) * sizeof(seed_type),
	       &mem_small, "seed");
  genomemap_len = (uint32_t **)
    //xrealloc_c(genomemap_len, sizeof(genomemap_len[0]) * n_seeds, sizeof(genomemap_len[0]) * (n_seeds - 1), &mem_genomemap);
    my_realloc(genomemap_len, sizeof(genomemap_len[0]) * n_seeds, sizeof(genomemap_len[0]) * (n_seeds - 1),
	       &mem_genomemap, "genomemap_len");
  genomemap = (uint32_t ***)
    //xrealloc_c(genomemap, sizeof(genomemap[0]) * n_seeds, sizeof(genomemap[0]) * (n_seeds - 1), &mem_genomemap);
    my_realloc(genomemap, sizeof(genomemap[0]) * n_seeds, sizeof(genomemap[0]) * (n_seeds - 1),
	       &mem_genomemap, "genomemap");
  genomemap_block = (ptr_and_sz *)
    my_realloc(genomemap_block, n_seeds * sizeof(genomemap_block[0]), (n_seeds - 1) * sizeof(genomemap_block[0]),
	       &mem_genomemap, "genomemap_block");

  xgzread(fp,seed + sn,sizeof(seed_type));
  max_seed_span = MAX(max_seed_span, seed[sn].span);
  min_seed_span = MIN(min_seed_span, seed[sn].span);
  avg_seed_span = 0;
  for(i = 0; i < n_seeds; i++) {
    avg_seed_span += seed[i].span;
  }
  avg_seed_span = avg_seed_span/n_seeds;

  // genomemap_len
  uint32_t capacity = (uint32_t)power4(Hflag? HASH_TABLE_POWER : seed[sn].weight);
  genomemap_len[sn] = (uint32_t *)
    //xmalloc_c(sizeof(genomemap_len[0][0]) * capacity, &mem_genomemap);
    my_malloc(sizeof(genomemap_len[0][0]) * capacity,
	      &mem_genomemap, "genomemap_len[%d]", sn);
  genomemap[sn] = (uint32_t **)
    //xmalloc_c(sizeof(genomemap[0][0]) * capacity, &mem_genomemap);
    my_malloc(sizeof(genomemap[0][0]) * capacity,
	      &mem_genomemap, "genomemap[%d]", sn);
  xgzread(fp, genomemap_len[sn], sizeof(uint32_t) * capacity);

  // total
  {
    uint32_t total;
    xgzread(fp, &total, sizeof(uint32_t));
    genomemap_block[sn].sz = (size_t)total * sizeof(uint32_t);
  }

  // genome_map
  //uint32_t * map;
  genomemap_block[sn].ptr =
    //xmalloc_c(sizeof(uint32_t) * total, &mem_genomemap);
    my_malloc(genomemap_block[sn].sz,
	      &mem_genomemap, "genomemap_block[%d].ptr", sn);
  xgzread(fp, genomemap_block[sn].ptr, genomemap_block[sn].sz);
  uint32_t * ptr;
  ptr = (uint32_t *)genomemap_block[sn].ptr;

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
  //char * name;
  //name = (char *)xmalloc(strlen(prefix) + n_seeds + 10);
  char name[strlen(prefix) + n_seeds + 10];

  int sn;
  for(sn = 0; sn < n_seeds; sn++) {
    sprintf(name,"%s.seed.%d", prefix, sn);
    save_genome_map_seed(name, sn);
  }

  sprintf(name, "%s.genome", prefix);
  gzFile fp = gzopen(name, "wb");
  if (fp == NULL){
    fprintf(stderr, "error: could not open genome file: %s.genome\n", prefix);
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
    /*
    xgzwrite(fp, (void *)genome_initbp, num_contigs * sizeof(uint32_t));
    */
  }

  gzclose(fp);
  return true;
}


static size_t
up_align(size_t size)
{
  return (size % MMAP_ALIGN == 0? size : ((size / MMAP_ALIGN) + 1) * MMAP_ALIGN);
}

static void
add_to_mmap(char * addr, char * * crt_end, size_t size, char * src = NULL)
{
  *(char **)addr = *crt_end;
  *crt_end += up_align(size);
  if (src != NULL) {
    memcpy(*(char **)addr, src, size);
  }
}


bool
genome_load_map_save_mmap(char * map_name, char const * mmap_name)
{
  map_header * h;
  int shm_fd;
  size_t map_size, capacity;

  gzFile genome_file;
  gzFile * seed_file = NULL;
  gzFile new_seed_file;
  int cn, sn;

  if (strchr(load_file, ',') == NULL) {
    // prefix: look for .genome, .seed.*
    int buf_size = strlen(map_name) + 20;
    char genome_file_name[buf_size];
    strncpy(genome_file_name, map_name, buf_size);
    strncat(genome_file_name, ".genome", buf_size);
    genome_file = gzopen(genome_file_name, "rb");
    if (genome_file == NULL) {
      crash(1, 1, "could not open genome file: %s", genome_file_name);
    } else {
      fprintf(stderr, "Loading genome from %s\n", genome_file_name);
    }

    n_seeds = 0;
    char seed_file_name[buf_size];
    snprintf(seed_file_name, buf_size, "%s.seed.%d", map_name, n_seeds);
    new_seed_file = gzopen(seed_file_name, "rb");
    while (new_seed_file != NULL) {
      n_seeds++;
      fprintf(stderr, "Loading seed from %s\n", seed_file_name);
      seed_file = (gzFile *)
	my_realloc(seed_file, n_seeds * sizeof(seed_file[0]), (n_seeds - 1) * sizeof(seed_file[0]),
		   &mem_small, "seed_file");
      seed_file[n_seeds - 1] =  new_seed_file;

      snprintf(seed_file_name, buf_size, "%s.seed.%d", map_name, n_seeds);
      new_seed_file = gzopen(seed_file_name, "rb");
    }
    if (n_seeds == 0) {
      crash(1, 0, "did not find seed file %s.seed.0", map_name);
    }
  } else {
    // comma-separated list: file.genome,file.seed.1,...
    char * c;
    c = strtok(map_name, ",");
    genome_file = gzopen(c, "rb");
    if (genome_file == NULL) {
      crash(1, 1, "could not open genome file: %s", c);
    } else {
      fprintf(stderr, "Loading genome from %s\n", c);
    }

    n_seeds = 0;
    while ((c = strtok(NULL, ",")) != NULL) {
      n_seeds++;
      seed_file = (gzFile *)
        my_realloc(seed_file, n_seeds * sizeof(seed_file[0]), (n_seeds - 1) * sizeof(seed_file[0]),
                   &mem_small, "seed_file");
      seed_file[n_seeds - 1] = gzopen(c, "rb");
      if (seed_file[n_seeds - 1] == NULL) {
	crash(1, 1, "could not open seed file: %s", c);
      } else {
	fprintf(stderr, "Loading seed from %s\n", c);
      }
    }
  }

  // at this point genome_file and seed_file are all open

  // need to load: num_contigs, contig_len, contig_names_len, seeds
  // from this, can estimate memory requirement directly
  uint32_t _shrimp_mode, _Hflag;
  xgzread(genome_file, &_shrimp_mode, sizeof(uint32_t));
  shrimp_mode = (shrimp_mode_t)_shrimp_mode;

  xgzread(genome_file, &_Hflag, sizeof(uint32_t));
  Hflag = _Hflag;
  
  xgzread(genome_file, &num_contigs, sizeof(uint32_t));

  // genome_len
  genome_len = (uint32_t *)
    my_malloc(num_contigs * sizeof(uint32_t),
              &mem_genomemap, "genome_len");
  xgzread(genome_file, genome_len, num_contigs * sizeof(uint32_t));

  // contig_offsets
  contig_offsets = (uint32_t *)
    my_malloc(num_contigs * sizeof(uint32_t),
              &mem_genomemap, "contig_offsets");
  xgzread(genome_file, contig_offsets, num_contigs * sizeof(uint32_t));

  // names / total
  contig_names = (char **)
    my_malloc(num_contigs * sizeof(contig_names[0]),
              &mem_genomemap, "contig_names");
  for (cn = 0; cn < num_contigs; cn++) {
    uint32_t len;
    xgzread(genome_file, &len, sizeof(uint32_t));
    contig_names[cn] = (char *)
      my_malloc(sizeof(char) * (len + 1),
                &mem_genomemap, "contig_names[%d]", cn);
    xgzread(genome_file, contig_names[cn], len + 1);
    assert(len == (uint32_t)strlen(contig_names[cn]));
  }

  seed = (seed_type *)
    my_malloc(n_seeds * sizeof(seed[0]),
	      &mem_small, "seed");
  for (sn = 0; sn < n_seeds; sn++) {
    xgzread(seed_file[sn], &_shrimp_mode, sizeof(uint32_t));
    if ((shrimp_mode_t)_shrimp_mode != shrimp_mode) {
      crash(1, 0, "shrimp_mode in seed file %d does not match shrimp mode from genome file", sn);
    }

    xgzread(seed_file[sn], &_Hflag, sizeof(uint32_t));
    if (Hflag != _Hflag) {
      crash(1, 0, "Hflag in seed file %d does not match Hlag from genome file", sn);
    }

    xgzread(seed_file[sn], &seed[sn], sizeof(seed[0]));
    max_seed_span = MAX(max_seed_span, seed[sn].span);
    min_seed_span = MIN(min_seed_span, seed[sn].span);
  }
  avg_seed_span = 0;
  for(sn = 0; sn < n_seeds; sn++) {
    avg_seed_span += seed[sn].span;
  }
  avg_seed_span = avg_seed_span/n_seeds;

  //
  // ready to estimate memory requirements
  //

  map_size = up_align(sizeof(struct map_header));

  // 1-dim arrays
  map_size += up_align(num_contigs * sizeof(genome_len[0]));
  map_size += up_align(num_contigs * sizeof(contig_offsets[0]));
  map_size += up_align(num_contigs * sizeof(contig_names[0]));
  map_size += up_align(num_contigs * sizeof(genome_contigs[0]));
  map_size += up_align(num_contigs * sizeof(genome_contigs_rc[0]));
  if (shrimp_mode == MODE_COLOUR_SPACE) {
    map_size += up_align(num_contigs * sizeof(genome_cs_contigs[0]));
    map_size += up_align(num_contigs * sizeof(genome_cs_contigs_rc[0]));
  }
  map_size += up_align(n_seeds * sizeof(seed[0]));
  map_size += up_align(n_seeds * sizeof(genomemap_len[0]));
  map_size += up_align(n_seeds * sizeof(genomemap[0]));
  if (Hflag) {
    map_size += up_align(n_seeds * sizeof(seed_hash_mask[0]));
  }

  // 2-dim arrays indexed by cn
  long long total_len = 0;
  for (cn = 0; cn < num_contigs; cn++) {
    map_size += up_align((strlen(contig_names[cn]) + 1) * sizeof(char));
    map_size += up_align(BPTO32BW(genome_len[cn]) * sizeof(uint32_t)); // genome_contigs
    map_size += up_align(BPTO32BW(genome_len[cn]) * sizeof(uint32_t)); // genome_contigs_rc
    if (shrimp_mode == MODE_COLOUR_SPACE) {
      map_size += up_align(BPTO32BW(genome_len[cn]) * sizeof(uint32_t)); // genome_cs_contigs
      map_size += up_align(BPTO32BW(genome_len[cn]) * sizeof(uint32_t)); // genome_cs_contigs_rc
    }
    total_len += genome_len[cn];
  }

  // 2-dim arrays indexed by sn
  for (sn = 0; sn < n_seeds; sn++) {
    if (Hflag) {
      map_size += up_align(BPTO32BW(max_seed_span) * sizeof(uint32_t));
    }
    capacity = power4(Hflag? HASH_TABLE_POWER : seed[sn].weight);
    map_size += up_align(capacity * sizeof(genomemap_len[0][0]));
    map_size += up_align(capacity * sizeof(genomemap[0][0]));
  }

  // for genomemap, in the worst case, each location appears once for every seed
  map_size += up_align((size_t)total_len * (size_t)n_seeds * sizeof(uint32_t));

  fprintf(stderr, "Allocating map of size: %.3gG\n", (double)map_size/(1024.0 * 1024.0 * 1024.0));

  if ((shm_fd = shm_open(mmap_name, O_CREAT | O_EXCL | O_RDWR, S_IREAD)) < 0 ) {
    crash(1, 1, "could not open mmap file %s for writing\n", mmap_name);
  }
  if (ftruncate(shm_fd, map_size) < 0) {
    crash(1, 1, "could not set size of mmap file %s to %lld", mmap_name, (long long)map_size);
  }
  if((h = (map_header *)mmap(0, map_size, (PROT_READ | PROT_WRITE), MAP_SHARED, shm_fd, 0)) == MAP_FAILED) {
    crash(1, 1, "could not mmap");
  }

  h->map_start = h;
  h->map_end = (char *)h + map_size;
  h->map_version = 1;

  h->shrimp_mode = shrimp_mode;
  h->Hflag = Hflag;
  h->num_contigs = num_contigs;
  h->n_seeds = n_seeds;
  h->min_seed_span = min_seed_span;
  h->max_seed_span = max_seed_span;
  h->avg_seed_span = avg_seed_span;

  // start copying pointers
  char * crt_end = (char *)h->map_start + up_align(sizeof(map_header));

  // genome_len: 1-dim; already loaded
  add_to_mmap((char*)&h->genome_len, &crt_end, num_contigs * sizeof(genome_len[0]), (char*)genome_len);

  // contig_offsets: 1-dim; already loaded
  add_to_mmap((char*)&h->contig_offsets, &crt_end, num_contigs * sizeof(contig_offsets[0]), (char*)contig_offsets);

  // contig_names: 2-dim, already loaded
  add_to_mmap((char*)&h->contig_names, &crt_end, num_contigs * sizeof(contig_names[0]));
  for (cn = 0; cn < num_contigs; cn++) {
    add_to_mmap((char*)&h->contig_names[cn], &crt_end, (strlen(contig_names[cn]) + 1) * sizeof(char), (char*)contig_names[cn]);
  }

  // genome_XX_contigs_YY: 2-dim; not loaded; block in file (except for _cs_rc)
  add_to_mmap((char*)&h->genome_contigs, &crt_end, num_contigs * sizeof(genome_contigs[0]));
  add_to_mmap((char*)&h->genome_contigs_rc, &crt_end, num_contigs * sizeof(genome_contigs_rc[0]));
  if (shrimp_mode == MODE_COLOUR_SPACE) {
    add_to_mmap((char*)&h->genome_cs_contigs, &crt_end, num_contigs * sizeof(genome_cs_contigs[0]));
    add_to_mmap((char*)&h->genome_cs_contigs_rc, &crt_end, num_contigs * sizeof(genome_cs_contigs_rc[0]));
  }

  // read blocks from file -- THESE ARE NOT ALIGNED on 8 bytes, only on 4!!
  uint32_t *ptr1, *ptr2, *ptr3 = NULL;
  uint32_t total;
  xgzread(genome_file, &total, sizeof(uint32_t));

  add_to_mmap((char*)&h->genome_contigs[0], &crt_end, (size_t)total * sizeof(uint32_t));
  xgzread(genome_file, h->genome_contigs[0], (size_t)total * sizeof(uint32_t));
  ptr1 = (uint32_t *)h->genome_contigs[0];

  add_to_mmap((char*)&h->genome_contigs_rc[0], &crt_end, (size_t)total * sizeof(uint32_t));
  xgzread(genome_file, h->genome_contigs_rc[0], (size_t)total * sizeof(uint32_t));
  ptr2 = (uint32_t *)h->genome_contigs_rc[0];

  if (shrimp_mode == MODE_COLOUR_SPACE) {
    add_to_mmap((char*)&h->genome_cs_contigs[0], &crt_end, (size_t)total * sizeof(uint32_t));
    xgzread(genome_file, h->genome_cs_contigs[0], (size_t)total * sizeof(uint32_t));
    ptr3 = (uint32_t *)h->genome_cs_contigs[0];
  }  
  
  for (cn = 0; cn < num_contigs; cn++) {
    h->genome_contigs[cn] = ptr1;
    ptr1 += BPTO32BW(genome_len[cn]);
    h->genome_contigs_rc[cn] = ptr2;
    ptr2 += BPTO32BW(genome_len[cn]);
    if (shrimp_mode == MODE_COLOUR_SPACE) {
      h->genome_cs_contigs[cn] = ptr3;
      ptr3 += BPTO32BW(genome_len[cn]);
    }
  }

  if (shrimp_mode == MODE_COLOUR_SPACE) {
    for (cn = 0; cn < num_contigs; cn++) {
      uint32_t * res = bitfield_to_colourspace(h->genome_contigs_rc[cn], genome_len[cn], false);
      add_to_mmap((char*)&h->genome_cs_contigs_rc[cn], &crt_end, BPTO32BW(genome_len[cn]) * sizeof(uint32_t), (char*)res);
    }
  }
  // done with per-contig data

  // next, seeds
  add_to_mmap((char*)&h->seed, &crt_end, n_seeds * sizeof(struct seed_type), (char*)seed);
  if (Hflag) {
    init_seed_hash_mask();
    add_to_mmap((char*)&h->seed_hash_mask, &crt_end, n_seeds * sizeof(seed_hash_mask[0]));
    for (sn = 0; sn < n_seeds; sn++) {
      add_to_mmap((char*)&h->seed_hash_mask[sn], &crt_end, BPTO32BW(max_seed_span) * sizeof(uint32_t), (char*)seed_hash_mask[sn]);
    }
  }

  // genomemap_len, genomemap: these are not loaded yet
  add_to_mmap((char*)&h->genomemap_len, &crt_end, n_seeds * sizeof(genomemap_len[0]));
  add_to_mmap((char*)&h->genomemap, &crt_end, n_seeds * sizeof(genomemap[0]));
  for (sn = 0; sn < n_seeds; sn++) {
    capacity = power4(Hflag? HASH_TABLE_POWER : seed[sn].weight);

    add_to_mmap((char*)&h->genomemap_len[sn], &crt_end, capacity * sizeof(genomemap_len[0][0]));
    xgzread(seed_file[sn], h->genomemap_len[sn], capacity * sizeof(genomemap_len[0][0]));

    add_to_mmap((char*)&h->genomemap[sn], &crt_end, capacity * sizeof(genomemap[0][0]));

    // genomemap is block-alloc-ed
    uint32_t total;
    xgzread(seed_file[sn], &total, sizeof(uint32_t));

    add_to_mmap((char*)&h->genomemap[sn][0], &crt_end, (size_t)total * sizeof(uint32_t));
    xgzread(seed_file[sn], h->genomemap[sn][0], (size_t)total * sizeof(uint32_t));

    uint32_t * ptr = (uint32_t *)h->genomemap[sn][0];
    for (size_t j = 0; j < capacity; j++) {
      h->genomemap[sn][j] = ptr;
      ptr += h->genomemap_len[sn][j];
    }
  }

  // DONE!!
  assert(crt_end < (char *)h->map_end);

  for (sn = 0; sn < n_seeds; sn++) {
    gzclose(seed_file[sn]);
  }
  gzclose(genome_file);
  close(shm_fd);

  fprintf(stderr, "Index successfully loaded index from [%s] to shared memory file [%s]\n", map_name, mmap_name);

  return true;
}


bool genome_load_mmap(char const * mmap_name)
{
  int shm_fd;
  map_header * h;
  void * map_start;
  void * map_end;

  if ((shm_fd = shm_open(mmap_name, O_RDONLY, S_IREAD)) < 0) {
    crash(1, 1, "could not open mmap file %s", mmap_name);
  }
  if ((h = (map_header *)mmap(NULL, sizeof(map_header), PROT_READ, MAP_PRIVATE, shm_fd, 0)) == MAP_FAILED) {
    crash(1, 1, "could not read map header from mmap file %s", mmap_name);
  }
  map_start = h->map_start;
  map_end = h->map_end;
  fprintf(stderr, "\nLoading shared memory index [%s] of size %.3gG\n", mmap_name,
	  (double)((char *)map_end - (char *)map_start)/(1024.0 * 1024.0 * 1024.0));
  munmap(h, sizeof(map_header));
  if ((h = (map_header *)mmap(map_start, (char *)map_end - (char *)map_start, PROT_READ, MAP_PRIVATE | MAP_FIXED, shm_fd, 0)) == MAP_FAILED
      || h != map_start) {
    crash(1, 1, "could not place mmap file %s at address %p", mmap_name, map_start);
  }
  close(shm_fd);

  shrimp_mode = h->shrimp_mode;
  Hflag = h->Hflag;
  num_contigs = h->num_contigs;
  n_seeds = h->n_seeds;
  min_seed_span = h->min_seed_span;
  max_seed_span = h->max_seed_span;
  avg_seed_span = h->avg_seed_span;

  genome_len = h->genome_len;
  contig_offsets = h->contig_offsets;
  contig_names = h->contig_names;

  genome_contigs = h->genome_contigs;
  genome_contigs_rc = h->genome_contigs_rc;
  genome_cs_contigs = h->genome_cs_contigs;
  genome_cs_contigs_rc = h->genome_cs_contigs_rc;

  seed = h->seed;
  seed_hash_mask = h->seed_hash_mask;

  genomemap_len = h->genomemap_len;
  genomemap = h->genomemap;

  fprintf(stderr, "Found %d contig%s:\n", num_contigs, num_contigs > 1? "s" : "");
  int cn;
  for (cn = 0; cn < num_contigs; cn++) {
    fprintf(stderr, "%s %u\n", contig_names[cn], genome_len[cn]);
  }
  fprintf(stderr, "\n");

#ifndef NDEBUG
  for (char * crt = (char *)h->map_start; crt < (char *)h->map_end; crt += 64) {
    not_used += *crt;
  }  
#endif

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

  //genome_len
  genome_len = (uint32_t *)
    //xmalloc(sizeof(uint32_t) * num_contigs);
    my_malloc(num_contigs * sizeof(uint32_t),
	      &mem_genomemap, "genome_len");
  xgzread(fp, genome_len, num_contigs * sizeof(uint32_t));

  // contig_offsets
  contig_offsets = (uint32_t *)
    //xmalloc(sizeof(uint32_t) * num_contigs);
    my_malloc(num_contigs * sizeof(uint32_t),
	      &mem_genomemap, "contig_offsets");
  xgzread(fp, contig_offsets, num_contigs * sizeof(uint32_t));

  // names / total
  contig_names = (char **)
    //xmalloc(sizeof(char *) * num_contigs);
    my_malloc(num_contigs * sizeof(contig_names[0]),
	      &mem_genomemap, "contig_names");
  for (i = 0; i < num_contigs; i++) {
    uint32_t len;
    xgzread(fp, &len, sizeof(uint32_t));
    contig_names[i] = (char *)
      //xmalloc(sizeof(char) * (len + 1));
      my_malloc(sizeof(char) * (len + 1),
		&mem_genomemap, "contig_names[%d]", i);
    xgzread(fp, contig_names[i], len + 1);
    assert(len == (uint32_t)strlen(contig_names[i]));
  }

  genome_contigs = (uint32_t **)
    //xmalloc(sizeof(uint32_t *) * num_contigs);
    my_malloc(num_contigs * sizeof(genome_contigs[0]),
	      &mem_genomemap, "genome_contigs");
  genome_contigs_rc = (uint32_t **)
    //xmalloc(sizeof(uint32_t *) * num_contigs);
    my_malloc(num_contigs * sizeof(genome_contigs_rc[0]),
	      &mem_genomemap, "genome_contigs_rc");
  if (shrimp_mode == MODE_COLOUR_SPACE) {
    genome_cs_contigs = (uint32_t **)
      //xmalloc(sizeof(genome_cs_contigs[0]) * num_contigs);
      my_malloc(num_contigs * sizeof(genome_cs_contigs[0]),
		&mem_genomemap, "genome_cs_contigs");
    genome_cs_contigs_rc = (uint32_t **)
      //xmalloc(sizeof(genome_cs_contigs_rc[0]) * num_contigs);
      my_malloc(num_contigs * sizeof(genome_cs_contigs_rc[0]),
		&mem_genomemap, "genome_cs_contigs_rc");
    /*
    genome_initbp = (int *)
      //xmalloc(sizeof(genome_initbp[0]) * num_contigs);
      my_malloc(num_contigs * sizeof(genome_initbp[0]),
	      &mem_genomemap, "genome_initbp[%d]", num_contigs);
    */
  }

  //genome_contigs / genome_contigs_rc / genome_cs_contigs / genome_initbp
  uint32_t *ptr1, *ptr2, *ptr3 = NULL;
  //total;
  {
    uint32_t total;
    xgzread(fp, &total, sizeof(uint32_t));
    genome_contigs_block.sz = (size_t)total * sizeof(uint32_t);
  }

  genome_contigs_block.ptr =
    //xmalloc(sizeof(uint32_t) * total);
    my_malloc(genome_contigs_block.sz,
	      &mem_genomemap, "genome_contigs_block.ptr");
  xgzread(fp, genome_contigs_block.ptr, genome_contigs_block.sz);
  ptr1 = (uint32_t *)genome_contigs_block.ptr;

  genome_contigs_rc_block.sz = genome_contigs_block.sz;
  genome_contigs_rc_block.ptr =
    //xmalloc(sizeof(uint32_t) * total);
    my_malloc(genome_contigs_rc_block.sz,
	      &mem_genomemap, "genome_contigs_rc_block.ptr");
  xgzread(fp, genome_contigs_rc_block.ptr, genome_contigs_rc_block.sz);
  ptr2 = (uint32_t *)genome_contigs_rc_block.ptr;

  if (shrimp_mode == MODE_COLOUR_SPACE) {
    genome_cs_contigs_block.sz = genome_contigs_block.sz;
    genome_cs_contigs_block.ptr =
      //xmalloc(sizeof(uint32_t) * total);
      my_malloc(genome_cs_contigs_block.sz,
		&mem_genomemap, "genome_cs_contigs_block.ptr");
    xgzread(fp, genome_cs_contigs_block.ptr, genome_cs_contigs_block.sz);
    ptr3 = (uint32_t *)genome_cs_contigs_block.ptr;

    /*
    xgzread(fp, genome_initbp, sizeof(uint32_t) * num_contigs);
    */
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
    capacity = (uint32_t)power4(Hflag? HASH_TABLE_POWER : seed[sn].weight);

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
  int i, sn;
  uint32_t j, capacity;

  for (sn = 0; sn < n_seeds; sn++){
    capacity = (uint32_t)power4(Hflag? HASH_TABLE_POWER : seed[sn].weight);
    //uint32_t mapidx = kmer_to_mapidx(kmerWindow, sn);
    if (load_file != NULL) {
      my_free(genomemap_block[sn].ptr, genomemap_block[sn].sz,
	      &mem_genomemap, "genomemap_block[%d].ptr", sn);
    } else {
      for (j = 0; j < capacity; j++) {
	if (genomemap[sn][j] != NULL) {
	  //free(genomemap[sn][j]);
	  my_free(genomemap[sn][j], genomemap_len[sn][j] * sizeof(genomemap[0][0][0]),
		  &mem_genomemap, "genomemap[%d][%u]", sn, j);
	}
      }
    }
    //free(genomemap[sn]);
    my_free(genomemap[sn], capacity * sizeof(genomemap[0][0]),
	    &mem_genomemap, "genomemap[%d]", sn);
    //free(genomemap_len[sn]);
    my_free(genomemap_len[sn], capacity * sizeof(genomemap_len[0][0]),
	    &mem_genomemap, "genomemap_len[%d]", sn);
  }
  //free(genomemap);
  my_free(genomemap, n_seeds * sizeof(genomemap[0]),
	  &mem_genomemap, "genomemap");
  //free(genomemap_len);
  my_free(genomemap_len, n_seeds * sizeof(genomemap_len[0]),
	  &mem_genomemap, "genomemap_len");
  if (load_file != NULL)
    my_free(genomemap_block, n_seeds * sizeof(genomemap_block[0]),
	    &mem_genomemap, "genomemap_block");

  if (load_file != NULL) {
    my_free(genome_contigs_block.ptr, genome_contigs_block.sz,
	    &mem_genomemap, "genome_contigs_block.ptr");
    my_free(genome_contigs_rc_block.ptr, genome_contigs_rc_block.sz,
	    &mem_genomemap, "genome_contigs_rc_block.ptr");
    if (shrimp_mode == MODE_COLOUR_SPACE) {
      my_free(genome_cs_contigs_block.ptr, genome_cs_contigs_block.sz,
	      &mem_genomemap, "genome_cs_contigs_block.ptr");
    }
  } else {
    for (i = 0; i < num_contigs; i++) {
      free(genome_contigs[i]);
      free(genome_contigs_rc[i]);
      if (shrimp_mode == MODE_COLOUR_SPACE) {
	free(genome_cs_contigs[i]);
      }
    }
  }
  if (shrimp_mode == MODE_COLOUR_SPACE) {
    for (i = 0; i < num_contigs; i++) {
      free(genome_cs_contigs_rc[i]);
    }
  }

  // genome_contigs, genome_contigs_rc, genome_cs_contigs, genome_cs_contigs_rc
  //free(genome_contigs);
  my_free(genome_contigs, num_contigs * sizeof(genome_contigs[0]),
	  &mem_genomemap, "genome_contigs");
  //free(genome_contigs_rc);
  my_free(genome_contigs_rc, num_contigs * sizeof(genome_contigs_rc[0]),
	  &mem_genomemap, "genome_contigs_rc");
  if (shrimp_mode==MODE_COLOUR_SPACE) {
    //free(genome_cs_contigs);
    my_free(genome_cs_contigs, num_contigs * sizeof(genome_cs_contigs[0]),
	    &mem_genomemap, "genome_cs_contigs");
    //free(genome_cs_contigs_rc);
    my_free(genome_cs_contigs_rc, num_contigs * sizeof(genome_cs_contigs_rc[0]),
	    &mem_genomemap, "genome_cs_contigs_rc");
    /*
    free(genome_initbp);
    */
  }

  // genome_len, contig_offsets
  //free(genome_len);
  my_free(genome_len, num_contigs * sizeof(uint32_t),
	  &mem_genomemap, "genome_len");
  //free(contig_offsets);
  my_free(contig_offsets, num_contigs * sizeof(uint32_t),
	  &mem_genomemap, "contig_offsets");

  // contig_names
  for (i = 0; i < num_contigs; i++) {
    if (load_file != NULL) {
      my_free(contig_names[i], (strlen(contig_names[i]) + 1) * sizeof(char),
	      &mem_genomemap, "contig_names[%d]", i);
    } else { // HACK: if no loadfile, names are alloc-ed from fasta_get_next..
      free(contig_names[i]);
    }
  }
  //free(contig_names);
  my_free(contig_names, num_contigs * sizeof(contig_names[0]),
	  &mem_genomemap, "contig_names");
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
  //uint32_t *kmerWindow;
  int sn;
  char *file;
  bool is_rna;

  //allocate memory for the genome map
  genomemap = (uint32_t ***)
    //xmalloc_c(n_seeds * sizeof(genomemap[0]), &mem_genomemap);
    my_malloc(n_seeds * sizeof(genomemap[0]),
	      &mem_genomemap, "genomemap");
  genomemap_len = (uint32_t **)
    //xmalloc_c(n_seeds * sizeof(genomemap_len[0]), &mem_genomemap);
    my_malloc(n_seeds * sizeof(genomemap_len[0]),
	      &mem_genomemap, "genomemap_len");

  for (sn = 0; sn < n_seeds; sn++) {
    capacity = (uint32_t)power4(Hflag? HASH_TABLE_POWER : seed[sn].weight);

    genomemap[sn] = (uint32_t **)
      //xcalloc_c(sizeof(uint32_t *) * capacity, &mem_genomemap);
      my_calloc(sizeof(uint32_t *) * capacity,
		&mem_genomemap, "genomemap[%d]", sn);
    //memset(genomemap[sn],0,sizeof(uint32_t *) * capacity); //???
    genomemap_len[sn] = (uint32_t *)
      //xcalloc_c(sizeof(uint32_t) * capacity, &mem_genomemap);
      my_calloc(sizeof(uint32_t) * capacity,
		&mem_genomemap, "genomemap_len[%d]", sn);

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
      contig_offsets = (uint32_t *)
	//xrealloc(contig_offsets,sizeof(uint32_t)*num_contigs);
	my_realloc(contig_offsets, num_contigs * sizeof(uint32_t), (num_contigs - 1) * sizeof(uint32_t),
		   &mem_genomemap, "contig_offsets");
      contig_offsets[num_contigs - 1] = i;
      contig_names = (char **)
	//xrealloc(contig_names,sizeof(char *)*num_contigs);
	my_realloc(contig_names, num_contigs * sizeof(contig_names[0]), (num_contigs - 1) * sizeof(contig_names[0]),
		   &mem_genomemap, "contig_names");
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
      genome_contigs = (uint32_t **)
	//xrealloc(genome_contigs,sizeof(uint32_t *)*num_contigs);
	my_realloc(genome_contigs, num_contigs * sizeof(genome_contigs[0]), (num_contigs - 1) * sizeof(genome_contigs[0]),
		   &mem_genomemap, "genome_contigs");
      genome_contigs[num_contigs-1] = read;
      genome_contigs_rc = (uint32_t **)
	//xrealloc(genome_contigs_rc,sizeof(uint32_t *)*num_contigs);
	my_realloc(genome_contigs_rc, sizeof(uint32_t *) * num_contigs, (num_contigs - 1) * sizeof(uint32_t *),
		   &mem_genomemap, "genome_contigs_rc");
      genome_contigs_rc[num_contigs-1] = reverse_complement_read_ls(read,seqlen,is_rna);
      if (shrimp_mode == MODE_COLOUR_SPACE){
	genome_cs_contigs = (uint32_t **)
	  //xrealloc(genome_cs_contigs,sizeof(uint32_t *)*num_contigs);
	  my_realloc(genome_cs_contigs, sizeof(uint32_t *) * num_contigs, (num_contigs - 1) * sizeof(uint32_t *),
		     &mem_genomemap, "genome_cs_contigs");
	genome_cs_contigs[num_contigs-1] = fasta_bitfield_to_colourspace(fasta,genome_contigs[num_contigs -1],seqlen,is_rna);
	genome_cs_contigs_rc = (uint32_t **)
	  //xrealloc(genome_cs_contigs_rc,sizeof(uint32_t *)*num_contigs);
	  my_realloc(genome_cs_contigs_rc, sizeof(uint32_t *) * num_contigs, (num_contigs - 1) * sizeof(uint32_t *),
		     &mem_genomemap, "genome_cs_contigs_rc");
	genome_cs_contigs_rc[num_contigs - 1] = fasta_bitfield_to_colourspace(fasta, genome_contigs_rc[num_contigs -1], seqlen, is_rna);
      }
      genome_len = (uint32_t *)
	//xrealloc(genome_len,sizeof(uint32_t)*num_contigs);
	my_realloc(genome_len, num_contigs * sizeof(uint32_t), (num_contigs - 1) * sizeof(uint32_t),
		   &mem_genomemap, "genome_len");
      genome_len[num_contigs - 1] = seqlen;

      if (shrimp_mode == MODE_COLOUR_SPACE){
	/*
	genome_initbp = (int *)
	  //xrealloc(genome_initbp,sizeof(genome_initbp[0])*num_contigs);
	  my_realloc(genome_initbp, sizeof(genome_initbp[0]) * num_contigs, (num_contigs - 1) * sizeof(genome_initbp[0]),
		     &mem_genomemap, "realloc genomeinitbp[%d]", num_contigs);
	genome_initbp[num_contigs - 1] = BASE_T; // EXTRACT(read,0); NOT REALLY NEEDED
	*/
	read = fasta_bitfield_to_colourspace(fasta,read,seqlen,is_rna);
      }

      //kmerWindow = (uint32_t *)xcalloc(sizeof(kmerWindow[0])*BPTO32BW(max_seed_span));
      uint32_t kmerWindow[BPTO32BW(max_seed_span)];

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
	  genomemap[sn][mapidx] = (uint32_t *)
	    //xrealloc_c(genomemap[sn][mapidx], sizeof(uint32_t) * (genomemap_len[sn][mapidx]), sizeof(uint32_t) * (genomemap_len[sn][mapidx] - 1), &mem_genomemap);
	    my_realloc(genomemap[sn][mapidx], sizeof(uint32_t) * (genomemap_len[sn][mapidx]), sizeof(uint32_t) * (genomemap_len[sn][mapidx] - 1),
		       &mem_genomemap, "genomemap[%d][%u]", sn, mapidx);
	  genomemap[sn][mapidx][genomemap_len[sn][mapidx] - 1] = i - seed[sn].span + 1;

	}
      }
      if (shrimp_mode==MODE_COLOUR_SPACE) {
	free(read);
      }
			
      free(seq);
      seq = NULL;
      name = NULL;

      //free(kmerWindow);
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
  int sn;
  uint32_t mapidx, capacity;

  for (sn = 0; sn < n_seeds; sn++) {
    capacity = (uint32_t)power4(Hflag? HASH_TABLE_POWER : seed[sn].weight);

    for (mapidx = 0; mapidx < capacity; mapidx++) {
      if (genomemap_len[sn][mapidx] > list_cutoff) {
	genomemap_len[sn][mapidx] = 0;
	if (load_file == NULL) {
	  //free(genomemap[sn][mapidx]);
	  my_free(genomemap[sn][mapidx], genomemap_len[sn][mapidx],
		  &mem_genomemap, "genomemap[%d][%u]", sn, mapidx);
	} // otherwise, this memory is block-allocated
	genomemap[sn][mapidx] = NULL;
      }
    }
  }
}
