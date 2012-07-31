/*
 * This file should contain extern declarations for global variables in gmapper.
 * This is the only file that should be included by other modules (seeds, mapping, etc).
 */

#ifndef _GMAPPER_H
#define _GMAPPER_H

#ifdef __cplusplus
//extern "C" {
#endif

#include "../gmapper/gmapper-definitions.h"
#include "../common/debug.h"
#include "../common/util.h"
#include "../common/time_counter.h"
#include "../common/gen-st.h"

#undef EXTERN
#undef STATIC
#ifdef _MODULE_GMAPPER
#include "../gmapper/gmapper-defaults.h"
#define EXTERN(_type, _id, _init_val) _type _id = _init_val
#define STATIC(_type, _id, _init_val) static _type _id = _init_val
#else
#define EXTERN(_type, _id, _init_val) extern _type _id
#define STATIC(_type, _id, _init_val)
#endif


/* shrimp mode */
EXTERN(shrimp_mode_t,	shrimp_mode,		DEF_SHRIMP_MODE);

EXTERN(shrimp_args_t,	shrimp_args,		{});


/* thread control */
EXTERN(int,			num_threads,		DEF_NUM_THREADS);
EXTERN(int,			chunk_size,		DEF_CHUNK_SIZE);
EXTERN(int,			not_used,		0);


/* parameters */
EXTERN(struct read_mapping_options_t *,		unpaired_mapping_options[2],	{});
EXTERN(int,					n_unpaired_mapping_options[2],	{});
EXTERN(struct readpair_mapping_options_t *,	paired_mapping_options,		NULL);
EXTERN(int,					n_paired_mapping_options,	0);

EXTERN(int,		mode_mirna,		false);
EXTERN(double,		window_len,		DEF_WINDOW_LEN);
EXTERN(double,		window_overlap,		DEF_WINDOW_OVERLAP);
EXTERN(int,		match_mode,		0);
EXTERN(int,		num_outputs,		DEF_NUM_OUTPUTS);
EXTERN(int,		max_alignments,		DEF_MAX_ALIGNMENTS);
EXTERN(int,		num_tmp_outputs,	20 + DEF_NUM_OUTPUTS);
EXTERN(int,		anchor_width,		DEF_ANCHOR_WIDTH);
EXTERN(int,		indel_taboo_len,	DEF_INDEL_TABOO_LEN);
EXTERN(uint32_t,	list_cutoff,		DEF_LIST_CUTOFF);
EXTERN(bool,		gapless_sw,		DEF_GAPLESS_SW);
EXTERN(bool,		hash_filter_calls,	DEF_HASH_FILTER_CALLS);
EXTERN(int,		longest_read_len,	DEF_LONGEST_READ_LENGTH);
EXTERN(bool,		trim,			false);
EXTERN(int,		trim_front,		0);
EXTERN(int,		trim_end,		0);
EXTERN(bool,		trim_first,		true);
EXTERN(bool,		trim_second,		true);
EXTERN(bool,		trim_illumina,		false);
EXTERN(char *,		save_file,		NULL);
EXTERN(char *,		load_file,		NULL);
EXTERN(char *,		save_mmap,		NULL);
EXTERN(char *,		load_mmap,		NULL);
EXTERN(unsigned int,	progress,		DEF_PROGRESS);

EXTERN(bool,		compute_mapping_qualities,	true);
EXTERN(bool,		no_qv_check,			false);
//EXTERN(int,		score_difference_mq_cutoff,	0);
EXTERN(bool,		all_contigs,			false);
EXTERN(bool,		use_sanger_qvs,			true);
EXTERN(int,		qual_vector_offset,		0);
EXTERN(int,		qual_delta,			33);
EXTERN(int,		min_avg_qv,			10);


/* Flags */
EXTERN(bool,		strata_flag,		false);		/* get only top scoring hits */
EXTERN(bool,		Cflag,			false);		/* do complement only */
EXTERN(bool,		Fflag,			false);		/* do positive (forward) only */
EXTERN(bool,		Hflag,			false);		/* use hash table, not lookup */
EXTERN(bool,		Pflag,			false);		/* pretty print results */
EXTERN(bool,		Rflag,			false);		/* add read sequence to output*/
EXTERN(bool,		Tflag,			true);		/* reverse sw full tie breaks */
EXTERN(bool,		Dflag,			false);		/* print statistics for each thread */
EXTERN(bool,		Eflag,			true);		/* output sam format */
EXTERN(bool,		Xflag,			false);		/* print insert histogram */
EXTERN(bool,		Yflag,			false);		/* print genome projection histogram */
EXTERN(bool,		Vflag,			true);		/* automatic genome index trimming */
EXTERN(bool,		Qflag,			true);		/* use fastq reads */
EXTERN(bool,		Gflag,			true);		/* global alignment flag ! */
EXTERN(bool,		Bflag,			false);		/* be like bfast - cs only! */
EXTERN(bool,		extra_sam_fields,	false);
EXTERN(bool,		single_best_mapping,	false);
EXTERN(bool,		improper_mappings,	true);
EXTERN(bool,		autodetect_input,	true);
EXTERN(bool,		ignore_qvs,		false);		/* if input is fastq, ignore qvs in analysis */
//EXTERN(bool,		hack,			false);

/* Scores */
EXTERN(int,		match_score,		DEF_LS_MATCH_SCORE);
EXTERN(int,		mismatch_score,		DEF_LS_MISMATCH_SCORE);
EXTERN(int,		a_gap_open_score,	DEF_LS_A_GAP_OPEN);
EXTERN(int,		a_gap_extend_score,	DEF_LS_A_GAP_EXTEND);
EXTERN(int,		b_gap_open_score,	DEF_LS_B_GAP_OPEN);
EXTERN(int,		b_gap_extend_score,	DEF_LS_B_GAP_EXTEND);
EXTERN(int,		crossover_score,	DEF_CS_XOVER_SCORE);

EXTERN(double,		score_alpha,		0.0);
EXTERN(double,		score_beta,		0.0);
EXTERN(double,		pr_mismatch,		0.0);
EXTERN(double,		pr_xover,		0.03);
EXTERN(double,		pr_del_open,		0.0);
EXTERN(double,		pr_del_extend,		0.0);
EXTERN(double,		pr_ins_open,		0.0);
EXTERN(double,		pr_ins_extend,		0.0);

EXTERN(double,		window_gen_threshold,	DEF_WINDOW_GEN_THRESHOLD);
EXTERN(double,		sw_vect_threshold,	DEF_SW_VECT_THRESHOLD);
EXTERN(double,		sw_full_threshold,	DEF_SW_FULL_THRESHOLD);


/* shrimp parameter/option parsing */
STATIC(struct option const,	standard_options[],	DEF_STANDARD_OPTIONS);
STATIC(struct option const,	colour_space_options[],	DEF_COLOUR_SPACE_OPTIONS);
STATIC(struct option const,	letter_space_options[], DEF_LETTER_SPACE_OPTIONS);
STATIC(size_t const,		standard_entries,	sizeof(standard_options)/sizeof(struct option));
STATIC(size_t const,		letter_entries,		sizeof(letter_space_options)/sizeof(struct option));
STATIC(size_t const,		colour_entries,		sizeof(colour_space_options)/sizeof(struct option));


/* pairing mode */
EXTERN(int,		pair_mode,			DEF_PAIR_MODE);
EXTERN(int,		min_insert_size,		DEF_MIN_INSERT_SIZE);
EXTERN(int,		max_insert_size,		DEF_MAX_INSERT_SIZE);
EXTERN(double,		insert_size_mean,		DEF_INSERT_SIZE_MEAN);
EXTERN(double,		insert_size_stddev,		DEF_INSERT_SIZE_STDDEV);
EXTERN(llint,		insert_histogram[100],		{});
EXTERN(int,		insert_histogram_bucket_size,	1);
EXTERN(int,		insert_histogram_load,		100);
EXTERN(char *,		reads_filename,			NULL);
EXTERN(char *,	 	left_reads_filename,		NULL);
EXTERN(char *,		right_reads_filename,		NULL);
EXTERN(bool,		single_reads_file,		true);

STATIC(char const * const,	pair_mode_string[5],		DEF_PAIR_MODE_STRING);
EXTERN(bool,			pair_reverse[5][2],		DEF_PAIR_REVERSE);


/* seed management */
EXTERN(int,			n_seeds,		0);
EXTERN(struct seed_type *,	seed,			NULL);
EXTERN(uint32_t * *,		seed_hash_mask,		NULL);
EXTERN(int,			max_seed_span,		0);
EXTERN(int,			min_seed_span,		MAX_SEED_SPAN);
EXTERN(int,			avg_seed_span,		0);


/* Thread output buffer */
EXTERN(char **,			thread_output_buffer,		NULL);
EXTERN(size_t *,		thread_output_buffer_sizes,	NULL);
EXTERN(char **,			thread_output_buffer_filled,	NULL);
EXTERN(unsigned int *,		thread_output_buffer_chunk,	NULL);
EXTERN(size_t,			thread_output_buffer_initial,	DEF_THREAD_OUTPUT_BUFFER_INITIAL);
EXTERN(size_t,			thread_output_buffer_increment,	DEF_THREAD_OUTPUT_BUFFER_INCREMENT);
EXTERN(size_t,			thread_output_buffer_safety,	DEF_THREAD_OUTPUT_BUFFER_SAFETY);
EXTERN(unsigned int,		thread_output_heap_capacity,	DEF_THREAD_OUTPUT_HEAP_CAPACITY);


/* SAM stuff */
EXTERN(FILE *,		unaligned_reads_file,		NULL);
EXTERN(FILE *,		aligned_reads_file,		NULL);
EXTERN(bool,		sam_unaligned,			false);
EXTERN(bool,		half_paired,			true); //output reads in paired mode that only have one mapping
EXTERN(bool,		sam_r2,				false);
EXTERN(char *,		sam_header_filename,		NULL);
EXTERN(char *,		sam_read_group_name,		NULL);
EXTERN(char *,		sam_sample_name,		NULL);
EXTERN(FILE *,		sam_header_hd,			NULL);
EXTERN(FILE *,		sam_header_sq,			NULL);
EXTERN(FILE *,		sam_header_rg,			NULL);
EXTERN(FILE *,		sam_header_pg,			NULL);


/* Statistics */
EXTERN(llint,			nreads,				0);
EXTERN(llint,			nreads_mod,			0);
EXTERN(llint,			total_reads_matched,		0);
EXTERN(llint,			total_pairs_matched,		0);
EXTERN(llint,			total_reads_matched_conf,	0);
EXTERN(llint,			total_pairs_matched_conf,	0);
EXTERN(llint,			total_reads_dropped,		0);
EXTERN(llint,			total_pairs_dropped,		0);
EXTERN(llint,			total_single_matches,		0);
EXTERN(llint,			total_paired_matches,		0);
EXTERN(llint,			total_dup_single_matches,	0);			/* number of duplicate hits */
EXTERN(llint,			total_dup_paired_matches,	0);

EXTERN(llint,			load_genome_usecs,		0);
EXTERN(llint,			mapping_wallclock_usecs,	0);

/* per-thread counts and statistics */
//EXTERN(llint,			read_handle_usecs,		0);
//EXTERN(llint,			wait_ticks,			0);
//EXTERN(llint,			anchor_list_ticks,		0);
//EXTERN(llint,			region_counts_ticks,		0);
//EXTERN(llint,			mp_region_counts_ticks,		0);
//EXTERN(llint,			hit_list_ticks,			0);
//EXTERN(llint,			pass1_ticks,			0);
//EXTERN(llint,			get_vector_hits_ticks,		0);
//EXTERN(llint,			pass2_ticks,			0);
//EXTERN(llint,			duplicate_removal_ticks,	0);
//EXTERN(stat_t,			anchor_list_init_size,		0);
//EXTERN(stat_t,			n_big_gaps_anchor_list,		0);
//EXTERN(stat_t,			n_anchors_discarded,		0);
EXTERN(int,			anchor_list_big_gap,		DEF_ANCHOR_LIST_BIG_GAP);

// thread-private globals
typedef struct tpg_t {
  llint read_handle_usecs;
  //llint wait_ticks;
  time_counter wait_tc;
  //llint anchor_list_ticks;
  time_counter anchor_list_tc;
  //llint region_counts_ticks;
  time_counter region_counts_tc;
  //llint mp_region_counts_ticks;
  time_counter mp_region_counts_tc;
  //llint hit_list_ticks;
  time_counter hit_list_tc;
  //llint pass1_ticks;
  time_counter pass1_tc;
  //llint get_vector_hits_ticks;
  time_counter get_vector_hits_tc;
  //llint pass2_ticks;
  time_counter pass2_tc;
  //llint duplicate_removal_ticks;
  time_counter duplicate_removal_tc;
  stat_t anchor_list_init_size;
  stat_t n_big_gaps_anchor_list;
  stat_t n_anchors_discarded;
} tpg_t;

EXTERN(tpg_t,	tpg,	{});
#pragma omp threadprivate(tpg)

EXTERN(count_t,			mem_genomemap,			{});
EXTERN(count_t,			mem_small,			{});
EXTERN(count_t,			mem_thread_buffer,		{});
EXTERN(count_t,			mem_mapping,			{});
EXTERN(count_t,			mem_sw,				{});


/* genome map */
EXTERN(uint32_t ***,		genomemap,			NULL);
EXTERN(uint32_t **,		genomemap_len,			NULL);
EXTERN(uint32_t *,		contig_offsets,			NULL);	/* offset info for genome contigs */
EXTERN(char **,			contig_names,			NULL);
EXTERN(int,			num_contigs,			0);
EXTERN(uint32_t **,		genome_contigs,			NULL);	/* genome -- always in letter */
EXTERN(uint32_t **,		genome_contigs_rc,		NULL);	/* reverse complemets */
EXTERN(uint32_t **,		genome_cs_contigs,		NULL);
EXTERN(uint32_t **,		genome_cs_contigs_rc,		NULL);
EXTERN(int *,			genome_initbp,			NULL);
EXTERN(uint32_t	*,		genome_len,			NULL);
EXTERN(bool,			genome_is_rna,			false);	/* is genome RNA (has uracil)?*/
EXTERN(long long int,		total_genome_size,		0);
EXTERN(gen_st,			contig_offsets_gen_st,		{});

EXTERN(ptr_and_sz *,		genomemap_block,		NULL);
EXTERN(ptr_and_sz,		genome_contigs_block,		{});
EXTERN(ptr_and_sz,		genome_contigs_rc_block,	{});
EXTERN(ptr_and_sz,		genome_cs_contigs_block,	{});


/* region handling */
EXTERN(bool,			use_regions,			DEF_USE_REGIONS);
EXTERN(int,			region_bits,			DEF_REGION_BITS);
EXTERN(int,			region_overlap,			DEF_REGION_OVERLAP);
EXTERN(int,			n_regions,			(1 << (32 - DEF_REGION_BITS)));

typedef uint16_t	region_map_t;
EXTERN(region_map_t *,		region_map[2][2],		{});
EXTERN(int,			region_map_id,			0);
EXTERN(int,			region_map_id_bits,		13);
//EXTERN(int,			region_map_max_count,		((1 << 8) - 1));
#pragma omp threadprivate(region_map, region_map_id)


/* contains inlined calls; uses gapless_sw and hash_filter_calls vars */
#include "../common/f1-wrapper.h"


//void		hit_free_sfrp(struct read_hit *);
void		read_free(struct read_entry *);
void		read_free_hit_list(struct read_entry *);
void		read_free_anchor_list(struct read_entry *);
void		read_free_full(struct read_entry *);


/* pulled off the web; this may or may not be any good */
static inline uint32_t
hash(uint32_t a)
{
  a = (a+0x7ed55d16) + (a<<12);
  a = (a^0xc761c23c) ^ (a>>19);
  a = (a+0x165667b1) + (a<<5);
  a = (a+0xd3a2646c) ^ (a<<9);
  a = (a+0xfd7046c5) + (a<<3);
  a = (a^0xb55a4f09) ^ (a>>16);
  return a;
}


/* hash-based version or kmer -> map index function for larger seeds */
static inline uint32_t
kmer_to_mapidx_hash(uint32_t *kmerWindow, int sn)
{
  static uint32_t maxidx = ((uint32_t)1 << 2*HASH_TABLE_POWER) - 1;
  uint32_t mapidx = 0;
  int i;

  assert(seed_hash_mask != NULL);

  for (i = 0; i < BPTO32BW(max_seed_span); i++)
    mapidx = hash((kmerWindow[i] & seed_hash_mask[sn][i]) ^ mapidx);

  return mapidx & maxidx;
}


/*
 * Compress the given kmer into an index in 'readmap' according to the seed.
 * While not optimal, this is only about 20% of the spaced seed scan time.
 *
 * This is the original version for smaller seeds.
 *
 * XXX- This algorithm only considers bases 0-3, which implies overlap
 *      when we have other bases (mainly uracil, but also wobble codes).
 *      This won't affect sensitivity, but may cause extra S-W calls.
 */
static inline uint32_t
kmer_to_mapidx_orig(uint32_t *kmerWindow, int sn)
{
  bitmap_type a = seed[sn].mask[0];
  uint32_t mapidx = 0;
  int i = 0;

  do {
    if ((a & 0x1) == 0x1) {
      mapidx <<= 2;
      mapidx |= ((kmerWindow[i/8] >> (i%8)*4) & 0x3);
    }
    a >>= 1;
    i++;
  } while (a != 0x0);

  assert(mapidx < power(4, seed[sn].weight));

  return mapidx;
}

#define KMER_TO_MAPIDX(kmer, sn) (Hflag? kmer_to_mapidx_hash((kmer), (sn)) : kmer_to_mapidx_orig((kmer), (sn)))

/* get contig number from absolute index */
static inline void
get_contig_num(uint32_t idx, int * cn) {

  if (num_contigs < 100)
    {
      *cn = 0;
      while (*cn < num_contigs - 1
	     && idx >= contig_offsets[*cn + 1])
	(*cn)++;
    }
  else
    {

  /*
  int l, r, m;

  l = 0;
  r = num_contigs;
  while (l + 1 < r) {
    m = (r + l)/2;
    if (idx < contig_offsets[m])
      r = m;
    else
      l = m;
  }
  *cn = l;
  */

      *cn = gen_st_search(&contig_offsets_gen_st, idx);
    }

  assert(contig_offsets[*cn] <= idx && idx < contig_offsets[*cn] + genome_len[*cn]);
}


#ifdef __cplusplus
//} /* extern "C" */
#endif

#endif
