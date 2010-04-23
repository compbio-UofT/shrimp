#ifndef _GMAPPER_H
#define _GMAPPER_H


#include "../common/sw-full-common.h"
#include "../common/debug.h"
#include "../common/anchors.h"
#include "../common/heap.h"


#define DEF_NUM_THREADS 1
#define DEF_CHUNK_SIZE 1000

#define DEF_HASH_FILTER_CALLS	true
#define DEF_GAPLESS_SW		false
#define DEF_LIST_CUTOFF		4294967295u // 2^32 - 1

#define PAIR_NONE	0
#define PAIR_OPP_IN	1
#define PAIR_OPP_OUT	2
#define PAIR_COL_FW	3
#define PAIR_COL_BW	4

static char const * const pair_mode_string[5] =
  { "None", "opposing strands; inwards", "opposing strands; outwards",
    "same strand; second is forward", "same strand; second is backward" };

static bool const pair_reverse[5][2] =
  { { 0, 0 }, // PAIR_NONE
    { 0, 0 }, // PAIR_OPP_IN
    { 1, 1 }, // PAIR_OPP_OUT
    { 0, 1 }, // PAIR_COL_FW
    { 1, 0 }  // PAIR_COL_BW
  };


#define DEF_PAIR_MODE		PAIR_NONE
#define DEF_MIN_INSERT_SIZE	50
#define DEF_MAX_INSERT_SIZE	2000


struct re_score {
  struct sw_full_results * sfrp;	/* alignment results (final pass) */
  struct anchor	anchor;

  uint		contig_num;	/* contig index (for filename)*/
  union {
    int		score;		/* doubles as heap cnt in [0] */
    uint	heap_elems;
  };
  union {
    uint	g_idx;		/* doubles as heap alloc in [0]*/
    uint	heap_capacity;
  };
  bool		rev_cmpl;	/* from contig's reverse cmpl */
};

struct range_restriction {
  uint	cn;
  uint	st;
  uint	g_start;
  uint	g_end;
};

struct read_entry {
  char *	name;
  char *	seq;
  re_score *	scores;
  uint32_t *	read[2];	/* the read as a bitstring */

  uint32_t *	mapidx[2];	/* per-seed list of mapidxs in read */
  bool *	mapidx_pos[2];	/* per-seed list of validity of mapidx positions in read; only if read has Ns */

  struct uw_anchor *	anchors[2];	/* list of anchors */

  struct read_hit *	hits[2];	/* list of hits */

  struct range_restriction * ranges;
  char *		range_string;

  uint		n_anchors[2];
  uint		n_hits[2];
  uint		n_ranges;

  uint32_t	sw_hits;
  uint32_t	final_matches;

  int8_t	initbp[2];		/* colour space init letter */
  uint8_t	read_len;
  uint8_t	window_len;

  uint8_t	max_n_kmers;	/* = read_len - min_seed_span + 1 */
  uint8_t	min_kmer_pos;	/* = 0 in LS; = 1 in CS */
  uint8_t	input_strand;
  bool		has_Ns;
  bool		is_rna;
};


struct read_hit {
  struct sw_full_results *	sfrp;
  struct anchor	anchor;
  uint		g_off;
  int		score_window_gen;
  int		score_vector;
  int		pct_score_vector;
  int		score_full;
  int		pct_score_full;
  int		score_max;
  uint		matches;
  uint		cn;
  int		pair_min;
  int		pair_max;
  uint16_t	w_len;
  uint8_t	st;
  uint8_t	gen_st;
};

struct read_hit_pair_holder {
  struct read_hit *	hit[2];
  int			insert_size;
};

struct read_hit_holder {
  struct read_hit *	hit;
};

DEF_HEAP(uint, struct read_hit_holder, unpaired)
DEF_HEAP(uint, struct read_hit_pair_holder, paired)


#endif
