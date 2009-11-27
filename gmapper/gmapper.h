#ifndef _GMAPPER_H
#define _GMAPPER_H

#include "../common/sw-full-common.h"
#include "../common/debug.h"

#define DEF_NUM_THREADS 1
#define DEF_CHUNK_SIZE 1000

struct re_score {
  struct sw_full_results * sfrp;	/* alignment results (final pass) */

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

struct read_entry {
  char *	name;
  re_score *	scores;
  uint32_t *	read;		/* the read as a bitstring */
  uint32_t *	read_rc;

  uint32_t *	mapidx;		/* per-seed list of mapidxs in read */
  uint32_t *	mapidx_rc;	/* per-seed list of mapidxs in read_rc */
  bool *	mapidx_pos;	/* per-seed list of validity of mapidx positions in read; only if read has Ns */
  bool *	mapidx_rc_pos;	/* per-seed list of validity of mapidx positions in read_rc; only if read has Ns */

  int8_t	initbp;		/* colour space init letter */
  int8_t	initbp_rc;
  uint8_t	read_len;
  uint8_t	window_len;
  uint32_t	sw_hits;
  uint32_t	final_matches;

  uint8_t	max_n_kmers;	/* = read_len - min_seed_span + 1 */
  uint8_t	min_kmer_pos;	/* = 0 in LS; = 1 in CS */
  bool		has_Ns;
};

#endif

