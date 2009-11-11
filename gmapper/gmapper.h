#ifndef _GMAPPER_H
#define _GMAPPER_H

#include "../common/sw-full-common.h"
#include "../common/debug.h"

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

  int8_t	initbp;		/* colour space init letter */
  int8_t	initbp_rc;
  uint8_t	read_len;
  uint8_t	window_len;
  uint32_t	sw_hits;
  uint32_t	final_matches;
};

#endif

