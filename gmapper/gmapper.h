//uncomment to debug
//#define DEBUG(x) fprintf(stderr,"DEBUG: " x "\n")
#ifndef DEBUG
#define DEBUG(X)
#endif



/*
 * Main entry for every read.
 */
struct read_entry {
  char *	name;
  uint32_t *	read;		/* the read as a bitstring */
  uint32_t * read_rev;

  int8_t	initbp;		/* colour space init letter */
  int8_t	initbp_rev;
  uint32_t	read_len;
  uint16_t  window_len;
  uint32_t	sw_hits;
};
