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

  int16_t	initbp;		/* colour space init letter */
  uint32_t	read1_len;
  uint16_t  window_len;
};
