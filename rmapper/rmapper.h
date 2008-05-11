/*	$Id$	*/

extern const bool use_colours;

/* default parameters - optimised for human */
#define DEF_SPACED_SEED_CS	"1111101111"
#define DEF_SPACED_SEED_LS	"111111011111"	/* longer for solexa/454 reads*/
#define	DEF_WINDOW_LEN		115.0		/* 115% of read length */
#define DEF_NUM_MATCHES		2
#define DEF_TABOO_LEN		4
#define DEF_NUM_OUTPUTS		100
#define DEF_MAX_READ_LEN	1000		/* high sanity mark */
#define DEF_KMER_STDDEV_LIMIT	-1		/* disabled by default */

#define DEF_MATCH_VALUE		100
#define DEF_MISMATCH_VALUE	-150
#define DEF_GAP_OPEN		-400
#define DEF_GAP_EXTEND		-70
#define DEF_XOVER_PENALTY	-140
#define DEF_SW_VECT_THRESHOLD	60.0	/* == DEF_SW_FULL_THRESHOLD in lspace */
#define DEF_SW_FULL_THRESHOLD	68.0	/* read_length x match_value x .68 */

/*
 * The maximum seed weight (maximum number of 1's in the seed) sets an
 * upper limit on our lookup table allocation size. The memory usage of
 * rmapper corresponds strongly to 4^MAX_SEED_WEIGHT * (sizeof(void *) +
 * sizeof(uint32_t)). At 16, this is 32GB on 32-bit and 48GB on 64-bit
 * architectures.
 */
#ifndef MAX_SEED_WEIGHT
#define MAX_SEED_WEIGHT		16
#endif

struct re_score {
	struct read_elem *parent;		/* associated read_elem */
	struct re_score  *next;			/* linked list */
	int		  contig_num;		/* config index (for filename)*/
	bool		  revcmpl;		/* from contig's reverse cmpl */
	int32_t		  score;		/* doubles as heap cnt in [0] */
	uint32_t	  index;
};

/*
 * Keep juicy bits of reads separate from the read structure itself. This significantly
 * reduces the cache footprint due to reads during scan() and hence speeds things up
 * enough to justify the ugliness.
 */
struct read_int {
	char		 *name;
	struct read_elem *next;			/* next in read list */
	uint32_t	 *read;			/* the read as a bitstring */
	uint32_t	  read_len;
	int		  initbp;		/* colour space init letter */

	int		  swhits;		/* num of hits with sw */
	struct re_score  *scores;		/* top 'num_ouputs' scores */
	uint32_t	  final_matches;	/* num of final output matches*/
};

struct read_elem {
	struct read_int  *ri;			/* pointer to the meat */

	/* the following are used during scan() */
	uint32_t	  last_swhit_idx;	/* index of last sw hit */
	uint16_t	  window_len;		/* per-read window length */
	uint8_t	  	  prev_hit;		/* prev index in 'hits' */
	uint8_t		  next_hit;		/* next index in 'hits' */
	struct {
		uint32_t  g_idx;		/* kmer index in genome */
		uint32_t  r_idx;		/* kmer index in read */
	} hits[0];				/* size depends on num_matches */
};

/*
 * Each index of `readmap' points to an array of `readmap_entry' structures.
 * These structures in turn point to all reads that contain the indexing kmer
 * as well as track where in the read the kmer exists for colinearity checking.
 */
struct readmap_entry {
	struct read_elem *re;		/* point to the read */
	int		  idx;		/* kmer's index in the read sequence */
};
