/*	$Id$	*/

extern const bool use_colours;
extern const bool use_dag;

/* default parameters - optimised for human */
#define DEF_SPACED_SEED_CS	"1111001111"	/* handle more adjacencies */
#define DEF_SPACED_SEED_LS	"111111011111"	/* longer for solexa/454 reads*/
#define DEF_SPACED_SEED_DAG	"11110111"	/* shorter for Helicos */ 
#define	DEF_WINDOW_LEN		115.0		/* 115% of read length */
#define DEF_NUM_MATCHES		2
#define DEF_HIT_TABOO_LEN	4
#define DEF_SEED_TABOO_LEN	0
#define DEF_NUM_OUTPUTS		100
#define DEF_MAX_READ_LEN	1000		/* high sanity mark */
#define DEF_KMER_STDDEV_LIMIT	-1		/* disabled by default */

/* DAG Scores/Parameters */
#define DEF_DAG_EPSILON		  0
#define DEF_DAG_READ_MATCH	  4
#define DEF_DAG_READ_GAP	 -2
#define DEF_DAG_READ_MISMATCH	 -4

#define DEF_DAG_REF_MATCH		11
#define DEF_DAG_REF_MISMATCH	       -10
#define DEF_DAG_REF_HALF_MATCH		4
#define DEF_DAG_REF_NEITHER_MATCH      -5
#define DEF_DAG_REF_MATCH_DELETION	5
#define DEF_DAG_REF_MISMATCH_DELETION  -6
#define DEF_DAG_REF_ERROR_INSERTION    -6
#define DEF_DAG_REF_WEIGHTED_THRESHOLD  8.0

/* SW Scores */
#define DEF_MATCH_VALUE		100
#define DEF_MATCH_VALUE_DAG	 5
#define DEF_MISMATCH_VALUE	-150
#define DEF_MISMATCH_VALUE_DAG	-6
#define DEF_A_GAP_OPEN		-400
#define DEF_A_GAP_OPEN_DAG	 0
#define DEF_B_GAP_OPEN		 DEF_A_GAP_OPEN
#define DEF_B_GAP_OPEN_DAG	 0
#define DEF_A_GAP_EXTEND	-70
#define DEF_A_GAP_EXTEND_DAG	-6
#define DEF_B_GAP_EXTEND	 DEF_A_GAP_EXTEND
#define DEF_B_GAP_EXTEND_DAG	-3
#define DEF_XOVER_PENALTY	-140	/* CS only */

#define DEF_SW_VECT_THRESHOLD	60.0	/* == DEF_SW_FULL_THRESHOLD in lspace */
#define DEF_SW_VECT_THRESHOLD_DAG 50.0	/* smaller for Helicos DAG */
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
	struct re_score  *next;			/* linked list (final pass) */
	struct sw_full_results *sfrp;		/* alignment results (final pass) */
	u_int		  contig_num;		/* contig index (for filename)*/
	bool		  revcmpl;		/* from contig's reverse cmpl */
	int32_t		  score;		/* doubles as heap cnt in [0] */
	uint32_t	  index;		/* doubles as heap alloc in [0]*/
};

/*
 * Keep juicy bits of reads separate from the read structure itself. This significantly
 * reduces the cache footprint due to reads during scan() and hence speeds things up
 * enough to justify the ugliness.
 */
struct read_int {
	uint32_t	  offset;		/* offset in read array */
	char		 *name;
	uint32_t	 *read1;		/* the read as a bitstring */
	uint32_t	  read1_len;
	uint32_t	 *read2;		/* second of Helicos pair */
	uint32_t	  read2_len;
	int		  initbp;		/* colour space init letter */

	dag_cookie_t	  dag_cookie;		/* kmer graph cookie for glue */

	int		  swhits;		/* num of hits with sw */
	struct re_score  *scores;		/* top 'num_ouputs' scores */
	uint32_t	  final_matches;	/* num of final output matches*/
};

struct read_hit {
	uint32_t  g_idx;		/* kmer index in genome */
	uint32_t  r_idx;		/* kmer index in read */
};

struct read_elem {
	struct read_int  *ri;			/* pointer to the meat */

	/* the following are used during scan() */
	uint32_t	  last_swhit_idx;	/* index of last sw hit */
	uint16_t	  window_len;		/* per-read window length */
	uint8_t	  	  prev_hit;		/* prev index in 'hits' */
	uint8_t		  next_hit;		/* next index in 'hits' */
	struct read_hit	  hits[0];		/* size depends on num_matches */
};

/*
 * Each index of `readmap' points to an array of `readmap_entry' structures.
 * These structures in turn point to all reads that contain the indexing kmer
 * as well as track where in the read the kmer exists for colinearity checking.
 */
struct readmap_entry {
	uint32_t	offset;		/* offset to the read in array */
	uint32_t	r_idx;		/* kmer's index in the read sequence */
};
