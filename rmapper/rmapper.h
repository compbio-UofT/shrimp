/*	$Id$	*/

#define MAX_KMERS_PER_READ	512

#define EXTRACT(_genome, _i) (((_genome)[(_i) >> 4] >> (2 * ((_i) & 15))) & 0x3)

#define DEF_SPACED_SEED		"11110111"
#define DEF_WINDOW_LEN		30
#define DEF_NUM_MATCHES		2
#define DEF_TABOO_LEN		4
#define DEF_NUM_OUTPUTS		100
#define DEF_MATCH_VALUE		100
#define DEF_MISMATCH_VALUE	-70
#define DEF_GAP_OPEN		-100
#define DEF_GAP_EXTEND		-70
#define DEF_SW_THRESHOLD	1875

struct re_score {
	int32_t  score;				/* doubles as heap cnt in [0] */
	uint32_t index;
};

struct read_elem {
	char		 *name;
	struct read_elem *next;			/* next in read list */
	uint32_t	  read[2];		/* XXX - max 32 bases */
	uint32_t	  read_len;

	int		  swhits;		/* num of hits with sw */
	uint32_t	  last_swhit_idx;	/* index of last sw hit */
	struct re_score  *scores;		/* top 'num_ouputs' scores */

	uint32_t	  prev_hit;		/* prev index in 'hits' */
	uint32_t	  next_hit;		/* next index in 'hits' */
	uint32_t	  hits[0];		/* size depends on num_matches*/
};

/* for kmer to read map */
struct read_node {
	struct read_elem *read;
	struct read_node *next;
};

/* for optarg (and to shut up icc) */
extern char *optarg;
extern int optind;
