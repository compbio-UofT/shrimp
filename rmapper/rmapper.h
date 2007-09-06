/*	$Id$	*/

#define RMAPPER_VERSION_STR	"0.42"

extern const int use_colours;

#define EXTRACT(_genome, _i) (((_genome)[(_i) / 8] >> (4 * ((_i) % 8))) & 0xf)
#define BPTO32BW(_x) (((_x) + 7) / 8)

#define DEF_SPACED_SEED		"11110111"
#define DEF_WINDOW_LEN		30
#define DEF_NUM_MATCHES		2
#define DEF_TABOO_LEN		4
#define DEF_NUM_OUTPUTS		100
#define DEF_MAX_READ_LEN	64
#define DEF_KMER_STDDEV_LIMIT	-1

#define DEF_MATCH_VALUE		100
#define DEF_MISMATCH_VALUE	-70
#define DEF_GAP_OPEN		-100
#define DEF_GAP_EXTEND		-70
#define DEF_XOVER_PENALTY	-200
#define DEF_SW_THRESHOLD	1875

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
	int32_t  score;				/* doubles as heap cnt in [0] */
	uint32_t index;
};

struct read_elem {
	char		 *name;
	struct read_elem *next;			/* next in read list */
	uint32_t	 *read;			/* the read as a bitstring */
	uint32_t	  read_len;

	int		  swhits;		/* num of hits with sw */
	uint32_t	  last_swhit_idx;	/* index of last sw hit */
	struct re_score  *scores;		/* top 'num_ouputs' scores */

	int		  initbp;		/* colour space init letter */

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
extern int   optind;

/*
 * Given the first letter of a pair corresponding to a colour, and the colour
 * itself, obtain the next letter.
 */
static inline int
cstols(int first_letter, int colour)
{
	assert(first_letter >= 0 && first_letter <= 3);
	assert(colour >= 0 && colour <= 3);

	if ((first_letter % 2) == 0)
		return ((4 + first_letter + colour) % 4);
	else
		return ((4 + first_letter - colour) % 4);
}

static inline int
lstocs(int first_letter, int second_letter)
{
	const int colourmat[5][5] = {
		{ 0, 1, 2, 3, 4 },
		{ 1, 0, 3, 2, 4 },
		{ 2, 3, 0, 1, 4 },
		{ 3, 2, 1, 0, 4 },
		{ 4, 4, 4, 4, 4 }
	};

	assert(first_letter  >= 0 && first_letter  <= 4);
	assert(second_letter >= 0 && second_letter <= 4);

	return (colourmat[first_letter][second_letter]);
}
