/*	$Id$	*/

struct sw_full_results {
	int read_start;					/* read index of map */
	int genome_start;				/* genome index of map*/
	int mapped;					/* read mapped length */
	int matches;					/* # of matches */
	int mismatches;					/* # of substitutions */
	int insertions;					/* # of insertions */
	int deletions;					/* # of deletions */
	int score;					/* final SW score */
};

int	sw_full_setup(int, int, int, int, int, int);
void	sw_full_stats(uint64_t *, uint64_t *, uint64_t *, double *);
void	sw_full(uint32_t *, int, int, uint32_t *, int, int, char **, char **,
	    struct sw_full_results *);
