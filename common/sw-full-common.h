/*	$Id$	*/

struct sw_full_results {
	/* Common fields */
	int read_start;					/* read index of map */
	int genome_start;				/* genome index of map*/
	int mapped;					/* read mapped length */
	int matches;					/* # of matches */
	int mismatches;					/* # of substitutions */
	int insertions;					/* # of insertions */
	int deletions;					/* # of deletions */
	int score;					/* final SW score */

	/* Colour space fields */
	int crossovers;					/* # of mat. xovers */
};
