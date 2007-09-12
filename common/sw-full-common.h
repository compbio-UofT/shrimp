/*	$Id$	*/

struct sw_full_results {
	/* Common fields */
	int read_start;					/* read index of map */
	int rmapped;					/* read mapped length */
	int genome_start;				/* genome index of map*/
	int gmapped;					/* genome mapped len */
	int matches;					/* # of matches */
	int mismatches;					/* # of substitutions */
	int insertions;					/* # of insertions */
	int deletions;					/* # of deletions */
	int score;					/* final SW score */

	/* Colour space fields */
	int crossovers;					/* # of mat. xovers */
};
