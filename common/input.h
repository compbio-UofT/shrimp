/*	$Id$	*/

struct input {
	char *read;				/* read name */
	char *genome;				/* genome/contig name */
	char *read_seq;				/* read sequence */
	bool  revcmpl;				/* strand */
	int   score;				/* alignment score */
	u_int genome_start;			/* start of alignment - genome*/
	u_int genome_end;			/* end of alignment - genome */
	u_int read_start;			/* start of alignment - read */
	u_int read_end;				/* end of alignment - read */ 
	u_int read_length;			/* length of the read */
	u_int matches;				/* number of matches */
	u_int mismatches;			/* number of mismatches */
	u_int insertions;			/* number of insertions */
	u_int deletions;			/* number of deletions */
	u_int crossovers;			/* number of crossovers (CS) */
};

bool	input_parseline(FILE *, struct input *);
