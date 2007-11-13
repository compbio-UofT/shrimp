/*	$Id$	*/

struct input {
	char    *read;				/* read name */
	char    *genome;			/* genome/contig name */
	char    *read_seq;			/* read sequence */
	uint8_t  flags;				/* various flags; see below */
	int32_t  score;				/* alignment score */
	uint32_t genome_start;			/* start of alignment - genome*/
	uint32_t genome_end;			/* end of alignment - genome */
	uint16_t read_start;			/* start of alignment - read */
	uint16_t read_end;			/* end of alignment - read */ 
	uint16_t read_length;			/* length of the read */
	uint16_t matches;			/* number of matches */
	uint16_t mismatches;			/* number of mismatches */
	uint16_t insertions;			/* number of insertions */
	uint16_t deletions;			/* number of deletions */
	uint16_t crossovers;			/* number of crossovers (CS) */
	double   pchance;			/* pchance from probcalc */
	double   pgenome;			/* pgenome from probcalc */
	double   normodds;			/* normodds from probcalc */
};

#define INPUT_FLAG_IS_REVCMPL	0x00000001
#define INPUT_FLAG_HAS_PCHANCE	0x00000002
#define INPUT_FLAG_HAS_PGENOME	0x00000004
#define INPUT_FLAG_HAS_NORMODDS	0x00000008

#define INPUT_IS_REVCMPL(_i)	((_i)->flags & INPUT_FLAG_IS_REVCMPL)
#define INPUT_HAS_PCHANCE(_i)	((_i)->flags & INPUT_FLAG_HAS_PCHANCE)
#define INPUT_HAS_PGENOME(_i)	((_i)->flags & INPUT_FLAG_HAS_PGENOME)
#define INPUT_HAS_NORMODDS(_i)	((_i)->flags & INPUT_FLAG_HAS_NORMODDS)

bool	input_parseline(FILE *, struct input *);
