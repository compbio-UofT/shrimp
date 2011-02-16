/*	$Id: shrimp2sam.c,v 1.23 2009/06/16 23:26:25 dlister Exp $	*/

#include <assert.h>
#include <dirent.h>
#include <errno.h>
#include <inttypes.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <zlib.h>

#include <sys/types.h>
#include <sys/stat.h>

#include "../common/dag_glue.h"
#include "../common/fasta.h"
#include "../common/dynhash.h"
#include "../common/input.h"
#include "../common/output.h"
#include "../common/util.h"
#include "../common/version.h"

#include "../gmapper/gmapper-definitions.h"

static shrimp_mode_t shrimp_mode;

static char ** ref_contig_names = NULL;
static uint32_t * ref_contig_lens = NULL;
static uint32_t n_ref_contigs = 0;

static dynhash_t read_list;		/* cache of reads we need */

/* probcalc/rmapper output cache */
struct fpo {
	struct input input;		/* input fields */
	struct fpo  *next_ordered;	/* global list in input order */
	struct fpo  *next_contig;	/* per-contig/revcmpl list */
	char        *output_normal;	/* alignment output */
	char        *output_pretty;	/* alignment output */
};

struct contig_ll {
	struct fpo  *head;		/* per-contig list */
	struct fpo  *head_revcmpl;
};

/* for read and contig lists */
struct sequence {
	char	 *name;
	char	 *sequence;
	uint32_t  sequence_len;
	int	  initbp;		/* initial base pair (colourspace) */
	bool	  revcmpl;
	bool	  is_rna;
};

#define MAX_READ_LEN	5000		/* ridiculously high upper bound */

static void
load_output_file(char *file)
{
	struct input inp;
	gzFile inf;
	struct sequence *val;
	bool freeread = false;

	inf = gzopen(file, "r");
	if (inf == NULL) {
		fprintf(stderr, "error: failed to open probcalc/rmapper output "
		    "file [%s]: %s\n", file, strerror(errno));
		exit(1);
	}

	// print the header:

	fprintf(stdout,"@HD\tVN:%i\tSO:%s\n",1,"unsorted");

	uint s;
	for(s = 0; s < n_ref_contigs; s++){
		fprintf(stdout,"@SQ\tSN:%s\tLN:%u\n",ref_contig_names[s],ref_contig_lens[s]);
	}
	fprintf(stdout,"@PG\tID:%s\tVN:%s\n","shrimp2sam",SHRIMP_VERSION_STRING);

	while (input_parseline(inf, &inp)) {
		char * read;
		if(dynhash_find(read_list, inp.read, NULL,(void **)&val)){
			if(inp.flags & INPUT_FLAG_IS_REVCMPL){
				read = (char *)xmalloc(sizeof(char)*(strlen(val->sequence)+1));
				freeread = true;
				reverse_complement_read_ls_text(val->sequence,read);
			} else {
				freeread = false;
				read = val->sequence;
			}
		} else {
			freeread = false;
			read = (char *)"*";
		}
		char * cigar;
		cigar = (char *)xmalloc(sizeof(char)*200);
		edit2cigar(inp.edit,inp.read_start,inp.read_end,inp.read_length,cigar);

		if(inp.flags & INPUT_FLAG_IS_REVCMPL){
		    	char * cigar_reverse = (char *)xmalloc(sizeof(char)*strlen(cigar)+1);
		    	*cigar_reverse = '\0';
		    	char * tmp = (char *)xmalloc(sizeof(char)*strlen(cigar)+1);
		    	*tmp = '\0';
		    	//char * tmp2 = (char *)xmalloc(sizeof(char)*strlen(cigar)+1);
		    	char * ptr = cigar;
		    	char * last = tmp;
		    	for (ptr = cigar; *ptr != '\0'; ptr++){
		    		*last = *ptr;
		    		last ++;
		    		*last = '\0';
		    		if (!(*ptr <= '9' && *ptr >= '0')){
		    			strcat(tmp,cigar_reverse);
		    			strcpy(cigar_reverse,tmp);
		    			last = tmp;
		    			*tmp = '\0';
		    		}
		    	}
		    	free(cigar);
		    	cigar = cigar_reverse;
		    	free(tmp);

		    }

		fprintf(stdout,"%s\t%i\t%s\t%u\t%i\t%s\t%s\t%u\t%i\t%s\t%s\tAS:i:%i\n",
				inp.read,
				((inp.flags & INPUT_FLAG_IS_REVCMPL) ? 16 :0),
				inp.genome,
				inp.genome_start + 1,
				255,
				cigar,
				"*",
				inp.read_start,
				0,
				read,
				"*",
				inp.score);

		free(cigar);
		if(freeread){
			free(read);
		}
	}

	gzclose(inf);
}

static void
load_ref_contigs(char **files,int nfiles){
	int c;
	fasta_t fasta;
	char *seq,*name;


	for (c = 0; c < nfiles; c++){
		fasta = fasta_open(files[c],MODE_LETTER_SPACE,false);
		while(1){
			if(!fasta_get_next_contig(fasta, &name, &seq, NULL)){
				break;
			}

			int n = n_ref_contigs;
			n_ref_contigs++;

			ref_contig_names = (char **)xrealloc(ref_contig_names,sizeof(char *)*n_ref_contigs);
			ref_contig_lens = (uint32_t *)xrealloc(ref_contig_lens,sizeof(uint32_t)*n_ref_contigs);

			ref_contig_names[n] = name;
			ref_contig_lens[n] = strlen(seq);
			free(seq);

		}
		fasta_close(fasta);
	}
}

static void
load_reads_file(char *fpath)
{
	fasta_t fasta;
	struct sequence *s;
	char *name, *seq;

	fasta = fasta_open(fpath, MODE_LETTER_SPACE, false);
	if (fasta == NULL) {
		fprintf(stderr, "error: failed to parse reads fasta file "
		    "[%s]\n", fpath);
		exit(1);
	}

	while (fasta_get_next_contig(fasta, &name, &seq, NULL)) {

		s = (struct sequence *)xmalloc(sizeof(*s));
		s->sequence = seq;
		s->sequence_len = strlen(seq);
		s->name = name;
		s->revcmpl = false;

		if (!dynhash_add(read_list, s->name, s)) {
			fprintf(stderr, "error: failed to add read to list - "
				"probably out of memory\n");
			exit(1);
		}
	}

	fasta_close(fasta);
}

static uint32_t
keyhasher(void *k)
{

	return (hash_string((char *)k));
}

static int
keycomparer(void *a, void *b)
{

	return (strcmp((char *)a, (char *)b));
}

static void
usage(char *progname)
{
	char *slash;

	slash = strrchr(progname, '/');
	if (slash != NULL)
		progname = slash + 1;

	fprintf(stderr, "usage: %s shrimp_output_file reads_file genome_file1 [genome_file2 ...]\n", progname);


	exit(1);
}

int
main(int argc, char **argv)
{
	char *progname;
	char const * optstr;
	int ch;

	fprintf(stderr, "--------------------------------------------------"
	    "------------------------------\n");
	fprintf(stderr, "shrimp2sam: %s SPACE.\nSHRiMP %s [%s]\n",
	    get_mode_string(shrimp_mode), SHRIMP_VERSION_STRING, get_compiler());
	fprintf(stderr, "--------------------------------------------------"
	    "------------------------------\n");

	progname = argv[0];

	if (shrimp_mode == MODE_COLOUR_SPACE){
		fprintf(stderr,"Colour space not currently supported");
		exit(1);
	}
	else {
		optstr = "";
	}

	while ((ch = getopt(argc, argv, optstr)) != -1) {
		switch (ch) {
		default:
			usage(progname);
		}
	}
	argc -= optind;
	argv += optind;

	if (argc < 3)
		usage(progname);


	read_list   = dynhash_create(keyhasher, keycomparer);

	fprintf(stderr,"Loading contigs\n");
	load_ref_contigs(argv + 2, argc - 2);
	fprintf(stderr,"Loading reads\n");
	load_reads_file(argv[1]);

	fprintf(stderr, "Loading shrimp output file...\n");
	load_output_file(argv[0]);
	return (0);
}
