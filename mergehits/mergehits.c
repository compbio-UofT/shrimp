/*	$Id$	*/

#include <assert.h>
#include <dirent.h>
#include <errno.h>
#include <inttypes.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
 
#include <sys/types.h>
#include <sys/stat.h>

#include "../common/fasta.h"
#include "../common/dynhash.h"
#include "../common/sw-full-common.h"
#include "../common/sw-full-cs.h"
#include "../common/sw-full-ls.h"
#include "../common/input.h"
#include "../common/output.h"
#include "../common/util.h"
#include "../common/version.h"

/* alignment list member */
struct alist_item {
	struct input       input;	/* input fields */
	struct alist_item *next;	/* next alignment in list */
};

/* cached data associated with a read */
struct sequence {
	char	 *name;
	uint32_t *sequence;
	uint32_t  sequence_len;
	int	  initbp;		      /* initial base pair (cs only) */
	struct    alist_item *alignments;     /* start of alignment list */
	struct    alist_item *last_alignment; /* end of alignment list */
	struct    sequence *next;	      /* next sequence in list */
};

static struct sequence *read_list;	/* list of read names and alignments */
static dynhash_t read_index;		/* name index for read_list */
static uint64_t nalignments;		/* number of alignments */

#ifdef USE_COLOURS
const bool use_colours = true;
#else
const bool use_colours = false;
#endif

static bool Rflag = false;		/* don't output read sequence */

#define MAX_READ_LEN	5000		/* ridiculously high upper bound */

static uint64_t
iter_reg_files(char *dir, void (*fh)(char *, void *), void *arg)
{
	char fpath[2048];
	struct stat sb;
	DIR *dp;
	struct dirent *de;
	uint64_t files;

	dp = opendir(dir);
	if (dp == NULL) {
		fprintf(stderr, "error: failed to open directory [%s]: %s\n",
		    dir, strerror(errno));
		exit(1);
	}

	files = 0;
	while (1) {
		de = readdir(dp);
		if (de == NULL)
			break;

		if (de->d_type != DT_REG && de->d_type != DT_LNK)
			continue;

		fpath[0] = '\0';
		strcpy(fpath, dir);
		strcat(fpath, "/");
		strcat(fpath, de->d_name);

		/* ensure it's a regular file or link to one */
		if (stat(fpath, &sb) == -1) {
			fprintf(stderr, "error: failed to stat file [%s]: %s\n",
			    fpath, strerror(errno));
			exit(1);
		}

		if (S_ISREG(sb.st_mode)) {
			fh(fpath, arg);
			files++;
		} else {
			fprintf(stderr, "warning: [%s] is neither a regular "
			    "file, nor a link to one; skipping...", fpath);
			continue;
		}
	}

	closedir(dp);

	return (files);
}

static void
load_output_file(char *file)
{
	struct input inp;
	FILE *fp;
	struct alist_item *alignment;
	struct sequence *seq, *last_read;

	fp = fopen(file, "r");
	if (fp == NULL) {
		fprintf(stderr, "error: failed to open probcalc/rmapper "
			"output file [%s]: %s\n", file, strerror(errno));
		exit(1);
	}

	assert(read_list == NULL);
	last_read = NULL;
	while (input_parseline(fp, &inp)) {
		alignment = xmalloc(sizeof(*alignment));
		memcpy(&alignment->input, &inp, sizeof(alignment->input));
		alignment->next = NULL;

		if (dynhash_find(read_index, inp.read, NULL, (void *)&seq)) {
			/*
			 * we've already seen an alignment for the same
			 * read, so add this one to that read's list
			 */
			assert(seq != NULL);
			assert(seq->last_alignment != NULL);
			assert(seq->last_alignment->next == NULL);
			seq->last_alignment->next = alignment;
			seq->last_alignment = alignment;
		} else {
			/*
		         * this alignment is for a read we haven't seen
			 * yet, so allocate a new list
			 */
			seq = xmalloc(sizeof(*seq));
			seq->name = xstrdup(inp.read);
			seq->sequence = NULL;
			seq->sequence_len = 0;
			seq->initbp = 0;
			seq->alignments = alignment;
			seq->last_alignment = alignment;
			seq->next = NULL;
			if (last_read == NULL) {
				read_list = seq;
			} else {
				last_read->next = seq;
			}
			last_read = seq;

			if (dynhash_add(read_index, seq->name, seq) == false) {
				fprintf(stderr, "error: failed to add "
					"read to index - "
					"probably out of memory\n");
				exit(1);
			}
			nalignments++;
		}
	}
}

static void
store_sequence(char *name, uint32_t *read, uint32_t read_len, int initbp)
{
	struct sequence *seq;

	if (!dynhash_find(read_index, name, NULL, (void *)&seq)) {
		/* this read has no alignments so skip it */
		return;
	}
	assert(seq != NULL);
	seq->sequence = xmalloc(sizeof(seq->sequence[0]) * BPTO32BW(read_len));
	memcpy(seq->sequence, read,
	    sizeof(seq->sequence[0]) * BPTO32BW(read_len));
	seq->sequence_len = read_len;
	seq->initbp = initbp;
}

static void
load_reads_file_helper(int base, ssize_t offset, int isnewentry, char *name,
    int initbp)
{
	static char     *seq_name;
	static uint32_t *read;
	static uint32_t  read_len;
	static int       read_initbp;

	/* shut up, icc */
	(void)offset;

	/* handle initial call to alloc resources */
	if (base == FASTA_ALLOC) {
		assert(seq_name == NULL);

		if (read == NULL) {
			read = xmalloc(sizeof(read[0]) *
			    BPTO32BW(MAX_READ_LEN));
			memset(read, 0, sizeof(read[0]) *
			    BPTO32BW(MAX_READ_LEN));
		}
		read_len = 0;
		read_initbp = 0;

		return;
	} else if (base == FASTA_DEALLOC) {
		assert(read != NULL);
		if (seq_name != NULL) {
			store_sequence(seq_name, read, read_len, read_initbp);
			free(seq_name);
			seq_name = NULL;
		}
		free(read);
		read = NULL;
		read_len = 0;
		read_initbp = 0;
		return;
	}

	if (isnewentry) {
		if (seq_name != NULL) {
			store_sequence(seq_name, read, read_len, read_initbp);
			free(seq_name);
			seq_name = NULL;
		}
		seq_name = xstrdup(name);
		read_initbp = initbp;
		memset(read, 0, sizeof(read[0]) * BPTO32BW(MAX_READ_LEN));
		read_len = 0;
	}

	assert(seq_name != NULL);
	assert(read_len < MAX_READ_LEN);
	assert(base >= 0 && base <= 15);

	bitfield_append(read, read_len++, base);
}

static void
load_reads_file(char *fpath, void *arg)
{
	ssize_t ret;

	/* shut up, icc */
	(void)arg;

	if (use_colours)
		ret = load_fasta(fpath, load_reads_file_helper, COLOUR_SPACE);
	else
		ret = load_fasta(fpath, load_reads_file_helper, LETTER_SPACE);

	if (ret == -1) {
		fprintf(stderr, "error: failed to reads parse fasta file "
		    "[%s]\n", fpath);
		exit(1);
	}
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
print_alignments()
{
	struct sequence   *seq;
	struct alist_item *alignment;
	char              *read_str;
	int32_t		   genome_align_start;
	uint32_t	   genome_align_size;
	uint32_t	   read_align_size;

	for (seq = read_list; seq != NULL; seq = seq->next) {
		printf(">%s", seq->name);
		for (alignment = seq->alignments; alignment != NULL;
		     alignment = alignment->next) {
			/*
			 * NB: internally 0 is first position, output uses 1.
			 *     adjust.
			 */
			genome_align_start =
			    INPUT_IS_REVCMPL(&alignment->input) ? 
			    -(alignment->input.genome_start + 1) :
			alignment->input.genome_start + 1;
			genome_align_size = alignment->input.genome_end -
			    alignment->input.genome_start + 1;
			read_align_size = alignment->input.read_end -
			    alignment->input.read_start + 1;
			printf(",%s.%d.%u.%hu.%u.%d", alignment->input.genome,
			    genome_align_start, genome_align_size,
			    alignment->input.read_start + 1, read_align_size,
			    alignment->input.score);
		}
		printf("\n");
		if (seq->sequence_len > 0) {
			read_str = readtostr(seq->sequence, seq->sequence_len,
			    use_colours, seq->initbp);
			puts(read_str);
		}
	}
}

static void
usage(char *progname)
{
	char *slash;

	slash = strrchr(progname, '/');
	if (slash != NULL)
		progname = slash + 1;

	fprintf(stderr, "usage: %s [-R reads_dir] shrimp_results_file\n",
	    progname);

	exit(1);
}

int
main(int argc, char **argv)
{
	char *fpout, *readsdir, *progname;
	int ch;

	/* shut up, gcc */
	readsdir = NULL;

	fprintf(stderr, "--------------------------------------------------"
	    "------------------------------\n");
	fprintf(stderr, "mergehits: %s SPACE. (SHRiMP version %s)\n",
	    (use_colours) ? "COLOUR" : "LETTER", SHRIMP_VERSION_STRING);
	fprintf(stderr, "--------------------------------------------------"
	    "------------------------------\n");

	progname = argv[0];

	while ((ch = getopt(argc, argv, "R:")) != -1) {
		switch (ch) {
		case 'R':
			Rflag = true;
			readsdir = optarg;
			break;
		default:
			usage(progname);
		}
	}
	argc -= optind;
	argv += optind;

	if (argc != 1)
		usage(progname);

	read_index   = dynhash_create(keyhasher, keycomparer);
	if (read_index == NULL) {
		fprintf(stderr, "error: failed to allocate read index\n");
		exit(1);
	}

	fpout = argv[0];

	fprintf(stderr, "Loading output file...\n");
	load_output_file(fpout);
	fprintf(stderr, "Loaded %" PRIu64 " alignments from rmapper/probcalc "
	    "output\n", nalignments);

	if (Rflag) {
		fprintf(stderr, "Loading read file(s)...\n");
		iter_reg_files(readsdir, load_reads_file, NULL);
	}

	print_alignments();

	return (0);
}
