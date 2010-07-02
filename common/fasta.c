/*	$Id: fasta.c,v 1.23 2009/06/16 23:26:21 rumble Exp $	*/

#include <assert.h>
#include <ctype.h>
#include <errno.h>
#include <fcntl.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <unistd.h>
#include <zlib.h>

#include <sys/types.h>
#include <sys/stat.h>
#include <sys/time.h>

#include "../common/fasta.h"
#include "../common/util.h"

static uint64_t total_ticks;

fasta_t
fasta_open(const char *file, shrimp_mode_t space, bool fastq)
{
	fasta_t fasta = NULL;
	struct stat sb;
	gzFile fp;
	uint64_t before = rdtsc();

	assert(space == MODE_COLOUR_SPACE || space == MODE_LETTER_SPACE);

	if (strcmp(file,"-")==0) {
		fp = gzdopen(fileno(stdin),"r");
	} else {
		if (stat(file, &sb)) {
			fprintf(stderr,"error: Syscall stat failed on file \"%s\"!\n",file);
			total_ticks += (rdtsc() - before);
			return NULL;
		}
		if (!S_ISREG(sb.st_mode)) {
			fprintf(stderr,"error: \"%s\" is not a regular file!\n",file);
			total_ticks += (rdtsc() - before);
			return NULL;
		}
		fp = gzopen(file, "r");
	}
	if (fp == NULL) {
		fprintf(stderr,"error: Failed to open \"%s\" for reading!\n",file);
		total_ticks += (rdtsc() - before);
		return NULL;		
	}

	fasta = (fasta_t)xmalloc(sizeof(*fasta));
	memset(fasta, 0, sizeof(*fasta));

	//if its fastq skip the header
	/*if (fastq) {
		z_off_t index=0;
		while (gzgets(fp, fasta->buffer, sizeof(fasta->buffer)) !=NULL) {
			if (fasta->buffer[0]!='#') {
				break;
			}
			index+=strlen(fasta->buffer);
		}
		//gzseek(fp,index,SEEK_SET);
		gzseek(fp,0,SEEK_SET);
	}*/

	fasta->fp = fp;
	fasta->file = xstrdup(file);
	fasta->space = space;
	fasta->fastq = fastq;
	fasta->parse_buffer_size = sizeof(fasta->buffer);
	fasta->parse_buffer = (char*)xmalloc(fasta->parse_buffer_size);
	fasta->header=true;
	memset(fasta->translate, -1, sizeof(fasta->translate));

	if (space == MODE_COLOUR_SPACE) {
		fasta->translate[(int)'0'] = BASE_0;
		fasta->translate[(int)'1'] = BASE_1;
		fasta->translate[(int)'2'] = BASE_2;
		fasta->translate[(int)'3'] = BASE_3;
		fasta->translate[(int)'4'] = BASE_N;
		fasta->translate[(int)'N'] = BASE_N;
		fasta->translate[(int)'n'] = BASE_N;
		fasta->translate[(int)'.'] = BASE_N;
		fasta->translate[(int)'X'] = BASE_X;
		fasta->translate[(int)'x'] = BASE_X;
	} else {	
		fasta->translate[(int)'A'] = BASE_A;
		fasta->translate[(int)'a'] = BASE_A;
		fasta->translate[(int)'C'] = BASE_C;
		fasta->translate[(int)'c'] = BASE_C;
		fasta->translate[(int)'G'] = BASE_G;
		fasta->translate[(int)'g'] = BASE_G;
		fasta->translate[(int)'T'] = BASE_T;
		fasta->translate[(int)'t'] = BASE_T;
		fasta->translate[(int)'U'] = BASE_U;
		fasta->translate[(int)'u'] = BASE_U;
		fasta->translate[(int)'M'] = BASE_M;
		fasta->translate[(int)'m'] = BASE_M;
		fasta->translate[(int)'R'] = BASE_R;
		fasta->translate[(int)'r'] = BASE_R;
		fasta->translate[(int)'W'] = BASE_W;
		fasta->translate[(int)'w'] = BASE_W;
		fasta->translate[(int)'S'] = BASE_S;
		fasta->translate[(int)'s'] = BASE_S;
		fasta->translate[(int)'Y'] = BASE_Y;
		fasta->translate[(int)'y'] = BASE_Y;
		fasta->translate[(int)'K'] = BASE_K;
		fasta->translate[(int)'k'] = BASE_K;
		fasta->translate[(int)'V'] = BASE_V;
		fasta->translate[(int)'v'] = BASE_V;
		fasta->translate[(int)'H'] = BASE_H;
		fasta->translate[(int)'h'] = BASE_H;
		fasta->translate[(int)'D'] = BASE_D;
		fasta->translate[(int)'d'] = BASE_D;
		fasta->translate[(int)'B'] = BASE_B;
		fasta->translate[(int)'b'] = BASE_B;
		fasta->translate[(int)'N'] = BASE_N;
		fasta->translate[(int)'n'] = BASE_N;
		fasta->translate[(int)'.'] = BASE_N;
		fasta->translate[(int)'X'] = BASE_X;
		fasta->translate[(int)'x'] = BASE_X;
	}
	return fasta;
}

void
fasta_close(fasta_t fasta)
{
	uint64_t before = rdtsc();

	gzclose(fasta->fp);
	free(fasta->file);
	free(fasta->parse_buffer);
	fasta->parse_buffer=0;
	free(fasta);

	total_ticks += (rdtsc() - before);
}

fasta_stats_t
fasta_stats()
{
	fasta_stats_t fs;

	fs = (fasta_stats_t)xmalloc(sizeof(*fs));
	fs->total_ticks = total_ticks;

	return (fs);
}

static char *
extract_name(char *buffer, char * * ranges)
{
	char *extracted;
	char *ret;
	int len;
	char * tok_save;

	assert(buffer[0] == '>' || buffer[0] == '@');

	/* Funny business for valgrind. See bottom of fasta_get_next. */
	/*
	extracted = strtrim(&buffer[1]);
	len = strlen(extracted);
	ret = (char *)xmalloc(len + 17);
	memcpy(ret, extracted, len);
	memset(ret + len, 0, 17);
	*/

	extracted = strtok_r(&buffer[1], "\t", &tok_save);
	extracted = strtrim(extracted);
	int full_length=strlen(extracted);
	for (len=0; len<full_length; len++) {
		if (extracted[len]==' ' || extracted[len]=='\t') {
			break;
		}
	}
	//len = strlen(extracted);
	ret = (char *)xmalloc(len + 17);
	memcpy(ret, extracted, len);
	memset(ret + len, 0, 17);

	if (ranges != NULL && (extracted = strtok_r(NULL, "\t", &tok_save)) != NULL) {
	  len = strlen(extracted);
	  *ranges = (char *)xmalloc(len + 17);
	  memcpy(*ranges, extracted, len);
	  memset(*ranges + len, 0, 17);
	}

	return (ret);
}


void fasta_write_fasta(FILE * file, char * seq) {
	assert(seq!=NULL);
	int index=0;
	int length=strlen(seq);
	char buffer[FASTA_PER_LINE];
	while (index<length) {
		int copy_over=MIN(FASTA_PER_LINE-1,length-index);
		memcpy(buffer,seq,copy_over);
		buffer[copy_over]='\0';
		fprintf(file,"%s\n",buffer);
		index+=copy_over;	
	}
}

void fasta_write_read(FILE * file, read_entry * re) {
		//TODO does not print out the ranges vars as well!
	if (re->qual==NULL) {
		fprintf(file,">%s\n",re->name);
		fasta_write_fasta(file,re->seq);
	} else {
		//its fastq
		fprintf(file,"@%s\n",re->name);
		fasta_write_fasta(file,re->seq);
		fprintf(file,"%s\n",re->plus_line);
		fasta_write_fasta(file,re->qual);
	}
};

bool
fasta_get_next_read_with_range(fasta_t fasta, read_entry * re )
{
	int i;
	uint32_t sequence_length = 0;
	uint32_t quality_length = 0;
	uint32_t plus_line_length = 0;
	uint32_t read_name_length = 0;
	uint64_t before = rdtsc();
	char c;
	assert(re!=NULL);
	if (fasta->fastq){
		c = '@';
	} else {
		c = '>';
	}
	re->name = re->seq = NULL;
	assert(fasta->parse_buffer!=NULL);
	assert(fasta->parse_buffer_size>0);

	/*
		The name of the read
	*/
	bool end_of_line=false;
	fasta->parse_buffer[0]='\0';
	while (fasta->leftover || fast_gzgets_safe(fasta) != NULL ) {
			//if was leftover, after this its gone
			fasta->leftover=false;
			while (fasta->parse_buffer_size <= (sizeof(fasta->buffer) + read_name_length)) {
					fasta->parse_buffer_size *= 2;
					fasta->parse_buffer = (char *)xrealloc(fasta->parse_buffer,fasta->parse_buffer_size);
			}
			for (i = 0; fasta->buffer[i] != '\0' && fasta->buffer[i]!='\n'; i++) {
				fasta->parse_buffer[read_name_length++] = fasta->buffer[i];
			}
			if (fasta->parse_buffer[0]!=c) {
				if (c=='>' && fasta->parse_buffer[0]=='@') {
					fprintf(stderr,"Expecting \">\" but got \"%c\" are you sure it's not FASTQ format?\n",fasta->parse_buffer[0]);
				} else if (c=='@' && fasta->parse_buffer[0]=='>') {
					fprintf(stderr,"Expecting \"@\" but got \"%c\" are you sure it's not FASTA format?\n",fasta->parse_buffer[0]);
				} else {
					fprintf(stderr,"Expecting \"%c\" but got \"%c\" are you sure it's right format?\n",c,fasta->parse_buffer[0]);
				}
				total_ticks += (rdtsc() - before);
				return (false);
			}
			if (fasta->buffer[i]=='\n') {
				end_of_line=true;
				fasta->parse_buffer[read_name_length]='\0';
				re->name = extract_name(fasta->parse_buffer, &(re->range_string));
				break;
			}

	}
	if (read_name_length==0) {
		total_ticks += (rdtsc() - before);
		return (false);
	}
	if (read_name_length<=1 || !end_of_line) {
		total_ticks += (rdtsc() - before);
		fprintf(stderr,"error: Invalid read name! Are you sure this is a FASTA or FASTQ file?\n%s\n",fasta->buffer);
		return (false);
	}	
	
	assert(re->name!=NULL);
	/*
		The read sequence
	*/
	fasta->parse_buffer[0]='\0';
	while (fast_gzgets_safe(fasta) != NULL ) {
			//Check if we have finished reading this sequence
			if (fasta->fastq && fasta->buffer[0]=='+') {
				break;
			} else if (!fasta->fastq && fasta->buffer[0]=='>') {
				fasta->leftover=true;
				break;
			}
			//otherwise keep reading the sequence
			while (fasta->parse_buffer_size <= (sizeof(fasta->buffer) + sequence_length)) {
					fasta->parse_buffer_size *= 2;
					fasta->parse_buffer = (char *)xrealloc(fasta->parse_buffer,fasta->parse_buffer_size);
			}
			for (i = 0; fasta->buffer[i] != '\0' && fasta->buffer[i]!='\n'; i++) {
				fasta->parse_buffer[sequence_length++] = fasta->buffer[i];
			}
	}
	if (sequence_length==0) {
		fprintf(stderr,"Read in sequence of length zero!\n");
		total_ticks += (rdtsc() - before);
		return (false);
	}
	assert(sequence_length>0);
	re->seq = (char *)xmalloc(sequence_length + 17);
	memcpy(re->seq, fasta->parse_buffer, sequence_length);
	memset(re->seq + sequence_length, 0, 17);
	if (fasta->fastq) {
		/*
			Read in the plus line
		*/
		//already got the first chunk above
		assert(fasta->buffer[0]=='+');
		do {
			//otherwise keep reading the sequence
			while (fasta->parse_buffer_size <= (sizeof(fasta->buffer) + plus_line_length)) {
					fasta->parse_buffer_size *= 2;
					fasta->parse_buffer = (char *)xrealloc(fasta->parse_buffer,fasta->parse_buffer_size);
			}
			for (i = 0; fasta->buffer[i] != '\0' && fasta->buffer[i]!='\n'; i++) {
				fasta->parse_buffer[plus_line_length++] = fasta->buffer[i];
			}
			if (fasta->buffer[i]=='\n') {
				fasta->parse_buffer[plus_line_length]='\0';
				break;
			}
		} while (fast_gzgets_safe(fasta) != NULL );
		if (plus_line_length<1 || fasta->buffer[strlen(fasta->buffer)-1]!='\n') {
			fprintf(stderr,"error: Error while readingin FASTQ entry!\n");
			total_ticks += (rdtsc() - before);
			return (false);	
		}
		re->plus_line = (char *)xmalloc(plus_line_length + 17);
		memcpy(re->plus_line, fasta->parse_buffer, plus_line_length);
		memset(re->plus_line + plus_line_length, 0, 17);
		
		/*
			Read in the qualities
		*/
		fasta->parse_buffer[0]='\0';
		while (fast_gzgets_safe(fasta) != NULL ) {
				while (fasta->parse_buffer_size <= (sizeof(fasta->buffer) + quality_length)) {
						fasta->parse_buffer_size *= 2;
						fasta->parse_buffer = (char *)xrealloc(fasta->parse_buffer,fasta->parse_buffer_size);
				}
				for (i = 0; fasta->buffer[i] != '\0' && fasta->buffer[i]!='\n'; i++) {
					fasta->parse_buffer[quality_length++] = fasta->buffer[i];
				}
				/*
					Want to make sure that we haven't written past buffer
					And want to make sure that we dont have more qual
					values then the length of the sequence
				*/
				assert(quality_length <= fasta->parse_buffer_size);
				assert(quality_length <= sequence_length);
			
				/*
					See if we are done reading quality values. This is a bit
					different for colour space and letter space. Since the first
					letter in the colour space read does not have a quality value	
				*/
				if (quality_length == sequence_length && shrimp_mode==MODE_LETTER_SPACE){
					break;
				} else if (quality_length == sequence_length-1 && shrimp_mode==MODE_COLOUR_SPACE) {
					break;
				} else {
					continue;
				}
		}
		if (quality_length != sequence_length && shrimp_mode==MODE_LETTER_SPACE){
			fprintf(stderr,"Read in quality string of wrong length!, %d vs %d\n",quality_length, sequence_length);
			free(re->seq);
			free(re->plus_line);
			total_ticks += (rdtsc() - before);
			return (false);
		} else if (quality_length != sequence_length-1 && shrimp_mode==MODE_COLOUR_SPACE) {
			fprintf(stderr,"Read in quality string of wrong length!, %d vs %d\n",quality_length, sequence_length);
			free(re->seq);
			free(re->plus_line);
			total_ticks += (rdtsc() - before);
			return (false);
		}
		re->qual = (char *)xmalloc(quality_length + 17);
		memcpy(re->qual, fasta->parse_buffer, quality_length);
		memset(re->qual + quality_length, 0, 17);
	}	


	/*
	 * Ensure nul-termination and allocate extra space to appease valgrind.
	 * I think strlen may do > char-sized loads and valgrind thinks it's
	 * conditionally jumping on uninitialised memory. That explanation
	 * doesn't seem right to me, but doing this makes it happy again.
	 */
	/**sequence = (char *)xmalloc(sequence_length + 17);
	memcpy(*sequence, readinseq, sequence_length);
	memset(*sequence + sequence_length, 0, 17);
	if (fasta->fastq){
		*qual = (char *)xmalloc(quality_length + 17);
		memcpy(*qual, readinqual, quality_length);
		memset(*qual + quality_length, 0, 17);
	}*/

	/* check if the sequence is rna (contains uracil and not thymine) */
	bool got_uracil = false;
	bool got_thymine = false;
	unsigned int j;

	for (j = 0; j < sequence_length; j++) {
		unsigned char chr = (re->seq)[j];
		got_thymine |= (chr == 'T' || chr == 't');
		got_uracil  |= (chr == 'U' || chr == 'u');
	}

	if (got_uracil && got_thymine)
		fprintf(stderr, "WARNING: sequence has both uracil and "
		    "thymine!?!\n");

	re->is_rna = (got_uracil && !got_thymine);

	assert(re->name != NULL);
	assert(re->seq != NULL);
	if (fasta->fastq) {
		assert(re->qual != NULL);
		assert(re->plus_line != NULL);
	}
	total_ticks += (rdtsc() - before);
	return (true);
}

int
fasta_get_initial_base(int space, char *sequence)
{

	assert(space == COLOUR_SPACE);

	switch (*sequence) {
	case 'A':
	case 'a':
		return (BASE_A);
	case 'C':
	case 'c':
		return (BASE_C);
	case 'G':
	case 'g':
		return (BASE_G);
	case 'T':
	case 't':
		return (BASE_T);
	}

	return (-1);
}

uint32_t *
fasta_bitfield_to_colourspace(fasta_t fasta, uint32_t *source, uint32_t length, bool is_rna)
{
	int a, lastbp = BASE_T;
	uint32_t *dst;
	uint32_t i;
	uint64_t before = rdtsc();

	assert(fasta->space == LETTER_SPACE);

	dst = (uint32_t *)xmalloc(BPTO32BW(length) * sizeof(uint32_t));
	memset(dst, 0, BPTO32BW(length) * sizeof(uint32_t));

	for (i = 0; i < length; i++) {
		a = EXTRACT(source, i);
		bitfield_insert(dst, i, lstocs(lastbp, a, is_rna));
		lastbp = a;
	}

	total_ticks += (rdtsc() - before);
	return (dst);
}

uint32_t *
fasta_sequence_to_bitfield(fasta_t fasta, char *sequence)
{
	uint32_t i, length, idx;
	uint32_t *bitfield;
	uint64_t before = rdtsc();
	int a;
	char c;

	length = strlen(sequence);
	bitfield = (uint32_t *)xmalloc(BPTO32BW(length) * sizeof(uint32_t));
	memset(bitfield, 0, BPTO32BW(length) * sizeof(uint32_t));

	i = 0;
	c = (char)sequence[0];
	if (fasta->space == COLOUR_SPACE) {
		if (c != 'A' && c != 'a' &&
		    c != 'C' && c != 'c' && 
		    c != 'G' && c != 'g' &&
		    c != 'T' && c != 't') {
			free(bitfield);
			total_ticks += (rdtsc() - before);
			return (NULL);
		}
		
		i = 1;
	}

	for (idx = 0; i < length; i++) {
		a = fasta->translate[(int)sequence[i]];
		if (a == -1) {
			fprintf(stderr, "error: invalid character ");
			if (isprint(a))
				fprintf(stderr, "(%c) ", a);
			else if (a != -1)
				fprintf(stderr, "(0x%x) ", a);
			fprintf(stderr, "in input file [%s]\n", fasta->file);
			fprintf(stderr, "       (Did you mix up letter "
			    "space and colour space programs?)\n");
			exit(1);
		}

		if (fasta->space == COLOUR_SPACE) {
			assert((a >= BASE_CS_MIN && a <= BASE_CS_MAX) ||
			    (a == BASE_N || a == BASE_X));
		} else {
			assert((a >= BASE_LS_MIN && a <= BASE_LS_MAX) ||
			    (a == BASE_N || a == BASE_X));
		}

		bitfield_append(bitfield, idx++, a);
	}

	if (fasta->space == COLOUR_SPACE)
		assert(idx == length - 1);
	else
		assert(idx == length);

	total_ticks += (rdtsc() - before);
	return (bitfield);
}

/*
 * Give BASE_x, return the appropriate character.
 *
 * NB: Since we're limited to 4-bits, BASE_X returns 'N'.
 */
char
base_translate(int base, bool use_colours)
{
	/*
	 * NB: colour-space only valid for 0-3 and BASE_N/BASE_X
	 *     BASE_N is reported as a skipped cycle: '.' in CS.
	 */
	char cstrans[] = { '0', '1', '2', '3', '!', '@', '#', '$',
			   '%', '^', '&', '*', '?', '~', ';', '.' };
	char lstrans[] = { 'A', 'C', 'G', 'T', 'U', 'M', 'R', 'W',
			   'S', 'Y', 'K', 'V', 'H', 'D', 'B', 'N' };

	if (use_colours) {
		assert((base >= BASE_CS_MIN && base <= BASE_CS_MAX) ||
		    (base == BASE_N || base == BASE_X));
		return (cstrans[base]);
	} else {
		assert((base >= BASE_LS_MIN && base <= BASE_LS_MAX) ||
		    (base == BASE_N || base == BASE_X));
		return (lstrans[base]);
	}
}
