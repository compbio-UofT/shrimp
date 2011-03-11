#define _MODULE_GMAPPER

#ifndef CXXFLAGS
#define CXXFLAGS "?"
#endif

#include <assert.h>
#include <ctype.h>
#include <errno.h>
#include <fcntl.h>
#include <inttypes.h>
#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <unistd.h>
#include <zlib.h>
#include <omp.h>	// OMP multithreading
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/time.h>
#include <getopt.h>

#include "../gmapper/gmapper.h"
#include "../gmapper/seeds.h"
#include "../gmapper/genome.h"
#include "../gmapper/mapping.h"

#include "../common/hash.h"
#include "../common/fasta.h"
#include "../common/util.h"
#include "../common/version.h"
#include "../common/sw-full-common.h"
#include "../common/sw-full-cs.h"
#include "../common/sw-full-ls.h"
#include "../common/sw-vector.h"
#include "../common/output.h"
#include "../common/input.h"
#include "../common/read_hit_heap.h"

/* heaps */
DEF_HEAP(uint32_t, char*, out)



/*
	Get hit stats
*/
int
hit_edit_distance(struct read_hit * rh) {
  //find how many perfect matches, off by 1 matches and off by 2 matches
  //       int matches;                            /* # of matches */
  //        int mismatches;                         /* # of substitutions */
  //      int insertions;                         /* # of insertions */
  //      int deletions;                          /* # of deletions */
	int edit_distance=0;
	edit_distance+=rh->sfrp->mismatches;
	edit_distance+=rh->sfrp->insertions;
	edit_distance+=rh->sfrp->deletions;
	assert(edit_distance>=0);	
	return edit_distance;	
}

//TODO move this to utils
void reverse(char* s, char* t) {
       int l=strlen(s);
       int i;
       for (i=0; i<l; i++) {
               switch (s[i]) {
                       case 'A': t[l-i-1]='T'; break;
                       case 'a': t[l-i-1]='t'; break;
                       case 'T': t[l-i-1]='A'; break;
                       case 't': t[l-i-1]='a'; break;

                       case 'C': t[l-i-1]='G'; break;
                       case 'c': t[l-i-1]='g'; break;
                       case 'G': t[l-i-1]='C'; break;
                       case 'g': t[l-i-1]='c'; break;

                       case '-': t[l-i-1]='-'; break;

                       case 'N': t[l-i-1]='N'; break;
                       case 'n': t[l-i-1]='n'; break;
		       case '.': t[l-i-1]='.'; break;

                       case 'R': t[l-i-1]='Y'; break;
                       case 'r': t[l-i-1]='y'; break;
                       case 'Y': t[l-i-1]='R'; break;
                       case 'y': t[l-i-1]='r'; break;

                       case 'S': t[l-i-1]='S'; break;
                       case 's': t[l-i-1]='s'; break;
                       case 'W': t[l-i-1]='W'; break;
                       case 'w': t[l-i-1]='w'; break;

                       case 'K': t[l-i-1]='M'; break;
                       case 'k': t[l-i-1]='m'; break;
                       case 'M': t[l-i-1]='K'; break;
                       case 'm': t[l-i-1]='k'; break;

                       case 'B': t[l-i-1]='V'; break;
                       case 'b': t[l-i-1]='v'; break;
                       case 'V': t[l-i-1]='B'; break;
                       case 'v': t[l-i-1]='b'; break;

                       case 'D': t[l-i-1]='H'; break;
                       case 'd': t[l-i-1]='h'; break;
                       case 'H': t[l-i-1]='D'; break;
                       case 'h': t[l-i-1]='d'; break;
	
                       default:
                               fprintf(stderr,"There has been a error in getting reverse complement of %s\n",s);
                               exit(1);
               }
       }
       t[l]='\0';
       //printf("%s vs %s\n", s , t);
       strcpy(s,t);
}

/*
	read_start and read_end are 1 based
*/
cigar_t * make_cigar(int read_start, int read_end , int read_length, char* qralign,char* dbalign) {
	cigar_t * cigar = (cigar_t*)xmalloc(sizeof(cigar_t));
	cigar->size=2; 
	assert(cigar->size>=2); //for start and end without loop
	cigar->ops=(char*)xmalloc(sizeof(char)*cigar->size);
	cigar->lengths=(uint32_t*)xmalloc(sizeof(uint32_t)*cigar->size);
	int used=0;
	if (read_start>1) {
		assert(cigar->size-used>0);
		cigar->ops[used]='S';
		cigar->lengths[used]=read_start-1;
		used++;
	}
	int i=0; int qralign_length=strlen(qralign);
	while (i<qralign_length) {
		int length; char op;
		if (qralign[i]=='-') {
			for (length=0; qralign[i+length]=='-' && i+length<qralign_length; length++);
			op='D';
		} else if (dbalign[i]=='-') {
			for (length=0; dbalign[i+length]=='-' && i+length<qralign_length; length++);
			op='I';
		} else {
			for (length=0; dbalign[i+length]!='-' && qralign[i+length]!='-' && i+length<qralign_length; length++);
			op='M';
		}
		while ((used+1)>=cigar->size) { //make it bigger, want to make sure have enough for one more after loop!
			assert(cigar->size!=0);
			cigar->size*=2;
			cigar->ops=(char*)xrealloc(cigar->ops,sizeof(char)*cigar->size);
			cigar->lengths=(uint32_t*)xrealloc(cigar->lengths,sizeof(uint32_t)*cigar->size);			
		}
		cigar->ops[used]=op;
		cigar->lengths[used]=length;
		i+=length;
		used++;		
	}
	if (read_end!=read_length) {
		assert(used<cigar->size); //by loop invariant and initial size >=2
		cigar->ops[used]='S';
		cigar->lengths[used]=read_length-read_end; 
		used++;
	}
	cigar->ops=(char*)xrealloc(cigar->ops,sizeof(char)*used);
	cigar->lengths=(uint32_t*)xrealloc(cigar->lengths,sizeof(uint32_t)*used);
	cigar->size=used;
	return cigar;
} 


void reverse_cigar(cigar_t * cigar) {
	char ops[cigar->size];
	uint32_t lengths[cigar->size];
	memcpy(ops,cigar->ops,sizeof(char)*cigar->size);
	memcpy(lengths,cigar->lengths,sizeof(uint32_t)*cigar->size);
	int i;	
	for (i=0; i<cigar->size; i++) {
		cigar->ops[i]=ops[cigar->size-i-1];
		cigar->lengths[i]=lengths[cigar->size-i-1];
	}
	return;
}

char* make_cigar_string(cigar_t * cigar) {
	int string_length=cigar->size; //1 char for each op
	int i;
	for (i=0; i<cigar->size; i++) {
		int j=0,length;
		for (length=cigar->lengths[i]; length>0; j++)
			length/=10;
		string_length+=j;
	}
	string_length++; // for null term
	char * ret = (char*)xmalloc(sizeof(char)*(string_length));
	int used=0;
	for (i=0; i<cigar->size; i++) {
		//printf("%d%c\n",cigar->lengths[i],cigar->ops[i]);
		used+=sprintf(ret+used,"%d%c",cigar->lengths[i],cigar->ops[i]);
	}
	//printf("%d vs %d, %d\n",used,string_length,cigar->size);
	assert(used+1==string_length);
	ret[used]='\0';
	return ret;
}

void free_cigar(cigar_t * cigar) {
	if (cigar->size>0) {
		assert(cigar->ops!=NULL);
		assert(cigar->lengths!=NULL);
		free(cigar->ops);
		free(cigar->lengths);
	}
	free(cigar);
}


/*
 * Print given hit.
 *
 */
void
hit_output(struct read_entry * re, struct read_hit * rh, struct read_hit * rh_mp,
	   char ** output1, char ** output2, bool first_in_pair, int* hits, int satisfying_alignments)
/*
 * This function sets the strings output1 and output2 to be the output for the current read and if in sam mode its matepair
 * It is capable of outputting regular shrimp output, pretty print output, and sam output
 *
 * re is the read_entry for the current read
 * rh is the read_hit for the current read
 * rh_mp is the read_hit for the current reads mate pair
 *
 * output1 is a pointer to a string to be used for the firest line of output
 * output1 is a pointer to a string to be used for remaining output lines (used for pretty print)
 *
 * paired is true if this read is paired
 * first is true if this is the first read in the pair
 *
 */
{
  assert(re !=NULL);
  if(!Eflag) {
  	assert(rh != NULL);
  	assert(rh->sfrp != NULL);
  	*output1 = output_normal(re->name, contig_names[rh->cn], rh->sfrp,
			   genome_len[rh->cn], shrimp_mode == MODE_COLOUR_SPACE, re->read[rh->st],
			   re->read_len, re->initbp[rh->st], rh->gen_st, Rflag);
	if (Pflag) {
		//pretty print output
		*output2 = output_pretty(re->name, contig_names[rh->cn], rh->sfrp,
				     genome_contigs[rh->cn], genome_len[rh->cn],
				     (shrimp_mode == MODE_COLOUR_SPACE), re->read[rh->st],
				     re->read_len, re->initbp[rh->st], rh->gen_st);
	}
  } else {
	int thread_id = omp_get_thread_num();
	char ** output_buffer = thread_output_buffer_filled+thread_id;
 	char * output_buffer_end = thread_output_buffer[thread_id] + thread_output_buffer_sizes[thread_id] - 1 + 1;
	while ( (size_t)(output_buffer_end - *output_buffer) < thread_output_buffer_safety) { 
		//fprintf(stderr,"%d incrementing buffer, free space %llu\n",thread_id,output_buffer_end - *output_buffer );
		size_t new_size = thread_output_buffer_sizes[thread_id]+thread_output_buffer_increment;
		size_t filled = thread_output_buffer_filled[thread_id]-thread_output_buffer[thread_id];
		//fprintf(stderr, "there are %llu bytes used\n",filled);
		thread_output_buffer[thread_id]=(char*)realloc(thread_output_buffer[thread_id],new_size);
		if (thread_output_buffer[thread_id]==NULL) {
			fprintf(stderr,"Hit output : realloc failed!\n");
			exit(1);
		}
		thread_output_buffer_sizes[thread_id]=new_size;
		thread_output_buffer_filled[thread_id]=thread_output_buffer[thread_id]+filled;
		output_buffer = thread_output_buffer_filled+thread_id;
		output_buffer_end = thread_output_buffer[thread_id] + thread_output_buffer_sizes[thread_id] - 1 + 1;
	}	
	//TODO change this size?
	//int buffer_size=MAX(longest_read_len,1000)*8;
	//*output1 = (char *)xmalloc(sizeof(char *)*buffer_size);
	//qname
	char * read_name = re->name;
	char qname[strlen(read_name)+1];
	strcpy(qname,read_name);
	//flag
	int flag;
	//rname
	char const * rname = "*";
	//pos
	int pos=0;
	//mapq
	int mapq=255;
	//cigar
	char * cigar=(char *)"*";
	cigar_t * cigar_binary=NULL;
	//mrnm
	const char * mrnm = "*"; //mate reference name
	//mpos
	int mpos=0;
	//isize
	int isize=0;
	//seq
	assert(shrimp_mode==MODE_COLOUR_SPACE || (signed int)strlen(re->seq)==re->read_len);
	assert(shrimp_mode==MODE_LETTER_SPACE || (signed int)strlen(re->seq)==re->read_len+1);
	char seq[re->read_len+1];
	if (shrimp_mode == MODE_LETTER_SPACE) {
		int i; 
		for (i=0; i<re->read_len; i++) {
			char c = re->seq[i];				
			switch(c) {
				case 'R':
				case 'Y':
				case 'S':
				case 'W':
				case 'K':
				case 'M':
				case 'B':
				case 'D':
				case 'H':
				case 'V':
					seq[i]='N';
					break;
				default:
					if (c>='a') {
						c-=32;
					}
					seq[i]=c;
					break;	
			}
		} 
		assert(i==re->read_len);
		seq[re->read_len]='\0';
	} else {
		seq[0]='*'; seq[1]='\0';
	}
	//qual
	int read_length = re->read_len;
	char qual[read_length+10];
	strcpy(qual,"*");
	//initialize flags	
	bool paired_read = re->paired;
	struct read_entry * re_mp = re->mate_pair;
	assert(!paired_read || re_mp!=NULL);
	bool paired_alignment = paired_read && (rh!=NULL && rh_mp!=NULL); //paired mapping, not paired read!
	//bool proper_pair = (paired_read && !query_unmapped && !mate_unmapped);
	//bool query_unmapped = (re->n_hits[0] + re->n_hits[1])>0 ? false : true;
	bool query_unmapped = (rh==NULL);
	bool mate_unmapped=false;
	bool reverse_strand = false;
	bool reverse_strand_mp = false;
	int genome_end_mp=0; int genome_start_mp=0;
	if (paired_read) {
		int min_read_name_length=MIN(strlen(read_name),strlen(re_mp->name));
		int i;
		for (i=0; i<min_read_name_length; i++) {
			if (read_name[i]==re_mp->name[i]) {
				qname[i]=read_name[i];
			} else {
				break;
			}
		}
		if (i>0 && (qname[i-1]==':' || qname[i-1]=='/')) {
			i--;
		}
		qname[i]='\0';
		//mate_unmapped=(re_mp->n_hits[0]+re_mp->n_hits[1])>0 ? false : true;
		mate_unmapped= (rh_mp==NULL);
		if (!mate_unmapped) {
			//char * read_name_mp = re->name;
			//char * rname_mp = contig_names[rh_mp->cn];
			int read_start_mp = rh_mp->sfrp->read_start+1; //1based
			//int read_length_mp = re_mp->read_len;
			int read_end_mp = read_start_mp + rh_mp->sfrp->rmapped -1; //1base
			int genome_length_mp = genome_len[rh_mp->cn];
			reverse_strand_mp = (rh_mp->gen_st ==1);
			if (!reverse_strand_mp) {
				genome_start_mp = rh_mp->sfrp->genome_start+1; // 0 based -> 1 based
			} else {
				int genome_right_most_coordinate = genome_length_mp - rh_mp->sfrp->genome_start;
				//rh->sfrp->deletions is deletions in the reference
				// This is when the read has extra characters that dont match into ref
				genome_start_mp = genome_right_most_coordinate - (read_end_mp - read_start_mp - rh_mp->sfrp->deletions + rh_mp->sfrp->insertions);
			}
			genome_end_mp=genome_start_mp+rh_mp->sfrp->gmapped-1;
			mpos=genome_start_mp;
			mrnm = contig_names[rh_mp->cn];
		}
	}
	bool second_in_pair = (paired_read && !first_in_pair);
	bool primary_alignment = false;
	bool platform_quality_fail = false;
	bool pcr_duplicate = false;
	int stored_alignments = MIN(num_outputs,satisfying_alignments); //IH
	//if the read has no mapping or if not in half_paired mode and the mate has no mapping
	if (query_unmapped || (!sam_half_paired && paired_read && mate_unmapped)) {
		mapq=0;
		if (Qflag && shrimp_mode == MODE_LETTER_SPACE ){
			strcpy(qual,re->qual);
		}
		flag = 
			( paired_read ?  0x0001 : 0) |
			( paired_alignment ? 0x0002 : 0) |
			( query_unmapped ? 0x0004 : 0) |
			( mate_unmapped ? 0x0008 : 0) |
			( reverse_strand ? 0x0010 : 0) |
			( reverse_strand_mp ? 0x0020 : 0) |
			( first_in_pair ? 0x0040 : 0) |
			( second_in_pair ? 0x0080 : 0) |
			( primary_alignment ? 0x0100 : 0) |
			( platform_quality_fail ? 0x0200 : 0) |
			( pcr_duplicate ? 0x0400 : 0);
		//char *extra = *output1 + sprintf(*output1,"%s\t%i\t%s\t%u\t%i\t%s\t%s\t%u\t%i\t%s\t%s",
		//	qname,flag,rname,pos,mapq,cigar,mrnm,mpos,
		//	isize,seq,qual);
		*output_buffer += snprintf(*output_buffer,output_buffer_end-*output_buffer,
			"%s\t%i\t%s\t%u\t%i\t%s\t%s\t%u\t%i\t%s\t%s",
			qname,flag,rname,pos,mapq,cigar,mrnm,mpos,
			isize,seq,qual);
		if (shrimp_mode == MODE_COLOUR_SPACE) {
			if (Qflag) {
				//extra = extra + sprintf(extra,"\tCQ:Z:%s",re->qual);
				*output_buffer += snprintf(*output_buffer,output_buffer_end-*output_buffer,"\tCQ:Z:%s",re->qual);
			} else {
				//extra = extra + sprintf(extra,"\tCQ:Z:%s",qual);
				*output_buffer += snprintf(*output_buffer,output_buffer_end-*output_buffer,"\tCQ:Z:%s",qual);
			}
			//extra = extra + sprintf(extra, "\tCS:Z:%s",re->seq);
			*output_buffer += snprintf(*output_buffer,output_buffer_end-*output_buffer, "\tCS:Z:%s",re->seq);
		}
		if (sam_r2) {
			if (shrimp_mode == MODE_COLOUR_SPACE) {
				//extra = extra + sprintf(extra, "\tX2:Z:%s",re_mp->seq);
				*output_buffer += snprintf(*output_buffer,output_buffer_end-*output_buffer, "\tX2:Z:%s",re_mp->seq);
			} else {
				//extra = extra + sprintf(extra, "\tR2:Z:%s",re_mp->seq);
				*output_buffer += snprintf(*output_buffer,output_buffer_end-*output_buffer, "\tR2:Z:%s",re_mp->seq);
			}
		}
		if (sam_read_group_name!=NULL ){
			//extra+=sprintf(extra,"\tRG:Z:%s",sam_read_group_name);
			*output_buffer += snprintf(*output_buffer,output_buffer_end-*output_buffer,"\tRG:Z:%s",sam_read_group_name);
		}
		*output_buffer += snprintf(*output_buffer,output_buffer_end-*output_buffer,"\n");	
		assert(satisfying_alignments==0);
		assert(stored_alignments==0);
		//assert(hits[0]==0);
		//assert(hits[1]==0);
		//assert(hits[2]==0);
		//assert(hits[3]==0);
		//extra = extra + sprintf(extra,"\tH0:i:%d\tH1:i:%d\tH2:i:%d\tNH:i:%d\tIH:i:%d",hits[0],hits[1],hits[2],found_alignments,stored_alignments);
		return;
	}
	assert(rh!=NULL);
	assert( !paired_read || (!query_unmapped && !mate_unmapped) || sam_half_paired);
	//start filling in the fields
	rname = contig_names[rh->cn];
	reverse_strand = (rh->gen_st == 1);
	int read_start = rh->sfrp->read_start+1; //1based
	int read_end = read_start + rh->sfrp->rmapped -1; //1base
	int genome_length = genome_len[rh->cn];
	cigar_binary = make_cigar(read_start,read_end,read_length,rh->sfrp->qralign,rh->sfrp->dbalign);

	int qralign_length=strlen(rh->sfrp->qralign);
	int i,j=0;
	int seq_length=0;
	if (shrimp_mode == MODE_LETTER_SPACE ) {
		j=read_start-1;
		seq_length=re->read_len;
	} else if (shrimp_mode == MODE_COLOUR_SPACE ) {
		seq_length=read_end-read_start+1;
	}
	assert(seq_length<=re->read_len);
	for(i=0;i<qralign_length;i++) {
		char c=rh->sfrp->qralign[i];
		if (c!='-') { 
			if (c>='a') {
				c-=32;
			}
			if (c!='A' && c!='a' &&
				c!='G' && c!='g' &&
				c!='C' && c!='c' &&
				c!='T' && c!='t' &&
				c!='N' && c!='n') {
				//see if we can figure out what its suppose to be
				c='N';
				if (rh->sfrp->dbalign[i]!='-') {
					char r = rh->sfrp->dbalign[i];
					if (r>='a') {
						r-=32;
					}
					switch (r) {
						case 'A':
							if (c=='R' || c=='W' || c=='M' || c=='D' || c=='H' || c=='V' )
								c='A';
							break;
						case 'C':
							if (c=='Y' || c=='S' || c=='M' || c=='B' || c=='H' || c=='V')
								c='C';
							break;
						case 'G':
							if (c=='R' || c=='S' || c=='K' || c=='B' || c=='D' || c=='V')
								c='G';
							break;
						case 'T':
							if (c=='Y' || c=='W' || c=='K' || c=='B' || c=='D' || c=='H') 
								c='T';
							break;
						default: 
							fprintf(stderr,"There has been an error in printing an alignment, %c\n",r);
							break;
					}
				}
			}
			seq[j++]=c;
		}
	}
	if (shrimp_mode == MODE_LETTER_SPACE ) {
		assert(j+(re->read_len-read_end)==seq_length);
		seq[j+(re->read_len-read_end)]='\0';
	} else if (shrimp_mode == MODE_COLOUR_SPACE) {
		assert(j==seq_length);
		seq[seq_length]='\0';
	}

	//if its letter space need to reverse the qual string if its backwards
	if (shrimp_mode == MODE_LETTER_SPACE) {
		//seq=re->seq;
		if (Qflag) {
			if (!reverse_strand) {
				strcpy(qual,re->qual);
			} else {
				int qual_len = strlen(re->qual); //not same as read_len, for color space reads... apperently.....
				assert((qual_len+1)<(re->read_len+10));
				int i;
				for (i=0; i<qual_len; i++) {
					qual[(qual_len-1)-i]=re->qual[i];
				}
				qual[qual_len]='\0';
			}
		}
	//else in colour space dont print a qual string
	//but get the seq differently and change 'S' to 'H' in cigar
	} else if (shrimp_mode == MODE_COLOUR_SPACE) {
		//also change 'S' in cigar to 'H'
		//clip the qual values
		for (i=0; i<cigar_binary->size; i++) {
			if (cigar_binary->ops[i]=='S') {
				cigar_binary->ops[i]='H';
			}
		}
		if (Bflag && Qflag) {
			int read_length=(read_end-read_start+1);
			for (i=0; i<read_length; i++) {
				qual[i]=re->qual[i+read_start-1];
			}
			qual[i]='\0';
			for (i=0; i<read_length-1; i++) {
				//this is different from bfast
				//qralign is already clipped! i.e. doesn't have clipped stuff and is
				//read orientation (not always on positive reference strand!)
				int first_position_mismatch = rh->sfrp->qralign[i] > 96;
				int second_position_mismatch = rh->sfrp->qralign[i+1] > 96;
				int base_qual=0;
				if (first_position_mismatch && second_position_mismatch ) {
					base_qual+=0;
				} else if (first_position_mismatch) {
					base_qual+=qual[i+1]-qual[i];
				} else if (second_position_mismatch) {
					base_qual+=qual[i]-qual[i+1]+33;
				} else {
					base_qual+=qual[i]+qual[i+1]+10-33;
				}
				base_qual=MIN('`',MAX(base_qual,'"'));
				qual[i]=base_qual;
			}
			if (reverse_strand) {
				for (i=0; i<read_length/2; i++) {
					char temp = qual[i];
					qual[i]=qual[read_length-i-1];
					qual[read_length-i-1]=temp;
				}
			}
		}
	}
	//get the pos
	int genome_start;
	if (!reverse_strand) {
		genome_start = rh->sfrp->genome_start+1; // 0 based -> 1 based
	} else {
		int genome_right_most_coordinate = genome_length - rh->sfrp->genome_start;
		//rh->sfrp->deletions is deletions in the reference
		// This is when the read has extra characters that dont match into ref
		genome_start = genome_right_most_coordinate - (read_end - read_start - rh->sfrp->deletions + rh->sfrp->insertions);
		char * tmp = (char*)xmalloc(sizeof(char)*(strlen(seq)+1));
		reverse(seq,tmp);
		free(tmp);
		reverse_cigar(cigar_binary);
	}
	int genome_end=genome_start+rh->sfrp->gmapped-1;
	pos=genome_start;
	//make the cigar string
	cigar = make_cigar_string(cigar_binary); 

	//do some stats using matepair
	if (paired_read && !mate_unmapped) {
		assert(rh_mp!=NULL && re_mp!=NULL);

		mrnm = (strcmp(rname,mrnm)==0) ? "=" : mrnm;
		//printf("%d %d %c, %d %d %c\n",genome_start, genome_end,reverse_strand ? '-' : '+' , genome_start_mp, genome_end_mp, reverse_strand_mp ? '-' : '+');
		int fivep = 0;
		int fivep_mp = 0;
		if (reverse_strand){
			fivep = genome_end;
		} else {
			fivep = genome_start - 1;
		}

		if (reverse_strand_mp){
			fivep_mp = genome_end_mp;
		} else {
			fivep_mp = genome_start_mp-1;
		}
		isize = (fivep_mp - fivep);
	}
	flag = 
		( paired_read ?  0x0001 : 0) |
		( paired_alignment ? 0x0002 : 0) |
		( query_unmapped ? 0x0004 : 0) |
		( mate_unmapped ? 0x0008 : 0) |
		( reverse_strand ? 0x0010 : 0) |
		( reverse_strand_mp ? 0x0020 : 0) |
		( first_in_pair ? 0x0040 : 0) |
		( second_in_pair ? 0x0080 : 0) |
		( primary_alignment ? 0x0100 : 0) |
		( platform_quality_fail ? 0x0200 : 0) |
		( pcr_duplicate ? 0x0400 : 0);
	//char *extra = *output1 + sprintf(*output1,"%s\t%i\t%s\t%u\t%i\t%s\t%s\t%u\t%i\t%s\t%s",
	//	qname,flag,rname,pos,mapq,cigar,mrnm,mpos,
	//	isize,seq,qual);
	*output_buffer += snprintf(*output_buffer,output_buffer_end-*output_buffer,"%s\t%i\t%s\t%u\t%i\t%s\t%s\t%u\t%i\t%s\t%s",
		qname,flag,rname,pos,mapq,cigar,mrnm,mpos,
		isize,seq,qual);
	//extra = extra + sprintf(extra,"\tAS:i:%d\tH0:i:%d\tH1:i:%d\tH2:i:%d\tNM:i:%d\tNH:i:%d\tIH:i:%d",rh->sfrp->score,hits[0],hits[1],hits[2],rh->sfrp->mismatches+rh->sfrp->deletions+rh->sfrp->insertions,found_alignments,stored_alignments);
		//MERGESAM DEPENDS ON SCORE BEING FIRST!
	*output_buffer += snprintf(*output_buffer,output_buffer_end-*output_buffer,
		"\tAS:i:%d\tH0:i:%d\tH1:i:%d\tH2:i:%d\tNM:i:%d\tNH:i:%d\tIH:i:%d",
		rh->sfrp->score,hits[0],hits[1],hits[2],rh->sfrp->mismatches+rh->sfrp->deletions+rh->sfrp->insertions,satisfying_alignments,stored_alignments);
	if (shrimp_mode == COLOUR_SPACE){
		//TODO
		//int first_bp = re->initbp[0];
		//printf("%s vs %s\n",readtostr(re->read[0],re->read_len,true,first_bp),re->seq);
		//assert(strcmp(readtostr(re->read[0],re->read_len,true,first_bp),re->seq)==0);
		if (Qflag) {
			//extra = extra + sprintf(extra,"\tCQ:Z:%s",re->qual);
			*output_buffer += snprintf(*output_buffer,output_buffer_end-*output_buffer,"\tCQ:Z:%s",re->qual);
		}
		//extra = extra + sprintf(extra, "\tCS:Z:%s\tCM:i:%d\tXX:Z:%s",re->seq,rh->sfrp->crossovers,rh->sfrp->qralign);
		*output_buffer += snprintf(*output_buffer,output_buffer_end-*output_buffer, "\tCS:Z:%s\tCM:i:%d\tXX:Z:%s",re->seq,rh->sfrp->crossovers,rh->sfrp->qralign);
	} 
	if (sam_r2) {
		if (shrimp_mode == MODE_COLOUR_SPACE) {
			//extra = extra + sprintf(extra, "\tX2:Z:%s",re_mp->seq);
			*output_buffer += snprintf(*output_buffer,output_buffer_end-*output_buffer, "\tX2:Z:%s",re_mp->seq);
		} else {
			//extra = extra + sprintf(extra, "\tR2:Z:%s",re_mp->seq);
			*output_buffer += snprintf(*output_buffer,output_buffer_end-*output_buffer, "\tR2:Z:%s",re_mp->seq);
		}
	}
	
	if (sam_read_group_name!=NULL) {
			//extra+=sprintf(extra,"\tRG:Z:%s",sam_read_group_name);
			*output_buffer+=snprintf(*output_buffer,output_buffer_end-*output_buffer,"\tRG:Z:%s",sam_read_group_name);
	}
	if (cigar_binary!=NULL) {
		free_cigar(cigar_binary);
		free(cigar);
	}
	*output_buffer += snprintf(*output_buffer,output_buffer_end-*output_buffer,"\n");	

    //to calculate the insert size we need to find the five' end of the reads
/*
	    inp.score);
    if(re_mp != NULL){
    	free(name);
    }
    free(read);
    free(cigar);
    format_free(fsp);*/

  }
}


/*
 * Free memory allocated by this read.
 */
void 
read_free(read_entry * re)
{
  free(re->name);
  free(re->seq);
  if (re->orig_seq!=re->seq) {
	free(re->orig_seq);
  }
  if (Qflag) {
        assert(re->qual!=NULL);
        free(re->qual);
	if (re->qual!=re->orig_qual) {
		free(re->orig_qual);
	}
        assert(re->plus_line!=NULL);
        free(re->plus_line);
  }
}


void
read_free_hit_list(struct read_entry * re)
{
  for (int st = 0; st < 2; st++) {
    if (re->hits[st] != NULL) {
      free(re->hits[st]);
      re->hits[st] = NULL;
      re->n_hits[st] = 0;
    }
  }
}


void
read_free_anchor_list(struct read_entry * re)
{
  for (int st = 0; st < 2; st++) {
    if (re->anchors[st] != NULL) {
      free(re->anchors[st]);
      re->anchors[st] = NULL;
      re->n_anchors[st] = 0;
    }
  }
}

void
read_free_full(struct read_entry * re)
{
  read_free(re);
  free(re->read[0]);
  free(re->read[1]);

  free(re->mapidx[0]);
  free(re->mapidx[1]);
  free(re->anchors[0]);
  free(re->anchors[1]);
  free(re->hits[0]);
  free(re->hits[1]);

  if (re->n_ranges > 0) {
    free(re->ranges);
  }
}


/*
 * Reverse read.
 *
 * This is useful for dealing with the various possibilities for matepair orientation in a unified way.
 */
static void
read_reverse(struct read_entry * re) {
  uint32_t * tmp1 = re->read[0];
  re->read[0] = re->read[1];
  re->read[1] = tmp1;
  
  int tmp2 = re->initbp[0];
  re->initbp[0] = re->initbp[1];
  re->initbp[1] = tmp2;

  re->input_strand = 1 - re->input_strand;
}

static uint
get_contig_number_from_name(char const * c)
{
  int cn;
  // hack: accept '>' for cn=0
  if (*c == '>')
    return 0;
  for (cn = 0; cn < num_contigs; cn++) {
    if (!strcmp(c, contig_names[cn]))
      break;
  }
  return cn;
}

/*
 * Compute range limitations for this read.
 */
static void
read_compute_ranges(struct read_entry * re)
{
  char * r, * r_save, * c, * c_save;
  uint g_start, g_end;
  int cn, st;

  assert(re->range_string != NULL);

  for (r = strtok_r(re->range_string, " ", &r_save); r != NULL; r = strtok_r(NULL, " ", &r_save))
    {
      c = strtok_r(r, ",", &c_save); // contig name
      if (c == NULL)
	continue;
      cn = get_contig_number_from_name(c);
      if (cn >= num_contigs)
	continue;

      c = strtok_r(NULL, ",", &c_save); // strand
      if (c == NULL)
	continue;
      if (*c == '+')
	st = 0;
      else if (*c == '-')
	st = 1;
      else
	continue;

      c = strtok_r(NULL, ",", &c_save); // g_start
      if (c == NULL)
	continue;
      g_start = (uint)atoll(c);
      if (g_start == 0)
	continue;
      g_start--;

      c = strtok_r(NULL, ",", &c_save); // g_end
      if (c == NULL)
	continue;
      g_end = (uint)atoll(c);
      if (g_end == 0)
	continue;
      g_end--;

      re->n_ranges++;
      re->ranges = (struct range_restriction *)xrealloc(re->ranges, re->n_ranges * sizeof(re->ranges[0]));
      re->ranges[re->n_ranges - 1].cn = cn;
      re->ranges[re->n_ranges - 1].st = st;
      re->ranges[re->n_ranges - 1].g_start = g_start;
      re->ranges[re->n_ranges - 1].g_end = g_end;
    }
}

/* Trim a read */
static void trim_read(struct read_entry * re) {
	re->orig_seq=strdup(re->seq);
	if (Qflag) {
		re->orig_qual=strdup(re->qual);
	}
	//handle colour space too!
	int length=strlen(re->seq);
	int i;
	for (i=0; i<length-trim_end-trim_front; i++) {
		re->seq[i]=re->seq[i+trim_front];
		if (Qflag) {
			re->qual[i]=re->qual[i+trim_front];
		}
	}
	memset(re->seq+i,0,17);
	if (Qflag) {
		memset(re->qual+i,0,17);
	}
	return;
}

/*
 * Launch the threads that will scan the reads
 */
static bool
launch_scan_threads(){
  fasta_t fasta = NULL, left_fasta = NULL, right_fasta = NULL;

  //open the fasta file and check for errors
  if (single_reads_file) {  
    fasta = fasta_open(reads_filename, shrimp_mode, Qflag);
    if (fasta == NULL) {
      fprintf(stderr, "error: failed to open read file [%s]\n", reads_filename);
      return (false);
    } else {
      fprintf(stderr, "- Processing read file [%s]\n", reads_filename);
    }
  } else {
    left_fasta = fasta_open(left_reads_filename,shrimp_mode,Qflag);
    if (left_fasta == NULL) {
      fprintf(stderr, "error: failed to open read file [%s]\n", left_reads_filename);
      return (false);
    }
    right_fasta = fasta_open(right_reads_filename,shrimp_mode,Qflag);
    if (right_fasta == NULL) {
      fprintf(stderr, "error: failed to open read file [%s]\n", right_reads_filename);
      return (false);
    }
    // WHY? the code above sets both ->space to shrimp_mode anyhow
    //if (right_fasta->space != left_fasta->space) {
    //  fprintf(stderr,"error: when using -1 and -2, both files must be either only colour space or only letter space!\n");
    //  return (false);
    //}
    fasta=left_fasta;
    fprintf(stderr, "- Processing read files [%s , %s]\n", left_reads_filename,right_reads_filename);
  }

  if ((pair_mode != PAIR_NONE || !single_reads_file) && (chunk_size%2)!=0) {
    fprintf(stderr,"error: when in paired mode or using options -1 and -2, then thread chunk size must be even\n"); 
  }

  bool read_more = true, more_in_left_file = true, more_in_right_file=true;

  /* initiate the thread buffers */
  thread_output_buffer_sizes = (size_t *)xcalloc_m(sizeof(size_t) * num_threads, "thread_output_buffer_sizes");
  thread_output_buffer_filled = (char * *)xcalloc_m(sizeof(char *) * num_threads, "thread_output_buffer_filled");
  thread_output_buffer = (char * *)xcalloc_m(sizeof(char *) * num_threads, "thread_output_buffer");
  thread_output_buffer_chunk = (unsigned int *)xcalloc_m(sizeof(unsigned int) * num_threads, "thread_output_buffer_chunk");
  
  unsigned int current_thread_chunk = 1;
  unsigned int next_chunk_to_print = 1;
  struct heap_out h; 
  heap_out_init(&h, thread_output_heap_capacity );

#pragma omp parallel shared(read_more,more_in_left_file,more_in_right_file, fasta) num_threads(num_threads)
  {
    int thread_id = omp_get_thread_num();
    struct read_entry * re_buffer;
    int load, i;
    uint64_t before;
    re_buffer = (struct read_entry *)xcalloc_m(chunk_size * sizeof(re_buffer[0]), "re_buffer");

    while (read_more) {
      before = rdtsc();

      //Read in this threads 'chunk'
#pragma omp critical (fill_reads_buffer)
      {
	thread_output_buffer_chunk[thread_id]=current_thread_chunk++;
	wait_ticks[omp_get_thread_num()] += rdtsc() - before;

	load = 0;
	assert(chunk_size>2);
	while (read_more && ((single_reads_file && load < chunk_size) || (!single_reads_file && load < chunk_size-1))) {
	  //if (!fasta_get_next_with_range(fasta, &re_buffer[load].name, &re_buffer[load].seq, &re_buffer[load].is_rna,
	  //				 &re_buffer[load].range_string, &re_buffer[load].qual))
	  if (single_reads_file) { 
	    if (fasta_get_next_read_with_range(fasta, &re_buffer[load])) {
	      load++;
	    } else { 
	      read_more = false;
	    }
	  } else {
	    //read from the left file
	    if (fasta_get_next_read_with_range(left_fasta, &re_buffer[load])) {
	      load++;
	    } else {
	      more_in_left_file = false;
	    }
	    //read from the right file
	    if (fasta_get_next_read_with_range(right_fasta, &re_buffer[load])) {
	      load++;
	    } else {
	      more_in_right_file = false;
	    }
	    //make sure that one is not smaller then the other
	    if (more_in_left_file != more_in_right_file) {
	      fprintf(stderr,"error: when using options -1 and -2, both files specified must have the same number of entries\n");
	      exit(1);
	    }
	    //keep reading?
	    read_more = more_in_left_file && more_in_right_file; 
	  }
	}

	nreads += load;

	// progress reporting
	if (progress > 0) {
	  nreads_mod += load;
	  if (nreads_mod >= progress)
	    fprintf(stderr, "\r%lld", nreads);
	  nreads_mod %= progress;
	}
      } // end critical section

      if (pair_mode != PAIR_NONE)
	assert(load % 2 == 0); // read even number of reads

      thread_output_buffer_sizes[thread_id] = thread_output_buffer_initial;
      thread_output_buffer[thread_id] = (char *)xmalloc_m(sizeof(char) * thread_output_buffer_sizes[thread_id], "thread_buffer");
      thread_output_buffer_filled[thread_id] = thread_output_buffer[thread_id];
      thread_output_buffer[thread_id][0] = '\0';

      for (i = 0; i < load; i++) {
	// if running in paired mode and first foot is ignored, ignore this one, too
	if (pair_mode != PAIR_NONE && i % 2 == 1 && re_buffer[i-1].ignore) {
	  read_free(&re_buffer[i-1]);
	  read_free(&re_buffer[i]);
	  continue;
	}

	//if (!(strcspn(re_buffer[i].seq, "nNxX.") == strlen(re_buffer[i].seq))) {
	//  if (pair_mode != PAIR_NONE && i % 2 == 1) {
	//    read_free(re_buffer+i-1);
	//    read_free(re_buffer+i);
	//  }
	//  continue;
	//}
	
	//Trim the reads
	if (trim) {
	  if (pair_mode != PAIR_NONE) {
	    if (trim_first) {
	      trim_read(&re_buffer[i-1]);
	    }
	    if (trim_second) {
	      trim_read(&re_buffer[i]);
	    }
	  } else {
	    trim_read(&re_buffer[i]);
	  }
	}	

	re_buffer[i].ignore = false;
	re_buffer[i].read[0] = fasta_sequence_to_bitfield(fasta, re_buffer[i].seq);
	re_buffer[i].read_len = strlen(re_buffer[i].seq);
	re_buffer[i].max_n_kmers = re_buffer[i].read_len - min_seed_span + 1;
	re_buffer[i].min_kmer_pos = 0;
	if (shrimp_mode == MODE_COLOUR_SPACE){
	  re_buffer[i].read_len--;
	  re_buffer[i].max_n_kmers -= 2; // 1st color always discarded from kmers
	  re_buffer[i].min_kmer_pos = 1;
	  re_buffer[i].initbp[0] = fasta_get_initial_base(shrimp_mode,re_buffer[i].seq);
	  re_buffer[i].initbp[1] = re_buffer[i].initbp[0];
	  re_buffer[i].read[1] = reverse_complement_read_cs(re_buffer[i].read[0], (int8_t)re_buffer[i].initbp[0],
							    (int8_t)re_buffer[i].initbp[1],
							    re_buffer[i].read_len, re_buffer[i].is_rna);
	} else {
	  re_buffer[i].read[1] = reverse_complement_read_ls(re_buffer[i].read[0], re_buffer[i].read_len, re_buffer[i].is_rna);
	}

	//Check if we can actually use this read
	if (re_buffer[i].max_n_kmers < 0
	    || re_buffer[i].read_len > longest_read_len) {
	  if (re_buffer[i].max_n_kmers < 0) {
	    fprintf(stderr, "warning: skipping read [%s]; smaller then any seed!\n",
		    re_buffer[i].name);
	    re_buffer[i].max_n_kmers=1;
	  } else {
	    fprintf(stderr, "warning: skipping read [%s]; it has length %d, maximum allowed is %d. Use --longest-read ?\n",
		    re_buffer[i].name, re_buffer[i].read_len, longest_read_len);
	  }
	  if (pair_mode == PAIR_NONE) {
	    read_free_full(&re_buffer[i]);
	  } else if (i%2 == 1) {
	    read_free_full(&re_buffer[i-1]);
	    read_free_full(&re_buffer[i]);
	  } else {
	    re_buffer[i].ignore=true;
	  }
	  continue;	
	}

	re_buffer[i].window_len = (uint16_t)abs_or_pct(window_len,re_buffer[i].read_len);
	re_buffer[i].input_strand = 0;

	re_buffer[i].mapidx[0] = NULL;
	re_buffer[i].mapidx[1] = NULL;
	re_buffer[i].anchors[0] = NULL;
	re_buffer[i].anchors[1] = NULL;
	re_buffer[i].hits[0] = NULL;
	re_buffer[i].hits[1] = NULL;
	re_buffer[i].ranges = NULL;
	re_buffer[i].n_ranges = 0;
	re_buffer[i].final_matches = 0;
	if (re_buffer[i].range_string != NULL) {
	  read_compute_ranges(&re_buffer[i]);
	  free(re_buffer[i].range_string);
	  re_buffer[i].range_string = NULL;
	}
	//free(re_buffer[i].seq);

	// time to do some mapping!
	if (pair_mode == PAIR_NONE)
	  {
	    handle_read(&re_buffer[i]);
	  }
	else if (i % 2 == 1)
	  {
	    if (pair_reverse[pair_mode][0])
	      read_reverse(&re_buffer[i-1]);
	    if (pair_reverse[pair_mode][1])
	      read_reverse(&re_buffer[i]);
	    re_buffer[i-1].paired=true;
	    re_buffer[i-1].first_in_pair=true;
	    re_buffer[i-1].mate_pair=&re_buffer[i];
	    re_buffer[i].paired=true;
	    re_buffer[i].first_in_pair=false;
	    re_buffer[i].mate_pair=&re_buffer[i-1];
	    handle_readpair(&re_buffer[i-1], &re_buffer[i]);
	  }
      }

      //fprintf(stdout,"%s",thread_output_buffer[thread_id]);
#pragma omp critical
      {
	struct heap_out_elem tmp;
	tmp.key = thread_output_buffer_chunk[thread_id];
	tmp.rest = thread_output_buffer[thread_id];
	thread_output_buffer[thread_id] = NULL;	
	heap_out_insert(&h, &tmp);
	heap_out_get_min(&h, &tmp);
	while (h.load > 0 && tmp.key == next_chunk_to_print) {
	  heap_out_extract_min(&h, &tmp);
	  fprintf(stdout, "%s", tmp.rest);
	  free(tmp.rest);
	  next_chunk_to_print++;
	}
      }
    }

    free(re_buffer);
  } // end parallel section

  if (progress > 0)
    fprintf(stderr, "\n");

  struct heap_out_elem tmp;
  while (h.load>0) {
    heap_out_extract_min(&h,&tmp);
    fprintf(stdout,"%s",tmp.rest);
    free(tmp.rest);
  }
  heap_out_destroy(&h);
  free(thread_output_buffer_sizes);
  free(thread_output_buffer_filled);
  free(thread_output_buffer);
  free(thread_output_buffer_chunk);
  if (single_reads_file) {
    fasta_close(fasta);
  } else {
    fasta_close(left_fasta);
    fasta_close(right_fasta);
  }
  return true;
}


static void
print_insert_histogram()
{
  int i;
  for (i = 0; i < 100; i++) {
    fprintf(stderr, "[%d-%d]: %.2f%%\n",
	    min_insert_size + i * insert_histogram_bucket_size,
	    min_insert_size + (i + 1) * insert_histogram_bucket_size - 1,
	    total_paired_matches == 0? 0 : ((double)insert_histogram[i] / (double)total_paired_matches) * 100);
  }
}


static void
print_statistics()
{
  static char const my_tab[] = "    ";

	uint64_t f1_invocs[num_threads], f1_cells[num_threads], f1_ticks[num_threads];
	double f1_secs[num_threads], f1_cellspersec[num_threads];
	uint64_t f1_total_invocs = 0, f1_total_cells = 0;
	double f1_total_secs = 0, f1_total_cellspersec = 0;
	uint64_t f1_calls_bypassed = 0;

	uint64_t f2_invocs[num_threads], f2_cells[num_threads], f2_ticks[num_threads];
	double f2_secs[num_threads], f2_cellspersec[num_threads];
	uint64_t f2_total_invocs = 0, f2_total_cells = 0;
	double f2_total_secs = 0, f2_total_cellspersec = 0;

	double scan_secs[num_threads], readload_secs[num_threads];
	double total_scan_secs = 0, total_wait_secs = 0, total_readload_secs = 0;

	double hz;
	uint64_t fasta_load_ticks;
	fasta_stats_t fs;

	fs = fasta_stats();
	fasta_load_ticks = fs->total_ticks;
	free(fs);
	hz = cpuhz();

#pragma omp parallel num_threads(num_threads) shared(hz)
	{
	  int tid = omp_get_thread_num();

	  f1_stats(&f1_invocs[tid], &f1_cells[tid], &f1_ticks[tid], NULL);

	  f1_secs[tid] = (double)f1_ticks[tid] / hz;
	  f1_cellspersec[tid] = (double)f1_cells[tid] / f1_secs[tid];
	  if (isnan(f1_cellspersec[tid]))
	    f1_cellspersec[tid] = 0;

	  if (shrimp_mode == MODE_COLOUR_SPACE)
	    sw_full_cs_stats(&f2_invocs[tid], &f2_cells[tid], &f2_ticks[tid]);
	  else
	    sw_full_ls_stats(&f2_invocs[tid], &f2_cells[tid], &f2_ticks[tid]);

	  f2_secs[tid] = (double)f2_ticks[tid] / hz;
	  f2_cellspersec[tid] = (double)f2_cells[tid] / f2_secs[tid];
	  if (isnan(f2_cellspersec[tid]))
	    f2_cellspersec[tid] = 0;

	  scan_secs[tid] = ((double)scan_ticks[tid] / hz) - f1_secs[tid] - f2_secs[tid];
	  scan_secs[tid] = MAX(0, scan_secs[tid]);
	  readload_secs[tid] = ((double)total_work_usecs / 1.0e6) - ((double)scan_ticks[tid] / hz) - ((double)wait_ticks[tid] / hz);
	}
	f1_stats(NULL, NULL, NULL, &f1_calls_bypassed);

	fprintf(stderr, "\nStatistics:\n");

	fprintf(stderr, "%sOverall:\n", my_tab);
	fprintf(stderr, "%s%s%-24s" "%.2f seconds\n", my_tab, my_tab,
		"Load Genome Time:", (double)map_usecs / 1.0e6);
	fprintf(stderr, "%s%s%-24s" "%.2f seconds\n", my_tab, my_tab,
		"Read Mapping Time:", (double)total_work_usecs / 1.0e6);

	fprintf(stderr, "\n");

	int i;
	for(i = 0; i < num_threads; i++){
	  total_scan_secs += scan_secs[i];
	  total_readload_secs += readload_secs[i];
	  total_wait_secs += (double)wait_ticks[i] / hz;

	  f1_total_secs += f1_secs[i];
	  f1_total_invocs += f1_invocs[i];
	  f1_total_cells += f1_cells[i];

	  f2_total_secs += f2_secs[i];
	  f2_total_invocs += f2_invocs[i];
	  f2_total_cells += f2_cells[i];
	}
	f1_total_cellspersec = f1_total_secs == 0? 0 : (double)f1_total_cells / f1_total_secs;
	f2_total_cellspersec = f2_total_secs == 0? 0 : (double)f2_total_cells / f2_total_secs;

	if (Dflag) {
	  fprintf(stderr, "%sPer-Thread Stats:\n", my_tab);
	  fprintf(stderr, "%s%s" "%11s %9s %9s %25s %25s %9s\n", my_tab, my_tab,
		  "", "Read Load", "Scan", "Vector SW", "Scalar SW", "Wait");
	  fprintf(stderr, "%s%s" "%11s %9s %9s %15s %9s %15s %9s %9s\n", my_tab, my_tab,
		  "", "Time", "Time", "Invocs", "Time", "Invocs", "Time", "Time");
	  fprintf(stderr, "\n");
	  for(i = 0; i < num_threads; i++) {
	    fprintf(stderr, "%s%s" "Thread %-4d %9.2f %9.2f %15s %9.2f %15s %9.2f %9.2f\n", my_tab, my_tab,
		    i, readload_secs[i], scan_secs[i], comma_integer(f1_invocs[i]), f1_secs[i],
		    comma_integer(f2_invocs[i]), f2_secs[i], (double)wait_ticks[i] / hz);
	  }
	  fprintf(stderr, "\n");
	}

	fprintf(stderr, "%sSpaced Seed Scan:\n", my_tab);
	fprintf(stderr, "%s%s%-24s" "%.2f seconds\n", my_tab, my_tab,
		"Run-time:", total_scan_secs);

	fprintf(stderr, "\n");

	fprintf(stderr, "%sVector Smith-Waterman:\n", my_tab);
	fprintf(stderr, "%s%s%-24s" "%.2f seconds\n", my_tab, my_tab,
		"Run-time:", f1_total_secs);
	fprintf(stderr, "%s%s%-24s" "%s\n", my_tab, my_tab,
		"Invocations:", comma_integer(f1_total_invocs));
	fprintf(stderr, "%s%s%-24s" "%s\n", my_tab, my_tab,
		"Bypassed Calls:", comma_integer(f1_calls_bypassed));
	fprintf(stderr, "%s%s%-24s" "%.2f million\n", my_tab, my_tab,
		"Cells Computed:", (double)f1_total_cells / 1.0e6);
	fprintf(stderr, "%s%s%-24s" "%.2f million\n", my_tab, my_tab,
		"Cells per Second:", f1_total_cellspersec / 1.0e6);

	fprintf(stderr, "\n");

	fprintf(stderr, "%sScalar Smith-Waterman:\n", my_tab);
	fprintf(stderr, "%s%s%-24s" "%.2f seconds\n", my_tab, my_tab,
		"Run-time:", f2_total_secs);
	fprintf(stderr, "%s%s%-24s" "%s\n", my_tab, my_tab,
		"Invocations:", comma_integer(f2_total_invocs));
	fprintf(stderr, "%s%s%-24s" "%.2f million\n", my_tab, my_tab,
		"Cells Computed:", (double)f2_total_cells / 1.0e6);
	fprintf(stderr, "%s%s%-24s" "%.2f million\n", my_tab, my_tab,
		"Cells per Second:", f2_total_cellspersec / 1.0e6);

	fprintf(stderr, "\n");

	fprintf(stderr, "%sMiscellaneous Totals:\n", my_tab);
	fprintf(stderr, "%s%s%-24s" "%.2f seconds\n", my_tab, my_tab,
		"Fasta Lib Time:", (double)fasta_load_ticks / hz);
	fprintf(stderr, "%s%s%-24s" "%.2f seconds\n", my_tab, my_tab,
		"Read Load Time:", total_readload_secs);
	fprintf(stderr, "%s%s%-24s" "%.2f seconds\n", my_tab, my_tab,
		"Wait Time:", total_wait_secs);

	fprintf(stderr, "\n");

	fprintf(stderr, "%sGeneral:\n", my_tab);
	if (pair_mode == PAIR_NONE)
	  {
	    fprintf(stderr, "%s%s%-24s" "%s    (%.4f%%)\n", my_tab, my_tab,
		    "Reads Matched:",
		    comma_integer(total_reads_matched),
		    (nreads == 0) ? 0 : ((double)total_reads_matched / (double)nreads) * 100);
	    fprintf(stderr, "%s%s%-24s" "%s    (%.4f%%)\n", my_tab, my_tab,
		    "Reads Dropped:",
		    comma_integer(total_reads_dropped),
		    (nreads == 0) ? 0 : ((double)total_reads_dropped / (double)nreads) * 100);
	    fprintf(stderr, "%s%s%-24s" "%s\n", my_tab, my_tab,
		    "Total Matches:",
		    comma_integer(total_single_matches));
	    fprintf(stderr, "%s%s%-24s" "%.2f\n", my_tab, my_tab,
		    "Avg Hits/Matched Read:",
		    (total_reads_matched == 0) ? 0 : ((double)total_single_matches / (double)total_reads_matched));
	    fprintf(stderr, "%s%s%-24s" "%s\n", my_tab, my_tab,
		    "Duplicate Hits Pruned:",
		    comma_integer(total_dup_single_matches));
	  }
	else // paired hits
	  {
	    fprintf(stderr, "%s%s%-40s" "%s    (%.4f%%)\n", my_tab, my_tab,
		    "Pairs Matched:",
		    comma_integer(total_pairs_matched),
		    (nreads == 0) ? 0 : ((double)total_pairs_matched / (double)(nreads/2)) * 100);
	    fprintf(stderr, "%s%s%-40s" "%s    (%.4f%%)\n", my_tab, my_tab,
		    "Pairs Dropped:",
		    comma_integer(total_pairs_dropped),
		    (nreads == 0) ? 0 : ((double)total_pairs_dropped / (double)(nreads/2)) * 100);
	    fprintf(stderr, "%s%s%-40s" "%s\n", my_tab, my_tab,
		    "Total Paired Matches:",
		    comma_integer(total_paired_matches));
	    fprintf(stderr, "%s%s%-40s" "%.2f\n", my_tab, my_tab,
		    "Avg Matches/Pair Matched:",
		    (total_pairs_matched == 0) ? 0 : ((double)total_paired_matches / (double)total_pairs_matched));
	    fprintf(stderr, "%s%s%-40s" "%s\n", my_tab, my_tab,
		    "Duplicate Paired Matches Pruned:",
		    comma_integer(total_dup_paired_matches));

	    /*
	    fprintf(stderr, "\n");

	    fprintf(stderr, "%s%s%-40s" "%s    (%.4f%%)\n", my_tab, my_tab,
		    "Additional Reads Matched Unpaired:",
		    comma_integer(total_reads_matched),
		    (nreads == 0) ? 0 : ((double)total_reads_matched / (double)nreads) * 100);
	    fprintf(stderr, "%s%s%-40s" "%s\n", my_tab, my_tab,
		    "Total Single Matches:",
		    comma_integer(total_single_matches));
	    fprintf(stderr, "%s%s%-40s" "%.2f\n", my_tab, my_tab,
		    "Avg Matches/Unpaired Matched Read:",
		    (total_reads_matched == 0) ? 0 : ((double)total_single_matches / (double)total_reads_matched));
	    fprintf(stderr, "%s%s%-40s" "%s\n", my_tab, my_tab,
		    "Duplicate Unpaired Matches Pruned:",
		    comma_integer(total_dup_single_matches));
	    */
	  }

	fprintf(stderr, "\n");

	fprintf(stderr, "%sMemory usage:\n", my_tab);
	fprintf(stderr, "%s%s%-24s" "%s\n", my_tab, my_tab,
		"Genomemap:",
		comma_integer(count_get_count(&mem_genomemap)));

	if (Xflag) {
	  print_insert_histogram();
	}
}

static void
usage(char * progname, bool full_usage){
  char *slash;
  int sn;

  if (n_seeds == 0)
    load_default_seeds(0);

  slash = strrchr(progname, '/');
  if (slash != NULL)
    progname = slash + 1;

  fprintf(stderr, 
	  "usage: %s [options/parameters] { <r> | -1 <r1> -2 <r2> } <g1> <g2>...\n", progname);
  fprintf(stderr,
	  "   <r>                  Reads filename, paired or unpaired\n");
  fprintf(stderr,
	  "   <r1>                 Upstream reads filename\n");
  fprintf(stderr,
	  "   <r2>                 Downstream reads filename\n");
  fprintf(stderr,
	  "   <g1> <g2>...         Space seperated list of genome filenames\n");
  fprintf(stderr,
	  "Parameters:\n");
  fprintf(stderr,
	  "   -s/--seeds           Spaced Seed(s)                (default: ");
  for (sn = 0; sn < n_seeds; sn++) {
    if (sn > 0)
      fprintf(stderr,
	  "                                                       ");
    fprintf(stderr, "%s%s\n", seed_to_string(sn), (sn == n_seeds - 1? ")" : ","));
  }
  fprintf(stderr,
	  "   -o/--report          Maximum Hits per Read         (default: %d)\n",
	  DEF_NUM_OUTPUTS);
  fprintf(stderr,
          "      --max-alignments  Max. align. per read  (0=all) (default: %d)\n",
	  DEF_MAX_ALIGNMENTS);
  fprintf(stderr,
	  "   -w/--match-window    Match Window Length           (default: %.02f%%)\n",
	  DEF_WINDOW_LEN);
  fprintf(stderr,
	  "   -n/--cmw-mode        Seed Matches per Window       (default: %d)\n",
	  DEF_NUM_MATCHES);
  if (full_usage) {
  fprintf(stderr,
	  "   -l/--cmw-overlap     Match Window Overlap Length   (default: %.02f%%)\n",
	  DEF_WINDOW_OVERLAP);
  fprintf(stderr,
	  "   -a/--anchor-width    Anchor Width Limiting Full SW (default: %d; disable: -1)\n",
	  DEF_ANCHOR_WIDTH);

  fprintf(stderr, "\n");
  fprintf(stderr,
	  "   -S/--save            Save Genome Proj. in File     (default: no)\n");
  fprintf(stderr,
	  "   -L/--load            Load Genome Proj. from File   (default: no)\n");
  fprintf(stderr,
	  "   -z/--cutoff          Projection List Cut-off Len.  (default: %u)\n",
	  DEF_LIST_CUTOFF);
  }

  fprintf(stderr, "\n");
  fprintf(stderr,
	  "   -m/--match           SW Match Score                (default: %d)\n",
	  DEF_MATCH_VALUE);
  fprintf(stderr,
	  "   -i/--mismatch        SW Mismatch Score             (default: %d)\n",
	  DEF_MISMATCH_VALUE);
  fprintf(stderr,
	  "   -g/--open-r          SW Gap Open Score (Reference) (default: %d)\n",
	  DEF_A_GAP_OPEN);
  fprintf(stderr,
	  "   -q/--open-q          SW Gap Open Score (Query)     (default: %d)\n",
	  DEF_B_GAP_OPEN);
  fprintf(stderr,
	  "   -e/--ext-r           SW Gap Extend Score(Reference)(default: %d)\n",
	  DEF_A_GAP_EXTEND);
  fprintf(stderr,
	  "   -f/--ext-q           SW Gap Extend Score (Query)   (default: %d)\n",
	  DEF_B_GAP_EXTEND);
  if (shrimp_mode == MODE_COLOUR_SPACE) {
  fprintf(stderr,
	  "   -x/--crossover       SW Crossover Score            (default: %d)\n",
	  DEF_XOVER_PENALTY);
  }
  fprintf(stderr,
	  "   -r/--cmw-threshold   Window Generation Threshold   (default: %.02f%%)\n",
	  DEF_WINDOW_GEN_THRESHOLD);
  if (shrimp_mode == MODE_COLOUR_SPACE) {
  fprintf(stderr,
	  "   -v/--vec-threshold   SW Vector Hit Threshold       (default: %.02f%%)\n",
	  DEF_SW_VECT_THRESHOLD);
  }
  fprintf(stderr,
	  "   -h/--full-threshold  SW Full Hit Threshold         (default: %.02f%%)\n",
	  DEF_SW_FULL_THRESHOLD);

  fprintf(stderr, "\n");

  fprintf(stderr,
	  "   -N/--threads         Number of Threads             (default: %d)\n",
	  DEF_NUM_THREADS);
  if (full_usage) {
  fprintf(stderr,
	  "   -K/--thread-chunk    Thread Chunk Size             (default: %d)\n",
	  DEF_CHUNK_SIZE);
  }

  fprintf(stderr, "\n");
  fprintf(stderr,
	  "   -p/--pair-mode       Paired Mode                   (default: %s)\n",
	  pair_mode_string[pair_mode]);
  fprintf(stderr,
	  "   -I/--isize           Min and Max Insert Size       (default: %d,%d)\n",
	  DEF_MIN_INSERT_SIZE, DEF_MAX_INSERT_SIZE);
  fprintf(stderr,
	  "      --longest-read    Maximum read length           (default: %d)\n",
	  DEF_LONGEST_READ_LENGTH);
  fprintf(stderr,
	  "   -1/--upstream        Upstream read pair file\n");
  fprintf(stderr,
	  "   -2/--downstream      Downstream read pair file\n");
  fprintf(stderr,
	  "      --un              Dump unaligned reads to file\n");
  fprintf(stderr,
	  "      --al              Dump aligned reads to file\n");
  fprintf(stderr,
	  "      --read-group      Attach SAM Read Group name\n");
  fprintf(stderr,
          "      --sam-header      Use file as SAM header\n");
  fprintf(stderr,
          "      --trim-front      Trim front of reads by this amount\n");
  fprintf(stderr,
          "      --trim-end        Trim end of reads by this amount\n");
  fprintf(stderr,
          "      --trim-first      Trim only first read in pair\n");
  fprintf(stderr,
          "      --trim-second     Trim only second read in pair\n");
  fprintf(stderr,
          "      --expected-isize  Use this to tie break for high scoring pairs\n");

  fprintf(stderr, "\n");
  fprintf(stderr, "Options:\n");

  fprintf(stderr,
	  "   -U/--ungapped        Perform Ungapped Alignment    (default: disabled)\n");
  fprintf(stderr,
          "      --global          Perform full global alignment (default: disabled)\n");
  if (shrimp_mode == MODE_COLOUR_SPACE) {
  fprintf(stderr,
          "      --bfast           Try to align like bfast       (default: disabled)\n");
  }
  fprintf(stderr,
	  "   -C/--negative        Negative Strand Aln. Only     (default: disabled)\n");
  fprintf(stderr,
	  "   -F/--positive        Positive Strand Aln. Only     (default: disabled)\n");
  fprintf(stderr,
	  "   -P/--pretty          Pretty Print Alignments       (default: disabled)\n");
  fprintf(stderr,
	  "   -E/--sam             Output SAM Format             (default: disabled)\n");
  fprintf(stderr,
  	  "   -Q/--fastq           Reads are in fastq format     (default: disabled)\n");
  if (full_usage) {
  fprintf(stderr,
	  "   -R/--print-reads     Print Reads in Output         (default: disabled)\n");
 // fprintf(stderr,
//	  "    -T    (does nothing since default) Reverse Tie-break on Negative Strand          (default: enabled)\n");
  fprintf(stderr,
	  "   -t/--tiebreak-off    Disable Reverse Tie-break\n");
  fprintf(stderr,
	  "                                  on Negative Strand  (default: enabled)\n");
  fprintf(stderr,
	  "   -X/--isize-hist      Print Insert Size Histogram   (default: disabled)\n");
  fprintf(stderr,
	  "   -Y/--proj-hist       Print Genome Proj. Histogram  (default: disabled)\n");
  fprintf(stderr,
	  "   -Z/--bypass-off      Disable Cache Bypass for SW\n");
  fprintf(stderr,
	  "                                    Vector Calls      (default: enabled)\n");
  fprintf(stderr,
	  "   -H/--spaced-kmers    Hash Spaced Kmers in Genome\n");
  fprintf(stderr,
	  "                                    Projection        (default: disabled)\n");
  fprintf(stderr,
	  "   -D/--thread-stats    Individual Thread Statistics  (default: disabled)\n");
  fprintf(stderr,
	  "   -V/--trim-off        Disable Automatic Genome\n");
  fprintf(stderr,
	  "                                 Index Trimming       (default: enabled)\n");
  }
  fprintf(stderr,
	  "      --sam-unaligned   Unaligned reads in SAM output (default: disabled)\n");
  fprintf(stderr,
	  "      --half-paired     Output half mapped read pairs (default: disabled)\n");
  fprintf(stderr,
	  "      --strata          Print only the best scoring hits\n");
  fprintf(stderr,
	  "   -?/--help            Full List of Parameters and Options\n");

  exit(1);
}

static void
print_settings() {
  static char const my_tab[] = "    ";
  int sn;

  fprintf(stderr, "Settings:\n");
  fprintf(stderr, "%s%-40s%s (%d/%d)\n", my_tab,
	  (n_seeds == 1) ? "Spaced Seed (weight/span)" : "Spaced Seeds (weight/span)",
	  seed_to_string(0), seed[0].weight, seed[0].span);
  for (sn = 1; sn < n_seeds; sn++) {
    fprintf(stderr, "%s%-40s%s (%d/%d)\n", my_tab, "",
	    seed_to_string(sn), seed[sn].weight, seed[sn].span);
  }

  fprintf(stderr, "%s%-40s%d\n", my_tab, "Number of Outputs per Read:", num_outputs);
  fprintf(stderr, "%s%-40s%d\n", my_tab, "Window Generation Mode:", num_matches);

  if (IS_ABSOLUTE(window_len)) {
    fprintf(stderr, "%s%-40s%u\n", my_tab, "Window Length:", (uint)-window_len);
  } else {
    fprintf(stderr, "%s%-40s%.02f%%\n", my_tab, "Window Length:", window_len);
  }

  if (IS_ABSOLUTE(window_overlap)) {
    fprintf(stderr, "%s%-40s%u\n", my_tab, "Window Overlap Length:", (uint)-window_overlap);
  } else {
    fprintf(stderr, "%s%-40s%.02f%%\n", my_tab, "Window Overlap Length:", window_overlap);
  }

  fprintf(stderr, "\n");
  fprintf(stderr, "%s%-40s%d\n", my_tab, "SW Match Score:", match_score);
  fprintf(stderr, "%s%-40s%d\n", my_tab, "SW Mismatch Score:", mismatch_score);
  fprintf(stderr, "%s%-40s%d\n", my_tab, "SW Gap Open Score (Ref):", a_gap_open_score);
  fprintf(stderr, "%s%-40s%d\n", my_tab, "SW Gap Open Score (Qry):", b_gap_open_score);
  fprintf(stderr, "%s%-40s%d\n", my_tab, "SW Gap Extend Score (Ref):", a_gap_extend_score);
  fprintf(stderr, "%s%-40s%d\n", my_tab, "SW Gap Extend Score (Qry):", b_gap_extend_score);
  if (shrimp_mode == MODE_COLOUR_SPACE) {
    fprintf(stderr, "%s%-40s%d\n", my_tab, "SW Crossover Score:", crossover_score);
  }

  fprintf(stderr, "\n");

  if (IS_ABSOLUTE(window_gen_threshold)) {
    fprintf(stderr, "%s%-40s%u\n", my_tab,
	    "Window Generation Threshold:", (uint)-window_gen_threshold);
  } else {
    fprintf(stderr, "%s%-40s%.02f%%\n", my_tab,
	    "Window Generation Threshold:", window_gen_threshold);
  }
  if (shrimp_mode == MODE_COLOUR_SPACE) {
    if (IS_ABSOLUTE(sw_vect_threshold)) {
      fprintf(stderr, "%s%-40s%u\n", my_tab, "SW Vector Hit Threshold:", (uint)-sw_vect_threshold);
    } else {
      fprintf(stderr, "%s%-40s%.02f%%\n", my_tab, "SW Vector Hit Threshold:", sw_vect_threshold);
    }
  }
  if (IS_ABSOLUTE(sw_full_threshold)) {
    fprintf(stderr, "%s%-40s%u\n", my_tab,
	    shrimp_mode == MODE_COLOUR_SPACE? "SW Full Hit Threshold:" : "SW Hit Threshold",
	    (uint)-sw_full_threshold);
  } else {
    fprintf(stderr, "%s%-40s%.02f%%\n", my_tab,
	    shrimp_mode == MODE_COLOUR_SPACE? "SW Full Hit Threshold:" : "SW Hit Threshold",
	    sw_full_threshold);
  }

  fprintf(stderr, "\n");

  fprintf(stderr, "%s%-40s%s\n", my_tab, "Paired mode:", pair_mode_string[pair_mode]);
  if (pair_mode != PAIR_NONE) {
    fprintf(stderr, "%s%-40smin:%d max:%d\n", my_tab, "Insert sizes:", min_insert_size, max_insert_size);
    if (Xflag) {
      fprintf(stderr, "%s%-40s%d\n", my_tab, "Bucket size:", insert_histogram_bucket_size);
    }
  }

  fprintf(stderr, "\n");

  fprintf(stderr, "%s%-40s%s\n", my_tab, "Gapless mode:", gapless_sw? "yes" : "no");
  fprintf(stderr, "%s%-40s%d\n", my_tab, "Number of threads:", num_threads);
  fprintf(stderr, "%s%-40s%d\n", my_tab, "Thread chunk size:", chunk_size);
  fprintf(stderr, "%s%-40s%s\n", my_tab, "Hash Filter Calls:", hash_filter_calls? "yes" : "no");
  fprintf(stderr, "%s%-40s%d%s\n", my_tab, "Anchor Width:", anchor_width,
	  anchor_width == -1? " (disabled)" : "");
  if (list_cutoff < DEF_LIST_CUTOFF) {
  fprintf(stderr, "%s%-40s%u\n", my_tab, "Index List Cutoff Length:", list_cutoff);
  }

}

static int
set_mode_from_string(char const * s) {
  if (!strcmp(s, "mirna")) {
    mode_mirna = true;

    //load_default_mirna_seeds();

    Hflag = true;
    gapless_sw = true;
    anchor_width = 0;
    a_gap_open_score = -255;
    b_gap_open_score = -255;
    hash_filter_calls = false;
    num_matches = 1;
    window_len = 100.0;

    return 1;
  } else {
    return 0;
  }
}


int main(int argc, char **argv){
	char **genome_files = NULL;
	int ngenome_files = 0;

	char *progname = argv[0];
	char const * optstr = NULL;
	char *c;
	int ch;

	bool a_gap_open_set, b_gap_open_set;
	bool a_gap_extend_set, b_gap_extend_set;
	bool num_matches_set = false;

	shrimp_args.argc=argc;
	shrimp_args.argv=argv;
	set_mode_from_argv(argv, &shrimp_mode);

	a_gap_open_set = b_gap_open_set = a_gap_extend_set = b_gap_extend_set = false;

	fprintf(stderr, "--------------------------------------------------"
			"------------------------------\n");
	fprintf(stderr, "gmapper: %s.\nSHRiMP %s\n[%s; CXXFLAGS=\"%s\"]\n",
		get_mode_string(shrimp_mode),
		SHRIMP_VERSION_STRING,
		get_compiler(), CXXFLAGS);
	fprintf(stderr, "--------------------------------------------------"
			"------------------------------\n");

	struct option getopt_long_options[standard_entries+MAX(colour_entries,letter_entries)];
	memcpy(getopt_long_options,standard_options,sizeof(standard_options));
	//TODO -t -9 -d -Z -D -Y
	switch(shrimp_mode){
	case MODE_COLOUR_SPACE:
		optstr = "?1:2:s:n:w:l:o:p:m:i:g:q:e:f:h:r:a:z:DCEFHI:K:L:M:N:PRS:TtUVXYZQx:v:";
		memcpy(getopt_long_options+standard_entries,colour_space_options,sizeof(colour_space_options));
		break;
	case MODE_LETTER_SPACE:
		optstr = "?1:2:s:n:w:l:o:p:m:i:g:q:e:f:h:r:a:z:DCEFHI:K:L:M:N:PRS:TtUVXYZQ";
		memcpy(getopt_long_options+standard_entries,letter_space_options,sizeof(letter_space_options));
		break;
	case MODE_HELICOS_SPACE:
		fprintf(stderr,"Helicose currently unsupported\n");
		exit(1);
		break;
	}

	

	while ((ch = getopt_long(argc,argv,optstr,getopt_long_options,NULL)) != -1){
		switch (ch) {
		case 9:
			strata_flag = true;
			break;
		case 10:
			unaligned_reads_file=fopen(optarg,"w");
			if (unaligned_reads_file==NULL) {
				fprintf(stderr,"error: cannot open file \"%s\" for writting\n",optarg);	
			}
			break;
		case 11:
			aligned_reads_file=fopen(optarg,"w");
			if (aligned_reads_file==NULL) {
				fprintf(stderr,"error: cannot open file \"%s\" for writting\n",optarg);	
			}
			break;
		case 12:
			sam_unaligned=true;
			break;
		case 13:
			longest_read_len=atoi(optarg);
			if (longest_read_len<200) {
				fprintf(stderr,"error: longest read length must be at least 200\n");
				exit(1);
			}
			break;
		case 14:
			max_alignments=atoi(optarg);
			break;
		case 15:
			Gflag = true;
			break;
		case 16:
			Bflag = true;
			Gflag = true;
			break;
		case 17:
			sam_read_group_name=strdup(optarg);
			sam_sample_name = strchr(sam_read_group_name,',');
			if (sam_sample_name==NULL) {
				fprintf(stderr,"error: sam read group needs to be two values, delimited by commas.\n");
				fprintf(stderr,"       the first value is unique read group identifier\n");
				fprintf(stderr,"       the second is the sample (use pool name where a pool is being sequence\n");
				exit(1);	
			}
			sam_sample_name[0]='\0';
			sam_sample_name++;
			break;
		case 18:
			sam_header_filename=optarg;
			break;	
		case 19:
			sam_half_paired=true;
			break;
		case 20:
			sam_r2=true;
			break;
		case '1':
			left_reads_filename = optarg;
			break;
		case '2':
			right_reads_filename = optarg;
			break;
		case 's':
			if (strchr(optarg, ',') == NULL) { // allow comma-separated seeds
				if (optarg[0] == 'w') {
					int weight = (int)atoi(&optarg[1]);
					if (!load_default_seeds(weight)) {
						fprintf(stderr, "error: invalid spaced seed weight (%d)\n", weight);
						exit(1);
					}
				} else {
					if (!add_spaced_seed(optarg)) {
						fprintf(stderr, "error: invalid spaced seed \"%s\"\n", optarg);
						exit (1);
					}
				}
			} else {
				c = strtok(optarg, ",");
				do {
                                	if (c[0] == 'w') {
	                                        int weight = (int)atoi(&c[1]);
        	                                if (!load_default_seeds(weight)) {
                	                                fprintf(stderr, "error: invalid spaced seed weight (%d)\n", weight);
                        	                        exit(1);
                                	        }
	                                } else {
        	                                if (!add_spaced_seed(c)) {
                	                                fprintf(stderr, "error: invalid spaced seed \"%s\"\n", c);
                        	                        exit (1);
                                	        }
                                	}
					c = strtok(NULL, ",");
				} while (c != NULL);
			}
			break;
		case 'n':
			num_matches = atoi(optarg);
			num_matches_set = true;
			break;
		case 'w':
			window_len = atof(optarg);
			if (window_len <= 0.0) {
				fprintf(stderr, "error: invalid window "
						"length\n");
				exit(1);
			}
			if (strcspn(optarg, "%.") == strlen(optarg))
				window_len = -window_len;		//absol.
			break;
		case 'o':
			num_outputs = atoi(optarg);
			num_tmp_outputs = 30 + num_outputs;
			break;
		case 'm':
			match_score = atoi(optarg);
			break;
		case 'i':
			mismatch_score = atoi(optarg);
			break;
		case 'g':
			a_gap_open_score = atoi(optarg);
			a_gap_open_set = true;
			break;
		case 'q':
			b_gap_open_score = atoi(optarg);
			b_gap_open_set = true;
			break;
		case 'e':
			a_gap_extend_score = atoi(optarg);
			a_gap_extend_set = true;
			break;
		case 'f':
			b_gap_extend_score = atoi(optarg);
			b_gap_extend_set = true;
			break;
		case 'x':
			assert(shrimp_mode == MODE_COLOUR_SPACE);
			crossover_score = atoi(optarg);
			break;
		case 'h':
			sw_full_threshold = atof(optarg);
			if (sw_full_threshold < 0.0) {
				fprintf(stderr, "error: invalid s-w full "
						"hit threshold\n");
				exit(1);
			}
			if (strcspn(optarg, "%.") == strlen(optarg))
				sw_full_threshold = -sw_full_threshold;	//absol.
			break;
		case 'v':
		  assert(shrimp_mode == MODE_COLOUR_SPACE);
			sw_vect_threshold = atof(optarg);
			if (sw_vect_threshold < 0.0) {
				fprintf(stderr, "error: invalid s-w vector "
						"hit threshold\n");
				exit(1);
			}
			if (strcspn(optarg, "%.") == strlen(optarg))
				sw_vect_threshold = -sw_vect_threshold;	//absol.
			break;
		case 'r':
		  window_gen_threshold = atof(optarg);
		  if (window_gen_threshold < 0.0) {
		    fprintf(stderr, "error: invalid window generation threshold [%s]\n", optarg);
		    exit(1);
		  }
		  if (strcspn(optarg, "%.") == strlen(optarg))
		    window_gen_threshold = -window_gen_threshold;	//absol.
		  break;
		case 'C':
			if (Fflag) {
				fprintf(stderr, "error: -C and -F are mutually "
						"exclusive\n");
				exit(1);
			}
			Cflag = true;
			break;
		case 'F':
			if (Cflag) {
				fprintf(stderr, "error: -C and -F are mutually "
						"exclusive\n");
				exit(1);
			}
			Fflag = true;
			break;
		case 'H':
			Hflag = true;
			break;
		case 'P':
			Pflag = true;
			break;
		case 'R':
			Rflag = true;
			break;
		case 't':
			Tflag= false;
			break;
		case 'T':
			Tflag = true;
			break;
			/*
			 * New options/parameters since SHRiMP 1.2.1
			 */
		case 'a':
			anchor_width = atoi(optarg);
			if (anchor_width < -1 || anchor_width >= 100) {
				fprintf(stderr, "error: anchor_width requested is invalid (%s)\n",
						optarg);
				exit(1);
			}
			break;
		case 'X':
			Xflag = true;
			break;
		case 'Y':
			Yflag = true;
			break;
		case 'l':
			window_overlap = atof(optarg);
			if (window_overlap <= 0.0) {
				fprintf(stderr, "error: invalid window overlap\n");
				exit(1);
			}
			if (strcspn(optarg, "%.") == strlen(optarg))
				window_overlap = -window_overlap;		//absol.
			break;
		case 'N':
			num_threads = atoi(optarg);
			break;
		case 'K':
			chunk_size = atoi(optarg);
			break;
		case 'S':
			save_file = optarg;
			break;
		case 'L':
			load_file = optarg;
			break;
		case 'D':
			Dflag = true;
			break;
		case '?':
			usage(progname, true);
			break;
		case 'Z':
		  hash_filter_calls = false;
		  break;
		case 'U':
		  gapless_sw = true;
		  anchor_width = 0;
		  a_gap_open_score = -255;
		  b_gap_open_score = -255;
		  hash_filter_calls = false;
		  break;
		case 'p':
		  if (!strcmp(optarg, "none")) {
		    pair_mode = PAIR_NONE;
		  } else if (!strcmp(optarg, "opp-in")) {
		    pair_mode = PAIR_OPP_IN;
		  } else if (!strcmp(optarg, "opp-out")) {
		    pair_mode = PAIR_OPP_OUT;
		  } else if (!strcmp(optarg, "col-fw")) {
		    pair_mode = PAIR_COL_FW;
		  } else if (!strcmp(optarg, "col-bw")) {
		    pair_mode = PAIR_COL_BW;
		  } else {
		    fprintf(stderr, "error: unrecognized pair mode (%s)\n", optarg);
		    exit(1);
		  }
		  break;
		case 'I':
		  c = strtok(strdup(optarg), ",");
		  if (c == NULL) {
		    fprintf(stderr, "error: format for insert sizes is \"-I 200,1000\"\n");
		    exit(1);
		  }
		  min_insert_size = atoi(c);
		  c = strtok(NULL, ",");
		  if (c == NULL) {
		    fprintf(stderr, "error: format for insert sizes is \"-I 200,1000\"\n");
		    exit(1);
		  }
		  max_insert_size = atoi(c);
		  if (min_insert_size > max_insert_size) {
		    fprintf(stderr, "error: invalid insert sizes (min:%d,max:%d)\n",
			    min_insert_size, max_insert_size);
		    exit(1);
		  }
		  break;
		case 'E':
			Eflag = true;
			break;
		case 'z':
		  list_cutoff = atoi(optarg);
		  if (list_cutoff == 0) {
		    fprintf(stderr, "error: invalid list cutoff (%s)\n", optarg);
		    exit(1);
		  }
		  break;
		case 'V':
			Vflag = false;
			break;
		case 'Q':
			Qflag = true;
			break;
                case 'M':
                        c = strtok(optarg, ",");
                        do {
                          if (!set_mode_from_string(c)) {
                            fprintf(stderr, "error: unrecognized mode (%s)\n", c);
                            exit(1);
                          } 
                          c = strtok(NULL, ",");
                        } while (c != NULL);
                        break;
		case 21:
			trim=true;
			trim_front=atoi(optarg);
			if (shrimp_mode == MODE_COLOUR_SPACE) {
				fprintf(stderr,"--trim-front cannot be used in colour space mode!\n");
				exit(1);
			}
			if (trim_front<0) {
				fprintf(stderr,"--trim-front value must be positive\n");
				exit(1);
			}
			break;
		case 22:
			trim=true;
			trim_end=atoi(optarg);
			if (trim_end<0) {
				fprintf(stderr,"--trim-end value must be positive\n");
				exit(1);
			}
			break;
		case 23:
			trim=true;
			trim_first=true;
			trim_second=false;
			break;
		case 24:
			trim=true;
			trim_second=true;
			trim_first=false;
			break;
		case 25:
			expected_isize=atoi(optarg);
			if (expected_isize<0) {
				fprintf(stderr,"Expected insert size needs to be positive!\n");
				exit(1);
			}
			break;
		default:
			usage(progname, false);
		}
	}

	argc -= optind;
	argv += optind;

	if (Gflag && gapless_sw) {
		fprintf(stderr,"error: cannot use global (or bfast) and ungapped mode at the same time!\n");
		usage(progname,false);
	}
	if (sam_unaligned && !Eflag) {
		fprintf(stderr,"error: when using flag --sam-unaligned must also use -E/--sam\n");
		usage(progname,false);
	}
	if (right_reads_filename != NULL || left_reads_filename !=NULL) {
		if (right_reads_filename == NULL || left_reads_filename == NULL ){
			fprintf(stderr,"error: when using \"%s\" must also specify \"%s\"\n",
				(left_reads_filename != NULL) ? "-1" : "-2",
				(left_reads_filename != NULL) ? "-2" : "-1");
			usage(progname,false);
		}
		single_reads_file=false;
		if (strcmp(right_reads_filename,"-")==0 && strcmp(left_reads_filename,"-")==0) {
			fprintf(stderr,"error: both -1 and -2 arguments cannot be stdin (\"-\")\n");
			usage(progname,false);
		}
	}
	if (pair_mode != PAIR_NONE && !num_matches_set) {
	  num_matches = 4;
	}
	if (pair_mode == PAIR_NONE && (!trim_first || !trim_second)) {
		fprintf(stderr,"error: cannot use --trim-first or --trim-second in unpaired mode\n");
		usage(progname,false);
	} 

	if (Xflag) {
	  if (pair_mode == PAIR_NONE) {
	    fprintf(stderr, "warning: insert histogram not available in unpaired mode; ignoring\n");
	    Xflag = false;
	  } else {
	    insert_histogram_bucket_size = ceil_div(max_insert_size - min_insert_size + 1, 100);
	  }
	}

	if(load_file != NULL && n_seeds != 0){
	  fprintf(stderr,"error: cannot specify seeds when loading genome map\n");
	  usage(progname,false);
	}

	if (n_seeds == 0 && load_file == NULL) {
          if (mode_mirna)
            load_default_mirna_seeds();
          else
            load_default_seeds(0);
	}

	if (Hflag){
	  init_seed_hash_mask();
	}

	if (save_file != NULL && load_file != NULL && list_cutoff == DEF_LIST_CUTOFF){
	  fprintf(stderr,"error: -L and -S allowed together only if list_cutoff is specified\n");
	  exit(1);
	}

	if (load_file != NULL && save_file != NULL)
	  { // args: none
	    if (argc != 0) {
	      fprintf(stderr, "error: when using both -L and -S, no extra files can be given\n");
	      usage(progname, false);
	    }
	  } 
	else if (load_file != NULL)
	  { // args: reads file
	    if (argc == 0 && single_reads_file) {
	      fprintf(stderr,"error: read_file not specified\n");
	      usage(progname, false);
	    } else if (argc == 1) {
	      if (single_reads_file) {
	      	reads_filename    = argv[0];
	      } else {
		fprintf(stderr,"error: cannot specify a reads file when using -L, -1 and -2\n");
		usage(progname,false);
	      }
	    }
	  }
	else if (save_file != NULL)
	  { // args: genome file(s)
	    if (argc == 0){
	      fprintf(stderr, "error: genome_file(s) not specified\n");
	      usage(progname,false);
	    }
	    genome_files  = &argv[0];
	    ngenome_files = argc;
	  }
	else if (single_reads_file)
	  { // args: reads file, genome file(s)
	    if (argc < 2) {
	      fprintf(stderr, "error: %sgenome_file(s) not specified\n",
		      (argc == 0) ? "reads_file, " : "");
	      usage(progname, false);
	    }
	    reads_filename    = argv[0];
	    genome_files  = &argv[1];
	    ngenome_files = argc - 1;
	  }
	else 
	  {
	   if( argc < 1) {
	      fprintf(stderr, "error: genome_file(s) not specified\n");
	      usage(progname, false);
	   }
	    genome_files  = &argv[0];
	    ngenome_files = argc;
	  }

	if (!Cflag && !Fflag) {
	  Cflag = Fflag = true;
	}

	if (pair_mode != PAIR_NONE && (!Cflag || !Fflag)) {
	  fprintf(stderr, "warning: in paired mode, both strands must be inspected; ignoring -C and -F\n");
	  Cflag = Fflag = true;
	}
	if (pair_mode == PAIR_NONE && sam_half_paired) {
	  fprintf(stderr, "error: cannot use option half-paired in non-paired mode!\n");
	  exit(1);
	}
	if (pair_mode == PAIR_NONE && sam_r2) {
	  fprintf(stderr, "error: cannot use option sam-r2 in non-paired mode!\n");
	  exit(1);
	}
	

	if (shrimp_mode == MODE_LETTER_SPACE) {
	  sw_vect_threshold = sw_full_threshold;
	}

	if (Eflag && Pflag) {
	  fprintf(stderr,"-E and -P are incompatable\n");
	  exit(1);
	}

	if (Eflag && Rflag) {
	  fprintf(stderr,"-E and -R are incompatable\n");
	  exit(1);
	}

	if (!valid_spaced_seeds()) {
	  fprintf(stderr, "error: invalid spaced seed\n");
	  if (!Hflag)
	    fprintf(stderr, "       for longer seeds, try using the -H flag\n");
	  exit(1);
	}

	if (!IS_ABSOLUTE(window_len) && window_len < 100.0) {
	  fprintf(stderr, "error: window length < 100%% of read length\n");
	  exit(1);
	}

	if (!IS_ABSOLUTE(window_overlap) && window_overlap > 100.0) {
	  fprintf(stderr, "warning: window overlap length > 100%% of window_length; resetting to 100%%\n");
	  window_overlap = 100.0;
	}

	if (num_matches < 1) {
	  fprintf(stderr, "error: invalid number of matches\n");
	  exit(1);
	}

	if (num_outputs < 1) {
	  fprintf(stderr, "error: invalid maximum hits per read\n");
	  exit(1);
	}

	if (a_gap_open_score > 0 || b_gap_open_score > 0) {
	  fprintf(stderr, "error: invalid gap open penalty\n");
	  exit(1);
	}

	if (a_gap_extend_score > 0 || b_gap_extend_score > 0) {
	  fprintf(stderr, "error: invalid gap extend penalty\n");
	  exit(1);
	}

	if (!IS_ABSOLUTE(sw_full_threshold) && sw_full_threshold > 100.0) {
	  fprintf(stderr, "error: invalid s-w full hit threshold\n");
	  exit(1);
	}

	if (shrimp_mode == MODE_COLOUR_SPACE && !IS_ABSOLUTE(sw_vect_threshold)
	    && sw_vect_threshold > 100.0) {
	  fprintf(stderr, "error: invalid s-w vector threshold\n");
	  exit(1);
	}

	if (!IS_ABSOLUTE(window_gen_threshold) && window_gen_threshold > 100.0) {
	  fprintf(stderr, "error: invalid window generation threshold\n");
	  exit(1);
	}

	if ((IS_ABSOLUTE(window_gen_threshold) && IS_ABSOLUTE(sw_full_threshold)
	     && -window_gen_threshold > -sw_full_threshold)
	    ||
	    (!IS_ABSOLUTE(window_gen_threshold) && !IS_ABSOLUTE(sw_full_threshold)
	     && window_gen_threshold > sw_full_threshold)) {
	  fprintf(stderr, "warning: window generation threshold is larger than sw threshold\n");
	}

	if ((a_gap_open_set && !b_gap_open_set) || (a_gap_extend_set && !b_gap_extend_set)) {
	  fputc('\n', stderr);
	}

	if (a_gap_open_set && !b_gap_open_set) {
	  fprintf(stderr, "Notice: Gap open penalty set for reference but not query; assuming symmetry.\n");
	  b_gap_open_score = a_gap_open_score;
	}
	if (a_gap_extend_set && !b_gap_extend_set) {
	  fprintf(stderr, "Notice: Gap extend penalty set for reference but not query; assuming symmetry.\n");
	  b_gap_extend_score = a_gap_extend_score;
	}

	if ((a_gap_open_set && !b_gap_open_set) || (a_gap_extend_set && !b_gap_extend_set)) {
	  fputc('\n', stderr);
	}

	if(load_file == NULL){
	  print_settings();
	}

	uint64_t before;
	before = gettimeinusecs();
	if (load_file != NULL){
		if (strchr(load_file, ',') == NULL) {
			//use prefix
			int buf_size = strlen(load_file) + 20;
			char * genome_name = (char *)xmalloc(sizeof(char)*buf_size);
			strncpy(genome_name,load_file,buf_size);
			strncat(genome_name,".genome",buf_size);
			fprintf(stderr,"Loading genome from %s\n",genome_name);
			if (!load_genome_map(genome_name)){
				fprintf(stderr, "error: loading from genome file \"%s\"\n", genome_name);
				exit (1);
			}
			free(genome_name);
			int seed_count = 0;
			char * seed_name = (char *)xmalloc(sizeof(char)*buf_size);
			char * buff = (char *)xmalloc(sizeof(char)*buf_size);
			strncpy(seed_name,load_file,buf_size);
			strncat(seed_name,".seed.",buf_size);
			sprintf(buff,"%d",seed_count);
			strncat(seed_name,buff,buf_size);
			FILE *f = fopen(seed_name,"r");
			while(f != NULL){
				fclose(f);
				fprintf(stderr,"Loading seed from %s\n",seed_name);
				if (!load_genome_map_seed(seed_name)) {
					fprintf(stderr, "error: loading from map file \"%s\"\n", seed_name);
					exit (1);
				}
				seed_count++;
				strncpy(seed_name,load_file,buf_size);
				strncat(seed_name,".seed.",buf_size);
				sprintf(buff,"%d",seed_count);
				strncat(seed_name,buff,buf_size);
				f = fopen(seed_name,"r");
			}
			free(seed_name);
			free(buff);

		} else {
			c = strtok(load_file, ",");
			fprintf(stderr,"Loading genome from %s\n",c);
			if (!load_genome_map(c)){
				fprintf(stderr, "error: loading from genome file \"%s\"\n", c);
				exit (1);
			}
			c = strtok(NULL, ",");
			do {
				fprintf(stderr,"Loading seed from %s\n",c);
				if (!load_genome_map_seed(c)) {
					fprintf(stderr, "error: loading from map file \"%s\"\n", c);
					exit (1);
				}
				c = strtok(NULL, ",");
			} while (c != NULL);
		}

		if (Hflag) {
		  init_seed_hash_mask();
		}

		print_settings();
	} else {
		if (!load_genome(genome_files,ngenome_files)){
			exit(1);
		}
	}

	map_usecs += (gettimeinusecs() - before);

	//
	// Automatic genome index trimming
	//
	if (Vflag && save_file == NULL && list_cutoff == DEF_LIST_CUTOFF) {
	  // this will be a mapping job; enable automatic trimming
	  int i, sn;
	  long long unsigned int total_genome_len = 0;
	  int max_seed_weight = 0;

	  for (i = 0; i < num_contigs; i++) {
	    total_genome_len += (long long unsigned int)genome_len[i];
	  }

	  if (Hflag) {
	    max_seed_weight = HASH_TABLE_POWER;
	  } else {
	    for (sn = 0; sn < n_seeds; sn++) {
	      if (seed[sn].weight > max_seed_weight) {
		max_seed_weight = seed[sn].weight;
	      }
	    }
	  }

	  // cutoff := max (1000, 100*(total_genome_len/4^max_seed_weight))
	  list_cutoff = 1000;
	  if ((uint32_t)((100llu * total_genome_len)/power(4, max_seed_weight)) > list_cutoff) {
	    list_cutoff = (uint32_t)((100llu * total_genome_len)/power(4, max_seed_weight));
	  }
	  fprintf(stderr, "Automatically trimming genome index lists longer than: %u\n", list_cutoff);
	}

	if (Yflag)
	  print_genomemap_stats();

	if (save_file != NULL) {
	  if (list_cutoff != DEF_LIST_CUTOFF) {
	    fprintf(stderr, "Trimming genome map lists longer than %u\n", list_cutoff);
	    trim_genome();
	  }

	  fprintf(stderr,"Saving genome map to %s\n",save_file);
	  if(save_genome_map(save_file)){
	    exit(0);
	  }
	  exit(1);
	}

	// set up new options structure
	// THIS SHOULD EVENTUALLY BE MERGED INTO OPTION READING
	assert(n_unpaired_mapping_options[0] == 0);
	assert(n_paired_mapping_options == 0);
	if (pair_mode == PAIR_NONE)
	  {
	    n_unpaired_mapping_options[0]++;
	    unpaired_mapping_options[0] = (struct read_mapping_options_t *)realloc((void *)unpaired_mapping_options[0],
										   n_unpaired_mapping_options[0] * sizeof(unpaired_mapping_options[0][0]));

	    unpaired_mapping_options[0][0].anchor_list.recompute = true;
	    unpaired_mapping_options[0][0].anchor_list.use_region_counts = false;
	    unpaired_mapping_options[0][0].anchor_list.use_pairing = false;
	    unpaired_mapping_options[0][0].anchor_list.min_count[0] = (num_matches == 2? 2 : 1);
	    unpaired_mapping_options[0][0].anchor_list.max_count[0] = 0;
	    unpaired_mapping_options[0][0].hit_list.recompute = true;
	    unpaired_mapping_options[0][0].hit_list.gapless = gapless_sw;
	    unpaired_mapping_options[0][0].hit_list.match_mode = num_matches;
	    unpaired_mapping_options[0][0].hit_list.threshold = window_gen_threshold;
	    unpaired_mapping_options[0][0].pass1.recompute =  true;
	    unpaired_mapping_options[0][0].pass1.only_paired = false;
	    unpaired_mapping_options[0][0].pass1.gapless = gapless_sw;
	    unpaired_mapping_options[0][0].pass1.num_outputs = num_tmp_outputs;
	    unpaired_mapping_options[0][0].pass1.threshold = sw_vect_threshold;
	    unpaired_mapping_options[0][0].pass1.window_overlap = window_overlap;
	    unpaired_mapping_options[0][0].pass2.recompute = true;
	    unpaired_mapping_options[0][0].pass2.strata = strata_flag;
	    unpaired_mapping_options[0][0].pass2.num_outputs = num_outputs;
	    unpaired_mapping_options[0][0].pass2.threshold = sw_full_threshold;
	    unpaired_mapping_options[0][0].stop_count = 0;
	  }
	else
	  {
	    n_paired_mapping_options++;
	    paired_mapping_options = (struct readpair_mapping_options_t *)realloc((void *)paired_mapping_options,
										  n_paired_mapping_options * sizeof(paired_mapping_options[0]));

	    paired_mapping_options[0].pairing.pair_mode = pair_mode;
	    paired_mapping_options[0].pairing.pair_up_hits = true;
	    paired_mapping_options[0].pairing.min_insert_size = min_insert_size;
	    paired_mapping_options[0].pairing.max_insert_size = max_insert_size;
	    paired_mapping_options[0].pairing.pass1_num_outputs = num_tmp_outputs;
	    paired_mapping_options[0].pairing.pass2_num_outputs = num_outputs;
	    paired_mapping_options[0].pairing.pass1_threshold = sw_vect_threshold;
	    paired_mapping_options[0].pairing.pass2_threshold = sw_full_threshold;

	    paired_mapping_options[0].read[0].anchor_list.recompute = true;
	    paired_mapping_options[0].read[0].anchor_list.use_region_counts = false;
	    paired_mapping_options[0].read[0].anchor_list.use_pairing = true;
	    paired_mapping_options[0].read[0].anchor_list.min_count[0] = (num_matches == 4? 2 : 1);
	    paired_mapping_options[0].read[0].anchor_list.min_count[1] = (num_matches == 4? 2 : 1);
	    paired_mapping_options[0].read[0].anchor_list.max_count[0] = 0;
	    paired_mapping_options[0].read[0].anchor_list.max_count[1] = 0;
	    paired_mapping_options[0].read[0].hit_list.recompute = true;
	    paired_mapping_options[0].read[0].hit_list.gapless = gapless_sw;
	    paired_mapping_options[0].read[0].hit_list.match_mode = (num_matches == 4? 2 : 1);
	    paired_mapping_options[0].read[0].hit_list.threshold = window_gen_threshold;
	    paired_mapping_options[0].read[0].pass1.recompute = true;
	    paired_mapping_options[0].read[0].pass1.only_paired = true;
	    paired_mapping_options[0].read[0].pass1.gapless = gapless_sw;
	    //paired_mapping_options[0].read[0].pass1.num_outputs = 0;
	    paired_mapping_options[0].read[0].pass1.threshold = sw_vect_threshold;
	    paired_mapping_options[0].read[0].pass1.window_overlap = window_overlap;
	    paired_mapping_options[0].read[0].pass2.recompute = true;
	    paired_mapping_options[0].read[0].pass2.strata = strata_flag;
	    //paired_mapping_options[0].read[0].pass2.num_outputs = 0;
	    paired_mapping_options[0].read[0].pass2.threshold = sw_full_threshold / 3;
	    paired_mapping_options[0].read[1] = paired_mapping_options[0].read[0];

	    if (!sam_half_paired)
	      {
		paired_mapping_options[0].pairing.stop_count = 0;
	      }
	    else // half_paired
	      {
		n_unpaired_mapping_options[0]++;
		n_unpaired_mapping_options[1]++;
		unpaired_mapping_options[0] = (struct read_mapping_options_t *)realloc((void *)unpaired_mapping_options[0],
										       n_unpaired_mapping_options[0] * sizeof(unpaired_mapping_options[0][0]));
		unpaired_mapping_options[1] = (struct read_mapping_options_t *)realloc((void *)unpaired_mapping_options[1],
										       n_unpaired_mapping_options[1] * sizeof(unpaired_mapping_options[1][0]));

		paired_mapping_options[0].pairing.stop_count = 1;
		paired_mapping_options[0].pairing.stop_threshold = paired_mapping_options[0].pairing.pass2_threshold;

		unpaired_mapping_options[0][0].anchor_list.recompute = false;
		unpaired_mapping_options[0][0].hit_list.recompute = false;
		unpaired_mapping_options[0][0].pass1.recompute = true; /// ??????????
		unpaired_mapping_options[0][0].pass1.gapless = gapless_sw;
		unpaired_mapping_options[0][0].pass1.only_paired = false;
		unpaired_mapping_options[0][0].pass1.num_outputs = num_tmp_outputs;
		unpaired_mapping_options[0][0].pass1.threshold = sw_vect_threshold;
		unpaired_mapping_options[0][0].pass1.window_overlap = window_overlap;
		unpaired_mapping_options[0][0].pass2.recompute = true;
		unpaired_mapping_options[0][0].pass2.strata = strata_flag;
		unpaired_mapping_options[0][0].pass2.num_outputs = num_outputs;
		unpaired_mapping_options[0][0].pass2.threshold = sw_full_threshold;
		unpaired_mapping_options[0][0].stop_count = 0;
		unpaired_mapping_options[1][0] = unpaired_mapping_options[0][0];

	      }
	  }

	//TODO setup need max window and max read len
	//int longest_read_len = 2000;
	int max_window_len = (int)abs_or_pct(window_len,longest_read_len);
#pragma omp parallel shared(longest_read_len,max_window_len,a_gap_open_score, a_gap_extend_score, b_gap_open_score, b_gap_extend_score,\
		match_score, mismatch_score,shrimp_mode,crossover_score,anchor_width) num_threads(num_threads)
	{
	  //hash_mark = 0;
	  //window_cache = (struct window_cache_entry *)xcalloc(1048576 * sizeof(window_cache[0]));

		if (f1_setup(max_window_len, longest_read_len,
			     a_gap_open_score, a_gap_extend_score, b_gap_open_score, b_gap_extend_score,
			     match_score, mismatch_score,
			     shrimp_mode == MODE_COLOUR_SPACE, false)) {
			fprintf(stderr, "failed to initialise vector "
					"Smith-Waterman (%s)\n", strerror(errno));
			exit(1);
		}

		int ret;
		if (shrimp_mode == MODE_COLOUR_SPACE) {
			/* XXX - a vs. b gap */
			ret = sw_full_cs_setup(max_window_len, longest_read_len,
					a_gap_open_score, a_gap_extend_score, match_score, mismatch_score,
					crossover_score, false, anchor_width);
		} else {
			ret = sw_full_ls_setup(max_window_len, longest_read_len,
					a_gap_open_score, a_gap_extend_score, b_gap_open_score, b_gap_extend_score,
					match_score, mismatch_score, false, anchor_width);
		}
		if (ret) {
			fprintf(stderr, "failed to initialise scalar "
					"Smith-Waterman (%s)\n", strerror(errno));
			exit(1);
		}
	}


	char * output;
	if (Eflag){
		int i;
		if (sam_header_filename!=NULL) {
			FILE * sam_header_file = fopen(sam_header_filename,"r");
			if (sam_header_file==NULL) {
				perror("Failed to open sam header file ");
				exit(1);
			}
			size_t buffer_size=2046;
			char buffer[buffer_size];
			size_t read; bool ends_in_newline=true;
			while ((read=fread(buffer,1,buffer_size-1,sam_header_file))) {
				buffer[read]='\0';
				fprintf(stdout,"%s",buffer);
				if (buffer[read-1]=='\n') {
					ends_in_newline=true;
				} else {
					ends_in_newline=false;
				}
			}
			if (!ends_in_newline) {
				fprintf(stdout,"\n");
			}
		} else {
			//Print sam header
			fprintf(stdout,"@HD\tVN:%s\tSO:%s\n","1.0","unsorted");

			for(i = 0; i < num_contigs; i++){
				fprintf(stdout,"@SQ\tSN:%s\tLN:%u\n",contig_names[i],genome_len[i]);
			}
		}
		//read group
		if (sam_read_group_name!=NULL) {
			fprintf(stdout,"@RG\tID:%s\tSM:%s\n",sam_read_group_name,sam_sample_name);
		}
		//print command line args used to invoke SHRiMP
		fprintf(stdout,"@PG\tID:%s\tVN:%s\tCL:","gmapper",SHRIMP_VERSION_STRING);
		for (i=0; i<(shrimp_args.argc-1); i++) {
			fprintf(stdout,"%s ",shrimp_args.argv[i]);
		}
		fprintf(stdout,"%s\n",shrimp_args.argv[i]);
	} else {
		output = output_format_line(Rflag);
		puts(output);
		free(output);
	}
	before = gettimeinusecs();
	bool launched = launch_scan_threads();
	if (!launched) {
		fprintf(stderr,"error: a fatal error occured while launching scan thread(s)!\n");
		exit(1);
	}
	total_work_usecs += (gettimeinusecs() - before);
	
	//if (load_file==NULL) {
		free_genome();
	//}
	print_statistics();
#pragma omp parallel shared(longest_read_len,max_window_len,a_gap_open_score, a_gap_extend_score, b_gap_open_score, b_gap_extend_score,\
		match_score, mismatch_score,shrimp_mode,crossover_score,anchor_width) num_threads(num_threads)
	{
		sw_vector_cleanup();
		if (shrimp_mode==MODE_COLOUR_SPACE) {
			sw_full_cs_cleanup();
		}
		sw_full_ls_cleanup();
		f1_free();	
	}
	int i;
	for (i=0; i<num_contigs; i++){
		free(contig_names[i]);
		if (i==0 || load_file==NULL) {
			free(genome_contigs[i]);
			free(genome_contigs_rc[i]);
		}
	}
	if (shrimp_mode==MODE_COLOUR_SPACE) {
		for (i=0; i<num_contigs && (i==0 || load_file==NULL); i++){
			free(genome_cs_contigs[i]);
		}
		free(genome_cs_contigs);
		free(genome_initbp);
	}
	free(genome_len);
	free(genome_contigs_rc);
	free(genome_contigs);
	free(contig_names);
	free(contig_offsets);
	free(seed);
	return 0;
}
