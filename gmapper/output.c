#define _MODULE_OUTPUT

#include <omp.h>

#include "output.h"
#include "../common/output.h"
#include "../common/sw-full-common.h"


/*
	read_start and read_end are 1 based
*/
static cigar_t *
make_cigar(int read_start, int read_end , int read_length, char* qralign,char* dbalign)
{
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


static void
reverse_cigar(cigar_t * cigar)
{
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


static char *
make_cigar_string(cigar_t * cigar)
{
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


static void
free_cigar(cigar_t * cigar)
{
	if (cigar->size>0) {
		assert(cigar->ops!=NULL);
		assert(cigar->lengths!=NULL);
		free(cigar->ops);
		free(cigar->lengths);
	}
	free(cigar);
}


//TODO move this to utils
static void
reverse(char* s, char* t)  // USING TWO ARRAYS FOR THIS IS BAD CODE!!!
{
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
		"\tAS:i:%d\tH0:i:%d\tH1:i:%d\tH2:i:%d\tNM:i:%d\tNH:i:%d\tIH:i:%d\tX0:i:%d",
				   rh->sfrp->score,hits[0],hits[1],hits[2],rh->sfrp->mismatches+rh->sfrp->deletions+rh->sfrp->insertions,satisfying_alignments,stored_alignments, rh->matches);
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


void
new_read_output(struct read_entry * re, struct read_hit * * hits_pass2, int * n_hits_pass2)
{
  int i;
  int hits[] = {0, 0, 0};

  for (i = 0; i < *n_hits_pass2; i++) {
    struct sw_full_results * sfrp = hits_pass2[i]->sfrp;
    int edit_distance = sfrp->mismatches + sfrp->insertions + sfrp->deletions;
    if (0 <= edit_distance && edit_distance < 3) {
      hits[edit_distance]++;
    }
  }

  /* Output sorted list, removing any duplicates. */
  for (i = 0; i < *n_hits_pass2; i++) {
    struct read_hit * rh = hits_pass2[i];
    char * output1 = NULL, * output2 = NULL;
    if (sam_half_paired) {
      int other_hits[] = {0, 0, 0};
      if (re->first_in_pair) {
	hit_output(re, rh, NULL, &output1, NULL, true, hits, *n_hits_pass2);
	hit_output(re->mate_pair, NULL, rh, &output2, NULL, false, other_hits, 0);
      } else {
	hit_output(re->mate_pair, NULL, rh, &output1, NULL, true, other_hits, 0);
	hit_output(re, rh, NULL, &output2, NULL, false, hits, *n_hits_pass2);
      }
    } else {
      hit_output(re, rh, NULL,  &output1, &output2, false, hits, *n_hits_pass2);
    }
    //legacy support for SHRiMP output
    if (!Eflag) {
      if (!Pflag) {
#pragma omp critical (stdout)
	{
	  fprintf(stdout, "%s\n", output1);
	}
      } else {
#pragma omp critical (stdout)
	{
	  fprintf(stdout, "%s\n\n%s\n", output1, output2);
	}
      }
      free(output1);
      if (Pflag || sam_half_paired) {
	free(output2);
      }
    }
  }
}


void
new_readpair_output(struct read_entry * re1, struct read_entry * re2,
		    struct read_hit_pair * hits_pass2, int * n_hits_pass2)
{
  int i;
  int hits1[] = {0, 0, 0};
  int hits2[] = {0, 0, 0};
  
  for (i = 0; i < *n_hits_pass2; i++) {
    struct sw_full_results * sfrp1 = hits_pass2[i].rh[0]->sfrp;
    struct sw_full_results * sfrp2 = hits_pass2[i].rh[1]->sfrp;
    int edit_distance1 = sfrp1->mismatches + sfrp1->insertions + sfrp1->deletions;
    int edit_distance2 = sfrp2->mismatches + sfrp2->insertions + sfrp2->deletions;
    if (0 <= edit_distance1 && edit_distance1 < 3) {
      hits1[edit_distance1]++;
    }
    if (0 <= edit_distance2 && edit_distance2 < 3) {
      hits2[edit_distance2]++;
    }
  }

  /* Output sorted list, removing any duplicates. */
  for (i = 0; i < *n_hits_pass2; i++) {
    struct read_hit * rh1 = hits_pass2[i].rh[0];
    struct read_hit * rh2 = hits_pass2[i].rh[1];
    uint bucket;

    char * output1 = NULL, * output2 = NULL, * output3 = NULL, * output4 = NULL;

    hit_output(re1, rh1,  rh2, &output1, &output2, true, hits1, *n_hits_pass2);
    hit_output(re2, rh2,  rh1, &output3, &output4, false, hits2, *n_hits_pass2);
    if (!Eflag) {
      if (!Pflag) {
#pragma omp critical (stdout)
	{
	  fprintf(stdout, "%s\n%s\n", output1, output3);
	}
      } else {
#pragma omp critical (stdout)
	{
	  fprintf(stdout, "%s\n%s\n%s\n%s\n", output1, output2, output3, output4);
	}
      }
    }

    if (Xflag) {
      if (hits_pass2[i].insert_size < min_insert_size)
	bucket = 0;
      else if (hits_pass2[i].insert_size > max_insert_size)
	bucket = 99;
      else
	bucket = (uint)((hits_pass2[i].insert_size - min_insert_size) / insert_histogram_bucket_size);
      if (bucket >= 100)
	bucket = 99;
#pragma omp atomic
      insert_histogram[bucket]++;
    }

    free(output1);
    free(output3);
    if (Pflag) {
      free(output2);
      free(output4);
    }
  }
}
