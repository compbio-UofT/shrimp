#define _MODULE_OUTPUT

#include <omp.h>

#include "gmapper.h"
#include "output.h"
#include "../common/output.h"
#include "mapping.h"
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
reverse_alignment_edit_string(char * editstr)
{
  int n = strlen(editstr);
  char * res = (char *)malloc((n + 1) * sizeof(char));
  int i = 0;
  while (i < n) {
    if (isdigit(editstr[n - 1 - i])) {
      int j = i + 1;
      while (j < n && isdigit(editstr[n - 1 - j])) j++;
      j--;
      memcpy(&res[i], &editstr[n - 1 - j], (j - i + 1) * sizeof(char));
      i = j + 1;
    } else if (editstr[n - 1 - i] == '-' || editstr[n - 1 - i] == 'x') {
      res[i] = editstr[n - 1 - i];
      i++;
    } else if (editstr[n - 1 - i] == ')') {
      res[i] = '(';
      i++;
    } else if (editstr[n - 1 - i] == '(') {
      res[i] = ')';
      i++;
    } else if (editstr[n - 1 - i] == 'A') {
      res[i] = 'T';
      i++;
    } else if (editstr[n - 1 - i] == 'C') {
      res[i] = 'G';
      i++;
    } else if (editstr[n - 1 - i] == 'G') {
      res[i] = 'C';
      i++;
    } else if (editstr[n - 1 - i] == 'T') {
      res[i] = 'A';
      i++;
    } else
      assert(0);
  }
  res[n] = 0;
  return res;
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
	   bool first_in_pair, int* hits, int satisfying_alignments, bool improper_mapping)
/*
 * This function sets the strings output1 and output2 to be the output for the current read and if in sam mode its matepair
 * It is capable of outputting regular shrimp output, pretty print output, and sam output
 *
 * re is the read_entry for the current read
 * rh is the read_hit for the current read
 * rh_mp is the read_hit for the current reads mate pair
 *
 * paired is true if this read is paired
 * first is true if this is the first read in the pair
 *
 */
{
  assert(re != NULL);
  //assert((rh != NULL && rh->sfrp != NULL) || (rh_mp != NULL && rh_mp->sfrp != NULL));

  int thread_id = omp_get_thread_num();
  char ** output_buffer = &thread_output_buffer_filled[thread_id];
  char * output_buffer_end = thread_output_buffer[thread_id] + thread_output_buffer_sizes[thread_id] - 1 + 1;
  while ( (size_t)(output_buffer_end - *output_buffer) < thread_output_buffer_safety) { 
    //fprintf(stderr,"%d incrementing buffer, free space %llu\n",thread_id,output_buffer_end - *output_buffer );
    size_t new_size = thread_output_buffer_sizes[thread_id]+thread_output_buffer_increment;
    size_t filled = thread_output_buffer_filled[thread_id]-thread_output_buffer[thread_id];
    //fprintf(stderr, "there are %llu bytes used\n",filled);
    //thread_output_buffer[thread_id]=(char*)realloc(thread_output_buffer[thread_id],new_size);
    thread_output_buffer[thread_id] = (char *)
      my_realloc(thread_output_buffer[thread_id], new_size, thread_output_buffer_sizes[thread_id],
		 &mem_thread_buffer, "realloc thread_output_buffer");
    /*
      if (thread_output_buffer[thread_id]==NULL) {
      fprintf(stderr,"Hit output : realloc failed!\n");
      exit(1);
      }
    */
    thread_output_buffer_sizes[thread_id]=new_size;
    thread_output_buffer_filled[thread_id]=thread_output_buffer[thread_id]+filled;
    output_buffer = thread_output_buffer_filled+thread_id;
    output_buffer_end = thread_output_buffer[thread_id] + thread_output_buffer_sizes[thread_id] - 1 + 1;
  }	

  if (!Eflag) // old shrimp output format
    {
      char * tmp_output;

      if (rh != NULL) {
	int score=rh->sfrp->score;
	rh->sfrp->score=rh->score_full;
	tmp_output = output_normal(re->name, contig_names[rh->cn], rh->sfrp,
				   genome_len[rh->cn], shrimp_mode == MODE_COLOUR_SPACE, re->read[rh->st],
				   re->read_len, re->initbp[rh->st], rh->gen_st, Rflag);
	*output_buffer += snprintf(*output_buffer, output_buffer_end - *output_buffer, "%s\n", tmp_output);
	free(tmp_output);

	if (Pflag) { //pretty print output
	  tmp_output = output_pretty(re->name, contig_names[rh->cn], rh->sfrp,
				     genome_contigs[rh->cn], genome_len[rh->cn],
				     (shrimp_mode == MODE_COLOUR_SPACE), re->read[rh->st],
				     re->read_len, re->initbp[rh->st], rh->gen_st);
	  *output_buffer += snprintf(*output_buffer, output_buffer_end - *output_buffer, "%s\n", tmp_output);
	  free(tmp_output);
	}
	rh->sfrp->score=score;
      } else { // this is an unmapped read (part of a pair)
	*output_buffer += snprintf(*output_buffer, output_buffer_end - *output_buffer, ">%s\n", re->name);
      }
    }
  else // SAM output
    {
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
	int mapq = (rh != NULL ? rh->sfrp->mqv : 0);
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
	bool paired_alignment = paired_read && (rh!=NULL && rh_mp!=NULL && !improper_mapping); //paired mapping, not paired read!
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
#ifndef NDEBUG
	int stored_alignments = MIN(num_outputs,satisfying_alignments); //IH
#endif
	//if the read has no mapping or if not in half_paired mode and the mate has no mapping
	if (query_unmapped || (!half_paired && paired_read && mate_unmapped)) {
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
	assert( !paired_read || (!query_unmapped && !mate_unmapped) || half_paired);
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
			if (qual_delta!=33) { 
				int qual_len = strlen(re->qual); //not same as read_len, for color space reads... apperently.....
				int i;
				for (i=0; i<qual_len; i++) {
					qual[i]=qual[i]-qual_delta+33;
				}
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
		if (Qflag) {
			if (Bflag) {
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
					for (i = 0; i < rh->sfrp->rmapped/2; i++) {
						char temp = qual[i];
						qual[i]=qual[rh->sfrp->rmapped-i-1];
						qual[rh->sfrp->rmapped-i-1]=temp;
					}
				}
			} else if (compute_mapping_qualities) {
		    		strcpy(qual, rh->sfrp->qual);
				if (reverse_strand) {
					for (i = 0; i < rh->sfrp->rmapped/2; i++) {
						char temp = qual[i];
						qual[i]=qual[rh->sfrp->rmapped-i-1];
						qual[rh->sfrp->rmapped-i-1]=temp;
					}
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

		if (strcmp(rname, mrnm) == 0) {
		  mrnm = "=";
		  //printf("%d %d %c, %d %d %c\n",genome_start, genome_end,reverse_strand ? '-' : '+' , genome_start_mp, genome_end_mp, reverse_strand_mp ? '-' : '+');
		  int fivep = 0;
		  int fivep_mp = 0;
		  if (reverse_strand)
		    fivep = genome_end;
		  else
		    fivep = genome_start - 1;

		  if (reverse_strand_mp)
		    fivep_mp = genome_end_mp;
		  else
		    fivep_mp = genome_start_mp-1;

		  isize = (fivep_mp - fivep);
		} else { // map to different chromosomes
		  isize = 0;
		}
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
	*output_buffer += snprintf(*output_buffer, output_buffer_end - *output_buffer, "\tAS:i:%d",
				   rh->score_full);

	if (compute_mapping_qualities && !all_contigs) {
	  if (pair_mode == PAIR_NONE)
	  {
	    *output_buffer += snprintf(*output_buffer, output_buffer_end - *output_buffer, "\tZ0:i:%d\tZ1:i:%d",
		double_to_neglog(rh->sfrp->z0), double_to_neglog(rh->sfrp->z1));
	  }
	  else // paired mode
	  {
	    if (rh != NULL && rh_mp != NULL && !improper_mapping) {
	      *output_buffer += snprintf(*output_buffer, output_buffer_end - *output_buffer, "\tZ2:i:%d\tZ3:i:%d\tZ4:i:%d\tZ6:i:%d",
		  double_to_neglog(rh->sfrp->z2), double_to_neglog(rh->sfrp->z3),
		  double_to_neglog(rh->sfrp->pr_top_random_at_location), double_to_neglog(rh->sfrp->insert_size_denom));
	    } else {
	      *output_buffer += snprintf(*output_buffer, output_buffer_end - *output_buffer, "\tZ0:i:%d\tZ1:i:%d\tZ4:i:%d\tZ5:i:%d",
		  double_to_neglog(rh->sfrp->z0), double_to_neglog(rh->sfrp->z1),
		  double_to_neglog(rh->sfrp->pr_top_random_at_location), double_to_neglog(rh->sfrp->pr_missed_mp));
	    }
	  }
	}

	//*output_buffer += snprintf(*output_buffer, output_buffer_end - *output_buffer, "\tH0:i:%d\tH1:i:%d\tH2:i:%d",
	//			   hits[0],hits[1],hits[2]);
	*output_buffer += snprintf(*output_buffer, output_buffer_end - *output_buffer, "\tNM:i:%d",
				   rh->sfrp->mismatches+rh->sfrp->deletions+rh->sfrp->insertions);
	//*output_buffer += snprintf(*output_buffer, output_buffer_end - *output_buffer, "\tNH:i:%d\tIH:i:%d",
	//			   satisfying_alignments, stored_alignments);
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
	if (extra_sam_fields) {
	  char * editstr = alignment_edit_string(rh->sfrp->dbalign, rh->sfrp->qralign);
	  if (reverse_strand) {
	    char * tmp = reverse_alignment_edit_string(editstr);
	    free(editstr);
	    editstr = tmp;
	  }
	  *output_buffer += snprintf(*output_buffer, output_buffer_end - *output_buffer,
				     "\tZM:i:%d\tZR:i:%d\tZV:i:%d\tZH:i:%d\tZE:Z:%s",
				     rh->matches, rh->score_window_gen,
				     rh->score_vector, rh->sfrp->score,
				     editstr);
	  free(editstr);
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


static void
compute_unpaired_mqv(read_entry * re, read_hit * * hits, int n_hits)
{
  int i;
  double z1 = 0.0;

  for (i = 0; i < n_hits; i++) {
    z1 += hits[i]->sfrp->posterior;
  }
  for (i = 0; i < n_hits; i++) {
    hits[i]->sfrp->z0 = hits[i]->sfrp->posterior;
    hits[i]->sfrp->z1 = z1;

    hits[i]->sfrp->mqv = qv_from_pr_corr(hits[i]->sfrp->posterior / z1);
    if (hits[i]->sfrp->mqv < 4) hits[i]->sfrp->mqv = 0;
  }
}


static inline double
get_pr_insert_size(double insert_size)
{
  double res;
  res = normal_cdf(insert_size + 10, insert_size_mean, insert_size_stddev) - normal_cdf(insert_size - 10, insert_size_mean, insert_size_stddev);
  //int bucket = (insert_size - min_insert_size) / insert_histogram_bucket_size;
  //if (bucket < 0) bucket = 0;
  //if (bucket > 99) bucket = 99;
  //res = (double)insert_histogram[bucket] / (double)insert_histogram_load;
  if (res < 1e-200) res = 1e-200;
  //fprintf(stderr,"INVOKED WITH %e, %e %e, res %e\n",insert_size,insert_size_mean, insert_size_stddev,res);
  return res;
}


static void
compute_paired_mqv(pair_entry * pe)
{
  int i, j, nip;
  double z1[2];
  double z3;
  double insert_size_denom;
  double pr_top_random[3] = { 1.0, 1.0, 1.0 };
  double pr_missed_mp[2];
  double class_select_denom;

  // first, compute z1s for unpaired hits, which give unpaired posteriors
  for (nip = 0; nip < 2; nip++) {
    z1[nip] = 0;
    for (i = 0; i < pe->re[nip]->n_final_unpaired_hits; i++)
      z1[nip] += pe->re[nip]->final_unpaired_hits[i].sfrp->posterior;
    // and write them back
    for (i = 0; i < pe->re[nip]->n_final_unpaired_hits; i++) {
      pe->re[nip]->final_unpaired_hits[i].sfrp->z0 = pe->re[nip]->final_unpaired_hits[i].sfrp->posterior;
      pe->re[nip]->final_unpaired_hits[i].sfrp->z1 = z1[nip];
    }
  }

  // renormalize insert sizes into a distribution
  insert_size_denom = 0.0;
  for (i = 0; i < pe->n_final_paired_hits; i++) {
    double pr_ins = get_pr_insert_size(pe->final_paired_hits[i].insert_size);
    insert_size_denom += pr_ins;
  }
  // write denom to all sfr structs of paired mappings
  for (nip = 0; nip < 2; nip++) {
    for (i = 0; i < pe->final_paired_hit_pool_size[nip]; i++) {
      pe->final_paired_hit_pool[nip][i].sfrp->insert_size_denom = insert_size_denom;
    }
  }

  // compute paired posteriors
  z3 = 0.0;
  for (nip = 0; nip < 2; nip++) {
    for (i = 0; i < pe->final_paired_hit_pool_size[nip]; i++) {
      read_hit * rhp = &pe->final_paired_hit_pool[nip][i];
      double tmp = 0.0;
      for (j = 0; j < rhp->n_paired_hit_idx; j++) {
	read_hit_pair * rhpp = &pe->final_paired_hits[rhp->paired_hit_idx[j]];
	read_hit * rh_mpP = &pe->final_paired_hit_pool[1 - nip][rhpp->rh_idx[1 - nip]];
	double pr_ins = get_pr_insert_size(rhpp->insert_size);
	tmp += pr_ins * rh_mpP->sfrp->posterior;
      }
      tmp *= rhp->sfrp->posterior;
      if (tmp < 1e-200) tmp = 1e-200;
      rhp->sfrp->z2 = tmp;
      if (nip == 0) z3 += tmp;
    }
  }
  for (nip = 0; nip < 2; nip++) {
    for (i = 0; i < pe->final_paired_hit_pool_size[nip]; i++) {
      pe->final_paired_hit_pool[nip][i].sfrp->z3 = z3;
    }
  }

  // compute the probability that the top mapping in each class is random
  // for unpaired, use the one with max posterior
  for (nip = 0; nip < 2; nip++) {
    // find the mapping with max score
    int max_idx = 0;
    if (pe->re[nip]->n_final_unpaired_hits == 0) continue;
    for (i = 1; i < pe->re[nip]->n_final_unpaired_hits; i++)
      if (pe->re[nip]->final_unpaired_hits[i].sfrp->z0 > pe->re[nip]->final_unpaired_hits[max_idx].sfrp->z0)
	max_idx = i;
    // find pr a mapping of that score arises by chance
    pr_top_random[nip] = pr_random_mapping_given_score(pe->re[nip], pe->re[nip]->final_unpaired_hits[max_idx].sfrp->posterior_score);
    // write the coefficient back
    for (i = 0; i < pe->re[nip]->n_final_unpaired_hits; i++)
      pe->re[nip]->final_unpaired_hits[i].sfrp->pr_top_random_at_location = pr_top_random[nip];
    pr_top_random[nip] *= (double)total_genome_size;
    if (pr_top_random[nip] > 1)
      pr_top_random[nip] = 1.0;
  }

  // for paired, both mappings must be random
  for (i = 0; i < pe->n_final_paired_hits; i++) {
    double tmp = pr_random_mapping_given_score(pe->re[0], pe->final_paired_hit_pool[0][pe->final_paired_hits[i].rh_idx[0]].sfrp->posterior_score);
    tmp *= pr_random_mapping_given_score(pe->re[1], pe->final_paired_hit_pool[1][pe->final_paired_hits[i].rh_idx[1]].sfrp->posterior_score);
    tmp *= 1000;
    if (tmp < pr_top_random[2]) pr_top_random[2] = tmp;
  }
  // write it back
  for (i = 0; i < pe->n_final_paired_hits; i++) {
    pe->final_paired_hit_pool[0][pe->final_paired_hits[i].rh_idx[0]].sfrp->pr_top_random_at_location = pr_top_random[2];
    pe->final_paired_hit_pool[1][pe->final_paired_hits[i].rh_idx[1]].sfrp->pr_top_random_at_location = pr_top_random[2];
  }
  pr_top_random[2] *= (double)total_genome_size;
  if (pr_top_random[2] > 1)
    pr_top_random[2] = 1.0;

  // finally, for unpaired mappings, compute the prob the algo missed the mate mapping
  for (nip = 0; nip < 2; nip++) {
    pr_missed_mp[nip] = get_pr_missed(pe->re[1-nip]);
    // write it to sfrp
    for (i = 0; i < pe->re[nip]->n_final_unpaired_hits; i++)
      pe->re[nip]->final_unpaired_hits[i].sfrp->pr_missed_mp = pr_missed_mp[nip];
  }

  // compute denominator selecting between classes of mappings
  class_select_denom = 0.0;
  if (pe->re[0]->n_final_unpaired_hits > 0)
    class_select_denom += pr_top_random[1] * pr_top_random[2] * pr_missed_mp[0];
  if (pe->re[1]->n_final_unpaired_hits > 0)
    class_select_denom += pr_top_random[0] * pr_top_random[2] * pr_missed_mp[1];
  if (pe->n_final_paired_hits > 0)
    class_select_denom += pr_top_random[0] * pr_top_random[1];

  // DONE! ready to compute mqvs
  // unpaired:
  for (nip = 0; nip < 2; nip++) {
    for (i = 0; i < pe->re[nip]->n_final_unpaired_hits; i++) {
      read_hit * rhp = &pe->re[nip]->final_unpaired_hits[i];
      double p_corr = (pr_top_random[1-nip] * pr_top_random[2] * pr_missed_mp[nip] / class_select_denom) * (rhp->sfrp->z0 / rhp->sfrp->z1);
      rhp->sfrp->mqv = qv_from_pr_corr(p_corr);
      if (rhp->sfrp->mqv < 4) rhp->sfrp->mqv = 0;
    }
  }
  // paired:
  for (i = 0; i < pe->n_final_paired_hits; i++) {
    for (nip = 0; nip < 2; nip++) {
      read_hit * rhp = &pe->final_paired_hit_pool[nip][pe->final_paired_hits[i].rh_idx[nip]];
      double p_corr = (pr_top_random[0] * pr_top_random[1] / class_select_denom) * (rhp->sfrp->z2 / (/*insert_size_denom * */rhp->sfrp->z3));
      rhp->sfrp->mqv = qv_from_pr_corr(p_corr);
      if (rhp->sfrp->mqv < 4) rhp->sfrp->mqv = 0;
    }
  }
}


static void
add_sam_hit_counts(struct read_hit * h, int * sam_hit_counts)
{
  int edit_distance = h->sfrp->mismatches + h->sfrp->insertions + h->sfrp->deletions;
  if (0 <= edit_distance && edit_distance < 3) {
    sam_hit_counts[edit_distance]++;
  }
}


void
read_output(struct read_entry * re, struct read_hit * * hits_pass2, int n_hits_pass2)
{
  int i;
  int sam_hit_counts[3] = {0, 0, 0};
  int first, last;

  if (n_hits_pass2 == 0)
    return;

  if (Eflag) {
    for (i = 0; i < n_hits_pass2; i++) {
      add_sam_hit_counts(hits_pass2[i], sam_hit_counts);
    }
  }

  // Compute mapping qualities only if in unpaired mode.
  // Note: this procedure might be called in paired mode by setting save_outputs to false; no mqv in that case
  first = 0;
  last = n_hits_pass2;
  if (pair_mode == PAIR_NONE && compute_mapping_qualities) {
    compute_unpaired_mqv(re, hits_pass2, n_hits_pass2);
    if (single_best_mapping) {
      int max_idx = 0;
      for (i = 1; i < n_hits_pass2; i++)
	if (hits_pass2[i]->sfrp->mqv > hits_pass2[max_idx]->sfrp->mqv)
	  max_idx = i;
      first = max_idx;
      last = max_idx + 1;
    }
  }

  bool good_unpair = false;
  for (i = first; i < last; i++) {
    struct read_hit * rh = hits_pass2[i];
    if (pair_mode == PAIR_NONE) {
      hit_output(re, rh, NULL, false, sam_hit_counts, n_hits_pass2, false);
      if (rh->sfrp->mqv >= 10)
	good_unpair = true;
    } else {
      if (re->first_in_pair) {
	hit_output(re, rh, NULL, true, sam_hit_counts, n_hits_pass2, false);
	hit_output(re->mate_pair, NULL, rh, false, NULL, 0, false);
     } else {
      	hit_output(re->mate_pair, NULL, rh, true, NULL, 0, false);
      	hit_output(re, rh, NULL, false, sam_hit_counts, n_hits_pass2, false);
     }
    }
  }
  if (pair_mode == PAIR_NONE && good_unpair) {
#pragma omp atomic
    total_reads_matched_conf++;
  }
}


void
readpair_output_no_mqv(pair_entry * pe, struct read_hit_pair * hits_pass2, int n_hits_pass2)
{
  int i;
  int sam_hit_counts[2][3] = { {0, 0, 0}, {0, 0, 0} };

  if (Eflag) {
    for (i = 0; i < n_hits_pass2; i++) {
      add_sam_hit_counts(hits_pass2[i].rh[0], sam_hit_counts[0]);
      add_sam_hit_counts(hits_pass2[i].rh[1], sam_hit_counts[1]);
    }
  }

  for (i = 0; i < n_hits_pass2; i++) {
    struct read_hit * rh1 = hits_pass2[i].rh[0];
    struct read_hit * rh2 = hits_pass2[i].rh[1];
    uint bucket;

    hit_output(pe->re[0], rh1,  rh2, true, sam_hit_counts[0], n_hits_pass2, false);
    hit_output(pe->re[1], rh2,  rh1, false, sam_hit_counts[1], n_hits_pass2, false);

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
  }
}


// Find paired mapping containing given mapping of one read
// and max mqv mapping of the other
int
get_idx_mp_max_mqv(pair_entry * pe, int nip, int idx)
{
  assert(pe != NULL && idx < pe->final_paired_hit_pool_size[nip]);
  assert(pe->final_paired_hit_pool[nip][idx].n_paired_hit_idx > 0);

  int best_other_mqv = -1;
  int i, idx_pair = -1; // init not needed
  read_hit * rhp = &pe->final_paired_hit_pool[nip][idx];
  for (i = 0; i < rhp->n_paired_hit_idx; i++) {
    read_hit * rh_mpP = &pe->final_paired_hit_pool[1 - nip][pe->final_paired_hits[rhp->paired_hit_idx[i]].rh_idx[1 - nip]];
    if (rh_mpP->sfrp->mqv > best_other_mqv) {
      best_other_mqv = rh_mpP->sfrp->mqv;
      idx_pair = rhp->paired_hit_idx[i];
    }
  }
  return idx_pair;
}


void
readpair_output(pair_entry * pe)
{
  int i, nip;
  int sam_hit_counts[2][3] = { {0, 0, 0}, {0, 0, 0} };
  int first[3], last[3];

  if (Eflag) {
    for (nip = 0; nip < 2; nip++) {
      for (i = 0; i < pe->re[nip]->n_final_unpaired_hits; i++)
	add_sam_hit_counts(&pe->re[nip]->final_unpaired_hits[i], sam_hit_counts[nip]);
      for (i = 0; i < pe->final_paired_hit_pool_size[nip]; i++)
	add_sam_hit_counts(&pe->final_paired_hit_pool[nip][i], sam_hit_counts[nip]);
    }
  }

  first[0] = 0;
  last[0] = pe->re[0]->n_final_unpaired_hits;
  first[1] = 0;
  last[1] = pe->re[1]->n_final_unpaired_hits;
  first[2] = 0;
  last[2] = pe->n_final_paired_hits;
  if (compute_mapping_qualities) {
    compute_paired_mqv(pe);
    if (single_best_mapping
	&& (pe->n_final_paired_hits > 0 || pe->re[0]->n_final_unpaired_hits > 0 || pe->re[1]->n_final_unpaired_hits > 0)) {
      // find out the max mqv for either read
      int max_is_paired[2];
      int max_idx_unpaired[2] = { -1, -1 };
      int max_idx_paired[2] = { -1, -1 };
      int max_idx[2];
      int max_mqv_unpaired[2] = { -1, -1 };
      int max_mqv_paired[2] = { -1, -1 };
      int max_mqv[2];
      int best_nip, idx_pair;

      for (nip = 0; nip < 2; nip++) {
	for (i = 0; i < pe->re[nip]->n_final_unpaired_hits; i++)
	  if (pe->re[nip]->final_unpaired_hits[i].sfrp->mqv > max_mqv_unpaired[nip]) {
	    max_mqv_unpaired[nip] = pe->re[nip]->final_unpaired_hits[i].sfrp->mqv;
	    max_idx_unpaired[nip] = i;
	  }
	for (i = 0; i < pe->final_paired_hit_pool_size[nip]; i++)
	  if (pe->final_paired_hit_pool[nip][i].sfrp->mqv > max_mqv_paired[nip]) {
	    max_mqv_paired[nip] = pe->final_paired_hit_pool[nip][i].sfrp->mqv;
	    max_idx_paired[nip] = i;
	  }
      }

      if (!all_contigs)
      {
	// output top mapping in each class

	// for unpaired, take max
	for (nip = 0; nip < 2; nip++)
	  if (max_idx_unpaired[nip] >= 0) {
	    first[nip] = max_idx_unpaired[nip];
	    last[nip] = max_idx_unpaired[nip] + 1;
	  }

	// for paired, take the max mqv, then the mate with the max mqv
	if (max_mqv_paired[0] > max_mqv_paired[1])
	  best_nip = 0;
	else
	  best_nip = 1;
	if (max_mqv_paired[best_nip] >= 0) {
	  idx_pair = get_idx_mp_max_mqv(pe, best_nip, max_idx_paired[best_nip]);
	  first[2] = idx_pair;
	  last[2] = idx_pair + 1;
	}
      }
      else // all-contigs
      {
	// output top mapping over all classes

	for (nip = 0; nip < 2; nip++) {
	  if (max_mqv_unpaired[nip] > max_mqv_paired[nip]) {
	    max_mqv[nip] = max_mqv_unpaired[nip];
	    max_is_paired[nip] = 0;
	    max_idx[nip] = max_idx_unpaired[nip];
	  } else {
	    max_mqv[nip] = max_mqv_paired[nip];
	    max_is_paired[nip] = 1;
	    max_idx[nip] = max_idx_paired[nip];
	  }
	}

	if (max_mqv[0] >= max_mqv[1])
	  // will output best hit for read0
	  best_nip = 0;
	else
	  best_nip = 1;

	if (max_is_paired[best_nip] == 1)
	{
	  // output a pair containing max
	  // choose max mq for the other read

	  idx_pair = get_idx_mp_max_mqv(pe, best_nip, max_idx[best_nip]);

	/*
	  // update the histogram of insert sizes
	  int bucket = (pe->final_paired_hits[idx_pair].insert_size - min_insert_size) / insert_histogram_bucket_size;
	  if (bucket < 0) bucket = 0;
	  if (bucket > 99) bucket = 99;
#pragma omp atomic
	  insert_histogram[bucket]++;
#pragma omp atomic
	  insert_histogram_load++;
	*/

	  last[0] = 0;
	  last[1] = 0;
	  first[2] = idx_pair;
	  last[2] = idx_pair + 1;
	}
	else
	{
	  // max mqv is in an unpaired hit
	  // see if it can be paired with best one from read1
	  int idx_best_other = -1;
	  double max_other_z0 = 0.0;
	  for (i = 0; i < pe->re[1 - best_nip]->n_final_unpaired_hits; i++) {
	    if (pe->re[1 - best_nip]->final_unpaired_hits[i].sfrp->z0 > max_other_z0) {
	      max_other_z0 = pe->re[1 - best_nip]->final_unpaired_hits[i].sfrp->z0;
	      idx_best_other = i;
	    }
	  }
	  int best_other_mqv = -1;
	  if (idx_best_other >= 0)
	    //best_other_mqv = qv_from_pr_corr(exp((-max_other_z0 + pe->re[1 - best_nip]->final_unpaired_hits[idx_best_other].sfrp->z1) / 1000));
	    best_other_mqv = qv_from_pr_corr(max_other_z0 / pe->re[1 - best_nip]->final_unpaired_hits[idx_best_other].sfrp->z1);
	  if (!improper_mappings || max_mqv_unpaired[best_nip] < 10 || best_other_mqv < 10) {
	    // output unpaired hit
	    last[2] = 0;
	    last[1-best_nip] = 0;
	    first[best_nip] = max_idx[best_nip];
	    last[best_nip] = max_idx[best_nip] + 1;
	  } else {
	    // unpaired qv >= 10:
	    // pair up hits, possibly across chromosomes
	    // for this, create a new read_hit_pair entry
	    read_hit_pair * rhpp;
	    pe->final_paired_hits = (read_hit_pair *)
		  my_realloc(pe->final_paired_hits,
		      (pe->n_final_paired_hits + 1) * sizeof(pe->final_paired_hits[0]),
		      pe->n_final_paired_hits * sizeof(pe->final_paired_hits[0]),
		      &mem_mapping, "final_paired_hits [%s,%s]", pe->re[0]->name, pe->re[1]->name);
	    pe->n_final_paired_hits++;
	    rhpp = &pe->final_paired_hits[pe->n_final_paired_hits - 1];
	    rhpp->rh[best_nip] = &pe->re[best_nip]->final_unpaired_hits[max_idx[best_nip]];
	    rhpp->rh[1 - best_nip] = &pe->re[1 - best_nip]->final_unpaired_hits[idx_best_other];
	    rhpp->score_max = rhpp->rh[0]->score_max + rhpp->rh[1]->score_max;
	    rhpp->insert_size = get_insert_size(rhpp->rh[best_nip], rhpp->rh[1 - best_nip]);
	    rhpp->improper_mapping = true;
	    last[0] = 0;
	    last[1] = 0;
	    first[2] = pe->n_final_paired_hits - 1;
	    last[2] = pe->n_final_paired_hits;
	  }
	}
      }

    }
  }

  // time to output stuff
  bool good_pair = false;
  bool good_unpair = false;
  for (i = first[2]; i < last[2]; i++) {
    read_hit * rh1 = pe->final_paired_hits[i].rh[0];
    read_hit * rh2 = pe->final_paired_hits[i].rh[1];

    if (rh1 == NULL) rh1 = &pe->final_paired_hit_pool[0][pe->final_paired_hits[i].rh_idx[0]];
    if (rh2 == NULL) rh2 = &pe->final_paired_hit_pool[1][pe->final_paired_hits[i].rh_idx[1]];

    hit_output(pe->re[0], rh1,  rh2, true, sam_hit_counts[0], pe->n_final_paired_hits, pe->final_paired_hits[i].improper_mapping);
    hit_output(pe->re[1], rh2,  rh1, false, sam_hit_counts[1], pe->n_final_paired_hits, pe->final_paired_hits[i].improper_mapping);

    if (!pe->final_paired_hits[i].improper_mapping && (rh1->sfrp->mqv >= 10 || rh2->sfrp->mqv >= 10)) {
      good_pair = true;
    } else if (pe->final_paired_hits[i].improper_mapping && (rh1->sfrp->mqv >= 10 || rh2->sfrp->mqv >= 10)) {
      good_unpair = true;
    }

    if (Xflag && !pe->final_paired_hits[i].improper_mapping) {
      // update the histogram of insert sizes
      int bucket = (pe->final_paired_hits[i].insert_size - min_insert_size) / insert_histogram_bucket_size;
      if (bucket < 0) bucket = 0;
      if (bucket > 99) bucket = 99;
#pragma omp atomic
      insert_histogram[bucket]++;
#pragma omp atomic
      insert_histogram_load++;
    }
  }
  for (nip = 0; nip < 2; nip++)
    for (i = first[nip]; i < last[nip]; i++) {
      read_entry * rep = pe->re[nip];
      read_hit * rhp = &rep->final_unpaired_hits[i];

      if (rep->first_in_pair) {
	hit_output(rep, rhp, NULL, true, sam_hit_counts[0], rep->n_final_unpaired_hits, false);
	hit_output(rep->mate_pair, NULL, rhp, false, NULL, 0, false);
      } else {
  	hit_output(rep->mate_pair, NULL, rhp, true, NULL, 0, false);
  	hit_output(rep, rhp, NULL, false, sam_hit_counts[1], rep->n_final_unpaired_hits, false);
      }

      if (rhp->sfrp->mqv >= 10) {
	good_unpair = true;
      }
    }

  if (good_pair) {
#pragma omp atomic
    total_pairs_matched_conf++;
  } else if (good_unpair) {
#pragma omp atomic
    total_reads_matched_conf++;
  }
}
