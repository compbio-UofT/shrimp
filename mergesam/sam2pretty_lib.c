//aug 10 2011
#include <stdio.h>
#include <ctype.h>
#include <assert.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

#include "sam2pretty_lib.h"
#include "render.h"

double inv_tnlog(int x) { 
        return exp(-(x/1000.0)); 
} 
 
int tnlog(double x) { 
        return (int)(1000.0 * -log(x)); 
}


bool sam2pretty_lib_verbose=false;

void pretty_stats(pretty * pa) {
	if (pa->pretty_length==0) {
		return;
	}
	assert(pa->pretty_read_string!=NULL);
	assert(pa->pretty_genome_string!=NULL);
	assert(pa->clipped==0);
	assert(pa->mismatches==0);
	assert(pa->matches==0);
	assert(pa->insertions==0);
	assert(pa->deletions==0);
	assert(pa->clipped_read_start!=0);
	assert(pa->clipped_read_end!=0);
	assert(pa->unclipped_read_length!=0);
	pa->clipped=(pa->clipped_read_start-1) + (pa->unclipped_read_length-pa->clipped_read_end);
	int32_t j, length;
	length=strlen(pa->pretty_read_string);
	//printf("%s,%d to %d\n",pa->cigar,pa->clipped_read_start-1,(length-(pa->unclipped_read_length-pa->clipped_read_end)));
	for (j=pa->clipped_read_start-1; j<pa->clipped_read_end; j++) {
		char g = pa->pretty_genome_string[j];
		char r = pa->pretty_read_string[j];
		char rc = r>90 ? r-32 : r;
		char gc = g>90 ? g-32 : g;
		if (gc=='-') {
			pa->insertions++;
		} else if (rc=='-') {
			pa->deletions++;
		} else if (rc=='.') {
			pa->skipped++;
		} else if (rc==gc) {
			pa->matches++;
		} else {
			pa->mismatches++;
		}
	}

}


int pretty_get_flag(pretty * pa ) {
	int flag=0;
	//0x0001 the read is paired in sequencing, no matter whether it is mapped in a pair
	flag|= (pa->paired_sequencing ? 0x0001 : 0);
	//0x0002 the read is mapped in a proper pair (depends on the protocol, normally inferred during alignment) 1
	flag|=(pa->proper_pair ? 0x0002 : 0);
	//0x0004 the query sequence itself is unmapped
	flag|=(pa->mapped ? 0 : 0x0004);
	//0x0008 the mate is unmapped 1
	flag|=(pa->mp_mapped ? 0 : 0x0008);
	//0x0010 strand of the query (0 for forward; 1 for reverse strand)
	flag|=(pa->reverse ? 0x0010 : 0);
	//0x0020 strand of the mate 1
	flag|=(pa->mp_reverse ? 0x0020 : 0);
	//0x0040 the read is the first read in a pair 1,2
	flag|=(pa->first_in_pair ? 0x0040 : 0);
	//0x0080 the read is the second read in a pair 1,2
	flag|=(pa->second_in_pair ? 0x0080 : 0);
	//0x0100 the alignment is not primary (a read having split hits may have multiple primary alignment records)
	flag|=(pa->primary_alignment ? 0x0100 : 0);
	//0x0200 the read fails platform/vendor quality checks
	flag|=(pa->platform_quality_fail ? 0x0200 : 0);
	//0x0400 the read is either a PCR duplicate or an optical duplicate
	flag|=(pa->pcr_duplicate ? 0x0400 : 0);
	return flag;
}

void pretty_match(pretty * pa) {
	if (pa->pretty_length==0) {
		return;
	}
	assert(pa->pretty_read_string!=NULL);
	assert(pa->pretty_genome_string!=NULL);
	char * g = pa->pretty_genome_string;
	char * r = pa->pretty_read_string;
	int32_t pretty_genome_string_length=strlen(g);
	int32_t pretty_read_string_length=strlen(r);
	//printf("|%s|\n|%s|\n",g,r);
	assert(pretty_genome_string_length==pretty_read_string_length);	
	char * m=(char*)malloc(sizeof(char)*(pretty_read_string_length+1));
	if (m==NULL){ 
		fprintf(stderr, "Malloc for pretty_match_string failed!\n");
		exit(1);
	}
	int32_t i;
	
	for (i=0; i<pretty_read_string_length; i++){
		char rc = r[i]>90 ? r[i]-32 : r[i];
		char gc = g[i]>90 ? g[i]-32 : g[i];
		if (r[i]>96) {
			if (gc==rc) {
				m[i]='X';
			} else {
				m[i]='x';
			}
		} else if (g[i]=='-' || r[i]=='-' || gc!=rc) {
			m[i]=' ';
		} else if (pa->has_cs_edit_string && pa->pretty_read_string[i]>96) {
			m[i]='X';
		} else {
			m[i]='|';
		} 
	}
	m[pretty_read_string_length]='\0';
	pa->pretty_match_string=m;
	return;
}


void make_pretty_edit_string(const char * pretty_match, const char * pretty_read, const char  * pretty_genome , int clipped_read_start, int clipped_read_end, char * target) {
	target[0]='\0';
	int strech=0;
	int i;
	bool read_gap=false;
	bool ref_gap=false;
	int read_gaps=0;
	for (i=clipped_read_start-1; i<read_gaps+clipped_read_end; i++) {
		//handle gaps in the genome
		if (pretty_genome[i]!='-' && ref_gap) {
			sprintf(target+strlen(target),")");
			ref_gap=false;
		}	
		if (pretty_match[i]=='|') {
			strech++;
		} else if (pretty_match[i]==' ') {
			if (strech!=0) {
				sprintf(target+strlen(target),"%d",strech);
			}
			if (pretty_genome[i]=='-' && ref_gap==false) {
				sprintf(target+strlen(target),"(");
				ref_gap=true;
			}
			sprintf(target+strlen(target),"%c",pretty_read[i]);
			strech=0;
		} else if (pretty_match[i]=='x' || pretty_match[i]=='X') {
			if (strech!=0) {
				sprintf(target+strlen(target),"%d",strech);
			}
			sprintf(target+strlen(target),"x");
			if (pretty_genome[i]=='-' && ref_gap==false) {
				sprintf(target+strlen(target),"(");
				ref_gap=true;
			}
			if (pretty_match[i]=='x') {
				sprintf(target+strlen(target),"%c",pretty_read[i]>96 ? pretty_read[i]-32 : pretty_read[i]);
				strech=0;
			} else {
				strech=1;
			}
		} else {
			fprintf(stderr,"Cannot find case for %c !!\n",pretty_match[i]);
			exit(1);
		}
		//handle gaps in the read
		if (pretty_read[i]=='-') {
			if (read_gap==false) {
				read_gap=true;
			} 
			read_gaps++;
		} else if (read_gap) {
			read_gap=false;
		}
	}
	if (strech!=0) {
		sprintf(target+strlen(target),"%d",strech);

	}
		if (ref_gap) {
			sprintf(target+strlen(target),")");
			ref_gap=false;
		}	
		if (read_gap) {
			read_gap=false;
		}
		
}

void pretty_print(pretty * pa) {
	char strand;
	int32_t genome_coordinate_1, genome_coordinate_2;
	if (!pa->reverse){ 
		strand='+';
		genome_coordinate_1=pa->genome_start_padded;
		genome_coordinate_2=pa->genome_end_padded;
	} else {
		strand='-';
		genome_coordinate_2=pa->genome_start_padded;
		genome_coordinate_1=pa->genome_end_padded;
	}
	char suffix[3]={'\0','\0','\0'};
	if (pa->first_in_pair) {
		strcpy(suffix,"/1");
	} else if (pa->second_in_pair) {
		strcpy(suffix,"/2");
	} else if (pa->paired_sequencing) {
		strcpy(suffix,"/?");
	}
	size_t pretty_edit_string_buffer_length=1;
	if (pa->pretty_read_string!=NULL && pa->pretty_match_string!=NULL) {
		pretty_edit_string_buffer_length=strlen(pa->read_string)*2;
	}
	char pretty_edit_string[pretty_edit_string_buffer_length];
	pretty_edit_string[0]='\0';
	if (pa->pretty_read_string!=NULL && pa->pretty_match_string!=NULL) {
		make_pretty_edit_string(pa->pretty_match_string,pa->pretty_read_string,pa->pretty_genome_string,pa->clipped_read_start,pa->clipped_read_end,pretty_edit_string);
	}
	if (pa->genome_start_padded!=0 && pa->genome_end_padded!=0 && pa->clipped_read_start !=0 && pa->clipped_read_end !=0 && pa->unclipped_read_length!=0)		
		printf(">%s%s\t%s\t%c\t%d\t%d\t%d\t%d\t%d\t%d\t%s\n",pa->read_name,suffix,pa->reference_name,strand,
			pa->genome_start_unpadded,pa->genome_end_unpadded,pa->clipped_read_start,pa->clipped_read_end,
			pa->unclipped_read_length,pa->score,pretty_edit_string);	
	if (pa->pretty_length>0) {
	if (sam2pretty_lib_verbose && pa->genome_start_padded!=0 && pa->genome_end_padded!=0 && pa->clipped_read_start !=0 && pa->clipped_read_end !=0 && pa->unclipped_read_length!=0)
		printf(">> Pretty SAM Alignment:\n");	
	if (pa->colour_space) {
		if (pa->pretty_genome_string!=NULL)
			printf("G:\t%d\t %s\t%d\n",genome_coordinate_1,pa->pretty_genome_string,genome_coordinate_2);
		if (pa->pretty_match_string!=NULL)
			printf("M:\t\t %s\n",pa->pretty_match_string);
		if (pa->pretty_read_string!=NULL)
			printf("T:\t%d\t %s\t%d\n",pa->clipped_read_start,pa->pretty_read_string,pa->clipped_read_end);
		if (pa->pretty_cs_string!=NULL)
			printf("R:\t%d\t%s\t%d\n",1,pa->pretty_cs_string,pa->unclipped_read_length);
		if (pa->has_cs_qualities) {
			printf("CQ:\t\t %s\t\n",pa->pretty_cs_qualities);
		}
		if (pa->pretty_read_qualities!=NULL)
			printf("Q:\t\t %s\t\n",pa->pretty_read_qualities);
	} else {
		if (pa->pretty_genome_string!=NULL)
			printf("G:\t%d\t%s\t%d\n",genome_coordinate_1,pa->pretty_genome_string,genome_coordinate_2);
		if (pa->pretty_match_string!=NULL)
			printf("M:\t\t%s\n",pa->pretty_match_string);
		if (pa->pretty_read_string!=NULL)
			printf("R:\t%d\t%s\t%d\n",pa->clipped_read_start,pa->pretty_read_string,pa->clipped_read_end);
		if (pa->pretty_read_qualities!=NULL)
			printf("Q:\t\t%s\t\n",pa->pretty_read_qualities);
	}
	}
	//pa->paired_sequencing = ((entry->core.flag & 0x0001) != 0) ? true : false ;
	//pa->proper_pair = ((entry->core.flag & 0x0002) != 0) ? true : false ;
	//pa->mapped = ((entry->core.flag & 0x0004) != 0) ? true : false ;
	//pa->mp_mapped = ((entry->core.flag & 0x0008) != 0) ? true : false ;
	//pa->reverse = ((entry->core.flag & 0x0010) != 0) ? true : false ;
	//pa->mp_reverse = ((entry->core.flag & 0x0020) != 0) ? true : false ;
	//pa->first_in_pair = ((entry->core.flag & 0x0040) != 0) ? true : false ;
	//pa->second_in_pair = ((entry->core.flag & 0x0080) != 0) ? true : false ;
	//pa->primary_alignment = ((entry->core.flag & 0x0100) != 0) ? true : false ;
	//pa->platform_quality_fail = ((entry->core.flag & 0x0200) != 0) ? true : false ;
	//pa->pcr_duplicate = ((entry->core.flag & 0x0400) != 0) ? true : false ;
	//print flags
	if (sam2pretty_lib_verbose) {
		printf("> SAM FLAG info:\n");
		printf(">>\tPaired Seqs:\t%s\tProper Pair:\t%s\n",pa->paired_sequencing ? "True" : "False",pa->proper_pair ? "True" : "False");
		printf(">>\tQuery Mapped:\t%s\tMate Mapped:\t%s\n",pa->mapped ? "True" : "False", pa->mp_mapped ? "True" : "False");
		printf(">>\tQuery Strand:\t%c\tMate Strand:\t%c\n",pa->reverse ? '-' : '+', pa->mp_reverse ? '-' : '+');
		printf(">>\tFirst in Pair:\t%s\tSecond in Pair:\t%s\n",pa->first_in_pair ? "True" : "False", pa->second_in_pair ? "True" : "False");
		printf(">>\tPrimary Algn:\t%s\tPltf. Q. Fail:\t%s\n",pa->primary_alignment ? "True" : "False", pa->platform_quality_fail ? "True" : "False");	
		printf(">>\tPCR duplicate:\t%s\n", pa->pcr_duplicate ? "True" : "False");
		//int32_t clipped;
		//int32_t mismatches;
		//int32_t deletions;
		//int32_t insertions;
		//int32_t matches;
		if (pa->pretty_length>0) {
			printf("> Internally Computed Stats:\n");
			printf(">>\tClipped:\t%d\tMatches:\t%d\tMisMatches:\t%d\n",pa->clipped,pa->matches,pa->mismatches);
			printf(">>\tInsertions:\t%d\tDeletions:\t%d\n",pa->insertions,pa->deletions);
			printf(">>\tEdit distance:\t%d\n",pa->mismatches+pa->insertions+pa->deletions+pa->clipped);
		}
		if (pa->edit_distance!=-1) {
			printf("> SAM Edit distance:\t%d\n",pa->edit_distance);
			if (pa->pretty_length>0){
				if (pa->edit_distance!=pa->mismatches+pa->insertions+pa->deletions) {
					fprintf(stderr,"read %s: SAM edit distance and calculated edit distance do not match!\n",pa->read_name);
				}
			}
		}
	}
}

pretty * pretty_new() {
	pretty * ret=(pretty*)malloc(sizeof(pretty));
	if (ret==NULL) {
		fprintf(stderr,"Failed to malloc new pretty\n");
		exit(1);
	}
	memset(ret,0,sizeof(pretty));
	return ret;
}

pretty * pretty_copy_base(pretty * parent) {
	pretty * child = pretty_new();
	*child = *parent;
	if (child->read_string!=NULL) {
		child->read_string=strdup(child->read_string);
	}
	if (child->read_qualities!=NULL) {
		child->read_qualities=strdup(child->read_qualities);
	}
	if (child->cs_qualities!=NULL) {
		child->cs_qualities=strdup(child->cs_qualities);
	}
	if (child->cs_string!=NULL) {	
		child->cs_string=strdup(child->cs_string);
	}
	if (child->read_group!=NULL) {
		child->read_group=strdup(child->read_group);
	}
	if (child->cs_edit_string!=NULL) {
		child->cs_edit_string=strdup(child->cs_edit_string);
	}	
	child->pretty_genome_string=NULL;
	child->pretty_read_string=NULL;
	child->pretty_read_qualities=NULL;
	child->pretty_cs_string=NULL;
	child->pretty_cs_qualities=NULL;
	child->cigar=NULL;
	if (child->reference_name!=NULL) {
		child->reference_name=strdup(child->reference_name);
	}
	if (child->mate_reference_name!=NULL) {
		child->mate_reference_name=strdup(child->mate_reference_name);
	}
	child->pretty_match_string=NULL;
	child->cigar_ops=NULL;
	child->cigar_lengths=NULL;
	if (child->read_name!=NULL) {
		child->read_name=strdup(child->read_name);
	}
	child->next=NULL;

	/*if (child->pretty_genome_string!=NULL) {
		child->pretty_genome_string=strdup(child->pretty_genome_string);
	}
	if (child->pretty_read_string!=NULL) {
		child->pretty_read_string=strdup(child->pretty_read_string);
	}
	if (child->pretty_read_qualities!=NULL) {
		child->pretty_read_qualities=strdup(child->pretty_read_qualities);
	}
	if (child->pretty_cs_string!=NULL) {
		child->pretty_cs_string=strdup(child->pretty_cs_string);
	}*/
	/*if (child->cigar!=NULL) {
		child->cigar=strdup(child->cigar);
	}*/
	/*if (child->num_cigar>0) {
		child->cigar_ops=(char*)malloc(sizeof(char)*child->num_cigar);
		child->cigar_lengths=(uint32_t*)malloc(sizeof(uint32_t)*child->num_cigar);
		int32_t i;
		for (i=0; i<child->num_cigar; i++) {
			child->cigar_ops[i]=parent->cigar_ops[i];
			child->cigar_lengths[i]=parent->cigar_lengths[i];
		}
	}*/
	/*if (child->pretty_match_string!=NULL) {
		child->pretty_match_string=strdup(child->pretty_match_string);
	}*/
	return child;
}

void pretty_free_fast(pretty * pa) {
	if (pa->read_name!=NULL) {
		free(pa->read_name);
	}
	if (pa->sam_string!=NULL) {
		free(pa->sam_string);
	}
	free(pa);
}

void pretty_free(pretty * pa) {
	//if (pa->genome_string != NULL) {
	//	free(pa->genome_string);
	//}
#ifdef SAM2PRETTY_DEBUG
	fprintf(stderr,"SAM2PRETTY %d FREE\n",pa->id);
#endif
	if (pa->pretty_genome_string != NULL) {
		free(pa->pretty_genome_string);
	}
	if (pa->read_string != NULL) {
		free(pa->read_string);
	}
	if (pa->pretty_read_string != NULL) {
		free(pa->pretty_read_string);
	}
	if (pa->read_qualities != NULL ) {
		free(pa->read_qualities);
	}
	if (pa->cs_qualities != NULL ) {
		free(pa->cs_qualities);
	}
	if (pa->pretty_read_qualities != NULL) {
		free(pa->pretty_read_qualities);
	}
	if (pa->pretty_cs_qualities != NULL) {
		free(pa->pretty_cs_qualities);
	}
	if (pa->cs_string != NULL) {
		free(pa->cs_string);
	}
	if (pa->pretty_cs_string != NULL) {
		free(pa->pretty_cs_string);
	}
	if (pa->cigar != NULL) {
		free(pa->cigar);
	}
	if (pa->reference_name != NULL) {
		free(pa->reference_name);
	}
	if (pa->mate_reference_name != NULL) {
		free(pa->mate_reference_name);
	}
	if (pa->read_name != NULL) {
		free(pa->read_name);
	}
	if (pa->cigar_ops != NULL) {
		free(pa->cigar_ops);
	}
	if (pa->cigar_lengths != NULL) {
		free(pa->cigar_lengths);
	}
	if (pa->pretty_match_string != NULL) {
		free(pa->pretty_match_string);
	}
	if (pa->has_read_group) {
		free(pa->read_group);
	}
	if (pa->has_cs_edit_string) {
		free(pa->cs_edit_string);
	}
	if (pa->sam_string!=NULL) {
		free(pa->sam_string);
	}
	if (pa->r2!=NULL) {
		free(pa->r2);
	}
	memset(pa,0,sizeof(pretty));
	free(pa);	
	return;				
}



void reverse_complement(char* sequence) {
	int32_t length=strlen(sequence);
	char buffer[length+1];
	int32_t i;
	for (i=0; i<length; i++) {
		char c = sequence[i];
		char d;
		switch (c) {
			case 'A': d='T'; break;		
			case 'a': d='t'; break;		
			case 'C': d='G'; break;		
			case 'c': d='g'; break;		
			case 'G': d='C'; break;		
			case 'g': d='c'; break;		
			case 'T': d='A'; break;		
			case 't': d='a'; break;		
			case 'N': d='N'; break;		
			case 'n': d='n'; break;		
			case '-': d='-'; break;		
			default: d=c; break;
		}
		buffer[length-1-i]=d;
	}
	buffer[length]='\0';
	strcpy(sequence,buffer);
}




//int32_t get_cigar_parse(char* cigar_ops,uint32_t * cigar_lengths, int32_t num_cigar, char* o_cigar) {
void pretty_cigar_parse(pretty * pa) {
	assert(pa->cigar!=NULL);
	assert(pa->cigar_ops!=NULL);
	assert(pa->cigar_lengths!=NULL);
	char cigar[strlen(pa->cigar)+1];
	strcpy(cigar,pa->cigar);
	int32_t last_seen_letter=-1;
	int32_t ops_so_far=0; int32_t i=0;
	while (cigar[i]!='\0') {
		assert(ops_so_far<pa->num_cigar);
		//check if its a letter
		if (cigar[i]>57) {
			pa->cigar_ops[ops_so_far]=cigar[i];
			cigar[i]='\0';
			pa->cigar_lengths[ops_so_far]=atoi(cigar+last_seen_letter+1);
			pa->pretty_length+=pa->cigar_lengths[ops_so_far]; //upper bound
			last_seen_letter=i;
			ops_so_far++;
		}
		i++;
	}
	return;
}


void calculate_genome_end(pretty * pa) {
	if (pa->num_cigar==0) {
		fprintf(stderr,"this read has no cigar string! %s\n",pa->read_name);	
		exit(1);
	}
	pa->genome_start_padded=pa->genome_start_unpadded;
	pa->genome_end_padded=pa->genome_start_unpadded;
	pa->genome_end_unpadded=pa->genome_start_unpadded;
	assert(pa->num_cigar>0);
	int32_t i;
	bool start = true;
	for (i=0; i<pa->num_cigar; i++) {
		int32_t length=pa->cigar_lengths[i];
		switch (pa->cigar_ops[i]) {
			case 'N':
			case 'D':
			case 'M':
				start=false;
				pa->genome_end_unpadded+=length;
				pa->genome_end_padded+=length;
				break;
			case 'S':
			case 'H':
				if (start) {
					pa->genome_start_padded-=length;
				} else {
					pa->genome_end_padded+=length;	
				}
				break;
			case 'P':
			case 'I':
				start=false;
				break;
			default:
				fprintf(stderr,"Error calculating pretty genome,%s\n",pa->cigar);
				exit(1);
		}
	}
	pa->genome_end_padded--;
	pa->genome_end_unpadded--;
	return;	
}
void calculate_insert_size(pretty * pa, pretty * pa_mp) {
	int isize=0;
	if (pa->reference_name==NULL || pa_mp->reference_name==NULL) {
		fprintf(stderr,"Cannot calculate insert size when reference names are not defined! %s, %s\n",pa->read_name,pa_mp->read_name);
		exit(1);
	}
	if (strcmp(pa->reference_name,pa_mp->reference_name) == 0) {
		pa->mate_reference_name= "=";
		pa_mp->mate_reference_name="=";
		int fivep = 0;
		int fivep_mp = 0;
		if (pa->reverse)
			fivep = pa->genome_end_unpadded;
		else
			fivep = pa->genome_start_unpadded - 1;

		if (pa_mp->reverse)
			fivep_mp = pa_mp->genome_end_unpadded;
		else
			fivep_mp = pa_mp->genome_start_unpadded-1;

		isize = (fivep_mp - fivep);
	} else { // map to different chromosomes
		isize = 0;
	}
	
	pa->isize=isize;
	pa_mp->isize=-isize;
	return;
}

//char* pretty_genome(char* cigar, int32_t num_cigar, int32_t strand, int32_t * genome_start, int32_t * genome_end) {
void pretty_genome(pretty * pa) {
	if (pa->pretty_length==0) {
		return;
	}
	assert(pa->genome_sequence!=NULL);
	assert(pa->genome_length!=0);
	//Allocate space for pretty string
	assert(pa->pretty_length>0);
	pa->pretty_genome_string = (char*)malloc(sizeof(char)*(pa->pretty_length+1));
	if (pa->pretty_genome_string==NULL) {
		fprintf(stderr,"Malloc failed for pretty_genome_string in pretty_genome!\n");
		exit(1);
	}
	//pa->genome_string = (char*)malloc(sizeof(char)*(pa->pretty_length+1));
	//if (pa->genome_string==NULL) {
	//	fprintf(stderr,"Malloc failed for genome_string in pretty_genome!\n");
	//	exit(1);
	//}
	//printf("allocated %d\n",pretty_seq_length+1);
	if (pa->num_cigar==0) {
		pa->genome_start_padded=pa->genome_start_unpadded;
		pa->genome_end_padded=pa->genome_start_unpadded+pa->pretty_length;
		pa->genome_end_unpadded=pa->genome_end_unpadded;
		strncpy(pa->pretty_genome_string,pa->genome_sequence+pa->genome_start_unpadded-1,pa->pretty_length);
		pa->pretty_genome_string[pa->pretty_length]='\0';	
		//strncpy(pa->genome_string,pa->genome_sequence+pa->genome_start_unpadded-1,pa->pretty_length);	
		//pa->genome_string[pa->pretty_length]='\0';	
		return;
	}
	pa->genome_start_padded=pa->genome_start_unpadded;
	pa->genome_end_padded=pa->genome_start_unpadded;
	pa->genome_end_unpadded=pa->genome_start_unpadded;
	assert(pa->num_cigar>0);
	int32_t j,i;
	bool start = true;
	int32_t current_pretty_genome_string_index=0;
	//int32_t current_genome_string_index=0;
	int32_t current_genome_index=pa->genome_start_unpadded-1; //1 based -> 0 based
	if (current_genome_index<0 || current_genome_index>=pa->genome_length) {
		fprintf(stderr,"The read \"%s\" has run off the genome!\n",pa->read_name);
		exit(1);
	}
	//assert(current_genome_index>=0);
	//assert(current_genome_index<pa->genome_length);
	for (i=0; i<pa->num_cigar; i++) {
		int32_t length=pa->cigar_lengths[i];
		switch (pa->cigar_ops[i]) {
			case 'N':
			case 'D':
			case 'M':
				start=false;
				pa->genome_end_unpadded+=length;
				pa->genome_end_padded+=length;
				for (j=0; j<length; j++) {
					if (current_genome_index>=pa->genome_length) {
						fprintf(stderr,"Read \"%s\" has run off the end of genome, index %d vs genome length %d!\n",pa->read_name,current_genome_index+1, pa->genome_length);
						exit(1);
					}
					pa->pretty_genome_string[current_pretty_genome_string_index++]=pa->genome_sequence[current_genome_index];
					//pa->genome_string[current_genome_string_index++]=pa->genome_sequence[current_genome_index];
					current_genome_index++;
				}
				break;
			case 'S':
			case 'H':
				if (start) {
					current_genome_index-=length;
					pa->genome_start_padded-=length;
				} else {
					pa->genome_end_padded+=length;	
				}
				for (j=0; j<length; j++) {
					if (current_genome_index<0 || current_genome_index>=pa->genome_length) {
						pa->pretty_genome_string[current_pretty_genome_string_index++]='-';
						//pa->genome_string[current_genome_string_index++]='-';
					} else {
						pa->pretty_genome_string[current_pretty_genome_string_index++]=pa->genome_sequence[current_genome_index];
						//pa->genome_string[current_genome_string_index++]=pa->genome_sequence[current_genome_index];
					}
					current_genome_index++;
				}
				break;
			case 'P':
			case 'I':
				start=false;
				for (j=0; j<length; j++) {
					pa->pretty_genome_string[current_pretty_genome_string_index++]='-';
				}
				break;
			default:
				fprintf(stderr,"Error in converting to pretty genome,%s\n",pa->cigar);
				exit(1);
				
						
		}
	}
	pa->genome_end_padded--;
	pa->genome_end_unpadded--;
	pa->pretty_genome_string[current_pretty_genome_string_index]='\0';
	//pa->genome_string[current_genome_string_index]='\0';
	if (pa->reverse){
		reverse_complement(pa->pretty_genome_string);
		//reverse_complement(pa->genome_string);
	}
	return;	
}


void remove_gaps(char * s) {
	size_t length=strlen(s);
	int i; 
	int gaps=0;
	for (i=0; i<length; i++) {
		if (s[i]=='-') {
			gaps++;
		} else {
			s[i-gaps]=s[i];
		}
	}
	s[i]='\0';
}


//char* pretty_ls_seq(char* seq, char* cigar, int32_t num_cigar, int32_t strand, int32_t * read_start, int32_t * read_end) {
void pretty_ls(pretty * pa) {
	if (pa->pretty_length==0) {
		pa->pretty_length=strlen(pa->read_string);
		pa->pretty_read_string = (char*)malloc(sizeof(char)*(pa->pretty_length+1));
		if (pa->pretty_read_string==NULL) {
			fprintf(stderr,"Malloc failed for ret in pretty_ls_seq!\n");
			exit(1);
		}
		strcpy(pa->pretty_read_string,pa->read_string);	
		pa->clipped_read_start=1;
		pa->clipped_read_end=pa->pretty_length;
		pa->clipped_read_length=pa->pretty_length;
		pa->unclipped_read_length=pa->clipped_read_length;
		return;
	}
	assert(pa->read_string!=NULL);
	assert(pa->pretty_read_string==NULL);
	pa->pretty_read_string = (char*)malloc(sizeof(char)*(pa->pretty_length+1));
	if (pa->pretty_read_string==NULL) {
		fprintf(stderr,"Malloc failed for ret in pretty_ls_seq!\n");
		exit(1);
	}

	int32_t current_pretty_read_string_index=0;
	int32_t current_read_string_index=0;
	pa->unclipped_read_length=pa->clipped_read_length;
	//printf("|%s| len %d, read_length %d\n",pa->read_string,strlen(pa->read_string),pa->clipped_read_length);
	//printf("start %d\nend %d\n",pa->clipped_read_start,pa->clipped_read_end);
	//printf("unclipped read length %d\n",pa->unclipped_read_length);
	int32_t j,i;
	int32_t deletions=0; int32_t insertions=0;
	//fprintf(stderr,"%s\n",pa->read_name);
	//for (i=0; i<pa->num_cigar; i++) {
	//	fprintf(stderr,"%d%c",pa->cigar_lengths[i],pa->cigar_ops[i]);
	//}
	//fprintf(stderr,"\n");
	char * source_string = strdup(pa->has_cs_edit_string ? pa->cs_edit_string : pa->read_string);
	if (pa->has_cs_edit_string && pa->reverse) {
		reverse_complement(source_string);
	}
	remove_gaps(source_string);
	for (i=0; i<pa->num_cigar; i++) {
		int32_t length=pa->cigar_lengths[i];
		switch (pa->cigar_ops[i]) {
			case 'S':
				current_read_string_index+=length;
				for (j=0; j<length; j++) {
					pa->pretty_read_string[current_pretty_read_string_index++]='-';
				}
				break;
			case 'H':
				pa->unclipped_read_length+=length;
				for (j=0; j<length; j++) {
					pa->pretty_read_string[current_pretty_read_string_index++]='-';
				}
				break;
			case 'D':
			case 'P':
				deletions+=length;
				for (j=0; j<length; j++) {
					pa->pretty_read_string[current_pretty_read_string_index++]='-';
				}
				break;
			case 'I':
				insertions+=length;
			case 'M':
				{
				//char * source_string=pa->has_cs_edit_string ? pa->cs_edit_string : pa->read_string;
				//char * source_string= pa->read_string;
				for (j=0; j<length; j++) {
					//if (current_read_string_index>=pa->clipped_read_length+deletions) {
					//	fprintf(stderr,"A5 offending read %s\n",pa->read_name);
					//	exit(1);
					//
					pa->pretty_read_string[current_pretty_read_string_index++]=source_string[current_read_string_index++];
				}
				}
				break;
			case 'N':
				for (j=0; j<length; j++) {
					pa->pretty_read_string[current_pretty_read_string_index++]='.';
				}
				break;
			default:
				fprintf(stderr,"Error in parsing cigar string in pretty_ls, %s\n",pa->cigar);
				exit(1);
				
						
		}
	}
	free(source_string);
	pa->pretty_read_string[current_pretty_read_string_index]='\0';
	//printf("ls : |%s| %d\n",pa->pretty_read_string,current_pretty_read_string_index);
	pa->clipped_read_start=1;
	pa->pretty_clipped_read_start=1;
	for (i=0; i<pa->num_cigar; i++) {
		if (pa->cigar_ops[i]=='S' || pa->cigar_ops[i]=='H') {
			pa->clipped_read_start+=pa->cigar_lengths[i];
			pa->pretty_clipped_read_start+=pa->cigar_lengths[i];
		} else {
			break;
		}
	}
	pa->clipped_read_end=pa->unclipped_read_length;
	pa->pretty_clipped_read_end=pa->unclipped_read_length+deletions;
	//printf("unclipped_read len %d\n",pa->unclipped_read_length);
	for (i=0; i<pa->num_cigar; i++) {
		if (pa->cigar_ops[pa->num_cigar-i-1]=='S' || pa->cigar_ops[pa->num_cigar-i-1]=='H') {
			pa->clipped_read_end-=pa->cigar_lengths[pa->num_cigar-i-1];
			pa->pretty_clipped_read_end-=pa->cigar_lengths[pa->num_cigar-i-1];
		} else {
			break;
		}
	}
	//printf("ls: %d and end %d\n",pa->clipped_read_start,pa->clipped_read_end);
	//printf("sd: %d\ned %d\n",start_deletions,end_deletions);	
	if (pa->reverse){
		//printf("about to flip/, %d, %d vs %d\n",pa->unclipped_read_length,
		//		pa->clipped_read_start,pa->clipped_read_end);
		int32_t tmp=pa->clipped_read_start;
		pa->clipped_read_start=pa->unclipped_read_length-pa->clipped_read_end+1;
		pa->clipped_read_end=pa->unclipped_read_length-tmp+1;
		int32_t pretty_tmp=pa->pretty_clipped_read_start;
		pa->pretty_clipped_read_start=deletions+pa->unclipped_read_length-pa->pretty_clipped_read_end+1;
		pa->pretty_clipped_read_end=deletions+pa->unclipped_read_length-pretty_tmp+1;
		reverse_complement(pa->pretty_read_string);
	}	
	//printf("ls : |%s| %d\n",pa->pretty_read_string,current_pretty_read_string_index);
	//printf("ls: readstart %d and readend %d\n",pa->clipped_read_start,pa->clipped_read_end);
	//printf("pr ls: readstart %d and readend %d\n",pa->pretty_clipped_read_start,pa->pretty_clipped_read_end);
	return;	
}

//char* pretty_cs_seq(char* seq, char* cigar, int32_t num_cigar, int32_t strand, int32_t * read_start, int32_t * read_end) {

void pretty_cs_qualities(pretty * pa) {
	if (pa->pretty_length==0) {
		return;
	}
	if (pa->cs_qualities==NULL) {
		fprintf(stderr,"pa->cs_qualities is NULL!\n");
		exit(1);
	}
	if (pa->cs_qualities[0]=='*' && pa->cs_qualities[1]=='\0'){ 
		pa->pretty_cs_qualities=strdup(pa->cs_qualities);
		return;
	}
	assert(pa->unclipped_read_length!=0);
	assert(pa->unclipped_read_length=strlen(pa->cs_qualities));
	assert(pa->pretty_cs_qualities==NULL);
	pa->pretty_cs_qualities = (char*)malloc(sizeof(char)*(pa->pretty_length+1));
	if (pa->pretty_cs_qualities==NULL) {
		fprintf(stderr,"Malloc failed for pretty_cs_qualities in pretty_qual_seq!\n");
		exit(1);
	}
	int32_t current_pretty_cs_qualities_index=0;
	int32_t current_cs_qualities_index=0;
	int32_t i,j;
	for (i=0; i<pa->num_cigar; i++) {
		int32_t index;
		if (!pa->reverse) {
			index=i;
		} else {
			index=pa->num_cigar-1-i;
		}
		int32_t length=pa->cigar_lengths[index];
		switch (pa->cigar_ops[index]) {
			case 'S':
			case 'I':
			case 'M':
				for(j=0; j<length; j++) {
					pa->pretty_cs_qualities[current_pretty_cs_qualities_index++]=pa->cs_qualities[current_cs_qualities_index++];
				}
				break;
			case 'N':
			case 'H':
			case 'D':
			case 'P':
				for (j=0; j<length; j++) {
					pa->pretty_cs_qualities[current_pretty_cs_qualities_index++]=' ';
				}
				break;
			default:
				fprintf(stderr,"Error parsing cigar string in pretty_qualities, %s\n",pa->cigar);
				exit(1);
		}
	}
	pa->pretty_cs_qualities[current_pretty_cs_qualities_index]='\0'; 
	return;	
}
void pretty_qualities(pretty * pa) {
	if (pa->pretty_length==0) {
		return;
	}
	assert(pa->read_qualities!=NULL);
	if (pa->has_cs_qualities) {
		pretty_cs_qualities(pa);
	}
	if (pa->read_qualities[0]=='*' && pa->read_qualities[1]=='\0'){ 
		pa->pretty_read_qualities=strdup(pa->read_qualities);
		return;
	}
	assert(pa->unclipped_read_length!=0);
	assert(pa->unclipped_read_length=strlen(pa->read_qualities));
	assert(pa->pretty_read_qualities==NULL);
	pa->pretty_read_qualities = (char*)malloc(sizeof(char)*(pa->pretty_length+1));
	if (pa->pretty_read_qualities==NULL) {
		fprintf(stderr,"Malloc failed for pretty_read_qualities in pretty_qual_seq!\n");
		exit(1);
	}
	int32_t current_pretty_read_qualities_index=0;
	int32_t current_read_qualities_index=0;
	int32_t i,j;
	for (i=0; i<pa->num_cigar; i++) {
		int32_t index;
		if (!pa->reverse) {
			index=i;
		} else {
			index=pa->num_cigar-1-i;
		}
		int32_t length=pa->cigar_lengths[index];
		switch (pa->cigar_ops[index]) {
			case 'S':
			case 'I':
			case 'M':
				for(j=0; j<length; j++) {
					pa->pretty_read_qualities[current_pretty_read_qualities_index++]=pa->read_qualities[current_read_qualities_index++];
				}
				break;
			case 'N':
			case 'D':
			case 'H':
			case 'P':
				for (j=0; j<length; j++) {
					pa->pretty_read_qualities[current_pretty_read_qualities_index++]=' ';
				}
				break;
			default:
				fprintf(stderr,"Error parsing cigar string in pretty_qualities, %s\n",pa->cigar);
				exit(1);
		}
	}
	pa->pretty_read_qualities[current_pretty_read_qualities_index]='\0'; 
	return;	
}

void pretty_cs(pretty * pa) {
	if (pa->pretty_length==0) {
		return;
	}
	assert(pa->cs_string!=NULL);
	assert(pa->pretty_cs_string==NULL);
	pa->pretty_cs_string = (char*)malloc(sizeof(char)*(pa->pretty_length+1+1)); //extra letter
	if (pa->pretty_cs_string==NULL) {
		fprintf(stderr,"Malloc failed for pretty_cs_string in pretty_cs_seq!\n");
		exit(1);
	}
	//copy over the letter at the start of cs read
	int32_t current_pretty_cs_string_index=0;
	int32_t current_cs_string_index=0;
	pa->pretty_cs_string[current_pretty_cs_string_index++]=pa->cs_string[current_cs_string_index++];
	int32_t i,j;
	for (i=0; i<pa->num_cigar; i++) {
		int32_t index;
		if (!pa->reverse) {
			index=i;
		} else {
			index=pa->num_cigar-1-i;
		}
		int32_t length=pa->cigar_lengths[index];
		switch (pa->cigar_ops[index]) {
			case 'S':
			case 'H':
			case 'I':
			case 'M':
				for(j=0; j<length; j++) {
					pa->pretty_cs_string[current_pretty_cs_string_index++]=pa->cs_string[current_cs_string_index++];
				}
				break;
			case 'D':
			case 'P':
				for (j=0; j<length; j++) {
					pa->pretty_cs_string[current_pretty_cs_string_index++]='-';
				}
				break;
			case 'N':
				for (j=0; j<length; j++) {
					pa->pretty_cs_string[current_pretty_cs_string_index++]='.';
				}
				break;
			default:
				fprintf(stderr,"Error parsing cigar string in pretty_cs, %s\n",pa->cigar);
				exit(1);
		}
	}
	assert(current_pretty_cs_string_index<=pa->pretty_length+1);
	pa->pretty_cs_string[current_pretty_cs_string_index]='\0'; 
	//printf("cs : |%s| %d\n",pa->pretty_cs_string,current_pretty_cs_string_index);
	//printf("cs: %d and end %d\n",pa->clipped_read_start,pa->clipped_read_end);
	return;	
}


char * pretty_reformat_cigar_string(pretty * pa) {
	int32_t cigar_string_length=pa->num_cigar+1;
	int32_t i;
	for (i=0; i<pa->num_cigar; i++) {
		cigar_string_length+=(pa->cigar_lengths[i]/10) +1;
	}
	char * temp_cigar_string = (char*)malloc(sizeof(char)*cigar_string_length);
	if (temp_cigar_string==NULL) {
		fprintf(stderr,"pretty_reformat_cigar_string : temp_cigar_string, failed to allocate memory!\n");
		exit(1);
	}
	char * ptr = temp_cigar_string;
	for (i=0; i<pa->num_cigar; i++) {
		ptr+=sprintf(ptr,"%d%c",pa->cigar_lengths[i],pa->cigar_ops[i]);
	}
	*ptr='\0';
	return temp_cigar_string;
} 


void shift_string(char* s, int32_t d) {
	int32_t i;
	for (i=0; i<(int32_t)strlen(s)-d; i++){
		s[i]=s[i+d];
	}
	s[i]='\0';
	return;
}


//return a linked list of pretty structures
//each of which represents one of the continuous exons
//in the alignment
pretty * pretty_sub_prettys(pretty * parent) {
	pretty * head=NULL; pretty * tail=NULL;
	int32_t i; int32_t last=1;
	for (i=0; i<parent->num_cigar; i++) {
		if (parent->cigar_ops[i]=='N') {
			if (last!=i+1) {
				pretty * new_node = pretty_sub_pretty(parent,last,i);
				if (tail==NULL) {
					head=new_node;
					tail=new_node;
				} else {
					tail->next=(struct pretty *)new_node;
					tail=new_node;
				}
			}
			last=i+2;
		}
	}
	if (last!=1 && last!=i+1) {
		pretty * new_node = pretty_sub_pretty(parent,last,i);
		if (tail==NULL) {
			head=new_node;
			tail=new_node;
		} else {
			tail->next=(struct pretty * )new_node;
			tail=new_node;
		}
	}
	/*pretty * temp = head;
	while (temp!=NULL) {
		pretty_print_sam(stdout,temp);
		temp=(pretty*)temp->next;
	}*/
	if (head==NULL) {
		head=parent;
	}
	return head;
}

pretty * pretty_sub_pretty(pretty * parent, int32_t from, int32_t to) {
	pretty * child = pretty_copy_base(parent);
	//allocate space for cigar array
	int32_t cigar_size = to - from + 1 + 2;
	child->cigar_ops = (char*)malloc(sizeof(char)*cigar_size);
	child->cigar_lengths = (uint32_t*)malloc(sizeof(uint32_t)*cigar_size);
	int32_t cigar_used=0;
	//Get the front out
	int32_t i; int32_t read_shift=0; int32_t genome_shift=0; int32_t hard_clipped=0;
	for (i=0; i<from-1; i++) {
		int32_t length=parent->cigar_lengths[i];
		switch (parent->cigar_ops[i]) {
			case 'M':
				genome_shift+=length;
				read_shift+=length;
				break;			
			case 'I':
				read_shift+=length;
				break;	
			case 'D':
			case 'N':
				genome_shift+=length;
				break;
			case 'P':
				break;	
			case 'S':
				read_shift+=length;
				break;
			case 'H':
				hard_clipped+=length;
				break;
			default:
				fprintf(stderr,"pretty_sub_pretty failed!\n");
				exit(1);
		}
	}
	if (hard_clipped>0 || read_shift>0) {
		child->genome_start_unpadded+=genome_shift;
		shift_string(child->read_string,read_shift);
		if (child->read_qualities[0]!='*' || child->read_qualities[1]!='\0') {
			shift_string(child->read_qualities,read_shift);
		}
		child->cigar_ops[cigar_used]='H';
		child->cigar_lengths[cigar_used++]=hard_clipped+read_shift;	
	}
	//do the middle
	for (i=from-1; i<to; i++) {
		child->cigar_ops[cigar_used]=parent->cigar_ops[i];
		child->cigar_lengths[cigar_used++]=parent->cigar_lengths[i];
	}
	//do the end	
	read_shift=0; genome_shift=0; hard_clipped=0; 
	for (i=to; i<parent->num_cigar; i++) {
		int32_t length=parent->cigar_lengths[i];
		switch (parent->cigar_ops[i]) {
			case 'M':
				genome_shift+=length;
				read_shift+=length;
				break;			
			case 'I':
				read_shift+=length;
				break;	
			case 'D':
			case 'N':
				genome_shift+=length;
				break;
			case 'P':
				break;	
			case 'S':
				read_shift+=length;
				break;
			case 'H':
				hard_clipped+=length;
				break;
			default:
				fprintf(stderr,"pretty_sub_pretty failed!\n");
				exit(1);
		}
	}
	if (hard_clipped>0 || read_shift>0) {
		child->read_string[strlen(child->read_string)-read_shift]='\0';
		if (child->read_qualities[0]!='*' || child->read_qualities[1]!='\0') {
			child->read_qualities[strlen(child->read_qualities)-read_shift]='\0';
		}
		child->cigar_ops[cigar_used]='H';
		child->cigar_lengths[cigar_used++]=read_shift+hard_clipped;
	}	
	child->num_cigar=cigar_used;
	child->cigar=pretty_reformat_cigar_string(child);
	//pretty_ls(child);
	//pretty_genome(child);
	//if (child->colour_space) {
	//	pretty_cs(child);
	//}
	//pretty_qualities(child);
	return child;
}



void pretty_print_sam_force_ls(FILE * f, pretty * pa) {
	pretty pa_dup=*pa;
	pa_dup.colour_space=false;
	pretty_print_sam(f,&pa_dup);	
}

void pretty_print_sam(FILE * f, pretty * pa) {
	fprintf(f,"%s\t%d\t%s\t%d\t%d\t%s\t%s\t%d\t%d\t%s\t%s",
		pa->read_name,
		pa->flags,
		pa->reference_name,
		pa->genome_start_unpadded,
		pa->mapq,
		pa->cigar,
		pa->mate_reference_name,
		pa->mate_genome_start_unpadded,
		pa->isize,
		pa->read_string,
		pa->read_qualities);
	if (pa->has_score) {
		fprintf(f,"\tAS:i:%d",pa->score);
	}	 
	if (pa->has_edit_distance) {
		fprintf(f,"\tNM:i:%d",pa->edit_distance);
	}	 
	if (pa->colour_space) {
		if (pa->has_cs_qualities) {
			fprintf(f,"\tCQ:Z:%s",pa->cs_qualities);
		}
		fprintf(f,"\tCS:Z:%s",pa->cs_string);
		if (pa->has_cs_mismatches) {
			fprintf(f,"\tCM:i:%d",pa->cs_mismatches);
		}
		if (pa->has_cs_edit_string) {
			fprintf(f,"\tXX:Z:%s",pa->cs_edit_string);
		}
	}
	if (pa->has_read_group) {
		fprintf(f,"\tRG:Z:%s",pa->read_group);
	}
	fprintf(f,"\n");
}

static inline int32_t field2int32_t(char* s) {
	return (s[0]<<8) + s[1]; 
}


static inline void not_sam(char * next_tab) {
	if (next_tab==NULL) {
		fprintf(stderr,"NOT SAM!\n");
		exit(1);
	}
}

static inline void switch_and_fill(int32_t k, char * data, pretty * pa) {
	switch (k) {
		case ('N'<<8)+'M':
			pa->edit_distance=atoi(data);			
			pa->has_edit_distance=true;
			break;
		case ('C'<<8)+'M':
			pa->cs_mismatches=atoi(data);
			pa->has_cs_mismatches=true;
			break;
		case ('R'<<8)+'2':
			pa->r2=data;
			pa->has_r2=true;
			break;
		case ('X'<<8)+'X':
			pa->cs_edit_string=data;
			pa->has_cs_edit_string=true;
			break;
		case ('R'<<8)+'G':
			pa->read_group=data;
			pa->has_read_group=true;
			break;
		case ('C'<<8)+'S':
			pa->colour_space=true;
			pa->cs_string=data;
			pa->has_cs_string=true;
			break;
		case ('C'<<8)+'Q':
			pa->cs_qualities=data;
			pa->has_cs_qualities=true;
			break;
		case ('I'<<8)+'H':
			pa->has_ih=true;
			pa->ih=atoi(data);
			break;
		case ('H'<<8)+'I':
			pa->has_hi=true;
			pa->hi=atoi(data);
			break;
		case ('H'<<8)+'0':
			pa->has_h0=true;
			pa->h0=atoi(data);
			break;
		case ('H'<<8)+'1':
			pa->has_h1=true;
			pa->h1=atoi(data);
			break;
		case ('H'<<8)+'2':
			pa->has_h2=true;
			pa->h2=atoi(data);
			break;
		case ('Z'<<8)+'0':
			pa->has_zs|=HAS_Z0;
			pa->z[0]=inv_tnlog(atoi(data));
			break;
		case ('Z'<<8)+'1':
			pa->has_zs|=HAS_Z1;	
			pa->z[1]=inv_tnlog(atoi(data));
			break;
		case ('Z'<<8)+'2':
			pa->has_zs|=HAS_Z2;
			pa->z[2]=inv_tnlog(atoi(data));
			break;
		case ('Z'<<8)+'3':
			pa->has_zs|=HAS_Z3;	
			pa->z[3]=inv_tnlog(atoi(data));
			break;
		case ('Z'<<8)+'4':
			pa->has_zs|=HAS_Z4;	
			pa->z[4]=inv_tnlog(atoi(data));
			break;
		case ('Z'<<8)+'5':
			pa->has_zs|=HAS_Z5;	
			pa->z[5]=inv_tnlog(atoi(data));
			break;
		case ('Z'<<8)+'6':
			pa->has_zs|=HAS_Z6;	
			pa->z[6]=inv_tnlog(atoi(data));
			break;
		default: 
			break;
	}
}


void pretty_from_aux_inplace(pretty * pa) {
	assert(pa!=NULL);
	if (pa->aux==NULL) {
		return;
	}
	char * start_of_string = pa->aux;
	//get the read name
	char * next_tab = (char*)memchr(start_of_string,'\t',pa->sam_string_length-(start_of_string-pa->sam_string)+1);
	while (next_tab!=NULL) {
		not_sam(next_tab); *next_tab='\0';
		int32_t k=field2int32_t(start_of_string);
		char * data=start_of_string+5;
		switch_and_fill(k,data,pa);
		start_of_string=++next_tab;
		next_tab = (char*)memchr(start_of_string,'\t',pa->sam_string_length-(start_of_string-pa->sam_string)+1);
	}
	//don't forget to process the last one
	int32_t k=field2int32_t(start_of_string);
	char * data=start_of_string+5;
	switch_and_fill(k,data,pa);
	pa->aux=NULL;
} 


pretty * pretty_from_string_inplace(char * sam_string,size_t length_of_string,pretty * pa) {
	memset(pa,0,sizeof(pretty));
	pa->sam_string=sam_string;
	pa->sam_string_length=length_of_string;
	char * start_of_string = sam_string;
	//get the read name
	char * next_tab = (char*)memchr(start_of_string,'\t',length_of_string);
	not_sam(next_tab); *next_tab='\0';
	pa->read_name=start_of_string;
	start_of_string=++next_tab;
	//get the flags
	next_tab=(char*)memchr(start_of_string,'\t',length_of_string-(next_tab-sam_string));
	not_sam(next_tab); *next_tab='\0';
	pa->flags=atoi(start_of_string);
	start_of_string=++next_tab;
	//get the reference name
	next_tab=(char*)memchr(start_of_string,'\t',length_of_string-(next_tab-sam_string));
	not_sam(next_tab); *next_tab='\0';
	pa->reference_name=start_of_string;
	start_of_string=++next_tab;
	//get the position
	next_tab=(char*)memchr(start_of_string,'\t',length_of_string-(next_tab-sam_string));
	not_sam(next_tab); *next_tab='\0';
	pa->genome_start_unpadded=atoi(start_of_string);
	start_of_string=++next_tab;
	//get the mapq
	next_tab=(char*)memchr(start_of_string,'\t',length_of_string-(next_tab-sam_string));
	not_sam(next_tab); *next_tab='\0';
	pa->mapq=atoi(start_of_string);
	start_of_string=++next_tab;
	//get the cigar strign
	next_tab=(char*)memchr(start_of_string,'\t',length_of_string-(next_tab-sam_string));
	not_sam(next_tab); *next_tab='\0';
	pa->cigar=start_of_string;
	start_of_string=++next_tab;
	//get the mate reference name
	next_tab=(char*)memchr(start_of_string,'\t',length_of_string-(next_tab-sam_string));
	not_sam(next_tab); *next_tab='\0';
	pa->mate_reference_name=start_of_string;
	start_of_string=++next_tab;
	//get the mate position
	next_tab=(char*)memchr(start_of_string,'\t',length_of_string-(next_tab-sam_string));
	not_sam(next_tab); *next_tab='\0';
	pa->mate_genome_start_unpadded=atoi(start_of_string);
	start_of_string=++next_tab;
	//get the isize
	next_tab=(char*)memchr(start_of_string,'\t',length_of_string-(next_tab-sam_string));
	not_sam(next_tab); *next_tab='\0';
	pa->isize=atoi(start_of_string);
	start_of_string=++next_tab;
	//get the sequence
	next_tab=(char*)memchr(start_of_string,'\t',length_of_string-(next_tab-sam_string));
	not_sam(next_tab); *next_tab='\0';
	pa->read_string=start_of_string;
	start_of_string=++next_tab;
	//get the qualities
	next_tab=(char*)memchr(start_of_string,'\t',length_of_string-(next_tab-sam_string));
	if (next_tab!=NULL) {
		*next_tab='\0';
	}
	pa->read_qualities=start_of_string;
	//get the score
	pa->aux=NULL;
	if (next_tab!=NULL && next_tab+6<length_of_string+sam_string) {
		start_of_string=++next_tab;
		next_tab=(char*)memchr(start_of_string,'\t',length_of_string-(next_tab-sam_string));
		if (start_of_string[0]=='A' && start_of_string[1]=='S') {
			if (next_tab!=NULL) {
				*next_tab='\0';
			}
			pa->score=atoi(start_of_string+5);
			pa->has_score=true;
			//TODO ERROR NEED TO FIX THIS!!!!!
			assert(pa->has_zs==0);
			if (next_tab!=NULL && (next_tab-sam_string)+6<length_of_string) {
				int i;
				for (i=0; next_tab!=NULL && i<SAM2PRETTY_NUM_ZS; i++) {
					if (next_tab[1]!='Z') {
						break;
					}
					int z_index=next_tab[2]-48;
					if ((pa->has_zs&(1<<z_index))!=0) {
						fprintf(stderr,"A fatal error occured parsing this sam line!\n");
						exit(1);
					}
					pa->has_zs|=(1<<z_index);
					
					start_of_string=++next_tab;
					
					if (length_of_string-(next_tab-sam_string)<6) {
						fprintf(stderr,"There has been an error in parsing ZX fields\n");
						exit(1);
					}
					next_tab=(char*)memchr(start_of_string,'\t',length_of_string-(next_tab-sam_string));
					if (next_tab!=NULL) {
						*next_tab='\0';											
					}
					//next field should be null terminated
					if (start_of_string[0]!='Z' && start_of_string[1]!=48+i) {
						fprintf(stderr,"there has been an error parsing ZX fields! there must be 5 of them!\n");
						exit(1);
					}
					pa->z[z_index]=inv_tnlog(atoi(start_of_string+5));	
				}
				if (next_tab!=NULL) {
					next_tab++;
				}
				pa->aux=next_tab;
			} else if (next_tab!=NULL) {
				pa->aux=next_tab;
			} else {
				pa->aux=NULL;
			}
		} else {
			pa->score=0;
			pa->has_score=false;
			pa->aux=start_of_string;
		}	
	}
	//done
	//0x0001 the read is paired in sequencing, no matter whether it is mapped in a pair
	pa->paired_sequencing = ((pa->flags & 0x0001) != 0) ? true : false ;
	//0x0002 the read is mapped in a proper pair (depends on the protocol, normally inferred during alignment) 1
	pa->proper_pair = ((pa->flags & 0x0002) != 0) ? true : false ;
	//0x0004 the query sequence itself is unmapped
	pa->mapped = ((pa->flags & 0x0004) != 0) ? false : true ;
	//0x0008 the mate is unmapped 1
	pa->mp_mapped = ((pa->flags & 0x0008) != 0) ? false : true ;
	//0x0010 strand of the query (0 for forward; 1 for reverse strand)
	pa->reverse = ((pa->flags & 0x0010) != 0) ? true : false ;
	//0x0020 strand of the mate 1
	pa->mp_reverse = ((pa->flags & 0x0020) != 0) ? true : false ;
	//0x0040 the read is the first read in a pair 1,2
	pa->first_in_pair = ((pa->flags & 0x0040) != 0) ? true : false ;
	//0x0080 the read is the second read in a pair 1,2
	pa->second_in_pair = ((pa->flags & 0x0080) != 0) ? true : false ;
	//0x0100 the alignment is not primary (a read having split hits may have multiple primary alignment records)
	pa->primary_alignment = ((pa->flags & 0x0100) != 0) ? true : false ;
	//0x0200 the read fails platform/vendor quality checks
	pa->platform_quality_fail = ((pa->flags & 0x0200) != 0) ? true : false ;
	//0x0400 the read is either a PCR duplicate or an optical duplicate
	pa->pcr_duplicate = ((pa->flags & 0x0400) != 0) ? true : false ;
	
	return pa;	
}

void fill_cigar_len(pretty * pa) { 
	int32_t i,cigar_len=strlen(pa->cigar);
	pa->num_cigar=0;
	for (i=0; i<cigar_len; i++) {
		if (pa->cigar[i]>=65) {
			pa->num_cigar++;	
		}
	}
}

pretty * pretty_from_string(char* sam_string) {
		if (sam_string[strlen(sam_string)-1]=='\n') {
			sam_string[strlen(sam_string)-1]='\0';
		}
		if (sam_string[0]=='\0') {
			return NULL;
		}
		pretty * pa = pretty_new();
		pa->sam_string_length=strlen(sam_string);
		pa->flags=0;

		char * save_ptr;
		int32_t parsed_tokens=0;
		char * tok;
		char * str = sam_string;
		while (  (parsed_tokens<11) && ( (tok = strtok_r(str, "\t", &save_ptr)) != NULL)) {
			switch (parsed_tokens) {
				case 0:
					pa->read_name=(char*)malloc(sizeof(char)*(strlen(tok)+1));
					if (pa->read_name==NULL) {
						fprintf(stderr,"Failed to allocate read read name string!\n");
						exit(1);
					}
					strcpy(pa->read_name,tok);
					pa->read_name_length=strlen(pa->read_name);
					break;
				case 1:
					pa->flags=atoi(tok);
					break;
				case 2:
					pa->reference_name=(char*)malloc(sizeof(char)*(strlen(tok)+1));
					if (pa->reference_name==NULL) {
						fprintf(stderr,"Failed to allocate reference name string!\n");
						exit(1);
					}
					strcpy(pa->reference_name,tok);
					break;
				case 3:
					pa->genome_start_unpadded=atoi(tok);
					/*if (pa->genome_start_unpadded==0) {
						pretty_free(pa);
						return NULL;
					}*/
					break;
				case 4:
					pa->mapq=atoi(tok);
					break;
				case 5:
					pa->cigar=(char*)malloc(sizeof(char)*(strlen(tok)+1));
					if (pa->cigar==NULL) {
						fprintf(stderr,"Failed to allocate cigar string!\n");
						exit(1);
					}
					strcpy(pa->cigar,tok);
					{
						fill_cigar_len(pa);
						/*int32_t i,cigar_len=strlen(pa->cigar);
						pa->num_cigar=0;
						for (i=0; i<cigar_len; i++) {
							if (pa->cigar[i]>=65) {
								pa->num_cigar++;	
							}
						}*/
					}
					break;
				case 6:
					pa->mate_reference_name=(char*)malloc(sizeof(char)*(strlen(tok)+1));
					if (pa->mate_reference_name==NULL) {
						fprintf(stderr,"Failed to allocate reference name string!\n");
						exit(1);
					}
					strcpy(pa->mate_reference_name,tok);
					break;
				case 7:
					pa->mate_genome_start_unpadded=atoi(tok);
					break;
				case 8:
					pa->isize=atoi(tok);
					break;
				case 9:
					pa->read_string=(char*)malloc(sizeof(char)*(strlen(tok)+1));
					if (pa->read_string==NULL) {
						fprintf(stderr,"Failed to allocate read string!\n");
						exit(1);
					}
					strcpy(pa->read_string,tok);
					pa->clipped_read_length=strlen(pa->read_string);
					//printf("Read in length is %d\n",pa->clipped_read_length);
					break;
				case 10:
					pa->read_qualities=(char*)malloc(sizeof(char)*(strlen(tok)+1));
					if (pa->read_qualities==NULL) {
						fprintf(stderr,"Failed to allocate read qualities string!\n");
						exit(1);
					}
					strcpy(pa->read_qualities,tok);
					break;
				default:
					break;
			}
			str=NULL;
			parsed_tokens++;			
		}
		
		assert(pa->read_qualities!=NULL);
		assert(parsed_tokens>=10);
		pa->edit_distance=-1;
		pa->skipped=-1;
		pa->colour_space=false;
		while ((tok = strtok_r(str, "\t", &save_ptr)) != NULL) {
			assert(tok[2]==':');
			assert(tok[4]==':');
			int32_t switch_on=field2int32_t(tok);
			tok+=5;
			switch (switch_on) {
				case ('A'<<8)+'S':
					pa->has_score=true;
					pa->score=atoi(tok);
					break;
				case ('N'<<8)+'M':
					pa->edit_distance=atoi(tok);			
					pa->has_edit_distance=true;
					break;
				case ('C'<<8)+'M':
					pa->cs_mismatches=atoi(tok);
					pa->has_cs_mismatches=true;
					break;
				case ('R'<<8)+'2':
					pa->r2=(char*)malloc(sizeof(char)*strlen(tok)+1);
					if (pa->r2==NULL) {
						fprintf(stderr,"Failed to allocate cs_edit_string!\n");
						exit(1);
					}
					strcpy(pa->r2,tok);
					pa->has_r2=true;
					break;
				case ('X'<<8)+'X':
					pa->cs_edit_string=(char*)malloc(sizeof(char)*strlen(tok)+1);
					if (pa->cs_edit_string==NULL) {
						fprintf(stderr,"Failed to allocate cs_edit_string!\n");
						exit(1);
					}
					strcpy(pa->cs_edit_string,tok);
					pa->has_cs_edit_string=true;
					break;
				case ('R'<<8)+'G':
					pa->read_group=(char*)malloc(sizeof(char)*strlen(tok)+1);
					if (pa->read_group==NULL) {
						fprintf(stderr,"Failed to allocate read_group!\n");
						exit(1);
					}
					strcpy(pa->read_group,tok);
					pa->has_read_group=true;
					break;
				case ('C'<<8)+'S':
					pa->colour_space=true;
					pa->cs_string=(char*)malloc(sizeof(char)*(strlen(tok)+1));
					if (pa->cs_string==NULL) {
						fprintf(stderr,"Failed to allocate cs_read_string!\n");
						exit(1);
					}
					strcpy(pa->cs_string,tok);
					pa->has_cs_string=true;
					break;
				case ('C'<<8)+'Q':
					/*if (pa->read_qualities==NULL || pa->read_qualities[0]!='*') {
						fprintf(stderr,"cs qualities unsupported in both fields!\n");
						exit(1);
					}*/
					//free(pa->read_qualities);
					pa->cs_qualities=(char*)malloc(sizeof(char)*(strlen(tok)+1));
					if (pa->cs_qualities==NULL) {
						fprintf(stderr,"Failed to allocate cs_qualities!\n");
						exit(1);
					}
					strcpy(pa->cs_qualities,tok);
					pa->has_cs_qualities=true;
					break;
				case ('I'<<8)+'H':
					pa->has_ih=true;
					pa->ih=atoi(tok);
					break;
				case ('H'<<8)+'I':
					pa->has_hi=true;
					pa->hi=atoi(tok);
					break;
				case ('H'<<8)+'0':
					pa->has_h0=true;
					pa->h0=atoi(tok);
					break;
				case ('H'<<8)+'1':
					pa->has_h1=true;
					pa->h1=atoi(tok);
					break;
				case ('H'<<8)+'2':
					pa->has_h2=true;
					pa->h2=atoi(tok);
					break;
				default:
					break;
			}
		}
		//Get all the flags	
		//0x0001 the read is paired in sequencing, no matter whether it is mapped in a pair
		pa->paired_sequencing = ((pa->flags & 0x0001) != 0) ? true : false ;
		//0x0002 the read is mapped in a proper pair (depends on the protocol, normally inferred during alignment) 1
		pa->proper_pair = ((pa->flags & 0x0002) != 0) ? true : false ;
		//0x0004 the query sequence itself is unmapped
		pa->mapped = ((pa->flags & 0x0004) != 0) ? false : true ;
		//0x0008 the mate is unmapped 1
		pa->mp_mapped = ((pa->flags & 0x0008) != 0) ? false : true ;
		//0x0010 strand of the query (0 for forward; 1 for reverse strand)
		pa->reverse = ((pa->flags & 0x0010) != 0) ? true : false ;
		//0x0020 strand of the mate 1
		pa->mp_reverse = ((pa->flags & 0x0020) != 0) ? true : false ;
		//0x0040 the read is the first read in a pair 1,2
		pa->first_in_pair = ((pa->flags & 0x0040) != 0) ? true : false ;
		//0x0080 the read is the second read in a pair 1,2
		pa->second_in_pair = ((pa->flags & 0x0080) != 0) ? true : false ;
		//0x0100 the alignment is not primary (a read having split hits may have multiple primary alignment records)
		pa->primary_alignment = ((pa->flags & 0x0100) != 0) ? true : false ;
		//0x0200 the read fails platform/vendor quality checks
		pa->platform_quality_fail = ((pa->flags & 0x0200) != 0) ? true : false ;
		//0x0400 the read is either a PCR duplicate or an optical duplicate
		pa->pcr_duplicate = ((pa->flags & 0x0400) != 0) ? true : false ;




		if (pa->num_cigar>0) {
			pa->cigar_ops=(char*)malloc(sizeof(char)*pa->num_cigar);
			if (pa->cigar_ops==NULL) {
				fprintf(stderr,"Failed to allocate cigar_ops string!\n");
				exit(1);
			}
			pa->cigar_lengths=(uint32_t*)malloc(sizeof(uint32_t)*pa->num_cigar);
			if (pa->cigar_lengths==NULL) {
				fprintf(stderr,"Failed to allocate cigar_lengths int32_t array!\n");
				exit(1);
			}
			pretty_cigar_parse(pa);
		}

		int32_t read_length;
		if (pa->colour_space) {
			read_length=strlen(pa->cs_string)-1;
		} else {
			read_length=strlen(pa->read_string);
		} 
		int32_t expected_length=0;
		int32_t i;
		for (i=0; i<pa->num_cigar; i++) {
			if ((pa->cigar_ops[i]!='H' && !pa->colour_space) && pa->cigar_ops[i]!='P' && pa->cigar_ops[i]!='D' && pa->cigar_ops[i]!='N') {
				expected_length+=pa->cigar_lengths[i];
			}
		}
		if (expected_length>0 && expected_length!=read_length) {
			fprintf(stderr,"Failed to parse, %s, expected length %d, read_length %d\n",pa->read_name,expected_length,read_length);
			return NULL;
		}
	
		return pa;
}

int32_t pretty_delete_cigar_field(pretty * pa, int32_t field) {
	if (pa->num_cigar<=field) {
		fprintf(stderr,"pretty_remap : A2 fatal error has occured!\n");
		exit(1);
	}
	int32_t k;
	for (k=field+1; k<pa->num_cigar; k++) {
		pa->cigar_ops[k-1]=pa->cigar_ops[k];
		pa->cigar_lengths[k-1]=pa->cigar_lengths[k];
	}
	pa->num_cigar--;
	return 0;
}


int32_t pretty_remap_header(FILE* f,char* contig_name,uint32_t offset_start, uint32_t offset_end)  {
	fprintf(f,"@HD\tVN:1\tSO:unsorted\n");
	fprintf(f,"@CO\tSUBSEQUENCE\tSTART:%u\tEND:%u\n",offset_start,offset_end);
	fprintf(f,"@SQ\tSN:%s\tLN:%u\n",contig_name,offset_end-offset_start+1);
	return 0;
}
int32_t pretty_header(FILE* f,char* contig_name,uint32_t sequence_size)  {
	fprintf(f,"@HD\tVN:1\tSO:unsorted\n");
	fprintf(f,"@SQ\tSN:%s\tLN:%u\n",contig_name,sequence_size);
	return 0;
}

//returns 0 if cannot remap!
int32_t pretty_remap(pretty * pa, uint32_t offset_start, uint32_t offset_end) {
	if (offset_start<=0 || offset_end<=0 || (offset_end<=offset_start)) {
		fprintf(stderr,"pretty_remap : offset_start and offset_end must be 1 based index!\n");
		exit(1);
	}
	int32_t genome_alignment_length=0;
	//Trim the begining if need be
	int32_t i;
	for (i=0; i<pa->num_cigar; i++) {
		switch (pa->cigar_ops[i]) {
			case 'M':
			case 'D':
			case 'N':
				genome_alignment_length+=pa->cigar_lengths[i];
				break;
			default:
				break;
		}
	}
	uint32_t genome_start = pa->genome_start_unpadded;
	uint32_t genome_end = pa->genome_start_unpadded+genome_alignment_length-1;
	if (genome_end < offset_start || genome_start > offset_end) {
		return 0;
	}
	//fprintf("rn: %s\n",pa->read_name);
	if (genome_start < offset_start) {
		int32_t trim_genome_front = offset_start - genome_start;
		genome_alignment_length -= trim_genome_front;
		//fix the CIGAR to account for this
		//insert trim of length 0 anyway!?
		if (pa->cigar_ops[0]!='S' && pa->cigar_ops[0]!='H') {
			pa->num_cigar++;
			char * temp_ops = (char*)malloc(sizeof(char)*pa->num_cigar);
			if (temp_ops==NULL) {
				fprintf(stderr,"pretty_remap : temp_ops, malloc failed!\n");
				exit(1);
			}
			uint32_t * temp_lengths = (uint32_t *)malloc(sizeof(uint32_t)*pa->num_cigar);			
			if (temp_lengths==NULL) {
				fprintf(stderr,"pretty_remap : temp_lengths, malloc failed!\n");
				exit(1);
			}
			temp_ops[0]='S'; temp_lengths[0]=0;
			memcpy(temp_ops+1,pa->cigar_ops,(pa->num_cigar-1)*sizeof(char));
			memcpy(temp_lengths+1,pa->cigar_lengths,(pa->num_cigar-1)*sizeof(uint32_t));
			free(pa->cigar_ops); free(pa->cigar_lengths);
			pa->cigar_ops=temp_ops; pa->cigar_lengths=temp_lengths;
		}
		//pa->cigar_lengths[0]+=trim_genome_front;
		uint32_t to_trim=trim_genome_front;
		//fprintf(stderr,"Trim %d from front\n",to_trim);
		int32_t j=1; int32_t deletions=0;
		int32_t inserts=0;
		while(to_trim>0) {
			if (pa->num_cigar<=j || pa->cigar_ops[j]=='S' || pa->cigar_ops[j]=='H') {
				fprintf(stderr,"pretty_remap : A1 fatal error has occured!\n");
				exit(1);
			}
			if ( pa->cigar_ops[j]=='P') {
				pretty_delete_cigar_field(pa,j);
			} else if (pa->cigar_ops[j]=='I') {
				pa->cigar_lengths[0]+=pa->cigar_lengths[j];
				inserts+=pa->cigar_lengths[j];
				pretty_delete_cigar_field(pa,j);
			} else if (pa->cigar_ops[j]=='D' || pa->cigar_ops[j]=='N') {
				if (pa->cigar_lengths[j]<=to_trim) {
					deletions+=pa->cigar_lengths[j];
					to_trim-=pa->cigar_lengths[j];
					pretty_delete_cigar_field(pa,j);
				} else {
					deletions+=to_trim;
					pa->cigar_lengths[j]-=to_trim;
					to_trim=0;
				}
			} else if (pa->cigar_lengths[j]<=to_trim) {
				pa->cigar_lengths[0]+=pa->cigar_lengths[j];
				to_trim-=pa->cigar_lengths[j];
				pretty_delete_cigar_field(pa,j);
			} else {
				pa->cigar_lengths[0]+=to_trim;
				pa->cigar_lengths[j]-=to_trim;
				to_trim=0;
			}
		}
		//check if its hard clip, then remove seq and qual if so
		if (pa->cigar_ops[0]=='H') {
			if (pa->read_string[0]!='*') {
				uint32_t i;
				uint32_t lim=strlen(pa->read_string)-trim_genome_front-inserts+deletions;
				for (i=0; i<lim; i++) {
					pa->read_string[i]=
						pa->read_string[i+trim_genome_front+inserts-deletions];
				}
				pa->read_string[i]='\0';
			}
			if (pa->read_qualities[0]!='*' || pa->read_qualities[1]!='\0') {
				uint32_t i;
				uint32_t lim=strlen(pa->read_qualities)-trim_genome_front-inserts+deletions;
				for (i=0; i<lim; i++) {
					pa->read_qualities[i]=
						pa->read_qualities[i+trim_genome_front+inserts-deletions];
				}
				pa->read_qualities[i]='\0';
			}
		}
		pa->genome_start_unpadded=1;
	} else {
		pa->genome_start_unpadded+=-offset_start+1;
	}

	if (genome_end > offset_end) {
		int32_t trim_genome_back = genome_end - offset_end;
		genome_alignment_length -= trim_genome_back;
		//fix the CIGAR to account for this
		//insert trim of length 0 anyway!?
		if (pa->cigar_ops[pa->num_cigar-1]!='S' && pa->cigar_ops[pa->num_cigar-1]!='H') {
			pa->num_cigar++;
			char * temp_ops = (char*)malloc(sizeof(char)*pa->num_cigar);
			if (temp_ops==NULL) {
				fprintf(stderr,"pretty_remap : temp_ops, malloc failed!\n");
				exit(1);
			}
			uint32_t * temp_lengths = (uint32_t *)malloc(sizeof(uint32_t)*pa->num_cigar);			
			if (temp_lengths==NULL) {
				fprintf(stderr,"pretty_remap : temp_lengths, malloc failed!\n");
				exit(1);
			}
			temp_ops[pa->num_cigar-1]='S'; temp_lengths[pa->num_cigar-1]=0;
			memcpy(temp_ops,pa->cigar_ops,(pa->num_cigar-1)*sizeof(char));
			memcpy(temp_lengths,pa->cigar_lengths,(pa->num_cigar-1)*sizeof(uint32_t));
			free(pa->cigar_ops); free(pa->cigar_lengths);
			pa->cigar_ops=temp_ops; pa->cigar_lengths=temp_lengths;
		}
		//printf("%d vs after %d\n",pa->cigar_lengths[pa->num_cigar-1],pa->cigar_lengths[pa->num_cigar-1]+trim_genome_back);
		uint32_t to_trim=trim_genome_back;
		//fprintf(stderr,"Trim %d from back\n",to_trim);
		int32_t j=pa->num_cigar-2; int32_t deletions=0;
		int32_t inserts=0;
		while(to_trim>0) {
			if (j<0 || pa->cigar_ops[j]=='S' || pa->cigar_ops[j]=='H') {
				fprintf(stderr,"pretty_remap : A3 fatal error has occured!, offending read %s, %d\n",pa->read_name,j);
				exit(1);
			}
			if (pa->cigar_ops[j]=='P') {
				pretty_delete_cigar_field(pa,j);
				j--;
			} else if (pa->cigar_ops[j]=='I') {
				pa->cigar_lengths[pa->num_cigar-1]+=pa->cigar_lengths[j];
				inserts+=pa->cigar_lengths[j];
				pretty_delete_cigar_field(pa,j);
				j--;
			} else if (pa->cigar_ops[j]=='D' || pa->cigar_ops[j]=='N') {
				if (pa->cigar_lengths[j]<=to_trim) {
					deletions+=pa->cigar_lengths[j];
					to_trim-=pa->cigar_lengths[j];
					pretty_delete_cigar_field(pa,j);
					j--;
				} else {
					deletions+=to_trim;
					pa->cigar_lengths[j]-=to_trim;
					to_trim=0;
				}
			} else if (pa->cigar_lengths[j]<=to_trim) {
				//printf("%d%c len %d\n",pa->cigar_lengths[j],pa->cigar_ops[j],pa->cigar_lengths[pa->num_cigar-1]);
				pa->cigar_lengths[pa->num_cigar-1]+=pa->cigar_lengths[j];
				to_trim-=pa->cigar_lengths[j];
				pretty_delete_cigar_field(pa,j);
				//printf("%d%c len %d\n",pa->cigar_lengths[j],pa->cigar_ops[j],pa->cigar_lengths[pa->num_cigar-1]);
				j--;
			} else {
				//printf("Trimming %d from -\n",to_trim);
				pa->cigar_lengths[pa->num_cigar-1]+=to_trim;
				pa->cigar_lengths[j]-=to_trim;
				to_trim=0;
			}
		//	printf("%d%c\n",pa->cigar_lengths[pa->num_cigar-1],pa->cigar_ops[pa->num_cigar-1]);
		}
		//check if its hard clip, then remove seq and qual if so
		if (pa->cigar_ops[pa->num_cigar-1]=='H') {
			if (pa->read_string[0]!='*') {
				pa->read_string[strlen(pa->read_string)-trim_genome_back-inserts+deletions]='\0';
			}
			if (pa->read_qualities[0]!='*' || pa->read_qualities[1]!='\0') {
				pa->read_qualities[strlen(pa->read_qualities)-trim_genome_back-inserts+deletions]='\0';
			}
		}
	}
	pa->mate_genome_start_unpadded+=-offset_start+1;
	if (pa->mate_genome_start_unpadded<1 || pa->mate_genome_start_unpadded>offset_end) {
		pa->mate_genome_start_unpadded=0;
	}
	if (genome_alignment_length<2) {
		return 0;
	}
	free(pa->cigar);
	pa->cigar=pretty_reformat_cigar_string(pa);
	return 1;	
}


/*void pretty_fill(fastas * fs, pretty * pa) {
		fasta * reference = fastas_get_reference(fs,pa->reference_name);
		if (reference!=NULL) {
			pa->genome_sequence=reference->sequence;
			pa->genome_length=reference->length;
			pretty_genome(pa);
			pretty_ls(pa);
			if (pa->colour_space) {
				pretty_cs(pa);
			}
			pretty_qualities(pa);
			pretty_stats(pa);
		} else {
			fprintf(stderr,"Could not find reference sequence!\n");
			exit(1);
		}
}*/


