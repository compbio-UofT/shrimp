#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "sam2pretty_lib.h"
#include "render.h"


size_t render_sam_unaligned_bound(pretty * pa ) { 
	size_t buffer_size=13*SIZE_TAB+SIZE_NEWLINE+SIZE_NULL+strlen(pa->read_name)+SIZE_FLAG+
		1+1+SIZE_MAPQ+1+1+1+1+1+1;
	if (pa->colour_space) {
		if (pa->has_cs_qualities) {
			buffer_size+=strlen(pa->cs_qualities)+SIZE_TAB+SIZE_SAM_AUX;
		}
		buffer_size+=strlen(pa->cs_string)+SIZE_TAB+SIZE_SAM_AUX;
	}
	if (pa->has_read_group) {
		buffer_size+=strlen(pa->read_group)+SIZE_TAB+SIZE_SAM_AUX;
	}
	if (pa->has_r2) {
		buffer_size+=strlen(pa->r2);
	}
	return buffer_size;
}

size_t render_sam_unaligned_string(pretty * pa, char * buffer, size_t buffer_size) {
	size_t position = snprintf(buffer,buffer_size,"%s\t%d\t%s\t%d\t%d\t%s\t%s\t%d\t%d\t%s\t%s",
		pa->read_name,
		pa->flags | 0x4 | 0x8,
		"*",
		0,
		0,
		"*",
		"*",
		0,
		0,
		pa->colour_space ? "*" : pa->read_string,
		pa->colour_space ? "*" : pa->read_qualities);
	if (pa->colour_space) {
		if (pa->has_cs_qualities) {
			position+=snprintf(buffer+position,buffer_size-position,"\tCQ:Z:%s",pa->cs_qualities);
		}
		position+=snprintf(buffer+position,buffer_size-position,"\tCS:Z:%s",pa->cs_string);
	}
	if (pa->has_read_group) {
		position+=snprintf(buffer+position,buffer_size-position,"\tRG:Z:%s",pa->read_group);
	}
	if (pa->has_r2) {
		position+=snprintf(buffer+position,buffer_size-position,"\tR2:Z:%s",pa->r2);
	}
	if (pa->aux!=NULL) {
		position+=snprintf(buffer+position,buffer_size-position,"\t%s",pa->aux);
	}
	return position;
}

void render_sam_unaligned(pretty * pa,bool inplace) {
	size_t buffer_size = render_sam_unaligned_bound(pa);
	char buffer[buffer_size];
	size_t position=render_sam_unaligned_string(pa,buffer,buffer_size);
	if (!inplace) {
		position+=snprintf(buffer+position,buffer_size-position,"\n");
		if (pa->sam_string!=NULL) {
			free(pa->sam_string);
		}
		if (buffer_size<position) {
			fprintf(stderr,"Failed to properly bound memory needed for updated sam alignment\n");
			exit(1);
		}
		pa->sam_string=(char*)malloc(sizeof(char)*(position+1));
		if (pa->sam_string==NULL) {
			fprintf(stderr,"Failed to allocate memory to store same string!\n");
			exit(1);
		}
		strcpy(pa->sam_string,buffer);
	} else {
		if (position>pa->sam_string_length) {
			fprintf(stderr,"failed in place update! %lu vs %lu\n",position,pa->sam_string_length);
			exit(1);
		}
		strcpy(pa->sam_string,buffer);

	}
}


size_t render_fastx_bound(pretty * pa )  {
	char * read; char * qualities=NULL;
	if (pa->colour_space) {
		read=pa->cs_string;
		if (pa->has_cs_qualities) {
			qualities=pa->cs_qualities;
			assert(strlen(pa->cs_string)-1==strlen(pa->cs_qualities));
		}
	} else {
		read=pa->read_string;
		if (pa->read_qualities[0]!='*') {
			qualities=pa->read_qualities;
		}
	}
	size_t buffer_size=1+strlen(pa->read_name)+1
		+strlen(read)+1;
	if (qualities!=NULL) {
		buffer_size+=1+1
			+strlen(qualities)+1;
	}	
	assert(buffer_size<1000);
	return buffer_size;
}

size_t render_fastx_string(pretty * pa, char * buffer, size_t buffer_size) {
	char * read; char * qualities=NULL;
	if (pa->colour_space) {
		read=pa->cs_string;
		if (pa->has_cs_qualities) {
			qualities=pa->cs_qualities;
			assert(strlen(pa->cs_string)-1==strlen(pa->cs_qualities));
		}
	} else {
		read=pa->read_string;
		if (pa->read_qualities[0]!='*') {
			qualities=pa->read_qualities;
		}
	}
	size_t position;
	if (read!=NULL && read[0]!='\0' && (read[0]!='*' || read[1]!='\0')) {
		if (qualities==NULL) {
			//fasta
			position=snprintf(buffer,buffer_size,">%s\n%s\n",pa->read_name,read);
		} else {
			//fastq
			position=snprintf(buffer,buffer_size,"@%s\n%s\n+\n%s",pa->read_name,read,qualities);
		} 
	} else {
		*buffer='\0';
	}
	return position;
}


void render_fastx(pretty* pa, bool inplace ) { 
	size_t buffer_size = render_fastx_bound(pa);
	char buffer[buffer_size];
	size_t position=render_fastx_string(pa,buffer,buffer_size);
	if (!inplace) {
		position+=snprintf(buffer+position,buffer_size-position,"\n");
		if (pa->sam_string!=NULL) {
			free(pa->sam_string);
		}
		if (buffer_size<position) {
			fprintf(stderr,"Failed to properly bound memory needed for updated sam alignment\n");
			exit(1);
		}
		pa->sam_string=(char*)malloc(sizeof(char)*(position+1));
		if (pa->sam_string==NULL) {
			fprintf(stderr,"Failed to allocate memory to store same string!\n");
			exit(1);
		}
		strcpy(pa->sam_string,buffer);
	} else {
		if (position>pa->sam_string_length) {
			fprintf(stderr,"Failed in place update! %lu vs %lu\n",position,pa->sam_string_length);
			exit(1);
		}
		strcpy(pa->sam_string,buffer);
	}
}

size_t render_sam_bound(pretty * pa ) {
	size_t buffer_size=13*SIZE_TAB+SIZE_NEWLINE+SIZE_NULL+strlen(pa->read_name)+SIZE_FLAG+
		strlen(pa->reference_name)+SIZE_POS+SIZE_MAPQ+strlen(pa->cigar)+
		strlen(pa->mate_reference_name)+SIZE_MPOS+SIZE_ISIZE+strlen(pa->read_string)+
		strlen(pa->read_qualities);
	assert(buffer_size<1000);
	if (pa->has_score) {
		buffer_size+=SIZE_TAB+SIZE_TAB+SIZE_SAM_AUX;
	}
	if (pa->has_edit_distance) {
		buffer_size+=SIZE_TAB+SIZE_TAB+SIZE_SAM_AUX;
	}
	if (pa->colour_space) {
		if (pa->has_cs_qualities) {
			buffer_size+=strlen(pa->cs_qualities)+SIZE_TAB+SIZE_SAM_AUX;
		}
		buffer_size+=strlen(pa->cs_string)+SIZE_TAB+SIZE_SAM_AUX;
	}
	if (pa->has_cs_mismatches) {
		buffer_size+=SIZE_32bit+SIZE_TAB+SIZE_SAM_AUX;
	}
	if (pa->has_cs_edit_string) {
		buffer_size+=strlen(pa->cs_edit_string)+SIZE_TAB+SIZE_SAM_AUX;
	}
	if (pa->has_read_group) {
		buffer_size+=strlen(pa->read_group)+SIZE_TAB+SIZE_SAM_AUX;
	}
	if (pa->has_r2) {
		buffer_size+=strlen(pa->r2)+SIZE_TAB+SIZE_SAM_AUX;
	}
	if (pa->has_ih) {
		buffer_size+=SIZE_32bit+SIZE_TAB+SIZE_SAM_AUX;
	}
	if (pa->has_hi) {
		buffer_size+=SIZE_32bit+SIZE_TAB+SIZE_SAM_AUX;
	}
	if (pa->has_h0) {
		buffer_size+=SIZE_32bit+SIZE_TAB+SIZE_SAM_AUX;
	}
	if (pa->has_h1) {
		buffer_size+=SIZE_32bit+SIZE_TAB+SIZE_SAM_AUX;
	}
	if (pa->has_h2) {
		buffer_size+=SIZE_32bit+SIZE_TAB+SIZE_SAM_AUX;
	}
	int i;
	for (i=0; i<SAM2PRETTY_NUM_ZS; i++) {
		if ((pa->has_zs&(1<<i))!=0) {
			buffer_size+=SIZE_32bit+SIZE_TAB+SIZE_SAM_AUX; //Z
		}	
	}
	if (pa->aux!=NULL) {
		buffer_size+=SIZE_TAB+strlen(pa->aux);
	}
	return buffer_size;
}

size_t render_sam_string(pretty * pa, char * buffer, size_t buffer_size) {
	if (!pa->mapped) {
		return render_sam_unaligned_string(pa,buffer,buffer_size);
	}	
	pa->flags=pretty_get_flag(pa);
	size_t position = snprintf(buffer,buffer_size,"%s\t%d\t%s\t%d\t%d\t%s\t%s\t%d\t%d\t%s\t%s",
		pa->read_name,
		pa->flags,
		pa->reference_name,
		pa->genome_start_unpadded,
		pa->mapq >= 4 ? pa->mapq : 0,
		pa->cigar,
		(strcmp(pa->reference_name,pa->mate_reference_name)==0 ? "=" : pa->mate_reference_name),
		pa->mate_genome_start_unpadded,
		pa->isize,
		pa->read_string,
		pa->read_qualities);
	if (pa->has_score) {
		position+=snprintf(buffer+position,buffer_size-position,"\tAS:i:%d",pa->score);
	}
	int i;	
	for (i=0; i<SAM2PRETTY_NUM_ZS; i++) {
		if ((pa->has_zs&(1<<i))!=0) {
			position+=snprintf(buffer+position,buffer_size-position,"\tZ%d:i:%d",i,tnlog(pa->z[i]));
		}	
	}
	if (pa->has_edit_distance) {
		position+=snprintf(buffer+position,buffer_size-position,"\tNM:i:%d",pa->edit_distance);
	}	 
	if (pa->colour_space) {
		if (pa->has_cs_qualities) {
			position+=snprintf(buffer+position,buffer_size-position,"\tCQ:Z:%s",pa->cs_qualities);
		}
		position+=snprintf(buffer+position,buffer_size-position,"\tCS:Z:%s",pa->cs_string);
	}
	if (pa->has_cs_mismatches) {
		position+=snprintf(buffer+position,buffer_size-position,"\tCM:i:%d",pa->cs_mismatches);
	}
	if (pa->has_cs_edit_string) {
		position+=snprintf(buffer+position,buffer_size-position,"\tXX:Z:%s",pa->cs_edit_string);
	}
	if (pa->has_r2) {
		position+=snprintf(buffer+position,buffer_size-position,"\tR2:Z:%s",pa->r2);
	}
	if (pa->has_hi) {
		position+=snprintf(buffer+position,buffer_size-position,"\tHI:i:%d",pa->hi);
	}
	if (pa->has_ih) {
		position+=snprintf(buffer+position,buffer_size-position,"\tIH:i:%d",pa->ih);
	}
	if (pa->has_h0) {
		position+=snprintf(buffer+position,buffer_size-position,"\tH0:i:%d",pa->h0);
	}
	if (pa->has_h1) {
		position+=snprintf(buffer+position,buffer_size-position,"\tH1:i:%d",pa->h1);
	}
	if (pa->has_h2) {
		position+=snprintf(buffer+position,buffer_size-position,"\tH2:i:%d",pa->h2);
	}
	if (pa->has_read_group) {
		position+=snprintf(buffer+position,buffer_size-position,"\tRG:Z:%s",pa->read_group);
	}
	if (pa->aux!=NULL) {
		position+=snprintf(buffer+position,buffer_size-position,"\t%s",pa->aux);
	}
	return position;
}

void render_sam(pretty * pa, bool inplace) {
	size_t buffer_size = render_sam_bound(pa);  
	char buffer[buffer_size];
	size_t position=render_sam_string(pa,buffer,buffer_size);
	if (!inplace) {
		position+=snprintf(buffer+position,buffer_size-position,"\n");
		if (pa->sam_string!=NULL) {
			free(pa->sam_string);
		}
		if (buffer_size<position) {
			fprintf(stderr,"Failed to properly bound memory needed for updated sam alignment\n");
			exit(1);
		}
		pa->sam_string=(char*)malloc(sizeof(char)*(position+1));
		if (pa->sam_string==NULL) {
			fprintf(stderr,"Failed to allocate memory to store same string!\n");
			exit(1);
		}
		strcpy(pa->sam_string,buffer);
	} else {
		if (position>pa->sam_string_length) {
			fprintf(stderr,"Failed in place update! %lu vs %lu\n",position,pa->sam_string_length);
			fprintf(stderr,"%s\n",pa->read_name);
			fprintf(stderr,"%s\n",pa->sam_string);
			fprintf(stderr,"|%s|\n",buffer);
			fprintf(stderr,"|%s|\n",pa->aux);
			fprintf(stderr,"%d\n",pa->has_zs);
			exit(1);
		}
		strcpy(pa->sam_string,buffer);

	}
}
