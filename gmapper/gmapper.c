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
#include "../common/sw-post.h"

/* heaps */
/*
typedef struct ptr_and_sz {
  void * ptr;
  size_t sz;
} ptr_and_sz;
*/
//DEF_HEAP(uint32_t, char *, out)
DEF_HEAP(uint32_t, struct ptr_and_sz, out)



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


/*
 * Free memory allocated by this read.
 */
void 
read_free(struct read_entry * re)
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
read_free_full(struct read_entry * re, count_t * counter)
{
  if (shrimp_mode == MODE_COLOUR_SPACE && Qflag) {
    if (re->crossover_score != NULL) {
      free(re->crossover_score);
      re->crossover_score = NULL;
    }
  }
  if (re->mapidx[0] != NULL)
    my_free(re->mapidx[0], n_seeds * re->max_n_kmers * sizeof(re->mapidx[0][0]), counter, "mapidx [%s]", re->name);
  if (re->mapidx[1] != NULL)
    my_free(re->mapidx[1], n_seeds * re->max_n_kmers * sizeof(re->mapidx[0][0]), counter, "mapidx [%s]", re->name);

  read_free_hit_list(re, counter);
  read_free_anchor_list(re, counter);

  if (re->n_ranges > 0)
    free(re->ranges);

  if (re->n_final_unpaired_hits > 0) {
    int i;
    for (i = 0; i < re->n_final_unpaired_hits; i++)
      free_sfrp(&re->final_unpaired_hits[i].sfrp, re, counter);
    my_free(re->final_unpaired_hits, re->n_final_unpaired_hits * sizeof(re->final_unpaired_hits[0]),
            counter, "final_unpaired_hits [%s]", re->name);
    re->n_final_unpaired_hits = 0;
    re->final_unpaired_hits = NULL;
  }

  free(re->read[0]);
  free(re->read[1]);
  read_free(re);
}


void
readpair_free_full(pair_entry * peP, count_t * counterP)
{
  int nip, i;

  if (peP->n_final_paired_hits > 0) {
    my_free(peP->final_paired_hits, peP->n_final_paired_hits * sizeof(peP->final_paired_hits[0]),
	counterP, "final_paired_hits [%s,%s]", peP->re[0]->name, peP->re[1]->name);
    peP->n_final_paired_hits = 0;
    peP->final_paired_hits = NULL;
  }

  for (nip = 0; nip < 2; nip++) {
    if (peP->final_paired_hit_pool_size[nip] > 0) {
      for (i = 0; i < peP->final_paired_hit_pool_size[nip]; i++) {
	free_sfrp(&peP->final_paired_hit_pool[nip][i].sfrp, peP->re[nip], counterP);
	if (peP->final_paired_hit_pool[nip][i].n_paired_hit_idx > 0) {
	  my_free(peP->final_paired_hit_pool[nip][i].paired_hit_idx, peP->final_paired_hit_pool[nip][i].n_paired_hit_idx * sizeof(peP->final_paired_hit_pool[nip][i].paired_hit_idx[0]),
	      counterP, "paired_hit_idx [%s]", peP->re[nip]->name);
	  peP->final_paired_hit_pool[nip][i].n_paired_hit_idx = 0;
	  peP->final_paired_hit_pool[nip][i].paired_hit_idx = NULL;
	}
      }
      my_free(peP->final_paired_hit_pool[nip], peP->final_paired_hit_pool_size[nip] * sizeof(peP->final_paired_hit_pool[nip][0]),
              counterP, "final_paired_hit_pool[%d] [%s,%s]", nip, peP->re[0]->name, peP->re[1]->name);
      peP->final_paired_hit_pool_size[nip] = 0;
      peP->final_paired_hit_pool[nip] = NULL;
    }
  }

  read_free_full(peP->re[0], counterP);
  read_free_full(peP->re[1], counterP);
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

/*
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
*/

/*
 * Compute range limitations for this read.
 */
/*
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
*/

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
	re->seq[i] = '\0';
	if (Qflag) {
		re->qual[i] = '\0';
	}
	return;
}

/*
 * Launch the threads that will scan the reads
 */
static bool
launch_scan_threads(fasta_t fasta, fasta_t left_fasta, fasta_t right_fasta)
{
  llint last_nreads, last_time_usecs;

  bool read_more = true, more_in_left_file = true, more_in_right_file=true;

  /* initiate the thread buffers */
  //thread_output_buffer_sizes = (size_t *)xcalloc_m(sizeof(size_t) * num_threads, "thread_output_buffer_sizes");
  thread_output_buffer_sizes = (size_t *)
    my_calloc(num_threads * sizeof(size_t),
	      &mem_thread_buffer, "thread_output_buffer_sizes");
  //thread_output_buffer_filled = (char * *)xcalloc_m(sizeof(char *) * num_threads, "thread_output_buffer_filled");
  thread_output_buffer_filled = (char * *)
    my_calloc(num_threads * sizeof(char *),
	      &mem_thread_buffer, "thread_output_buffer_filled");
  //thread_output_buffer = (char * *)xcalloc_m(sizeof(char *) * num_threads, "thread_output_buffer");
  thread_output_buffer = (char * *)
    my_calloc(num_threads * sizeof(char *),
	      &mem_thread_buffer, "thread_output_buffer");
  //thread_output_buffer_chunk = (unsigned int *)xcalloc_m(sizeof(unsigned int) * num_threads, "thread_output_buffer_chunk");
  thread_output_buffer_chunk = (unsigned int *)
    my_calloc(num_threads * sizeof(unsigned int),
	      &mem_thread_buffer, "thread_output_buffer_chunk");
  
  unsigned int current_thread_chunk = 1;
  unsigned int next_chunk_to_print = 1;
  struct heap_out h; 
  heap_out_init(&h, thread_output_heap_capacity );

  if (progress > 0) {
    fprintf(stderr, "done r/hr r/core-hr\n");
    last_nreads = 0;
    last_time_usecs = gettimeinusecs();
  }

#pragma omp parallel shared(read_more,more_in_left_file,more_in_right_file, fasta) num_threads(num_threads)
  {
    int thread_id = omp_get_thread_num();
    struct read_entry * re_buffer;
    int load, i;
    //llint before, after;
    //re_buffer = (struct read_entry *)xmalloc_m(chunk_size * sizeof(re_buffer[0]), "re_buffer");
    re_buffer = (read_entry *)my_malloc(chunk_size * sizeof(re_buffer[0]), &mem_thread_buffer, "re_buffer");

    while (read_more) {
      memset(re_buffer, 0, chunk_size * sizeof(re_buffer[0]));

      //before = rdtsc();
      TIME_COUNTER_START(tpg.wait_tc);

      //Read in this threads 'chunk'
#pragma omp critical (fill_reads_buffer)
      {
	//after = rdtsc();
	//tpg.wait_ticks += MAX(after - before, 0);
	TIME_COUNTER_STOP(tpg.wait_tc);

	thread_output_buffer_chunk[thread_id]=current_thread_chunk++;

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
	  if (nreads_mod >= progress) {
	    llint time_usecs = gettimeinusecs();
#pragma omp critical (cs_stderr)
	    {
	    fprintf(stderr, "%lld %d %d.\r", nreads,
		    (int)(((double)(nreads - last_nreads)/(double)(time_usecs - last_time_usecs)) * 3600.0 * 1.0e6),
		    (int)(((double)(nreads - last_nreads)/(double)(time_usecs - last_time_usecs)) * 3600.0 * 1.0e6 * (1/(double)num_threads)) );
	    }
	    last_nreads = nreads;
	    last_time_usecs = time_usecs;
	  }
	  nreads_mod %= progress;
	}
      } // end critical section

      if (pair_mode != PAIR_NONE)
	assert(load % 2 == 0); // read even number of reads

      thread_output_buffer_sizes[thread_id] = thread_output_buffer_initial;
      //thread_output_buffer[thread_id] = (char *)xmalloc_m(sizeof(char) * thread_output_buffer_sizes[thread_id], "thread_buffer");
      thread_output_buffer[thread_id] = (char *)
	my_malloc(thread_output_buffer_sizes[thread_id] * sizeof(char),
		  &mem_thread_buffer, "thread_output_buffer[]");
      thread_output_buffer_filled[thread_id] = thread_output_buffer[thread_id];
      thread_output_buffer[thread_id][0] = '\0';

      for (i = 0; i < load; i++) {
	// if running in paired mode and first foot is ignored, ignore this one, too
	if (pair_mode != PAIR_NONE && i % 2 == 1 && re_buffer[i-1].ignore) {
	  //read_free(&re_buffer[i-1]);
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
	  if (pair_mode == PAIR_NONE) {
	    trim_read(&re_buffer[i]);
	  } else if (i % 2 == 1){
	    if (trim_first) {
	      trim_read(&re_buffer[i-1]);
	    }
	    if (trim_second) {
	      trim_read(&re_buffer[i]);
	    }
	  }
	}
        if (shrimp_mode == MODE_LETTER_SPACE && trim_illumina) { 
		int trailing_Bs=0;
		for (int j=0; j<(int)strlen(re_buffer[i].qual); j++) {
			if (re_buffer[i].qual[j]!='B') {
				trailing_Bs=0;
			} else {
				trailing_Bs+=1;
			}
		}
		if (trailing_Bs>0) {
			re_buffer[i].seq[strlen(re_buffer[i].seq)-trailing_Bs]='\0';
			re_buffer[i].qual[strlen(re_buffer[i].qual)-trailing_Bs]='\0';
		}
	}
	//compute average quality value
	re_buffer[i].read_len = strlen(re_buffer[i].seq);
	if (Qflag && !ignore_qvs && min_avg_qv >= 0) {
	  //fprintf(stderr, "read:[%s] qual:[%s]", re_buffer[i].name, re_buffer[i].qual);
	  re_buffer[i].avg_qv = 0;
	  for (char * c = re_buffer[i].qual; *c != 0; c++) {
	    re_buffer[i].avg_qv += (*c - qual_delta);
	  }
	  //fprintf(stderr, " avg_qv:%d\n", re_buffer[i].avg_qv);
	}
	if (Qflag && !ignore_qvs && !no_qv_check) {
		for (char * c =re_buffer[i].qual; *c !=0; c++) {
			int qual_value=(*c-qual_delta);
			if (qual_value<-10 || qual_value>50) {
				fprintf(stderr,"The qv-offset might be set incorrectly! Currenty qvs are interpreted as PHRED+%d\
 and a qv of %d was observed. To disable this error, etiher set the offset correctly or disable this check (see README).\n",qual_delta,qual_value);
				exit(1);
			}
		}
	}

	re_buffer[i].max_n_kmers = re_buffer[i].read_len - min_seed_span + 1;
	re_buffer[i].read[0] = fasta_sequence_to_bitfield(fasta, re_buffer[i].seq);
	if (shrimp_mode == MODE_COLOUR_SPACE) {
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
	if (re_buffer[i].max_n_kmers < 0)
	  re_buffer[i].max_n_kmers = 0;
	if (re_buffer[i].read_len > 0)
	  re_buffer[i].avg_qv /= re_buffer[i].read_len;

	//Check if we can actually use this read
	if (//re_buffer[i].max_n_kmers <= 0
	    //|| 
	    re_buffer[i].read_len > longest_read_len
	    || (Qflag && !ignore_qvs && re_buffer[i].avg_qv < min_avg_qv) // ignore reads with low avg qv
	    ) {
	  if (re_buffer[i].max_n_kmers <= 0) {
	    fprintf(stderr, "warning: skipping read [%s]; smaller then any seed!\n",
		    re_buffer[i].name);
	    re_buffer[i].max_n_kmers=1;
	  } else if (re_buffer[i].read_len > longest_read_len) {
	    fprintf(stderr, "warning: skipping read [%s]; it has length %d, maximum allowed is %d. Use --longest-read ?\n",
		    re_buffer[i].name, re_buffer[i].read_len, longest_read_len);
	  } else {
	    //fprintf(stderr, "skipping read [%s] with avg_qv:%g\n", re_buffer[i].name, (double)re_buffer[i].avg_qv);
	  }
	  if (pair_mode == PAIR_NONE) {
	    #pragma omp atomic
	    total_reads_dropped++;

	    read_free_full(&re_buffer[i], &mem_mapping);
	  } else {
	    #pragma omp atomic
	    total_pairs_dropped++;

	    if (i%2 == 1) {
	      read_free_full(&re_buffer[i-1], &mem_mapping);
	      read_free_full(&re_buffer[i], &mem_mapping);
	    } else {
	      read_free_full(&re_buffer[i], &mem_mapping);
	      re_buffer[i].ignore = true;
	    }
	  }
	  continue;	
	}

	re_buffer[i].window_len = (uint16_t)abs_or_pct(window_len,re_buffer[i].read_len);
	// compute position-based crossover scores based on qvs
	if (shrimp_mode == MODE_COLOUR_SPACE && Qflag && !ignore_qvs) {
	  int j;

	  re_buffer[i].crossover_score = (int *)xmalloc(re_buffer[i].read_len * sizeof(re_buffer[i].crossover_score[0]));
	  for (j = 0; j < re_buffer[i].read_len; j++) {
	    re_buffer[i].crossover_score[j] = (int)(score_alpha * log(pr_err_from_qv(re_buffer[i].qual[j] - qual_delta) / 3.0) / log(2.0));
	    if (re_buffer[i].crossover_score[j] > -1) {
	      re_buffer[i].crossover_score[j] = -1;
	    } else if (re_buffer[i].crossover_score[j] < 2*crossover_score) {
	      re_buffer[i].crossover_score[j] = 2*crossover_score;
	    }
	  }
	}

	if (re_buffer[i].range_string != NULL) {
	  assert(0); // not maintained
	  //read_compute_ranges(&re_buffer[i]);
	  free(re_buffer[i].range_string);
	  re_buffer[i].range_string = NULL;
	}
	//free(re_buffer[i].seq);

	// time to do some mapping!
	if (pair_mode == PAIR_NONE)
	  {
	    handle_read(&re_buffer[i], unpaired_mapping_options[0], n_unpaired_mapping_options[0]);
	    read_free_full(&re_buffer[i], &mem_mapping);
	  }
	else if (i % 2 == 1)
	  {
	    pair_entry pe;
	    memset(&pe, 0, sizeof(pe));
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
	    pe.re[0] = &re_buffer[i-1];
	    pe.re[1] = &re_buffer[i];
	    handle_readpair(&pe, paired_mapping_options, n_paired_mapping_options);
	    readpair_free_full(&pe, &mem_mapping);
	  }
      }

      // free unused memory while the buffer waits in the output heap
      size_t new_size = (thread_output_buffer_filled[thread_id] - thread_output_buffer[thread_id]) + 1;
      thread_output_buffer[thread_id] = (char *)
	my_realloc(thread_output_buffer[thread_id], new_size, thread_output_buffer_sizes[thread_id],
		   &mem_thread_buffer, "thread_output_buffer[]");

      //fprintf(stdout,"%s",thread_output_buffer[thread_id]);
#pragma omp critical
      {
	struct heap_out_elem tmp;
	tmp.key = thread_output_buffer_chunk[thread_id];
	//tmp.rest = thread_output_buffer[thread_id];
	tmp.rest.ptr = thread_output_buffer[thread_id];
	tmp.rest.sz = new_size; //thread_output_buffer_sizes[thread_id];
	thread_output_buffer[thread_id] = NULL;	
	heap_out_insert(&h, &tmp);
	heap_out_get_min(&h, &tmp);
	while (h.load > 0 && tmp.key == next_chunk_to_print) {
	  heap_out_extract_min(&h, &tmp);
	  //fprintf(stdout, "%s", tmp.rest);
	  fprintf(stdout, "%s", (char *)tmp.rest.ptr);
	  //free(tmp.rest);
	  my_free(tmp.rest.ptr, tmp.rest.sz,
		  &mem_thread_buffer, "thread_output_buffer[]");
	  next_chunk_to_print++;
	}
      }
    }

    //free(re_buffer);
    my_free(re_buffer, chunk_size * sizeof(re_buffer[0]),
	    &mem_thread_buffer, "re_buffer");
  } // end parallel section

  if (progress > 0)
    fprintf(stderr, "\n");

  //assert(h.load == 0);
  
  struct heap_out_elem tmp;
  while (h.load>0) {
    heap_out_extract_min(&h,&tmp);
    fprintf(stdout,"%s",(char *)tmp.rest.ptr);
    //free(tmp.rest);
    my_free(tmp.rest.ptr, tmp.rest.sz,
	    &mem_thread_buffer, "thread_output_buffer[]");
  }
  
  heap_out_destroy(&h);

  //free(thread_output_buffer_sizes);
  my_free(thread_output_buffer_sizes, sizeof(size_t) * num_threads,
	  &mem_thread_buffer, "thread_output_buffer_sizes");
  //free(thread_output_buffer_filled);
  my_free(thread_output_buffer_filled, sizeof(char *) * num_threads,
	  &mem_thread_buffer, "thread_output_buffer_filled");
  //free(thread_output_buffer);
  my_free(thread_output_buffer, sizeof(char *) * num_threads,
	  &mem_thread_buffer, "thread_output_buffer");
  //free(thread_output_buffer_chunk);
  my_free(thread_output_buffer_chunk, sizeof(unsigned int) * num_threads,
	  &mem_thread_buffer, "thread_output_buffer_chunk");

  return true;
}


static char *
thres_to_buff(char * buff, double * thres)
{
  if (IS_ABSOLUTE(*thres)) {
    sprintf(buff, "%u", -(uint)(*thres));
  } else {
    sprintf(buff, "%.02f%%", *thres);
  }
  return buff;
}

static char *
bool_to_buff(char * buff, bool * val)
{
  sprintf(buff, "%s", *val ? "true" : "false");
  return buff;
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

typedef struct tp_stats_t {
  uint64_t f1_invocs, f1_cells, f1_ticks;
  double f1_secs, f1_cellspersec;
  uint64_t f2_invocs, f2_cells, f2_ticks;
  double f2_secs, f2_cellspersec;
  uint64_t fwbw_invocs, fwbw_cells, fwbw_ticks;
  double fwbw_secs, fwbw_cellspersec;
  double scan_secs, readparse_secs, read_handle_overhead_secs;
  double anchor_list_secs, hit_list_secs;
  double region_counts_secs, mp_region_counts_secs, duplicate_removal_secs;
  double pass1_secs, get_vector_hits_secs, pass2_secs;
} tp_stats_t;

static void
print_statistics()
{
  static char const my_tab[] = "    ";

  //uint64_t f1_invocs[num_threads], f1_cells[num_threads], f1_ticks[num_threads];
  //double f1_secs[num_threads], f1_cellspersec[num_threads];
  uint64_t f1_total_invocs = 0, f1_total_cells = 0;
  double f1_total_secs = 0, f1_total_cellspersec = 0;
  uint64_t f1_calls_bypassed = 0;
  //uint64_t f2_invocs[num_threads], f2_cells[num_threads], f2_ticks[num_threads];
  //double f2_secs[num_threads], f2_cellspersec[num_threads];
  uint64_t f2_total_invocs = 0, f2_total_cells = 0;
  double f2_total_secs = 0, f2_total_cellspersec = 0;
  //uint64_t fwbw_invocs[num_threads], fwbw_cells[num_threads], fwbw_ticks[num_threads];
  //double fwbw_secs[num_threads], fwbw_cellspersec[num_threads];
  uint64_t fwbw_total_invocs = 0, fwbw_total_cells = 0;
  double fwbw_total_secs = 0, fwbw_total_cellspersec = 0;
  //double scan_secs[num_threads], readparse_secs[num_threads], read_handle_overhead_secs[num_threads];
  //double anchor_list_secs[num_threads], hit_list_secs[num_threads];
  //double region_counts_secs[num_threads], mp_region_counts_secs[num_threads], duplicate_removal_secs[num_threads];
  //double pass1_secs[num_threads], get_vector_hits_secs[num_threads], pass2_secs[num_threads];
  double total_scan_secs = 0, total_wait_secs = 0, total_readparse_secs = 0;

  tp_stats_t * tps = (tp_stats_t *)malloc(num_threads * sizeof(tps[0]));
  tpg_t * tpgA = (tpg_t *)malloc(num_threads * sizeof(tpgA[0]));

  double hz;
  double fasta_load_secs;
  fasta_stats_t fs;

  fs = fasta_stats();
  fasta_load_secs = fs->total_secs;
  free(fs);
  hz = cpuhz();

#pragma omp parallel num_threads(num_threads) shared(hz)
  {
    int tid = omp_get_thread_num();

    memcpy(&tpgA[tid], &tpg, sizeof(tpg_t));

    f1_stats(&tps[tid].f1_invocs, &tps[tid].f1_cells, &tps[tid].f1_secs, NULL);

    //tps[tid].f1_secs = (double)tps[tid].f1_ticks / hz;
    tps[tid].f1_cellspersec = (double)tps[tid].f1_cells / tps[tid].f1_secs;
    if (isnan(tps[tid].f1_cellspersec))
      tps[tid].f1_cellspersec = 0;

    if (shrimp_mode == MODE_COLOUR_SPACE) {
      sw_full_cs_stats(&tps[tid].f2_invocs, &tps[tid].f2_cells, &tps[tid].f2_secs);
      post_sw_stats(&tps[tid].fwbw_invocs, &tps[tid].fwbw_cells, &tps[tid].fwbw_secs);
    } else {
      sw_full_ls_stats(&tps[tid].f2_invocs, &tps[tid].f2_cells, &tps[tid].f2_secs);
      tps[tid].fwbw_secs = 0;
    }

    //tps[tid].f2_secs = (double)tps[tid].f2_ticks / hz;
    tps[tid].f2_cellspersec = (double)tps[tid].f2_cells / tps[tid].f2_secs;
    if (isnan(tps[tid].f2_cellspersec))
      tps[tid].f2_cellspersec = 0;

    //tps[tid].fwbw_secs = (double)tps[tid].fwbw_ticks / hz;
    if (shrimp_mode == MODE_COLOUR_SPACE) {
      tps[tid].fwbw_cellspersec = (double)tps[tid].fwbw_cells / tps[tid].fwbw_secs;
      if (isnan(tps[tid].fwbw_cellspersec))
        tps[tid].fwbw_cellspersec = 0;
    }

    tps[tid].readparse_secs = ((double)mapping_wallclock_usecs / 1.0e6) - ((double)tpg.read_handle_usecs / 1.0e6) - time_counter_get_secs(&tpg.wait_tc);

    tps[tid].scan_secs = ((double)tpg.read_handle_usecs / 1.0e6) - tps[tid].f1_secs - tps[tid].f2_secs - tps[tid].fwbw_secs;
    tps[tid].scan_secs = MAX(0, tps[tid].scan_secs);

    //tps[tid].anchor_list_secs = (double)tpg.anchor_list_ticks / hz;
    tps[tid].anchor_list_secs = time_counter_get_secs(&tpg.anchor_list_tc);
    //tps[tid].hit_list_secs = (double)tpg.hit_list_ticks / hz;
    tps[tid].hit_list_secs = time_counter_get_secs(&tpg.hit_list_tc);
    //tps[tid].duplicate_removal_secs = (double)tpg.duplicate_removal_ticks / hz;
    tps[tid].duplicate_removal_secs = time_counter_get_secs(&tpg.duplicate_removal_tc);
    //tps[tid].region_counts_secs = (double)tpg.region_counts_ticks / hz;
    tps[tid].region_counts_secs = time_counter_get_secs(&tpg.region_counts_tc);
    //tps[tid].mp_region_counts_secs = (double)tpg.mp_region_counts_ticks / hz;
    tps[tid].mp_region_counts_secs = time_counter_get_secs(&tpg.mp_region_counts_tc);
    //tps[tid].pass1_secs = (double)tpg.pass1_ticks / hz;
    tps[tid].pass1_secs = time_counter_get_secs(&tpg.pass1_tc);
    //tps[tid].get_vector_hits_secs = (double)tpg.get_vector_hits_ticks / hz;
    tps[tid].get_vector_hits_secs = time_counter_get_secs(&tpg.get_vector_hits_tc);
    //tps[tid].pass2_secs = (double)tpg.pass2_ticks / hz;
    tps[tid].pass2_secs = time_counter_get_secs(&tpg.pass2_tc);
    /*
	  tps[tid].anchor_list_secs = (double)anchor_list_usecs[tid] / 1.0e6;
          tps[tid].hit_list_secs = (double)hit_list_usecs[tid] / 1.0e6;
          tps[tid].duplicate_removal_secs = (double)duplicate_removal_usecs[tid] / 1.0e6;
          tps[tid].region_counts_secs = (double)region_counts_usecs[tid] / 1.0e6;
     */

    tps[tid].read_handle_overhead_secs = tps[tid].scan_secs
      - tps[tid].region_counts_secs - tps[tid].anchor_list_secs - tps[tid].hit_list_secs - tps[tid].duplicate_removal_secs;
  }
  f1_stats(NULL, NULL, NULL, &f1_calls_bypassed);

  fprintf(stderr, "\nStatistics:\n");

  fprintf(stderr, "%sOverall:\n", my_tab);
  fprintf(stderr, "%s%s%-24s" "%.2f seconds\n", my_tab, my_tab,
          "Load Genome Time:", (double)load_genome_usecs / 1.0e6);
  fprintf(stderr, "%s%s%-24s" "%.2f seconds\n", my_tab, my_tab,
          "Read Mapping Time:", (double)mapping_wallclock_usecs / 1.0e6);
  fprintf(stderr, "%s%s%-24s" "%s\n", my_tab, my_tab,
          "Reads per hour:", comma_integer((int)(((double)nreads/(double)mapping_wallclock_usecs) * 3600.0 * 1.0e6)));
  fprintf(stderr, "%s%s%-24s" "%s\n", my_tab, my_tab,
          "Reads per core-hour:", comma_integer((int)(((double)nreads/(double)mapping_wallclock_usecs) * 3600.0 * 1.0e6 * (1/(double)num_threads))));

  fprintf(stderr, "\n");

  int i;
  for(i = 0; i < num_threads; i++){
    total_scan_secs += tps[i].scan_secs;
    total_readparse_secs += tps[i].readparse_secs;
    total_wait_secs += time_counter_get_secs(&tpgA[i].wait_tc);

    f1_total_secs += tps[i].f1_secs;
    f1_total_invocs += tps[i].f1_invocs;
    f1_total_cells += tps[i].f1_cells;

    f2_total_secs += tps[i].f2_secs;
    f2_total_invocs += tps[i].f2_invocs;
    f2_total_cells += tps[i].f2_cells;

    if (shrimp_mode == MODE_COLOUR_SPACE) {
      fwbw_total_secs += tps[i].fwbw_secs;
      fwbw_total_invocs += tps[i].fwbw_invocs;
      fwbw_total_cells += tps[i].fwbw_cells;
    } else {
      fwbw_total_secs = 0;
      fwbw_total_invocs = 0;
      fwbw_total_cells = 0;
    }
  }
  f1_total_cellspersec = f1_total_secs == 0? 0 : (double)f1_total_cells / f1_total_secs;
  f2_total_cellspersec = f2_total_secs == 0? 0 : (double)f2_total_cells / f2_total_secs;
  fwbw_total_cellspersec = fwbw_total_secs == 0? 0 : (double)fwbw_total_cells / fwbw_total_secs;

  if (Dflag) {
    fprintf(stderr, "%sPer-Thread Stats:\n", my_tab);
    fprintf(stderr, "%s%s" "%11s %9s %9s %9s %9s %9s %9s %9s %9s %9s %9s %25s %25s %25s %9s\n", my_tab, my_tab,
            "", "ReadParse", "Scan", "Reg Cnts", "MPRegCnt", "Anch List", "Hit List", "Pass1", "Vect Hits", "Pass2", "Dup Remv",
            "Vector SW", "Scalar SW", "Post SW", "Wait");
    fprintf(stderr, "%s%s" "%11s %9s %9s %9s %9s %9s %9s %9s %9s %9s %9s %15s %9s %15s %9s %15s %9s %9s\n", my_tab, my_tab,
            "", "Time", "Time", "Time", "Time", "Time", "Time", "Time", "Time", "Time", "Time",
            "Invocs", "Time", "Invocs", "Time", "Invocs", "Time", "Time");
    fprintf(stderr, "\n");
    for(i = 0; i < num_threads; i++) {
      fprintf(stderr, "%s%s" "Thread %-4d %9.2f %9.2f %9.2f %9.2f %9.2f %9.2f %9.2f %9.2f %9.2f %9.2f %15s %9.2f %15s %9.2f %15s %9.2f %9.2f\n", my_tab, my_tab,
	  i, tps[i].readparse_secs, tps[i].scan_secs,
	  tps[i].region_counts_secs, tps[i].mp_region_counts_secs, tps[i].anchor_list_secs, tps[i].hit_list_secs,
	      tps[i].pass1_secs, tps[i].get_vector_hits_secs, tps[i].pass2_secs, tps[i].duplicate_removal_secs,
	  comma_integer(tps[i].f1_invocs), tps[i].f1_secs,
	  comma_integer(tps[i].f2_invocs), tps[i].f2_secs,
	  comma_integer(tps[i].fwbw_invocs), tps[i].fwbw_secs,
	      time_counter_get_secs(&tpgA[i].wait_tc));
    }
    for (i = 0; i < num_threads; i++) {
      fprintf (stderr, "thrd:%d anchor_list_init_size:(%.2f, %.2f) anchors_discarded:(%.2f, %.2f) big_gaps:(%.2f, %.2f)\n",
	  i, stat_get_mean(&tpgA[i].anchor_list_init_size), stat_get_sample_stddev(&tpgA[i].anchor_list_init_size),
	  stat_get_mean(&tpgA[i].n_anchors_discarded), stat_get_sample_stddev(&tpgA[i].n_anchors_discarded),
	  stat_get_mean(&tpgA[i].n_big_gaps_anchor_list), stat_get_sample_stddev(&tpgA[i].n_big_gaps_anchor_list));
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

  if (shrimp_mode == MODE_COLOUR_SPACE) {
  fprintf(stderr, "%sForward-Backward:\n", my_tab);
  fprintf(stderr, "%s%s%-24s" "%.2f seconds\n", my_tab, my_tab,
          "Run-time:", fwbw_total_secs);
  fprintf(stderr, "%s%s%-24s" "%s\n", my_tab, my_tab,
          "Invocations:", comma_integer(fwbw_total_invocs));
  fprintf(stderr, "%s%s%-24s" "%.2f million\n", my_tab, my_tab,
          "Cells Computed:", (double)fwbw_total_cells / 1.0e6);
  fprintf(stderr, "%s%s%-24s" "%.2f million\n", my_tab, my_tab,
          "Cells per Second:", fwbw_total_cellspersec / 1.0e6);
  fprintf(stderr, "\n");
  }

  fprintf(stderr, "%sMiscellaneous Totals:\n", my_tab);
  fprintf(stderr, "%s%s%-24s" "%.2f seconds\n", my_tab, my_tab,
          "Fasta Lib Time:", fasta_load_secs);
  fprintf(stderr, "%s%s%-24s" "%.2f seconds\n", my_tab, my_tab,
          "Read Load Time:", total_readparse_secs);
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
            "... with QV >= 10:",
            comma_integer(total_reads_matched_conf),
            (nreads == 0) ? 0 : ((double)total_reads_matched_conf / (double)nreads) * 100);
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
            "... with QV >= 10:",
            comma_integer(total_pairs_matched_conf),
            (nreads == 0) ? 0 : ((double)total_pairs_matched_conf / (double)(nreads/2)) * 100);
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

    if (half_paired) {
      fprintf(stderr, "\n");

      fprintf(stderr, "%s%s%-40s" "%s    (%.4f%%)\n", my_tab, my_tab,
              "Additional Reads Matched Unpaired:",
              comma_integer(total_reads_matched),
              (nreads == 0) ? 0 : ((double)total_reads_matched / (double)nreads) * 100);
      fprintf(stderr, "%s%s%-40s" "%s    (%.4f%%)\n", my_tab, my_tab,
              "... with QV >= 10:",
              comma_integer(total_reads_matched_conf),
              (nreads == 0) ? 0 : ((double)total_reads_matched_conf / (double)nreads) * 100);
      fprintf(stderr, "%s%s%-40s" "%s\n", my_tab, my_tab,
              "Total Unpaired Matches:",
              comma_integer(total_single_matches));
      fprintf(stderr, "%s%s%-40s" "%.2f\n", my_tab, my_tab,
              "Avg Matches/Unpaired Matched Read:",
              (total_reads_matched == 0) ? 0 : ((double)total_single_matches / (double)total_reads_matched));
      fprintf(stderr, "%s%s%-40s" "%s\n", my_tab, my_tab,
              "Duplicate Unpaired Matches Pruned:",
              comma_integer(total_dup_single_matches));
    }
  }

  fprintf(stderr, "\n");

  fprintf(stderr, "%sMemory usage:\n", my_tab);
  fprintf(stderr, "%s%s%-24s" "%s\n", my_tab, my_tab,
          "Genomemap:",
          comma_integer(count_get_count(&mem_genomemap)));

  if (Xflag) {
    print_insert_histogram();
  }

  free(tps);
  free(tpgA);
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
	  "   -n/--cmw-mode        Match Mode                    (default: unpaired:%d paired:%d)\n",
	  DEF_MATCH_MODE_UNPAIRED, DEF_MATCH_MODE_PAIRED);
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
	  shrimp_mode == MODE_LETTER_SPACE? DEF_LS_MATCH_SCORE : DEF_CS_MATCH_SCORE);
  fprintf(stderr,
	  "   -i/--mismatch        SW Mismatch Score             (default: %d)\n",
	  shrimp_mode == MODE_LETTER_SPACE? DEF_LS_MISMATCH_SCORE : DEF_CS_MISMATCH_SCORE);
  fprintf(stderr,
	  "   -g/--open-r          SW Gap Open Score (Reference) (default: %d)\n",
	  shrimp_mode == MODE_LETTER_SPACE? DEF_LS_A_GAP_OPEN : DEF_CS_A_GAP_OPEN);
  fprintf(stderr,
	  "   -q/--open-q          SW Gap Open Score (Query)     (default: %d)\n",
	  shrimp_mode == MODE_LETTER_SPACE? DEF_LS_B_GAP_OPEN : DEF_CS_B_GAP_OPEN);
  fprintf(stderr,
	  "   -e/--ext-r           SW Gap Extend Score(Reference)(default: %d)\n",
	  shrimp_mode == MODE_LETTER_SPACE? DEF_LS_A_GAP_EXTEND : DEF_CS_A_GAP_EXTEND);
  fprintf(stderr,
	  "   -f/--ext-q           SW Gap Extend Score (Query)   (default: %d)\n",
	  shrimp_mode == MODE_LETTER_SPACE? DEF_LS_B_GAP_EXTEND : DEF_CS_B_GAP_EXTEND);
  if (shrimp_mode == MODE_COLOUR_SPACE) {
  fprintf(stderr,
	  "   -x/--crossover       SW Crossover Score            (default: %d)\n",
	  DEF_CS_XOVER_SCORE);
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
	  "      --single-best-mapping Report only the best mapping(s), this is not strata (see README)\n");
  fprintf(stderr,
          "      --all-contigs     Report a maximum of 1 mapping for each read.\n");
  fprintf(stderr,
          "      --no-mapping-qualities Do not compute mapping qualities\n");
  fprintf(stderr,
	  "      --insert-size-dist Specifies the mean and stddev of the insert sizes\n");
  fprintf(stderr,
	  "      --no-improper-mappings (see README)\n");
  if (full_usage) {
  fprintf(stderr,
          "      --trim-front      Trim front of reads by this amount\n");
  fprintf(stderr,
          "      --trim-end        Trim end of reads by this amount\n");
  fprintf(stderr,
          "      --trim-first      Trim only first read in pair\n");
  fprintf(stderr,
          "      --trim-second     Trim only second read in pair\n");
  fprintf(stderr,
          "      --min-avg-qv      The minimum average quality value of a read\n");
  fprintf(stderr,
          "      --progress        Display a progress line each <value> reads. (default %d)\n",progress);
  fprintf(stderr,
 	  "      --save-mmap       Save genome projection to shared memory\n");
  fprintf(stderr,
          "      --load-mmap       Load genome projection from shared memory\n");
  fprintf(stderr,
          "      --indel-taboo-len Prevent indels from starting or ending in the tail\n");
  fprintf(stderr,
          "      --shrimp-format   Output mappings in SHRiMP format (default: %s)\n",Eflag ? "disabled" : "enabled");
  fprintf(stderr,
          "      --qv-offset       (see README)\n");
  fprintf(stderr,
          "      --sam-header-hd   (see README)\n");
  fprintf(stderr,
          "      --sam-header-sq   (see README)\n");
  fprintf(stderr,
          "      --sam-header-rg   (see README)\n");
  fprintf(stderr,
          "      --sam-header-pg   (see README)\n");
  fprintf(stderr,
          "      --no-autodetect-input (see README)\n");
  }
  fprintf(stderr, "\n");
  fprintf(stderr, "Options:\n");

  fprintf(stderr,
	  "   -U/--ungapped        Perform Ungapped Alignment    (default: %s)\n", gapless_sw ? "enabled" : "disabled");
  fprintf(stderr,
          "      --global          Perform full global alignment (default: %s)\n", Gflag ? "enabled" : "disabled");
  fprintf(stderr,
          "      --local           Perform local alignment       (default: %s)\n", Gflag ? "disabled" : "enabled");
  if (shrimp_mode == MODE_COLOUR_SPACE) {
  fprintf(stderr,
          "      --bfast           Try to align like bfast       (default: %s)\n", Bflag ? "enabled" : "disabled");
  } else {
  fprintf(stderr,
          "      --trim-illumina   Trim trailing B qual values   (default: %s)\n", trim_illumina ? "enabled" : "disabled");
  }
  fprintf(stderr,
	  "   -C/--negative        Negative Strand Aln. Only     (default: %s)\n", Cflag ? "enabled" : "disabled");
  fprintf(stderr,
	  "   -F/--positive        Positive Strand Aln. Only     (default: %s)\n", Fflag ? "enabled" : "disabled");
  fprintf(stderr,
	  "   -P/--pretty          Pretty Print Alignments       (default: %s)\n", Pflag ? "enabled" : "disabled");
  fprintf(stderr,
	  "   -E/--sam             Output SAM Format             (default: %s)\n", Eflag ? "enabled" : "disabled");
  fprintf(stderr,
  	  "   -Q/--fastq           Reads are in fastq format     (default: %s)\n", Qflag ? "enabled" : "disabled");
  if (full_usage) {
  fprintf(stderr,
	  "   -R/--print-reads     Print Reads in Output         (default: %s)\n", Rflag ? "enabled" : "disabled");
 // fprintf(stderr,
//	  "    -T    (does nothing since default) Reverse Tie-break on Negative Strand          (default: enabled)\n");
  fprintf(stderr,
	  "   -t/--tiebreak-off    Disable Reverse Tie-break\n");
  fprintf(stderr,
	  "                                  on Negative Strand  (default: %s)\n",Tflag ? "enabled" : "disabled");
  fprintf(stderr,
	  "   -X/--isize-hist      Print Insert Size Histogram   (default: %s)\n", Xflag ? "enabled" : "disabled");
  fprintf(stderr,
	  "   -Y/--proj-hist       Print Genome Proj. Histogram  (default: %s)\n", Yflag ? "enabled" : "disabled");
  fprintf(stderr,
	  "   -Z/--bypass-off      Disable Cache Bypass for SW\n");
  fprintf(stderr,
	  "                                    Vector Calls      (default: %s)\n", hash_filter_calls ? "enabled" : "disabled");
  fprintf(stderr,
	  "   -H/--spaced-kmers    Hash Spaced Kmers in Genome\n");
  fprintf(stderr,
	  "                                    Projection        (default: %s)\n", Hflag ? "enabled" : "disabled");
  fprintf(stderr,
	  "   -D/--thread-stats    Individual Thread Statistics  (default: %s)\n", Dflag ? "enabled" : "disabled");
  fprintf(stderr,
	  "   -V/--trim-off        Disable Automatic Genome\n");
  fprintf(stderr,
	  "                                 Index Trimming       (default: %s)\n", Vflag ? "enabled" : "disabled");
  fprintf(stderr,
	  "      --pr-xover        Set the Default Crossover Probability\n");
  fprintf(stderr,
	  "                                                      (default: %.1e)\n", pr_xover);
  }
  fprintf(stderr,
	  "      --sam-unaligned   Unaligned reads in SAM output (default: %s)\n", sam_unaligned ? "enabled" : "disabled");
  fprintf(stderr,
	  "      --half-paired     Output half mapped read pairs (default: %s)\n", half_paired ? "enabled" : "disabled");
  fprintf(stderr,
	  "      --strata          Print only the best scoring hits\n");
  fprintf(stderr,
	  "   -?/--help            Full List of Parameters and Options\n");
  exit(1);
}


static void
print_pairing_options(struct pairing_options * options)
{
  char buff[2][100];

  fprintf(stderr, "[\n");
  fprintf(stderr, "  pairing:%s\n", pair_mode_string[options->pair_mode]);
  fprintf(stderr, "  min_insert:%d\n", options->min_insert_size);
  fprintf(stderr, "  max_insert:%d\n", options->max_insert_size);
  fprintf(stderr, "  pass1_num_outputs:%d\n", options->pass1_num_outputs);
  fprintf(stderr, "  pass1_threshold:%s\n", thres_to_buff(buff[0], &options->pass1_threshold));
  fprintf(stderr, "  pass2_num_outputs:%d\n", options->pass2_num_outputs);
  fprintf(stderr, "  pass2_threshold:%s\n", thres_to_buff(buff[0], &options->pass2_threshold));
  fprintf(stderr, "  strata:%s\n", bool_to_buff(buff[0], &options->strata));
  fprintf(stderr, "  save_outputs:%s\n", bool_to_buff(buff[0], &options->save_outputs));
  fprintf(stderr, "  stop_count:%d\n", options->stop_count);
  fprintf(stderr, "  stop_threshold:%s\n", thres_to_buff(buff[0], &options->stop_threshold));
  fprintf(stderr, "]\n");
}


static void
print_read_mapping_options(struct read_mapping_options_t * options, bool is_paired)
{
  char buff[2][100];

  fprintf(stderr, "[\n");

  fprintf(stderr, "  regions:\n");
  //fprintf(stderr, "  [\n");
  fprintf(stderr, "    recompute:%s\n", bool_to_buff(buff[0], &options->regions.recompute));
  //if (!options->regions.recompute)
  //fprintf(stderr, "  ]\n");
  //else
  //  fprintf(stderr, ", min_seed:%d, max_seed:%d]\n",
  //    options->regions.min_seed, options->regions.max_seed);

  fprintf(stderr, "  anchor_list:\n");
  //fprintf(stderr, "  [\n");
  fprintf(stderr, "    recompute:%s\n", bool_to_buff(buff[0], &options->anchor_list.recompute));
  if (options->anchor_list.recompute) {
    fprintf(stderr, "    collapse:%s\n", bool_to_buff(buff[0], &options->anchor_list.collapse));
    fprintf(stderr, "    use_region_counts:%s\n", bool_to_buff(buff[0], &options->anchor_list.use_region_counts));
    fprintf(stderr, "    use_mp_region_counts:%d\n", options->anchor_list.use_mp_region_counts);
  }
  //fprintf(stderr, "  ]\n");

  fprintf(stderr, "  hit_list:\n");
  //fprintf(stderr, "  [\n");
  fprintf(stderr, "    recompute:%s\n", bool_to_buff(buff[0], &options->hit_list.recompute));
  if (options->hit_list.recompute) {
    fprintf(stderr, "    gapless:%s\n", bool_to_buff(buff[0], &options->hit_list.gapless));
    fprintf(stderr, "    match_mode:%d\n", options->hit_list.match_mode);
    fprintf(stderr, "    threshold:%s\n", thres_to_buff(buff[0], &options->hit_list.threshold));
  }
  //fprintf(stderr, "  ]\n");

  fprintf(stderr, "  pass1:\n");
  //fprintf(stderr, "  [\n");
  fprintf(stderr, "    recompute:%s\n", bool_to_buff(buff[0], &options->pass1.recompute));
  if (options->pass1.recompute) {
    fprintf(stderr, "    threshold:%s\n", thres_to_buff(buff[0], &options->pass1.threshold));
    fprintf(stderr, "    window_overlap:%s\n", thres_to_buff(buff[1], &options->pass1.window_overlap));
    fprintf(stderr, "    min_matches:%d\n", options->pass1.min_matches);
    fprintf(stderr, "    gapless:%s\n", bool_to_buff(buff[0], &options->pass1.gapless));
    if (is_paired) {
      fprintf(stderr, "    only_paired:%s\n", bool_to_buff(buff[0], &options->pass1.only_paired));
    }
    if (!is_paired) {
      fprintf(stderr, "    num_outputs:%d\n", options->pass1.num_outputs);
    }
  }
  //fprintf(stderr, "  ]\n");

  fprintf(stderr, "  pass2:\n");
  //fprintf(stderr, "  [\n");
  fprintf(stderr, "    threshold:%s\n", thres_to_buff(buff[0], &options->pass2.threshold));
  if (!is_paired) {
    fprintf(stderr, "    strata:%s\n", bool_to_buff(buff[0], &options->pass2.strata));
    fprintf(stderr, "    save_outputs:%s\n", bool_to_buff(buff[0], &options->pass2.save_outputs));
    fprintf(stderr, "    num_outputs:%d\n", options->pass2.num_outputs);
  }
  //fprintf(stderr, "  ]\n");

  if (!is_paired) {
    fprintf(stderr, "  stop:\n");
    //fprintf(stderr, "  [\n");
    fprintf(stderr, "    stop_count:%d\n", options->pass2.stop_count);
    if (options->pass2.stop_count > 0) {
      fprintf(stderr, "    stop_threshold:%s\n", thres_to_buff(buff[0], &options->pass2.stop_threshold));
    }
    //fprintf(stderr, "  ]\n");
  }

  fprintf(stderr, "]\n");
}


static void
print_settings() {
  char buff[100];
  static char const my_tab[] = "    ";
  int sn, i;

  fprintf(stderr, "Settings:\n");

  // Seeds
  fprintf(stderr, "%s%-40s%s (%d/%d)\n", my_tab,
	  (n_seeds == 1) ? "Spaced Seed (weight/span)" : "Spaced Seeds (weight/span)",
	  seed_to_string(0), seed[0].weight, seed[0].span);
  for (sn = 1; sn < n_seeds; sn++) {
    fprintf(stderr, "%s%-40s%s (%d/%d)\n", my_tab, "",
	    seed_to_string(sn), seed[sn].weight, seed[sn].span);
  }

  // Global settings
  fprintf(stderr, "\n");
  fprintf(stderr, "%s%-40s%d\n", my_tab, "Number of threads:", num_threads);
  fprintf(stderr, "%s%-40s%d\n", my_tab, "Thread chunk size:", chunk_size);
  fprintf(stderr, "%s%-40s%s\n", my_tab, "Window length:", thres_to_buff(buff, &window_len));

  fprintf(stderr, "%s%-40s%s\n", my_tab, "Hash filter calls:", hash_filter_calls? "yes" : "no");
  fprintf(stderr, "%s%-40s%d%s\n", my_tab, "Anchor width:", anchor_width,
	  anchor_width == -1? " (disabled)" : "");
  fprintf(stderr, "%s%-40s%d%s\n", my_tab, "Indel taboo Len:", indel_taboo_len,
	  indel_taboo_len == 0? " (disabled)" : "");
  if (list_cutoff < DEF_LIST_CUTOFF) {
  fprintf(stderr, "%s%-40s%u\n", my_tab, "Index list cutoff length:", list_cutoff);
  }
  fprintf(stderr, "%s%-40s%s\n", my_tab, "Gapless mode:", gapless_sw? "yes" : "no");
  fprintf(stderr, "%s%-40s%s\n", my_tab, "Global alignment:", Gflag? "yes" : "no");
  fprintf(stderr, "%s%-40s%s\n", my_tab, "Region filter:", use_regions? "yes" : "no");
  if (use_regions) {
  fprintf(stderr, "%s%-40s%d\n", my_tab, "Region size:", (1 << region_bits));
  fprintf(stderr, "%s%-40s%d\n", my_tab, "Region overlap:", region_overlap);
  }
  if (Qflag) {
  fprintf(stderr, "%s%-40s%s\n", my_tab, "Ignore QVs:", ignore_qvs? "yes" : "no");
  }
  if (Qflag && !ignore_qvs) {
  fprintf(stderr, "%s%-40s%d%s\n", my_tab, "Minimum average qv:", min_avg_qv, min_avg_qv < 0? " (none)" : "");
  //fprintf(stderr, "%s%-40sPHRED+%d\n", my_tab, "QV input encoding:", qual_delta);
  }
  fprintf(stderr, "%s%-40s%s\n", my_tab, "Compute mapping qualities:", compute_mapping_qualities? "yes" : "no");
  if (compute_mapping_qualities)
  {
  fprintf(stderr, "%s%-40s%s\n", my_tab, "All contigs:", all_contigs? "yes" : "no");
  fprintf(stderr, "%s%-40s%s\n", my_tab, "Single best mapping:", single_best_mapping? "yes" : "no");
  }
  //fprintf(stderr, "%s%-40s%s\n", my_tab, "Hack:", hack? "yes" : "no");

  // Scores
  fprintf(stderr, "\n");
  fprintf(stderr, "%s%-40s%-10d\n", my_tab, "SW Match Score:", match_score);
  fprintf(stderr, "%s%-40s%-10d\t[%.1e]\n", my_tab, "SW Mismatch Score [Prob]:", mismatch_score, pr_mismatch);
  fprintf(stderr, "%s%-40s%-10d\t[%.1e]\n", my_tab, "SW Del Open Score [Prob]:", a_gap_open_score, pr_del_open);
  fprintf(stderr, "%s%-40s%-10d\t[%.1e]\n", my_tab, "SW Ins Open Score [Prob]:", b_gap_open_score, pr_ins_open);
  fprintf(stderr, "%s%-40s%-10d\t[%.1e]\n", my_tab, "SW Del Extend Score [Prob]:", a_gap_extend_score, pr_del_extend);
  fprintf(stderr, "%s%-40s%-10d\t[%.1e]\n", my_tab, "SW Ins Extend Score [Prob]:", b_gap_extend_score, pr_ins_extend);
  if (shrimp_mode == MODE_COLOUR_SPACE) {
    fprintf(stderr, "%s%-40s%-10d\t[%.1e]\n", my_tab, "SW Crossover Score [Prob]:", crossover_score, pr_xover);
  }

  fprintf(stderr, "\n");
  if (n_paired_mapping_options > 0) { // paired mapping
    for (i = 0; i < n_paired_mapping_options; i++) {
      fprintf(stderr, "Paired mapping options, set [%d]\n", i);
      print_pairing_options(&paired_mapping_options[i].pairing);
      print_read_mapping_options(&paired_mapping_options[i].read[0], true);
      print_read_mapping_options(&paired_mapping_options[i].read[1], true);
    }

    if (n_unpaired_mapping_options[0] > 0) {
      fprintf(stderr, "\n");
      for (i = 0; i < n_unpaired_mapping_options[0]; i++) {
	fprintf(stderr, "Unpaired mapping options for first read in a pair, set [%d]\n", i);
	print_read_mapping_options(&unpaired_mapping_options[0][i], false);
      }
    }

    if (n_unpaired_mapping_options[1] > 0) {
      fprintf(stderr, "\n");
      for (i = 0; i < n_unpaired_mapping_options[1]; i++) {
	fprintf(stderr, "Unpaired mapping options for second read in a pair, set [%d]\n", i);
	print_read_mapping_options(&unpaired_mapping_options[1][i], false);
      }
    }
  } else {
    for (i = 0; i < n_unpaired_mapping_options[0]; i++) {
      fprintf(stderr, "Unpaired mapping options, set [%d]\n", i);
      print_read_mapping_options(&unpaired_mapping_options[0][i], false);
    }
  }
  fprintf(stderr, "\n");

  return;

  fprintf(stderr, "%s%-40s%d\n", my_tab, "Number of Outputs per Read:", num_outputs);
  fprintf(stderr, "%s%-40s%d\n", my_tab, "Window Generation Mode:", match_mode);

  if (IS_ABSOLUTE(window_overlap)) {
    fprintf(stderr, "%s%-40s%u\n", my_tab, "Window Overlap Length:", (uint)-window_overlap);
  } else {
    fprintf(stderr, "%s%-40s%.02f%%\n", my_tab, "Window Overlap Length:", window_overlap);
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
    match_mode = 1;
    window_len = 100.0;
    Gflag = false;
    compute_mapping_qualities = false;

    return 1;
  } else {
    return 0;
  }
}


static void
get_pair_mode(char * c, int * pair_mode)
{
  if (c == NULL) {
    fprintf(stderr, "error: invalid pair mode\n");
    exit(1);
  }
  if (!strcmp(c, "none")) {
    *pair_mode = PAIR_NONE;
  } else if (!strcmp(c, "opp-in")) {
    *pair_mode = PAIR_OPP_IN;
  } else if (!strcmp(c, "opp-out")) {
    *pair_mode = PAIR_OPP_OUT;
  } else if (!strcmp(c, "col-fw")) {
    *pair_mode = PAIR_COL_FW;
  } else if (!strcmp(c, "col-bw")) {
    *pair_mode = PAIR_COL_BW;
  } else {
    fprintf(stderr, "error: unrecognized pair mode (%s)\n", c);
    exit(1);
  }
}


static void
get_int(char * c, int * d)
{
  if (c == NULL) {
    fprintf(stderr, "error: invalid integer\n");
    exit(1);
  }
  *d = atoi(c);
}


static void
get_threshold(char * c, double * t)
{
  if (c == NULL) {
    fprintf(stderr, "error: invalid threshold\n");
    exit(1);
  }
  *t = atof(c);
  if (*t < 0.0) {
    fprintf(stderr, "error: invalid threshold [%s]\n", c);
    exit(1);
  }
  if (strcspn(c, "%.") == strlen(c))
    *t = -(*t);	//absol.
}


static void
get_bool(char * c, bool * b)
{
  if (!strcmp(c, "true") || !strcmp(c, "1")) {
    *b = true;
  } else if (!strcmp(c, "false") || !strcmp(c, "0")) {
    *b = false;
  } else {
    fprintf(stderr, "error: invalid bool\n");
    exit(1);
  }
}


static void
get_pairing_options(char * c, struct pairing_options * options)
{
  char * p;
  if (c == NULL) {
    fprintf(stderr, "error: invalid pairing options\n");
    exit(1);
  }
  p = strtok(c, ",");
  get_pair_mode(p, &options->pair_mode);
  p = strtok(NULL, ",");
  get_int(p, &options->min_insert_size);
  p = strtok(NULL, ",");
  get_int(p, &options->max_insert_size);
  p = strtok(NULL, ",");
  get_int(p, &options->pass1_num_outputs);
  p = strtok(NULL, ",");
  get_threshold(p, &options->pass1_threshold);
  p = strtok(NULL, ",");
  get_int(p, &options->pass2_num_outputs);
  p = strtok(NULL, ",");
  get_threshold(p, &options->pass2_threshold);
  p = strtok(NULL, ",");
  get_int(p, &options->stop_count);
  p = strtok(NULL, ",");
  get_threshold(p, &options->stop_threshold);
  p = strtok(NULL, ",");
  get_bool(p, &options->strata);  
  p = strtok(NULL, ",");
  get_bool(p, &options->save_outputs);
}


static void
get_read_mapping_options(char * c, struct read_mapping_options_t * options, bool is_paired)
{
  char * p, * q;
  char * save_ptr;
  if (c == NULL) {
    fprintf(stderr, "error: invalid read mapping options\n");
    exit(1);
  }
  fprintf(stderr, "parsing read_mapping_options [%s] (%s)\n", c, is_paired? "paired" : "unpaired");
  // regions
  q = strtok_r(c, "/", &save_ptr);
  logit(0, "parsing region options: %s", q);
  p = strtok(q, ",");
  get_bool(p, &options->regions.recompute);
  //if (options->regions.recompute) {
  //p = strtok(NULL, ",");
  //get_int(p, &options->regions.min_seed);
  //p = strtok(NULL, ",");
  //get_int(p, &options->regions.max_seed);
  //}
  // anchor_list
  q = strtok_r(NULL, "/", &save_ptr);
  logit(0, "parsing anchor_list options: %s", q);
  p = strtok(q, ",");
  get_bool(p, &options->anchor_list.recompute);
  if (options->anchor_list.recompute) {
    p = strtok(NULL, ",");
    get_bool(p, &options->anchor_list.collapse);
    p = strtok(NULL, ",");
    get_bool(p, &options->anchor_list.use_region_counts);
    if (is_paired) {
      p = strtok(NULL, ",");
      get_int(p, &options->anchor_list.use_mp_region_counts);
    }
  }
  // hit_list
  q = strtok_r(NULL, "/", &save_ptr);
  logit(0, "parsing hit_list options: %s", q);
  p = strtok(q, ",");
  get_bool(p, &options->hit_list.recompute);
  if (options->hit_list.recompute) {
    p = strtok(NULL, ",");
    get_bool(p, &options->hit_list.gapless);
    p = strtok(NULL, ",");
    get_int(p, &options->hit_list.match_mode);
    p = strtok(NULL, ",");
    get_threshold(p, &options->hit_list.threshold);
  }
  // pass1
  q = strtok_r(NULL, "/", &save_ptr);
  logit(0, "parsing pass1 options: %s", q);
  p = strtok(q, ",");
  get_bool(p, &options->pass1.recompute);
  if (options->pass1.recompute) {
    p = strtok(NULL, ",");
    get_threshold(p, &options->pass1.threshold);
    p = strtok(NULL, ",");
    get_threshold(p, &options->pass1.window_overlap);
    p = strtok(NULL, ",");
    get_int(p, &options->pass1.min_matches);
    p = strtok(NULL, ",");
    get_bool(p, &options->pass1.gapless);
    if (is_paired) {
      p = strtok(NULL, ",");
      get_bool(p, &options->pass1.only_paired);
    }
    if (!is_paired) {
      p = strtok(NULL, ",");
      get_int(p, &options->pass1.num_outputs);
    }
  }
  // pass2
  q = strtok_r(NULL, "/", &save_ptr);
  logit(0, "parsing pass2 options: %s", q);
  p = strtok(q, ",");
  get_threshold(p, &options->pass2.threshold);
  if (!is_paired) {
    p = strtok(NULL, ",");
    get_bool(p, &options->pass2.strata);
    p = strtok(NULL, ",");
    get_bool(p, &options->pass2.save_outputs);
    p = strtok(NULL, ",");
    get_int(p, &options->pass2.num_outputs);
  }
  // stop
  if (!is_paired) {
    q = strtok_r(NULL, "/", &save_ptr);
    logit(0, "parsing stop options: %s", q);
    p = strtok(q, ",");
    get_int(p, &options->pass2.stop_count);
    if (options->pass2.stop_count > 0) {
      p = strtok(NULL, ",");
      get_threshold(p, &options->pass2.stop_threshold);
    }
  }
}


int main(int argc, char **argv){
	char **genome_files = NULL;
	int ngenome_files = 0;

	char *progname = argv[0];
	char const * optstr = NULL;
	char *c, *save_c;
	int ch;
	int i, sn, cn;

	llint before;

	bool a_gap_open_set, b_gap_open_set;
	bool a_gap_extend_set, b_gap_extend_set;
	bool match_score_set, mismatch_score_set, xover_score_set;
	bool match_mode_set = false;
	bool qual_delta_set = false;

	fasta_t fasta = NULL, left_fasta = NULL, right_fasta = NULL;

	my_alloc_init(64l*1024l*1024l*1024l, 64l*1024l*1024l*1024l);

	shrimp_args.argc=argc;
	shrimp_args.argv=argv;
	set_mode_from_argv(argv, &shrimp_mode);

	a_gap_open_set = b_gap_open_set = a_gap_extend_set = b_gap_extend_set = false;
	match_score_set = mismatch_score_set = xover_score_set = false;
	if (shrimp_mode == MODE_COLOUR_SPACE) {
	  match_score = DEF_CS_MATCH_SCORE;
	  mismatch_score = DEF_CS_MISMATCH_SCORE;
	  a_gap_open_score = DEF_CS_A_GAP_OPEN;
	  b_gap_open_score = DEF_CS_B_GAP_OPEN;
	  a_gap_extend_score = DEF_CS_A_GAP_EXTEND;
	  b_gap_extend_score = DEF_CS_B_GAP_EXTEND;
	}

	fasta_reset_stats();

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

	
	//get a copy of the command line
	size_t size_of_command_line=0;
	for (i=0; i<argc; i++) {
		size_of_command_line+=strlen(argv[i])+1;
	}
	size_of_command_line++;
	char command_line[size_of_command_line];
	size_t offset=0;
	for (i=0; i<argc; i++) {
		offset+=sprintf(command_line+offset,"%s",argv[i]);
		if (i+1!=argc) {
			offset+=sprintf(command_line+offset," ");
		}
	}
	assert(offset+1<=size_of_command_line);

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
			logit(0, "as of v2.2.0, --global is on by default");
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
			half_paired = false;
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
			match_mode = atoi(optarg);
			match_mode_set = true;
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
			num_tmp_outputs = 20 + num_outputs;
			break;
		case 'm':
			match_score = atoi(optarg);
			match_score_set = true;
			break;
		case 'i':
			mismatch_score = atoi(optarg);
			mismatch_score_set = true;
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
			xover_score_set = true;
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
			Eflag = false;
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
		  c = strtok(optarg, ",");
		  if (c == NULL)
		    crash(1, 0, "format for insert sizes is \"-I 200,1000\"\n");

		  min_insert_size = atoi(c);
		  if (min_insert_size < 0) {
		    logit(0, "insert sizes must be nonnegative; check README. resetting min_insert_size from [%s] to 0", optarg);
		    min_insert_size = 0;
		  }

		  c = strtok(NULL, ",");
		  if (c == NULL)
		    crash(1, 0, "format for insert sizes is \"-I 200,1000\"\n");

		  max_insert_size = atoi(c);
		  if (max_insert_size < 0) {
		    logit(0, "insert sizes must be nonnegative; check README. resetting max_insert_size from [%s] to 0", optarg);
		    max_insert_size = 0;
		  }

		  if (min_insert_size > max_insert_size)
		    crash(1, 0, "invalid insert sizes (min:%d,max:%d)\n", min_insert_size, max_insert_size);
		  break;
		case 'E': // on by default; accept option for bw compatibility
		  logit(0, "as of v2.2.0, -E/--sam is on by default");
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
			  match_mode_set = true;
                          c = strtok(NULL, ",");
                        } while (c != NULL);
                        break;
		case 21:
			trim_front=atoi(optarg);
			if (shrimp_mode == MODE_COLOUR_SPACE) {
				fprintf(stderr,"--trim-front cannot be used in colour space mode!\n");
				exit(1);
			}
			if (trim_front<0) {
				fprintf(stderr,"--trim-front value must be positive\n");
				exit(1);
			}
			if (trim_front>0)
				trim=true;
			break;
		case 22:
			trim_end=atoi(optarg);
			if (trim_end<0) {
				fprintf(stderr,"--trim-end value must be positive\n");
				exit(1);
			}
			if (trim_end>0)
				trim=true;
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
			c = strtok(optarg, ",");
			if (c == NULL)
			  crash(1, 0, "argmuent for insert-size-dist should be \"mean,stddev\" [%s]", optarg);
			insert_size_mean = atof(c);
			c = strtok(NULL, ",");
			if (c == NULL)
			  crash(1, 0, "argmuent for insert-size-dist should be \"mean,stddev\" [%s]", optarg);
			insert_size_stddev = atof(c);
			break;
		case 26:
		  use_regions = !use_regions;
		  break;
		case 27:
		  region_overlap = atoi(optarg);
		  if (region_overlap < 0) {
		    fprintf(stderr, "region overlap must be non-negative!\n");
		    exit(1);
		  }
		  break;
		case 28:
		  if (n_unpaired_mapping_options[0] > 0 || n_unpaired_mapping_options[1] > 0) {
		    fprintf(stderr, "warning: unpaired mapping options set before paired mapping options! the latter take precedence.\n");
		    half_paired = true;
		  }
		  n_paired_mapping_options++;
		  paired_mapping_options = (struct readpair_mapping_options_t *)
		    //xrealloc(paired_mapping_options, n_paired_mapping_options * sizeof(paired_mapping_options[0]));
		    my_realloc(paired_mapping_options, n_paired_mapping_options * sizeof(paired_mapping_options[0]), (n_paired_mapping_options - 1) * sizeof(paired_mapping_options[0]),
			       &mem_small, "paired_mapping_options");
		  c = strtok_r(optarg, ";", &save_c);
		  get_pairing_options(c, &paired_mapping_options[n_paired_mapping_options - 1].pairing);
		  c = strtok_r(NULL, ";", &save_c);
		  get_read_mapping_options(c, &paired_mapping_options[n_paired_mapping_options - 1].read[0], true);
		  c = strtok_r(NULL, ";", &save_c);
		  get_read_mapping_options(c, &paired_mapping_options[n_paired_mapping_options - 1].read[1], true);
		  // HACK SETTINGS
		  pair_mode = paired_mapping_options[0].pairing.pair_mode;
		  break;
		case 29:
		  int nip;
		  c = strtok(optarg, ";");
		  if (c == NULL || (*c != '0' && *c != '1')) {
		    fprintf(stderr, "error: invalid unpaired mapping options:[%s]\n", optarg);
		    exit(1);
		  }
		  if (n_paired_mapping_options > 0)
		    half_paired = true;
		  nip = (*c == '0'? 0 : 1);
		  n_unpaired_mapping_options[nip]++;
		  unpaired_mapping_options[nip] = (struct read_mapping_options_t *)
		    //xrealloc(unpaired_mapping_options[nip], n_unpaired_mapping_options[nip] * sizeof(unpaired_mapping_options[nip][0]));
		    my_realloc(unpaired_mapping_options[nip], n_unpaired_mapping_options[nip] * sizeof(unpaired_mapping_options[nip][0]), (n_unpaired_mapping_options[nip] - 1) * sizeof(unpaired_mapping_options[nip][0]),
			       &mem_small, "unpaired_mapping_options[%d]", nip);
		  c = strtok(NULL, ";");
		  get_read_mapping_options(c, &unpaired_mapping_options[nip][n_unpaired_mapping_options[nip] - 1], false);
		  break;
		case 30:
		  min_avg_qv = atoi(optarg);
		  if (min_avg_qv < -2 || min_avg_qv > 40) {
		    fprintf(stderr, "error: invalid minimum average quality value (%s)\n", optarg);
		  }
		  break;
		case 31:
		  extra_sam_fields = true;
		  break;
		case 32:
		  region_bits = atoi(optarg);
		  if (region_bits < 8 || region_bits > 20) {
		    crash(1, 0, "invalid number of region bits: %s; must be between 8 and 20", optarg);
		  }
		  n_regions = (1 << (32 - region_bits));
		  break;
		case 33:
		  progress = atoi(optarg);
		  break;
		case 34:
		  save_mmap = optarg;
		  break;
		case 35:
		  load_mmap = optarg;
		  break;
		case 36:
		  indel_taboo_len = atoi(optarg);
		  if (indel_taboo_len < 0)
		    crash(1, 0, "invalid indel taboo len: [%s]", optarg);
		  break;
		case 37:
		  single_best_mapping = true;
		  break;
		case 38:
		  all_contigs = true;
		  break;
		case 39:
		  compute_mapping_qualities = false;
		  break;
		case 40:
		  Eflag = false;
		  break;
		case 41: // half-paired: accept option for bw compatibility
		  logit(0, "as of v2.2.0, --half-paired is on by default");
		  half_paired = true;
		  break;
		case 42: // no-improper-mappings
		  improper_mappings = false;
		  break;
		case 43: // qual value offset
		  qual_delta = atoi(optarg);
		  qual_delta_set = true;
		  break;
		case 44: // sam-header-hd
		  sam_header_hd = fopen(optarg, "r");
		  if (sam_header_hd == NULL)
		    crash(1, 1, "cannot open sam header file with HD lines [%s]", optarg);
		  break;
		case 45: // sam-header-sq
		  sam_header_sq = fopen(optarg, "r");
		  if (sam_header_sq == NULL)
		    crash(1, 1, "cannot open sam header file with SQ lines [%s]", optarg);
		  break;
		case 46: // sam-header-rg
		  sam_header_rg = fopen(optarg, "r");
		  if (sam_header_rg == NULL)
		    crash(1, 1, "cannot open sam header file with RG lines [%s]", optarg);
		  break;
		case 47: // sam-header-pg
		  sam_header_pg = fopen(optarg, "r");
		  if (sam_header_pg == NULL)
		    crash(1, 1, "cannot open sam header file with PG lines [%s]", optarg);
		  break;
		case 48: // no-autodetect-input
		  autodetect_input = false;
		  break;
		case 3:
		  trim_illumina=true;
		  break;
		case 123:
		  no_qv_check=true;
		  break;
		case 124: // local alignment
		  Gflag = false;
		  break;
		case 125: // --ignore-qvs
		  ignore_qvs = true;
		  break;
		case 126:
		  pr_xover = atof(optarg);
		  break;
		default:
			usage(progname, false);
		}
	}

	argc -= optind;
	argv += optind;

	if ((pair_mode != PAIR_NONE || !single_reads_file) && (chunk_size % 2) != 0) {
	  logit(0, "in paired mode or if using options -1 and -2, the thread chunk size must be even; adjusting it to [%d]", chunk_size + 1);
	  chunk_size++;
	}

	if (!Gflag && compute_mapping_qualities) {
	  logit(0, "mapping qualities are not available in local alignment mode; disabling them");
	  compute_mapping_qualities = false;
	}

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
	if (!match_mode_set) {
	  match_mode = (pair_mode == PAIR_NONE? DEF_MATCH_MODE_UNPAIRED : DEF_MATCH_MODE_PAIRED);
	}
	if (pair_mode == PAIR_NONE && (!trim_first || !trim_second)) {
		fprintf(stderr,"error: cannot use --trim-first or --trim-second in unpaired mode\n");
		usage(progname,false);
	} 

	// set up insert size histogram
	if (Xflag && pair_mode == PAIR_NONE) {
	  fprintf(stderr, "warning: insert histogram not available in unpaired mode; ignoring\n");
	  Xflag = false;
	}
	if (pair_mode != PAIR_NONE) {
	  insert_histogram_bucket_size = ceil_div(max_insert_size - min_insert_size + 1, 100);
	  for (i = 0; i < 100; i++) {
	    insert_histogram[i] = 1;
	    insert_histogram_load += insert_histogram[i];
	  }
	}

	if(load_file != NULL && n_seeds != 0){
	  fprintf(stderr,"error: cannot specify seeds when loading genome map\n");
	  usage(progname,false);
	}

	if (n_seeds == 0 && load_file == NULL && load_mmap == NULL) {
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

	if (load_file != NULL && (save_file != NULL || save_mmap != NULL))
	  { // args: none
	    if (argc != 0) {
	      fprintf(stderr, "error: when using both -L and -S, no extra files can be given\n");
	      usage(progname, false);
	    }
	  }
	else if (load_file != NULL || load_mmap != NULL)
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
	/*
	if (pair_mode == PAIR_NONE && half_paired) {
	  fprintf(stderr, "error: cannot use option half-paired in non-paired mode!\n");
	  exit(1);
	}
	*/
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

	if ((pair_mode == PAIR_NONE && (match_mode < 1 || match_mode > 2))
	    || (pair_mode != PAIR_NONE && (match_mode < 2 || match_mode > 4))) {
	  fprintf(stderr, "error: invalid match mode [pair_mode=%d;match_mode=%d]\n", pair_mode, match_mode);
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
	  //fprintf(stderr, "warning: window generation threshold is larger than sw threshold\n");
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

	/* Set probabilities from scores */
	if (shrimp_mode == MODE_COLOUR_SPACE) {
	  // CS: pr_xover ~= .03 => alpha => pr_mismatch => rest
	  //pr_xover = .03;
	  score_alpha = (double)crossover_score / (log(pr_xover/3)/log(2.0));
	  pr_mismatch = 1.0/(1.0 + 1.0/3.0 * pow(2.0, ((double)match_score - (double)mismatch_score)/score_alpha));
	} else {
	  // LS: pr_mismatch ~= .01 => alpha => rest
	  pr_mismatch = .01;
	  score_alpha = ((double)match_score - (double)mismatch_score)/(log((1 - pr_mismatch)/(pr_mismatch/3.0))/log(2.0));
	}
	score_beta = (double)match_score - 2 * score_alpha - score_alpha * log(1 - pr_mismatch)/log(2.0);
	pr_del_open = pow(2.0, (double)a_gap_open_score/score_alpha);
	pr_ins_open = pow(2.0, (double)b_gap_open_score/score_alpha);
	pr_del_extend = pow(2.0, (double)a_gap_extend_score/score_alpha);
	pr_ins_extend = pow(2.0, ((double)b_gap_extend_score - score_beta)/score_alpha);

	//score_difference_mq_cutoff = (int)rint(10.0 * score_alpha);

#ifdef DEBUG_SCORES
	fprintf(stderr, "probabilities from scores:\talpha=%.9g\tbeta=%.9g\n", score_alpha, score_beta);
	if (shrimp_mode == MODE_COLOUR_SPACE)
	  fprintf(stderr, "pr_xover=%.9g\n", pr_xover);
	fprintf(stderr, "pr_mismatch=%.9g\npr_del_open=%.9g\tpr_del_extend=%.9g\t(pr_del1=%.9g)\npr_ins_open=%.9g\tpr_ins_extend=%.9g\t(pr_ins1=%.9g)\n",
		pr_mismatch,
		pr_del_open, pr_del_extend, pr_del_open*pr_del_extend, 
		pr_ins_open, pr_ins_extend, pr_ins_open*pr_ins_extend);

	// sanity check:
	fprintf(stderr, "scores from probabilities:\n");
	fprintf(stderr, "match_score=%g\nmismatch_score=%g\n",
		2 * score_alpha + score_beta + score_alpha * log(1 - pr_mismatch) / log(2.0),
		2 * score_alpha + score_beta + score_alpha * log(pr_mismatch/3) / log(2.0));
	if (shrimp_mode == MODE_COLOUR_SPACE)
	  fprintf(stderr, "crossover_score=%g\n", score_alpha * log(pr_xover/3) / log(2.0));
	fprintf(stderr, "a_gap_open_score=%g\ta_gap_extend_score=%g\nb_gap_open_score=%g\tb_gap_extend_score=%g\n",
		score_alpha * log(pr_del_open) / log(2.0),
		score_alpha * log(pr_del_extend) / log(2.0),
		score_alpha * log(pr_ins_open) / log(2.0),
		score_alpha * log(pr_ins_extend) / log(2.0) + score_beta);
#endif

	// set up new options structure
	// THIS SHOULD EVENTUALLY BE MERGED INTO OPTION READING
	if (n_unpaired_mapping_options[0] == 0 && n_paired_mapping_options == 0) {
	  if (pair_mode == PAIR_NONE)
	    {
	      n_unpaired_mapping_options[0]++;
	      //unpaired_mapping_options[0] = (struct read_mapping_options_t *)xcalloc(n_unpaired_mapping_options[0] * sizeof(unpaired_mapping_options[0][0]));
	      unpaired_mapping_options[0] = (struct read_mapping_options_t *)
		my_calloc(n_unpaired_mapping_options[0] * sizeof(unpaired_mapping_options[0][0]),
			  &mem_small, "unpaired_mapping_options[0]");

	      unpaired_mapping_options[0][0].regions.recompute = (match_mode == 2 && use_regions);
	      //unpaired_mapping_options[0][0].regions.min_seed = -1;
	      //unpaired_mapping_options[0][0].regions.max_seed = -1;
	      unpaired_mapping_options[0][0].anchor_list.recompute = true;
	      unpaired_mapping_options[0][0].anchor_list.collapse = true;
	      unpaired_mapping_options[0][0].anchor_list.use_region_counts = (match_mode == 2 && use_regions);
	      unpaired_mapping_options[0][0].anchor_list.use_mp_region_counts = 0;
	      unpaired_mapping_options[0][0].hit_list.recompute = true;
	      unpaired_mapping_options[0][0].hit_list.gapless = gapless_sw;
	      unpaired_mapping_options[0][0].hit_list.match_mode = match_mode;
	      unpaired_mapping_options[0][0].hit_list.threshold = window_gen_threshold;
	      unpaired_mapping_options[0][0].pass1.recompute =  true;
	      unpaired_mapping_options[0][0].pass1.only_paired = false;
	      unpaired_mapping_options[0][0].pass1.gapless = gapless_sw;
	      unpaired_mapping_options[0][0].pass1.min_matches = match_mode; // this is 1 or 2 in unpaired mode
	      unpaired_mapping_options[0][0].pass1.num_outputs = num_tmp_outputs;
	      unpaired_mapping_options[0][0].pass1.threshold = sw_vect_threshold;
	      unpaired_mapping_options[0][0].pass1.window_overlap = window_overlap;
	      unpaired_mapping_options[0][0].pass2.strata = strata_flag;
	      unpaired_mapping_options[0][0].pass2.save_outputs = false;
	      unpaired_mapping_options[0][0].pass2.num_outputs = num_outputs;
	      unpaired_mapping_options[0][0].pass2.threshold = sw_full_threshold;
	      unpaired_mapping_options[0][0].pass2.stop_count = 0;
	    }
	  else
	    {
	      n_paired_mapping_options++;
	      //paired_mapping_options = (struct readpair_mapping_options_t *)xcalloc(n_paired_mapping_options * sizeof(paired_mapping_options[0]));
	      paired_mapping_options = (struct readpair_mapping_options_t *)
		my_calloc(n_paired_mapping_options * sizeof(paired_mapping_options[0]),
			  &mem_small, "paired_mapping_options");

	      paired_mapping_options[0].pairing.pair_mode = pair_mode;
	      paired_mapping_options[0].pairing.min_insert_size = min_insert_size;
	      paired_mapping_options[0].pairing.max_insert_size = max_insert_size;
	      paired_mapping_options[0].pairing.strata = strata_flag;
	      paired_mapping_options[0].pairing.save_outputs = compute_mapping_qualities;
	      paired_mapping_options[0].pairing.pass1_num_outputs = num_tmp_outputs;
	      paired_mapping_options[0].pairing.pass2_num_outputs = num_outputs;
	      paired_mapping_options[0].pairing.pass1_threshold = sw_vect_threshold;
	      paired_mapping_options[0].pairing.pass2_threshold = sw_full_threshold;

	      paired_mapping_options[0].read[0].regions.recompute = use_regions && match_mode != 2;
	      //paired_mapping_options[0].read[0].regions.min_seed = -1;
	      //paired_mapping_options[0].read[0].regions.max_seed = -1;
	      paired_mapping_options[0].read[0].anchor_list.recompute = true;
	      paired_mapping_options[0].read[0].anchor_list.collapse = true;
	      paired_mapping_options[0].read[0].anchor_list.use_region_counts = use_regions && match_mode != 2;
	      if (use_regions) {
		paired_mapping_options[0].read[0].anchor_list.use_mp_region_counts = (match_mode == 4 && !half_paired? 1
										      : match_mode == 3 && half_paired? 2
										      : match_mode == 3 && !half_paired? 3
										      : 0);
	      }
	      paired_mapping_options[0].read[0].hit_list.recompute = true;
	      paired_mapping_options[0].read[0].hit_list.gapless = gapless_sw;
	      paired_mapping_options[0].read[0].hit_list.match_mode = (match_mode == 4? 2
								       : match_mode == 3? 3
								       : 1);
	      paired_mapping_options[0].read[0].hit_list.threshold = window_gen_threshold;
	      paired_mapping_options[0].read[0].pass1.recompute = true;
	      paired_mapping_options[0].read[0].pass1.only_paired = true;
	      paired_mapping_options[0].read[0].pass1.gapless = gapless_sw;
	      paired_mapping_options[0].read[0].pass1.min_matches = (match_mode == 4? 2 : 1);
	      paired_mapping_options[0].read[0].pass1.threshold = sw_vect_threshold;
	      paired_mapping_options[0].read[0].pass1.window_overlap = window_overlap;
	      paired_mapping_options[0].read[0].pass2.strata = strata_flag;
	      paired_mapping_options[0].read[0].pass2.threshold = sw_full_threshold * 0.5;
	      paired_mapping_options[0].read[1] = paired_mapping_options[0].read[0];

	      if (!half_paired)
		{
		  paired_mapping_options[0].pairing.stop_count = 0;
		}
	      else // half_paired
		{
		  paired_mapping_options[0].pairing.stop_count = 1;
		  paired_mapping_options[0].pairing.stop_threshold = 101.0; //paired_mapping_options[0].pairing.pass2_threshold;

		  n_unpaired_mapping_options[0]++;
		  n_unpaired_mapping_options[1]++;
		  //unpaired_mapping_options[0] = (struct read_mapping_options_t *)xcalloc(n_unpaired_mapping_options[0] * sizeof(unpaired_mapping_options[0][0]));
		  unpaired_mapping_options[0] = (struct read_mapping_options_t *)
		    my_calloc(n_unpaired_mapping_options[0] * sizeof(unpaired_mapping_options[0][0]),
			      &mem_small, "unpaired_mapping_options[0]");
		  //unpaired_mapping_options[1] = (struct read_mapping_options_t *)xcalloc(n_unpaired_mapping_options[1] * sizeof(unpaired_mapping_options[1][0]));
		  unpaired_mapping_options[1] = (struct read_mapping_options_t *)
		    my_calloc(n_unpaired_mapping_options[1] * sizeof(unpaired_mapping_options[1][0]),
			      &mem_small, "unpaired_mapping_options[1]");

		  unpaired_mapping_options[0][0].regions.recompute = false;
		  unpaired_mapping_options[0][0].anchor_list.recompute = false;
		  unpaired_mapping_options[0][0].hit_list.recompute = false;
		  unpaired_mapping_options[0][0].pass1.recompute = true;
		  unpaired_mapping_options[0][0].pass1.gapless = gapless_sw;
		  unpaired_mapping_options[0][0].pass1.min_matches = 2;
		  unpaired_mapping_options[0][0].pass1.only_paired = false;
		  unpaired_mapping_options[0][0].pass1.num_outputs = num_tmp_outputs;
		  unpaired_mapping_options[0][0].pass1.threshold = sw_vect_threshold;
		  unpaired_mapping_options[0][0].pass1.window_overlap = window_overlap;
		  unpaired_mapping_options[0][0].pass2.strata = strata_flag;
		  unpaired_mapping_options[0][0].pass2.save_outputs = compute_mapping_qualities;
		  unpaired_mapping_options[0][0].pass2.num_outputs = num_outputs;
		  unpaired_mapping_options[0][0].pass2.threshold = sw_full_threshold;
		  unpaired_mapping_options[0][0].pass2.stop_count = 0;

		  unpaired_mapping_options[1][0] = unpaired_mapping_options[0][0];

		}
	    }
	}

	if(load_file == NULL && load_mmap == NULL) {
	  print_settings();
	}

	if (load_file != NULL && save_mmap != NULL) {
	  exit(genome_load_map_save_mmap(load_file, save_mmap) == true ? 0 : 1);
	}

	before = gettimeinusecs();
	if (load_mmap != NULL) {
	  genome_load_mmap(load_mmap);
	} else if (load_file != NULL){
		if (strchr(load_file, ',') == NULL) {
			//use prefix
			int buf_size = strlen(load_file) + 20;
			//char * genome_name = (char *)xmalloc(sizeof(char)*buf_size);
			char genome_name[buf_size];
			strncpy(genome_name,load_file,buf_size);
			strncat(genome_name,".genome",buf_size);
			fprintf(stderr,"Loading genome from %s\n",genome_name);
			if (!load_genome_map(genome_name)){
				fprintf(stderr, "error: loading from genome file \"%s\"\n", genome_name);
				exit (1);
			}
			//free(genome_name);
			int seed_count = 0;
			//char * seed_name = (char *)xmalloc(sizeof(char)*buf_size);
			char seed_name[buf_size];
			//char * buff = (char *)xmalloc(sizeof(char)*buf_size);
			char buff[buf_size];
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
			//free(seed_name);
			//free(buff);

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

		//print_settings();
	} else {
		if (!load_genome(genome_files,ngenome_files)){
			exit(1);
		}
	}

	load_genome_usecs += (gettimeinusecs() - before);

	// initialize general search tree for contig offsets
	gen_st_init(&contig_offsets_gen_st, 17, contig_offsets, num_contigs);

	//
	// Automatic genome index trimming
	//
	if (Vflag && save_file == NULL && list_cutoff == DEF_LIST_CUTOFF) {
	  // this will be a mapping job; enable automatic trimming
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
	  //fprintf(stderr, "    %-40s%d\n", "Trimming index lists longer than:", list_cutoff);
	  //fprintf(stderr, "\n");
	}

	if (load_file != NULL || load_mmap != NULL) {
	  print_settings();
	}

	if (Yflag)
	  print_genomemap_stats();

	if (save_file != NULL) {
	  if (list_cutoff != DEF_LIST_CUTOFF) {
	    fprintf(stderr, "\nTrimming index lists longer than: %u\n", list_cutoff);
	    trim_genome();
	  }

	  fprintf(stderr,"Saving genome map to %s\n",save_file);
	  if(!save_genome_map(save_file)){
	    exit(1);
	  }
	  exit(0);
	}

	// compute total genome size
	for (cn = 0; cn < num_contigs; cn++)
	  total_genome_size += genome_len[cn];

	//TODO setup need max window and max read len
	//int longest_read_len = 2000;
	int max_window_len = (int)abs_or_pct(window_len,longest_read_len);

	// open input files, and set Qflag accordingly
	if (single_reads_file) {  
	  fasta = fasta_open(reads_filename, shrimp_mode, Qflag, autodetect_input? &Qflag : NULL);
	  if (fasta == NULL) {
	    crash(1, 1, "failed to open read file [%s]", reads_filename);
	  } else {
	    fprintf(stderr, "- Processing read file [%s]\n", reads_filename);
	  }
	} else {
	  left_fasta = fasta_open(left_reads_filename, shrimp_mode, Qflag, autodetect_input? &Qflag : NULL);
	  if (left_fasta == NULL) {
	    crash(1, 1, "failed to open read file [%s]", left_reads_filename);
	  }
	  right_fasta = fasta_open(right_reads_filename, shrimp_mode, Qflag);
	  if (right_fasta == NULL) {
	    crash(1, 1, "failed to open read file [%s]", right_reads_filename);
	  }
	  // WHY? the code above sets both ->space to shrimp_mode anyhow
	  //if (right_fasta->space != left_fasta->space) {
	  //  fprintf(stderr,"error: when using -1 and -2, both files must be either only colour space or only letter space!\n");
	  //  return (false);
	  //}
	  fasta = left_fasta;
	  fprintf(stderr, "- Processing read files [%s , %s]\n", left_reads_filename, right_reads_filename);
	}

	// set default quality value delta
	if (Qflag) {
	  if (!qual_delta_set) {
	    if (shrimp_mode == MODE_LETTER_SPACE)
	      qual_delta = DEF_LS_QUAL_DELTA;
	    else
	      qual_delta = DEF_CS_QUAL_DELTA;
	    logit(0, "quality value format not set explicitly; using PHRED+%d", qual_delta);
	  } else {
	    logit(0, "quality value format set to PHRED+%d", qual_delta);
	  }
	}


#pragma omp parallel shared(longest_read_len,max_window_len,a_gap_open_score, a_gap_extend_score, b_gap_open_score, b_gap_extend_score,\
		match_score, mismatch_score,shrimp_mode,crossover_score,anchor_width) num_threads(num_threads)
	{
	  // init thread-private globals
	  memset(&tpg, 0, sizeof(tpg_t));
	  tpg.wait_tc.type = DEF_FAST_TIME_COUNTER;
	  tpg.region_counts_tc.type = DEF_FAST_TIME_COUNTER;
	  tpg.mp_region_counts_tc.type = DEF_FAST_TIME_COUNTER;
	  tpg.anchor_list_tc.type = DEF_FAST_TIME_COUNTER;
	  tpg.hit_list_tc.type = DEF_FAST_TIME_COUNTER;
	  tpg.get_vector_hits_tc.type = DEF_FAST_TIME_COUNTER;
	  tpg.pass1_tc.type = DEF_FAST_TIME_COUNTER;
	  tpg.pass2_tc.type = DEF_FAST_TIME_COUNTER;
	  tpg.duplicate_removal_tc.type = DEF_FAST_TIME_COUNTER;

	  /* region handling */
	  if (use_regions) {
	    region_map_id = 0;
	    for (int number_in_pair = 0; number_in_pair < 2; number_in_pair++)
	      for (int st = 0; st < 2; st++)
		//region_map[number_in_pair][st] = (int32_t *)xcalloc(n_regions * sizeof(region_map[0][0][0]));
		region_map[number_in_pair][st] = (region_map_t *)
		  my_calloc(n_regions * sizeof(region_map[0][0][0]),
			    &mem_mapping, "region_map");
	  }

	  if (f1_setup(max_window_len, longest_read_len,
		       a_gap_open_score, a_gap_extend_score, b_gap_open_score, b_gap_extend_score,
		       match_score, shrimp_mode == MODE_LETTER_SPACE? mismatch_score : match_score + crossover_score,
		       shrimp_mode == MODE_COLOUR_SPACE, true)) {
	    fprintf(stderr, "failed to initialise vector Smith-Waterman (%s)\n", strerror(errno));
	    exit(1);
	  }

	  int ret;
	  if (shrimp_mode == MODE_COLOUR_SPACE) {
	    /* XXX - a vs. b gap */
	    ret = sw_full_cs_setup(max_window_len, longest_read_len,
				   a_gap_open_score, a_gap_extend_score, b_gap_open_score, b_gap_extend_score,
				   match_score, mismatch_score,
				   crossover_score, true, anchor_width, indel_taboo_len);
	  } else {
	    ret = sw_full_ls_setup(max_window_len, longest_read_len,
				   a_gap_open_score, a_gap_extend_score, b_gap_open_score, b_gap_extend_score,
				   match_score, mismatch_score, true, anchor_width);
	  }
	  if (ret) {
	    fprintf(stderr, "failed to initialise scalar Smith-Waterman (%s)\n", strerror(errno));
	    exit(1);
	  }

	  /* post_sw */
	  if (shrimp_mode == MODE_COLOUR_SPACE) {
	    post_sw_setup(max_window_len + longest_read_len,
			  pr_mismatch, pr_xover, pr_del_open, pr_del_extend, pr_ins_open, pr_ins_extend,
			  Qflag && !ignore_qvs, use_sanger_qvs, qual_vector_offset, qual_delta, true);
	  }

	}


	char * output;
	if (Eflag){
	  int i;
	  if (sam_header_filename != NULL) {
	    FILE * sam_header_file = fopen(sam_header_filename, "r");
	    if (sam_header_file == NULL) {
	      perror("Failed to open sam header file ");
	      exit(1);
	    }
	    cat(sam_header_file, stdout);
	    fclose(sam_header_file);
	  } else {
	    // HD line
	    if (sam_header_hd != NULL) {
	      cat(sam_header_hd, stdout);
	    } else {
	      fprintf(stdout,"@HD\tVN:%s\tSO:%s\n","1.0","unsorted");
	    }

	    // SQ lines
	    if (sam_header_sq != NULL) {
	      cat(sam_header_sq, stdout);
	    } else {
	      for(i = 0; i < num_contigs; i++){
		fprintf(stdout,"@SQ\tSN:%s\tLN:%u\n",contig_names[i],genome_len[i]);
	      }
	    }

	    // RG lines
	    if (sam_header_rg != NULL) {
	      cat(sam_header_rg, stdout);
	    } else if (sam_read_group_name != NULL) {
	      fprintf(stdout, "@RG\tID:%s\tSM:%s\n", sam_read_group_name, sam_sample_name);
	    }

	    // PG lines
	    if (sam_header_pg != NULL) {
	      cat(sam_header_pg, stdout);
	    } else {
	      fprintf(stdout, "@PG\tID:%s\tVN:%s\tCL:%s\n", "gmapper", SHRIMP_VERSION_STRING, command_line);
	    }
	  }
	} else {
	  output = output_format_line(Rflag);
	  puts(output);
	  free(output);
	}
	before = gettimeinusecs();
	bool launched = launch_scan_threads(fasta, left_fasta, right_fasta);
	if (!launched) {
	  fprintf(stderr,"error: a fatal error occured while launching scan thread(s)!\n");
	  exit(1);
	}
	mapping_wallclock_usecs += (gettimeinusecs() - before);

	if (single_reads_file) {
	  fasta_close(fasta);
	} else {
	  fasta_close(left_fasta);
	  fasta_close(right_fasta);
	}
	
	print_statistics();
#pragma omp parallel shared(longest_read_len,max_window_len,a_gap_open_score, a_gap_extend_score, b_gap_open_score, b_gap_extend_score,	\
			    match_score, mismatch_score,shrimp_mode,crossover_score,anchor_width) num_threads(num_threads)
	{
	  sw_vector_cleanup();
	  if (shrimp_mode==MODE_COLOUR_SPACE) {
	    sw_full_cs_cleanup();
	    post_sw_cleanup();
	  }
	  sw_full_ls_cleanup();
	  f1_free();

	  if (use_regions) {
	    for (int number_in_pair = 0; number_in_pair < 2; number_in_pair++)
	      for (int st = 0; st < 2; st++)
		//free(region_map[number_in_pair][st]);
		my_free(region_map[number_in_pair][st], n_regions * sizeof(region_map[0][0][0]),
			&mem_mapping, "region_map");
	  }
	}

	gen_st_delete(&contig_offsets_gen_st);

	if (load_mmap != NULL) {
	  // munmap?
	} else {
	  free_genome();

	  //free(seed);
	  if (Hflag) {
	    int sn;
	    for (sn = 0; sn < n_seeds; sn++) {
	      my_free(seed_hash_mask[sn], BPTO32BW(max_seed_span) * sizeof(seed_hash_mask[sn][0]),
		      &mem_small, "seed_hash_mask[%d]", sn);
	    }
	    my_free(seed_hash_mask, n_seeds * sizeof(seed_hash_mask[0]),
		    &mem_small, "seed_hash_mask");
	  }
	  my_free(seed, n_seeds * sizeof(seed[0]),
		  &mem_small, "seed");
	}

	//free(paired_mapping_options);
	my_free(paired_mapping_options, n_paired_mapping_options * sizeof(paired_mapping_options[0]),
		&mem_small, "paired_mapping_options");
	//free(unpaired_mapping_options[0]);
	my_free(unpaired_mapping_options[0], n_unpaired_mapping_options[0] * sizeof(unpaired_mapping_options[0][0]),
		&mem_small, "unpaired_mapping_options[0]");
	//free(unpaired_mapping_options[1]);
	my_free(unpaired_mapping_options[1], n_unpaired_mapping_options[1] * sizeof(unpaired_mapping_options[0][0]),
		&mem_small, "unpaired_mapping_options[1]");

	if (sam_read_group_name != NULL)
	  free(sam_read_group_name);

	// close some files
	if (aligned_reads_file != NULL)
	  fclose(aligned_reads_file);
	if (unaligned_reads_file != NULL)
	  fclose(unaligned_reads_file);
	if (sam_header_hd != NULL)
	  fclose(sam_header_hd);
	if (sam_header_sq != NULL)
	  fclose(sam_header_sq);
	if (sam_header_rg != NULL)
	  fclose(sam_header_rg);
	if (sam_header_pg != NULL)
	  fclose(sam_header_pg);

#ifdef MYALLOC_ENABLE_CRT
	fprintf(stderr, "crt_mem: %lld\n", (long long)crt_mem);
#endif
#ifndef NDEBUG
	fprintf(stderr, "mem_genomemap: max=%lld crt=%lld\n", (long long)count_get_max(&mem_genomemap), (long long)count_get_count(&mem_genomemap));
	fprintf(stderr, "mem_mapping: max=%lld crt=%lld\n", (long long)count_get_max(&mem_mapping), (long long)count_get_count(&mem_mapping));
	fprintf(stderr, "mem_thread_buffer: max=%lld crt=%lld\n", (long long)count_get_max(&mem_thread_buffer), (long long)count_get_count(&mem_thread_buffer));
	fprintf(stderr, "mem_small: max=%lld crt=%lld\n", (long long)count_get_max(&mem_small), (long long)count_get_count(&mem_small));
	fprintf(stderr, "mem_sw: max=%lld crt=%lld\n", (long long)count_get_max(&mem_sw), (long long)count_get_count(&mem_sw));
#endif
	return 0;
}
