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
DEF_HEAP(uint32_t, char *, out)



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
read_free_full(struct read_entry * re)
{
  read_free(re);
  free(re->read[0]);
  free(re->read[1]);

  free(re->mapidx[0]);
  free(re->mapidx[1]);
  read_free_anchor_list(re);
  read_free_hit_list(re);

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
    re_buffer = (struct read_entry *)xmalloc_m(chunk_size * sizeof(re_buffer[0]), "re_buffer");

    while (read_more) {
      memset(re_buffer, 0, chunk_size * sizeof(re_buffer[0]));

      before = rdtsc();

      //Read in this threads 'chunk'
#pragma omp critical (fill_reads_buffer)
      {
	wait_ticks[omp_get_thread_num()] += rdtsc() - before;

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

	re_buffer[i].read[0] = fasta_sequence_to_bitfield(fasta, re_buffer[i].seq);
	re_buffer[i].read_len = strlen(re_buffer[i].seq);
	re_buffer[i].max_n_kmers = re_buffer[i].read_len - min_seed_span + 1;
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
	    re_buffer[i].ignore = true;
	  }
	  continue;	
	}

	re_buffer[i].window_len = (uint16_t)abs_or_pct(window_len,re_buffer[i].read_len);

	if (re_buffer[i].range_string != NULL) {
	  read_compute_ranges(&re_buffer[i]);
	  free(re_buffer[i].range_string);
	  re_buffer[i].range_string = NULL;
	}
	//free(re_buffer[i].seq);

	// time to do some mapping!
	if (pair_mode == PAIR_NONE)
	  {
	    new_handle_read(&re_buffer[i], unpaired_mapping_options[0], n_unpaired_mapping_options[0]);
	    read_free_full(&re_buffer[i]);
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
	    new_handle_readpair(&re_buffer[i-1], &re_buffer[i], paired_mapping_options, n_paired_mapping_options);
	    read_free_full(&re_buffer[i-1]);
	    read_free_full(&re_buffer[i]);	    
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
	double anchor_list_secs[num_threads], hit_list_secs[num_threads];
	double region_counts_secs[num_threads], duplicate_removal_secs[num_threads];
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
	  anchor_list_secs[tid] = (double)anchor_list_ticks[tid] / hz;
          hit_list_secs[tid] = (double)hit_list_ticks[tid] / hz;
          duplicate_removal_secs[tid] = (double)duplicate_removal_ticks[tid] / hz;
          region_counts_secs[tid] = (double)region_counts_ticks[tid] / hz;
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
	  fprintf(stderr, "%s%s" "%11s %9s %9s %9s %9s %9s %9s %25s %25s %9s\n", my_tab, my_tab,
		  "", "Read Load", "Scan", "Reg Cnts", "Anch List", "Hit List", "Dup Remv",
		  "Vector SW", "Scalar SW", "Wait");
	  fprintf(stderr, "%s%s" "%11s %9s %9s %9s %9s %9s %9s %15s %9s %15s %9s %9s\n", my_tab, my_tab,
		  "", "Time", "Time", "Time", "Time", "Time", "Time",
		  "Invocs", "Time", "Invocs", "Time", "Time");
	  fprintf(stderr, "\n");
	  for(i = 0; i < num_threads; i++) {
	    fprintf(stderr, "%s%s" "Thread %-4d %9.2f %9.2f %9.2f %9.2f %9.2f %9.2f %15s %9.2f %15s %9.2f %9.2f\n", my_tab, my_tab,
		    i, readload_secs[i], scan_secs[i],
		    region_counts_secs[i], anchor_list_secs[i], hit_list_secs[i], duplicate_removal_secs[i],
		    comma_integer(f1_invocs[i]), f1_secs[i],
		    comma_integer(f2_invocs[i]), f2_secs[i], (double)wait_ticks[i] / hz);
	  }
          for (i = 0; i < num_threads; i++) {
            fprintf (stderr, "thrd:%d anchor_list_init_size:(%.2f, %.2f) anchors_discarded:(%.2f, %.2f) big_gaps:(%.2f, %.2f)\n",
              i, stat_get_mean(&anchor_list_init_size[i]), stat_get_sample_stddev(&anchor_list_init_size[i]),
              stat_get_mean(&n_anchors_discarded[i]), stat_get_sample_stddev(&n_anchors_discarded[i]),
              stat_get_mean(&n_big_gaps_anchor_list[i]), stat_get_sample_stddev(&n_big_gaps_anchor_list[i]));
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

	    if (sam_half_paired) {
	      fprintf(stderr, "\n");

	      fprintf(stderr, "%s%s%-40s" "%s    (%.4f%%)\n", my_tab, my_tab,
		      "Additional Reads Matched Unpaired:",
		      comma_integer(total_reads_matched),
		      (nreads == 0) ? 0 : ((double)total_reads_matched / (double)nreads) * 100);
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
  fprintf(stderr, "%s%-40s%s\n", my_tab, "Region Filter:", use_regions? "yes" : "no");
  if (use_regions) {
  fprintf(stderr, "%s%-40s%d\n", my_tab, "Region Overlap:", region_overlap);
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
			num_tmp_outputs = 20 + num_outputs;
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
	    unpaired_mapping_options[0] = (struct read_mapping_options_t *)malloc(n_unpaired_mapping_options[0] * sizeof(unpaired_mapping_options[0][0]));

	    unpaired_mapping_options[0][0].regions.recompute = use_regions;
	    unpaired_mapping_options[0][0].regions.min_seed = 0;
	    unpaired_mapping_options[0][0].regions.max_seed = n_seeds - 1;
	    unpaired_mapping_options[0][0].anchor_list.recompute = true;
	    unpaired_mapping_options[0][0].anchor_list.collapse = true;
	    unpaired_mapping_options[0][0].anchor_list.use_region_counts = use_regions;
	    unpaired_mapping_options[0][0].anchor_list.min_count[0] = (num_matches == 2? 2 : 1);
	    unpaired_mapping_options[0][0].anchor_list.max_count[0] = 0;
	    unpaired_mapping_options[0][0].anchor_list.min_count[1] = 0;
	    unpaired_mapping_options[0][0].anchor_list.max_count[1] = 0;
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
	    unpaired_mapping_options[0][0].pass2.strata = strata_flag;
	    unpaired_mapping_options[0][0].pass2.num_outputs = num_outputs;
	    unpaired_mapping_options[0][0].pass2.threshold = sw_full_threshold;
	    unpaired_mapping_options[0][0].pass2.stop_count = 0;
	  }
	else
	  {
	    n_paired_mapping_options++;
	    paired_mapping_options = (struct readpair_mapping_options_t *)malloc(n_paired_mapping_options * sizeof(paired_mapping_options[0]));

	    paired_mapping_options[0].pairing.pair_mode = pair_mode;
	    paired_mapping_options[0].pairing.pair_up_hits = true;
	    paired_mapping_options[0].pairing.min_insert_size = min_insert_size;
	    paired_mapping_options[0].pairing.max_insert_size = max_insert_size;
	    paired_mapping_options[0].pairing.strata = strata_flag;
	    paired_mapping_options[0].pairing.min_num_matches = 3;
	    paired_mapping_options[0].pairing.pass1_num_outputs = num_tmp_outputs;
	    paired_mapping_options[0].pairing.pass2_num_outputs = num_outputs;
	    paired_mapping_options[0].pairing.pass1_threshold = sw_vect_threshold;
	    paired_mapping_options[0].pairing.pass2_threshold = sw_full_threshold;

	    paired_mapping_options[0].read[0].regions.recompute = use_regions;
	    paired_mapping_options[0].read[0].regions.min_seed = 0;
	    paired_mapping_options[0].read[0].regions.max_seed = n_seeds - 1;

	    paired_mapping_options[0].read[0].anchor_list.recompute = true;
	    paired_mapping_options[0].read[0].anchor_list.collapse = true;
	    paired_mapping_options[0].read[0].anchor_list.use_region_counts = use_regions;
	    paired_mapping_options[0].read[0].anchor_list.min_count[0] = (num_matches == 4? 2 : 1);
	    paired_mapping_options[0].read[0].anchor_list.min_count[1] = 0; //(num_matches == 4? 2 : 1);
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
		paired_mapping_options[0].pairing.stop_count = 1;
		paired_mapping_options[0].pairing.stop_threshold = paired_mapping_options[0].pairing.pass2_threshold;

		n_unpaired_mapping_options[0]++;
		n_unpaired_mapping_options[1]++;
		unpaired_mapping_options[0] = (struct read_mapping_options_t *)malloc(n_unpaired_mapping_options[0] * sizeof(unpaired_mapping_options[0][0]));
		unpaired_mapping_options[1] = (struct read_mapping_options_t *)malloc(n_unpaired_mapping_options[1] * sizeof(unpaired_mapping_options[1][0]));

		unpaired_mapping_options[0][0].regions.recompute = false;
		//if (!use_regions) {
		  unpaired_mapping_options[0][0].anchor_list.recompute = false;
		  unpaired_mapping_options[0][0].hit_list.recompute = false;
		/*
		} else {
		  unpaired_mapping_options[0][0].anchor_list.recompute = true;
		  unpaired_mapping_options[0][0].anchor_list.collapse = true;
		  unpaired_mapping_options[0][0].anchor_list.use_region_counts = use_regions;
		  unpaired_mapping_options[0][0].anchor_list.min_count[0] = (num_matches == 2? 2 : 1);
		  unpaired_mapping_options[0][0].anchor_list.max_count[0] = 0;
		  unpaired_mapping_options[0][0].anchor_list.min_count[1] = 0;
		  unpaired_mapping_options[0][0].anchor_list.max_count[1] = 0;
		  unpaired_mapping_options[0][0].hit_list.recompute = true;
		  unpaired_mapping_options[0][0].hit_list.gapless = gapless_sw;
		  unpaired_mapping_options[0][0].hit_list.match_mode = num_matches;
		  unpaired_mapping_options[0][0].hit_list.threshold = window_gen_threshold;
		}
		*/
		unpaired_mapping_options[0][0].pass1.recompute = true;
		unpaired_mapping_options[0][0].pass1.gapless = gapless_sw;
		unpaired_mapping_options[0][0].pass1.only_paired = false;
		unpaired_mapping_options[0][0].pass1.num_outputs = num_tmp_outputs;
		unpaired_mapping_options[0][0].pass1.threshold = sw_vect_threshold;
		unpaired_mapping_options[0][0].pass1.window_overlap = window_overlap;
		unpaired_mapping_options[0][0].pass2.strata = strata_flag;
		unpaired_mapping_options[0][0].pass2.num_outputs = num_outputs;
		unpaired_mapping_options[0][0].pass2.threshold = sw_full_threshold;
		unpaired_mapping_options[0][0].pass2.stop_count = 0;
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

	  /* region handling */
	  if (use_regions) {
	    region_map_id = 0;
	    for (int number_in_pair = 0; number_in_pair < 2; number_in_pair++)
	      for (int st = 0; st < 2; st++)
		region_map[number_in_pair][st] = (int32_t *)xcalloc(n_regions * sizeof(region_map[0][0][0]));
	  }

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
#pragma omp parallel shared(longest_read_len,max_window_len,a_gap_open_score, a_gap_extend_score, b_gap_open_score, b_gap_extend_score,	\
			    match_score, mismatch_score,shrimp_mode,crossover_score,anchor_width) num_threads(num_threads)
	{
	  sw_vector_cleanup();
	  if (shrimp_mode==MODE_COLOUR_SPACE) {
	    sw_full_cs_cleanup();
	  }
	  sw_full_ls_cleanup();
	  f1_free();

	  if (use_regions) {
	    for (int number_in_pair = 0; number_in_pair < 2; number_in_pair++)
	      for (int st = 0; st < 2; st++)
		free(region_map[number_in_pair][st]);
	  }
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
