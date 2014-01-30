/*
 * This file should contain only #defines with default values.
 */

#ifndef _GMAPPER_DEFAULTS_H
#define _GMAPPER_DEFAULTS_H


#define DEF_SHRIMP_MODE		MODE_LETTER_SPACE

#define DEF_NUM_THREADS		1
#define DEF_MAX_THREADS		100
#define DEF_CHUNK_SIZE		1000
#define DEF_PROGRESS		100000
#define USE_PREFETCH

#define DEF_HASH_FILTER_CALLS	true
#define DEF_GAPLESS_SW		false
#define DEF_LIST_CUTOFF		4294967295u // 2^32 - 1

#define DEF_USE_REGIONS		true
#define DEF_REGION_BITS		11
#define DEF_REGION_OVERLAP	50

#define DEF_ANCHOR_LIST_BIG_GAP	1024

#define DEF_PAIR_MODE		PAIR_NONE
#define DEF_MIN_INSERT_SIZE	0
#define DEF_MAX_INSERT_SIZE	1000
#define DEF_INSERT_SIZE_MEAN	200
#define DEF_INSERT_SIZE_STDDEV	100

#define DEF_WINDOW_LEN		140.0
#define DEF_WINDOW_OVERLAP	90.0
#define DEF_MATCH_MODE_PAIRED	4
#define DEF_MATCH_MODE_UNPAIRED	2
#define DEF_NUM_OUTPUTS		10
#define DEF_ANCHOR_WIDTH	8	/* width around anchors in full SW */
#define DEF_INDEL_TABOO_LEN	0

#define DEF_LS_QUAL_DELTA	64
#define DEF_CS_QUAL_DELTA	33

/* SW Scores */
#define DEF_LS_MATCH_SCORE	10
#define DEF_LS_MISMATCH_SCORE	-15
#define DEF_LS_A_GAP_OPEN	-33
#define DEF_LS_B_GAP_OPEN	DEF_LS_A_GAP_OPEN
#define DEF_LS_A_GAP_EXTEND	-7
#define DEF_LS_B_GAP_EXTEND	-3

#define DEF_CS_MATCH_SCORE	10
#define DEF_CS_MISMATCH_SCORE	-24
#define DEF_CS_XOVER_SCORE	-20
#define DEF_CS_A_GAP_OPEN	-33
#define DEF_CS_B_GAP_OPEN	DEF_CS_A_GAP_OPEN
#define DEF_CS_A_GAP_EXTEND	-7
#define DEF_CS_B_GAP_EXTEND	-3


/* Score Thresholds */
#define DEF_WINDOW_GEN_THRESHOLD	(55.0)	/* Min required to generate match window */
//SHRiMP v 2.0.1
//#define DEF_SW_VECT_THRESHOLD	60.0	/* == DEF_SW_FULL_THRESHOLD in lspace */
//#define DEF_SW_FULL_THRESHOLD	68.0	/* read_length x match_value x .68 */
//SHRiMP v 2.0.2
#define DEF_SW_VECT_THRESHOLD	(47.0)	/* == DEF_SW_FULL_THRESHOLD in lspace */
#define DEF_SW_FULL_THRESHOLD	(50.0)	/* read_length x match_value x .55 */


#define DEF_MAX_ALIGNMENTS 0
#define DEF_LONGEST_READ_LENGTH	1000

#define DEF_STANDARD_OPTIONS \
{\
	{"un",1,0,10},\
	{"al",1,0,11},\
	{"upstream",1,0,'1'},\
	{"downstream",1,0,'2'},\
	{"sam-unaligned",0,0,12},\
	{"longest-read",1,0,13},\
	{"seeds",1,0,'s'},\
	{"report",1,0,'o'},\
	{"match-window",1,0,'w'},\
	{"cmw-mode",1,0,'n'},\
	{"cmw-overlap",1,0,'l'},\
	{"anchor-width",1,0,'a'},\
	{"save",1,0,'S'},\
	{"load",1,0,'L'},\
	{"cutoff",1,0,'z'},\
	{"match",1,0,'m'},\
	{"mismatch",1,0,'i'},\
	{"open-r",1,0,'g'},\
	{"open-q",1,0,'q'},\
	{"ext-r",1,0,'e'},\
	{"ext-q",1,0,'f'},\
	{"cmw-threshold",1,0,'r'},\
	{"full-threshold",1,0,'h'},\
	{"threads",1,0,'N'},\
	{"thread-chunk",1,0,'K'},\
	{"pair-mode",1,0,'p'},\
	{"isize",1,0,'I'},\
	{"ungapped",0,0,'U'},\
	{"negative",0,0,'C'},\
	{"positive",0,0,'F'},\
	{"pretty",0,0,'P'},\
	{"sam",0,0,'E'},\
	{"fastq",0,0,'Q'},\
	{"print-reads",0,0,'R'},\
	{"rev-tiebreak",0,0,'T'},\
	{"tiebreak-off",0,0,'t'},\
	{"isize-histogram",0,0,'X'},\
	{"proj-histogram",0,0,'Y'},\
	{"cachebypass-off",0,0,'Z'},\
	{"help",0,0,'?'},\
	{"hash-spaced-kmers",0,0,'H'},\
	{"thread-stats",0,0,'D'},\
	{"trim-off",0,0,'V'},\
	{"strata",0,0,9},\
	{"max-alignments",1,0,14},\
	{"global",0,0,15},\
	{"read-group",1,0,17},\
	{"sam-header",1,0,18},\
	{"no-half-paired",0,0,19},\
	{"sam-r2",0,0,20},\
	{"mode",1,0,'M'},\
	{"trim-front",1,0,21},\
	{"trim-end",1,0,22},\
	{"trim-first",0,0,23},\
	{"trim-second",0,0,24},\
	{"insert-size-dist",1,0,25},\
	{"use-regions",0,0,26},\
	{"region-overlap",1,0,27},\
	{"paired-options",1,0,28},\
	{"unpaired-options",1,0,29},\
	{"min-avg-qv",1,0,30},\
	{"extra-sam-fields",0,0,31},\
	{"region-bits",1,0,32},\
	{"progress",1,0,33},\
	{"save-mmap",1,0,34},\
	{"load-mmap",1,0,35},\
	{"indel-taboo-len",1,0,36},\
	{"single-best-mapping",0,0,37},\
	{"all-contigs",0,0,38},\
	{"no-mapping-qualities",0,0,39},\
	{"shrimp-format",0,0,40},\
	{"half-paired",0,0,41},\
	{"no-improper-mappings",0,0,42},\
	{"qv-offset",1,0,43},\
	{"sam-header-hd",1,0,44},\
	{"sam-header-sq",1,0,45},\
	{"sam-header-rg",1,0,46},\
	{"sam-header-pg",1,0,47},\
	{"no-autodetect-input",0,0,48},\
	{"local",0,0,124},\
	{"no-qv-check",0,0,123},\
	{"ignore-qvs",0,0,125},\
	{"pr-xover",1,0,126}\
}

#define DEF_COLOUR_SPACE_OPTIONS \
{\
	{"crossover",1,0,'x'},\
	{"vec-threshold",1,0,'v'},\
	{"bfast",0,0,16},\
	{0,0,0,0}\
}

#define DEF_LETTER_SPACE_OPTIONS \
{\
	{"trim-illumina",0,0,3},\
	{0,0,0,0}\
}

#define DEF_PAIR_MODE_STRING \
{\
  "none",\
  "opposing strands; inwards",\
  "opposing strands; outwards",\
  "same strand; second is forward",\
  "same strand; second is backward"\
}

#define DEF_PAIR_REVERSE \
{\
  { 0, 0 },\
  { 0, 0 },\
  { 1, 1 },\
  { 0, 1 },\
  { 1, 0 }\
}


#define DEF_DEF_SEEDS_CS_MIN_WEIGHT 10
#define DEF_DEF_SEEDS_CS_MAX_WEIGHT 18
#define DEF_DEF_SEEDS_CS_WEIGHT 12
#define DEF_DEF_SEEDS_CS_CNT { 4, 4, 3, 0, 0, 0, 4, 0, 4 }
#define DEF_DEF_SEEDS_CS \
{\
  { "111110011111", "111100110001111", "111100100100100111", "111001000100001001111" }, \
  { "1111001111111", "1111100110001111", "11110010010001001111", "11100110010000100100111" },\
  { "11110111101111", "1111011100100001111", "1111000011001101111"}, \
  { },\
  { },\
  { },\
  { "111111101110111111", "1111100101101101011111", "11110011001010100011011111", "111101001100000100110011010111" },\
  { },\
  { "11111011111110111111", "11110111011010111011111", "11111100110101101001011111", "11111010101100100010011101111" }\
}
// { "111110001111111", "111100111001001111", "111001001000111001111", "1111001000010001001001111" },

#define DEF_DEF_SEEDS_LS_MIN_WEIGHT 10
#define DEF_DEF_SEEDS_LS_MAX_WEIGHT 18
#define DEF_DEF_SEEDS_LS_WEIGHT 12
#define DEF_DEF_SEEDS_LS_CNT { 4, 4, 3, 0, 0, 0, 4, 0, 4 }
#define DEF_DEF_SEEDS_LS \
{\
  { "111110011111", "111100110001111", "111100100100100111", "111001000100001001111" },\
  { "1111001111111", "1111100110001111", "11110010010001001111", "11100110010000100100111" },\
  { "11110111101111", "1111011100100001111", "1111000011001101111"}, \
  { },\
  { },\
  { },\
  { "111111101110111111", "1111100101101101011111", "11110011001010100011011111", "111101001100000100110011010111" },\
  { },\
  { "11111011111110111111", "11110111011010111011111", "11111100110101101001011111", "11111010101100100010011101111" }\
}
//  { "111101101011111", "111010110011001111", "1110110001011010111", "11110010100000100110111" },

#define DEF_DEF_SEEDS_MIRNA_CNT 5
#define DEF_DEF_SEEDS_MIRNA \
{\
  "00111111001111111100",\
  "00111111110011111100",\
  "00111111111100111100",\
  "00111111111111001100",\
  "00111111111111110000"\
}


#define DEF_THREAD_OUTPUT_BUFFER_INITIAL	1024*1024*10
#define DEF_THREAD_OUTPUT_BUFFER_INCREMENT	1024*1024*10
#define DEF_THREAD_OUTPUT_BUFFER_SAFETY		1024*500
#define DEF_THREAD_OUTPUT_HEAP_CAPACITY		16384


#endif
