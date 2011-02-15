/*
 * This file should contain only #defines with default values.
 */

#ifndef _GMAPPER_DEFAULTS_H
#define _GMAPPER_DEFAULTS_H


#define DEF_SHRIMP_MODE		MODE_LETTER_SPACE

#define DEF_NUM_THREADS		1
#define DEF_CHUNK_SIZE		1000

#define DEF_HASH_FILTER_CALLS	true
#define DEF_GAPLESS_SW		false
#define DEF_LIST_CUTOFF		4294967295u // 2^32 - 1

#define DEF_PAIR_MODE		PAIR_NONE
#define DEF_MIN_INSERT_SIZE	50
#define DEF_MAX_INSERT_SIZE	2000

#define DEF_WINDOW_LEN		140.0
#define DEF_WINDOW_OVERLAP	20.0
#define DEF_NUM_MATCHES		2
#define DEF_NUM_OUTPUTS		10
#define DEF_ANCHOR_WIDTH	8	/* width around anchors in full SW */

/* SW Scores */
#define DEF_MATCH_VALUE		10
#define DEF_MISMATCH_VALUE	-15
#define DEF_A_GAP_OPEN		-40
#define DEF_B_GAP_OPEN		 DEF_A_GAP_OPEN
#define DEF_A_GAP_EXTEND	-7
#define DEF_B_GAP_EXTEND	 DEF_A_GAP_EXTEND
#define DEF_XOVER_PENALTY	-14	/* CS only */

/* Score Thresholds */
#define DEF_WINDOW_GEN_THRESHOLD	55.0	/* Min required to generate match window */
//SHRiMP v 2.0.1
//#define DEF_SW_VECT_THRESHOLD	60.0	/* == DEF_SW_FULL_THRESHOLD in lspace */
//#define DEF_SW_FULL_THRESHOLD	68.0	/* read_length x match_value x .68 */
//SHRiMP v 2.0.2
#define DEF_SW_VECT_THRESHOLD	50.0	/* == DEF_SW_FULL_THRESHOLD in lspace */
#define DEF_SW_FULL_THRESHOLD	55.0	/* read_length x match_value x .55 */


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
	{"cmv-threshold",1,0,'r'},\
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
	{"spaced-kmers",0,0,'H'},\
	{"thread-stats",0,0,'D'},\
	{"trim-off",0,0,'V'},\
	{"strata",0,0,9},\
	{"max-alignments",1,0,14},\
	{"global",0,0,15},\
	{"read-group",1,0,17},\
	{"sam-header",1,0,18},\
	{"half-paired",0,0,19},\
	{"sam-r2",0,0,20},\
	{"mode",1,0,'M'},\
	{"trim-front",1,0,200},\
	{"trim-end",1,0,201},\
	{"trim-first",1,0,202},\
	{"trim-second",1,0,203},\
	{"expected-isize",1,0,204}\
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
#define DEF_DEF_SEEDS_CS_CNT { 4, 4, 4, 0, 0, 0, 4, 0, 4 }
#define DEF_DEF_SEEDS_CS \
{\
  { "111110011111", "111100110001111", "111100100100100111", "111001000100001001111" }, \
  { "1111001111111", "1111100110001111", "11110010010001001111", "11100110010000100100111" },\
  { "111110001111111", "111100111001001111", "111001001000111001111", "1111001000010001001001111" },\
  { },\
  { },\
  { },\
  { "111111101110111111", "1111100101101101011111", "11110011001010100011011111", "111101001100000100110011010111" },\
  { },\
  { "11111011111110111111", "11110111011010111011111", "11111100110101101001011111", "11111010101100100010011101111" }\
}

#define DEF_DEF_SEEDS_LS_MIN_WEIGHT 10
#define DEF_DEF_SEEDS_LS_MAX_WEIGHT 18
#define DEF_DEF_SEEDS_LS_WEIGHT 12
#define DEF_DEF_SEEDS_LS_CNT { 4, 4, 4, 0, 0, 0, 4, 0, 4 }
#define DEF_DEF_SEEDS_LS \
{\
  { "111110011111", "111100110001111", "111100100100100111", "111001000100001001111" },\
  { "1111001111111", "1111100110001111", "11110010010001001111", "11100110010000100100111" },\
  { "111101101011111", "111010110011001111", "1110110001011010111", "11110010100000100110111" },\
  { },\
  { },\
  { },\
  { "111111101110111111", "1111100101101101011111", "11110011001010100011011111", "111101001100000100110011010111" },\
  { },\
  { "11111011111110111111", "11110111011010111011111", "11111100110101101001011111", "11111010101100100010011101111" }\
}

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
#define DEF_THREAD_OUTPUT_HEAP_CAPACITY		1024


#endif
