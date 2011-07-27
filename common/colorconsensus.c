#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <ctype.h>
#include <libgen.h>
#include <getopt.h>

#define left(i) ( ((i) >> 2) & 3)
#define right(i) ( (i) & 3)
#define MIN2(i,j) ( ((i)< (j))? (i):(j))

/******************************************************************************
 * CONSTANTS
 ******************************************************************************/
static const char * CMDLINE_FLAGS = "e:c:?";
static const struct option CMDLINE_OPTIONS[] = {
            { "error-rate", required_argument, NULL, 'e' },
            { "conf-level", required_argument, NULL, 'c' },
            { "help", no_argument, NULL, '?' },
            { 0, 0, 0, 0}
};

#define DEFAULT_ERROR_RATE          0.04
#define DEFAULT_CONF_LEVEL          0.9


double glob_errrate;
double confrate;

double neglogsixteenth;
double neglogfourth;

int vitorfb = 0;

int baseoffset = 1; //this is 0 if 0-based coords, 1 if 1-based. not that 0 
                    // is possibly buggy due to -0 being misinterpreted.


typedef struct column{
  double forwards[16]; //we'll misuse this for the viterbi
  double backwards[16];
  double forwscale; //adding numerical stability
  double backscale; //adding numerical stability
  int ncols;
  int nlets;
  int letssize;
  int colssize;
  int* lets;
  int* cols;
  double* letserrrate;
  double* colserrrate;
  char backpointer[16]; //previous state for viterbi
} states;

char bbmap[256]; /* byte to bit map; from chars to [0,1,2,3,4] */
char letmap[4] = {'T', 'G', 'C', 'A'}; /* reverse of this map */


/* In order to understand any of the code below, you need to understand color-space;
   Specifically that LETTER ^ LETTER = COLOR : T (00) ^ C (10) = 2 (10), etc. 
   And that LETTER ^ COLOR = NEXTLETTER: T (00) ^ 3 (11) = A (11). */

/* compute prior probability of the node given the emissions. letters are thought to be at the
   left "side" of the pair emitted by the node */

double nodePrior(states* allstates, int i, int j) { //i  is state, j is node in the state
  double val = 0;
  double errrate;
  int let, col, k;
  for (k = 0; k < allstates[i].nlets; k++) {
    let = allstates[i].lets[k];
    errrate = allstates[i].letserrrate[k];
    if (left(j) == let) {
      val = val - log(1-errrate);
    }
    else {
      val = val - log(errrate/3.0);
    }
  }
  //fprintf(stdout, "%g\n", val);
  for (k = 0; k < allstates[i].ncols; k++) {
    col = allstates[i].cols[k];
    errrate = allstates[i].colserrrate[k];
    if ((left(j) ^ right(j)) == col) {
      val = val - log(1-errrate);
    }
    else {
      val = val - log(errrate/3.0);
    }
    //fprintf(stdout, "%d %g\n", col, val);
  }
  return val;
}

/* Little helper for debugging */

void printStates(states* allstates, int stateslen, char* result, double* distrib, FILE* stream) {
  int i,j,k;
  fprintf(stream, "\nCONTIG %d", stateslen);

  for (i=0; i< stateslen; i++) {
    fprintf(stream, "\nCOLORS[%d] ",i);
    for (k = 0; k < allstates[i].ncols; k++) {
            fprintf(stream, "%d ",allstates[i].cols[k]);
    }
  }
  for (i=0; i< stateslen; i++) {
    fprintf(stream, "\nFORWARDSS[%d] ",i);
    for (j=0; j< 16; j++) {
      fprintf(stream, "%.5g ",allstates[i].forwards[j]);
    }    
  }
  for (i=0; i< stateslen; i++) {
    fprintf(stream, "\nBACKWARDSS[%d] ",i);
    for (j=0; j< 16; j++) {
      fprintf(stream, "%.5g ",allstates[i].backwards[j]);
    }    
  }

  for (i=0; i< stateslen; i++) {
    fprintf(stream, "\nLETS[%d] ",i);
    for (k = 0; k < allstates[i].nlets; k++) {
      fprintf(stream, "%d ",allstates[i].lets[k]);
    }
  
    fprintf(stream, "%c",result[i]);
    fprintf(stream, " %.5g",distrib[i]);
  }
  
  fprintf(stream, "\n");
}

/*maximum posterior traceback */

char* post_traceback (states* allstates, int stateslen, double norm_px, double* poster) {
  char* result = (char*) calloc (stateslen + 1, 1);
  int i = 0, j, maxval;
  double distrib[4];

  for (i = 0; i < stateslen; i++) {
    for (j=0; j< 4; j++) distrib[j] = 0;
    for (j = 0; j < 16; j++) {
      //     fprintf(stderr, "%g %g %g\n",allstates[i].forwards[j], allstates[i].backwards[j], norm_px); 
      distrib[left(j)] += exp(-1 * (allstates[i].forwards[j] + allstates[i].backwards[j] + allstates[i].forwscale + allstates[i].backscale - norm_px)); 
      //      fprintf(stderr, "distrib[%d,%d] = %g\n", i,j, exp(-1 * (allstates[i].forwards[j] + allstates[i].backwards[j] + allstates[i].forwscale + allstates[i].backscale - norm_px)));
    }
    maxval = 0;
    for (j=0; j< 4; j++)  {
      //      fprintf(stderr, "let_distrib[%d,%d] = %g\n", i,j, distrib[j]);

      if (distrib[j] >distrib[maxval]) 
	maxval = j;
    }
    //    fprintf (stderr, "\n");
    poster[i] = distrib[maxval];
    if (poster[i] > confrate) {
      result[i] = letmap[maxval];
    }
    else {
      result[i] = 'N';
    }
  }

  return result;
}


/*viterbi traceback */

char* vit_traceback (states* allstates, int stateslen) {
  char* result = (char*) calloc (stateslen + 1, 1);
  int i,j;
  int minval, prev;

  for (i = stateslen -1; i >= 0; i--) {
    minval = 0;
    for (j = 0; j< 16; j++) {
      if (allstates[i].forwards[j] < allstates[i].forwards[minval]) {
	minval = j;
      }
    }
    prev = allstates[i].backpointer[minval];
    if (i && (left(minval) != right (prev))) {
      fprintf (stderr, "BACKTRACE error %d %d %d\n", i, minval, prev);
      exit(2);
    }
    result[i] = letmap[left(minval)];
  }
  return result;
}


void viterbi (states* allstates, int stateslen) {
  int i,j,k,let,col;
  int minback;
  double valback;
  double val;
  i = 0;
  for (j = 0; j < 16; j++) {
    allstates[i].forwards[j] = nodePrior(allstates,i,j);
  }

  for (i=1; i < stateslen; i++) {    
    for (j = 0; j < 16; j++) {
      allstates[i].forwards[j] = nodePrior(allstates,i,j);

      minback = left(j);
      for (k = 1; k < 16; k++) {
	if (left(j) == right(k)) {
	  if (allstates[i-1].forwards[k] < allstates[i-1].forwards[minback]) {
	    minback = k;
	  }
	}
      }
      valback = allstates[i-1].forwards[minback];
      allstates[i].forwards[j] += valback;
      allstates[i].backpointer[j] = minback;
    }
  }
}

double do_backwards (states* allstates, int stateslen) {
  int i,j,k,let,col;
  double val;
  
  i = stateslen-1;
  allstates[i].backscale = 999999999;
  for (j = 0; j < 16; j++) {
    allstates[i].backwards[j] = neglogsixteenth;
    allstates[i].backscale = MIN2 (allstates[i].backscale, allstates[i].backwards[j]);
  }
  for (j = 0; j < 16; j++) {
    allstates[i].backwards[j] -= allstates[i].backscale;
  }

  for (i = stateslen-2; i >=0; i--) {    
    allstates[i].backscale = 999999999;
    for (j = 0; j < 16; j++) {
      for (k = 0; k < 16; k++) {
	if (right(j) == left(k)) {
	  val = nodePrior(allstates,i+1,k);
	  allstates[i].backwards[j] += exp(-1*(val + allstates[i+1].backwards[k]));
	}
      }
      //      fprintf(stdout, "bw was [%d, %d] = %g\n", i, j, allstates[i].backwards[j]); 

      allstates[i].backwards[j] = -log(allstates[i].backwards[j]) + neglogfourth;
      allstates[i].backscale = MIN2 (allstates[i].backscale, allstates[i].backwards[j]);
    }
    for (j = 0; j < 16; j++) {

      allstates[i].backwards[j] -= allstates[i].backscale;
      //      fprintf(stdout, "bw is [%d, %d] = %g\n", i, j, allstates[i].backwards[j]); 
    }
    allstates[i].backscale += allstates[i+1].backscale;

  }
  val = 0;
  i = 0;
  for (j = 0; j < 16; j++) {
    val += exp(-1*(allstates[i].backwards[j] + nodePrior(allstates,i,j) + neglogsixteenth));
  }
  return -log(val) + allstates[0].backscale;
}

double do_forwards (states* allstates, int stateslen) {
  int i,j,k,let,col;
  double val;
  
  i = 0; j = 0;
  allstates[i].forwscale = 999999999;
  for (j = 0; j < 16; j++) {
    allstates[i].forwards[j] = nodePrior(allstates,i,j) + neglogsixteenth;
    allstates[i].forwscale = MIN2 (allstates[i].forwscale, allstates[i].forwards[j]);
  }
  for (j = 0; j < 16; j++) {
    allstates[i].forwards[j] -= allstates[i].forwscale;
  }

  for (i=1; i < stateslen; i++) {    
    allstates[i].forwscale = 999999999;
    for (j = 0; j < 16; j++) {
      val = nodePrior(allstates,i,j);
      for (k = 0; k < 16; k++) {
	if (left(j) == right(k)) {
	  allstates[i].forwards[j] += exp(-1*(allstates[i-1].forwards[k]));
	}
      }
      allstates[i].forwards[j] = val + neglogfourth - log(allstates[i].forwards[j]);
      allstates[i].forwscale = MIN2 (allstates[i].forwscale, allstates[i].forwards[j]);
    }
    for (j = 0; j < 16; j++) {
      allstates[i].forwards[j] -= allstates[i].forwscale;
    }
    allstates[i].forwscale += allstates[i-1].forwscale;
  }

  val = 0;
  i = stateslen-1;
  for (j = 0; j < 16; j++) {
    val += exp(-1*(allstates[i].forwards[j] + neglogsixteenth));
  }
  return -log(val)+ allstates[i].forwscale;
}

double forward_backward (states* allstates, int stateslen) {
  double no1, no2;
  no1 = do_forwards(allstates, stateslen);
  no2 = do_backwards(allstates, stateslen);
  //fprintf (stderr, "SANITY CHECK: no1 == no2 %g %g\n", no1, no2);
  // don't really want a hard assert due to precision issues
  return no1;
}


/* updates emissions of states of hmm for a given read */

void updateStates (states* allstates, int pos, char* read, char* qual) {
  int readlen = strlen(read)-1;
  int first = bbmap[read[0]] ^ bbmap[read[1]];
  int i=0, base, rdpos;

  if (!isalpha(read[i]) || bbmap[read[i]] < 0) {
    fprintf(stderr, "Parse error read %s, pos %d\n",read, i);
    exit (2);
  }
  for (i=1; i <= readlen; i++) {
    if (!isdigit(read[i]) || bbmap[read[i]] < 0) {
      fprintf(stderr, "Parse error read %s, pos %d\n",read, i);
      exit (2);
    }
  }
  /* pos < 0 means read on negative strand */

  if (pos < 0) {
    first = first ^ 3; /* revcomp the letter*/
    base = -pos + readlen - baseoffset - 1;
  }
  else {
    base = pos-baseoffset;
  }

  if (allstates[base].nlets ==  allstates[base].letssize) {
    allstates[base].lets = realloc(allstates[base].lets, sizeof(int) * (allstates[base].letssize * 2 + 10));
    allstates[base].letserrrate = realloc(allstates[base].letserrrate, sizeof(double) * (allstates[base].letssize * 2 + 10));
    allstates[base].letssize = allstates[base].letssize * 2 + 10;
  }
  allstates[base].lets[allstates[base].nlets] = first;
  if (qual == NULL)
  {
      allstates[base].letserrrate[allstates[base].nlets] = glob_errrate;
  }
  else
  {
      allstates[base].letserrrate[allstates[base].nlets] = 
              pow(10.0, -((double)(qual[base] - '!'))/10.0);
  }
  allstates[base].nlets++;

  for (i = 2;  i <= readlen; i++) {
    if (pos < 0) {
      base = -pos + i - baseoffset - 2;
      rdpos = readlen - i + 2;
    }
    else {
      base = pos+i-2 - baseoffset;
      rdpos = i;
    }
    if (allstates[base].ncols ==  allstates[base].colssize) {
      allstates[base].cols = realloc(allstates[base].cols, sizeof(int) * (allstates[base].colssize * 2 + 10));
      allstates[base].colserrrate = realloc(allstates[base].colserrrate, sizeof(double) * (allstates[base].colssize * 2 + 10));
      allstates[base].colssize = allstates[base].colssize * 2 + 10;
    }
    allstates[base].cols[allstates[base].ncols] = bbmap[read[rdpos]];
    if (qual == NULL)
    {
        allstates[base].colserrrate[allstates[base].ncols] = glob_errrate;
    }
    else
    {
        allstates[base].colserrrate[allstates[base].ncols] = 
                pow(10.0, -((double)(qual[base] - '!'))/10.0);
    }
    allstates[base].ncols++;
  }  
}

char* parseContig(char* buffer) {
  int len, pos;
  char eval, contigname[1024], read[1024], qual[1024];
  int i,t = sscanf(buffer, "Contig %s %d", contigname, &len);
  states* allstates;
  char* result;
  double norm_px;
  double logerr;
  double* posteriors;
  if (t != 2) {
    fprintf(stderr, "Parse error %s\n", buffer);
    exit (2);
  }
  if (len < 0 || len > 10000000) {
    fprintf(stderr, "Length out of range %d\n", len);
    exit (2);
  }
  allstates = (states*) calloc (len, sizeof(states));
  posteriors = (double*) calloc (len, sizeof(double));
  buffer = fgets(buffer, 1024, stdin);
  while (buffer && !strstr(buffer, "Contig")) {
    t= sscanf (buffer, "%d %s %s", &pos, read, qual);
    if (t < 2) {
      fprintf(stderr, "Parse error %s\n", buffer);
      exit (2);
    }
    if (t == 2)
        updateStates(allstates, pos, read, NULL);
    else
        updateStates(allstates, pos, read, qual);
    buffer = fgets(buffer, 1024, stdin);
  }
  /* Viterbi*/
  if (vitorfb) {
    viterbi(allstates,len);
    result = vit_traceback(allstates, len);
  }
  else {
  /* Forward-Backward */
    norm_px = forward_backward(allstates,len);
    result = post_traceback(allstates, len, norm_px, posteriors);
  }

  fprintf(stdout, "%s %d\n", contigname, len);
  fprintf(stdout, "%s\n", result);
  if (!vitorfb) {
    /* print confidence values -- log(errorlikelihoods) -- if we ran forw-backw */

    for (i = 0; i < len; i++) {
      if (1-posteriors[i] > .0000000001) {
	logerr = -log(1-posteriors[i]);
	t = MIN2(((int)logerr),9);
      }
      else { t = 9; }
      eval = '0' + t;
      fprintf(stdout, "%c", eval);
    }
    fprintf(stdout,"\n");
  }
  for (i=0; i< len; i++) {
    free(allstates[i].lets);
    free(allstates[i].cols);
    free(allstates[i].letserrrate);
    free(allstates[i].colserrrate);
  }
  free (allstates);
  free (posteriors);
  free (result);
  return buffer;
}

void init() {
  int i;
  for (i = 0; i < 256; i++) {
    bbmap[i] = -1;
  }
  bbmap['T'] = 0;
  bbmap['G'] = 1;    
  bbmap['C'] = 2;    
  bbmap['A'] = 3;
  bbmap['0'] = 0;
  bbmap['1'] = 1;    
  bbmap['2'] = 2;    
  bbmap['3'] = 3;
  neglogsixteenth = - log (1.0/16.0);
  neglogfourth = -log (1.0/4.0);
}

void usage(const char * progname) {
    fprintf(stderr,
            "Usage: %s [OPTIONS]\n"
            "\n"
            "Options:\n"
            "   -e --error-rate=<#>     The expected per-color error rate.\n"
            "                           (default: %f)\n"
            "   -c --conf-level=<#>     The confidence level required for a base\n"
            "                           to be called. (default: %f)\n"
            "   -? --help               Prints this message and exits.\n"
            "\n"
            "Notes:\n"
            "   Input is read from stdin, and output is written to stdout.\n"
            "\n", 
                progname,
                DEFAULT_ERROR_RATE,
                DEFAULT_CONF_LEVEL);
}

int main (int argc, char** argv) {

    char * progname;
    char buffer[1024];
    char* test;
    int opt;

    progname = basename(argv[0]);

    /* Viterbi disabled; use FB -tsmith */
    vitorfb = 0;
    glob_errrate = DEFAULT_ERROR_RATE;
    confrate = DEFAULT_CONF_LEVEL;

    /* Parse arguments */
    for (;;)
    {
        opt = getopt_long(argc, argv, CMDLINE_FLAGS, CMDLINE_OPTIONS, NULL);
        if (opt == -1)
            break;
        
        switch (opt)
        {
            /* Per-color error rate */
            case 'e':
                glob_errrate = atof(optarg);
                if (glob_errrate <= 0.0)
                {
                    fprintf(stderr, "error-rate must be positive\n");
                    goto print_usage;
                }
                break;
            case 'c':
                confrate = atof(optarg);
                if (confrate <= 0.0)
                {
                    fprintf(stderr, "conf-level must be positive\n");
                    goto print_usage;
                }
                break;
            /* Help, or unknown option */
            case '?':
            default:
                goto print_usage;
        }
    }
    
    init();
    
    test = fgets(buffer, 1024, stdin);
    while (test && strstr(test, "Contig")) {
      test = parseContig(test);
    }
    if (!feof(stdin)) {
      fprintf(stderr, "Parse Error!!! %s\n", test);
      exit (1);
    }
    return 0;
    
print_usage:
    usage(progname);
}



