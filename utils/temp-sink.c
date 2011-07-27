#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <getopt.h>
#include <errno.h>
#include <unistd.h>
#include <string.h>

#define DEF_BLOCK_SIZE (1048576)
#define MIN(a, b) ((a) < (b) ? (a) : (b))

struct option getopt_long_options[] = {
  {"block-size", 1, 0, 'b'},
  {"help", 0, 0, '?'}
};

char const getopt_short_options[] = "b:?";

FILE * tmp_file, * dest_file;
char * prog_name;
char * buffer, * buffer2;
size_t block_size = DEF_BLOCK_SIZE;
size_t total_tmp, total_dest, total_in, total_out, total_out2;


void
logit(FILE * f, int stop, int exit_code, char const * msg, ...)
{
  va_list fmtargs;

  va_start(fmtargs, msg);
  vfprintf(f, msg, fmtargs);
  va_end(fmtargs);

  if (stop)
    exit(exit_code);
}


void
flush_tmp_to_dest()
{
  size_t crt, total_buff2_in, total_buff2_out;

  rewind(tmp_file);
  crt = 0;

  total_buff2_in = fread(buffer2, sizeof(char), MIN(total_tmp - crt, block_size), tmp_file);
  if (ferror(tmp_file))
    logit(stderr, 1, 1, "error reading tmp_file (%s)\n", strerror(errno));

  while (total_buff2_in > 0) {
    total_buff2_out = fwrite(buffer2, sizeof(char), total_buff2_in, dest_file);
    crt += total_buff2_out;
    if (total_buff2_out != total_buff2_in)
      logit(stderr, 1, 1, "error writing dest_file (%s)\n", strerror(errno));

    if (feof(tmp_file))
      break;

    total_buff2_in = fread(buffer2, sizeof(char), block_size, tmp_file);
    if (ferror(tmp_file))
      logit(stderr, 1, 1, "error reading tmp_file (%s)\n", strerror(errno));
  }

  if (crt != total_tmp)
    logit(stderr, 1, 1, "error: something went wrong (total_tmp=%lld crt=%lld)\n", (long long)total_tmp, (long long)crt);

  total_dest += crt;

  rewind(tmp_file);
  if (ftruncate(fileno(tmp_file), 0))
    logit(stderr, 1, 1, "error truncating tmp_file (%s)\n", strerror(errno));
  total_tmp = 0;
}


void
usage()
{
  fprintf(stderr, "\nSynopsis:\n\n");
  fprintf(stderr, "  Copy stdin to dest_file via tmp_file.\n\n");
  fprintf(stderr, "  If tmp_file runs out of space, the program flushes its contents to dest_file,\n");
  fprintf(stderr, "  truncates tmp_file to zero length, and continues as before.\n\n");
  fprintf(stderr, "  If successful, all of stdin is in dest_file, and tmp_file is empty.\n\n");
  fprintf(stderr, "Options:\n\n");
  fprintf(stderr, "  -b | --block-size <size>\n");
  fprintf(stderr, "    Set block size for read/write operations.\n");
}


int
main(int argc, char * argv[])
{
  char ch;
  int bypass_tmp;
  size_t test_out;

  prog_name = argv[0];

  while ((ch = getopt_long(argc, argv, getopt_short_options, getopt_long_options, NULL)) != -1) {
    switch (ch) {
    case 'b':
      block_size = atoi(optarg);
      if (block_size < 64)
	logit(stderr, 1, 1, "error: block size to small (\"%s\" => %d)\n", optarg, block_size);
      break;
    case '?':
    default:
      logit(stderr, 0, 0, "Use: %s [tmp_file] [dest_file]\n", prog_name);
      usage();
      exit(0);
    }
  }
  argc -= optind;
  argv += optind;

  if (argc < 2) {
    logit(stderr, 1, 1, "Use: %s [tmp_file] [dest_file]\n", prog_name);
  }

  buffer = (char *)malloc(block_size * sizeof(char));
  if (buffer == NULL) 
    logit(stderr, 1, 1, "error: could not allocate buffer of size %d (%s)\n", block_size, strerror(errno));
  buffer2 = (char *)malloc(block_size * sizeof(char));
  if (buffer2 == NULL) 
    logit(stderr, 1, 1, "error: could not allocate buffer of size %d (%s)\n", block_size, strerror(errno));

  tmp_file = fopen(argv[0], "w+");
  if (tmp_file == NULL)
    logit(stderr, 1, 1, "error: could not open file %s (%s)\n", argv[0], strerror(errno));

  dest_file = fopen(argv[1], "w+");
  if (dest_file == NULL)
    logit(stderr, 1, 1, "error: could not open file %s (%s)\n", argv[1], strerror(errno));

  // test whether tmp_file has room for one block
  test_out = fwrite(buffer, sizeof(char), block_size, tmp_file);
  if (test_out != block_size) {
    if (errno != ENOSPC)
      logit(stderr, 1, 1, "error writing tmp_file (%s)\n", strerror(errno));

    bypass_tmp = 1;
  } else {
    bypass_tmp = 0;
  }
  
  rewind(tmp_file);
  if (ftruncate(fileno(tmp_file), 0))
    logit(stderr, 1, 1, "error truncating tmp_file (%s)\n", strerror(errno));

  if (!bypass_tmp)
    {
      total_tmp = total_dest = 0;

      total_in = fread(buffer, sizeof(char), block_size, stdin);
      if (ferror(stdin))
	logit(stderr, 1, 1, "error reading stdin (%s)\n", strerror(errno));

      while (total_in > 0) {
	total_out = fwrite(buffer, sizeof(char), total_in, tmp_file);
	total_tmp += total_out;
	if (total_out != total_in) {
	  if (errno != ENOSPC)
	    logit(stderr, 1, 1, "error writing tmp_file (%s)\n", strerror(errno));

	  logit(stderr, 0, 0, "warning: temp disk full (total_tmp:%lld)\n", (long long int)total_tmp);

	  // move tmp_file to dest_file
	  flush_tmp_to_dest();

	  // write remaining bytes in buffer directly to dest
	  total_out2 = fwrite(buffer + total_out, sizeof(char), total_in - total_out, dest_file);
	  total_dest += total_out2;
	  if (total_out2 != total_in - total_out)
	    logit(stderr, 1, 1, "error writing dest_file (%s)\n", strerror(errno));
	}

	if (feof(tmp_file))
	  break;

	total_in = fread(buffer, sizeof(char), block_size, stdin);
	if (ferror(stdin))
	  logit(stderr, 1, 1, "error reading stdin (%s)\n", strerror(errno));
      }
      if (total_tmp > 0)
	flush_tmp_to_dest();
    }
  else // bypass tmp_file entirely
    {
      logit(stderr, 0, 0, "warning: tmp_file cannot hold even a single block; bypassing it entirely\n");

      total_in = fread(buffer, sizeof(char), block_size, stdin);
      if (ferror(stdin))
	logit(stderr, 1, 1, "error reading stdin (%s)\n", strerror(errno));

      while (total_in > 0) {
      	total_out = fwrite(buffer, sizeof(char), total_in, dest_file);
	total_tmp += total_out;
	if (total_out != total_in) {
	  logit(stderr, 1, 1, "error writing dest_file (%s)\n", strerror(errno));
	}

	if (feof(tmp_file))
	  break;

	total_in = fread(buffer, sizeof(char), block_size, stdin);
	if (ferror(stdin))
	  logit(stderr, 1, 1, "error reading stdin (%s)\n", strerror(errno));
      }
    }

  fclose(tmp_file);
  fclose(dest_file);
  free(buffer);
  free(buffer2);

  return 0;
}
