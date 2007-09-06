# $Id$
CFLAGS=-Wall -Werror -O3 -static
LDFLAGS=-lm

all: rmapper-cs rmapper-ls

rmapper-cs: rmapper.c fasta-cs.o sw-vector.o sw-full-cs.o sw-full-ls.o util.o
	$(CC) $(CFLAGS) -DUSE_COLOURS -o $@ $+ $(LDFLAGS)

rmapper-ls: rmapper.c fasta-ls.o sw-vector.o sw-full-cs.o sw-full-ls.o util.o
	$(CC) $(CFLAGS) -o $@ $+ $(LDFLAGS)

fasta-cs.o: fasta.c fasta.h
	$(CC) $(CFLAGS) -c -o $@ $<

fasta-ls.o: fasta.c fasta.h
	$(CC) $(CFLAGS) -c -o $@ $<

sw-vector.o: sw-vector.c rmapper.h
	$(CC) $(CFLAGS) -c -o $@ $<

sw-full-cs.o: sw-full-cs.c rmapper.h
	$(CC) $(CFLAGS) -c -o $@ $<

sw-full-ls.o: sw-full-ls.c rmapper.h
	$(CC) $(CFLAGS) -c -o $@ $<

util.o: util.c util.h
	$(CC) $(CFLAGS) -c -o $@ $<

clean:
	rm -f *.o *.core rmapper-cs rmapper-ls
