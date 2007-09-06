# $Id$
CFLAGS=-Wall -Werror -O3 -static
LDFLAGS=-lm

all: rmapper-cs rmapper-ls

rmapper-cs: rmapper/rmapper.c common/fasta-cs.o rmapper/sw-vector.o \
    rmapper/sw-full-cs.o rmapper/sw-full-ls.o common/util.o
	$(CC) $(CFLAGS) -DUSE_COLOURS -o $@ $+ $(LDFLAGS)

rmapper-ls: rmapper/rmapper.c common/fasta-cs.o rmapper/sw-vector.o \
    rmapper/sw-full-cs.o rmapper/sw-full-ls.o common/util.o
	$(CC) $(CFLAGS) -o $@ $+ $(LDFLAGS)

rmapper/sw-vector.o: rmapper/sw-vector.c rmapper/rmapper.h
	$(CC) $(CFLAGS) -c -o $@ $<

rmapper/sw-full-cs.o: rmapper/sw-full-cs.c rmapper/rmapper.h
	$(CC) $(CFLAGS) -c -o $@ $<

rmapper/sw-full-ls.o: rmapper/sw-full-ls.c rmapper/rmapper.h
	$(CC) $(CFLAGS) -c -o $@ $<

common/fasta-cs.o: common/fasta.c common/fasta.h
	$(CC) $(CFLAGS) -c -o $@ $<

common/fasta-ls.o: common/fasta.c common/fasta.h
	$(CC) $(CFLAGS) -c -o $@ $<

common/util.o: common/util.c common/util.h
	$(CC) $(CFLAGS) -c -o $@ $<

clean:
	rm -f rmapper-cs rmapper-ls
	find . -name '*.o' |xargs rm -f
	find . -name  '*.core' |xargs rm -f
