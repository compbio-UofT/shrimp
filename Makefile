# $Id$
CFLAGS=-Wall -Werror -O3 -static
LDFLAGS=-lm

all: bin/rmapper-cs bin/rmapper-ls bin/colourise bin/finalpass bin/revcmpl \
    bin/splitreads bin/splittigs

#
# rmapper/
#

bin/rmapper-cs: rmapper/rmapper.c common/fasta-cs.o rmapper/sw-vector.o \
    rmapper/sw-full-cs.o rmapper/sw-full-ls.o common/util.o
	$(CC) $(CFLAGS) -DUSE_COLOURS -o $@ $+ $(LDFLAGS)

bin/rmapper-ls: rmapper/rmapper.c common/fasta-cs.o rmapper/sw-vector.o \
    rmapper/sw-full-cs.o rmapper/sw-full-ls.o common/util.o
	$(CC) $(CFLAGS) -o $@ $+ $(LDFLAGS)

rmapper/sw-vector.o: rmapper/sw-vector.c rmapper/rmapper.h
	$(CC) $(CFLAGS) -c -o $@ $<

rmapper/sw-full-cs.o: rmapper/sw-full-cs.c rmapper/rmapper.h
	$(CC) $(CFLAGS) -c -o $@ $<

rmapper/sw-full-ls.o: rmapper/sw-full-ls.c rmapper/rmapper.h
	$(CC) $(CFLAGS) -c -o $@ $<

#
# colourise
#

bin/colourise: colourise/colourise.c colourise/fasta.o
	$(CC) $(CFLAGS) -o $@ $+

bin/colourise/fasta.o: colourise/fasta.c colourise/fasta.h
	$(CC) $(CFLAGS) -c -o $@ $<

#
# finalpass 
#

bin/finalpass: finalpass/finalpass.c common/lookup.o common/red_black_tree.o \
    common/util.o
	$(CC) $(CFLAGS) -o $@ $+ $(LDFLAGS)

#
# finalprint
#

#
# revcmpl
#

bin/revcmpl: revcmpl/revcmpl.c
	$(CC) $(CFLAGS) -o $@ $+

#
# splitreads
#

bin/splitreads: splitreads/splitreads.c
	$(CC) $(CFLAGS) -o $@ $+

#
# splittigs
#

bin/splittigs: splittigs/splittigs.c
	$(CC) $(CFLAGS) -o $@ $+

#
# common/
#

common/fasta-cs.o: common/fasta.c common/fasta.h
	$(CC) $(CFLAGS) -c -o $@ $<

common/fasta-ls.o: common/fasta.c common/fasta.h
	$(CC) $(CFLAGS) -c -o $@ $<

common/lookup.o: common/lookup.c common/lookup.h
	$(CC) $(CFLAGS) -c -o $@ $<

common/red_black_tree.o: common/red_black_tree.c common/red_black_tree.h
	$(CC) $(CFLAGS) -c -o $@ $<

common/util.o: common/util.c common/util.h
	$(CC) $(CFLAGS) -c -o $@ $<

#
# cleanup
#

clean:
	rm -f bin/rmapper-cs bin/rmapper-ls bin/colourise bin/finalpass \
	    bin/revcmpl bin/splitreads bin/splittigs
	find . -name '*.o' |xargs rm -f
	find . -name  '*.core' |xargs rm -f
