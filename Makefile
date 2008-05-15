# $Id$
CC=/home/rumble/intel/cce/10.0.025/bin/icc
CFLAGS=-Kc++ -wd981,1572 -axP -O3 -ipo -D__STDC_FORMAT_MACROS -D__STDC_LIMIT_MACROS -static

#CC=g++
#CFLAGS=-O3 -mmmx -msse -msse2 -DSHRIMP_BUGS=OFF -D__STDC_FORMAT_MACROS -D__STDC_LIMIT_MACROS
#CFLAGS=-pg -mmmx -msse -msse2 -DSHRIMP_BUGS=OFF -D__STDC_FORMAT_MACROS -D__STDC_LIMIT_MACROS
#CFLAGS=-g -mmmx -msse -msse2 -DSHRIMP_BUGS=OFF -D__STDC_FORMAT_MACROS -D__STDC_LIMIT_MACROS
LDFLAGS=-lm
CP=cp
MV=mv

#all: bin/rmapper bin/colourise bin/probcalc \
#    bin/prettyprint-cs bin/prettyprint-ls bin/mergehits-cs bin/mergehits-ls \
#    bin/revcmpl bin/splitreads bin/splittigs

#
# rmapper/
#

bin/rmapper: rmapper/rmapper.c common/fasta.o common/dag_align.o \
     common/dag_glue.o common/dag_kmers.o common/sw-vector.o \
     common/sw-full-cs.o common/sw-full-ls.o common/output.o common/util.o
	$(CC) $(CFLAGS) -o $@ $+ $(LDFLAGS)
	$(CP) bin/rmapper bin/rmapper-cs
	$(CP) bin/rmapper bin/rmapper-dag
	$(MV) bin/rmapper bin/rmapper-ls

#
# colourise
#

#bin/colourise: colourise/colourise.c colourise/fasta.o
#	$(CC) $(CFLAGS) -o $@ $+

#bin/colourise/fasta.o: colourise/fasta.c colourise/fasta.h
#	$(CC) $(CFLAGS) -c -o $@ $<

#
# probcalc 
#

#bin/probcalc: probcalc/probcalc.c common/fasta.o common/dynhash.o \
#    common/input.o common/output.o common/util.o
#	$(CC) $(CFLAGS) -o $@ $+ $(LDFLAGS)

#
# prettyprint
#

#bin/prettyprint-cs: prettyprint/prettyprint.c common/fasta.o common/dynhash.o \
#    common/sw-full-cs.o common/sw-full-ls.o common/input.o common/output.o \
#    common/util.o
#	$(CC) $(CFLAGS) -DUSE_COLOURS -o $@ $+ $(LDFLAGS)

#bin/prettyprint-ls: prettyprint/prettyprint.c common/fasta.o common/dynhash.o \
#    common/sw-full-cs.o common/sw-full-ls.o common/input.o common/output.o \
#    common/util.o
#	$(CC) $(CFLAGS) -o $@ $+ $(LDFLAGS)

#
# mergehits
#

#bin/mergehits-cs: mergehits/mergehits.c common/fasta.o common/dynhash.o \
#    common/input.o common/output.o common/util.o
#	$(CC) $(CFLAGS) -DUSE_COLOURS -o $@ $+ $(LDFLAGS)

#bin/mergehits-ls: mergehits/mergehits.c common/fasta.o common/dynhash.o \
#    common/input.o common/output.o common/util.o
#	$(CC) $(CFLAGS) -o $@ $+ $(LDFLAGS)

#
# revcmpl
#

#bin/revcmpl: revcmpl/revcmpl.c
#	$(CC) $(CFLAGS) -o $@ $+

#
# splitreads
#

#bin/splitreads: splitreads/splitreads.c
#	$(CC) $(CFLAGS) -o $@ $+

#
# splittigs
#

#bin/splittigs: splittigs/splittigs.c
#	$(CC) $(CFLAGS) -o $@ $+

#
# common/
#

common/fasta.o: common/fasta.c common/fasta.h
	$(CC) $(CFLAGS) -c -o $@ $<

common/dag_align.o: common/dag_align.cpp common/dag_align.h
	$(CC) $(CFLAGS) -c -o $@ $<

common/dag_glue.o: common/dag_glue.cpp common/dag_glue.h
	$(CC) $(CFLAGS) -c -o $@ $<

common/dag_kmers.o: common/dag_kmers.cpp
	$(CC) $(CFLAGS) -c -o $@ $<

common/dynhash.o: common/dynhash.c common/dynhash.h
	$(CC) $(CFLAGS) -c -o $@ $<

common/input.o: common/input.c common/input.h
	$(CC) $(CFLAGS) -c -o $@ $<

common/output.o: common/output.c common/output.h
	$(CC) $(CFLAGS) -c -o $@ $<

common/sw-full-cs.o: common/sw-full-cs.c common/sw-full-cs.h common/sw-full-common.h
	$(CC) $(CFLAGS) -c -o $@ $<

common/sw-full-ls.o: common/sw-full-ls.c common/sw-full-ls.h common/sw-full-common.h
	$(CC) $(CFLAGS) -c -o $@ $<

common/sw-vector.o: common/sw-vector.c common/sw-vector.h
	$(CC) $(CFLAGS) -c -o $@ $<

common/util.o: common/util.c common/util.h
	$(CC) $(CFLAGS) -c -o $@ $<

#
# cleanup
#

clean:
	rm -f bin/rmapper* bin/colourise bin/probcalc \
	    bin/prettyprint-cs bin/prettyprint-ls bin/mergehits-cs \
	    bin/mergehits-ls bin/revcmpl bin/splitreads bin/splittigs
	find . -name '*.o' |xargs rm -f
	find . -name  '*.core' |xargs rm -f
