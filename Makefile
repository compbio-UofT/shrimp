# $Id$
ifndef CXXFLAGS
CXXFLAGS=-O3 -mmmx -msse -msse2 -Wall -Werror
endif
override CXXFLAGS+=-D__STDC_FORMAT_MACROS -D__STDC_LIMIT_MACROS -DSHRIMP_BUGS=OFF
LDFLAGS=-lm -lz
LN=ln

# statically compile on all but OS X.
UNAME=$(shell uname -s)
ifneq ($(UNAME),Darwin)
override CXXFLAGS+=-static
endif

all: bin/rmapper bin/probcalc bin/prettyprint bin/mergehits

#
# rmapper/
#

bin/rmapper: rmapper/rmapper.c common/fasta.o common/dag_align.o \
    common/dag_glue.o common/dag_kmers.o common/sw-vector.o \
    common/sw-full-cs.o common/sw-full-ls.o common/input.o \
    common/output.o common/util.o
	$(CXX) $(CXXFLAGS) -o $@ $+ $(LDFLAGS)
	$(LN) -sf rmapper bin/rmapper-cs
	$(LN) -sf rmapper bin/rmapper-hs
	$(LN) -sf rmapper bin/rmapper-ls

#
# probcalc 
#

bin/probcalc: probcalc/probcalc.c common/fasta.o common/dynhash.o \
    common/input.o common/output.o common/util.o
	$(CXX) $(CXXFLAGS) -o $@ $+ $(LDFLAGS)

#
# prettyprint
#

bin/prettyprint: prettyprint/prettyprint.c common/fasta.o common/dynhash.o \
    common/sw-full-cs.o common/sw-full-ls.o common/input.o common/output.o \
    common/util.o
	$(CXX) $(CXXFLAGS) -o $@ $+ $(LDFLAGS)
	$(LN) -sf prettyprint bin/prettyprint-cs
	$(LN) -sf prettyprint bin/prettyprint-hs
	$(LN) -sf prettyprint bin/prettyprint-ls

#
# mergehits
#

bin/mergehits: mergehits/mergehits.c common/fasta.o common/dynhash.o \
    common/input.o common/output.o common/util.o
	$(CXX) $(CXXFLAGS) -o $@ $+ $(LDFLAGS)
	$(LN) -sf mergehits bin/mergehits-cs
	$(LN) -sf mergehits bin/mergehits-hs
	$(LN) -sf mergehits bin/mergehits-ls

#
# common/
#

common/fasta.o: common/fasta.c common/fasta.h
	$(CXX) $(CXXFLAGS) -c -o $@ $<

common/dag_align.o: common/dag_align.cpp common/dag_align.h
	$(CXX) $(CXXFLAGS) -c -o $@ $<

common/dag_glue.o: common/dag_glue.cpp common/dag_glue.h
	$(CXX) $(CXXFLAGS) -c -o $@ $<

common/dag_kmers.o: common/dag_kmers.cpp
	$(CXX) $(CXXFLAGS) -c -o $@ $<

common/dynhash.o: common/dynhash.c common/dynhash.h
	$(CXX) $(CXXFLAGS) -c -o $@ $<

common/input.o: common/input.c common/input.h
	$(CXX) $(CXXFLAGS) -c -o $@ $<

common/output.o: common/output.c common/output.h
	$(CXX) $(CXXFLAGS) -c -o $@ $<

common/sw-full-cs.o: common/sw-full-cs.c common/sw-full-cs.h common/sw-full-common.h
	$(CXX) $(CXXFLAGS) -c -o $@ $<

common/sw-full-ls.o: common/sw-full-ls.c common/sw-full-ls.h common/sw-full-common.h
	$(CXX) $(CXXFLAGS) -c -o $@ $<

common/sw-vector.o: common/sw-vector.c common/sw-vector.h
	$(CXX) $(CXXFLAGS) -c -o $@ $<

common/util.o: common/util.c common/util.h
	$(CXX) $(CXXFLAGS) -c -o $@ $<

#
# cleanup
#

clean:
	rm -f bin/rmapper* bin/colourise bin/probcalc \
	    bin/prettyprint* bin/mergehits*
	find . -name '*.o' |xargs rm -f
	find . -name  '*.core' |xargs rm -f
	find . -name '*.pyc' |xargs rm -f
