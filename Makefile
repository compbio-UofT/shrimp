# $Id: Makefile,v 1.23 2009/06/16 23:26:20 rumble Exp $
#CXXFLAGS=-Kc++ -wd383,981,1572 -axP -O3 -ipo -openmp -DNDEBUG -static-intel
#CXXFLAGS=-Kc++ -wd383,981,1572 -axP -O3 -ipo -openmp -static-intel
#CXXFLAGS=-Kc++ -O2 -openmp -DNDEBUG -static-intel -g 

#CXXFLAGS=-fopenmp -Wall -Wno-deprecated -g -DDEBUG_KMERS -DDEBUG_HIT_LIST_CREATION -DDEBUG_HIT_LIST_PASS1 -DDEBUG_SW_FULL_CALLS -DDEBUG_ANCHOR_LIST 
ifndef CXXFLAGS
#CXXFLAGS=-g -p -mmmx -msse -msse2 -fopenmp -Wall -Wno-deprecated -DNDEBUG
#CXXFLAGS=-g -O3 -mmmx -msse -msse2 -fopenmp -Wall -Wno-deprecated -DNDEBUG
#CXXFLAGS=-mmmx -msse -msse2 -fopenmp -Wall -Wno-deprecated -DNDEBUG -DDEBUG_KMERS -DDEBUG_HIT_LIST_PASS1
#CXXFLAGS=-g -p -mmmx -msse -msse2 -fopenmp -Wall -Wno-deprecated  -DNDEBUG
CXXFLAGS=-g -fopenmp -Wall -Wno-deprecated 
#CXXFLAGS=-fopenmp -g -O1 -Wall -DNDEBUG
endif
override CXXFLAGS+=-D__STDC_FORMAT_MACROS -D__STDC_LIMIT_MACROS
#CXX=/opt/intel/cce/10.1.015/bin/icc
#CXX=/filer/misko/intel_suite/compilerpro-12.0.0.084/bin/intel64/icc

LD=$(CXX)
LDFLAGS=-lm -lz -lstdc++
LN=ln

all: bin/gmapper bin/probcalc bin/prettyprint bin/probcalc_mp \
    bin/shrimp_var bin/shrimp2sam utils/split-contigs bin/mergesam

#
# mapper /
#
mapper/mapper.o: mapper/mapper.c mapper/mapper.h
	 $(CXX) $(CXXFLAGS) -c -o $@ $<

#
# gmapper /
#
bin/gmapper: gmapper/gmapper.o common/fasta.o common/util.o \
    common/bitmap.o common/sw-vector.o common/sw-gapless.o common/sw-full-cs.o \
    common/sw-full-ls.o common/output.o common/anchors.o common/input.o \
    common/read_hit_heap.o
	$(LD) $(CXXFLAGS) -o $@ $+ $(LDFLAGS)
	$(LN) -sf gmapper bin/gmapper-cs
	$(LN) -sf gmapper bin/gmapper-ls

gmapper/gmapper.o: gmapper/gmapper.c common/bitmap.h gmapper/gmapper.h \
    common/debug.h common/f1-wrapper.h
	$(CXX) $(CXXFLAGS) -c -o $@ $<

bin/mergesam: mergesam/merge_sam.o mergesam/sam2pretty_lib.o mergesam/merge_sam_main.o mergesam/merge_sam_heap.o
	$(LD) $(CXXFLAGS) -o $@ $+ $(LDFLAGS)

mergesam/merge_sam.o: mergesam/merge_sam.c mergesam/merge_sam.h 
	$(CXX) $(CXXFLAGS) -c -o $@ $<

mergesam/merge_sam_main.o: mergesam/merge_sam_main.c  
	$(CXX) $(CXXFLAGS) -c -o $@ $<

mergesam/merge_sam_heap.o: mergesam/merge_sam_heap.c mergesam/merge_sam_heap.h
	$(CXX) $(CXXFLAGS) -c -o $@ $<

mergesam/sam2pretty_lib.o: mergesam/sam2pretty_lib.c mergesam/sam2pretty_lib.h
	$(CXX) $(CXXFLAGS) -c -o $@ $<
	

#
# probcalc 
#
bin/probcalc: probcalc/probcalc.o common/fasta.o common/dynhash.o \
    common/input.o common/output.o common/util.o
	$(LD) $(CXXFLAGS) -o $@ $+ $(LDFLAGS)

probcalc/probcalc.o: probcalc/probcalc.c
	$(CXX) $(CXXFLAGS) -c -o $@ $<

#
# probcalc_mp 
#
bin/probcalc_mp: probcalc_mp/probcalc_mp.o
	$(LD) $(CXXFLAGS) -o $@ $+ $(LDFLAGS)

probcalc_mp/probcalc_mp.o: probcalc_mp/probcalc_mp.c
	$(CXX) $(CXXFLAGS) -c -o $@ $<

#
# shrimp_var 
#
bin/shrimp_var: shrimp_var/shrimp_var.o
	$(LD) $(CXXFLAGS) -o $@ $+ $(LDFLAGS)

shrimp_var/shrimp_var.o: shrimp_var/shrimp_var.c
	$(CXX) $(CXXFLAGS) -c -o $@ $<


#
# prettyprint
#
bin/prettyprint: prettyprint/prettyprint.o common/fasta.o common/dynhash.o \
    common/sw-full-cs.o common/sw-full-ls.o common/input.o common/output.o \
    common/util.o common/anchors.o
	$(LD) $(CXXFLAGS) -o $@ $+ $(LDFLAGS)
	$(LN) -sf prettyprint bin/prettyprint-cs
	$(LN) -sf prettyprint bin/prettyprint-ls

prettyprint/prettyprint.o: prettyprint/prettyprint.c
	$(CXX) $(CXXFLAGS) -c -o $@ $<

#
# shrimp2sam
#
bin/shrimp2sam: shrimp2sam/shrimp2sam.o common/fasta.o common/dynhash.o \
    common/input.o common/output.o \
    common/util.o common/anchors.o
	$(LD) $(CXXFLAGS) -o $@ $+ $(LDFLAGS)
	$(LN) -sf prettyprint bin/prettyprint-ls

shrimp2sam/shrimp2sam.o: shrimp2sam/shrimp2sam.c
	$(CXX) $(CXXFLAGS) -c -o $@ $<

#
# bin/split-contigs
#
utils/split-contigs: utils/split-contigs.o common/fasta.o common/util.o
	$(LD) $(CXXFLAGS) -o $@ $+ $(LDFLAGS)

utils/split-contigs.o: utils/split-contigs.c
	$(CXX) $(CXXFLAGS) -c -o $@ $<

#
# common/
#
common/read_hit_heap.o: common/read_hit_heap.c common/read_hit_heap.h gmapper/gmapper.h
	$(CXX) $(CXXFLAGS) -c -o $@ $<

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

common/sw-full-cs.o: common/sw-full-cs.c common/sw-full-cs.h \
    common/sw-full-common.h common/util.h common/anchors.h
	$(CXX) $(CXXFLAGS) -c -o $@ $<

common/sw-full-ls.o: common/sw-full-ls.c common/sw-full-ls.h \
    common/sw-full-common.h common/util.h
	$(CXX) $(CXXFLAGS) -c -o $@ $<

common/sw-vector.o: common/sw-vector.c common/sw-vector.h common/util.h
	$(CXX) $(CXXFLAGS) -c -o $@ $<

common/sw-gapless.o: common/sw-gapless.c common/sw-gapless.h common/util.h
	$(CXX) $(CXXFLAGS) -c -o $@ $<

common/util.o: common/util.c common/util.h
	$(CXX) $(CXXFLAGS) -c -o $@ $<

common/anchors.o: common/anchors.c common/anchors.h common/util.h
	$(CXX) $(CXXFLAGS) -c -o $@ $<

common/bitmap.o: common/bitmap.c common/bitmap.h
	$(CXX) $(CXXFLAGS) -c -o $@ $<

#
# cleanup
#
clean:
	rm -f bin/colourise bin/probcalc bin/gmapper* \
	    bin/prettyprint* bin/probcalc_mp bin/shrimp_var \
	    bin/shrimp2sam utils/split-contigs bin/mergesam
	find . -name '*.o' |xargs rm -f
	find . -name  '*.core' |xargs rm -f
	find . -name '*.pyc' |xargs rm -f
