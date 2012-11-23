ifndef BUILD_TYPE
  BUILD_TYPE=production
endif

ifdef USE_ICC
  CXX=/opt/intel/cce/10.1.015/bin/icc
endif

ifndef CXXFLAGS
  ifeq ($(BUILD_TYPE), production)
    CXXFLAGS=-g -O2 -DNDEBUG
  else
    ifeq ($(BUILD_TYPE), testing)
      CXXFLAGS=-g -O2
    else
      CXXFLAGS=-g
    endif
  endif
  ifdef USE_ICC
    CXXFLAGS+=-Kc++ -wd383,981,1572 -axP -ipo -openmp -static-intel
  else
    CXXFLAGS+=-mmmx -msse -msse2 -fopenmp -Wall -Wno-deprecated
  endif
endif

GIT_VERSION=$(shell ./get_git_version)
override CXXFLAGS+=-D__STDC_FORMAT_MACROS -D__STDC_LIMIT_MACROS -DGIT_VERSION=$(GIT_VERSION)

LD=$(CXX)

UNAME := $(shell uname)
ifeq ($(UNAME), Darwin)
LDFLAGS=-lm -lz -lstdc++
else
LDFLAGS=-lm -lz -lstdc++ -lrt
endif 

LN=ln

all: bin/gmapper bin/probcalc bin/prettyprint bin/probcalc_mp \
    bin/shrimp_var bin/shrimp2sam utils/split-contigs bin/fasta2fastq bin/mergesam utils/temp-sink 

#
# mapper /
#
mapper/mapper.o: mapper/mapper.c mapper/mapper.h
	 $(CXX) $(CXXFLAGS) -c -o $@ $<

#
# gmapper /
#
bin/gmapper: gmapper/gmapper.o gmapper/seeds.o gmapper/genome.o gmapper/mapping.o gmapper/output.o \
    common/fasta.o common/util.o \
    common/bitmap.o common/sw-vector.o common/sw-gapless.o common/sw-full-cs.o \
    common/sw-full-ls.o common/output.o common/anchors.o common/input.o \
    common/read_hit_heap.o common/sw-post.o common/my-alloc.o common/gen-st.o
	$(LD) $(CXXFLAGS) -o $@ $+ $(LDFLAGS)
	$(LN) -sf gmapper bin/gmapper-cs
	$(LN) -sf gmapper bin/gmapper-ls

gmapper/gmapper.o: gmapper/gmapper.c common/bitmap.h gmapper/gmapper.h gmapper/gmapper-defaults.h \
    common/debug.h common/f1-wrapper.h common/version.h
	$(CXX) $(CXXFLAGS) -DCXXFLAGS="\"$(CXXFLAGS)\"" -c -o $@ $<

bin/lineindex: mergesam/lineindex.o mergesam/lineindex_lib.o mergesam/file_buffer.o
	$(LD) $(CXXFLAGS) -o $@ $+ $(LDFLAGS)

bin/fasta2fastq: mergesam/file_buffer.o mergesam/fasta_reader.o mergesam/fasta2fastq.o mergesam/lineindex_lib.o
	$(LD) $(CXXFLAGS) -o $@ $+ $(LDFLAGS)



mergesam/fasta_reader.o: mergesam/fasta_reader.c mergesam/fasta_reader.h
	$(CXX) $(CXXFLAGS) -c -o $@ $<

bin/mergesam: mergesam/file_buffer.o mergesam/sam2pretty_lib.o mergesam/mergesam_heap.o mergesam/mergesam.o mergesam/fastx_readnames.o mergesam/sam_reader.o mergesam/render.o
	$(LD) $(CXXFLAGS) -o $@ $+ $(LDFLAGS)

mergesam/mergesam.o: mergesam/mergesam.c
	$(CXX) $(CXXFLAGS) -c -o $@ $<

mergesam/mergesam_heap.o: mergesam/mergesam_heap.c  mergesam/mergesam_heap.h
	$(CXX) $(CXXFLAGS) -c -o $@ $<

mergesam/file_buffer.o: mergesam/file_buffer.c mergesam/file_buffer.h
	$(CXX) $(CXXFLAGS) -c -o $@ $<

mergesam/sam_reader.o: mergesam/sam_reader.c mergesam/sam_reader.h
	$(CXX) $(CXXFLAGS) -c -o $@ $<

mergesam/fastx_readnames.o: mergesam/fastx_readnames.c mergesam/fastx_readnames.h
	$(CXX) $(CXXFLAGS) -c -o $@ $<

mergesam/sam2pretty_lib.o: mergesam/sam2pretty_lib.c mergesam/sam2pretty_lib.h
	$(CXX) $(CXXFLAGS) -c -o $@ $<

mergesam/fasta2fastq.o: mergesam/fasta2fastq.c
	$(CXX) $(CXXFLAGS) -c -o $@ $<

mergesam/lineindex.o: mergesam/lineindex.c
	$(CXX) $(CXXFLAGS) -c -o $@ $<

mergesam/lineindex_lib.o: mergesam/lineindex_lib.c
	$(CXX) $(CXXFLAGS) -c -o $@ $<

mergesam/render.o: mergesam/render.c mergesam/render.h mergesam/sam2pretty_lib.h
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

prettyprint/prettyprint.o: prettyprint/prettyprint.c common/version.h
	$(CXX) $(CXXFLAGS) -c -o $@ $<

#
# shrimp2sam
#
bin/shrimp2sam: shrimp2sam/shrimp2sam.o common/fasta.o common/dynhash.o \
    common/input.o common/output.o \
    common/util.o common/anchors.o
	$(LD) $(CXXFLAGS) -o $@ $+ $(LDFLAGS)
	$(LN) -sf prettyprint bin/prettyprint-ls

shrimp2sam/shrimp2sam.o: shrimp2sam/shrimp2sam.c common/version.h
	$(CXX) $(CXXFLAGS) -c -o $@ $<

#
# utils/split-contigs
#
utils/split-contigs: utils/split-contigs.o common/fasta.o common/util.o
	$(LD) $(CXXFLAGS) -o $@ $+ $(LDFLAGS)

utils/split-contigs.o: utils/split-contigs.c
	$(CXX) $(CXXFLAGS) -c -o $@ $<

#
# utils/temp-sink
#
utils/temp-sink: utils/temp-sink.o
	$(LD) $(CXXFLAGS) -o $@ $+ $(LDFLAGS)

utils/temp-sink.o: utils/temp-sink.c
	$(CXX) $(CXXFLAGS) -c -o $@ $<

#
# gmapper/
#
gmapper/seeds.o: gmapper/seeds.c gmapper/seeds.h gmapper/gmapper.h
	$(LD) $(CXXFLAGS) -c -o $@ $<

gmapper/genome.o: gmapper/genome.c gmapper/genome.h gmapper/gmapper.h
	$(LD) $(CXXFLAGS) -c -o $@ $<

gmapper/mapping.o: gmapper/mapping.c gmapper/mapping.h gmapper/gmapper.h
	$(LD) $(CXXFLAGS) -c -o $@ $<

gmapper/output.o: gmapper/output.c gmapper/output.h gmapper/gmapper.h
	$(LD) $(CXXFLAGS) -c -o $@ $<

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

common/sw-post.o: common/sw-post.c common/sw-post.h common/util.h common/fasta.h common/sw-full-common.h
	$(CXX) $(CXXFLAGS) -c -o $@ $<

common/my-alloc.o: common/my-alloc.c common/my-alloc.h
	$(CXX) $(CXXFLAGS) -c -o $@ $<

common/gen-st.o: common/gen-st.c common/gen-st.h
	$(CXX) $(CXXFLAGS) -c -o $@ $<

#
# cleanup
#
clean:
	rm -f bin/colourise bin/probcalc bin/gmapper* \
	    bin/prettyprint* bin/probcalc_mp bin/shrimp_var \
	    bin/shrimp2sam utils/split-contigs bin/mergesam utils/temp-sink bin/fasta2fastq
	find . -name '*.o' |xargs rm -f
	find . -name  '*.core' |xargs rm -f
	find . -name '*.pyc' |xargs rm -f
