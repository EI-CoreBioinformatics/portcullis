AUTOMAKE_OPTIONS = subdir-objects
ACLOCAL_AMFLAGS = -I m4 ${ACLOCAL_FLAGS}

EXTRA_DIST = \
	resources/artha_chr4.fa \
	resources/bam1.bam \
	resources/bam2.bam \
	resources/clipped3.bam \
	resources/clipped3.bam.bai \
	resources/ecoli.bam.depth \
	resources/ecoli.fa \
	resources/sorted.bam \
	resources/sorted.bam.bai \
	resources/unsorted.bam
	
TEST_EXTENSIONS = .sh
SH_LOG_COMPILER = $(SHELL)
AM_SH_LOG_FLAGS =

TESTS = check_unit_tests test_substeps.sh test_full.sh
check_PROGRAMS = check_unit_tests

noinst_HEADERS =    gtest/gtest.h \
		    gtest/src/gtest-all.cc \
		    gtest/src/gtest_main.cc

check_unit_tests_SOURCES =	bam_tests.cpp \
				seq_utils_tests.cpp \
				kmer_tests.cpp \
				smote_tests.cpp \
				intron_tests.cpp \
				junction_tests.cpp \
				check_portcullis.cc

check_unit_tests_CXXFLAGS = -O0 @AM_CXXFLAGS@ @CXXFLAGS@

check_unit_tests_CPPFLAGS =	-isystem $(top_srcdir)/deps/htslib-1.3 \
				-isystem $(top_srcdir)/deps/ranger-0.3.8/include \
				-I$(top_srcdir)/lib/include \
				-DRESOURCESDIR=\"$(top_srcdir)/tests/resources\" \
				-DDATADIR=\"$(datadir)\" \
				@AM_CPPFLAGS@
				@CPPFLAGS@

check_unit_tests_LDFLAGS =	-L../lib \
				-L../deps/htslib-1.3 \
				-L../deps/ranger-0.3.8 \
				@AM_LDFLAGS@ \
				@LDFLAGS@
	
check_unit_tests_LDADD =	libgtest.la \
				-lz \
				-lportcullis \
				-lphts \
				-lranger \
				@AM_LIBS@	

clean-local: clean-local-check
.PHONY: clean-local-check
clean-local-check:
	-rm -rf temp
	
include gtest.mk
	
