AUTOMAKE_OPTIONS = subdir-objects
ACLOCAL_AMFLAGS = -I m4 ${ACLOCAL_FLAGS}

AM_CXXFLAGS = -g -g3 -gdwarf-4 -O0 @AM_CXXFLAGS@
 
AM_LDFLAGS = -L$(top_srcdir)/src @AM_LDFLAGS@

TESTS = check_unit_tests
check_PROGRAMS = check_unit_tests

check_unit_tests_SOURCES = ../src/bam/bam_master.cc \
		    ../src/bam/bam_alignment.cc \
		    ../src/bam/bam_reader.cc \
		    ../src/bam/genome_mapper.cc \
		    ../src/bam/depth_parser.cc \
		    ../src/intron.cc \
		    ../src/junction.cc \
		    bam_tests.cpp \
		    bam_filter_tests.cpp \
		    seq_utils_tests.cpp \
		    intron_tests.cpp \
		    junction_tests.cpp \
		    check_portcullis.cc

check_unit_tests_CPPFLAGS = -isystem $(top_srcdir)/deps/htslib-1.2.1 \
		    -isystem $(top_srcdir)/deps/samtools-1.2 \
		    -isystem $(top_srcdir)/deps/dlib-18.16 \
		    -I$(top_srcdir)/src \
		    @AM_CPPFLAGS@ 

check_unit_tests_LDADD = $(top_srcdir)/deps/htslib-1.2.1/libhts.a \
		    libgtest.la \
		    @AM_LIBS@	
	
include gtest.mk
	