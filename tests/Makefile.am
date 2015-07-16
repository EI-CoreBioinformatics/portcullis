AUTOMAKE_OPTIONS = subdir-objects
ACLOCAL_AMFLAGS = -I m4 ${ACLOCAL_FLAGS}

AM_CXXFLAGS = -g -g3 -gdwarf-4 -O0 @AM_CXXFLAGS@
AM_CPPFLAGS = -I$(top_srcdir)/deps/htslib-1.2.1 -I$(top_srcdir)/deps/samtools-1.2 -I$(top_srcdir)/deps/dlib-18.16 -I$(top_srcdir)/src @AM_CPPFLAGS@ 
AM_LDFLAGS = -L$(top_srcdir)/src @AM_LDFLAGS@
LDADD = $(top_srcdir)/deps/htslib-1.2.1/libhts.a @BOOST_UNIT_TEST_FRAMEWORK_LIB@ @AM_LIBS@

TESTS_ENVIRONMENT=BOOST_TEST_REPORT_LEVEL=detailed BOOST_TEST_LOG_LEVEL=all
TESTS = check_bam check_seq check_intron check_junction
check_PROGRAMS = check_bam check_seq check_intron check_junction

check_bam_SOURCES = ../src/bam/bam_master.cc \
		    ../src/bam/bam_alignment.cc \
		    ../src/bam/bam_reader.cc \
		    ../src/bam/depth_parser.cc \
		    ../src/bam/genome_mapper.cc \
		    bam_tests.cpp
	
check_seq_SOURCES = seq_utils_tests.cpp

check_intron_SOURCES =  ../src/bam/bam_master.cc \
			../src/intron.cc \
			intron_tests.cpp
	
check_junction_SOURCES =    ../src/bam/bam_master.cc \
			    ../src/bam/bam_alignment.cc \
			    ../src/bam/bam_reader.cc \
			    ../src/bam/genome_mapper.cc \
			    ../src/intron.cc \
			    ../src/junction.cc \
			    junction_tests.cpp


	