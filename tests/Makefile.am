AUTOMAKE_OPTIONS = subdir-objects
ACLOCAL_AMFLAGS = -I m4 ${ACLOCAL_FLAGS}

TESTS_ENVIRONMENT=BOOST_TEST_REPORT_LEVEL=detailed BOOST_TEST_LOG_LEVEL=all
TESTS = check_portcullis
check_PROGRAMS = check_portcullis
check_portcullis_CXXFLAGS = -g3 -O0 @AM_CXXFLAGS@
check_portcullis_CPPFLAGS = -I$(top_srcdir)/deps/htslib-1.2.1 -I$(top_srcdir)/deps/samtools-1.2 -I$(top_srcdir)/deps/dlib-18.16 -I$(top_srcdir)/src @AM_CPPFLAGS@ 
check_portcullis_LDFLAGS = -L$(top_srcdir)/src @AM_LDFLAGS@
check_portcullis_LDADD = $(top_srcdir)/deps/htslib-1.2.1/libhts.a @AM_LIBS@ -lboost_unit_test_framework 
check_portcullis_SOURCES =  ../src/htslib_helper.cc \
			    ../src/samtools_helper.cc
			    bam_utils_tests.cpp \
			    bam_filter_tests.cpp \
			    seq_utils_tests.cpp \
			    intron_tests.cpp \
			    junction_tests.cpp \
			    portcullis_tests.cpp \
			    depth_parser_tests.cpp \
			    genome_mapper_tests.cpp \
			    check_portcullis.cc
	