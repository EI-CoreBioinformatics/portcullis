AUTOMAKE_OPTIONS = subdir-objects
ACLOCAL_AMFLAGS = -I m4 ${ACLOCAL_FLAGS}

TESTS = check_portcullis
check_PROGRAMS = check_portcullis
check_portcullis_CXXFLAGS = -g3 @AM_CXXFLAGS@
check_portcullis_CPPFLAGS = -I$(top_srcdir)/deps/htslib-1.2.1 -I$(top_srcdir)/deps/samtools-1.2 -I$(top_srcdir)/include -I$(top_srcdir)/src -DPREFIX="\"$(prefix)\"" @AM_CPPFLAGS@ 
check_portcullis_LDFLAGS = -L$(top_srcdir)/deps/htslib-1.2.1 -L$(top_srcdir)/src @AM_LDFLAGS@
check_portcullis_LDADD = @AM_LIBS@ -lboost_unit_test_framework 
check_portcullis_SOURCES =  bam_utils_tests.cpp \
			    seq_utils_tests.cpp \
			    intron_tests.cpp \
			    junction_tests.cpp \
			    portcullis_tests.cpp \
			    depth_parser_tests.cpp \
			    genome_mapper_tests.cpp \
			    check_portcullis.cc
	