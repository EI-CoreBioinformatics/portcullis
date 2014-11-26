AUTOMAKE_OPTIONS = subdir-objects
ACLOCAL_AMFLAGS = -I m4 ${ACLOCAL_FLAGS}

TESTS = check_portcullis
check_PROGRAMS = check_portcullis
check_portcullis_CXXFLAGS = -g3 @AM_CXXFLAGS@
check_portcullis_CPPFLAGS = @AM_CPPFLAGS@ -I$(top_srcdir)/src -I$(top_srcdir)/src/.deps
check_portcullis_LDFLAGS = @AM_LDFLAGS@ -L$(top_srcdir)/src/.deps @AM_LDFLAGS@
check_portcullis_LDADD = @AM_LIBS@ -lboost_unit_test_framework 
check_portcullis_SOURCES =   bam_utils_tests.cpp \
			    seq_utils_tests.cpp \
			    intron_tests.cpp \
			    junction_tests.cpp \
			    portculis_tests.cpp \
			    depth_parser_tests.cpp \
			    genome_mapper_tests.cpp \
			    check_portcullis.cc
	