AUTOMAKE_OPTIONS = subdir-objects
ACLOCAL_AMFLAGS = -I m4 ${ACLOCAL_FLAGS}

TESTS = check_unit_tests
check_PROGRAMS = check_unit_tests

check_unit_tests_SOURCES =	bam_tests.cpp \
				seq_utils_tests.cpp \
				intron_tests.cpp \
				junction_tests.cpp \
				check_portcullis.cc

check_unit_tests_CXXFLAGS = -g -g3 -gdwarf-4 -O0 @AM_CXXFLAGS@

check_unit_tests_CPPFLAGS =	-isystem $(top_srcdir)/deps/htslib-1.3 \
				-isystem $(top_srcdir)/deps/dlib-18.16 \
				-isystem $(top_srcdir)/lib/include \
				-DRESOURCESDIR=\"$(srcdir)/resources\" \
				-DDATADIR=\"$(datadir)\" \
				@AM_CPPFLAGS@ 

check_unit_tests_LDADD =	$(top_srcdir)/lib/libportcullis-0.11.la \
				$(top_srcdir)/deps/htslib-1.3/libhts.a \
				libgtest.la \
				@AM_LIBS@	
	
include gtest.mk
	