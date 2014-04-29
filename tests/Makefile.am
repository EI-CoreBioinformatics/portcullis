AUTOMAKE_OPTIONS = subdir-objects
ACLOCAL_AMFLAGS = -I m4 ${ACLOCAL_FLAGS}

TESTS = check_portculis
check_PROGRAMS = check_portculis
check_portculis_SOURCES = check_portculis.cc portculis_tests.cpp genome_mapper_tests.cpp
check_portculis_CPPFLAGS = -DSEQAN_HAS_ZLIB=1 @SAMTOOLS_CPPFLAGS@ @SEQAN_CPPFLAGS@ @BOOST_CPPFLAGS@ @AM_CPPFLAGS@ -I$(top_srcdir)/include -I$(top_srcdir)/src
check_portculis_CXXFLAGS = -Wno-conversion-null -Wno-unused-result -Wno-unused-function -Wno-sign-compare -Wnon-virtual-dtor @AM_CXXFLAGS@
check_portculis_LDFLAGS = @SAMTOOLS_LDFLAGS@ @BOOST_LDFLAGS@ @AM_LDFLAGS@
check_portculis_LDADD = -lboost_unit_test_framework @BOOST_LIBS@ -lbam -lpthread -lz @AM_LIBS@

