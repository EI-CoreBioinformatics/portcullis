AUTOMAKE_OPTIONS = subdir-objects
ACLOCAL_AMFLAGS = -I m4 ${ACLOCAL_FLAGS}

# additional include paths necessary to compile the C++ library
AM_CPPFLAGS = -DSEQAN_HAS_ZLIB=1 -I$(top_srcdir)/src/inc -I$(top_srcdir)/src @SEQAN_CPPFLAGS@ @BOOST_CPPFLAGS@ @AM_CPPFLAGS@
AM_CXXFLAGS = -Wno-conversion-null -Wno-unused-result -Wno-unused-function -Wno-sign-compare -Wnon-virtual-dtor @AM_CXXFLAGS@
AM_LDFLAGS = @BOOST_SYSTEM_LDFLAGS@ @BOOST_FILESYSTEM_LDFLAGS@ @BOOST_TEST_LDFLAGS@
LIBS = @BOOST_SYSTEM_LIBS@ @BOOST_FILESYSTEM_LIBS@ @BOOST_TEST_LIBS@ -lpthread -lz @AM_LIBS@

TESTS = #check_comp
#check_PROGRAMS = #check_comp

#check_comp_SOURCES = #check_comp.cc

