AUTOMAKE_OPTIONS = subdir-objects
ACLOCAL_AMFLAGS = -I m4 ${ACLOCAL_FLAGS}

# Create portculis binary
bin_PROGRAMS = portculis

portculis_CXXFLAGS = -Wno-conversion-null -Wno-unused-result -Wno-unused-function -Wno-sign-compare -Wnon-virtual-dtor @AM_CXXFLAGS@
portculis_CPPFLAGS = -DSEQAN_HAS_ZLIB=1 -I$(top_srcdir)/include -I$(top_srcdir)/src @SEQAN_CPPFLAGS@ @BOOST_CPPFLAGS@ @AM_CPPFLAGS@
portculis_LDFLAGS = @BOOST_LDFLAGS@ @AM_LDFLAGS@
portculis_LDADD = @BOOST_LIBS@ -lpthread -lz @AM_LIBS@
portculis_SOURCES = \
                portculis.cc \
                portculis_args.hpp
