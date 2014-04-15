AUTOMAKE_OPTIONS = subdir-objects
ACLOCAL_AMFLAGS = -I m4 ${ACLOCAL_FLAGS}

# Create portculis binary
bin_PROGRAMS = portculis

portculis_CXXFLAGS = -Wno-conversion-null -Wno-unused-result -Wno-unused-function -Wno-sign-compare -Wnon-virtual-dtor @AM_CXXFLAGS@
portculis_CPPFLAGS = -DSEQAN_HAS_ZLIB=1 -I$(top_srcdir)/src/inc -I$(top_srcdir)/src @BOOST_CPPFLAGS@ @AM_CPPFLAGS@
portculis_LDFLAGS = @BOOST_SYSTEM_LDFLAGS@ @BOOST_FILESYSTEM_LDFLAGS@
portculis_LDADD = @BOOST_SYSTEM_LIBS@ @BOOST_FILESYSTEM_LIBS@ -lpthread -lz @AM_LIBS@
portculis_SOURCES = \
                portculis.cc \
                portculis_args.hpp
