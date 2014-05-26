AUTOMAKE_OPTIONS = subdir-objects
ACLOCAL_AMFLAGS = -I m4 ${ACLOCAL_FLAGS}

# Create portculis binary
bin_PROGRAMS = portculis

#portculis_CXXFLAGS = -g3 -Wall -Werror -Wl,-rpath,@BAMTOOLS_PATH@/lib @AM_CXXFLAGS@
portculis_CXXFLAGS = -g3 -Wl,-rpath,@BAMTOOLS_PATH@/lib @AM_CXXFLAGS@
portculis_CPPFLAGS = -DCPLUSPLUS @BAMTOOLS_CPPFLAGS@ @SAMTOOLS_CPPFLAGS@ @BOOST_CPPFLAGS@ @AM_CPPFLAGS@ -I$(top_srcdir)/include -I$(top_srcdir)/src 
portculis_LDFLAGS = @BAMTOOLS_LDFLAGS@ @SAMTOOLS_LDFLAGS@ @BOOST_LDFLAGS@ @AM_LDFLAGS@
portculis_LDADD = -lbam -lbamtools -lbamtools-utils @BOOST_LIBS@ -lpthread -lz
portculis_SOURCES = portculis.cc
                
