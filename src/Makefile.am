AUTOMAKE_OPTIONS = subdir-objects
ACLOCAL_AMFLAGS = -I m4 ${ACLOCAL_FLAGS}

# Create portculis binary
bin_PROGRAMS = portculis

portculis_CXXFLAGS = -Wno-conversion-null -Wno-unused-result -Wno-unused-function -Wno-sign-compare -Wnon-virtual-dtor -Wl,-rpath,@BAMTOOLS_PATH@/lib @AM_CXXFLAGS@
portculis_CPPFLAGS = -DCPLUSPLUS @BAMTOOLS_CPPFLAGS@ @SAMTOOLS_CPPFLAGS@ @BOOST_CPPFLAGS@ @AM_CPPFLAGS@ -I$(top_srcdir)/include -I$(top_srcdir)/src 
portculis_LDFLAGS = @BAMTOOLS_LDFLAGS@ @SAMTOOLS_LDFLAGS@ @BOOST_LDFLAGS@ @AM_LDFLAGS@
portculis_LDADD = -lbam -lbamtools @BOOST_LIBS@ -lpthread -lz
portculis_SOURCES = \
	bamtools_sort.h \
	bamtools_sort.cpp \    
	portculis.cc
                
