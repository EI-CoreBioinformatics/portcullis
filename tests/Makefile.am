AUTOMAKE_OPTIONS = subdir-objects
ACLOCAL_AMFLAGS = -I m4 ${ACLOCAL_FLAGS}

TESTS = check_portculis
check_PROGRAMS = check_portculis
check_portculis_CPPFLAGS = @BAMTOOLS_CPPFLAGS@ @SAMTOOLS_CPPFLAGS@ @DLIB_CPPFLAGS@ @BOOST_CPPFLAGS@ @AM_CPPFLAGS@ -I$(top_srcdir)/src -I$(top_srcdir)/src/.deps
#check_portculis_CXXFLAGS = -g3 -Wall -Werror -Wl,-rpath,@BAMTOOLS_PATH@/lib @AM_CXXFLAGS@
#check_portculis_CXXFLAGS = -g3 -Wl,-rpath,@BAMTOOLS_PATH@/lib @AM_CXXFLAGS@
check_portculis_CXXFLAGS = -g3 @AM_CXXFLAGS@
check_portculis_LDFLAGS = @BAMTOOLS_LDFLAGS@ @SAMTOOLS_LDFLAGS@ @BOOST_LDFLAGS@ -L$(top_srcdir)/src/.deps @AM_LDFLAGS@
check_portculis_LDADD = -lboost_unit_test_framework @BOOST_LIBS@ -lbam -lbamtools -lpthread -lz
check_portculis_SOURCES =   bam_utils_tests.cpp \
			    seq_utils_tests.cpp \
			    intron_tests.cpp \
			    junction_tests.cpp \
			    portculis_tests.cpp \
			    genome_mapper_tests.cpp \
			    check_portculis.cc

			    
			    


