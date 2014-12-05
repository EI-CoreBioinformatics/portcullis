AUTOMAKE_OPTIONS = subdir-objects

# Create portculis binary
bin_PROGRAMS = portcullis

portcullis_CXXFLAGS = -g3 -O2 @AM_CXXFLAGS@
portcullis_CPPFLAGS = @AM_CPPFLAGS@ -I$(top_srcdir)/src -I$(top_srcdir)/include 
portcullis_LDFLAGS = @AM_LDFLAGS@
portcullis_LDADD = @AM_LIBS@
portcullis_SOURCES =	bam_utils.hpp \
			depth_parser.hpp \
			genome_mapper.hpp \
			intron.hpp \
			junction.hpp \
			junction_builder.hpp \
			junction_system.hpp \
			prepare.hpp \
			filter.hpp \
			cluster.hpp \
			seq_utils.hpp \
			portcullis.cc
                
