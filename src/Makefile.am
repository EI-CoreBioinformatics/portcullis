AUTOMAKE_OPTIONS = subdir-objects

# Create portculis binary
bin_PROGRAMS = portculis

portculis_CXXFLAGS = -g3 -O2
portculis_CPPFLAGS = @AM_CPPFLAGS@ -I$(top_srcdir)/src 
portculis_LDFLAGS = @AM_LDFLAGS@
portculis_LDADD = @AM_LIBS@
portculis_SOURCES =	bam_utils.hpp \
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
			portculis.cc
                
