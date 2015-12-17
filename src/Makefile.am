AUTOMAKE_OPTIONS = subdir-objects

# Create portculis binary
bin_PROGRAMS = portcullis

portcullis_noinstdir = $(top_srcdir)/src

portcullis_CXXFLAGS = -g -g3 -gdwarf-4 -O3 -pedantic @AM_CXXFLAGS@

portcullis_CPPFLAGS =	-isystem $(top_srcdir)/deps/htslib-1.2.1 \
			-isystem $(top_srcdir)/deps/samtools-1.2 \
			-isystem $(top_srcdir)/deps/dlib-18.16 \
			-I$(top_srcdir)/src \
			@AM_CPPFLAGS@  

portcullis_LDFLAGS = @AM_LDFLAGS@

portcullis_LDADD = $(top_srcdir)/deps/htslib-1.2.1/libhts.a @AM_LIBS@


noinst_HEADERS =	bam/bam_master.hpp \
			bam/bam_alignment.hpp \
			bam/bam_reader.hpp \
			bam/bam_writer.hpp \
			bam/depth_parser.hpp \
			bam/genome_mapper.hpp \
			intron.hpp \
			junction.hpp \
			junction_system.hpp \
			prepare.hpp \
			junction_builder.hpp \
			junction_filter.hpp \
			bam_filter.hpp

portcullis_SOURCES =	bam/bam_master.cc \
			bam/bam_alignment.cc \
			bam/bam_reader.cc \
			bam/bam_writer.cc \
			bam/depth_parser.cc \
			bam/genome_mapper.cc \
			portcullis_fs.hpp \
			intron.cc \
			junction.cc \
			junction_system.cc \
			prepare.cc \
			junction_builder.cc \
			bam_filter.cc \
			junction_filter.cc \
			cluster.hpp \
			test.hpp \
			train.hpp \
			seq_utils.hpp \
			portcullis.cc
                
