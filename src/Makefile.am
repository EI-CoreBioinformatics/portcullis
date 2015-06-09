AUTOMAKE_OPTIONS = subdir-objects

# Create portculis binary
bin_PROGRAMS = portcullis

portcullisdir = $(top_srcdir)/src

portcullis_CXXFLAGS = -g3 -O3 @AM_CXXFLAGS@

portcullis_CPPFLAGS =	-I$(top_srcdir)/deps/htslib-1.2.1 \
			-I$(top_srcdir)/deps/samtools-1.2 \
			-I$(top_srcdir)/deps/dlib-18.16 \
			-I$(top_srcdir)/src \
			@AM_CPPFLAGS@  

portcullis_LDFLAGS = @AM_LDFLAGS@

portcullis_LDADD = $(top_srcdir)/deps/htslib-1.2.1/libhts.a @AM_LIBS@

portcullis_HEADERS =	samtools_helper.hpp \
			junction.hpp \
			prepare.hpp \
			junction_builder.hpp \
			bam_filter.hpp

portcullis_SOURCES =	portcullis_fs.hpp \
			samtools_helper.cc \
			intron.hpp \
			junction.cc \
			junction_system.hpp \
			prepare.cc \
			junction_builder.cc \
			bam_filter.cc \
			junction_filter.hpp \
			cluster.hpp \
			test.hpp \
			train.hpp \
			seq_utils.hpp \
			portcullis.cc
                
