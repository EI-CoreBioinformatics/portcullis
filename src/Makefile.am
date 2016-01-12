AUTOMAKE_OPTIONS = subdir-objects

# Create portculis binary
bin_PROGRAMS = portcullis

portcullis_noinstdir = $(top_srcdir)/src

portcullis_CXXFLAGS = -O3 -pedantic @AM_CXXFLAGS@ @CXXFLAGS@

portcullis_CPPFLAGS =	-isystem $(top_srcdir)/deps/htslib-1.3 \
			-isystem $(top_srcdir)/lib/include \
			-I$(top_srcdir)/src \
			-DDATADIR=\"$(datadir)\" \
			@AM_CPPFLAGS@ \
			@CPPFLAGS@

portcullis_LDFLAGS =	-all-static \
			-L../lib \
			-L../deps/htslib-1.3 \
			@AM_LDFLAGS@ \
			@LDFLAGS@

portcullis_LDADD =	-lportcullis \
			-lphts \
			@AM_LIBS@


noinst_HEADERS =	prepare.hpp \
			junction_builder.hpp \
			junction_filter.hpp \
			bam_filter.hpp

portcullis_SOURCES =	prepare.cc \
			junction_builder.cc \
			bam_filter.cc \
			junction_filter.cc \
			cluster.hpp \
			test.hpp \
			train.hpp \
			portcullis.cc
                
