AUTOMAKE_OPTIONS = subdir-objects

# Create portculis binary
bin_PROGRAMS = portcullis

portcullis_noinstdir = $(top_srcdir)/src

portcullis_CXXFLAGS =	-O3 -pedantic -ansi -Wall -Wextra \
			-Wno-unused-function -Wno-unused-const-variable \
			-std=c++11 \
			@AM_CXXFLAGS@

portcullis_CPPFLAGS =	-isystem $(top_srcdir)/deps/htslib-1.3 \
			-isystem $(top_srcdir)/deps/ranger-0.3.8/include \
			-isystem $(top_srcdir)/lib/include \
			-I$(top_srcdir)/src \
			-DDATADIR=\"$(datadir)\" \
			@AM_CPPFLAGS@

portcullis_LDFLAGS =	-static \
			-L../lib \
			-L../deps/htslib-1.3 \
			-L../deps/ranger-0.3.8 \
			@AM_LDFLAGS@

portcullis_LDADD =	@AM_LIBS@ \
			-lportcullis \
			-lranger \
			-lphts \
			-lz
						

noinst_HEADERS =	prepare.hpp \
			junction_builder.hpp \
			junction_filter.hpp \
			bam_filter.hpp \
			train.hpp 

portcullis_SOURCES =	prepare.cc \
			junction_builder.cc \
			bam_filter.cc \
			junction_filter.cc \
			train.cc \
			portcullis.cc
