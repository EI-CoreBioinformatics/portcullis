AUTOMAKE_OPTIONS = subdir-objects

# Create portculis binary
bin_PROGRAMS = portcullis

portcullis_noinstdir = $(top_srcdir)/src

portcullis_CXXFLAGS = -g -g3 -gdwarf-4 -O3 -pedantic @AM_CXXFLAGS@

portcullis_CPPFLAGS =	-isystem $(top_srcdir)/deps/htslib-1.3 \
			-isystem $(top_srcdir)/lib/include \
			-I$(top_srcdir)/src \
			-DDATADIR=\"$(datadir)\" \
			@AM_CPPFLAGS@

portcullis_LDFLAGS = -all-static @AM_LDFLAGS@

portcullis_LDADD =	../lib/libportcullis-0.11.la \
			../deps/htslib-1.3/libphts.la \
			-ldl \
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
                
