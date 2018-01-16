AUTOMAKE_OPTIONS = subdir-objects

# Create portculis binary
bin_PROGRAMS = portcullis

portcullis_noinstdir = $(top_srcdir)/src

portcullis_CXXFLAGS =	\
			-O3 -pedantic -ansi -Wall -Wextra \
			-Wno-unused-function -Wno-unused-variable \
			-std=c++11 \
			@AM_CXXFLAGS@

portcullis_CPPFLAGS = \
			-isystem $(top_srcdir)/deps/htslib-1.3 \
			-isystem $(top_srcdir)/deps/ranger-0.3.8/include \
			-isystem $(top_srcdir)/deps/boost \
      -isystem $(top_srcdir)/lib/include \
			-DDATADIR=\"$(datadir)\" \
			-DPORTCULLIS_SITE_PKGS='"$(libdir)/python$(PYTHON_VERSION)/site-packages"' \
			@AM_CPPFLAGS@

portcullis_LDFLAGS = \
			-static \
			-L../lib \
			-L../deps/htslib-1.3 \
			-L../deps/ranger-0.3.8 \
			@AM_LDFLAGS@

portcullis_LDADD = \
			-lportcullis \
			-lranger \
			-lphts \
			../deps/boost/stage/lib/libboost_timer.a \
      ../deps/boost/stage/lib/libboost_chrono.a \
      ../deps/boost/stage/lib/libboost_filesystem.a \
      ../deps/boost/stage/lib/libboost_program_options.a \
      ../deps/boost/stage/lib/libboost_system.a \
			@AM_LIBS@


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
