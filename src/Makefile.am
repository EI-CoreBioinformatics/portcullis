AUTOMAKE_OPTIONS = subdir-objects

# Create portculis binary
bin_PROGRAMS = portcullis

portcullis_noinstdir = $(top_srcdir)/src

portcullis_CXXFLAGS =	-m64\
			-Wno-unused-function -Wno-unused-variable -Wno-ignored-qualifiers \
                        -Wall \
			@AM_CXXFLAGS@

portcullis_CPPFLAGS = \
			-isystem $(top_srcdir)/deps/htslib-1.3 \
			-isystem $(top_srcdir)/deps/ranger-0.3.8/include \
			-isystem $(top_srcdir)/deps/boost/build/include \
			-isystem $(top_srcdir)/lib/include \
			-DDATADIR=\"$(datadir)\" \
			-DPORTCULLIS_SITE_PKGS='"$(libdir)/python$(PYTHON_VERSION)/site-packages"' \
			@AM_CPPFLAGS@

portcullis_LDFLAGS = \
			-static \
			-L../lib \
			-L../deps/boost/build/lib \
			-L../deps/htslib-1.3 \
			-L../deps/ranger-0.3.8 \
			@AM_LDFLAGS@

portcullis_LDADD = \
			../deps/boost/build/lib/libboost_timer.a \
			../deps/boost/build/lib/libboost_chrono.a \
			../deps/boost/build/lib/libboost_filesystem.a \
			../deps/boost/build/lib/libboost_program_options.a \
			../deps/boost/build/lib/libboost_system.a \
			-lportcullis \
			-lranger \
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
			portcullis.cc
