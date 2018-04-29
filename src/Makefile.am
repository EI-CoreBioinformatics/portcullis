AUTOMAKE_OPTIONS = subdir-objects

# Create portculis binary
bin_PROGRAMS = portcullis

portcullis_noinstdir = $(top_srcdir)/src

portcullis_CXXFLAGS = \
			-m64 \
			-Wno-unused-function -Wno-unused-variable -Wno-ignored-qualifiers \
			-Wall \
			@AM_CXXFLAGS@

portcullis_CPPFLAGS = \
			-I$(top_srcdir)/deps/htslib-1.3 \
			-I$(top_srcdir)/deps/ranger-0.3.8/include \
			-I$(top_srcdir)/deps/boost/build/include \
			-I$(top_srcdir)/lib/include \
			-DDATADIR=\"$(datadir)\" \
			@AM_CPPFLAGS@

portcullis_LDFLAGS = \
			-static \
			@AM_LDFLAGS@

portcullis_LDADD = \
			$(top_builddir)/deps/boost/build/lib/libboost_timer.a \
			$(top_builddir)/deps/boost/build/lib/libboost_chrono.a \
			$(top_builddir)/deps/boost/build/lib/libboost_filesystem.a \
			$(top_builddir)/deps/boost/build/lib/libboost_program_options.a \
			$(top_builddir)/deps/boost/build/lib/libboost_system.a \
			$(top_builddir)/deps/htslib-1.3/libphts.la \
			$(top_builddir)/deps/ranger-0.3.8/libranger.la \
			$(top_builddir)/lib/libportcullis.la \
			@AM_LIBS@


noinst_HEADERS = \
			prepare.hpp \
			junction_builder.hpp \
			junction_filter.hpp \
			bam_filter.hpp

portcullis_SOURCES = \
			prepare.cc \
			junction_builder.cc \
			bam_filter.cc \
			junction_filter.cc \
			portcullis.cc
