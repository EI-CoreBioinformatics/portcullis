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
			-I$(top_srcdir)/lib/include \
			-DDATADIR=\"$(datadir)\" \
			-DEXECPREFIX=\"$(exec_prefix)\" \
			@AM_CPPFLAGS@

portcullis_LDFLAGS = \
			-static \
			@AM_LDFLAGS@

portcullis_LDADD = \
			$(top_builddir)/deps/htslib-1.3/libphts.la \
			$(top_builddir)/deps/ranger-0.3.8/libranger.la \
			$(top_builddir)/lib/libportcullis.la \
			-lboost_timer \
			-lboost_chrono \
			-lboost_filesystem \
			-lboost_program_options \
			-lboost_system \
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
