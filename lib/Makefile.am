AUTOMAKE_OPTIONS = subdir-objects

pkgconfigdir = $(libdir)/pkgconfig
pkgconfig_DATA = portcullis-0.11.pc

# Create portcullis library
lib_LTLIBRARIES = libportcullis-0.11.la

# Library inputs
libportcullis_0_11_la_LDFLAGS = -version-info 0:11:0 -L$(top_srcdir)/deps/htslib-1.3
libportcullis_0_11_la_SOURCES =	\
	src/bam_master.cc \
	src/bam_alignment.cc \
	src/bam_reader.cc \
	src/bam_writer.cc \
	src/depth_parser.cc \
	src/genome_mapper.cc \
	src/intron.cc \
	src/junction.cc \
	src/junction_system.cc

library_includedir=$(includedir)/portcullis-@PACKAGE_VERSION@/portcullis
PI = include/portcullis
library_include_HEADERS = \
	$(PI)/bam/bam_master.hpp \
	$(PI)/bam/bam_alignment.hpp \
	$(PI)/bam/bam_reader.hpp \
	$(PI)/bam/bam_writer.hpp \
	$(PI)/bam/depth_parser.hpp \
	$(PI)/bam/genome_mapper.hpp \
	$(PI)/intron.hpp \
	$(PI)/junction.hpp \
	$(PI)/junction_system.hpp \
	$(PI)/portcullis_fs.hpp \
	$(PI)/seq_utils.hpp

libportcullis_0_11_la_CPPFLAGS = \
	-isystem $(top_srcdir)/deps/htslib-1.3 \
	-isystem $(top_srcdir)/lib/include \
	-DDATADIR=\"$(datadir)\" \
	@AM_CPPFLAGS@

libportcullis_0_11_la_CXXFLAGS = -g -g3 -gdwarf-4 -O3 -pedantic @AM_CXXFLAGS@
