AUTOMAKE_OPTIONS = subdir-objects

pkgconfigdir = $(libdir)/pkgconfig
pkgconfig_DATA = portcullis.pc

# Create portcullis library
lib_LTLIBRARIES = libportcullis.la

# Library inputs
libportcullis_la_LDFLAGS = \
	-version-info 19:0:19 \
	-static \
	@AM_LDFLAGS@

EXTRA_libportcullis_la_DEPENDENCIES = \
	$(top_builddir)/deps/boost/build/lib/libboost_timer.a \
	$(top_builddir)/deps/boost/build/lib/libboost_chrono.a \
	$(top_builddir)/deps/boost/build/lib/libboost_filesystem.a \
	$(top_builddir)/deps/boost/build/lib/libboost_system.a
	

libportcullis_la_LIBADD = \
	$(top_builddir)/deps/htslib-1.3/libphts.la \
	$(top_builddir)/deps/ranger-0.3.8/libranger.la \
	-lz \
	@AM_LIBS@

libportcullis_la_SOURCES =	\
	src/bam_master.cc \
	src/bam_alignment.cc \
	src/bam_reader.cc \
	src/bam_writer.cc \
	src/depth_parser.cc \
	src/genome_mapper.cc \
	src/markov_model.cc \
	src/model_features.cc \
	src/intron.cc \
	src/junction.cc \
	src/junction_system.cc \
    src/performance.cc \
    src/knn.cc \
	src/enn.cc \
	src/smote.cc

library_includedir=$(includedir)/portcullis-@PACKAGE_VERSION@/portcullis
PI = include/portcullis
library_include_HEADERS = \
	$(PI)/bam/bam_master.hpp \
	$(PI)/bam/bam_alignment.hpp \
	$(PI)/bam/bam_reader.hpp \
	$(PI)/bam/bam_writer.hpp \
	$(PI)/bam/depth_parser.hpp \
	$(PI)/bam/genome_mapper.hpp \
	$(PI)/ml/markov_model.hpp \
	$(PI)/ml/model_features.hpp \
	$(PI)/ml/performance.hpp \
	$(PI)/ml/k_fold.hpp \
	$(PI)/ml/knn.hpp \
	$(PI)/ml/enn.hpp \
	$(PI)/ml/smote.hpp \
	$(PI)/kmer.hpp \
	$(PI)/python_helper.hpp \
	$(PI)/intron.hpp \
	$(PI)/junction.hpp \
	$(PI)/junction_system.hpp \
	$(PI)/portcullis_fs.hpp \
	$(PI)/seq_utils.hpp


libportcullis_la_CPPFLAGS = \
	-I$(top_srcdir)/deps/htslib-1.3 \
	-I$(top_srcdir)/deps/ranger-0.3.8/include \
	-I$(top_srcdir)/deps/boost/build/include \
  	-I$(top_srcdir)/lib/include \
	-DDATADIR=\"$(datadir)\" \
	@AM_CPPFLAGS@

libportcullis_la_CXXFLAGS = \
	-m64 \
	-Wno-unused-function -Wno-unused-variable -Wno-ignored-qualifiers \
	-Wall \
	-fPIC \
	@AM_CXXFLAGS@
