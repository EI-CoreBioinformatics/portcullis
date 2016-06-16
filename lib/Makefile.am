AUTOMAKE_OPTIONS = subdir-objects

pkgconfigdir = $(libdir)/pkgconfig
pkgconfig_DATA = portcullis.pc

# Create portcullis library
lib_LTLIBRARIES = libportcullis.la

# Library inputs
libportcullis_la_LDFLAGS = -version-info 11:0:11
libportcullis_la_SOURCES =	\
	src/bam_master.cc \
	src/bam_alignment.cc \
	src/bam_reader.cc \
	src/bam_writer.cc \
	src/depth_parser.cc \
	src/genome_mapper.cc \
	src/markov_model.cc \
	src/model_features.cc \
	src/python_exe.cc \
	src/intron.cc \
	src/junction.cc \
	src/junction_system.cc \
	src/rule_parser.cc \
	src/performance.cc \
	src/ss_forest.cc \
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
	$(PI)/ml/ss_forest.hpp \
	$(PI)/kmer.hpp \
	$(PI)/python_exe.hpp \
	$(PI)/intron.hpp \
	$(PI)/junction.hpp \
	$(PI)/junction_system.hpp \
	$(PI)/portcullis_fs.hpp \
	$(PI)/seq_utils.hpp \
	$(PI)/rule_parser.hpp
	

libportcullis_la_CPPFLAGS = \
	-isystem $(top_srcdir)/deps/htslib-1.3 \
	-isystem $(top_srcdir)/deps/ranger-0.3.8/include \
	-isystem $(top_srcdir)/lib/include \
	-DDATADIR=\"$(datadir)\" \
	@AM_CPPFLAGS@ \
	@CPPFLAGS@

libportcullis_la_CXXFLAGS = -g -g3 -gdwarf-4 -O3 -pedantic @AM_CXXFLAGS@
