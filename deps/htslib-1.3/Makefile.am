AUTOMAKE_OPTIONS = subdir-objects
ACLOCAL_AMFLAGS = -I m4

pkgconfigdir = $(libdir)/pkgconfig
pkgconfig_DATA = phtslib.pc

HTSLIB_CPPFLAGS = -I$(top_srcdir) -I$(srcdir) @CPPFLAGS@
HTSLIB_LDFLAGS = -pthread @LDFLAGS@
HTSLIB_LIBS = libphts.la -lz -lm -ldl

# Create portcullis library
lib_LTLIBRARIES = libphts.la

libphts_la_LDFLAGS = -version-info 1:3:0 $(HTSLIB_LDFLAGS)
libphts_la_CPPFLAGS = $(HTSLIB_CPPFLAGS)
libphts_la_CXXFLAGS = -g -g3 -gdwarf-4 -O3 @AM_CXXFLAGS@
libphts_la_SOURCES = \
	bgzf.c \
	faidx.c \
	hfile.c \
	hfile_net.c \
	hts.c \
	htsfile.c \
	kfunc.c \
	knetfile.c \
	kstring.c \
	md5.c \
	plugin.c \
	regidx.c \
	sam.c \
	synced_bcf_reader.c \
	tbx.c \
	vcf.c \
	vcf_sweep.c \
	vcfutils.c \
	cram/cram_codecs.c \
	cram/cram_decode.c \
	cram/cram_encode.c \
	cram/cram_external.c \
	cram/cram_index.c \
	cram/cram_io.c \
	cram/cram_samtools.c \
	cram/cram_stats.c \
	cram/files.c \
	cram/mFILE.c \
	cram/open_trace_file.c \
	cram/pooled_alloc.c \
	cram/rANS_static.c \
	cram/sam_header.c \
	cram/string_alloc.c \
	cram/thread_pool.c \
	cram/vlen.c \
	cram/zfio.c \
	hfile_internal.h \
	hts_internal.h

htslib_library_includedir=$(includedir)/phtslib-@PACKAGE_VERSION@/htslib
HI = $(top_srcdir)/htslib
htslib_library_include_HEADERS = \
	$(HI)/bgzf.h \
	$(HI)/cram.h \
	$(HI)/faidx.h \
	$(HI)/hfile.h \
	$(HI)/hts_defs.h \
	$(HI)/hts.h \
	$(HI)/kbitset.h \
	$(HI)/kfunc.h \
	$(HI)/khash.h \
	$(HI)/khash_str2int.h \
	$(HI)/klist.h \
	$(HI)/knetfile.h \
	$(HI)/kseq.h \
	$(HI)/ksort.h \
	$(HI)/kstring.h \
	$(HI)/regidx.h \
	$(HI)/sam.h \
	$(HI)/synced_bcf_reader.h \
	$(HI)/tbx.h \
	$(HI)/vcf.h \
	$(HI)/vcf_sweep.h \
	$(HI)/vcfutils.h
	
cram_library_includedir=$(includedir)/phtslib-@PACKAGE_VERSION@/cram
CI = $(top_srcdir)/cram
cram_library_include_HEADERS = \
	$(CI)/cram_codecs.h \
	$(CI)/cram_decode.h \
	$(CI)/cram_encode.h \
	$(CI)/cram.h \
	$(CI)/cram_index.h \
	$(CI)/cram_io.h \
	$(CI)/cram_samtools.h \
	$(CI)/cram_stats.h \
	$(CI)/cram_structs.h \
	$(CI)/mFILE.h \
	$(CI)/misc.h \
	$(CI)/open_trace_file.h \
	$(CI)/os.h \
	$(CI)/pooled_alloc.h \
	$(CI)/rANS_byte.h \
	$(CI)/rANS_static.h \
	$(CI)/sam_header.h \
	$(CI)/string_alloc.h \
	$(CI)/thread_pool.h \
	$(CI)/vlen.h \
	$(CI)/zfio.h

TESTS = \
	test/fieldarith \
	test/sam
	
check_PROGRAMS = \
	test/fieldarith \
	test/sam
	
HTSLIB_TEST_CXXFLAGS = -g -g3 -gdwarf-4 -O0 -Wno-deprecated-declarations @AM_CXXFLAGS@
	
test_fieldarith_SOURCES = test/fieldarith.c
test_fieldarith_CXXFLAGS = $(HTSLIB_TEST_CXXFLAGS)
test_fieldarith_CPPFLAGS = -I.. $(HTSLIB_CPPFLAGS)
test_fieldarith_LDFLAGS = $(HTSLIB_LDFLAGS)
test_fieldarith_LDADD = $(HTSLIB_LIBS)

test_sam_SOURCES = test/sam.c
test_sam_CXXFLAGS = $(HTSLIB_TEST_CXXFLAGS)
test_sam_CPPFLAGS = $(HTSLIB_CPPFLAGS)
test_sam_LDFLAGS = $(HTSLIB_LDFLAGS)
test_sam_LDADD = $(HTSLIB_LIBS)



clean-local: clean-local-check
.PHONY: clean-local-check
clean-local-check:
	-rm -rf test/*.tmp.* test/*.tmp
