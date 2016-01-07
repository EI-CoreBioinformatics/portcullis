AUTOMAKE_OPTIONS = subdir-objects

pkgconfigdir = $(libdir)/pkgconfig
pkgconfig_DATA = htslib.pc

# Create portcullis library
lib_LTLIBRARIES = htslib-@PACKAGE_VERSION@.la

BUILT_TEST_PROGRAMS = \
	test/fieldarith \
	test/hfile \
	test/sam \
	test/test-regidx \
	test/test_view \
	test/test-vcf-api \
	test/test-vcf-sweep