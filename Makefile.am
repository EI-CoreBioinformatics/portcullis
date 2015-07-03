
# Tells aclocal to put 'aclocal.m4' into the 'm4' directory
ACLOCAL_AMFLAGS = -I m4 ${ACLOCAL_FLAGS}

# ADDITIONAL FILES TO INSTALL
EXTRA_DIST = etc \
README.md \
COPYING \
AUTHORS \
tests/resources \
deps/dlib-18.16/dlib \
deps/htslib-1.2.1/*.c \
deps/htslib-1.2.1/*.h \
deps/htslib-1.2.1/*.mk \
deps/htslib-1.2.1/htslib.pc.in \
deps/htslib-1.2.1/*.1 \
deps/htslib-1.2.1/*.5 \
deps/htslib-1.2.1/cram/*.c \
deps/htslib-1.2.1/cram/*.h \
deps/htslib-1.2.1/htslib/*.h \
deps/htslib-1.2.1/test/*.sam \
deps/htslib-1.2.1/test/*.fa \
deps/htslib-1.2.1/test/*.fai \
deps/htslib-1.2.1/test/*.pl \
deps/htslib-1.2.1/test/*.c \
deps/htslib-1.2.1/README \
deps/htslib-1.2.1/LICENSE \
deps/htslib-1.2.1/INSTALL \
deps/htslib-1.2.1/NEWS \
deps/samtools-1.2/*.c \
deps/samtools-1.2/*.h \
deps/samtools-1.2/*.1 \
deps/samtools-1.2/*/*.c \
deps/samtools-1.2/*/*.h \
deps/samtools-1.2/*/*.pl \
deps/samtools-1.2/*/*.py \
deps/samtools-1.2/*/*.sam \
deps/samtools-1.2/test/*/*.c \
deps/samtools-1.2/test/*/*.sam \
deps/samtools-1.2/test/*/*.bam \
deps/samtools-1.2/test/*/*.bai \
deps/samtools-1.2/test/*/*.cram \
deps/samtools-1.2/test/*/*.crai \
deps/samtools-1.2/test/*/*.fa \
deps/samtools-1.2/test/*/*.fai \
deps/samtools-1.2/test/*/*.expected \
deps/samtools-1.2/test/*/*.expected.err \
deps/samtools-1.2/misc/*.lua \
deps/samtools-1.2/misc/HmmGlocal.java \
deps/samtools-1.2/misc/plot-bamstats \
deps/samtools-1.2/test/dat/cram_md5/* \
deps/samtools-1.2/test/mpileup/md5/* \
deps/samtools-1.2/test/mpileup/expected/* \
deps/samtools-1.2/README \
deps/samtools-1.2/LICENSE \
deps/samtools-1.2/INSTALL \
deps/samtools-1.2/AUTHORS \
deps/samtools-1.2/NEWS

#dist-hook:
#	rm -rf `find $(distdir)/deps/dlib-18.16 -type d -name .git`
#	rm -rf $(distdir)/deps/dlib-18.16/.gitignore
#	rm -rf `find $(distdir)/deps/htslib-1.2.1 -type d -name .git`
#	rm -rf $(distdir)/deps/htslib-1.2.1/.gitignore
#	rm -rf `find $(distdir)/deps/samtools-1.2 -type d -name .git`
#	rm -rf $(distdir)/deps/samtools-1.2/.gitignore
#	rm -rf `find $(distdir)/tests/resources -type d -name .git`
#	rm -rf $(distdir)/tests/resources/.gitignore


# FILES NOT TO INSTALL
dist_noinst_SCRIPTS = autogen.sh clean.sh

dist_bin_SCRIPTS =  scripts/filter_junctions.py \
		    scripts/fs2bed.py \
		    scripts/fs2gff.py \
		    scripts/venn.py \
		    scripts/gff2bed.py

# 
dist_sysconf_DATA = etc/default_filter.json

# SRC DIRS
SUBDIRS = deps/htslib-1.2.1 deps/samtools-1.2 src tests


