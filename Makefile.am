
# Tells aclocal to put 'aclocal.m4' into the 'm4' directory
ACLOCAL_AMFLAGS = -I m4 ${ACLOCAL_FLAGS}

# ADDITIONAL FILES TO INSTALL
EXTRA_DIST = etc \
README.md \
COPYING \
AUTHORS \
tests/resources \
deps/dlib-18.16/dlib

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
dist_noinst_SCRIPTS = autogen.sh antigen.sh

dist_bin_SCRIPTS =  scripts/bed12.py \
		    scripts/bed12_set.py \
		    scripts/cleanup_spanki.py \
		    scripts/cufflinks2spankicov.py \
		    scripts/filter_junctions.py \
		    scripts/fs2bed.py \
		    scripts/fs2gff.py \
		    scripts/venn.py \
		    scripts/gff2bed.py \
		    scripts/star_tab2bed.py \
		    scripts/spanki2bed.py \
		    scripts/spanki_filter.py \
		    scripts/tab_bed12_set.py \
		    scripts/hpc.json \
		    scripts/read_gen.py \
		    scripts/read_gen.snakefile \
		    scripts/example_pp_ath.yaml

# 
dist_sysconf_DATA = etc/default_filter.json

# SRC DIRS
SUBDIRS = deps/htslib-1.3 src tests


