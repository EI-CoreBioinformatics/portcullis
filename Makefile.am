
# Tells aclocal to put 'aclocal.m4' into the 'm4' directory
ACLOCAL_AMFLAGS = -I m4 ${ACLOCAL_FLAGS}

# ADDITIONAL FILES TO INSTALL
EXTRA_DIST = data/default_filter.json \
	data/default_filter.0.9.2.json
	README.md \
	COPYING \
	AUTHORS \
	tests/resources \
	deps/dlib-18.16/dlib


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

# Config files
configdir = $(datadir)/portcullis
dist_config_DATA = data/default_filter.json

# SRC DIRS
make_dirs=deps/htslib-1.3 lib src tests

if test @MAKE_DOCS@ = 1; then
    make_dirs +=  doc

SUBDIRS = $(make_dirs)


