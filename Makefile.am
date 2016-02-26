
# Tells aclocal to put 'aclocal.m4' into the 'm4' directory
ACLOCAL_AMFLAGS = -I m4 ${ACLOCAL_FLAGS}

# Blank man target (used in doc)
man: 
	cd doc && $(MAKE) $(AM_MAKEFLAGS) man
.PHONY: man
	
# ADDITIONAL FILES TO INSTALL
EXTRA_DIST = \
	README.md \
	COPYING \
	AUTHORS \
	NEWS

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
		    scripts/predict.py \
		    scripts/predict.snakefile \
		    scripts/example_pp_ath.yaml \
		    scripts/analyse_alignments.py \
		    scripts/analyse_junctions.py \
		    scripts/bed_merge.py \
		    scripts/tab.py \
		    scripts/tab_bed12_set.py \
		    scripts/tab_in_bedref.py \
		    scripts/tab_combine.py

# Config files
configdir = $(datadir)/portcullis
dist_config_DATA = data/default_filter.json

# SRC DIRS
make_dirs=deps/htslib-1.3 deps/ranger-0.3.8 lib src tests

if MAKE_DOCS
    make_dirs += doc
endif

SUBDIRS = $(make_dirs)


