
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

dist_bin_SCRIPTS =  scripts/analyse_alignments.py \
		    scripts/analyse_junctions.py \
		    scripts/bed12.py \
		    scripts/bed12_set.py \
		    scripts/bed_merge.py \
		    scripts/bed_v_ref.py \
		    scripts/cleanup_spanki.py \
		    scripts/create_ref.py \
		    scripts/cufflinks2spankicov.py \
		    scripts/example_pp_ath.yaml \
		    scripts/filter_junctions.py \
		    scripts/fs2bed.py \
		    scripts/fs2gff.py \
		    scripts/gff2bed.py \
		    scripts/gtf2bed.py \
		    scripts/hpc.json \
		    scripts/ml.py \
		    scripts/pca.py \
		    scripts/performance.py \
		    scripts/predict.snakefile \
		    scripts/read_gen.snakefile \
		    scripts/snakey.py \
		    scripts/spanki2bed.py \
		    scripts/spanki_filter.py \
		    scripts/split_error_classes.py \
		    scripts/star_tab2bed.py \
		    scripts/tab.py \
		    scripts/tab_bed12_set.py \
		    scripts/tab_in_bedref.py \
		    scripts/tab_combine.py \
		    scripts/tophat2bed12.py \
		    scripts/truesight2bed.py \
		    scripts/venn.py
		    
# Config files
configdir = $(datadir)/portcullis
dist_config_DATA = \
		    data/default_filter.json \
		    data/default_model.forest

# SRC DIRS
make_dirs=deps/htslib-1.3 deps/ranger-0.3.8 lib src tests

if MAKE_DOCS
    make_dirs += doc
endif

SUBDIRS = $(make_dirs)
