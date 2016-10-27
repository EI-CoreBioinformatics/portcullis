
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
		    scripts/bed_merge.py \
		    scripts/bed_v_ref.py \
		    scripts/create_ref.py \
		    scripts/performance.py \
		    scripts/portcullis_convert.py \
		    scripts/portcullis_rulefilt.py \
		    scripts/split_error_classes.py \
		    scripts/tab.py \
		    scripts/tab_bed12_set.py \
		    scripts/tab_in_bedref.py \
		    scripts/tab_combine.py
		    
		    
# Config files
configdir = $(datadir)/portcullis
dist_config_DATA = \
		    data/default_filter.json \
		    data/selftrain_initial_neg.layer1.json \
		    data/selftrain_initial_neg.layer2.json \
		    data/selftrain_initial_neg.layer3.json \
		    data/selftrain_initial_neg.layer4.json \
		    data/selftrain_initial_neg.layer5.json \
		    data/selftrain_initial_neg.layer6.json \
		    data/selftrain_initial_neg.layer7.json \
		    data/selftrain_initial_pos.layer1.json \
		    data/selftrain_initial_pos.layer2.json \
		    data/selftrain_initial_pos.layer3.json

# SRC DIRS
make_dirs=deps/htslib-1.3 deps/ranger-0.3.8 lib src tests

if MAKE_DOCS
    make_dirs += doc
endif

SUBDIRS = $(make_dirs)
