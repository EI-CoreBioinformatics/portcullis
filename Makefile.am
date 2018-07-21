
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
dist_noinst_SCRIPTS = autogen.sh antigen.sh build_boost.sh
		    
# Config files
configdir = $(datadir)/portcullis
dist_config_DATA = \
		    data/default_filter.json \
		    data/low_juncs_filter.json
balanceddir = $(datadir)/portcullis/balanced
dist_balanced_DATA = \
		    data/balanced/selftrain_initial_neg.layer1.json \
		    data/balanced/selftrain_initial_neg.layer2.json \
		    data/balanced/selftrain_initial_neg.layer3.json \
		    data/balanced/selftrain_initial_neg.layer4.json \
		    data/balanced/selftrain_initial_neg.layer5.json \
		    data/balanced/selftrain_initial_neg.layer6.json \
		    data/balanced/selftrain_initial_neg.layer7.json \
		    data/balanced/selftrain_initial_pos.layer1.json \
		    data/balanced/selftrain_initial_pos.layer2.json \
		    data/balanced/selftrain_initial_pos.layer3.json

precisedir = $(datadir)/portcullis/precise
dist_precise_DATA = \
		    data/precise/selftrain_initial_neg.layer1.json \
		    data/precise/selftrain_initial_neg.layer2.json \
		    data/precise/selftrain_initial_neg.layer3.json \
		    data/precise/selftrain_initial_neg.layer4.json \
		    data/precise/selftrain_initial_neg.layer5.json \
		    data/precise/selftrain_initial_neg.layer6.json \
		    data/precise/selftrain_initial_neg.layer7.json \
		    data/precise/selftrain_initial_pos.layer1.json \
		    data/precise/selftrain_initial_pos.layer2.json \
		    data/precise/selftrain_initial_pos.layer3.json

# SRC DIRS
make_dirs=deps/htslib-1.3 deps/ranger-0.3.8 lib src tests scripts

if MAKE_DOCS
    make_dirs += doc
endif

SUBDIRS = $(make_dirs)
