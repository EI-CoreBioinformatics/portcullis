#!/bin/sh
#
# An utility script to remove all generated files.
#
# Running autogen.sh will be required after running this script since
# the 'configure' script will also be removed.
#
# This script is mainly useful when testing autoconf/automake changes
# and as a part of their development process.

# If there's a Makefile, then run the 'distclean' target first (which
# will also remove the Makefile).
if test -f Makefile; then
  make distclean
fi

# Also clean docs
cd doc
if test -f Makefile; then
    make clean
fi
cd ..

# Remove all tar-files (assuming there are some packages), and any expanded directories.
rm -f *.tar.* *.tgz
rm -Rf portcullis-*/

# Remove the build_aux directory
rm -Rf build-aux

# Also remove the autotools cache directory.
rm -Rf autom4te.cache deps/htslib-1.3/autom4te.cache deps/samtools-1.3/autom4te.cache

# Remove rest of the generated files.
rm -f Makefile.in tests/Makefile.in config.* aclocal.m4 configure depcomp install-sh missing libtool stamp-h1 ltmain.sh 

