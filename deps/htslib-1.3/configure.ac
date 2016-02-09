# Configure script for htslib, a C library for high-throughput sequencing data.
#
#    Copyright (C) 2015 Genome Research Ltd.
#
#    Author: John Marshall <jm18@sanger.ac.uk>
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
# THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
# FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
# DEALINGS IN THE SOFTWARE.

# Autoconf initialisation
AC_PREREQ([2.68])
AC_INIT([HTSlib],[1.3],[samtools-help@lists.sourceforge.net],[],[http://www.htslib.org/])
AC_CONFIG_SRCDIR([hts.c])
AC_CONFIG_AUX_DIR([build-aux])
AC_CONFIG_MACRO_DIR([m4])
AC_PROG_CC
AC_PROG_INSTALL
AC_COPYRIGHT([Portions copyright (C) 2015 Genome Research Ltd.

This configure script is free software: you are free to change and
redistribute it.  There is NO WARRANTY, to the extent permitted by law.])

# Automake initialisation.
AM_INIT_AUTOMAKE([1.11 subdir-objects foreign no-define tar-ustar])
AM_SILENT_RULES([yes])

LT_INIT

dnl FIXME This pulls in dozens of standard header checks
AC_FUNC_MMAP
AC_CHECK_FUNCS(gmtime_r)

SAVE_LIBS=${LIBS}

AC_SEARCH_LIBS([dlopen], [dl], [DL_LIB="-ldl"])
AC_SEARCH_LIBS([sin], [m], [M_LIB="-lm"])

AX_PTHREAD
AX_CHECK_ZLIB

LIBS=${SAVE_LIBS}

AM_CXXFLAGS="-g -g3 -gdwarf-4 -Wno-deprecated-declarations"
AC_SUBST([AM_CXXFLAGS])

AM_CPPFLAGS=${ZLIB_CPPFLAGS}
AC_SUBST([AM_CPPFLAGS])

AM_LDFLAGS="${PTHREAD_CFLAGS} ${ZLIB_LDFLAGS}"
AC_SUBST([AM_LDFLAGS])

AM_LIBS="${M_LIB} ${DL_LIB} ${ZLIB_LIB} ${PTHREAD_LIBS}"
AC_SUBST([AM_LIBS])

AC_CONFIG_HEADERS([config.h])
AC_CONFIG_FILES(phtslib.pc Makefile)
AC_OUTPUT
