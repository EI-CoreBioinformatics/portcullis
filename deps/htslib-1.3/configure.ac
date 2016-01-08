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

dnl Process this file with autoconf to produce a configure script
AC_INIT([HTSlib],[1.3],[samtools-help@lists.sourceforge.net],[],[http://www.htslib.org/])
AC_PREREQ([2.69])  dnl This version introduced 4-argument AC_CHECK_HEADER
AC_CONFIG_SRCDIR([hts.c])
AC_CONFIG_AUX_DIR([build-aux])
AC_CONFIG_MACRO_DIR([m4])

# Automake initialisation.
AM_INIT_AUTOMAKE([1.11 subdir-objects foreign no-define tar-ustar])
AM_SILENT_RULES([yes])

dnl Copyright notice to be copied into the generated configure script
AC_COPYRIGHT([Portions copyright (C) 2015 Genome Research Ltd.

This configure script is free software: you are free to change and
redistribute it.  There is NO WARRANTY, to the extent permitted by law.])

AC_PROG_CC
AC_PROG_INSTALL
LT_INIT

dnl FIXME This pulls in dozens of standard header checks
AC_FUNC_MMAP
AC_CHECK_FUNCS(gmtime_r)

AX_CHECK_ZLIB

AC_CONFIG_HEADERS([config.h])
AC_CONFIG_FILES(phtslib.pc Makefile)
AC_OUTPUT
