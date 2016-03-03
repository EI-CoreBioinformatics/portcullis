##########################################################################
# Setup
##########################################################################

# Autoconf setup
AC_PREREQ([2.68])
AC_INIT([ranger],[0.3.8],[wright@imbs.uni-luebeck.de],[ranger],[])
AC_CONFIG_SRCDIR([src/Tree.cpp])
AC_CONFIG_AUX_DIR([build-aux])
AC_CONFIG_MACRO_DIR([m4])
AC_CONFIG_HEADERS([config.h])
AC_CONFIG_FILES([ranger.pc Makefile])
AC_LANG([C++])
AC_PROG_INSTALL
AC_PROG_CXX
AC_PROG_CC

# Ensure there's C++11 support
AX_CXX_COMPILE_STDCXX_11(,[mandatory])

# Automake setup
AM_INIT_AUTOMAKE([1.11 subdir-objects foreign -Wall no-define tar-ustar])
m4_pattern_allow([AM_PROG_AR])
AM_PROG_AR
AM_SILENT_RULES([yes])

# Libtool setup
LT_INIT

######################################################################
# Dependency checks
######################################################################

# Check some (almost) standard functionality is present that we require
#AC_CHECK_HEADER_STDBOOL
AC_FUNC_ERROR_AT_LINE
AC_TYPE_INT16_T
AC_TYPE_INT32_T
AC_TYPE_INT64_T
AC_TYPE_UINT8_T
AC_TYPE_UINT16_T
AC_TYPE_UINT32_T
AC_TYPE_UINT64_T
AC_TYPE_SIZE_T
AC_CHECK_FUNCS([floor])
AC_CHECK_FUNCS([memset])
AC_CHECK_FUNCS([sqrt])
AC_C_INLINE
AC_CHECK_HEADERS([math.h])
AC_CHECK_HEADERS([string.h])
AC_CHECK_HEADERS([iostream])
AC_CHECK_HEADERS([stdio.h])
AC_CHECK_HEADERS([fstream])
AC_CHECK_HEADERS([vector])

#########

AM_CXXFLAGS="-g -g3 -gdwarf-4"
AC_SUBST([AM_CXXFLAGS])

AM_CPPFLAGS=""
AC_SUBST([AM_CPPFLAGS])

AM_LDFLAGS=""
AC_SUBST([AM_LDFLAGS])

AM_LIBS=""
AC_SUBST([AM_LIBS])

AC_OUTPUT
