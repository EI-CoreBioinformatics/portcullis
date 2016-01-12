##########################################################################
# Setup
##########################################################################

# Autoconf setup
AC_PREREQ([2.68])
AC_INIT([portcullis],[0.11.0],[daniel.mapleson@tgac.ac.uk],[portcullis],[http://www.tgac.ac.uk])
AC_CONFIG_SRCDIR([src/portcullis.cc])
AC_CONFIG_AUX_DIR([build-aux])
AC_CONFIG_MACRO_DIR([m4])
AC_CONFIG_HEADERS([config.h])
AC_CONFIG_FILES([Makefile doc/Makefile doc/source/conf.py lib/portcullis-0.11.pc lib/Makefile src/Makefile tests/compat.sh tests/Makefile])
AC_CONFIG_SUBDIRS([deps/htslib-1.3])
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
AC_CHECK_HEADERS([glob.h])

# Save libs (restore as and when needed)
SAVED_LIBS="${LIBS}"

# Check pthreads, zlib exist and are available
AX_PTHREAD(, [AC_MSG_ERROR([pthreads lib not found.  Please ensure that pthreads is properly built and configured.])])
AX_CHECK_ZLIB

# Check for dlib (required)
AX_DLIB

AC_CHECK_FUNCS(clock_gettime, [], [
  AC_CHECK_LIB(rt, clock_gettime, [
    AC_DEFINE(HAVE_CLOCK_GETTIME, 1)
    RT_LIB="-lrt"    
  ])
])

# Check for samtools (reqired)
AC_CHECK_PROG([with_samtools],[samtools],[yes],[no])
if test "$with_samtools" == no ; then
  AC_MSG_ERROR([samtools not found])
fi 

# Check for sphinx to build documentation (optional)
AC_CHECK_PROG([sphinx], [sphinx-build], [yes], [no])
if test x"${sphinx}" = x"no"; then
    AC_MSG_WARN([Sphinx not detected, won't create documentation.  You can still find the documentation online at: https://portcullis.readthedocs.org/en/latest/])
    MAKE_DOCS=0
else
    MAKE_DOCS=1
fi

AC_SUBST([MAKE_DOCS])

## Check for boost (required)
AX_BOOST_BASE([1.52],, [AC_MSG_ERROR([Boost not found.  Please ensure that boost is properly built and the BOOST_ROOT environment variable is set.  Alternatively you can override BOOST_ROOT with the --with-boost option.])])
AX_BOOST_SYSTEM
AX_BOOST_FILESYSTEM
AX_BOOST_CHRONO
AX_BOOST_PROGRAM_OPTIONS
AX_BOOST_TIMER

# Combine BOOST variables (apart for BOOST_TEST)
BOOST_LIBS="${BOOST_TIMER_STATIC_LIB} ${BOOST_CHRONO_STATIC_LIB} ${BOOST_FILESYSTEM_STATIC_LIB} ${BOOST_PROGRAM_OPTIONS_STATIC_LIB} ${BOOST_SYSTEM_STATIC_LIB}"
AC_SUBST([BOOST_LIBS])

LIBS="${SAVED_LIBS}"


#########

AM_CXXFLAGS="-DCPLUSPLUS"
AC_SUBST([AM_CXXFLAGS])

AM_CPPFLAGS="${BOOST_CPPFLAGS} ${DLIB_CPPFLAGS} ${ZLIB_CPPFLAGS}"
AC_SUBST([AM_CPPFLAGS])

AM_LDFLAGS="${BOOST_LDFLAGS} ${PTHREAD_CFLAGS} ${ZLIB_LDFLAGS}"
AC_SUBST([AM_LDFLAGS])

AM_LIBS="${BOOST_LIBS} ${ZLIB_LIB} ${RT_LIB}"
AC_SUBST([AM_LIBS])

AC_OUTPUT
