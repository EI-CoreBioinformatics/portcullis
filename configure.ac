##########################################################################
# Setup
##########################################################################

# Autoconf setup
AC_PREREQ([2.68])
AC_INIT([portcullis],[1.1.0],[daniel.mapleson@earlham.ac.uk],[portcullis],[http://www.earlham.ac.uk])
AC_CONFIG_SRCDIR([src/portcullis.cc])
AC_CONFIG_AUX_DIR([build-aux])
AC_CONFIG_MACRO_DIR([m4])
AC_CONFIG_HEADERS([config.h])
AC_CONFIG_FILES([Makefile doc/Makefile doc/source/conf.py lib/portcullis.pc lib/Makefile src/Makefile tests/compat.sh tests/Makefile
                scripts/junctools/setup.py scripts/junctools/junctools/__init__.py
                scripts/portcullis/setup.py scripts/portcullis/portcullis/__init__.py
                scripts/Makefile])
AC_CONFIG_SUBDIRS([deps/htslib-1.3 deps/ranger-0.3.8])
: ${CXXFLAGS=""}
: ${CFLAGS=""}
AC_LANG([C++])
AC_LANG_CPLUSPLUS
AC_PROG_CXXCPP
AC_PROG_CXX
AC_PROG_INSTALL

# Ensure there's C++11 support
AX_CXX_COMPILE_STDCXX_11([noext],[mandatory])
CXXFLAGS=""
CFLAGS=""

# Automake setup
AM_INIT_AUTOMAKE([1.11 subdir-objects foreign -Wall no-define tar-ustar])
m4_pattern_allow([AM_PROG_AR])
AM_PROG_AR
AM_SILENT_RULES([yes])

# Libtool setup
LT_INIT

# Detect OS
UNAME_S=`uname -s`
if test x"${UNAME_S}" == x"Linux"; then
    AC_DEFINE([OS_LINUX], [1], [Linux detected])
fi
if test x"${UNAME_S}" == x"Darwin"; then
    AC_DEFINE([OS_MAC], [1], [Mac detected])
fi

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

# Required for HTS lib
AC_SEARCH_LIBS([dlopen], [dl],
    [DL_LIB="-ldl"],
    [AC_MSG_ERROR([dl lib not found.  Please ensure that dl is properly built and configured.])])

AC_SEARCH_LIBS([sin], [m],
    [M_LIB="-lm"],
    [AC_MSG_ERROR([libm not found.  Please ensure that maths lib is properly built and configured.])])

# Required for Boost.  Doesn't exist on Mac
if test x"${UNAME_S}" != x"Darwin"; then
    AC_SEARCH_LIBS([clock_gettime], [rt],
        [RT_LIB="-lrt"],
        [AC_MSG_ERROR([rt lib not found.  Please ensure that rt is properly built and configured.])])
else
    RT_LIB=""
fi

# Check pthreads, zlib exist and are available
AX_PTHREAD(, [AC_MSG_ERROR([pthreads lib not found.  Please ensure that pthreads is properly built and configured.])])
AX_CHECK_ZLIB(, [AC_MSG_ERROR([zlib not found.  Please ensure that zlib is properly built and configured.])])

# Check for samtools (reqired)
AC_CHECK_PROG([with_samtools],[samtools],[yes],[no])
if test "$with_samtools" == no ; then
    AC_MSG_ERROR([samtools not found])
fi

AX_PYTHON_DEVEL([>= '3.5'])

if [[ -z "${PYTHON_VERSION}" ]]; then
    AC_MSG_ERROR([Python3 not detected. If running on a debian system try "sudo apt-get install python3-dev".])
fi

if [[ -z "${PYTHON_EXTRA_LIBS}" ]]; then
    AC_MSG_ERROR([Python3 detected but Python3 development library was not found. If running on a debian system try "sudo apt-get install python3-dev".])
fi

pybin=python${PYTHON_VERSION}

# Check for required python libraries
AX_PYTHON_MODULE(pandas, 1, ${pybin})
#AX_PYTHON_MODULE(matplotlib, [AC_MSG_ERROR([matplotlib not found])], ${pybin})
#AX_PYTHON_MODULE(scipy, [AC_MSG_ERROR([scipy not found])], ${pybin})
#AX_PYTHON_MODULE(sphinx, [AC_MSG_ERROR([sphinx not found])], ${pybin})
#AX_PYTHON_MODULE(sklearn, [AC_MSG_ERROR([sklearn not found])], ${pybin})

AC_ARG_ENABLE([py-install],
    AS_HELP_STRING([--disable-py-install], [Disable installation of python packages.  If this is disabled, portcullis will not install the python portcullis and junctools packages to <prefix>/lib/python<version>/site-packages.  The user can manually install the portcullis and junctools packages using setuptools.]), make_py_install="no", make_py_install="yes")

AM_CONDITIONAL([MAKE_PY_INSTALL], [test x$make_py_install = xyes])
AC_SUBST([MAKE_PY_INSTALL])


# Check for sphinx executable to build documentation
AC_CHECK_PROG([sphinx], [sphinx-build], [yes], [no])
if [[ "${sphinx}" == "yes" ]]; then
    #Now check if the version is correct
    sphinx_major_version=`sphinx-build --version 2>&1 | cut -d' ' -f3 | cut -d '.' -f 1 | tr -d " \t\n\r"`
    sphinx_minor_version=`sphinx-build --version 2>&1 | cut -d' ' -f3 | cut -d '.' -f 2 | tr -d " \t\n\r"`
    if test x"${sphinx_major_version}" = x""; then sphinx_major_version=0; fi
    if test x"${sphinx_minor_version}" = x""; then sphinx_minor_version=0; fi
    major_size=${#sphinx_major_version}
    minor_size=${#sphinx_minor_version}
    if (( ("${major_size}" > 2) || ("${minor_size}" > 2) )); then sphinx_major_version=0; sphinx_minor_version=0; fi
    if (( ("${sphinx_major_version}" < 1 ) || (  ( "${sphinx_major_version}" >= 1 ) && ( "${sphinx_minor_version}" < 3 ) ) )); then
        AC_MSG_WARN([Sphinx detected, but version is too old.  Cant create documentation.])
        sphinx="no"
    fi
else
    AC_MSG_WARN([Sphinx not detected, cant create documentation.])
    sphinx="no"
fi
AM_CONDITIONAL([MAKE_DOCS], [test x$sphinx = xyes])
AC_SUBST([MAKE_DOCS])

AM_CXXFLAGS="-g -O2 -std=c++11"
AC_SUBST([AM_CXXFLAGS])

AM_CPPFLAGS="${ZLIB_CPPFLAGS} ${PYTHON_CPPFLAGS}"
AC_SUBST([AM_CPPFLAGS])

AM_LDFLAGS="-static-libstdc++ ${ZLIB_LDFLAGS} ${PYTHON_RPATH}"
AC_SUBST([AM_LDFLAGS])

AM_LIBS="${PTHREAD_CFLAGS} ${ZLIB_LIB} ${RT_LIB} ${PYTHON_EXTRA_LIBS} ${PYTHON_LIBS}"
AC_SUBST([AM_LIBS])

AC_OUTPUT

if [[ "${sphinx}" == "no" ]]; then
    AC_MSG_WARN([Sphinx not detected, or version is too old.  Wont create documentation.  You can still find the documentation online at: https://portcullis.readthedocs.io/en/latest/])
fi
