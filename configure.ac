##########################################################################
# Setup
##########################################################################

# Autoconf setup
AC_PREREQ([2.68])
AC_INIT([portcullis],[1.0.2],[daniel.mapleson@earlham.ac.uk],[portcullis],[http://www.earlham.ac.uk])
AC_CONFIG_SRCDIR([src/portcullis.cc])
AC_CONFIG_AUX_DIR([build-aux])
AC_CONFIG_MACRO_DIR([m4])
AC_CONFIG_HEADERS([config.h])
AC_CONFIG_FILES([Makefile doc/Makefile doc/source/conf.py lib/portcullis.pc lib/Makefile src/Makefile tests/compat.sh tests/Makefile scripts/setup.py scripts/junctools/__init__.py scripts/Makefile])
AC_CONFIG_SUBDIRS([deps/htslib-1.3 deps/ranger-0.3.8])
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

AC_ARG_ENABLE([junctools],
    AS_HELP_STRING([--disable-junctools], [Disable build of python package junctools]), make_junctools="no", make_junctools="yes")

if [[ "${make_junctools}" == "yes" ]]; then
    
    AX_PYTHON_DEVEL([>= '3.1'])

    if [[ -z "${PYTHON_VERSION}" ]]; then
        AC_MSG_ERROR([Python3 not detected. If running on a debian system try "sudo apt-get install python3-dev".])
    fi    

    if [[ -z "${PYTHON_EXTRA_LIBS}" ]]; then
        AC_MSG_ERROR([Python3 detected but Python3 development library was not found. If running on a debian system try "sudo apt-get install python3-dev".])
    fi

    pybin=python${PYTHON_VERSION}

    # Check for required python libraries
    #AX_PYTHON_MODULE(numpy, [AC_MSG_ERROR([numpy not found])], ${pybin})
    #AX_PYTHON_MODULE(matplotlib, [AC_MSG_ERROR([matplotlib not found])], ${pybin})
    #AX_PYTHON_MODULE(scipy, [AC_MSG_ERROR([scipy not found])], ${pybin})
    #AX_PYTHON_MODULE(sphinx, [AC_MSG_ERROR([sphinx not found])], ${pybin})
    #AX_PYTHON_MODULE(sklearn, [AC_MSG_ERROR([sklearn not found])], ${pybin})

    
fi
AM_CONDITIONAL([MAKE_JUNCTOOLS], [test x$make_junctools = xyes])
AC_SUBST([MAKE_JUNCTOOLS])


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



## Check for boost (required)
AC_ARG_ENABLE([dynamic-boost], AS_HELP_STRING([--enable-dynamic-boost], [Boost is statically linked by default.  Use this option if you would prefer to dynamically link boost]), dynboost="yes", dynboost="no")
AX_BOOST_BASE([1.52],, [AC_MSG_ERROR([Boost not found.  Please ensure that boost is properly built and the BOOST_ROOT environment variable is set.  Alternatively you can override BOOST_ROOT with the --with-boost option.])])
AX_BOOST_SYSTEM
AX_BOOST_FILESYSTEM
AX_BOOST_CHRONO
AX_BOOST_PROGRAM_OPTIONS
AX_BOOST_TIMER


# Combine BOOST variables (apart for BOOST_TEST)
BOOST_DYN_LIBS="${BOOST_TIMER_LIB} ${BOOST_CHRONO_LIB} ${BOOST_FILESYSTEM_LIB} ${BOOST_PROGRAM_OPTIONS_LIB} ${BOOST_SYSTEM_LIB}"
AC_SUBST([BOOST_DYN_LIBS])

BOOST_STATIC_LIBS="${BOOST_TIMER_STATIC_LIB} ${BOOST_CHRONO_STATIC_LIB} ${BOOST_FILESYSTEM_STATIC_LIB} ${BOOST_PROGRAM_OPTIONS_STATIC_LIB} ${BOOST_SYSTEM_STATIC_LIB} "
AC_SUBST([BOOST_STATIC_LIBS])

# Prefer static linking of boost.  If one or more static libs are not present 
# (this happens with some default boost installations on some environments) then
# use dynamic libs instead.
if [[ "${dynboost}" == "yes" ]]; then
    BOOST_LIBS="${BOOST_DYN_LIBS}"
else
    if [[ -z "${BOOST_TIMER_STATIC_LIB}" ]] || [[ -z "${BOOST_CHRONO_STATIC_LIB}" ]] || [[ -z "${BOOST_FILESYSTEM_STATIC_LIB}" ]] || [[ -z "${BOOST_PROGRAM_OPTIONS_STATIC_LIB}" ]] || [[ -z "${BOOST_SYSTEM_STATIC_LIB}" ]]; then
        AC_MSG_WARN([Not all static boost libraries could be found.  Will use dynamic libraries instead.])
        BOOST_LIBS="${BOOST_DYN_LIBS}"
        dynboost="yes"
    else
        BOOST_LIBS="${BOOST_STATIC_LIBS}"
    fi
fi
AC_SUBST([BOOST_LIBS])

LIBS="${SAVED_LIBS}"


#########

AM_CXXFLAGS="-g -g3 -gdwarf-4"
AC_SUBST([AM_CXXFLAGS])

AM_CPPFLAGS="${BOOST_CPPFLAGS} ${ZLIB_CPPFLAGS}"
AC_SUBST([AM_CPPFLAGS])

AM_LDFLAGS="${BOOST_LDFLAGS} ${PTHREAD_CFLAGS} ${ZLIB_LDFLAGS}"
AC_SUBST([AM_LDFLAGS])

# Order is important!
AM_LIBS="${M_LIB} ${DL_LIB} ${BOOST_LIBS} ${RT_LIB} ${ZLIB_LIB} ${PTHREAD_LIBS}"
AC_SUBST([AM_LIBS])

AC_OUTPUT

if [[ "$make_junctools" == "no" ]]; then
    AC_MSG_NOTICE([Junctools not created at user request])    
fi

if [[ "${sphinx}" == "no" ]]; then
    AC_MSG_WARN([Sphinx not detected, or version is too old.  Wont create documentation.  You can still find the documentation online at: https://portcullis.readthedocs.io/en/latest/])    
fi

if [[ "${dynboost}" == "yes" ]]; then
    AC_MSG_NOTICE([Boost will be dynamically linked to your executable])
else
    AC_MSG_NOTICE([Boost will be statically linked to your executable])
fi
