##########################################################################
# Setup
##########################################################################

# Autoconf initialistion.  Sets package name version and contact details
AC_INIT([Portcullis],[0.8.7],[daniel.mapleson@tgac.ac.uk],[portcullis],[http://www.tgac.ac.uk])

# I built this with autoconf 2.53 and have only tested it on 2.69 systems.  You can change this but I don't
# make any guarentees
AC_PREREQ([2.53])

# Check that the configure script is being run from the right directory
AC_CONFIG_SRCDIR([src/portcullis.cc])

# Put auxillary build scripts in 'build-aux' directory
AC_CONFIG_AUX_DIR([build-aux])

# Set the macro directory
AC_CONFIG_MACRO_DIR([m4])

# Automake initialisation.  Sets some default flags.  We use 'foreign' because we have a markdown README file.
# I've only tried this with automake 1.11 and above, so you require at least that version installed on your 
# machine if you want to regenerate 'the Makefile.in's.  If you're looking at this comment in the 'configure' 
#script then this the am__api_version, will be set to whatever I used to create this file.  Or it might be 
#whatever you used if you made it yourself with autoconf (or autoreconf, or the autogen.sh script)
AM_INIT_AUTOMAKE([1.11 subdir-objects foreign -Wall no-define tar-ustar])

# Most of portcullis is written in C++.  This is what this macro means: "Do compilation tests using CXX and CXXCPP and use extension .C for test programs. Use compilation flags: CPPFLAGS with CXXCPP, and both CPPFLAGS and CXXFLAGS with CXX."
AC_LANG_CPLUSPLUS

# Set output variable CXXCPP to a command that runs the C++ preprocessor. If ‘$CXX -E’ doesn't work, /lib/cpp is used. It is portable to run CXXCPP only on files with a .c, .C, .cc, or .cpp extension.
AC_PROG_CXXCPP

# AC_PROG_CXX determine a C++ compiler to use. Check whether the environment variable CXX or CCC (in that order) is set; if so, then set output variable CXX to its value.  We ignore the CXXFLAGS already set in the environment. This is necessary to avoid AC_PROG_XX automatically adding "-g -O2" to CXXFLAGS
CXXFLAGS=""
AC_PROG_CXX

# Determine a C compiler to use. If CC is not already set in the environment, check for gcc and cc, then for other C compilers. Set output variable CC to the name of the compiler found.
# NOTE: This probably isn't necessary.  I don't think we have any C code in portcullis.
CFLAGS=""
AC_PROG_CC

# Ensure there's C++11 support
AX_CXX_COMPILE_STDCXX_11(,[mandatory])

# Delegates building of libraries to libtool
LT_PREREQ([2.4.2])
LT_INIT

# Create the 'install-sh' script, which is used as the 'install' make target
AC_PROG_INSTALL



######################################################################
# Dependency checks
######################################################################

# Check some (almost) standard functionality is present that we require to run portcullis.
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

# Check pthreads and zlib exist and are available
AX_PTHREAD(, [AC_MSG_ERROR([pthreads lib not found.  Please ensure that pthreads is properly built and configured.])])
AX_CHECK_ZLIB(, [AC_MSG_ERROR([Zlib not found.  Please ensure that zlib is properly built and configured.])])

## Check for boost
AX_BOOST_BASE([1.52],, [AC_MSG_ERROR([Boost not found.  Please ensure that boost is properly built and the BOOST_ROOT environment variable is set.  Alternatively you can override BOOST_ROOT with the --with-boost option.])])
AX_BOOST_SYSTEM
AX_BOOST_FILESYSTEM
AX_BOOST_CHRONO
AX_BOOST_PROGRAM_OPTIONS
AX_BOOST_TIMER
AX_BOOST_UNIT_TEST_FRAMEWORK

# Combine BOOST variables (apart for BOOST_TEST)
BOOST_LIBS="${BOOST_TIMER_STATIC_LIB} ${BOOST_CHRONO_STATIC_LIB} ${BOOST_FILESYSTEM_STATIC_LIB} ${BOOST_PROGRAM_OPTIONS_STATIC_LIB} ${BOOST_SYSTEM_STATIC_LIB}"
AC_SUBST([BOOST_LIBS])

LIBS="${SAVED_LIBS}"


#########

AM_CXXFLAGS="-DCPLUSPLUS"
AC_SUBST([AM_CXXFLAGS])

AM_CPPFLAGS="${BOOST_CPPFLAGS} ${ZLIB_CPPFLAGS}"
AC_SUBST([AM_CPPFLAGS])

AM_LDFLAGS="${BOOST_LDFLAGS} ${ZLIB_LDFLAGS}"
AC_SUBST([AM_LDFLAGS])

AM_LIBS="${BOOST_LIBS} ${PTHREAD_CFLAGS} ${ZLIB_LIB}"
AC_SUBST([AM_LIBS])

AC_CONFIG_HEADERS([config.h])
AC_CONFIG_FILES([Makefile deps/htslib-1.2.1/Makefile deps/samtools-1.2/Makefile src/Makefile tests/Makefile])
AC_OUTPUT
