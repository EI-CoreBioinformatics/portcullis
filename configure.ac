##########################################################################
# Setup
##########################################################################

# Autoconf initialistion.  Sets package name version and contact details
AC_INIT([Portculis],[0.1.1],[daniel.mapleson@tgac.ac.uk],[portculis],[http://www.tgac.ac.uk])

# I built this with autoconf 2.53 and have only tested it on 2.69 systems.  You can change this but I don't
# make any guarentees
AC_PREREQ([2.53])

# Check that the configure script is being run from the right directory
AC_CONFIG_SRCDIR([src/portculis.cc])

# Put auxillary build scripts in 'build-aux' directory
AC_CONFIG_AUX_DIR([build-aux])

# Set the macro directory
AC_CONFIG_MACRO_DIR([m4])

# Automake initialisation.  Sets some default flags.  We use 'foreign' because we have a markdown README file.
# I've only tried this with automake 1.11 and above, so you require at least that version installed on your machine if you want to regenerate 'the Makefile.in's.  If you're looking at this comment in the 'configure' script then this the am__api_version, will be set to whatever I used to create this file.  Or it might be whatever you used if you made it yourself with autoconf (or autoreconf, or the autogen.sh script)
AM_INIT_AUTOMAKE([1.11 foreign -Wall no-define])

# Most of portculis is written in C++.  This is what this macro means: "Do compilation tests using CXX and CXXCPP and use extension .C for test programs. Use compilation flags: CPPFLAGS with CXXCPP, and both CPPFLAGS and CXXFLAGS with CXX."
AC_LANG_CPLUSPLUS

# Set output variable CXXCPP to a command that runs the C++ preprocessor. If ‘$CXX -E’ doesn't work, /lib/cpp is used. It is portable to run CXXCPP only on files with a .c, .C, .cc, or .cpp extension.
AC_PROG_CXXCPP

# AC_PROG_CXX determine a C++ compiler to use. Check whether the environment variable CXX or CCC (in that order) is set; if so, then set output variable CXX to its value.  We ignore the CXXFLAGS already set in the environment. This is necessary to avoid AC_PROG_XX automatically adding "-g -O2" to CXXFLAGS
CXXFLAGS=""
AC_PROG_CXX

# Determine a C compiler to use. If CC is not already set in the environment, check for gcc and cc, then for other C compilers. Set output variable CC to the name of the compiler found.
# NOTE: This probably isn't necessary.  I don't think we have any C code in portculis.
CFLAGS=""
AC_PROG_CC

# Delegates building of libraries to libtool
LT_PREREQ([2.4.2])
LT_INIT

# Create the 'install-sh' script, which is used as the 'install' make target
AC_PROG_INSTALL



######################################################################
# Debug compilation support.
######################################################################

AC_MSG_CHECKING([whether to build with debug information])
AC_ARG_ENABLE([debug],
    [AS_HELP_STRING([--enable-debug],
        [enable debug data generation (def=no)])],
    [debugit="$enableval"],
    [debugit=no])
AC_MSG_RESULT([$debugit])

# If in debug mode make a really slow binary with debug symbols, otherwise make a fast one without debug symbols.
# Generally CPPFLAGS is reserved for include paths and use CXXFLAGS is for c++ compiler settings.
if test x"$debugit" = x"yes"; then
    AC_DEFINE([DEBUG],[],[Debug Mode])
    AM_CXXFLAGS="-g -Wno-uninitialized -O0"
else
    AC_DEFINE([NDEBUG],[],[No-debug Mode])
    AM_CXXFLAGS="-O2"
fi

######################################################################


######################################################################
# DOXYGEN SUPPORT
######################################################################

DX_HTML_FEATURE(ON)
DX_CHM_FEATURE(OFF)
DX_CHI_FEATURE(OFF)
DX_MAN_FEATURE(ON)
DX_RTF_FEATURE(OFF)
DX_XML_FEATURE(OFF)
DX_PDF_FEATURE(ON)
DX_PS_FEATURE(OFF)

DX_INIT_DOXYGEN([$PACKAGE_NAME],[doxygen.cfg])

######################################################################


######################################################################
# Dependency checks
######################################################################


# Check some (almost) standard functionality is present that we require to run KAT.
AC_CHECK_HEADER_STDBOOL
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
AC_CHECK_HEADERS([string.h])
AC_CHECK_HEADERS([iostream])
AC_CHECK_HEADERS([stdio.h])
AC_CHECK_HEADERS([fstream])
AC_CHECK_HEADERS([vector])
AC_CHECK_HEADERS([getopt.h])



# Check pthreads and zlib exist and are available
AC_CHECK_LIB(pthread, pthread_exit)
AC_CHECK_LIB(z, gzread)

## Check for boost
BOOST_REQUIRE([1.52])

# Header only
BOOST_STRING_ALGO
BOOST_FOREACH

# Header + lib
BOOST_FILESYSTEM
BOOST_THREADS
BOOST_TEST
BOOST_SYSTEM


# Allow users to specify custom directories for seqan library and doxygen.  This might save some users from having to set the environment variables manually.

AC_ARG_WITH([seqan],
        [AS_HELP_STRING([--with-seqan=PREFIX],
                [Use this if your seqan library installation prefix was set to a non-standard location.])],
        [SEQAN_PATH=$withval], [SEQAN_PATH=0])

AC_ARG_WITH([doxygen],
        [AS_HELP_STRING([--with-doxygen=PREFIX],
		[Use this if your doxygen installation prefix was set to a non-standard location. Doxygen is only required for building source code documentation.  i.e. if doxygen is not available you will not be able to run [make doxygen-doc] but KAT will still work fine without this.])],
        [DOXYGEN_PATH=$withval], [DOXYGEN_PATH=0])


# Update vars for seqan if required
if test ${SEQAN_PATH} != 0; then
        AM_CPPFLAGS="-I${SEQAN_PATH}/include $AM_CPPFLAGS"
fi

# Update vars for doxygen if required
if test ${DOXYGEN_PATH} != 0; then
	AM_PATH="${DOXYGEN_PATH}/bin:${AM_PATH}"
fi

######### Header checks (need to use the actual command line vars here so we need to store and restore these)

OLD_CPPFLAGS=${CPPFLAGS}
OLD_LDFLAGS=${LDFLAGS}
CPPFLAGS="${AM_CPPFLAGS} ${CPPFLAGS}"
LDFLAGS="${AM_LDFLAGS} ${LDFLAGS}"


#### Put any header or lib checks here ####

AC_CHECK_HEADER([seqan/sequence.h], [], [
         echo "Seqan-library is required to compile and run portculis.  Please make sure seqan-library is installed on your system."
         echo "If seqan-library is installed but the headers are not currently available through the CPPFLAGS environment"
         echo "variable then please re-run \"./configure\".  You can optionally use the \"--with-seqan=prefix\" argument"
         echo "to specify the root directory of your seqan-library installation."
         exit -1],
         [])



# Restore the previous CPPFLAGS
CPPFLAGS=${OLD_CPPFLAGS}
LDFLAGS=${OLD_LDFLAGS}



#########

echo "Transfering variables:"
echo " BOOST_CPPFLAGS: ${BOOST_CPPFLAGS}"
echo " BOOST_FILESYSTEM_LDFLAGS : ${BOOST_FILESYSTEM_LDFLAGS}"
echo " BOOST_FILESYSTEM_LIBS : ${BOOST_FILESYSTEM_LIBS}"


#########

AC_SUBST([AM_CPPFLAGS])
AC_SUBST([AM_CXXFLAGS])
AC_SUBST([AM_CFLAGS])
AC_SUBST([AM_LDFLAGS])
AC_SUBST([AM_PATH])
AC_SUBST([AM_LIBS])
AC_SUBST([BOOST_CPPFLAGS])
AC_SUBST([BOOST_FILESYSTEM_LDFLAGS])
AC_SUBST([BOOST_FILESYSTEM_LIBS])
AC_SUBST([BOOST_SYSTEM_LDFLAGS])
AC_SUBST([BOOST_SYSTEM_LIBS])

AC_CONFIG_HEADERS([config.h])
AC_CONFIG_FILES([Makefile src/Makefile tests/Makefile])
AC_OUTPUT
