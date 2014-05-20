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
# Dependency checks
######################################################################

# Check some (almost) standard functionality is present that we require to run portculis.
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


# Check pthreads and zlib exist and are available
AC_CHECK_LIB(pthread, pthread_exit)
AC_CHECK_LIB(z, gzread)


## Check for boost
AX_BOOST_BASE([1.52],, [AC_MSG_ERROR([Boost not found.  Please ensure that boost is properly built and the BOOST_ROOT environment variable is set.  Alternatively you can override BOOST_ROOT with the --with-boost option.])])
AX_BOOST_FILESYSTEM
AX_BOOST_PROGRAM_OPTIONS
AX_BOOST_SYSTEM
AX_BOOST_THREAD
AX_BOOST_UNIT_TEST_FRAMEWORK

# Combine BOOST variables (apart for BOOST_TEST)
#BOOST_LDFLAGS="${BOOST_FILESYSTEM_LDFLAGS} ${BOOST_PROGRAM_OPTIONS_LDFLAGS} ${BOOST_SYSTEM_LDFLAGS} ${BOOST_THREAD_LDFLAGS}"
BOOST_LIBS="${BOOST_FILESYSTEM_LIB} ${BOOST_PROGRAM_OPTIONS_LIB} ${BOOST_SYSTEM_LIB} ${BOOST_THREAD_LIB}"



# Allow users to specify custom directories for samtools and bamtools.  
# This might save some users from having to set the environment variables manually.

AC_ARG_WITH([samtools],
        [AS_HELP_STRING([--with-samtools=PREFIX],
                [Use this if your samtools installation prefix was set to a non-standard location.])],
        [SAMTOOLS_PATH=$withval], [SAMTOOLS_PATH=0])

AC_ARG_WITH([bamtools],
        [AS_HELP_STRING([--with-bamtools=PREFIX],
                [Use this if your bamtools installation prefix was set to a non-standard location.])],
        [BAMTOOLS_PATH=$withval], [BAMTOOLS_PATH=0])


# Update vars for samtools if required
if test ${SAMTOOLS_PATH} != 0; then
        SAMTOOLS_CPPFLAGS="-I${SAMTOOLS_PATH}/include"
        SAMTOOLS_LDFLAGS="-L${SAMTOOLS_PATH}/lib"
fi

# Update vars for samtools if required
if test ${BAMTOOLS_PATH} != 0; then
        BAMTOOLS_CPPFLAGS="-I${BAMTOOLS_PATH}/include"
        BAMTOOLS_LDFLAGS="-L${BAMTOOLS_PATH}/lib"
fi


######### Header checks (need to use the actual command line vars here so we need to store and restore these)

OLD_CPPFLAGS=${CPPFLAGS}
OLD_LDFLAGS=${LDFLAGS}
CPPFLAGS="${BAMTOOLS_CPPFLAGS} ${SAMTOOLS_CPPFLAGS} ${SEQAN_CPPFLAGS} ${AM_CPPFLAGS} ${CPPFLAGS}"
LDFLAGS="${BAMTOOLS_LDFLAGS} ${SAMTOOLS_LDFLAGS} ${AM_LDFLAGS} ${LDFLAGS}"


#### Put any header or lib checks here ####

AC_CHECK_HEADER([api/BamMultiReader.h], [], [
    AC_MSG_ERROR([Bamtools not found. Please ensure that the bamtools include directory can be found on the CPPFLAGS env var. Also check that the Bamtools lib directory can be found on the LDFLAGS env var.  Alternatively, you can try the --with-bamtools option, which assumes you have an include and lib dir under a bamtools build directory.])
])

AC_CHECK_LIB(bam, fai_build, [], [
    AC_MSG_ERROR([Samtools not found. Please ensure that the Samtools include directory can be found on the CPPFLAGS env var. Also check that the Samtools lib directory can be found on the LDFLAGS env var.  Alternatively, you can try the --with-samtools option, which assumes you have an include and lib dir under a samtools build directory.])
])

# Restore the previous CPPFLAGS
CPPFLAGS=${OLD_CPPFLAGS}
LDFLAGS=${OLD_LDFLAGS}



#########

AC_SUBST([AM_CPPFLAGS])
AC_SUBST([AM_CXXFLAGS])
AC_SUBST([AM_CFLAGS])
AC_SUBST([AM_LDFLAGS])
AC_SUBST([AM_PATH])
AC_SUBST([AM_LIBS])
AC_SUBST([SAMTOOLS_CPPFLAGS])
AC_SUBST([SAMTOOLS_LDFLAGS])
AC_SUBST([BAMTOOLS_PATH])
AC_SUBST([BAMTOOLS_CPPFLAGS])
AC_SUBST([BAMTOOLS_LDFLAGS])
AC_SUBST([SEQAN_CPPFLAGS])
AC_SUBST([BOOST_CPPFLAGS])
AC_SUBST([BOOST_LDFLAGS])
AC_SUBST([BOOST_LIBS])
AC_SUBST([BOOST_TEST_LDFLAGS])
AC_SUBST([BOOST_TEST_LIBS])

AC_CONFIG_HEADERS([config.h])
AC_CONFIG_FILES([Makefile src/Makefile tests/Makefile])
AC_OUTPUT
