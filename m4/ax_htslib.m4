#  ********************************************************************
#  This file is part of Portculis.
#
#  Portculis is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  Portculis is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with Portculis.  If not, see <http://www.gnu.org/licenses/>.
#  *******************************************************************

AC_DEFUN([AX_HTSLIB],
[
    AC_ARG_WITH([htslib],
        [AS_HELP_STRING([--with-htslib=PREFIX],
                [Use this if your htslib installation prefix was set to a non-standard location.])],
        [HTSLIB_PATH=$withval], [HTSLIB_PATH=0])

    # Update vars for htslib if required
    if test ${HTSLIB_PATH} != 0; then
        HTSLIB_CPPFLAGS="-I${HTSLIB_PATH}/include"
        HTSLIB_LDFLAGS="-L${HTSLIB_PATH}/lib"
        HTSLIB_LIB="-lhts"
        
        # We need to do this in case we are using samtools
        export HTSLIB_CPPFLAGS
        export HTSLIB_LDFLAGS
        
        OLD_CPPFLAGS=${CPPFLAGS}
        OLD_LDFLAGS=${LDFLAGS}
        CPPFLAGS="${AM_CPPFLAGS} ${CPPFLAGS} ${HTSLIB_CPPFLAGS}"
        LDFLAGS="${AM_LDFLAGS} ${LDFLAGS} ${HTSLIB_LDFLAGS}"
    fi

    AC_LANG_PUSH([C])

    # Check header exists
    AC_CHECK_HEADER([htslib/faidx.h], [], [
        AC_MSG_ERROR([HTSlib headers not found. Please ensure that the htslib headers can be found on the CPPFLAGS env var. Alternatively, you can try the --with-htslib option, which assumes you have the following setup <htslib_dir>/include/htslib/<headers>.])
    ])

    # Check lib and function exists
    AC_CHECK_LIB(hts, fai_build, [], [
        AC_MSG_ERROR([HTSlib library not found. Please ensure that the HTSlib libs can be found on the LDFLAGS env var.  Alternatively, you can try the --with-htslib option, which assumes you have the following setup <htslib_dir>/lib/<libs>.])
    ])

    AC_LANG_POP([C])

    AC_DEFINE(HAVE_HTSLIB,,[define if the hts library is available])
    
    # Restore the previous environment variables if required
    if test ${HTSLIB_PATH} != 0; then
        CPPFLAGS=${OLD_CPPFLAGS}
        LDFLAGS=${OLD_LDFLAGS}
    fi
])