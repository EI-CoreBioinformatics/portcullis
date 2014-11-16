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

AC_DEFUN([SAMTOOLS],
[
    AC_ARG_WITH([samtools],
        [AS_HELP_STRING([--with-samtools=PREFIX],
                [Use this if your samtools installation prefix was set to a non-standard location.])],
        [SAMTOOLS_PATH=$withval], [SAMTOOLS_PATH=0])

    # Update vars for samtools if required
    if test ${SAMTOOLS_PATH} != 0; then
        SAMTOOLS_CPPFLAGS="-I${SAMTOOLS_PATH}/include -I${SAMTOOLS_PATH}/include/bam"
        SAMTOOLS_LDFLAGS="-L${SAMTOOLS_PATH}/lib -L${SAMTOOLS_PATH}/lib64"
        SAMTOOLS_BIN="${SAMTOOLS_BIN}/bin"

        AC_SUBST([SAMTOOLS_CPPFLAGS])
        AC_SUBST([SAMTOOLS_LDFLAGS])
        AC_SUBST([SAMTOOLS_BIN])

        OLD_CPPFLAGS=${CPPFLAGS}
        OLD_LDFLAGS=${LDFLAGS}
        CPPFLAGS="${AM_CPPFLAGS} ${CPPFLAGS} ${SAMTOOLS_CPPFLAGS}"
        LDFLAGS="${AM_LDFLAGS} ${LDFLAGS} ${SAMTOOLS_LDFLAGS}"
    fi

    # Check header exists
    AC_CHECK_HEADER([faidx.h], [], [
        AC_MSG_ERROR([Samtools headers not found. Please ensure that the samtools include directory can be found on the CPPFLAGS env var. Alternatively, you can try the --with-samtools option, which assumes you have an include and lib dir under a samtools build directory.])
    ])

    # Check lib and function exists
    AC_CHECK_LIB(bam, fai_build, [], [
        AC_MSG_ERROR([Samtools library not found. Please ensure that the Samtools lib directory can be found on the LDFLAGS env var.  Alternatively, you can try the --with-samtools option, which assumes you have an include and lib dir under a samtools build directory.])
    ])
    
    # Restore the previous environment variables if required
    if test ${SAMTOOLS_PATH} != 0; then
        CPPFLAGS=${OLD_CPPFLAGS}
        LDFLAGS=${OLD_LDFLAGS}
    fi
])