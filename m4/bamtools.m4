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

AC_DEFUN([BAMTOOLS],
[
    AC_ARG_WITH([bamtools],
        [AS_HELP_STRING([--with-bamtools=PREFIX],
                [Use this if your bamtools installation prefix was set to a non-standard location.])],
        [BAMTOOLS_PATH=$withval], [BAMTOOLS_PATH=0])

    # Update vars for bamtools if required
    if test ${BAMTOOLS_PATH} != 0; then
        BAMTOOLS_CPPFLAGS="-I${BAMTOOLS_PATH}/include -I${BAMTOOLS_PATH}/include/bamtools"
        BAMTOOLS_LDFLAGS="-L${BAMTOOLS_PATH}/lib -L${BAMTOOLS_PATH}/lib64"

        AC_SUBST([BAMTOOLS_CPPFLAGS])
        AC_SUBST([BAMTOOLS_LDFLAGS])
        
        OLD_CPPFLAGS=${CPPFLAGS}
        OLD_LDFLAGS=${LDFLAGS}
        CPPFLAGS="${AM_CPPFLAGS} ${CPPFLAGS} ${BAMTOOLS_CPPFLAGS}"
        LDFLAGS="${AM_LDFLAGS} ${LDFLAGS} ${BAMTOOLS_LDFLAGS}"
    fi

    # Check header exists
    AC_CHECK_HEADER([api/BamMultiReader.h], [], [
        AC_MSG_ERROR([Bamtools not found. Please ensure that the bamtools include directory can be found on the CPPFLAGS env var. Also check that the Bamtools lib directory can be found on the LDFLAGS env var.  Alternatively, you can try the --with-bamtools option, which assumes you have an include and lib dir under a bamtools build directory.])
    ])
    
    # Restore the previous environment variables if required
    if test ${BAMTOOLS_PATH} != 0; then
        CPPFLAGS=${OLD_CPPFLAGS}
        LDFLAGS=${OLD_LDFLAGS}
    fi
])