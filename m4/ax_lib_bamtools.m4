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

AC_DEFUN([AX_LIB_BAMTOOLS],
[
    AC_MSG_CHECKING([if bamtools is wanted])
    AC_ARG_WITH([bamtools],
        [AS_HELP_STRING([--with-bamtools=DIR],
                [Search for bamtools in DIR/include and DIR/lib.  Use if your bamtools was installed to a non-standard location.])],
        [if test "$withval" != "no" ; then
            AC_MSG_RESULT([yes])
            if test -d "$withval" ; then
                BAMTOOLS_HOME="$withval"
            else
                AC_MSG_WARN([Sorry, $withval does not exist, checking usual places])
            fi
        else
            AC_MSG_RESULT([no])
        fi], 
        [AC_MSG_RESULT([yes])]
    )

    if test -f "${BAMTOOLS_HOME}/include/bamtools/api/BamReader.h" ; then
        BAMTOOLS_CPPFLAGS="-I${BAMTOOLS_HOME}/include/bamtools"
        BAMTOOLS_LDFLAGS="-L${BAMTOOLS_HOME}/lib"
    elif test -f "${BAMTOOLS_HOME}/include/api/BamReader.h" ; then
        BAMTOOLS_CPPFLAGS="-I${BAMTOOLS_HOME}/include"
        BAMTOOLS_LDFLAGS="-L${BAMTOOLS_HOME}/lib"
    elif test -f "${BAMTOOLS_HOME}/api/BamReader.h" ; then
        BAMTOOLS_CPPFLAGS="-I${BAMTOOLS_HOME}"
        BAMTOOLS_LDFLAGS="-L${BAMTOOLS_HOME}"
    elif test -f "/usr/local/include/bamtools/api/BamReader.h" ; then
        BAMTOOLS_HOME="/usr/local"
        BAMTOOLS_CPPFLAGS="-I${BAMTOOLS_HOME}/include/bamtools"
        BAMTOOLS_LDFLAGS="-L${BAMTOOLS_HOME}/lib"
    else
        BAMTOOLS_HOME="/usr"
        BAMTOOLS_CPPFLAGS="-I${BAMTOOLS_HOME}/include/bamtools"
        BAMTOOLS_LDFLAGS=""
    fi

    # Update vars for bamtools if required
    if test -n "${BAMTOOLS_HOME}" ; then
        
        OLD_CPPFLAGS=${CPPFLAGS}
        OLD_LDFLAGS=${LDFLAGS}
        OLD_LIBS=${LIBS}
        CPPFLAGS="${AM_CPPFLAGS} ${CPPFLAGS} ${BAMTOOLS_CPPFLAGS}"
        LDFLAGS="${AM_LDFLAGS} ${LDFLAGS} ${BAMTOOLS_LDFLAGS}"
    
        AC_LANG_PUSH(C++)

        # Check header exists
        AC_CHECK_HEADER([api/BamMultiReader.h], [], [
            AC_MSG_ERROR([Bamtools headers not found. Please ensure that the bamtools include directory can be found on the CPPFLAGS env var. Alternatively, you can try the --with-bamtools option, which assumes the bamtools headers can be found at <bamtools_dir>/include/bamtools/<headers>.])
        ])

        # Check lib and function exists
        AX_CPP_CHECK_LIB(bamtools,
            [#include <api/BamMultiReader.h>],
            [],
            [],
            [AC_MSG_ERROR([Bamtools library not found. Please ensure that the bamtools libraries can be found on the LDFLAGS env var.  Alternatively, you can try the --with-bamtools option, which assumes the bamtools libraries can be found at <bamtools_dir>/lib/<libs>.])]
        )

        AC_LANG_POP(C++)
        
        BAMTOOLS_LIB="-lbamtools"

        # Restore the previous environment variables if required
        CPPFLAGS=${OLD_CPPFLAGS}
        LDFLAGS=${OLD_LDFLAGS}
        LIBS=${OLD_LIBS}
    fi
])