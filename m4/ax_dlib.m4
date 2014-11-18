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

AC_DEFUN([AX_DLIB],
[
    AC_MSG_CHECKING([if dlib is wanted])
    AC_ARG_WITH([dlib],
       [AS_HELP_STRING([--with-dlib=PREFIX],
	        [Search for dlib in DIR/include.  Use if your dlib library is installed in a non-standard location.])],
        [if test "$withval" != "no" ; then
            AC_MSG_RESULT([yes])
            if test -d "$withval" ; then
                DLIB_PATH="$withval"
            else
                AC_MSG_WARN([Sorry, $withval does not exist, checking usual places])
            fi
        else
            AC_MSG_RESULT([no])
        fi], 
        [AC_MSG_RESULT([yes])]
    )

    if test -f "${DLIB_PATH}/include/dlib/svm.h" ; then
        DLIB_CPPFLAGS="-I${DLIB_PATH}/include"        
    elif test -f "${DLIB_PATH}/dlib/svm.h" ; then
        DLIB_CPPFLAGS="-I${DLIB_PATH}"
    elif test -f "/usr/local/include/dlib/svm.h" ; then
        DLIB_PATH="/usr/local"
        DLIB_CPPFLAGS="-I${DLIB_PATH}/include"
    else
        DLIB_PATH="/usr"
        DLIB_CPPFLAGS="-I${DLIB_PATH}/include"
    fi

    # Update vars for bamtools if required
    if test -n "${DLIB_PATH}" ; then
        
        OLD_CPPFLAGS=${CPPFLAGS}
        CPPFLAGS="${AM_CPPFLAGS} ${CPPFLAGS} ${DLIB_CPPFLAGS}"

        #echo "CPPFLAGS=${CPPFLAGS}"

        AC_LANG_PUSH(C++)

        # Check header exists
        AC_CHECK_HEADER([dlib/svm.h], [], [
            AC_MSG_ERROR([dlib headers not found.  Please ensure that the dlib headers directory can be found on the CPPFLAGS env var. Alternatively, you can try the --with-dlib option, which expects to find the dlib headers at <dlib_dir>/include/dlib/<headers>.])
        ])

        AC_DEFINE(HAVE_DLIB,,[define if the dlib is available])

        AC_LANG_POP(C++)

        # Restore the previous environment variables if required
        CPPFLAGS=${OLD_CPPFLAGS}
    fi
])