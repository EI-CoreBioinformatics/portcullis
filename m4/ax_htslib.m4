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
    AC_MSG_CHECKING([if htslib is wanted])
    AC_ARG_WITH(
        [htslib],
        AS_HELP_STRING(
            [--with-htslib=DIR],
            [search for htslib in DIR/include and DIR/lib]
        ),
        [if test "$withval" != no ; then
            AC_MSG_RESULT([yes])
            if test -d "$withval" ; then
                HTSLIB_PATH="$withval"
            else
                AC_MSG_WARN([Sorry, $withval does not exist, checking usual places])
            fi
        else
            AC_MSG_RESULT([no])
        fi],
        [AC_MSG_RESULT([yes])]
    )

    HTSLIB_LIB="-lhts"

    if test -f "${HTSLIB_PATH}/include/htslib/hts.h" ; then
        HTSLIB_CPPFLAGS="-I${HTSLIB_PATH}/include"
        HTSLIB_LDFLAGS="-L${HTSLIB_PATH}/lib"
    elif test -f "/usr/local/include/htslib/hts.h" ; then
        HTSLIB_PATH="/usr/local"
        HTSLIB_CPPFLAGS="-I${HTSLIB_PATH}/include"
        HTSLIB_LDFLAGS="-L${HTSLIB_PATH}/lib"
    else
        HTSLIB_PATH="/usr"
        HTSLIB_CPPFLAGS="-I${HTSLIB_PATH}/include"
        HTSLIB_LDFLAGS=""
    fi


    #
    # Locate samtools, if wanted
    #
    if test -n "${HTSLIB_PATH}" ; then
        
        OLD_CPPFLAGS=${CPPFLAGS}
        OLD_LDFLAGS=${LDFLAGS}
        CPPFLAGS="${AM_CPPFLAGS} ${CPPFLAGS} ${HTSLIB_CPPFLAGS}"
        LDFLAGS="${AM_LDFLAGS} ${LDFLAGS} ${HTSLIB_LDFLAGS}"
    
        AC_LANG_PUSH([C])

        # Check header exists
        AC_CHECK_HEADER([htslib/faidx.h], [ac_cv_faidx_h=yes], [ac_cv_faidx_h=no])

        # Check lib and function exists
        AC_CHECK_LIB(hts, fai_build, [ac_cv_libhts=yes], [ac_cv_libhts=no])

        AC_LANG_POP([C])

        #echo "${ac_cv_faidx_h}"
        #echo "${ac_cv_libhts}"

        CPPFLAGS=${OLD_CPPFLAGS}
        LDFLAGS=${OLD_LDFLAGS}

        AC_MSG_CHECKING([htslib])

        if test "${ac_cv_libhts}" = "yes" && test "${ac_cv_faidx_h}" = "yes" ; then
            AC_MSG_RESULT([ok])
            HTSLIB_OK=1
            AC_DEFINE(HAVE_HTSLIB, [1],[define if the hts library is available])
        else
            AC_MSG_RESULT([failed])    
            AC_MSG_ERROR([HTSlib headers not found. Please ensure that the htslib headers can be found on the CPPFLAGS env var. Alternatively, you can try the --with-htslib option, which assumes you have the following setup <htslib_dir>/include/htslib/<headers>.])
        fi
    fi
])