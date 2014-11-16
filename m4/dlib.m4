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

AC_DEFUN([DLIB],
[
    AC_ARG_WITH([dlib],
       [AS_HELP_STRING([--with-dlib=PREFIX],
		[Use this if your dlib installation prefix was set to a non-standard location.])],
       [DLIB_PATH=$withval], [DLIB_PATH=0])

    # Update vars for dlib if required (cover all bases for the header files)
    if test ${DLIB_PATH} != 0; then
   	DLIB_CPPFLAGS="-I${DLIB_PATH} -I${DLIB_PATH}/include"

        AC_SUBST([DLIB_CPPFLAGS])
        
        OLD_CPPFLAGS=${CPPFLAGS}
        CPPFLAGS="${AM_CPPFLAGS} ${CPPFLAGS} ${DLIB_CPPFLAGS}"
    fi

    # Check header exists
    AC_CHECK_HEADER([dlib/svm.h], [], [
        AC_MSG_ERROR([dlib not found.  Please ensure that the dlib directory can be found on the CPPFLAGS env var. Alternatively, you can try the --with-dlib option.])
    ])
    
    # Restore the previous environment variables if required
    if test ${DLIB_PATH} != 0; then
        CPPFLAGS=${OLD_CPPFLAGS}
    fi
])