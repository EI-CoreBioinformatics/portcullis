# ===========================================================================
#       dlib check
# ===========================================================================
#
# SYNOPSIS
#
#   AX_DLIB([action-if-found], [action-if-not-found])
#
# DESCRIPTION
#
#   This macro searches for an installed dlib library. If nothing was
#   specified when calling configure, it searches first in /usr/local and
#   then in /usr, /opt/local and /sw. If the --with-dlib=DIR is specified,
#   it will try to find it in DIR/dlib/clustering.h. If
#   --without-dlib is specified, the library is not searched at all.
#
#   If the header file (dlib.h) is not found,
#   shell commands 'action-if-not-found' is run. If 'action-if-not-found' is
#   not specified, the configuration exits on error, asking for a valid dlib
#   installation directory or --without-dlib.
#
#   If both header file and library are found, shell commands
#   'action-if-found' is run. If 'action-if-found' is not specified, the
#   default action appends '-isystem ${DLIB_HOME}' to CPPFLAGS, and calls
#   AC_DEFINE(HAVE_DLIB). You should use autoheader to include a definition
#   for this symbol in a config.h file. Sample usage in a C/C++ source is as
#   follows:
#
#     #ifdef HAVE_DLIB
#     #include <dlib/clustering.h>
#     #endif /* HAVE_DLIB */
#
# LICENSE
#


#serial 14

AU_ALIAS([DLIB], [AX_DLIB])
AC_DEFUN([AX_DLIB],
#
# Handle user hints
#
[AC_MSG_CHECKING(if dlib is wanted)
dlib_places="/usr/local/include /usr/include /opt/local/include"
AC_ARG_WITH([dlib],
[  --with-dlib=DIR        root directory path of dlib installation @<:@defaults to
                          /usr/local or /usr if not found in /usr/local@:>@],
[if test "$withval" != no ; then
  AC_MSG_RESULT(yes)
  if test -d "$withval"
  then
    dlib_places="$withval $dlib_places"
  else
    AC_MSG_WARN([Sorry, $withval does not exist, checking usual places])
  fi
else
  dlib_places=
  AC_MSG_RESULT(no)
fi],
[AC_MSG_RESULT(yes)])

#
# Locate dlib, if wanted
#
if test -n "${dlib_places}"
then
	# check the user supplied or any other more or less 'standard' place:
	#   Most UNIX systems      : /usr/local and /usr
	#   MacPorts / Fink on OSX : /opt/local respectively /sw
	for DLIB_HOME in ${dlib_places} ; do
          if test -f "${DLIB_HOME}/dlib/clustering.h"; then break; fi
	  DLIB_HOME=""
	done

  DLIB_OLD_CPPFLAGS=$CPPFLAGS
  if test -n "${DLIB_HOME}"; then
        CPPFLAGS="$CPPFLAGS -isystem ${DLIB_HOME}"
  fi
  
  AC_LANG_SAVE
  AC_LANG([C++])
  AC_CHECK_HEADER([dlib/clustering.h], [dlib_h=yes], [dlib_h=no])
  AC_LANG_RESTORE

  if test "$dlib_h" = "yes"
  then
    #
    # If header were found, action-if-found
    #
    m4_ifblank([$1],[
                
                DLIB_CPPFLAGS="-isystem ${DLIB_HOME}"
                AC_SUBST([DLIB_CPPFLAGS])
                CPPFLAGS="$DLIB_OLD_CPPFLAGS"
                
                AC_DEFINE([HAVE_DLIB], [1],
                          [Define to 1 if you have dlib])
               ],[
                # Restore variables
                CPPFLAGS="$DLIB_OLD_CPPFLAGS"
                $1
               ])
  else
    #
    # If either header was not found, action-if-not-found
    #
    m4_default([$2],[
                AC_MSG_ERROR([either specify a valid dlib installation with --with-dlib=DIR or disable dlib usage with --without-dlib])
                ])
  fi
fi
])