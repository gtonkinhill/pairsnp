AC_INIT([amhello], [1.0], [bug-automake@gnu.org])
AM_INIT_AUTOMAKE([-Wall -Werror foreign])
AC_PROG_CC
AC_PROG_CXX

# Armadillo 
AC_LANG_PUSH([C++])
AC_CHECK_HEADERS([armadillo],[], [AC_MSG_ERROR("Please install armadillo")])
AC_LANG_PUSH([C++])

AC_CONFIG_HEADERS([config.h])
AC_CONFIG_FILES([
 Makefile
 src/Makefile
])
AC_OUTPUT
