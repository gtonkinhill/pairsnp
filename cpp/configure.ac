AC_INIT([pairsnp], [0.0.1], [g.tonkinhill@gmail.com])
AM_INIT_AUTOMAKE([-Wall -Werror foreign])
AC_PROG_CC

CXXFLAGS="$CXXFLAGS -std=c++1y -g -O3"
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
