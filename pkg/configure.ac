## 
## This file is part of fftw
##
## Process this file with autoconf to produce a configure script.
##

AC_PREREQ(2.61)

## ---[ VERSION ]-------------------------------------------------------------

AC_INIT()
AC_CONFIG_AUX_DIR(src/build/autoconf)

## pick up compiler as will be used by R CMD INSTALL/SHLIB
if test -n "${R_HOME}"; then
  CC=`${R_HOME}/bin/R CMD config CC`
  CPP=`${R_HOME}/bin/R CMD config CPP`
  CFLAGS=`${R_HOME}/bin/R CMD config CFLAGS`
  CXX=`${R_HOME}/bin/R CMD config CXX`
  CXXFLAGS=`${R_HOME}/bin/R CMD config CXXFLAGS`
fi

AC_PROG_CPP
AC_PROG_CC

if test `uname` = "Darwin" ; then
  darwin="yes"
else
  darwin="no"
fi

## --- FFTW 3 ----------------------------------------------------------------
AC_ARG_WITH(fftw3-prefix,
[  --with-ffw3-prefix=DIR     specify location of FFTW3 headers / libs (defaults to "/usr/local")],
[fftw3_prefix="${withval}"],
[fftw3_prefix="/usr/local/"]
)

AC_ARG_WITH(fftw3-includes,
[  --with-fftw3-includes=DIR  specify location of FFTW 3 headers],
[CPPFLAGS="${CPPFLAGS} -I${withval}"],
[CPPFLAGS="${CPPFLAGS} -I${fftw3_prefix}/include"]
)

if test $darwin != yes; then
  AC_CHECK_HEADERS(fftw3.h)
  if test "x$ac_cv_header_fftw3_h" = xno; then
    AC_ERROR([missing required header fftw3.h])
  fi
fi

AC_ARG_WITH(fftw3-libs,
[  --with-fftw3-libs=DIR      specify location of FFTW 3 libs],
[LIBS="${LIBS} -L${withval}"],
[LIBS="${LIBS} -L${fftw3_prefix}/lib"]
)

AC_ARG_WITH(fftw3-libname,
[  --with-fftw3-libname=NAME  specify Library name (defaults to "fftw3")], 
[lfftw3=${withval}], 
[lfftw3=fftw3]
)

AC_CHECK_LIB($lfftw3, fftw_import_system_wisdom)
this=`eval echo '${'$as_ac_Lib'}'`
if test "x$this" != xyes; then
  AC_ERROR([missing required library ${lfftw3}])
fi

## --- Output ----------------------------------------------------------------

AC_SUBST(CPPFLAGS)
AC_SUBST(CXXFLAGS)
AC_SUBST(LIBS)
AC_SUBST(AGLCPPFLAGS)
AC_SUBST(AGLLIBS)
AC_SUBST(BUILD_AGL_TRUE)

AC_OUTPUT(src/Makevars)
