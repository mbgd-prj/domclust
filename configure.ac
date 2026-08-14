dnl Process this file with autoconf to produce a configure script.
AC_PREREQ([2.69])
AC_INIT([domclust], [1.2], [uchiyama@nibb.ac.jp])
AC_CONFIG_SRCDIR([src/domclust.c])

dnl Preserve the historical -O3 default (autoconf's AC_PROG_CC would
dnl otherwise default to "-g -O2" when CFLAGS is unset).
: ${CFLAGS="-O3"}

AC_PROG_CC
AC_PROG_INSTALL

dnl Detect OpenMP support; OPENMP_CFLAGS is left empty (build proceeds
dnl without OpenMP) when the compiler/environment doesn't support it.
AC_OPENMP

AC_ARG_ENABLE([dev],
  [AS_HELP_STRING([--enable-dev],
    [build development version with gene-neighbor / domain-cut / new-alignment support (equivalent to the old "make all2")])],
  [enable_dev=$enableval], [enable_dev=no])

AC_ARG_ENABLE([large],
  [AS_HELP_STRING([--enable-large],
    [enable large-scale mode with nucleotide reverse-match direction tracking (equivalent to the old "make all3")])],
  [enable_large=$enableval], [enable_large=no])

FEATURE_CFLAGS=
EXTRA_OBJS=

AS_IF([test "x$enable_dev" = xyes], [
  FEATURE_CFLAGS="$FEATURE_CFLAGS -DWITH_NEIGHBOR -DWITH_DOMCUT -DWITH_CALNEWALI"
  EXTRA_OBJS="$EXTRA_OBJS genedp.o domcut.o calnewali.o mergehit.o poslist.o hitlist.o"
])

AS_IF([test "x$enable_large" = xyes], [
  FEATURE_CFLAGS="$FEATURE_CFLAGS -DLARGE"
])

AC_SUBST([FEATURE_CFLAGS])
AC_SUBST([EXTRA_OBJS])
AC_SUBST([OPENMP_CFLAGS])

AC_CHECK_LIB([m], [sqrt])

AC_CONFIG_FILES([Makefile])
AC_OUTPUT
