dnl Process this file with autoconf to produce a configure script.

AC_INIT(configure.ac)
AC_CONFIG_MACRO_DIR([m4])


AC_CANONICAL_HOST

AM_INIT_AUTOMAKE(gfsview, 0.5.0)
AM_CONFIG_HEADER(config.h)
AM_MAINTAINER_MODE

AC_ARG_WITH(mpicc,
              [  --with-mpicc=path       pick a specific path to mpicc],
              [with_mpicc=$withval],
              [with_mpicc=''])

# check if MPI support is disabled
AC_ARG_ENABLE(mpi,
[  --disable-mpi           do not compile MPI support],
[ case "${enableval}" in
	yes) ;;
	*) with_mpicc=no ;;
  esac])

if test "x$with_mpicc" != "xno" ; then
  if test "x$with_mpicc" != "x" ; then
     if test -x "$with_mpicc" ; then
       use_mpicc=yes
       CC="$with_mpicc"
     fi
  else
     AC_CHECK_PROG(use_mpicc, mpicc, yes)
     if test "x$use_mpicc" = "xyes" ; then
       CC=mpicc
     fi
  fi
  if test "x$use_mpicc" = "xyes" ; then
     AC_DEFINE(HAVE_MPI, 1, [Define to 1 if you have MPI])
  fi
fi

AC_ISC_POSIX
AC_PROG_CC
AM_PROG_CC_STDC
AM_PROG_CC_C_O
AC_HEADER_STDC

dnl Initialize libtool
AC_LIBTOOL_WIN32_DLL
AM_PROG_LIBTOOL
AM_PROG_CC_C_O

if test x$GCC = xyes ; then
  CFLAGS="$CFLAGS -Wall -Werror-implicit-function-declaration -Wmissing-prototypes -Wmissing-declarations -pipe -std=c99 -D_XOPEN_SOURCE=500"
fi

GERRIS_VERSION="1.2.0"
PKG_CHECK_MODULES(GERRIS2D,  [gerris2D  >= $GERRIS_VERSION])
PKG_CHECK_MODULES(GERRIS3D,  [gerris3D  >= $GERRIS_VERSION])
PKG_CHECK_MODULES(OSMESA, [ osmesa ], [ have_osmesa=yes ], [ have_osmesa=no ])

dnl alternative check for systems without pkg-config for OSMesa
if test "x$have_osmesa" = "xno"; then
   AC_CHECK_LIB(OSMesa, OSMesaMakeCurrent, [
      have_osmesa=yes
      OSMESA_LIBS="-lOSMesa -lGL -lGLU "
   ], [
      FLAGS=$CFLAGS
      CFLAGS="$CFLAGS -L/usr/X11R6/lib"
      AC_CHECK_LIB(OSMesa, OSMesaCreateContext, [
         have_osmesa=yes
         OSMESA_LIBS="-lOSMesa -lGL -lGLU "
         LDFLAGS="$LDFLAGS -L/usr/X11R6/lib"
      ], [
         AC_MSG_WARN([OSMesa not found. Batch mode will not be available.])
         have_osmesa=no
         OSMESA_LIBS=""
      ])
      CFLAGS=$FLAGS
   ])
fi

AM_CONDITIONAL(HAVE_OSMESA, test "x$have_osmesa" = "xyes")

dnl Fix for MacOSX and cygwin libtool
case "$build" in
  *-apple-darwin*) 
  		   NO_UNDEFINED="" 
		   CFLAGS="$CFLAGS -D_DARWIN_C_SOURCE"
		   ;;
  *-pc-cygwin*)    NO_UNDEFINED="" ;;
  *)               NO_UNDEFINED="-no-undefined" ;;
esac
AC_SUBST(NO_UNDEFINED)

GTK_PACKAGES="gthread-2.0 gtk+-2.0 >= 2.4.0 gtkglext-1.0 >= 1.0.6 libstartup-notification-1.0"
if test "x$have_osmesa" = "xyes"; then
        PKG_CHECK_MODULES(GTK, [ $GTK_PACKAGES ],
	[
	  have_gtk=yes
        ],[
	  AC_MSG_WARN([GTK dependencies not met.])
	  have_gtk=no
	])
else
	PKG_CHECK_MODULES(GTK, [ $GTK_PACKAGES ])
	have_gtk=yes
fi

AM_CONDITIONAL(HAVE_GTK, test "x$have_gtk" = "xyes")

PKG_CHECK_MODULES(FTGL, [ ftgl >= 2.1.3 ],
	[
	  AC_DEFINE(HAVE_FTGL, 1, [Define to 1 if you have FTGL])
        ],[
	  AC_MSG_WARN([FTGL dependencies not met.])
	])

AC_SUBST(OSMESA_CFLAGS)
AC_SUBST(OSMESA_LIBS)
AC_SUBST(GTK_CFLAGS)
AC_SUBST(GTK_LIBS)
AC_SUBST(FTGL_CFLAGS)
AC_SUBST(FTGL_LIBS)
AC_SUBST(SN_CFLAGS)
AC_SUBST(SN_LIBS)
AC_SUBST(GERRIS2D_CFLAGS)
AC_SUBST(GERRIS2D_LIBS)
AC_SUBST(GERRIS3D_CFLAGS)
AC_SUBST(GERRIS3D_LIBS)

# check whether the code is under darcs

AM_CONDITIONAL(DARCS_CONTROLLED, test -d _darcs)

AC_OUTPUT([
Makefile
view/Makefile
view/desktop/Makefile
gl2ps/Makefile
gl/Makefile
batch/Makefile
debian/Makefile
])

