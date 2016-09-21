## Process this file with automake to produce Makefile.in

INCLUDES = \
	-DPACKAGE_DATA_DIR=\""$(datadir)/gfsview"\" \
	-DPACKAGE_LOCALE_DIR=\""$(prefix)/$(DATADIRNAME)/locale"\" \
        -DG_LOG_DOMAIN=\"GfsGl\" \
	-I$(top_srcdir)

lib_LTLIBRARIES = libgfsgl2D.la libgfsgl3D.la

SRC =   gfsgl.c gfsgl.h \
	trackball.c trackball.h

SRC3D = $(SRC) gfsgl3D.h \
        iso.c

libgfsgl2D_la_CFLAGS = $(AM_CFLAGS) @GERRIS2D_CFLAGS@ @FTGL_CFLAGS@
libgfsgl2D_la_LDFLAGS = $(NO_UNDEFINED) -export-dynamic
libgfsgl2D_la_LIBADD = @GERRIS2D_LIBS@ @FTGL_LIBS@ -lGLU
libgfsgl2D_la_SOURCES = $(SRC) gfsgl2D.h


libgfsgl3D_la_CFLAGS = $(AM_CFLAGS) @GERRIS3D_CFLAGS@ @FTGL_CFLAGS@
libgfsgl3D_la_LDFLAGS = $(NO_UNDEFINED) -export-dynamic
libgfsgl3D_la_LIBADD = @GERRIS3D_LIBS@ @FTGL_LIBS@ -lGLU
libgfsgl3D_la_SOURCES = $(SRC3D)
