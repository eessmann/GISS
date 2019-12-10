## Process this file with automake to produce Makefile.in

INCLUDES = \
	-DPACKAGE_DATA_DIR=\""$(datadir)/gfsview"\" \
	-DPACKAGE_LOCALE_DIR=\""$(prefix)/$(DATADIRNAME)/locale"\" \
        -DG_LOG_DOMAIN=\"GfsView\" \
	-I$(top_srcdir)

noinst_LTLIBRARIES = librender2D.la librender3D.la

librender2D_la_SOURCES = render.c render.h
librender2D_la_CFLAGS = @GERRIS2D_CFLAGS@ @OSMESA_CFLAGS@
librender2D_la_LIBADD = -L$(top_builddir)/gl2ps -lgl2ps \
                        -L$(top_builddir)/gl -lgfsgl2D \
                        @OSMESA_LIBS@
librender2D_la_DEPENDENCIES = $(top_builddir)/gl2ps/libgl2ps.la $(top_builddir)/gl/libgfsgl2D.la

librender3D_la_SOURCES = render.c render.h
librender3D_la_CFLAGS = @GERRIS3D_CFLAGS@ @OSMESA_CFLAGS@
librender3D_la_LIBADD = -L$(top_builddir)/gl2ps -lgl2ps \
                        -L$(top_builddir)/gl -lgfsgl3D \
                        @OSMESA_LIBS@
librender3D_la_DEPENDENCIES = $(top_builddir)/gl2ps/libgl2ps.la $(top_builddir)/gl/libgfsgl3D.la

bin_PROGRAMS = gfsview-batch2D gfsview-batch3D

gfsview_batch2D_SOURCES = gfsview-batch.c
gfsview_batch2D_CFLAGS = @GERRIS2D_CFLAGS@
gfsview_batch2D_LDADD = -lrender2D
gfsview_batch2D_DEPENDENCIES = librender2D.la

gfsview_batch3D_SOURCES = gfsview-batch.c
gfsview_batch3D_CFLAGS = @GERRIS3D_CFLAGS@
gfsview_batch3D_LDADD = -lrender3D
gfsview_batch3D_DEPENDENCIES = librender3D.la

pkglib_LTLIBRARIES = \
	libgfsview2D.la \
	libgfsview3D.la

pkglibdir = $(libdir)/gerris

libgfsview2D_la_SOURCES = module.c
libgfsview2D_la_CFLAGS = @GERRIS2D_CFLAGS@
libgfsview2D_la_LDFLAGS = $(NO_UNDEFINED) -release 0.0.1 -export-dynamic -module
libgfsview2D_la_LIBADD = -lrender2D
libgfsview2D_la_DEPENDENCIES = librender2D.la

libgfsview3D_la_SOURCES = module.c
libgfsview3D_la_CFLAGS = @GERRIS3D_CFLAGS@
libgfsview3D_la_LDFLAGS = $(NO_UNDEFINED) -release 0.0.1 -export-dynamic -module
libgfsview3D_la_LIBADD = -lrender3D
libgfsview3D_la_DEPENDENCIES = librender3D.la
