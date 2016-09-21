## Process this file with automake to produce Makefile.in

SUBDIRS = desktop

EXTRA_DIST = \
	mangle \
	glade/interface.c \
	glade/gfsview.glade \
	glade/gfsview.gladep \
	gfsview

INCLUDES = \
	-DPACKAGE_DATA_DIR=\""$(datadir)/gfsview"\" \
	-DPACKAGE_LOCALE_DIR=\""$(prefix)/$(DATADIRNAME)/locale"\" \
	-I$(top_srcdir)

bin_PROGRAMS = gfsview2D gfsview3D

bin_SCRIPTS = gfsview

BUILT_SOURCES = glade/mangled_interface.c

CLEANFILES = $(BUILT_SOURCES)

SRC = \
	glade/mangled_interface.c glade/interface.h \
	glade/callbacks.c glade/callbacks.h \
	glade/support.c glade/support.h \
	$(top_srcdir)/gl2ps/gl2ps.c $(top_srcdir)/gl2ps/gl2ps.h \
        main.c  \
	gfkgl.c gfkgl.h

gfsview2D_SOURCES = $(SRC) gfkgl2D.h
gfsview2D_CFLAGS = @SN_CFLAGS@ @GTK_CFLAGS@ @GERRIS2D_CFLAGS@
gfsview2D_LDADD = -L$(top_builddir)/gl2ps -lgl2ps \
                  -L$(top_builddir)/gl -lgfsgl2D \
	          @SN_LIBS@ @GTK_LIBS@ @GERRIS2D_LIBS@
gfsview2D_DEPENDENCIES = $(top_builddir)/gl2ps/libgl2ps.la $(top_builddir)/gl/libgfsgl2D.la

gfsview3D_SOURCES = $(SRC) gfkgl3D.h
gfsview3D_CFLAGS = @SN_CFLAGS@ @GTK_CFLAGS@ @GERRIS3D_CFLAGS@
gfsview3D_LDADD = -L$(top_builddir)/gl2ps -lgl2ps \
                  -L$(top_builddir)/gl -lgfsgl3D \
                  @SN_LIBS@ @GTK_LIBS@ @GERRIS3D_LIBS@
gfsview3D_DEPENDENCIES = $(top_builddir)/gl2ps/libgl2ps.la $(top_builddir)/gl/libgfsgl3D.la

glade/mangled_interface.c: $(srcdir)/glade/interface.c $(srcdir)/glade/interface.h
	sh $(srcdir)/mangle < $(srcdir)/glade/interface.c > $(srcdir)/glade/mangled_interface.c
