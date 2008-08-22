## Process this file with automake to produce Makefile.in

SUBDIRS = RStarTree fes2004

AM_CPPFLAGS = -DGFS_MODULES_DIR=\"$(libdir)/gerris\"

INCLUDES = -I$(top_srcdir)/src -I$(includedir) \
           -DG_LOG_DOMAIN=\"Gfs-modules\" $(GTS_CFLAGS)

if HAS_LIBPROJ
MAP = libmap2D.la libmap3D.la libmap2D3.la
endif
if BUILD_TIDE
TIDE = libtide2D.la libtide3D.la libtide2D3.la
endif

pkglib_LTLIBRARIES = \
	$(MAP) \
	libterrain2D.la \
	libterrain3D.la \
	libterrain2D3.la \
	$(TIDE) \
	libstokes2D.la \
	libstokes3D.la

bin_PROGRAMS = \
	xyz2rsurface \
	rsurfacequery \
	rsurfacedraw

EXTRA_DIST = \
	map.mod \
	tide.mod \
	terrain.mod \
	stokes.mod

BUILT_SOURCES = \
	map.c \
	tide.c \
	terrain.c \
	stokes.c

AM_LDFLAGS = $(NO_UNDEFINED)\
        -version-info $(LT_CURRENT):$(LT_REVISION):$(LT_AGE)\
	-release $(LT_RELEASE) -export-dynamic

libmap3D_la_SOURCES = map.c
libmap3D_la_LIBADD = -lproj
libmap2D_la_SOURCES = map.c
libmap2D_la_CFLAGS = $(AM_CFLAGS) -DFTT_2D=1
libmap2D_la_LIBADD = -lproj
libmap2D3_la_SOURCES = map.c
libmap2D3_la_CFLAGS = $(AM_CFLAGS) -DFTT_2D3=1
libmap2D3_la_LIBADD = -lproj

libterrain3D_la_SOURCES = terrain.c rsurface.c rsurface.h
libterrain3D_la_LIBADD = -LRStarTree -lcSmRST

libterrain2D_la_CFLAGS = $(AM_CFLAGS) -DFTT_2D=1
libterrain2D_la_SOURCES = terrain.c rsurface.c rsurface.h
libterrain2D_la_LIBADD = -LRStarTree -lcSmRST

libterrain2D3_la_SOURCES = terrain.c rsurface.c rsurface.h
libterrain2D3_la_CFLAGS = $(AM_CFLAGS) -DFTT_2D3=1
libterrain2D3_la_LIBADD = -LRStarTree -lcSmRST

libtide3D_la_SOURCES = tide.c
libtide3D_la_CFLAGS = $(AM_CFLAGS) $(GSL_CFLAGS)
libtide3D_la_LIBADD = -Lfes2004 -lfes -lnetcdf $(GSL_LIBS)
libtide2D_la_SOURCES = tide.c
libtide2D_la_CFLAGS = $(AM_CFLAGS) $(GSL_CFLAGS) -DFTT_2D=1
libtide2D_la_LIBADD = -Lfes2004 -lfes -lnetcdf $(GSL_LIBS)
libtide2D3_la_SOURCES = tide.c
libtide2D3_la_CFLAGS = $(AM_CFLAGS) $(GSL_CFLAGS) -DFTT_2D3=1
libtide2D3_la_LIBADD = -Lfes2004 -lfes -lnetcdf $(GSL_LIBS)

libstokes3D_la_SOURCES = stokes.c CW263.f
libstokes3D_la_CFLAGS = $(AM_CFLAGS)
libstokes2D_la_SOURCES = stokes.c CW263.f
libstokes2D_la_CFLAGS = $(AM_CFLAGS) -DFTT_2D=1

xyz2rsurface_SOURCES = xyz2rsurface.c rsurface.c rsurface.h
xyz2rsurface_LDADD = -LRStarTree -lcSmRST
xyz2rsurface_CFLAGS = $(AM_CFLAGS)

rsurfacequery_SOURCES = rsurfacequery.c rsurface.c rsurface.h
rsurfacequery_LDADD = -LRStarTree -lcSmRST
rsurfacequery_CFLAGS = $(AM_CFLAGS)

rsurfacedraw_SOURCES = rsurfacedraw.c rsurface.c rsurface.h
rsurfacedraw_LDADD = -LRStarTree -lcSmRST
rsurfacedraw_CFLAGS = $(AM_CFLAGS)

if HAVE_MODULES
%.c : %.mod
	@echo "/* $@" > $@
	@echo " * This is a generated file.  Please modify \`$<'" >> $@
	@echo " */" >> $@
	@echo "#line 1 \""$<"\"" >> $@
	@cat $< >> $@
else
%.c : %.mod
	@echo "/* $@" > $@
	@echo " * This is a generated file.  Please modify \`$<'" >> $@
	@echo " */" >> $@
	@echo "#line 1 \""$<"\"" >> $@
	@replace="s/g_module_check_init/gfs_init_`basename $< .mod`/g";\
	 sed $$replace < $< >> $@
endif

clean-generic:
	$(RM) $(BUILT_SOURCES)

dist-hook:
	for file in $(BUILT_SOURCES); do $(RM) $(distdir)/$$file; done
