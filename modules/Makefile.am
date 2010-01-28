## Process this file with automake to produce Makefile.in

SUBDIRS = RStarTree fes2004 wavewatch

AM_CPPFLAGS = -DGFS_MODULES_DIR=\"$(libdir)/gerris\"

INCLUDES = -I$(top_srcdir)/src -I$(includedir) \
           -DG_LOG_DOMAIN=\"Gfs-modules\" $(GTS_CFLAGS)

if HAS_LIBPROJ
MAP = libmap2D.la libmap3D.la libmap2D3.la
endif
if BUILD_TIDE
TIDE = libtide2D.la libtide3D.la libtide2D3.la
endif
if BUILD_STOKES
STOKES = libstokes2D.la libstokes3D.la
endif
if BUILD_WAVEWATCH
WAVEWATCH = libwavewatch2D.la
endif

pkglib_LTLIBRARIES = \
	$(MAP) \
	libterrain2D.la \
	libterrain3D.la \
	libterrain2D3.la \
	libparticulates2D.la \
	libparticulates3D.la \
	libparticulates2D3.la \
	libokada2D.la \
	libokada3D.la \
	$(TIDE) \
	$(STOKES) \
	$(WAVEWATCH)

bin_PROGRAMS = \
	xyz2rsurface \
	rsurfacequery \
	rsurfacedraw \
	rsurfaceinfo

EXTRA_DIST = \
	map.mod \
	tide.mod \
	terrain.mod \
	stokes.mod \
	wavewatch.mod \
	particulates.mod \
	okada.mod

BUILT_SOURCES = \
	map.c \
	tide.c \
	terrain.c \
	stokes.c \
	wavewatch.c \
	particulates.c \
	okada.c

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
libterrain3D_la_DEPENDENCIES = RStarTree/libcSmRST.la

libterrain2D_la_CFLAGS = $(AM_CFLAGS) -DFTT_2D=1
libterrain2D_la_SOURCES = terrain.c rsurface.c rsurface.h
libterrain2D_la_LIBADD = -LRStarTree -lcSmRST
libterrain2D_la_DEPENDENCIES = RStarTree/libcSmRST.la

libterrain2D3_la_SOURCES = terrain.c rsurface.c rsurface.h
libterrain2D3_la_CFLAGS = $(AM_CFLAGS) -DFTT_2D3=1
libterrain2D3_la_LIBADD = -LRStarTree -lcSmRST
libterrain2D3_la_DEPENDENCIES = RStarTree/libcSmRST.la

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

libwavewatch2D_la_SOURCES = wavewatch.c
libwavewatch2D_la_LIBADD = -Lwavewatch -lwavewatch -L/usr/lib/gcc/i486-linux-gnu/4.2 -lgfortran
libwavewatch2D_la_CFLAGS = $(AM_CFLAGS) -DFTT_2D=1 \
  -DW3SRCE=`nm wavewatch/libwavewatch.a | grep w3srce | awk '{ if ($$2 == "T") print $$3; }'` \
  -DGFSW3INIT=`nm wavewatch/libwavewatch.a | grep gfsw3_init | awk '{ if ($$2 == "T") print $$3; }'`
libwavewatch2D_la_DEPENDENCIES = wavewatch/libwavewatch.a

xyz2rsurface_SOURCES = xyz2rsurface.c rsurface.c rsurface.h
xyz2rsurface_LDADD = -LRStarTree -lcSmRST
xyz2rsurface_CFLAGS = $(AM_CFLAGS)
xyz2rsurface_DEPENDENCIES = RStarTree/libcSmRST.la

rsurfacequery_SOURCES = rsurfacequery.c rsurface.c rsurface.h
rsurfacequery_LDADD = -LRStarTree -lcSmRST
rsurfacequery_CFLAGS = $(AM_CFLAGS)
rsurfacequery_DEPENDENCIES = RStarTree/libcSmRST.la

rsurfacedraw_SOURCES = rsurfacedraw.c rsurface.c rsurface.h
rsurfacedraw_LDADD = -LRStarTree -lcSmRST
rsurfacedraw_CFLAGS = $(AM_CFLAGS)
rsurfacedraw_DEPENDENCIES = RStarTree/libcSmRST.la

rsurfaceinfo_SOURCES = rsurfaceinfo.c rsurface.c rsurface.h
rsurfaceinfo_LDADD = -LRStarTree -lcSmRST
rsurfaceinfo_CFLAGS = $(AM_CFLAGS)
rsurfaceinfo_DEPENDENCIES = RStarTree/libcSmRST.la

libparticulates3D_la_SOURCES = particulates.c
libparticulates3D_la_CFLAGS = $(AM_CFLAGS) $(GSL_CFLAGS)
libparticulates2D_la_SOURCES = particulates.c
libparticulates2D_la_CFLAGS = $(AM_CFLAGS) $(GSL_CFLAGS) -DFTT_2D=1
libparticulates2D3_la_SOURCES = particulates.c
libparticulates2D3_la_CFLAGS = $(AM_CFLAGS) $(GSL_CFLAGS) -DFTT_2D3=1

libokada3D_la_SOURCES = okada.c
libokada3D_la_CFLAGS = $(AM_CFLAGS)
libokada2D_la_SOURCES = okada.c
libokada2D_la_CFLAGS = $(AM_CFLAGS) -DFTT_2D=1

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
