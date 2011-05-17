## Process this file with automake to produce Makefile.in

SUBDIRS = RStarTree kdt fes2004 wavewatch

AM_CPPFLAGS = -DGFS_MODULES_DIR=\"$(libdir)/gerris\"

INCLUDES = -I$(top_srcdir)/src -I$(includedir) \
           -DG_LOG_DOMAIN=\"Gfs-modules\" $(GTS_CFLAGS)

if HAS_LIBPROJ
MAP = libmap2D.la libmap3D.la
endif
if BUILD_TIDE
TIDE = libtide2D.la libtide3D.la
endif
if BUILD_STOKES
STOKES = libstokes2D.la libstokes3D.la
endif
if BUILD_WAVEWATCH
WAVEWATCH = libwavewatch2D.la
endif
if BUILD_HYPRE
HYPRE = libhypre2D.la libhypre3D.la
endif
if BUILD_LIS
LIS = liblis2D.la liblis3D.la
endif
if BUILD_AGMG
AGMG = libagmg2D.la libagmg3D.la
endif
if BUILD_ODE
ODE = libode2D.la libode3D.la
endif
if HAS_GSL
BUBBLES = libbubbles2D.la libbubbles3D.la
endif

pkglib_LTLIBRARIES = \
	$(MAP) \
	libterrain2D.la \
	libterrain3D.la \
	libparticulates2D.la \
	libparticulates3D.la \
	libokada2D.la \
	libokada3D.la \
	libdf33D.la \
	libelectrohydro2D.la \
	libelectrohydro3D.la \
	libtopics2D.la \
	libtopics3D.la \
	$(TIDE) \
	$(STOKES) \
	$(WAVEWATCH) \
	$(HYPRE) \
	$(LIS) \
	$(AGMG) \
	$(ODE) \
        $(BUBBLES)

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
	bubbles.mod \
	df3.mod \
	electrohydro.mod \
	okada.mod \
	hypre.mod \
	lis.mod \
	agmg.mod \
	ode.mod \
	topics.mod

BUILT_SOURCES = \
	map.c \
	tide.c \
	terrain.c \
	stokes.c \
	wavewatch.c \
	particulates.c \
	bubbles.c \
	df3.c \
	electrohydro.c \
	okada.c \
	hypre.c \
	lis.c \
	agmg.c \
	ode.c \
	topics.c

AM_LDFLAGS = $(NO_UNDEFINED)\
        -version-info $(LT_CURRENT):$(LT_REVISION):$(LT_AGE)\
	-release $(LT_RELEASE) -export-dynamic -module

libmap3D_la_SOURCES = map.c
libmap3D_la_LIBADD = -lproj
libmap2D_la_SOURCES = map.c
libmap2D_la_CFLAGS = $(AM_CFLAGS) -DFTT_2D=1
libmap2D_la_LIBADD = -lproj

KDTLIBS = -LRStarTree -lcSmRST -Lkdt -lkdt
KDTDEPS = RStarTree/libcSmRST.la kdt/libkdt.la

libterrain3D_la_SOURCES = terrain.c rsurface.c rsurface.h
libterrain3D_la_LIBADD = $(KDTLIBS)
libterrain3D_la_DEPENDENCIES = $(KDTDEPS)

libterrain2D_la_CFLAGS = $(AM_CFLAGS) -DFTT_2D=1
libterrain2D_la_SOURCES = terrain.c rsurface.c rsurface.h
libterrain2D_la_LIBADD = $(KDTLIBS)
libterrain2D_la_DEPENDENCIES = $(KDTDEPS)

libtide3D_la_SOURCES = tide.c
libtide3D_la_CFLAGS = $(AM_CFLAGS) $(GSL_CFLAGS)
libtide3D_la_LIBADD = -Lfes2004 -lfes -lnetcdf $(GSL_LIBS)
libtide2D_la_SOURCES = tide.c
libtide2D_la_CFLAGS = $(AM_CFLAGS) $(GSL_CFLAGS) -DFTT_2D=1
libtide2D_la_LIBADD = -Lfes2004 -lfes -lnetcdf $(GSL_LIBS)

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

libode3D_la_SOURCES = ode.c
libode3D_la_CFLAGS = $(AM_CFLAGS) $(ODE_CFLAGS)
libode3D_la_LIBADD = $(ODE_LIBS)
libode2D_la_SOURCES = ode.c
libode2D_la_CFLAGS = $(AM_CFLAGS) $(ODE_CFLAGS) -DFTT_2D=1
libode2D_la_LIBADD = $(ODE_LIBS)

xyz2rsurface_SOURCES = xyz2rsurface.c rsurface.c rsurface.h
xyz2rsurface_LDADD = $(KDTLIBS)
xyz2rsurface_CFLAGS = $(AM_CFLAGS)
xyz2rsurface_DEPENDENCIES = $(KDTDEPS)

rsurfacequery_SOURCES = rsurfacequery.c rsurface.c rsurface.h
rsurfacequery_LDADD = -LRStarTree -lcSmRST -Lkdt -lkdt
rsurfacequery_CFLAGS = $(AM_CFLAGS)
rsurfacequery_DEPENDENCIES = $(KDTDEPS)

rsurfacedraw_SOURCES = rsurfacedraw.c rsurface.c rsurface.h
rsurfacedraw_LDADD = $(KDTLIBS)
rsurfacedraw_CFLAGS = $(AM_CFLAGS)
rsurfacedraw_DEPENDENCIES = $(KDTDEPS)

rsurfaceinfo_SOURCES = rsurfaceinfo.c rsurface.c rsurface.h
rsurfaceinfo_LDADD = $(KDTLIBS)
rsurfaceinfo_CFLAGS = $(AM_CFLAGS)
rsurfaceinfo_DEPENDENCIES = $(KDTDEPS)

libparticulates3D_la_SOURCES = particulates.c particulates.h
libparticulates3D_la_CFLAGS = $(AM_CFLAGS)
libparticulates2D_la_SOURCES = particulates.c particulates.h
libparticulates2D_la_CFLAGS = $(AM_CFLAGS) -DFTT_2D=1

libbubbles3D_la_SOURCES = bubbles.c
libbubbles3D_la_CFLAGS = $(AM_CFLAGS) $(GSL_CFLAGS)
libbubbles3D_la_LIBADD = $(GSL_LIBS)
libbubbles2D_la_SOURCES = bubbles.c
libbubbles2D_la_CFLAGS = $(AM_CFLAGS) $(GSL_CFLAGS) -DFTT_2D=1
libbubbles2D_la_LIBADD = $(GSL_LIBS)

libhypre3D_la_SOURCES = hypre.c
libhypre3D_la_CFLAGS = $(AM_CFLAGS)
libhypre3D_la_LIBADD = -lHYPRE -lm
libhypre2D_la_SOURCES = hypre.c
libhypre2D_la_CFLAGS = $(AM_CFLAGS) -DFTT_2D=1
libhypre2D_la_LIBADD = -lHYPRE -lm

liblis3D_la_SOURCES = lis.c
liblis3D_la_CFLAGS = $(AM_CFLAGS)
liblis3D_la_LIBADD = -llis -lgfortran -lm
liblis2D_la_SOURCES = lis.c
liblis2D_la_CFLAGS = $(AM_CFLAGS) -DFTT_2D=1
liblis2D_la_LIBADD = -llis -lgfortran -lm

libagmg3D_la_SOURCES = agmg.c dagmg.f90
libagmg3D_la_CFLAGS = $(AM_CFLAGS)
libagmg3D_la_LIBADD = -llapack -lblas -lgfortran -lm
libagmg2D_la_SOURCES = agmg.c dagmg.f90
libagmg2D_la_CFLAGS = $(AM_CFLAGS) -DFTT_2D=1
libagmg2D_la_LIBADD = -llapack -lblas -lgfortran -lm

libokada3D_la_SOURCES = okada.c
libokada3D_la_CFLAGS = $(AM_CFLAGS)
libokada2D_la_SOURCES = okada.c
libokada2D_la_CFLAGS = $(AM_CFLAGS) -DFTT_2D=1

libdf33D_la_SOURCES = df3.c
libdf33D_la_CFLAGS = $(AM_CFLAGS)

libelectrohydro3D_la_SOURCES = electrohydro.c
libelectrohydro3D_la_CFLAGS = $(AM_CFLAGS)
libelectrohydro2D_la_SOURCES = electrohydro.c
libelectrohydro2D_la_CFLAGS = $(AM_CFLAGS) -DFTT_2D=1

libtopics3D_la_SOURCES = topics.c
libtopics3D_la_CFLAGS = $(AM_CFLAGS)
libtopics2D_la_SOURCES = topics.c
libtopics2D_la_CFLAGS = $(AM_CFLAGS) -DFTT_2D=1

if HAVE_MODULES
.mod.c:
	@echo "/* $@" > $@
	@echo " * This is a generated file.  Please modify \`$<'" >> $@
	@echo " */" >> $@
	@echo "#line 1 \""$<"\"" >> $@
	@cat $< >> $@
else
.mod.c:
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
