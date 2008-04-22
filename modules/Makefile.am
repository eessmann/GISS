## Process this file with automake to produce Makefile.in

SUBDIRS = RStarTree

INCLUDES = -I$(top_srcdir)/src -I$(includedir) \
           -DG_LOG_DOMAIN=\"Gfs-modules\" $(GTS_CFLAGS)

if HAS_LIBPROJ
MAP = libmap2D.la libmap3D.la libmap2D3.la
endif

pkglib_LTLIBRARIES = \
	$(MAP) \
	libterrain2D.la \
	libterrain3D.la \
	libterrain2D3.la \
	libtide2D.la \
	libtide3D.la \
	libtide2D3.la

bin_PROGRAMS = \
	xyz2rsurface \
	rsurfacequery

EXTRA_DIST = \
	map.mod \
	tide.mod \
	terrain.mod

BUILT_SOURCES = \
	map.c \
	tide.c \
	terrain.c

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

libterrain3D_la_SOURCES = terrain.c rsurface.c
libterrain3D_la_LIBADD = -LRStarTree -lcSmRST

libterrain2D_la_CFLAGS = $(AM_CFLAGS) -DFTT_2D=1
libterrain2D_la_SOURCES = terrain.c rsurface.c
libterrain2D_la_LIBADD = -LRStarTree -lcSmRST

libterrain2D3_la_SOURCES = terrain.c rsurface.c
libterrain2D3_la_CFLAGS = $(AM_CFLAGS) -DFTT_2D3=1
libterrain2D3_la_LIBADD = -LRStarTree -lcSmRST

libtide3D_la_SOURCES = tide.c
libtide3D_la_LIBADD = -L/home/popinet/local/lib -lfes -lnetcdf -lgsl -lgslcblas -lm -lpthread
libtide2D_la_SOURCES = tide.c
libtide2D_la_CFLAGS = $(AM_CFLAGS) -DFTT_2D=1
libtide2D_la_LIBADD = -L/home/popinet/local/lib -lfes -lnetcdf -lgsl -lgslcblas -lm -lpthread
libtide2D3_la_SOURCES = tide.c
libtide2D3_la_CFLAGS = $(AM_CFLAGS) -DFTT_2D3=1
libtide2D3_la_LIBADD = -L/home/popinet/local/lib -lfes -lnetcdf -lgsl -lgslcblas -lm -lpthread

xyz2rsurface_SOURCES = xyz2rsurface.c rsurface.c
xyz2rsurface_LDADD = -LRStarTree -lcSmRST

rsurfacequery_SOURCES = rsurfacequery.c rsurface.c
rsurfacequery_LDADD = -LRStarTree -lcSmRST

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
