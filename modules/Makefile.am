## Process this file with automake to produce Makefile.in

INCLUDES = -I$(top_srcdir)/src -I$(includedir) \
           -DG_LOG_DOMAIN=\"Gfs-modules\" $(GTS_CFLAGS)

if HAS_LIBPROJ
MAP = libmap2D.la libmap3D.la libmap2D3.la
endif

pkglib_LTLIBRARIES = \
	$(MAP) \
	libtide2D.la \
	libtide3D.la \
	libtide2D3.la

EXTRA_DIST = \
	map.mod \
	tide.mod

BUILT_SOURCES = \
	map.c \
	tide.c

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

libtide3D_la_SOURCES = tide.c
libtide3D_la_LIBADD = -L/home/popinet/local/lib -lfes -lnetcdf -lgsl -lgslcblas -lm -lpthread
libtide2D_la_SOURCES = tide.c
libtide2D_la_CFLAGS = $(AM_CFLAGS) -DFTT_2D=1
libtide2D_la_LIBADD = -L/home/popinet/local/lib -lfes -lnetcdf -lgsl -lgslcblas -lm -lpthread
libtide2D3_la_SOURCES = tide.c
libtide2D3_la_CFLAGS = $(AM_CFLAGS) -DFTT_2D3=1
libtide2D3_la_LIBADD = -L/home/popinet/local/lib -lfes -lnetcdf -lgsl -lgslcblas -lm -lpthread

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
