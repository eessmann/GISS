## Process this file with automake to produce Makefile.in

INCLUDES = -I$(top_srcdir) -I$(includedir) \
           -DG_LOG_DOMAIN=\"Gfs\" $(GTS_CFLAGS)

bin_SCRIPTS=gfs-config

gfs-config: gfs-config.in
gerris2D.pc: gerris2D.pc.in
gerris3D.pc: gerris3D.pc.in
dx2D.mdf: dx2D.mdf.in
dx3D.mdf: dx3D.mdf.in

pkgconfigdir = $(libdir)/pkgconfig
pkgconfig_DATA = gerris2D.pc gerris3D.pc

lib_LTLIBRARIES = libgfs2D.la libgfs3D.la libgfs2D3.la

if HAVE_DX
DX_MODS = libdx2D.la libdx3D.la 
DX_MDF = dx2D.mdf dx3D.mdf
DX_2DSRC = dx.c dx2D.c
DX_3DSRC = dx.c dx3D.c
DX_BUILT_SRC = dx2D.c dx3D.c
MDF2C = $(DX_PATH)/bin/mdf2c
else
DX_MODS = 
DX_MDF =
DX_2DSRC =
DX_3DSRC =
DX_BUILT_SRC =
MDF2C = cat
endif

BUILT_SOURCES= \
	gfs-config  \
	gerris2D.pc \
	gerris3D.pc \
	dx2D.mdf    \
	dx3D.mdf    \
	$(DX_BUILT_SRC)

pkglib_LTLIBRARIES = $(DX_MODS)
dist_pkglib_DATA = $(DX_MDF)

if HAVE_MPI
MPI_SRC = mpi_boundary.c
MPI_HDS = mpi_boundary.h
else
MPI_SRC = 
MPI_HDS = 
endif

GFS_HDS = \
	ftt.h \
	fluid.h \
	variable.h \
	output.h \
	solid.h \
	poisson.h \
	advection.h \
	boundary.h \
	timestep.h \
	domain.h \
	init.h \
	refine.h \
	event.h \
	simulation.h \
	graphic.h \
	adaptive.h \
	source.h \
	tension.h \
	vof.h \
	utils.h \
	ocean.h \
	$(MPI_HDS)

pkginclude_HEADERS = \
	$(GFS_HDS)

include_HEADERS = \
	gfs.h \
	gfsconfig.h

SRC = \
	ftt.c \
	fluid.c \
	variable.c \
	output.c \
	solid.c \
	poisson.c \
	advection.c \
	boundary.c \
	timestep.c \
	domain.c \
	$(MPI_SRC) \
	init.c \
	refine.c \
	event.c \
	simulation.c \
	graphic.c \
	adaptive.c \
	source.c \
	tension.c \
	vof.c \
	utils.c \
	ocean.c \
	$(GFS_HDS)

libgfs3D_la_LDFLAGS = $(NO_UNDEFINED)\
        -version-info $(LT_CURRENT):$(LT_REVISION):$(LT_AGE)\
	-release $(LT_RELEASE) -export-dynamic
libgfs3D_la_SOURCES = $(SRC)

libgfs2D_la_LDFLAGS = $(NO_UNDEFINED)\
        -version-info $(LT_CURRENT):$(LT_REVISION):$(LT_AGE)\
	-release $(LT_RELEASE) -export-dynamic
libgfs2D_la_SOURCES = $(SRC)
libgfs2D_la_CFLAGS = $(AM_CFLAGS) -DFTT_2D=1

libgfs2D3_la_LDFLAGS = $(NO_UNDEFINED)\
        -version-info $(LT_CURRENT):$(LT_REVISION):$(LT_AGE)\
	-release $(LT_RELEASE) -export-dynamic
libgfs2D3_la_SOURCES = $(SRC)
libgfs2D3_la_CFLAGS = $(AM_CFLAGS) -DFTT_2D3=1

CLEANFILES = $(BUILT_SOURCES)

EXTRA_DIST = \
	mpi_boundary.c \
	mpi_boundary.h \
	ftt_internal.c \
	dx.c

bin_PROGRAMS = gerris2D gerris3D gerris2D3 gerris2D_dir gerris3D_dir

gerris2D_SOURCES = gerris.c
gerris2D_CFLAGS = $(AM_CFLAGS) -DFTT_2D=1
gerris2D_LDADD = $(GFS2D_LIBS)

gerris3D_SOURCES = gerris.c
gerris3D_LDADD = $(GFS3D_LIBS)

gerris2D3_SOURCES = gerris.c
gerris2D3_CFLAGS = $(AM_CFLAGS) -DFTT_2D3=1
gerris2D3_LDADD = $(GFS2D3_LIBS)

gerris2D_dir_SOURCES = gerris_dir.c
gerris2D_dir_CFLAGS = $(AM_CFLAGS) -DFTT_2D=1
gerris2D_dir_LDADD = $(GFS2D_LIBS)

gerris3D_dir_SOURCES = gerris_dir.c
gerris3D_dir_LDADD = $(GFS3D_LIBS)

.mdf.c:
	$(MDF2C) $^ > $@

libdx3D_la_LDFLAGS = -version-info $(LT_CURRENT):$(LT_REVISION):$(LT_AGE)\
	-release $(LT_RELEASE)
libdx3D_la_SOURCES = $(DX_3DSRC)
libdx3D_la_CFLAGS = $(AM_CFLAGS) $(DX_CFLAGS)
libdx3D_la_LIBADD = $(GFS3D_LIBS)

libdx2D_la_LDFLAGS = -version-info $(LT_CURRENT):$(LT_REVISION):$(LT_AGE)\
	-release $(LT_RELEASE)
libdx2D_la_SOURCES = $(DX_2DSRC)
libdx2D_la_CFLAGS = $(AM_CFLAGS) $(DX_CFLAGS) -DFTT_2D=1
libdx2D_la_LIBADD = $(GFS2D_LIBS)
