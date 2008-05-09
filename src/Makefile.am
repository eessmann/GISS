## Process this file with automake to produce Makefile.in

AUTOMAKE_OPTIONS = 1.6

AM_CPPFLAGS = -DGFS_MODULES_DIR=\"$(libdir)/gerris\"

INCLUDES = -I$(top_srcdir) -I$(includedir) \
           -DG_LOG_DOMAIN=\"Gfs\" $(GTS_CFLAGS)

bin_SCRIPTS=gfs-config

gfs-config: gfs-config.in
gerris2D.pc: gerris2D.pc.in
gerris2D3.pc: gerris2D3.pc.in
gerris3D.pc: gerris3D.pc.in

pkgconfigdir = $(libdir)/pkgconfig
pkgconfig_DATA = gerris2D.pc gerris3D.pc gerris2D3.pc
pkglib_DATA = m4.awk

lib_LTLIBRARIES = libgfs2D.la libgfs3D.la libgfs2D3.la

BUILT_SOURCES= \
	gfs-config  \
	gerris2D.pc \
	gerris2D3.pc \
	gerris3D.pc

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
	wave.h \
	levelset.h \
	isocube.h \
	cartesian.h \
	surface.h \
	unstructured.h \
	map.h \
	version.h \
	$(MPI_HDS)

pkginclude_HEADERS = \
	$(GFS_HDS) \
	spatial.h \
	function.h

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
	wave.c \
	levelset.c \
	myc.h \
	myc2d.h \
	cartesian.c \
	surface.c \
	unstructured.c \
	map.c \
	$(GFS_HDS)

simulation.c: version.h

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
	m4.awk

bin_PROGRAMS = gerris2D gerris3D gerris2D3

gerris2D_SOURCES = gerris.c
gerris2D_CFLAGS = $(AM_CFLAGS) -DFTT_2D=1
gerris2D_LDADD = $(GFS2D_LIBS)

gerris3D_SOURCES = gerris.c
gerris3D_LDADD = $(GFS3D_LIBS)

gerris2D3_SOURCES = gerris.c
gerris2D3_CFLAGS = $(AM_CFLAGS) -DFTT_2D3=1
gerris2D3_LDADD = $(GFS2D3_LIBS)

if DARCS_CONTROLLED
DARCS_VERSION = darcsversion
else
DARCS_VERSION = 
endif

version.h: $(DARCS_VERSION)

darcsversion:
	sh darcsversion.sh
