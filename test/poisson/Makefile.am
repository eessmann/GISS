## Process this file with automake to produce Makefile.in

SUBDIRS = convergence order

INCLUDES = -I$(top_srcdir)/src -I$(includedir) -DG_LOG_DOMAIN=\"Gfs-test\"\
            $(GTS_CFLAGS)
LDADD = $(GFS2D_LIBS)
AM_CFLAGS = -DFTT_2D=1

noinst_PROGRAMS = \
	simple

EXTRA_DIST = \
	command.awk \
	options.awk \
	models
