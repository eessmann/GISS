## Process this file with automake to produce Makefile.in

INCLUDES = -I$(top_srcdir)/src -I$(includedir) -DG_LOG_DOMAIN=\"Gfs-test\"\
            $(GTS_CFLAGS)
LDADD = $(GFS2D_LIBS)
AM_CFLAGS = -DFTT_2D=1

EXTRA_DIST = match.c read.c test.dom stationary.sh stationary1.sh restart.sh \
	stationary2.sh restartbox.sh

if HAVE_MPI

noinst_PROGRAMS = match read

TESTS = match.sh read.sh restart.sh restartbox.sh \
	stationary.sh stationary1.sh stationary2.sh
else

TESTS = restart.sh restartbox.sh stationary1.sh stationary2.sh

endif

match.sh:
	@echo "#!/bin/sh" > match.sh
	@echo "mpirun -np 2 ./match" >> match.sh
	@chmod +x match.sh

read.sh:
	@echo "#!/bin/sh" > read.sh
	@echo "mpirun -np 3 ./read test.dom" >> read.sh
	@chmod +x read.sh
