## Process this file with automake to produce Makefile.in

INCLUDES = -I$(top_srcdir)/src -I$(includedir) -DG_LOG_DOMAIN=\"Gfs-test\"\
            $(GTS_CFLAGS)
LDADD = $(GFS2D_LIBS)
CFLAGS += -DFTT_2D=1

EXTRA_DIST = \
	reference2D \
	parameters \
	scripts

REFERENCE = reference2D

test.sh:
	@echo "#! /bin/sh" > test.sh
	@echo "rm -f report.ps" >> test.sh
	@echo "./scripts/check.sh \"$(REFERENCE)/*\"" >> test.sh
	@chmod +x test.sh

report.ps: scripts/report.sh timestamp
	./scripts/report.sh "$(REFERENCE)/*"

clean-generic:
	./scripts/cleanup.sh

dist-hook:
	./scripts/cleanup.sh

if HAVE_MPI
  TESTS = test.sh ./scripts/mpi.sh
else
  TESTS = test.sh
endif
