## Process this file with automake to produce Makefile.in

if BUILD_WAVEWATCH
WAVEWATCH = libwavewatch.a
endif

noinst_LIBRARIES = $(WAVEWATCH)

EXTRA_DIST = gfsad3 gfsw3initmd.ftn wavewatch.h ww3_grid.inp

libwavewatch.a: gfsw3initmd.o
	rm -f libwavewatch.a
	ar cr libwavewatch.a `grep WWATCH3_DIR $(HOME)/.wwatch3.env | awk '{print $$2 "/obj/w3*.o"}'` gfsw3initmd.o
	ranlib libwavewatch.a

gfsw3initmd.o: gfsw3initmd.ftn
	sh gfsad3 gfsw3initmd.ftn
