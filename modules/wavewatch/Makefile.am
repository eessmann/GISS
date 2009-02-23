## Process this file with automake to produce Makefile.in

if BUILD_WAVEWATCH
WAVEWATCH = libwavewatch.a
CHECKTYPES = checktypes
BUILT = cfortrantypes.h
endif

noinst_LIBRARIES = $(WAVEWATCH)
noinst_PROGRAMS = $(CHECKTYPES)

EXTRA_DIST = gfsad3 gfsw3initmd.ftn wavewatch.h fchecktypes.ftn

BUILT_SOURCES = $(BUILT)

checktypes_SOURCES = cchecktypes.c
checktypes_LDADD = fchecktypes.o

libwavewatch_a_SOURCES = cfortrantypes.h

libwavewatch.a: gfsw3initmd.o $(HOME)/.wwatch3.env
	rm -f libwavewatch.a
	$(AR) cru libwavewatch.a `grep WWATCH3_DIR $(HOME)/.wwatch3.env | awk '{print $$2 "/obj/w3*.o"}'` gfsw3initmd.o
	$(RANLIB) libwavewatch.a

gfsw3initmd.o: gfsw3initmd.ftn $(HOME)/.wwatch3.env gfsad3
	sh gfsad3 gfsw3initmd.ftn

fchecktypes.o: fchecktypes.ftn $(HOME)/.wwatch3.env gfsad3
	sh gfsad3 fchecktypes.ftn

cfortrantypes.h: checktypes
	./checktypes > cfortrantypes.h
