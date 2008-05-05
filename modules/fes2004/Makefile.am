if BUILD_TIDE
TIDE = libfes.a
endif

noinst_LIBRARIES = $(TIDE)

libfes_a_CFLAGS = $(AM_CFLAGS) $(GSL_CFLAGS)
libfes_a_SOURCES = \
	fes2004_alloc.c \
	fes2004_error.c \
	fes2004_extraction.c \
	fes2004_init.c \
	fes2004_io.c \
	fes2004_prediction.c \
	fes2004_kernel.c \
	fes.h \
	fes2004_lib.h
