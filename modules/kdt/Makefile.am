## Process this file with automake to produce Makefile.in

noinst_LTLIBRARIES = libkdt.la

bin_PROGRAMS = \
	xyz2kdt \
	kdtquery

libkdt_la_CFLAGS = $(AM_CFLAGS) -D_FILE_OFFSET_BITS=64
libkdt_la_SOURCES = \
	kdt.c \
	kdt.h

xyz2kdt_SOURCES = xyz2kdt.c kdt.h
xyz2kdt_LDADD = -lkdt
xyz2kdt_CFLAGS = $(AM_CFLAGS)
xyz2kdt_DEPENDENCIES = libkdt.la

kdtquery_SOURCES = kdtquery.c kdt.h
kdtquery_LDADD = -lkdt
kdtquery_CFLAGS = $(AM_CFLAGS)
kdtquery_DEPENDENCIES = libkdt.la
