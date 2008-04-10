# Make RStarTree C-Version R2.0

# Cflags = -g -ansi -pedantic -Wall
# Cflags = -g -ansi -pedantic
Cflags = -O -ansi -pedantic -Wall
# Cflags = -O -ansi -pedantic

CC = gcc


compile = $(CC) $(Cflags) -c
link = $(CC) $(Cflags) -o


RSTos = RStarTree.o RSTInterUtil.o RSTInOut.o RSTInstDel.o RSTUtil.o \
        RSTQuery.o  RSTJoin.o


libraries = libcSmRST.a
programs = trst


  RANLIB = ranlib $(libraries)          # Berkely
# RANLIB = @ echo 'skipping ranlib'     # SVR4


all: $(libraries) $(programs)


#
# Libraries
#
libcSmRST.a: $(RSTos)
	ar rcv libcSmRST.a $(RSTos)
	$(RANLIB)
	mv README ..
	mv libcSmRST.a ..
	cp RStarTree.h ..
	cp RSTBase.h ..

#
# Main Programs
#
trst: trst.o \
      libcSmRST.a \
      Set.o
	$(link) trst trst.o Set.o -L.. -lcSmRST

trst.o: trst.c \
        RStarTree.h \
        Set.h
	$(compile) trst.c

#
# compile the modules.
#
RStarTree.o: RStarTree.c \
             RStarTree.h \
             RSTInterUtil.h \
             RSTInOut.h \
             RSTInstDel.h \
             RSTQuery.h \
             RSTJoin.h
	$(compile) RStarTree.c

RSTInterUtil.o: RSTInterUtil.c \
                RSTInterUtil.h \
                RStarTree.h
	$(compile) RSTInterUtil.c

RSTInOut.o: RSTInOut.c \
            RSTInOut.h \
            RStarTree.h
	$(compile) RSTInOut.c

RSTInstDel.o: RSTInstDel.c \
              RSTInstDel.h \
              RStarTree.h \
              RSTUtil.h \
              RSTInOut.h
	$(compile) RSTInstDel.c

RSTUtil.o: RSTUtil.c \
           RSTUtil.h \
           RStarTree.h
	$(compile) RSTUtil.c

RSTQuery.o: RSTQuery.c \
            RSTQuery.h \
            RStarTree.h \
            RSTUtil.h \
            RSTInOut.h
	$(compile) RSTQuery.c

RSTJoin.o: RSTJoin.c \
           RSTJoin.h \
           RStarTree.h \
           RSTUtil.h \
           RSTInOut.h
	$(compile) RSTJoin.c

Set.o: Set.c \
       Set.h
	$(compile) Set.c

#
# interface dependencies
#
RStarTree.h: RSTBase.h
	touch RStarTree.h

RSTInterUtil.h: RStarTree.h
	touch RSTInterUtil.h

RSTInOut.h: RStarTree.h
	touch RSTInOut.h

RSTInstDel.h: RStarTree.h
	touch RSTInstDel.h

RSTUtil.h: RStarTree.h
	touch RSTUtil.h

RSTQuery.h: RStarTree.h
	touch RSTQuery.h

RSTJoin.h: RStarTree.h
	touch RSTJoin.h

Set.h: SetDef.h
	touch Set.h
