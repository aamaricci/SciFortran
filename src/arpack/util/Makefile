include ../../../make.inc
include ../ARmake.inc

DEST=../../


OBJS  = icnteq.o icopy.o iset.o iswap.o ivout.o second.o

SOBJ  = svout.o  smout.o

DOBJ  = dvout.o  dmout.o

COBJ  = cvout.o  cmout.o 

ZOBJ  = zvout.o  zmout.o

.SUFFIXES:.f90 .f

.f.o:	
	$(FC) -c $(FCOPT) $< -o $(DEST)$@

all: double complex16

single: $(SOBJ) $(OBJS)

double: $(DOBJ) $(OBJS)

complex: $(COBJ) $(OBJS)

complex16: $(ZOBJ) $(OBJS)


