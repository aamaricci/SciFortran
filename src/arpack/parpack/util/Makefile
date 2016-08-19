include ../../../../make.inc
include ../../ARmake.inc

DEST=../../../

OBJS  = pivout.o

SOBJ  = psvout.o  psmout.o

DOBJ  = pdvout.o  pdmout.o

COBJ  = pcvout.o  pcmout.o 

ZOBJ  = pzvout.o  pzmout.o

.SUFFIXES:.f90 .f

.f.o:	
	$(MPIFC) -c $(FCOPT) $< -o $(DEST)$@

all: double complex16

single: $(SOBJ) $(OBJS)

double: $(DOBJ) $(OBJS)

complex: $(COBJ) $(OBJS)

complex16: $(ZOBJ) $(OBJS)

