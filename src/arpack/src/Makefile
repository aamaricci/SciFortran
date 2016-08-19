include ../../../make.inc
include ../ARmake.inc

DEST=../../

SOBJ  = sgetv0.o slaqrb.o sstqrb.o ssortc.o ssortr.o sstatn.o sstats.o \
	snaitr.o snapps.o snaup2.o snaupd.o snconv.o sneigh.o sngets.o \
	ssaitr.o ssapps.o ssaup2.o ssaupd.o ssconv.o sseigt.o ssgets.o \
	sneupd.o sseupd.o ssesrt.o

DOBJ  = dgetv0.o dlaqrb.o dstqrb.o dsortc.o dsortr.o dstatn.o dstats.o \
	dnaitr.o dnapps.o dnaup2.o dnaupd.o dnconv.o dneigh.o dngets.o \
	dsaitr.o dsapps.o dsaup2.o dsaupd.o dsconv.o dseigt.o dsgets.o \
	dneupd.o dseupd.o dsesrt.o

COBJ  = cnaitr.o cnapps.o cnaup2.o cnaupd.o cneigh.o cneupd.o cngets.o \
        cgetv0.o csortc.o cstatn.o

ZOBJ  = znaitr.o znapps.o znaup2.o znaupd.o zneigh.o zneupd.o zngets.o \
        zgetv0.o zsortc.o zstatn.o



.SUFFIXES:.f90 .f

all: double complex16

single: $(SOBJ)

double: $(DOBJ)

complex: $(COBJ)

complex16: $(ZOBJ)

.f90.o:	
	$(FC) -c $(FCOPT) $(MOPT)$(DEST) $< -o $(DEST)$@

.f.o:	
	$(FC) -c $(FCOPT) $< -o $(DEST)$@

clean:
	rm -fv *.o *.mod

