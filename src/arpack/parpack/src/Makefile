include ../../../../make.inc
include ../../ARmake.inc

DEST=../../../

SOBJ  = psgetv0.o  \
	psnaitr.o psnapps.o psnaup2.o psnaupd.o psneigh.o psngets.o \
	pssaitr.o pssapps.o pssaup2.o pssaupd.o psseigt.o pssgets.o \
	psneupd.o psseupd.o pslarnv.o pslamch.o psnorm2.o

DOBJ  = pdgetv0.o \
	pdnaitr.o pdnapps.o pdnaup2.o pdnaupd.o pdneigh.o pdngets.o \
	pdsaitr.o pdsapps.o pdsaup2.o pdsaupd.o pdseigt.o pdsgets.o \
	pdneupd.o pdseupd.o pdlarnv.o pdlamch.o pdnorm2.o

COBJ  = pcnaitr.o pcnapps.o pcnaup2.o pcnaupd.o pcneigh.o \
		pcneupd.o pcngets.o pcgetv0.o pscnorm2.o pclarnv.o

ZOBJ  = pznaitr.o pznapps.o pznaup2.o pznaupd.o pzneigh.o \
		pzneupd.o pzngets.o pzgetv0.o pdznorm2.o pzlarnv.o

.SUFFIXES:.f90 .f

.f.o:
	$(MPIFC) -c $(FCOPT) $< -o $(DEST)$@

all: double complex16

single: $(SOBJ)

double: $(DOBJ)

complex: $(COBJ)

complex16: $(ZOBJ)


