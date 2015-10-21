include ../make.inc

.SUFFIXES:.f90 .f

all: allmod lib

allmod: chkder.o dogleg.o enorm.o enorm2.o fdjac1.o fdjac2.o hybrd.o hybrd1.o hybrj.o hybrj1.o lmder.o lmder1.o lmdif.o lmdif1.o lmpar.o lmstr.o lmstr1.o qform.o qrfac.o qrsolv.o r1mpyq.o r1updt.o r8vec_print.o rwupdt.o timestamp.o

lib:
	ar cvq $(LIBMINPACK) `ls *.o | sort | uniq`  
	${RANLIB} $(LIBMINPACK)
	@echo "$(LIBMINPACK) is done..."
	@echo ' '
	rsync -avP *.o $(OBJ_INSTALL)/

.f90.o:	
	$(FC) -c $(FFLAG) $<

.f.o:	
	$(FC) -c $(FFLAG) $< 

clean:
	rm -fv  *.o
