include ../../make.inc

OBJS=chkder.o dogleg.o enorm.o enorm2.o fdjac1.o fdjac2.o hybrd.o hybrd1.o hybrj.o hybrj1.o lmder.o lmder1.o lmdif.o lmdif1.o lmpar.o lmstr.o lmstr1.o qform.o qrfac.o qrsolv.o r1mpyq.o r1updt.o r8vec_print.o rwupdt.o timestamp.o

.SUFFIXES:.f90 .f


lib: ${OBJS}

.f90.o:	
	$(FC) -c $(FCOPT) $(MOPT)../ $< -o ../$@

.f.o:	
	$(FC) -c $(FCOPT) $< -o ../$@

clean:
	rm -fv *.o *.mod
