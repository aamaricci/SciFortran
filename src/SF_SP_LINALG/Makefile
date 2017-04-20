include ../../make.inc

OBJS=SF_SP_LINALG.o


.SUFFIXES:.f90 .f


lib: ${OBJS}

.f90.o:	
	$(FC) -c $(FCOPT) $(MOPT)../ $< -o ../$@

.f.o:	
	$(FC) -c $(FCOPT) $< -o ../$@

clean:
	rm -fv *.o *.mod
