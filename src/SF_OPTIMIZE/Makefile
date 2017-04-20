include ../../make.inc

OBJS=optimize_broyden_routines.o OPTIMIZE_ROOT_FINDING.o minimize.o cgfam.o cgsearch.o cgblas.o optimize_cgfit_routines.o OPTIMIZE_ROOT_FINDING.o OPTIMIZE_MINIMIZE.o SF_OPTIMIZE.o


.SUFFIXES:.f90 .f


lib: ${OBJS}

.f90.o:	
	$(FC) -c $(FCOPT) $(MOPT)../ $< -o ../$@

.f.o:	
	$(FC) -c $(FCOPT) $< -o ../$@

clean:
	rm -fv *.o *.mod
