include ../../make.inc

OBJS=SF_COLORS.o SF_PAULI.o SF_FONTS.o SF_CONSTANTS.o SF_TIMER.o SF_MISC.o SF_ARRAYS.o


.SUFFIXES:.f90 .f


lib: ${OBJS}

.f90.o:	
	$(FC) -c $(FCOPT) $(MOPT)../ $< -o ../$@

.f.o:	
	$(FC) -c $(FCOPT) $< -o ../$@

clean:
	rm -fv *.o *.mod
