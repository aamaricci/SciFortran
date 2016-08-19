include ../../make.inc


OBJS=quadpack_qagi.o quadpack_qags.o quadpack_qawf.o quadpack_qaws.o quadpack_qag.o  quadpack_qagp.o quadpack_qawc.o quadpack_qawo.o quadpack_qng.o quadpack_aux.o MOD_QUADPACK.o

.SUFFIXES:.f90 .f


lib: ${OBJS}

.f90.o:	
	$(FC) -c $(FCOPT) $(MOPT)../ $< -o ../$@

.f.o:	
	$(FC) -c $(FCOPT) $< -o ../$@

clean:
	rm -fv *.o *.mod
