include $(HOME)/lib/flag.mk

ROUTINES= fftpack.f
OBJ  =   fftpack.o

all: comp libfftpack.a

comp:  ${ROUTINES}
	@rm -fv *.o *.mod 
	echo 'files erased...'
	@rm -f $(MYLIB)/libfftpack.a
	echo 'old lib erased...'
	@echo "compile"
	$(FC) -c $(FFLAG) ${ROUTINES} 

libfftpack.a:
	ar cvq $@ `ls *.o | sort | uniq`
	${RANLIB} $@
	@mv -fv *.a $(MYLIB)/
	@echo 'LIBRARY IS DONE.........'
	echo 'files erased....'
	@rm -f *.o
	@echo ' '
clean:
	rm -fv $(OBJ
	rm -fv $(MYLIB)/libfftpack.a
