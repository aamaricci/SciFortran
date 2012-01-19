include $(HOME)/lib/opt.mk

ROUTINES  = CHRPACK.f90

all: libchrpack.a libchrpack_deb.a

libchrpack.a:
	@rm -fv *.o *.mod $(MYLIB)/$@	  
	$(FC) -c $(FFLAG) $(ROUTINES) 
	ar cvq $@ `ls *.o | sort | uniq`  
	ranlib $@			  
	@mv -fv *.a $(MYLIB)/
	@rm -f *.o			  
	@echo 'LIBRARY IS DONE.........'
	@echo ' '

libchrpack_deb.a:
	@rm -fv *.o *.mod $(MYLIB)/$@	  
	$(FC) -c $(DFLAG) $(ROUTINES)
	ar cvq $@ `ls *.o | sort | uniq`  
	ranlib $@			  
	@mv -fv *.a $(MYLIB)/		  
	@rm -f *.o			  
	@echo 'LIBRARY IS DONE.........'
	@echo ' '




test: chrpak_test.f90
	$(FC) $(FFLAG) -o test chrpak_test.f90 -L./lib -lchrpack

