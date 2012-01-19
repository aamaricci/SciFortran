include $(HOME)/lib/opt.mk

ROUTINES=MATRIX.f90 

all: libmatrix.a libmatrix_deb.a

libmatrix.a:
	@rm -fv *.o *.mod $(MYLIB)/$@	  
	$(FC) -c $(FFLAG) $(ROUTINES)
	ar cvq $@ `ls *.o | sort | uniq`  
	ranlib $@			  
	@mv -fv *.a $(MYLIB)/		  
	@mv -fv *.mod $(MYINCLUDE)/	  
	@rm -f *.o			  
	@echo 'LIBRARY IS DONE.........'
	@echo ' '

libmatrix_deb.a:
	@rm -fv *.o *.mod $(MYLIB)/$@	  
	$(FC) -c $(DFLAG) $(ROUTINES)
	ar cvq $@ `ls *.o | sort | uniq`  
	ranlib $@			  
	@mv -fv *.a $(MYLIB)/		  
	@mv -fv *.mod $(MYINCLUDE)/debug/	  
	@rm -f *.o			  
	@echo 'LIBRARY IS DONE.........'
	@echo ' '


