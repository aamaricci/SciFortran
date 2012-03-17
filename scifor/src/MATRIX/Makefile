include $(SFDIR)/etc/opt.mk

ROUTINES=MATRIX.f90 

all: libmatrix.a libmatrix_deb.a

libmatrix.a:
	@rm -fv *.o *.mod $(SFLIB)/$@	  
	$(FC) -c $(FFLAG) $(ROUTINES)
	ar cvq $@ `ls *.o | sort | uniq`  
	ranlib $@			  
	@mv -fv *.a $(SFLIB)/		  
	@mv -fv *.mod $(SFINCLUDE)/	  
	@rm -f *.o			  
	@echo 'LIBRARY IS DONE.........'
	@echo ' '

libmatrix_deb.a:
	@rm -fv *.o *.mod $(SFLIB)/$@	  
	$(FC) -c $(DFLAG) $(ROUTINES)
	ar cvq $@ `ls *.o | sort | uniq`  
	ranlib $@			  
	@mv -fv *.a $(SFLIB)/		  
	@mv -fv *.mod $(SFINCLUDE)/debug/	  
	@rm -f *.o			  
	@echo 'LIBRARY IS DONE.........'
	@echo ' '


