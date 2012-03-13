include $(SFDIR)/etc/opt.mk

ROUTINES=SQUARE_LATTICE.f90 

all: liblattice.a liblattice_deb.a

liblattice.a:
	@rm -fv *.o *.mod $(SFLIB)/$@	  
	$(FC) -c $(FFLAG) $(ROUTINES) -I$(SFINCLUDE)
	ar cvq $@ `ls *.o | sort | uniq`  
	ranlib $@			  
	@mv -fv *.a $(SFLIB)/		  
	@mv -fv *.mod $(SFINCLUDE)/	  
	@rm -f *.o			  
	@echo 'LIBRARY IS DONE.........'
	@echo ' '

liblattice_deb.a:
	@rm -fv *.o *.mod $(SFLIB)/$@	  
	$(FC) -c $(DFLAG) $(ROUTINES) -I$(SFINCLUDE)/debug
	ar cvq $@ `ls *.o | sort | uniq`  
	ranlib $@			  
	@mv -fv *.a $(SFLIB)/		  
	@mv -fv *.mod $(SFINCLUDE)/debug/	  
	@rm -f *.o			  
	@echo 'LIBRARY IS DONE.........'
	@echo ' '


