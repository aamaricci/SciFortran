include $(SFDIR)/etc/opt.mk

ROUTINES = BRENT.f90 BROYDEN.f90 ROOTFIND.f90


all: libroot.a libroot_deb.a

libroot.a:
	@rm -fv *.o *.mod $(SFLIB)/$@	  
	$(FC) -c $(FFLAG) $(ROUTINES) 
	ar cvq $@ `ls *.o | sort | uniq`  
	ranlib $@			  
	@mv -fv *.a $(SFLIB)/		  
	@mv -fv *.mod $(SFINCLUDE)/	  
	@rm -f *.o			  
	@echo 'LIBRARY IS DONE.........'
	@echo ' '

libroot_deb.a:
	@rm -fv *.o *.mod $(SFLIB)/$@	  
	$(FC) -c $(DFLAG) $(ROUTINES)
	ar cvq $@ `ls *.o | sort | uniq`  
	ranlib $@			  
	@mv -fv *.a $(SFLIB)/		  
	@mv -fv *.mod $(SFINCLUDE)/debug/	  
	@rm -f *.o			  
	@echo 'LIBRARY IS DONE.........'
	@echo ' '



