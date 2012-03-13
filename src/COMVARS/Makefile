include $(SFDIR)/etc/opt.mk

ROUTINES=COMVARS.f90
REV = $(shell git rev-parse HEAD)

all: libcomvars.a libcomvars_deb.a

libcomvars.a:
	@rm -fv *.o *.mod $(SFLIB)/$@
	@echo 'character(len=41),parameter,public :: sf_version = "$(REV)"' > scifor_version.inc
	$(FC) -c $(FFLAG) $(ROUTINES) 
	ar cvq $@ `ls *.o | sort | uniq`  
	ranlib $@			  
	@mv -fv *.a $(SFLIB)/		  
	@mv -fv *.mod $(SFINCLUDE)/	  
	@rm -f *.o			  
	@echo 'LIBRARY IS DONE.........'
	@echo ' '

libcomvars_deb.a:
	@rm -fv *.o *.mod $(SFLIB)/$@	  
	@echo 'character(len=41),parameter :: sf_version = "$(REV)"' > scifor_version.inc
	$(FC) -c $(DFLAG) $(ROUTINES)
	ar cvq $@ `ls *.o | sort | uniq`  
	ranlib $@			  
	@mv -fv *.a $(SFLIB)/		  
	@mv -fv *.mod $(SFINCLUDE)/debug/	  
	@rm -f *.o			  
	@echo 'LIBRARY IS DONE.........'
	@echo ' '




