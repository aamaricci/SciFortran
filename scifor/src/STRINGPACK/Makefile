include $(SFDIR)/etc/opt.mk

ROUTINES  = CHRPACK.f90 STRING_FUNX_1.f90 STRING_FUNX_2.f90 STRINGPACK.f90

all: libstringpack.a libstringpack_deb.a

libstringpack.a:
	@rm -fv *.o *.mod $(SFLIB)/$@	  
	$(FC) -c $(FFLAG) $(ROUTINES) 
	ar cvq $@ `ls *.o | sort | uniq`  
	ranlib $@
	@mv -fv *.a $(SFLIB)/
	@mv -fv *.mod $(SFINCLUDE)/
	@rm -f *.o			  
	@echo 'LIBRARY IS DONE.........'
	@echo ' '

libstringpack_deb.a:
	@rm -fv *.o *.mod $(SFLIB)/$@	  
	$(FC) -c $(DFLAG) $(ROUTINES)
	ar cvq $@ `ls *.o | sort | uniq`  
	ranlib $@			  
	@mv -fv *.a $(SFLIB)/		  
	@mv -fv *.mod $(SFINCLUDE)/debug/
	@rm -f *.o			  
	@echo 'LIBRARY IS DONE.........'
	@echo ' '




test: chrpak_test.f90
	$(FC) $(FFLAG) -o test chrpak_test.f90 -L./lib -lstringpack

