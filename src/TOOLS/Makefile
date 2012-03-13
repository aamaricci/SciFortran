include $(SFDIR)/etc/opt.mk
ROUTINES=TOOLS.f90 

all: libtools.a libtools_deb.a

libtools.a:
	@rm -fv *.o *.mod $(SFLIB)/$@	  
	$(FC) -c $(FFLAG) $(ROUTINES) -I$(SFINCLUDE)
	ar cvq $@ `ls *.o | sort | uniq`  
	ranlib $@			  
	@mv -fv *.a $(SFLIB)/		  
	@mv -fv *.mod $(SFINCLUDE)/	  
	@rm -f *.o			  
	@echo 'LIBRARY IS DONE.........'
	@echo ' '

libtools_deb.a:
	@rm -fv *.o *.mod $(SFLIB)/$@	  
	$(FC) -c $(DFLAG) $(ROUTINES) -I$(SFINCLUDE)/debug
	ar cvq $@ `ls *.o | sort | uniq`  
	ranlib $@			  
	@mv -fv *.a $(SFLIB)/		  
	@mv -fv *.mod $(SFINCLUDE)/debug/	  
	@rm -f *.o			  
	@echo 'LIBRARY IS DONE.........'
	@echo ' '
