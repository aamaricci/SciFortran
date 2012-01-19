include $(HOME)/lib/opt.mk

ROUTINES=FFTGF.f90
all: libfftgf.a libfftgf_deb.a

libfftgf.a:
	@rm -fv *.o *.mod $(MYLIB)/$@	  
	$(FC) -c $(FFLAG) $(ROUTINES) -I$(FGSLDIR)/include/ -L$(FGSLDIR)/lib/ -L$(LIBOPT)/gsl -lfgsl_ifort -lgsl -lgslcblas
	ar cvq $@ `ls *.o | sort | uniq`  
	ranlib $@			  
	@mv -fv *.a $(MYLIB)/		  
	@mv -fv *.mod $(MYINCLUDE)/	  
	@rm -f *.o			  
	@echo 'LIBRARY IS DONE.........'
	@echo ' '

libfftgf_deb.a:
	@rm -fv *.o *.mod $(MYLIB)/$@	  
	$(FC) -c $(DFLAG) $(ROUTINES) -I$(FGSLDIR)/include/ -L$(FGSLDIR)/lib/ -L$(LIBOPT)/gsl -lfgsl_ifort -lgsl -lgslcblas
	ar cvq $@ `ls *.o | sort | uniq`  
	ranlib $@			  
	@mv -fv *.a $(MYLIB)/		  
	@mv -fv *.mod $(MYINCLUDE)/debug/	  
	@rm -f *.o			  
	@echo 'LIBRARY IS DONE.........'
	@echo ' '
