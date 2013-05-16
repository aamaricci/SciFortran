include options.inc
REV = $(shell git rev-parse HEAD)

all: version allmod libscifor.a libscifor_deb.a

allmod: COMVARS PARSECMD TIMER IOTOOLS TOOLS OPTIMIZE FFTGF FUNCTIONS GREENFUNX INTEGRATE  SQUARE_LATTICE LIST MATRIX SPLINE RANDOM STATISTICS PADE SCIFOR

libscifor.a: 
	@rsync -avP mods/*.mod $(SFINCLUDE)/
	ar cvq $@ `ls objs/*.o | sort | uniq`  
	ranlib $@
	@mv -fv *.a $(SFLIB)/
	@echo 'LIBRARY IS DONE.........'
	@echo ' '

libscifor_deb.a: 
	@rsync -avP debug_mods/*.mod $(SFINCLUDE)/debug/
	ar cvq $@ `ls debug_objs/*.o | sort | uniq`  
	ranlib $@
	@mv -fv *.a $(SFLIB)/
	@echo 'LIBRARY IS DONE.........'
	@echo ' '

version:
	@echo 'character(len=41),parameter,public :: sf_version = "$(REV)"' > scifor_version.inc


clean: 
	@rm -fv *.o *.mod *.a mods/* objs/* debug_*/* $(SFINCLUDE)/*.mod $(SFINCLUDE)/debug/*.mod $(SFLIB)/*scifor*.a


#STD
COMVARS: 
	$(FC) -c $(FFLAG) $@.f90 -I./mods/ -o objs/$@.o $(MOPT)./mods/
	$(FC) -c $(DFLAG) $@.f90 -I./debug_mods/ -o debug_objs/$@.o $(MOPT)./debug_mods/
mod_COMVARS: COMVARS
	@rsync -avP mods/*.mod $(SFINCLUDE)/
	@rsync -avP debug_mods/*.mod $(SFINCLUDE)/debug/
	ar cvq libscifor.a `ls objs/*.o | sort | uniq`  
	ar cvq libscifor_deb.a `ls debug_objs/*.o | sort | uniq`  
	ranlib libscifor.a
	ranlib libscifor_deb.a
	@mv -fv *.a $(SFLIB)/
	@echo 'LIBRARY IS UPDATED.........'
	@echo ' '

#STD
PARSECMD: 
	$(FC) -c $(FFLAG) $@.f90 -I./mods/ -o objs/$@.o $(MOPT)./mods/
	$(FC) -c $(DFLAG) $@.f90 -I./debug_mods/ -o debug_objs/$@.o $(MOPT)./debug_mods/
mod_PARSECMD: PARSECMD
	@rsync -avP mods/*.mod $(SFINCLUDE)/
	@rsync -avP debug_mods/*.mod $(SFINCLUDE)/debug/
	ar cvq libscifor.a `ls objs/*.o | sort | uniq`  
	ar cvq libscifor_deb.a `ls debug_objs/*.o | sort | uniq`  
	ranlib libscifor.a
	ranlib libscifor_deb.a
	@mv -fv *.a $(SFLIB)/
	@echo 'LIBRARY IS UPDATED.........'
	@echo ' '


BROYDEN: 
	$(FC) -c $(FFLAG) $@.f90 -I./mods/ -o objs/$@.o $(MOPT)./mods/
	$(FC) -c $(DFLAG) $@.f90 -I./debug_mods/ -o debug_objs/$@.o $(MOPT)./debug_mods/


BRENT: 
	$(FC) -c $(FFLAG) $@.f90 -I./mods/ -o objs/$@.o $(MOPT)./mods/
	$(FC) -c $(DFLAG) $@.f90 -I./debug_mods/ -o debug_objs/$@.o $(MOPT)./debug_mods/

WRAP_MINPACK: 
	$(FC) -c $(FFLAG) $@.f90 -I./mods/ -o objs/$@.o $(MOPT)./mods/
	$(FC) -c $(DFLAG) $@.f90 -I./debug_mods/ -o debug_objs/$@.o $(MOPT)./debug_mods/


OPTIMIZE: BROYDEN BRENT WRAP_MINPACK
	$(FC) -c $(FFLAG) $@.f90 -I./mods/ -o objs/$@.o $(MOPT)./mods/
	$(FC) -c $(DFLAG) $@.f90 -I./debug_mods/ -o debug_objs/$@.o $(MOPT)./debug_mods/
mod_OPTIMIZE: OPTIMIZE 
	@rsync -avP mods/*.mod $(SFINCLUDE)/
	@rsync -avP debug_mods/*.mod $(SFINCLUDE)/debug/
	ar cvq libscifor.a `ls objs/*.o | sort | uniq`  
	ar cvq libscifor_deb.a `ls debug_objs/*.o | sort | uniq`  
	ranlib libscifor.a
	ranlib libscifor_deb.a
	@mv -fv *.a $(SFLIB)/
	@echo 'LIBRARY IS UPDATED.........'
	@echo ' '


TIMER: 
	$(FC) -c $(FFLAG) $@.f90 -I./mods/ -o objs/$@.o $(MOPT)./mods/
	$(FC) -c $(DFLAG) $@.f90 -I./debug_mods/ -o debug_objs/$@.o $(MOPT)./debug_mods/
mod_TIMER: TIMER 
	@rsync -avP mods/*.mod $(SFINCLUDE)/
	@rsync -avP debug_mods/*.mod $(SFINCLUDE)/debug/
	ar cvq libscifor.a `ls objs/*.o | sort | uniq`  
	ar cvq libscifor_deb.a `ls debug_objs/*.o | sort | uniq`  
	ranlib libscifor.a
	ranlib libscifor_deb.a
	@mv -fv *.a $(SFLIB)/
	@echo 'LIBRARY IS UPDATED.........'
	@echo ' '

FFTGF:  COMVARS SPLINE
	$(FC) -c $(FFLAG) $@.f90 -I./mods/ -o objs/$@.o $(MOPT)./mods/
	$(FC) -c $(DFLAG) $@.f90 -I./debug_mods/ -o debug_objs/$@.o $(MOPT)./debug_mods/
mod_FFTGF: FFTGF 
	@rsync -avP mods/*.mod $(SFINCLUDE)/
	@rsync -avP debug_mods/*.mod $(SFINCLUDE)/debug/
	ar cvq libscifor.a `ls objs/*.o | sort | uniq`  
	ar cvq libscifor_deb.a `ls debug_objs/*.o | sort | uniq`  
	ranlib libscifor.a
	ranlib libscifor_deb.a
	@mv -fv *.a $(SFLIB)/
	@echo 'LIBRARY IS UPDATED.........'
	@echo ' '

FUNCTIONS: 
	$(FC) -c $(FFLAG) $@.f90 -I./mods/ -o objs/$@.o $(MOPT)./mods/
	$(FC) -c $(DFLAG) $@.f90 -I./debug_mods/ -o debug_objs/$@.o $(MOPT)./debug_mods/
mod_FUNCTIONS: FUNCTIONS 
	@rsync -avP mods/*.mod $(SFINCLUDE)/
	@rsync -avP debug_mods/*.mod $(SFINCLUDE)/debug/
	ar cvq libscifor.a `ls objs/*.o | sort | uniq`  
	ar cvq libscifor_deb.a `ls debug_objs/*.o | sort | uniq`  
	ranlib libscifor.a
	ranlib libscifor_deb.a
	@mv -fv *.a $(SFLIB)/
	@echo 'LIBRARY IS UPDATED.........'
	@echo ' '


GREENFUNX: 
	$(FC) -c $(FFLAG) $@.f90 -I./mods/ -o objs/$@.o $(MOPT)./mods/
	$(FC) -c $(DFLAG) $@.f90 -I./debug_mods/ -o debug_objs/$@.o $(MOPT)./debug_mods/
mod_GREENFUNX: GREENFUNX 
	@rsync -avP mods/*.mod $(SFINCLUDE)/
	@rsync -avP debug_mods/*.mod $(SFINCLUDE)/debug/
	ar cvq libscifor.a `ls objs/*.o | sort | uniq`  
	ar cvq libscifor_deb.a `ls debug_objs/*.o | sort | uniq`  
	ranlib libscifor.a
	ranlib libscifor_deb.a
	@mv -fv *.a $(SFLIB)/
	@echo 'LIBRARY IS UPDATED.........'
	@echo ' '

INTEGRATE: 
	$(FC) -c $(FFLAG) integrate_d_quadpack.f90 -I./mods/ -o objs/integrate_d_quadpack.o $(MOPT)./mods/
	$(FC) -c $(FFLAG)  $@.f90 -I./mods/ -o objs/$@.o $(MOPT)./mods/
	$(FC) -c $(DFLAG) integrate_d_quadpack.f90 -I./debug_mods/ -o debug_objs/integrate_d_quadpack.o $(MOPT)./debug_mods/
	$(FC) -c $(DFLAG)  $@.f90 -I./debug_mods/ -o debug_objs/$@.o $(MOPT)./debug_mods/
mod_INTEGRATE: INTEGRATE 
	@rsync -avP mods/*.mod $(SFINCLUDE)/
	@rsync -avP debug_mods/*.mod $(SFINCLUDE)/debug/
	ar cvq libscifor.a `ls objs/*.o | sort | uniq`  
	ar cvq libscifor_deb.a `ls debug_objs/*.o | sort | uniq`  
	ranlib libscifor.a
	ranlib libscifor_deb.a
	@mv -fv *.a $(SFLIB)/
	@echo 'LIBRARY IS UPDATED.........'
	@echo ' '

SQUARE_LATTICE: 
	$(FC) -c $(FFLAG) $@.f90 -I./mods/ -o objs/$@.o $(MOPT)./mods/
	$(FC) -c $(DFLAG) $@.f90 -I./debug_mods/ -o debug_objs/$@.o $(MOPT)./debug_mods/
mod_SQUARE_LATTICE: SQUARE_LATTICE 
	@rsync -avP mods/*.mod $(SFINCLUDE)/
	@rsync -avP debug_mods/*.mod $(SFINCLUDE)/debug/
	ar cvq libscifor.a `ls objs/*.o | sort | uniq`  
	ar cvq libscifor_deb.a `ls debug_objs/*.o | sort | uniq`  
	ranlib libscifor.a
	ranlib libscifor_deb.a
	@mv -fv *.a $(SFLIB)/
	@echo 'LIBRARY IS UPDATED.........'
	@echo ' '

MATRIX: 
	$(FC) -c $(FFLAG) $@.f90 -I./mods/ -o objs/$@.o $(MOPT)./mods/
	$(FC) -c $(DFLAG) $@.f90 -I./debug_mods/ -o debug_objs/$@.o $(MOPT)./debug_mods/
mod_MATRIX: MATRIX 
	@rsync -avP mods/*.mod $(SFINCLUDE)/
	@rsync -avP debug_mods/*.mod $(SFINCLUDE)/debug/
	ar cvq libscifor.a `ls objs/*.o | sort | uniq`  
	ar cvq libscifor_deb.a `ls debug_objs/*.o | sort | uniq`  
	ranlib libscifor.a
	ranlib libscifor_deb.a
	@mv -fv *.a $(SFLIB)/
	@echo 'LIBRARY IS UPDATED.........'
	@echo ' '

SPLINE: 
	$(FC) -c $(FFLAG) $@.f90 -I./mods/ -o objs/$@.o $(MOPT)./mods/
	$(FC) -c $(DFLAG) $@.f90 -I./debug_mods/ -o debug_objs/$@.o $(MOPT)./debug_mods/
mod_SPLINE: SPLINE 
	@rsync -avP mods/*.mod $(SFINCLUDE)/
	@rsync -avP debug_mods/*.mod $(SFINCLUDE)/debug/
	ar cvq libscifor.a `ls objs/*.o | sort | uniq`  
	ar cvq libscifor_deb.a `ls debug_objs/*.o | sort | uniq`  
	ranlib libscifor.a
	ranlib libscifor_deb.a
	@mv -fv *.a $(SFLIB)/
	@echo 'LIBRARY IS UPDATED.........'
	@echo ' '

RANDOM: 
	$(FC) -c $(FFLAG) $@.f90 -I./mods/ -o objs/$@.o $(MOPT)./mods/
	$(FC) -c $(DFLAG) $@.f90 -I./debug_mods/ -o debug_objs/$@.o $(MOPT)./debug_mods/
mod_RANDOM: RANDOM 
	@rsync -avP mods/*.mod $(SFINCLUDE)/
	@rsync -avP debug_mods/*.mod $(SFINCLUDE)/debug/
	ar cvq libscifor.a `ls objs/*.o | sort | uniq`  
	ar cvq libscifor_deb.a `ls debug_objs/*.o | sort | uniq`  
	ranlib libscifor.a
	ranlib libscifor_deb.a
	@mv -fv *.a $(SFLIB)/
	@echo 'LIBRARY IS UPDATED.........'
	@echo ' '

STATISTICS: 
	$(FC) -c $(FFLAG) $@.f90 -I./mods/ -o objs/$@.o $(MOPT)./mods/
	$(FC) -c $(DFLAG) $@.f90 -I./debug_mods/ -o debug_objs/$@.o $(MOPT)./debug_mods/
mod_STATISTICS: STATISTICS 
	@rsync -avP mods/*.mod $(SFINCLUDE)/
	@rsync -avP debug_mods/*.mod $(SFINCLUDE)/debug/
	ar cvq libscifor.a `ls objs/*.o | sort | uniq`  
	ar cvq libscifor_deb.a `ls debug_objs/*.o | sort | uniq`  
	ranlib libscifor.a
	ranlib libscifor_deb.a
	@mv -fv *.a $(SFLIB)/
	@echo 'LIBRARY IS UPDATED.........'
	@echo ' '

TOOLS:  
	$(FC) -c $(FFLAG) $@.f90 -I./mods/ -o objs/$@.o $(MOPT)./mods/
	$(FC) -c $(DFLAG) $@.f90 -I./debug_mods/ -o debug_objs/$@.o $(MOPT)./debug_mods/
mod_TOOLS: COMVARS PARSECMD TIMER IOTOOLS TOOLS 
	@rsync -avP mods/*.mod $(SFINCLUDE)/
	@rsync -avP debug_mods/*.mod $(SFINCLUDE)/debug/
	ar cvq libscifor.a `ls objs/*.o | sort | uniq`  
	ar cvq libscifor_deb.a `ls debug_objs/*.o | sort | uniq`  
	ranlib libscifor.a
	ranlib libscifor_deb.a
	@mv -fv *.a $(SFLIB)/
	@echo 'LIBRARY IS UPDATED.........'
	@echo ' '

PADE:
	$(FC) -c $(FFLAG) $@.f90 -I./mods/ -o objs/$@.o $(MOPT)./mods/
	$(FC) -c $(DFLAG) $@.f90 -I./debug_mods/ -o debug_objs/$@.o $(MOPT)./debug_mods/
mod_PADE: PADE 
	@rsync -avP mods/*.mod $(SFINCLUDE)/
	@rsync -avP debug_mods/*.mod $(SFINCLUDE)/debug/
	ar cvq libscifor.a `ls objs/*.o | sort | uniq`  
	ar cvq libscifor_deb.a `ls debug_objs/*.o | sort | uniq`  
	ranlib libscifor.a
	ranlib libscifor_deb.a
	@mv -fv *.a $(SFLIB)/
	@echo 'LIBRARY IS UPDATED.........'
	@echo ' '

IOTOOLS: 
	$(FC) -c $(FFLAG) IOFILE.f90 -I./mods/ -o objs/IOFILE.o $(MOPT)./mods/	
	$(FC) -c $(FFLAG) SLPLOT.f90 -I./mods/ -o objs/SLPLOT.o $(MOPT)./mods/
	$(FC) -c $(FFLAG) SLREAD.f90 -I./mods/ -o objs/SLREAD.o $(MOPT)./mods/
	$(FC) -c $(FFLAG) $@.f90 -I./mods/ -o objs/$@.o $(MOPT)./mods/
	$(FC) -c $(DFLAG) IOFILE.f90 -I./debug_mods/ -o debug_objs/IOFILE.o $(MOPT)./debug_mods/	
	$(FC) -c $(DFLAG) SLPLOT.f90 -I./debug_mods/ -o debug_objs/SLPLOT.o $(MOPT)./debug_mods/
	$(FC) -c $(DFLAG) SLREAD.f90 -I./debug_mods/ -o debug_objs/SLREAD.o $(MOPT)./debug_mods/
	$(FC) -c $(DFLAG) $@.f90 -I./debug_mods/ -o debug_objs/$@.o $(MOPT)./debug_mods/
mod_IOTOOLS: COMVARS IOTOOLS 
	@rsync -avP mods/*.mod $(SFINCLUDE)/
	@rsync -avP debug_mods/*.mod $(SFINCLUDE)/debug/
	ar cvq libscifor.a `ls objs/*.o | sort | uniq`  
	ar cvq libscifor_deb.a `ls debug_objs/*.o | sort | uniq`  
	ranlib libscifor.a
	ranlib libscifor_deb.a
	@mv -fv *.a $(SFLIB)/
	@echo 'LIBRARY IS UPDATED.........'
	@echo ' '


LIST: 
	$(FC) -c $(FFLAG) LIST_D_ORDERED.f90 -I./mods/ -o objs/LIST_D_ORDERED.o $(MOPT)./mods/	
	$(FC) -c $(FFLAG) LIST_D_UNORDERED.f90 -I./mods/ -o objs/LIST_D_UNORDERED.o $(MOPT)./mods/	
	$(FC) -c $(FFLAG) LIST_Z_UNORDERED.f90 -I./mods/ -o objs/LIST_Z_UNORDERED.o $(MOPT)./mods/	
	$(FC) -c $(DFLAG) LIST_D_ORDERED.f90 -I./debug_mods/ -o debug_objs/LIST_D_ORDERED.o $(MOPT)./debug_mods/	
	$(FC) -c $(DFLAG) LIST_D_UNORDERED.f90 -I./debug_mods/ -o debug_objs/LIST_D_UNORDERED.o $(MOPT)./debug_mods/	
	$(FC) -c $(DFLAG) LIST_Z_UNORDERED.f90 -I./debug_mods/ -o debug_objs/LIST_Z_UNORDERED.o $(MOPT)./debug_mods/	
mod_LIST: LIST 
	@rsync -avP mods/*.mod $(SFINCLUDE)/
	@rsync -avP debug_mods/*.mod $(SFINCLUDE)/debug/
	ar cvq libscifor.a `ls objs/*.o | sort | uniq`  
	ar cvq libscifor_deb.a `ls debug_objs/*.o | sort | uniq`  
	ranlib libscifor.a
	ranlib libscifor_deb.a
	@mv -fv *.a $(SFLIB)/
	@echo 'LIBRARY IS UPDATED.........'
	@echo ' '


SCIFOR: 
	$(FC) -c $(FFLAG) $@.f90 -I./mods/ -o objs/$@.o $(MOPT)./mods/
	$(FC) -c $(DFLAG) $@.f90 -I./debug_mods/ -o debug_objs/$@.o $(MOPT)./debug_mods/
mod_SCIFOR:  
	@rsync -avP mods/*.mod $(SFINCLUDE)/
	@rsync -avP debug_mods/*.mod $(SFINCLUDE)/debug/
	ar cvq libscifor.a `ls objs/*.o | sort | uniq`  
	ar cvq libscifor_deb.a `ls debug_objs/*.o | sort | uniq`  
	ranlib libscifor.a
	ranlib libscifor_deb.a
	@mv -fv *.a $(SFLIB)/
	@echo 'LIBRARY IS UPDATED.........'
	@echo ' '
