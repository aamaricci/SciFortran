include options.inc
REV = $(shell git rev-parse HEAD)
DIR = ./scifor_portable

all: version allmod libscifor.a libscifor_deb.a

allmod: SCIFOR_VERSION COMVARS PARSE_INPUT TIMER IOTOOLS FFTGF TOOLS ERROR ARRAYS  DERIVATE OPTIMIZE FUNCTIONS GREENFUNX INTEGRATE VECTORS SQUARE_LATTICE LIST MATRIX INTERPOLATE RANDOM STATISTICS SCIFOR


libscifor.a: 
	@rsync -avP mods/*.mod $(SFMODS)/
	ar cvq $@ `ls objs/*.o | sort | uniq`  
	ranlib $@
	@mv -fv *.a $(SFLIB)/
	@echo 'LIBRARY IS DONE.........'
	@echo ' '

libscifor_deb.a: 
	@rsync -avP debug_mods/*.mod $(SFMODS)/debug/
	ar cvq $@ `ls debug_objs/*.o | sort | uniq`  
	ranlib $@
	@mv -fv *.a $(SFLIB)/
	@echo 'LIBRARY IS DONE.........'
	@echo ' '

version:
	@echo "Actual version is: $(REV)"
	@echo 'character(len=41),parameter,public :: sf_version = "$(REV)"' > scifor_version.inc


clean: 
	@rm -fv *.o *.mod *.a mods/* objs/* debug_*/* $(SFMODS)/*.mod $(SFMODS)/debug/*.mod $(SFLIB)/*scifor*.a

portable:
	./create_portable.sh $(DIR)

libupdate:
	@rsync -avP mods/*.mod $(SFMODS)/
	@rsync -avP debug_mods/*.mod $(SFMODS)/debug/
	ar cvq libscifor.a `ls objs/*.o | sort | uniq`  
	ar cvq libscifor_deb.a `ls debug_objs/*.o | sort | uniq`  
	ranlib libscifor.a
	ranlib libscifor_deb.a
	@mv -fv *.a $(SFLIB)/
	@echo 'LIBRARY IS UPDATED.........'
	@echo ' '



SCIFOR_VERSION: version
	$(FC) -c $(FFLAG) $@.f90 -I./mods/ -o objs/$@.o $(MOPT)./mods/
	$(FC) -c $(DFLAG) $@.f90 -I./debug_mods/ -o debug_objs/$@.o $(MOPT)./debug_mods/
mod_SCIFOR_VERSION: SCIFOR_VERSION libupdate


COMVARS: 
	$(FC) -c $(FFLAG) $@.f90 -I./mods/ -o objs/$@.o $(MOPT)./mods/
	$(FC) -c $(DFLAG) $@.f90 -I./debug_mods/ -o debug_objs/$@.o $(MOPT)./debug_mods/
mod_COMVARS: COMVARS libupdate

PARSE_INPUT: 
	$(FC) -c $(FFLAG)  PARSE_LIST_INPUT.f90 -o objs/PARSE_LIST_INPUT.o $(MOPT)./mods/
	$(FC) -c $(DFLAG)  PARSE_LIST_INPUT.f90 -o debug_objs/PARSE_LIST_INPUT.o.o $(MOPT)./debug_mods/
	$(FC) -c $(FFLAG) $@.f90 -I./mods/ -o objs/$@.o $(MOPT)./mods/
	$(FC) -c $(DFLAG) $@.f90 -I./debug_mods/ -o debug_objs/$@.o $(MOPT)./debug_mods/
mod_PARSE_INPUT: PARSE_INPUT libupdate


TIMER: 
	$(FC) -c $(FFLAG) $@.f90 -I./mods/ -o objs/$@.o $(MOPT)./mods/
	$(FC) -c $(DFLAG) $@.f90 -I./debug_mods/ -o debug_objs/$@.o $(MOPT)./debug_mods/
mod_TIMER: TIMER libupdate 

IOTOOLS: 
	$(FC) -c $(FFLAG) IOFILE.f90 -I./mods/ -o objs/IOFILE.o $(MOPT)./mods/	
	$(FC) -c $(FFLAG) IOPLOT.f90 -I./mods/ -o objs/IOPLOT.o $(MOPT)./mods/
	$(FC) -c $(FFLAG) IOREAD.f90 -I./mods/ -o objs/IOREAD.o $(MOPT)./mods/
	$(FC) -c $(FFLAG) $@.f90 -I./mods/ -o objs/$@.o $(MOPT)./mods/
	$(FC) -c $(DFLAG) IOFILE.f90 -I./debug_mods/ -o debug_objs/IOFILE.o $(MOPT)./debug_mods/	
	$(FC) -c $(DFLAG) IOPLOT.f90 -I./debug_mods/ -o debug_objs/IOPLOT.o $(MOPT)./debug_mods/
	$(FC) -c $(DFLAG) IOREAD.f90 -I./debug_mods/ -o debug_objs/IOREAD.o $(MOPT)./debug_mods/
	$(FC) -c $(DFLAG) $@.f90 -I./debug_mods/ -o debug_objs/$@.o $(MOPT)./debug_mods/
mod_IOTOOLS: COMVARS IOTOOLS libupdate


FFTGF:  COMVARS INTERPOLATE
	$(FC) -c $(FFLAG) $@.f90 -I./mods/ -o objs/$@.o $(MOPT)./mods/
	$(FC) -c $(DFLAG) $@.f90 -I./debug_mods/ -o debug_objs/$@.o $(MOPT)./debug_mods/
mod_FFTGF: FFTGF libupdate 


TOOLS:  
	$(FC) -c $(FFLAG) $@.f90 -I./mods/ -o objs/$@.o $(MOPT)./mods/
	$(FC) -c $(DFLAG) $@.f90 -I./debug_mods/ -o debug_objs/$@.o $(MOPT)./debug_mods/
mod_TOOLS: COMVARS TIMER FFTGF TOOLS libupdate


ERROR:  
	$(FC) -c $(FFLAG) $@.f90 -I./mods/ -o objs/$@.o $(MOPT)./mods/
	$(FC) -c $(DFLAG) $@.f90 -I./debug_mods/ -o debug_objs/$@.o $(MOPT)./debug_mods/
mod_ERROR: ERROR libupdate


ARRAYS:  
	$(FC) -c $(FFLAG) $@.f90 -I./mods/ -o objs/$@.o $(MOPT)./mods/
	$(FC) -c $(DFLAG) $@.f90 -I./debug_mods/ -o debug_objs/$@.o $(MOPT)./debug_mods/
mod_ARRAYS: ARRAYS libupdate


DERIVATE: 
	$(FC) -c $(FFLAG) $@.f90 -I./mods/ -o objs/$@.o $(MOPT)./mods/
	$(FC) -c $(DFLAG)  $@.f90 -I./debug_mods/ -o debug_objs/$@.o $(MOPT)./debug_mods/
mod_DERIVATE: DERIVATE libupdate


OPTIMIZE_ROOT_FINDING: 
	$(FC) -c $(FFLAG)  optimize_broyden_routines.f90 -o objs/optimize_broyden_routines.o $(MOPT)./mods/
	$(FC) -c $(DFLAG)  optimize_broyden_routines.f90 -o debug_objs/optimize_broyden_routines.o $(MOPT)./debug_mods/
# $(FC) -c $(FFLAG)  minpack.f90 -o objs/minpack.o $(MOPT)./mods/
# $(FC) -c $(DFLAG)  minpack.f90 -o debug_objs/minpack.o $(MOPT)./debug_mods/
	$(FC) -c $(FFLAG) -fpp -D_$(PRECOMP_MINPACK) $@.f90 -I./mods/ -o objs/$@.o $(MOPT)./mods/
	$(FC) -c $(DFLAG) -fpp -D_$(PRECOMP_MINPACK) $@.f90 -I./debug_mods/ -o debug_objs/$@.o $(MOPT)./debug_mods/
OPTIMIZE_MINIMIZE: 
	$(FC) -c $(FFLAG)  minimize.f -o objs/minimize.o $(MOPT)./mods/
	$(FC) -c $(FFLAG)  CGFAM.f -o objs/CGFAM.o $(MOPT)./mods/
	$(FC) -c $(FFLAG)  CGSEARCH.f -o objs/CGSEARCH.o $(MOPT)./mods/
	$(FC) -c $(FFLAG)  CGBLAS.f -o objs/CGBLAS.o $(MOPT)./mods/
	$(FC) -c $(FFLAG)  optimize_cgfit_routines.f90 -I./mods/ -o objs/optimize_cgfit_routines.o $(MOPT)./mods/ 
	$(FC) -c $(FFLAG)  $@.f90 -I./mods/  -o objs/$@.o $(MOPT)./mods/
	$(FC) -c $(DFLAG)  minimize.f -o debug_objs/minimize.o $(MOPT)./debug_mods/
	$(FC) -c $(DFLAG)  CGFAM.f -o debug_objs/CGFAM.o $(MOPT)./debug_mods/
	$(FC) -c $(DFLAG)  CGSEARCH.f -o debug_objs/CGSEARCH.o $(MOPT)./debug_mods/
	$(FC) -c $(DFLAG)  CGBLAS.f -o debug_objs/CGBLAS.o $(MOPT)./debug_mods/
	$(FC) -c $(DFLAG)  optimize_cgfit_routines.f90 -I./debug_mods/ -o debug_objs/optimize_cgfit_routines.o $(MOPT)./debug_mods/
	$(FC) -c $(DFLAG)  $@.f90 -I./debug_mods/ -o debug_objs/$@.o $(MOPT)./debug_mods/
OPTIMIZE: OPTIMIZE_ROOT_FINDING OPTIMIZE_MINIMIZE 
	$(FC) -c $(FFLAG) $@.f90 -I./mods/ -o objs/$@.o $(MOPT)./mods/
	$(FC) -c $(DFLAG) $@.f90 -I./debug_mods/ -o debug_objs/$@.o $(MOPT)./debug_mods/
mod_OPTIMIZE: OPTIMIZE libupdate


FUNCTIONS: 
	$(FC) -c $(FFLAG) $@.f90 -I./mods/ -o objs/$@.o $(MOPT)./mods/
	$(FC) -c $(DFLAG) $@.f90 -I./debug_mods/ -o debug_objs/$@.o $(MOPT)./debug_mods/
mod_FUNCTIONS: FUNCTIONS libupdate


GREENFUNX: 
	$(FC) -c $(FFLAG) $@.f90 -I./mods/ -o objs/$@.o $(MOPT)./mods/
	$(FC) -c $(DFLAG) $@.f90 -I./debug_mods/ -o debug_objs/$@.o $(MOPT)./debug_mods/
mod_GREENFUNX: GREENFUNX libupdate


INTEGRATE: 
# $(FC) -c $(FFLAG) integrate_quadpack.f90 -I./mods/ -o objs/integrate_quadpack.o $(MOPT)./mods/
# $(FC) -c $(DFLAG) integrate_quadpack.f90 -I./debug_mods/ -o debug_objs/integrate_quadpack.o $(MOPT)./debug_mods/
	$(FC) -c $(FFLAG)  $@.f90 -I./mods/ -o objs/$@.o $(MOPT)./mods/
	$(FC) -c $(DFLAG)  $@.f90 -I./debug_mods/ -o debug_objs/$@.o $(MOPT)./debug_mods/
mod_INTEGRATE: INTEGRATE libupdate


VECTORS: 
	$(FC) -c $(FFLAG) $@.f90 -I./mods/ -o objs/$@.o $(MOPT)./mods/
	$(FC) -c $(DFLAG) $@.f90 -I./debug_mods/ -o debug_objs/$@.o $(MOPT)./debug_mods/
mod_VECTORS: VECTORS libupdate


SQUARE_LATTICE: 
	$(FC) -c $(FFLAG) $@.f90 -I./mods/ -o objs/$@.o $(MOPT)./mods/
	$(FC) -c $(DFLAG) $@.f90 -I./debug_mods/ -o debug_objs/$@.o $(MOPT)./debug_mods/
mod_SQUARE_LATTICE: VECTORS SQUARE_LATTICE libupdate


LIST: 
	$(FC) -c $(FFLAG) LIST_D_ORDERED.f90 -I./mods/ -o objs/LIST_D_ORDERED.o $(MOPT)./mods/	
	$(FC) -c $(FFLAG) LIST_D_UNORDERED.f90 -I./mods/ -o objs/LIST_D_UNORDERED.o $(MOPT)./mods/	
	$(FC) -c $(FFLAG) LIST_Z_UNORDERED.f90 -I./mods/ -o objs/LIST_Z_UNORDERED.o $(MOPT)./mods/	
	$(FC) -c $(DFLAG) LIST_D_ORDERED.f90 -I./debug_mods/ -o debug_objs/LIST_D_ORDERED.o $(MOPT)./debug_mods/	
	$(FC) -c $(DFLAG) LIST_D_UNORDERED.f90 -I./debug_mods/ -o debug_objs/LIST_D_UNORDERED.o $(MOPT)./debug_mods/	
	$(FC) -c $(DFLAG) LIST_Z_UNORDERED.f90 -I./debug_mods/ -o debug_objs/LIST_Z_UNORDERED.o $(MOPT)./debug_mods/	
mod_LIST: LIST libupdate


MATRIX: 
	$(FC) -c $(FFLAG) $@.f90 -I./mods/ -o objs/$@.o $(MOPT)./mods/
	$(FC) -c $(DFLAG) $@.f90 -I./debug_mods/ -o debug_objs/$@.o $(MOPT)./debug_mods/
mod_MATRIX: MATRIX libupdate


INTERPOLATE: 
	$(FC) -c $(FFLAG) $@.f90 -I./mods/ -o objs/$@.o $(MOPT)./mods/
	$(FC) -c $(DFLAG) $@.f90 -I./debug_mods/ -o debug_objs/$@.o $(MOPT)./debug_mods/
mod_INTERPOLATE: INTERPOLATE libupdate


RANDOM: 
	$(FC) -c $(FFLAG) $@.f90 -I./mods/ -o objs/$@.o $(MOPT)./mods/
	$(FC) -c $(DFLAG) $@.f90 -I./debug_mods/ -o debug_objs/$@.o $(MOPT)./debug_mods/
mod_RANDOM: RANDOM libupdate

STATISTICS: 
	$(FC) -c $(FFLAG) $@.f90 -I./mods/ -o objs/$@.o $(MOPT)./mods/
	$(FC) -c $(DFLAG) $@.f90 -I./debug_mods/ -o debug_objs/$@.o $(MOPT)./debug_mods/
mod_STATISTICS: STATISTICS libupdate


SCIFOR: 
	$(FC) -c $(FFLAG) $@.f90 -I./mods/ -o objs/$@.o $(MOPT)./mods/
	$(FC) -c $(DFLAG) $@.f90 -I./debug_mods/ -o debug_objs/$@.o $(MOPT)./debug_mods/
mod_SCIFOR:  SCIFOR libupdate
