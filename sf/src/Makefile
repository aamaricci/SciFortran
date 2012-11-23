include opt.mk
REV = $(shell git rev-parse HEAD)

all: allmod libscifor.a libscifor_deb.a version

allmod: COMVARS PARSECMD TIMER TOOLS BROYDEN BRENT FFTGF GREENFUNX INTEGRATE IOTOOLS SQUARE_LATTICE LIST MATRIX SPLINE RANDOM STATISTICS PADE

allmod_deb: COMVARS_DEB PARSECMD_DEB TIMER_DEB TOOLS_DEB BROYDEN_DEB BRENT_DEB  GREENFUNX_DEB  FFTGF_DEB INTEGRATE_DEB IOTOOLS_DEB SQUARE_LATTICE_DEB LIST_DEB MATRIX_DEB SPLINE_DEB RANDOM_DEB STATISTICS_DEB  PADE_DEB

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




#STD
COMVARS: 
	$(FC) -c $(FFLAG) $@.f90 -I./mods/ -o objs/$@.o $(MOPT) mods/
	$(FC) -c $(DFLAG) $@.f90 -I./debug_mods/ -o debug_objs/$@.o $(MOPT) ./debug_mods/
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
	$(FC) -c $(FFLAG) $@.f90 -I./mods/ -o objs/$@.o $(MOPT) mods/
	$(FC) -c $(DFLAG) $@.f90 -I./debug_mods/ -o debug_objs/$@.o $(MOPT) ./debug_mods/
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
	$(FC) -c $(FFLAG) $@.f90 -I./mods/ -o objs/$@.o $(MOPT) mods/
	$(FC) -c $(DFLAG) $@.f90 -I./debug_mods/ -o debug_objs/$@.o $(MOPT) ./debug_mods/
mod_BROYDEN: BROYDEN 
	@rsync -avP mods/*.mod $(SFINCLUDE)/
	@rsync -avP debug_mods/*.mod $(SFINCLUDE)/debug/
	ar cvq libscifor.a `ls objs/*.o | sort | uniq`  
	ar cvq libscifor_deb.a `ls debug_objs/*.o | sort | uniq`  
	ranlib libscifor.a
	ranlib libscifor_deb.a
	@mv -fv *.a $(SFLIB)/
	@echo 'LIBRARY IS UPDATED.........'
	@echo ' '


BRENT: 
	$(FC) -c $(FFLAG) $@.f90 -I./mods/ -o objs/$@.o $(MOPT) mods/
	$(FC) -c $(DFLAG) $@.f90 -I./debug_mods/ -o debug_objs/$@.o $(MOPT) ./debug_mods/
mod_BRENT: BRENT 
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
	$(FC) -c $(FFLAG) $@.f90 -I./mods/ -o objs/$@.o $(MOPT) mods/
	$(FC) -c $(DFLAG) $@.f90 -I./debug_mods/ -o debug_objs/$@.o $(MOPT) ./debug_mods/
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
	$(FC) -c $(FFLAG) $@.f90 -I./mods/ -o objs/$@.o $(MOPT) mods/
	$(FC) -c $(DFLAG) $@.f90 -I./debug_mods/ -o debug_objs/$@.o $(MOPT) ./debug_mods/
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


GREENFUNX: 
	$(FC) -c $(FFLAG) $@.f90 -I./mods/ -o objs/$@.o $(MOPT) mods/
	$(FC) -c $(DFLAG) $@.f90 -I./debug_mods/ -o debug_objs/$@.o $(MOPT) ./debug_mods/
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
	$(FC) -c $(FFLAG) $@.f90 -I./mods/ -o objs/$@.o $(MOPT) mods/
	$(FC) -c $(DFLAG) $@.f90 -I./debug_mods/ -o debug_objs/$@.o $(MOPT) ./debug_mods/
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
	$(FC) -c $(FFLAG) $@.f90 -I./mods/ -o objs/$@.o $(MOPT) mods/
	$(FC) -c $(DFLAG) $@.f90 -I./debug_mods/ -o debug_objs/$@.o $(MOPT) ./debug_mods/
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
	$(FC) -c $(FFLAG) $@.f90 -I./mods/ -o objs/$@.o $(MOPT) mods/
	$(FC) -c $(DFLAG) $@.f90 -I./debug_mods/ -o debug_objs/$@.o $(MOPT) ./debug_mods/
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
	$(FC) -c $(FFLAG) $@.f90 -I./mods/ -o objs/$@.o $(MOPT) mods/
	$(FC) -c $(DFLAG) $@.f90 -I./debug_mods/ -o debug_objs/$@.o $(MOPT) ./debug_mods/
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
	$(FC) -c $(FFLAG) $@.f90 -I./mods/ -o objs/$@.o $(MOPT) mods/
	$(FC) -c $(DFLAG) $@.f90 -I./debug_mods/ -o debug_objs/$@.o $(MOPT) ./debug_mods/
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
	$(FC) -c $(FFLAG) $@.f90 -I./mods/ -o objs/$@.o $(MOPT) mods/
	$(FC) -c $(DFLAG) $@.f90 -I./debug_mods/ -o debug_objs/$@.o $(MOPT) ./debug_mods/
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
	$(FC) -c $(FFLAG) $@.f90 -I./mods/ -o objs/$@.o $(MOPT) mods/
	$(FC) -c $(DFLAG) $@.f90 -I./debug_mods/ -o debug_objs/$@.o $(MOPT) ./debug_mods/
mod_TOOLS: TOOLS 
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
	$(FC) -c $(FFLAG) $@.f90 -I./mods/ -o objs/$@.o $(MOPT) mods/
	$(FC) -c $(DFLAG) $@.f90 -I./debug_mods/ -o debug_objs/$@.o $(MOPT) ./debug_mods/
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
	$(FC) -c $(FFLAG) IOFILE.f90 -I./mods/ -o objs/IOFILE.o $(MOPT) mods/	
	$(FC) -c $(FFLAG) SLPLOT.f90 -I./mods/ -o objs/SLPLOT.o $(MOPT) mods/
	$(FC) -c $(FFLAG) SLREAD.f90 -I./mods/ -o objs/SLREAD.o $(MOPT) mods/
	$(FC) -c $(FFLAG) $@.f90 -I./mods/ -o objs/$@.o $(MOPT) mods/
	$(FC) -c $(DFLAG) IOFILE.f90 -I./debug_mods/ -o debug_objs/IOFILE.o $(MOPT) ./debug_mods/	
	$(FC) -c $(DFLAG) SLPLOT.f90 -I./debug_mods/ -o debug_objs/SLPLOT.o $(MOPT) ./debug_mods/
	$(FC) -c $(DFLAG) SLREAD.f90 -I./debug_mods/ -o debug_objs/SLREAD.o $(MOPT) ./debug_mods/
	$(FC) -c $(DFLAG) $@.f90 -I./debug_mods/ -o debug_objs/$@.o $(MOPT) ./debug_mods/
mod_IOTOOLS: IOTOOLS 
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
	$(FC) -c $(FFLAG) D_ORDERED_LIST.f90 -I./mods/ -o objs/D_ORDERED_LIST.o $(MOPT) mods/	
	$(FC) -c $(FFLAG) D_UNORDERED_LIST.f90 -I./mods/ -o objs/D_UNORDERED_LIST.o $(MOPT) mods/	
	$(FC) -c $(FFLAG) Z_UNORDERED_LIST.f90 -I./mods/ -o objs/Z_UNORDERED_LIST.o $(MOPT) mods/	
	$(FC) -c $(DFLAG) D_ORDERED_LIST.f90 -I./debug_mods/ -o debug_objs/D_ORDERED_LIST.o $(MOPT) ./debug_mods/	
	$(FC) -c $(DFLAG) D_UNORDERED_LIST.f90 -I./debug_mods/ -o debug_objs/D_UNORDERED_LIST.o $(MOPT) ./debug_mods/	
	$(FC) -c $(DFLAG) Z_UNORDERED_LIST.f90 -I./debug_mods/ -o debug_objs/Z_UNORDERED_LIST.o $(MOPT) ./debug_mods/	
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

# STRINGPACK: 
# 	$(FC) -c $(FFLAG) STRING_FUNX_1.f90 -I./mods/ -o objs/STRING_FUNX_1.o $(MOPT) mods/	
# 	$(FC) -c $(FFLAG) STRING_FUNX_2.f90 -I./mods/ -o objs/STRING_FUNX_2.o $(MOPT) mods/	
# 	$(FC) -c $(FFLAG) $@.f90 -I./mods/ -o objs/$@.o $(MOPT) mods/	
# 	$(FC) -c $(DFLAG) STRING_FUNX_1.f90 -I./debug_mods/ -o debug_objs/STRING_FUNX_1.o $(MOPT) ./debug_mods/	
# 	$(FC) -c $(DFLAG) STRING_FUNX_2.f90 -I./debug_mods/ -o debug_objs/STRING_FUNX_2.o $(MOPT) ./debug_mods/	
# 	$(FC) -c $(DFLAG) STRINGPACK.f90 -I./debug_mods/ -o debug_objs/$@.o $(MOPT) ./debug_mods/
# mod_STRINGPACK: STRINGPACK 
# 	@rsync -avP mods/*.mod $(SFINCLUDE)/
# 	@rsync -avP debug_mods/*.mod $(SFINCLUDE)/debug/
# 	ar cvq libscifor.a `ls objs/*.o | sort | uniq`  
# 	ar cvq libscifor_deb.a `ls debug_objs/*.o | sort | uniq`  
# 	ranlib libscifor.a
# 	ranlib libscifor_deb.a
# 	@mv -fv *.a $(SFLIB)/
# 	@echo 'LIBRARY IS UPDATED.........'
# 	@echo ' '

clean: 
	@rm -fv *.o *.mod *.a mods/* objs/* debug_*/* $(SFINCLUDE)/*.mod $(SFINCLUDE)/debug/*.mod
