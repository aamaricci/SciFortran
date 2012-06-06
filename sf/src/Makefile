include $(SFDIR)/etc/opt.mk
REV = $(shell git rev-parse HEAD)

all: allmod libscifor.a libscifor_deb.a version

allmod: \
COMVARS\
BROYDEN\
BRENT\
CGFIT\
CHRONOBAR\
FFTGF\
GREENFUNX\
INTEGRATE\
IOTOOLS\
SQUARE_LATTICE\
LIST\
MATRIX\
SPLINE\
RANDOM\
STATISTICS\
TOOLS\
STRINGPACK\
PADE

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
	$(FC) -c $(FFLAG) $@.f90 -I$(SFINCLUDE) -o objs/$@.o -module mods/
	$(FC) -c $(DFLAG) $@.f90 -I$(SFINCLUDE) -o debug_objs/$@.o -module ./debug_mods/
mod_comvars: COMVARS
	@rsync -avP mods/*.mod $(SFINCLUDE)/
	@rsync -avP debug_mods/*.mod $(SFINCLUDE)/debug/
	ar cvq libscifor.a `ls objs/*.o | sort | uniq`  
	ar cvq libscifor_deb.a `ls debug_objs/*.o | sort | uniq`  
	ranlib libscifor.a
	ranlib libscifor_deb.a
	@mv -fv *.a $(SFLIB)/
	@echo 'LIBRARY IS UPDATED.........'
	@echo ' '

CGFIT: 
	$(FC) -c $(FFLAG) $@.f90 -I$(SFINCLUDE) -o objs/$@.o -module mods/
	$(FC) -c $(DFLAG) $@.f90 -I$(SFINCLUDE) -o debug_objs/$@.o -module ./debug_mods/
mod_CGFIT: CGFIT
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
	$(FC) -c $(FFLAG) $@.f90 -I$(SFINCLUDE) -o objs/$@.o -module mods/
	$(FC) -c $(DFLAG) $@.f90 -I$(SFINCLUDE) -o debug_objs/$@.o -module ./debug_mods/
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
	$(FC) -c $(FFLAG) $@.f90 -I$(SFINCLUDE) -o objs/$@.o -module mods/
	$(FC) -c $(DFLAG) $@.f90 -I$(SFINCLUDE) -o debug_objs/$@.o -module ./debug_mods/
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

CHRONOBAR: 
	$(FC) -c $(FFLAG) $@.f90 -I$(SFINCLUDE) -o objs/$@.o -module mods/
	$(FC) -c $(DFLAG) $@.f90 -I$(SFINCLUDE) -o debug_objs/$@.o -module ./debug_mods/
mod_CHRONOBAR: CHRONOBAR 
	@rsync -avP mods/*.mod $(SFINCLUDE)/
	@rsync -avP debug_mods/*.mod $(SFINCLUDE)/debug/
	ar cvq libscifor.a `ls objs/*.o | sort | uniq`  
	ar cvq libscifor_deb.a `ls debug_objs/*.o | sort | uniq`  
	ranlib libscifor.a
	ranlib libscifor_deb.a
	@mv -fv *.a $(SFLIB)/
	@echo 'LIBRARY IS UPDATED.........'
	@echo ' '

FFTGF: 
	$(FC) -c $(FFLAG) $@.f90 -I$(SFINCLUDE) -o objs/$@.o -module mods/
	$(FC) -c $(DFLAG) $@.f90 -I$(SFINCLUDE) -o debug_objs/$@.o -module ./debug_mods/
mod_fftgf: FFTGF 
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
	$(FC) -c $(FFLAG) $@.f90 -I$(SFINCLUDE) -o objs/$@.o -module mods/
	$(FC) -c $(DFLAG) $@.f90 -I$(SFINCLUDE) -o debug_objs/$@.o -module ./debug_mods/
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
	$(FC) -c $(FFLAG) $@.f90 -I$(SFINCLUDE) -o objs/$@.o -module mods/
	$(FC) -c $(DFLAG) $@.f90 -I$(SFINCLUDE) -o debug_objs/$@.o -module ./debug_mods/
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
	$(FC) -c $(FFLAG) $@.f90 -I$(SFINCLUDE) -o objs/$@.o -module mods/
	$(FC) -c $(DFLAG) $@.f90 -I$(SFINCLUDE) -o debug_objs/$@.o -module ./debug_mods/
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
	$(FC) -c $(FFLAG) $@.f90 -I$(SFINCLUDE) -o objs/$@.o -module mods/
	$(FC) -c $(DFLAG) $@.f90 -I$(SFINCLUDE) -o debug_objs/$@.o -module ./debug_mods/
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
	$(FC) -c $(FFLAG) $@.f90 -I$(SFINCLUDE) -o objs/$@.o -module mods/
	$(FC) -c $(DFLAG) $@.f90 -I$(SFINCLUDE) -o debug_objs/$@.o -module ./debug_mods/
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
	$(FC) -c $(FFLAG) $@.f90 -I$(SFINCLUDE) -o objs/$@.o -module mods/
	$(FC) -c $(DFLAG) $@.f90 -I$(SFINCLUDE) -o debug_objs/$@.o -module ./debug_mods/
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
	$(FC) -c $(FFLAG) $@.f90 -I$(SFINCLUDE) -o objs/$@.o -module mods/
	$(FC) -c $(DFLAG) $@.f90 -I$(SFINCLUDE) -o debug_objs/$@.o -module ./debug_mods/
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
	$(FC) -c $(FFLAG) $@.f90 -I$(SFINCLUDE) -o objs/$@.o -module mods/
	$(FC) -c $(DFLAG) $@.f90 -I$(SFINCLUDE) -o debug_objs/$@.o -module ./debug_mods/
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
	$(FC) -c $(FFLAG) $@.f90 -I$(SFINCLUDE) -o objs/$@.o -module mods/
	$(FC) -c $(DFLAG) $@.f90 -I$(SFINCLUDE) -o debug_objs/$@.o -module ./debug_mods/
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
	$(FC) -c $(FFLAG) IOFILE.f90 -I$(SFINCLUDE) -o objs/IOFILE.o -module mods/	
	$(FC) -c $(FFLAG) SLPLOT.f90 -I$(SFINCLUDE) -o objs/SLPLOT.o -module mods/
	$(FC) -c $(FFLAG) SLREAD.f90 -I$(SFINCLUDE) -o objs/SLREAD.o -module mods/
	$(FC) -c $(FFLAG) $@.f90 -I$(SFINCLUDE) -o objs/$@.o -module mods/
	$(FC) -c $(DFLAG) IOFILE.f90 -I$(SFINCLUDE) -o debug_objs/IOFILE.o -module ./debug_mods/	
	$(FC) -c $(DFLAG) SLPLOT.f90 -I$(SFINCLUDE) -o debug_objs/SLPLOT.o -module ./debug_mods/
	$(FC) -c $(DFLAG) SLREAD.f90 -I$(SFINCLUDE) -o debug_objs/SLREAD.o -module ./debug_mods/
	$(FC) -c $(DFLAG) $@.f90 -I$(SFINCLUDE) -o debug_objs/$@.o -module ./debug_mods/
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
	$(FC) -c $(FFLAG) D_ORDERED_LIST.f90 -I$(SFINCLUDE) -o objs/D_ORDERED_LIST.o -module mods/	
	$(FC) -c $(FFLAG) D_UNORDERED_LIST.f90 -I$(SFINCLUDE) -o objs/D_UNORDERED_LIST.o -module mods/	
	$(FC) -c $(FFLAG) Z_UNORDERED_LIST.f90 -I$(SFINCLUDE) -o objs/Z_UNORDERED_LIST.o -module mods/	
	$(FC) -c $(DFLAG) D_ORDERED_LIST.f90 -I$(SFINCLUDE) -o debug_objs/D_ORDERED_LIST.o -module ./debug_mods/	
	$(FC) -c $(DFLAG) D_UNORDERED_LIST.f90 -I$(SFINCLUDE) -o debug_objs/D_UNORDERED_LIST.o -module ./debug_mods/	
	$(FC) -c $(DFLAG) Z_UNORDERED_LIST.f90 -I$(SFINCLUDE) -o debug_objs/Z_UNORDERED_LIST.o -module ./debug_mods/	
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

STRINGPACK: 
	$(FC) -c $(FFLAG) STRING_FUNX_1.f90 -I$(SFINCLUDE) -o objs/STRING_FUNX_1.o -module mods/	
	$(FC) -c $(FFLAG) STRING_FUNX_2.f90 -I$(SFINCLUDE) -o objs/STRING_FUNX_2.o -module mods/	
	$(FC) -c $(FFLAG) $@.f90 -I$(SFINCLUDE) -o objs/$@.o -module mods/	
	$(FC) -c $(DFLAG) STRING_FUNX_1.f90 -I$(SFINCLUDE) -o debug_objs/STRING_FUNX_1.o -module ./debug_mods/	
	$(FC) -c $(DFLAG) STRING_FUNX_2.f90 -I$(SFINCLUDE) -o debug_objs/STRING_FUNX_2.o -module ./debug_mods/	
	$(FC) -c $(DFLAG) STRINGPACK.f90 -I$(SFINCLUDE) -o debug_objs/$@.o -module ./debug_mods/
mod_STRINGPACK: STRINGPACK 
	@rsync -avP mods/*.mod $(SFINCLUDE)/
	@rsync -avP debug_mods/*.mod $(SFINCLUDE)/debug/
	ar cvq libscifor.a `ls objs/*.o | sort | uniq`  
	ar cvq libscifor_deb.a `ls debug_objs/*.o | sort | uniq`  
	ranlib libscifor.a
	ranlib libscifor_deb.a
	@mv -fv *.a $(SFLIB)/
	@echo 'LIBRARY IS UPDATED.........'
	@echo ' '


#DEB

clean: 
	@rm -fv *.o *.mod *.a mods/* objs/* debug_*/*
