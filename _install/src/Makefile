include ../make.inc

REV = $(shell git rev-parse HEAD)

DIR = ./scifor_portable

all: version mkdir allmod libscifor.a


mkdir:
	@mkdir -p $(OBJ_DIR)
	@mkdir -p $(MOD_DIR)

allmod: SCIFOR_VERSION CONSTANTS FONTS TIMER IOTOOLS FFT TOOLS ARRAYS DERIVATE OPTIMIZE FUNCTIONS INTEGRATE LIST MATRIX INTERPOLATE RANDOM STATISTICS SCIFOR


libscifor.a: 
	@rsync -avP $(MOD_DIR)/*.mod $(MOD_SCIFOR)/
	ar crv $(LIB_SCIFOR) `ls $(OBJ_DIR)/*.o | sort | uniq`  
	ranlib $(LIB_SCIFOR)
	@echo 'LIBRARY IS DONE.........'
	@echo ' '


version:
	@echo "Actual version is: $(REV)"
	@echo 'character(len=41),parameter,public :: version_scifor = "$(REV)"' > scifor_version.inc

clean: 
	@rm -fv *.o *.mod *.a

portable:
	./create_portable.sh $(DIR)

libupdate:
	@echo "	rsync -avP $(MOD_DIR)/*.mod $(MOD_SCIFOR)/"
	@rsync -avP $(MOD_DIR)/*.mod $(MOD_SCIFOR)/
	ar crv $(LIB_SCIFOR) `ls $(OBJ_DIR)/*.o | sort | uniq`  
	ranlib $(LIB_SCIFOR)
	@echo 'LIBRARY IS UPDATED.........'
	@echo ' '



SCIFOR_VERSION: version
	$(FC) -c $(FFLAGS) $@.f90 -I$(MOD_DIR)/ -o $(OBJ_DIR)/$@.o $(MOPT)$(MOD_DIR)/
mod_SCIFOR_VERSION: SCIFOR_VERSION libupdate



CONSTANTS: 
	$(FC) -c $(FFLAGS) $@.f90 -I$(MOD_DIR)/ -o $(OBJ_DIR)/$@.o $(MOPT)$(MOD_DIR)/
mod_CONSTANTS: CONSTANTS libupdate



FONTS: 
	$(FC) -c $(FFLAGS) $@.f90 -I$(MOD_DIR)/ -o $(OBJ_DIR)/$@.o $(MOPT)$(MOD_DIR)/
mod_FONTS: FONTS libupdate



TIMER: 
	$(FC) -c $(FFLAGS) $@.f90 -I$(MOD_DIR)/ -o $(OBJ_DIR)/$@.o $(MOPT)$(MOD_DIR)/
mod_TIMER: TIMER libupdate 



IOTOOLS: 
	$(FC) -c $(FFLAGS) IOFILE.f90 -I$(MOD_DIR)/ -o $(OBJ_DIR)/IOFILE.o $(MOPT)$(MOD_DIR)/	
	$(FC) -c $(FFLAGS) IOPLOT.f90 -I$(MOD_DIR)/ -o $(OBJ_DIR)/IOPLOT.o $(MOPT)$(MOD_DIR)/
	$(FC) -c $(FFLAGS) IOREAD.f90 -I$(MOD_DIR)/ -o $(OBJ_DIR)/IOREAD.o $(MOPT)$(MOD_DIR)/
	$(FC) -c $(FFLAGS) $@.f90 -I$(MOD_DIR)/ -o $(OBJ_DIR)/$@.o $(MOPT)$(MOD_DIR)/
mod_IOTOOLS: COMVARS IOTOOLS libupdate



FFT:  CONSTANTS INTERPOLATE
	$(FC) -c $(FFLAGS) $@.f90 -I$(MOD_DIR)/ -o $(OBJ_DIR)/$@.o $(MOPT)$(MOD_DIR)/
mod_FFT: FFT libupdate 



TOOLS:  
	$(FC) -c $(FFLAGS) $@.f90 -I$(MOD_DIR)/ -o $(OBJ_DIR)/$@.o $(MOPT)$(MOD_DIR)/
mod_TOOLS: TOOLS libupdate



ARRAYS:  
	$(FC) -c $(FFLAGS) $@.f90 -I$(MOD_DIR)/ -o $(OBJ_DIR)/$@.o $(MOPT)$(MOD_DIR)/
mod_ARRAYS: ARRAYS libupdate


DERIVATE: 
	$(FC) -c $(FFLAGS) $@.f90 -I$(MOD_DIR)/ -o $(OBJ_DIR)/$@.o $(MOPT)$(MOD_DIR)/
mod_DERIVATE: DERIVATE libupdate



OPTIMIZE_ROOT_FINDING: 
	$(FC) -c $(FFLAGS) optimize_broyden_routines.f90 -o $(OBJ_DIR)/optimize_broyden_routines.o $(MOPT)$(MOD_DIR)/
	$(FC) -c $(FFLAGS) $@.f90 -I$(MOD_DIR)/ -o $(OBJ_DIR)/$@.o $(MOPT)$(MOD_DIR)/
OPTIMIZE_MINIMIZE: 
	$(FC) -c $(FFLAGS) minimize.f -o $(OBJ_DIR)/minimize.o $(MOPT)$(MOD_DIR)/
	$(FC) -c $(FFLAGS) CGFAM.f -o $(OBJ_DIR)/CGFAM.o $(MOPT)$(MOD_DIR)/
	$(FC) -c $(FFLAGS) CGSEARCH.f -o $(OBJ_DIR)/CGSEARCH.o $(MOPT)$(MOD_DIR)/
	$(FC) -c $(FFLAGS) CGBLAS.f -o $(OBJ_DIR)/CGBLAS.o $(MOPT)$(MOD_DIR)/
	$(FC) -c $(FFLAGS) optimize_cgfit_routines.f90 -I$(MOD_DIR)/ -o $(OBJ_DIR)/optimize_cgfit_routines.o $(MOPT)$(MOD_DIR)/ 
	$(FC) -c $(FFLAGS) $@.f90 -I$(MOD_DIR)/  -o $(OBJ_DIR)/$@.o $(MOPT)$(MOD_DIR)/
OPTIMIZE: OPTIMIZE_ROOT_FINDING OPTIMIZE_MINIMIZE 
	$(FC) -c $(FFLAGS) $@.f90 -I$(MOD_DIR)/ -o $(OBJ_DIR)/$@.o $(MOPT)$(MOD_DIR)/
mod_OPTIMIZE: OPTIMIZE libupdate


FUNCTIONS: 
	$(FC) -c $(FFLAGS) $@.f90 -I$(MOD_DIR)/ -o $(OBJ_DIR)/$@.o $(MOPT)$(MOD_DIR)/
mod_FUNCTIONS: FUNCTIONS libupdate



INTEGRATE: 
	$(FC) -c $(FFLAGS) $@.f90 -I$(MOD_DIR)/ -o $(OBJ_DIR)/$@.o $(MOPT)$(MOD_DIR)/
mod_INTEGRATE: INTEGRATE libupdate



LIST: 
	$(FC) -c $(FFLAGS) LIST_D_ORDERED.f90 -I$(MOD_DIR)/ -o $(OBJ_DIR)/LIST_D_ORDERED.o $(MOPT)$(MOD_DIR)/	
	$(FC) -c $(FFLAGS) LIST_D_UNORDERED.f90 -I$(MOD_DIR)/ -o $(OBJ_DIR)/LIST_D_UNORDERED.o $(MOPT)$(MOD_DIR)/	
	$(FC) -c $(FFLAGS) LIST_Z_UNORDERED.f90 -I$(MOD_DIR)/ -o $(OBJ_DIR)/LIST_Z_UNORDERED.o $(MOPT)$(MOD_DIR)/	
mod_LIST: LIST libupdate



MATRIX: 
	$(FC) -c $(FFLAGS) $@.f90 -I$(MOD_DIR)/ -o $(OBJ_DIR)/$@.o $(MOPT)$(MOD_DIR)/
mod_MATRIX: MATRIX libupdate



INTERPOLATE: 
	$(FC) -c $(FFLAGS) $@.f90 -I$(MOD_DIR)/ -o $(OBJ_DIR)/$@.o $(MOPT)$(MOD_DIR)/
mod_INTERPOLATE: INTERPOLATE libupdate


RANDOM: 
	$(FC) -c $(FFLAGS) $@.f90 -I$(MOD_DIR)/ -o $(OBJ_DIR)/$@.o $(MOPT)$(MOD_DIR)/
mod_RANDOM: RANDOM libupdate


STATISTICS: 
	$(FC) -c $(FFLAGS) $@.f90 -I$(MOD_DIR)/ -o $(OBJ_DIR)/$@.o $(MOPT)$(MOD_DIR)/
mod_STATISTICS: STATISTICS libupdate


SCIFOR: 
	$(FC) -c $(FFLAGS) $@.f90 -I$(MOD_DIR)/ -o $(OBJ_DIR)/$@.o $(MOPT)$(MOD_DIR)/
mod_SCIFOR:  SCIFOR libupdate
