include ../make.inc

REV = $(shell git rev-parse HEAD)

DIR = ./scifor_portable

all: version mkdir allmod libscifor.a rmdir


mkdir:
	@mkdir -p $(OBJ_DIR)
	@mkdir -p $(MOD_DIR)

rmdir:
	@rm -rfv $(OBJ_DIR)
	@rm -rfv $(MOD_DIR)

allmod: SCIFOR_VERSION SF_COLORS SF_PAULI SF_FONTS SF_CONSTANTS SF_TIMER SF_IOTOOLS SF_FFT SF_MISC SF_ARRAYS SF_DERIVATE SF_OPTIMIZE SF_SPECIAL SF_INTEGRATE SF_LINALG SF_INTERPOLATE SF_RANDOM SF_STAT SCIFOR


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

cmdexit:
	exit 1


SCIFOR_VERSION: version
	$(FC) -c $(FFLAGS) $@.f90 -I$(MOD_DIR)/ -o $(OBJ_DIR)/$@.o $(MOPT)$(MOD_DIR)/

SF_COLORS: 
	$(FC) -c $(FFLAGS) $@.f90 -I$(MOD_DIR)/ -o $(OBJ_DIR)/$@.o $(MOPT)$(MOD_DIR)/

SF_PAULI: 
	$(FC) -c $(FFLAGS) $@.f90 -I$(MOD_DIR)/ -o $(OBJ_DIR)/$@.o $(MOPT)$(MOD_DIR)/

SF_FONTS: 
	$(FC) -c $(FFLAGS) $@.f90 -I$(MOD_DIR)/ -o $(OBJ_DIR)/$@.o $(MOPT)$(MOD_DIR)/

SF_CONSTANTS: 
	$(FC) -c $(FFLAGS) $@.f90 -I$(MOD_DIR)/ -o $(OBJ_DIR)/$@.o $(MOPT)$(MOD_DIR)/



SF_TIMER: 
	$(FC) -c $(FFLAGS) $@.f90 -I$(MOD_DIR)/ -o $(OBJ_DIR)/$@.o $(MOPT)$(MOD_DIR)/

SF_IOTOOLS: 
	$(FC) -c $(FFLAGS) IOFILE.f90 -I$(MOD_DIR)/ -o $(OBJ_DIR)/IOFILE.o $(MOPT)$(MOD_DIR)/	
	$(FC) -c $(FFLAGS) IOPLOT.f90 -I$(MOD_DIR)/ -o $(OBJ_DIR)/IOPLOT.o $(MOPT)$(MOD_DIR)/
	$(FC) -c $(FFLAGS) IOREAD.f90 -I$(MOD_DIR)/ -o $(OBJ_DIR)/IOREAD.o $(MOPT)$(MOD_DIR)/
	$(FC) -c $(FFLAGS) $@.f90 -I$(MOD_DIR)/ -o $(OBJ_DIR)/$@.o $(MOPT)$(MOD_DIR)/

SF_FFT: 
	$(FC) -c $(FFLAGS) $@.f90 -I$(MOD_DIR)/ -o $(OBJ_DIR)/$@.o $(MOPT)$(MOD_DIR)/

SF_MISC:  
	$(FC) -c $(FFLAGS) $@.f90 -I$(MOD_DIR)/ -o $(OBJ_DIR)/$@.o $(MOPT)$(MOD_DIR)/

SF_ARRAYS:  
	$(FC) -c $(FFLAGS) $@.f90 -I$(MOD_DIR)/ -o $(OBJ_DIR)/$@.o $(MOPT)$(MOD_DIR)/

SF_DERIVATE: 
	$(FC) -c $(FFLAGS) $@.f90 -I$(MOD_DIR)/ -o $(OBJ_DIR)/$@.o $(MOPT)$(MOD_DIR)/

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
SF_OPTIMIZE: OPTIMIZE_ROOT_FINDING OPTIMIZE_MINIMIZE 
	$(FC) -c $(FFLAGS) $@.f90 -I$(MOD_DIR)/ -o $(OBJ_DIR)/$@.o $(MOPT)$(MOD_DIR)/

SF_SPECIAL: 
	$(FC) -c $(FFLAGS) $@.f90 -I$(MOD_DIR)/ -o $(OBJ_DIR)/$@.o $(MOPT)$(MOD_DIR)/

SF_INTEGRATE: 
	$(FC) -c $(FFLAGS) $@.f90 -I$(MOD_DIR)/ -o $(OBJ_DIR)/$@.o $(MOPT)$(MOD_DIR)/

SF_LINALG: 
	$(FC) -c $(FFLAGS) $@.f90 -I$(MOD_DIR)/ -o $(OBJ_DIR)/$@.o $(MOPT)$(MOD_DIR)/

SF_INTERPOLATE: 
	$(FC) -c $(FFLAGS) $@.f90 -I$(MOD_DIR)/ -o $(OBJ_DIR)/$@.o $(MOPT)$(MOD_DIR)/

SF_RANDOM: 
	$(FC) -c $(FFLAGS) $@.f90 -I$(MOD_DIR)/ -o $(OBJ_DIR)/$@.o $(MOPT)$(MOD_DIR)/

SF_STAT: 
	$(FC) -c $(FFLAGS) $@.f90 -I$(MOD_DIR)/ -o $(OBJ_DIR)/$@.o $(MOPT)$(MOD_DIR)/

SCIFOR: 
	$(FC) -c $(FFLAGS) $@.f90 -I$(MOD_DIR)/ -o $(OBJ_DIR)/$@.o $(MOPT)$(MOD_DIR)/



# LIST: 
# 	$(FC) -c $(FFLAGS) LIST_D_ORDERED.f90 -I$(MOD_DIR)/ -o $(OBJ_DIR)/LIST_D_ORDERED.o $(MOPT)$(MOD_DIR)/	
# 	$(FC) -c $(FFLAGS) LIST_D_UNORDERED.f90 -I$(MOD_DIR)/ -o $(OBJ_DIR)/LIST_D_UNORDERED.o $(MOPT)$(MOD_DIR)/	
# 	$(FC) -c $(FFLAGS) LIST_Z_UNORDERED.f90 -I$(MOD_DIR)/ -o $(OBJ_DIR)/LIST_Z_UNORDERED.o $(MOPT)$(MOD_DIR)/	

