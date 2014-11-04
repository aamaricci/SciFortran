include ../make.inc

REV = $(shell git rev-parse HEAD)

DIR = ./dmft_tools_portable

all: version mkdir allmod libdmftt.a


mkdir:
	@mkdir -p $(OBJ_DIR)
	@mkdir -p $(MOD_DIR)

allmod: DMFT_TOOLS_VERSION PARSE_INPUT FFTGF TOOLS ERROR VECTORS SQUARE_LATTICE DMFT_TOOLS


libdmftt.a: 
	@rsync -avP $(MOD_DIR)/*.mod $(MOD_DMFTT)/
	ar crv $(LIB_DMFTT) `ls $(OBJ_DIR)/*.o | sort | uniq`  
	ranlib $(LIB_DMFTT)
	@echo 'LIBRARY IS DONE.........'
	@echo ' '


version:
	@echo "Actual version is: $(REV)"
	@echo 'character(len=41),parameter,public :: sf_version = "$(REV)"' > dmft_tools_version.inc

clean: 
	@rm -fv *.o *.mod *.a

portable:
	./create_portable.sh $(DIR)

libupdate:
	@echo "	rsync -avP $(MOD_DIR)/*.mod $(MOD_DMFTT)/"
	@rsync -avP $(MOD_DIR)/*.mod $(MOD_DMFTT)/
	ar crv $(LIB_DMFTT) `ls $(OBJ_DIR)/*.o | sort | uniq`  
	ranlib $(LIB_DMFTT)
	@echo 'LIBRARY IS UPDATED.........'
	@echo ' '



DMFT_TOOLS_VERSION: version
	$(FC) -c $(FFLAGS) $@.f90 -I$(MOD_DIR)/ -o $(OBJ_DIR)/$@.o $(MOPT)$(MOD_DIR)/
mod_DMFT_TOOLS_VERSION: DMFT_TOOLS_VERSION libupdate


PARSE_INPUT: 
	$(FC) -c $(FFLAGS)  PARSE_LIST_INPUT.f90 -o $(OBJ_DIR)/PARSE_LIST_INPUT.o $(MOPT)$(MOD_DIR)/
	$(FC) -c $(FFLAGS) $@.f90 -I$(MOD_DIR)/ -o $(OBJ_DIR)/$@.o $(MOPT)$(MOD_DIR)/
mod_PARSE_INPUT: PARSE_INPUT libupdate


FFTGF:  
	$(FC) -c $(FFLAGS) $@.f90 -I$(MOD_DIR)/ -o $(OBJ_DIR)/$@.o $(MOPT)$(MOD_DIR)/
mod_FFTGF: FFTGF libupdate 


TOOLS:  
	$(FC) -c $(FFLAGS) $@.f90 -I$(MOD_DIR)/ -o $(OBJ_DIR)/$@.o $(MOPT)$(MOD_DIR)/
mod_TOOLS:  TOOLS libupdate


ERROR:  
	$(FC) -c $(FFLAGS) $@.f90 -I$(MOD_DIR)/ -o $(OBJ_DIR)/$@.o $(MOPT)$(MOD_DIR)/
mod_ERROR: ERROR libupdate


VECTORS: 
	$(FC) -c $(FFLAGS) $@.f90 -I$(MOD_DIR)/ -o $(OBJ_DIR)/$@.o $(MOPT)$(MOD_DIR)/
mod_VECTORS: VECTORS libupdate


SQUARE_LATTICE: 
	$(FC) -c $(FFLAGS) $@.f90 -I$(MOD_DIR)/ -o $(OBJ_DIR)/$@.o $(MOPT)$(MOD_DIR)/
mod_SQUARE_LATTICE: SQUARE_LATTICE libupdate


DMFT_TOOLS: 
	$(FC) -c $(FFLAGS) $@.f90 -I$(MOD_DIR)/ -o $(OBJ_DIR)/$@.o $(MOPT)$(MOD_DIR)/
mod_DMFT_TOOLS:  DMFT_TOOLS libupdate
