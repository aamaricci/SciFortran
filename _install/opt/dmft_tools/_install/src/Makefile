include ../make.inc

REV = $(shell git rev-parse HEAD)

DIR = ./dmft_tools_portable

all: version mkdir allmod libdmftt.a rmdir


mkdir:
	@mkdir -p $(OBJ_DIR)
	@mkdir -p $(MOD_DIR)

rmdir:
	@rm -rfv $(OBJ_DIR)
	@rm -rfv $(MOD_DIR)

allmod: DMFT_TOOLS_VERSION DMFT_PARSE_INPUT DMFT_FFTGF DMFT_FFTAUX DMFT_MISC DMFT_CONVERGENCE DMFT_VECTORS DMFT_TIGHT_BINDING DMFT_TOOLS


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
	@rsync -avP --del $(MOD_DIR)/*.mod $(MOD_DMFTT)/
	ar crv $(LIB_DMFTT) `ls $(OBJ_DIR)/*.o | sort | uniq`  
	ranlib $(LIB_DMFTT)
	@echo 'LIBRARY IS UPDATED.........'
	@echo ' '

cmdexit:
	exit 1


DMFT_TOOLS_VERSION: version
	$(FC) -c $(FFLAGS) $@.f90 -I$(MOD_DIR)/ -I$(MOD_SCIFOR)/ -o $(OBJ_DIR)/$@.o $(MOPT)$(MOD_DIR)/

DMFT_PARSE_INPUT: 
	$(FC) -c $(FFLAGS)  PARSE_LIST_INPUT.f90 -o $(OBJ_DIR)/PARSE_LIST_INPUT.o $(MOPT)$(MOD_DIR)/
	$(FC) -c $(FFLAGS) $@.f90 -I$(MOD_DIR)/ -I$(MOD_SCIFOR)/ -o $(OBJ_DIR)/$@.o $(MOPT)$(MOD_DIR)/

DMFT_FFTGF:  
	$(FC) -c $(FFLAGS) $@.f90 -I$(MOD_DIR)/ -I$(MOD_SCIFOR)/ -o $(OBJ_DIR)/$@.o $(MOPT)$(MOD_DIR)/

DMFT_FFTAUX:  
	$(FC) -c $(FFLAGS) $@.f90 -I$(MOD_DIR)/ -I$(MOD_SCIFOR)/ -o $(OBJ_DIR)/$@.o $(MOPT)$(MOD_DIR)/

DMFT_MISC:  
	$(FC) -c $(FFLAGS) $@.f90 -I$(MOD_DIR)/ -I$(MOD_SCIFOR)/ -o $(OBJ_DIR)/$@.o $(MOPT)$(MOD_DIR)/

DMFT_CONVERGENCE:  
	$(FC) -c $(FFLAGS) $@.f90 -I$(MOD_DIR)/ -I$(MOD_SCIFOR)/ -o $(OBJ_DIR)/$@.o $(MOPT)$(MOD_DIR)/

DMFT_VECTORS: 
	$(FC) -c $(FFLAGS) $@.f90 -I$(MOD_DIR)/ -I$(MOD_SCIFOR)/ -o $(OBJ_DIR)/$@.o $(MOPT)$(MOD_DIR)/

DMFT_TIGHT_BINDING: 
	$(FC) -c $(FFLAGS) $@.f90 -I$(MOD_DIR)/ -I$(MOD_SCIFOR)/ -o $(OBJ_DIR)/$@.o $(MOPT)$(MOD_DIR)/

DMFT_TOOLS: 
	$(FC) -c $(FFLAGS) $@.f90 -I$(MOD_DIR)/ -I$(MOD_SCIFOR)/ -o $(OBJ_DIR)/$@.o $(MOPT)$(MOD_DIR)/

mod_DMFT_TOOLS:  DMFT_TOOLS libupdate
