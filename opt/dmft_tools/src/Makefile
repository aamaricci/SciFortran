include ../make.inc

REV = $(shell git rev-parse HEAD)

DIR = ./dmft_tools_portable

.SUFFIXES:.f90 .f

all: version allmod libdmftt.a


allmod: DMFT_TOOLS_VERSION.o PARSE_LIST_INPUT.o DMFT_PARSE_INPUT.o DMFT_FFTGF.o DMFT_FFTAUX.o DMFT_MISC.o DMFT_CONVERGENCE.o DMFT_VECTORS.o DMFT_TIGHT_BINDING.o DMFT_TOOLS.o


libdmftt.a: 
	ar crv $(LIB_DMFTT) `ls *.o | sort | uniq` 
	ranlib $(LIB_DMFTT)
	@echo 'LIBRARY IS DONE.........'
	@echo ' '

version:
	@echo "Actual version is: $(REV)"
	@echo 'character(len=41),parameter,public :: sf_version = "$(REV)"' > dmft_tools_version.inc

clean:
	@echo ""
	@echo "Cleaning ..."
	rm -fv *.o
	@echo ""
	@echo "done."
	@echo ""

portable:
	./create_portable.sh $(DIR)


cmdexit:
	exit 1



.f90.o:
	$(FC) -c $(FFLAGS) -I$(INC_SCIFOR)/ -I$(INC_TARGET)/ $(MOPT)$(INC_TARGET)/ $<

.f.o:
	$(FC) -c $(FFLAGS) $< 


# DMFT_TOOLS_VERSION: version
# 	$(FC) -c $(FFLAGS) $@.f90 -I$(MOD_DIR)/ -I$(MOD_SCIFOR)/ -o $(OBJ_DIR)/$@.o $(MOPT)$(MOD_DIR)/

# DMFT_PARSE_INPUT: 
# 	$(FC) -c $(FFLAGS)  PARSE_LIST_INPUT.f90 -o $(OBJ_DIR)/PARSE_LIST_INPUT.o $(MOPT)$(MOD_DIR)/
# 	$(FC) -c $(FFLAGS) $@.f90 -I$(MOD_DIR)/ -I$(MOD_SCIFOR)/ -o $(OBJ_DIR)/$@.o $(MOPT)$(MOD_DIR)/

# DMFT_FFTGF:  
# 	$(FC) -c $(FFLAGS) $@.f90 -I$(MOD_DIR)/ -I$(MOD_SCIFOR)/ -o $(OBJ_DIR)/$@.o $(MOPT)$(MOD_DIR)/

# DMFT_FFTAUX:  
# 	$(FC) -c $(FFLAGS) $@.f90 -I$(MOD_DIR)/ -I$(MOD_SCIFOR)/ -o $(OBJ_DIR)/$@.o $(MOPT)$(MOD_DIR)/

# DMFT_MISC:  
# 	$(FC) -c $(FFLAGS) $@.f90 -I$(MOD_DIR)/ -I$(MOD_SCIFOR)/ -o $(OBJ_DIR)/$@.o $(MOPT)$(MOD_DIR)/

# DMFT_CONVERGENCE:  
# 	$(FC) -c $(FFLAGS) $@.f90 -I$(MOD_DIR)/ -I$(MOD_SCIFOR)/ -o $(OBJ_DIR)/$@.o $(MOPT)$(MOD_DIR)/

# DMFT_VECTORS: 
# 	$(FC) -c $(FFLAGS) $@.f90 -I$(MOD_DIR)/ -I$(MOD_SCIFOR)/ -o $(OBJ_DIR)/$@.o $(MOPT)$(MOD_DIR)/

# DMFT_TIGHT_BINDING: 
# 	$(FC) -c $(FFLAGS) $@.f90 -I$(MOD_DIR)/ -I$(MOD_SCIFOR)/ -o $(OBJ_DIR)/$@.o $(MOPT)$(MOD_DIR)/

# DMFT_TOOLS: 
# 	$(FC) -c $(FFLAGS) $@.f90 -I$(MOD_DIR)/ -I$(MOD_SCIFOR)/ -o $(OBJ_DIR)/$@.o $(MOPT)$(MOD_DIR)/

# mod_DMFT_TOOLS:  DMFT_TOOLS libupdate
