include ../make.inc

REV = $(shell git rev-parse HEAD)

DIR = ./scifor_portable

.SUFFIXES:.f90 .f

.PHONY: arpack blas fftpack lapack minpack parpack quadpack sf_init sf_mpi sf_iotools sf_fft  sf_derivate sf_optimize sf_special sf_integrate sf_interpolate sf_random sf_linalg sf_sp_linalg  sf_parse_input scifor 


serial:   version  arpack blas fftpack lapack minpack quadpack scifor lib
parallel: version parpack blas fftpack lapack minpack quadpack scifor lib
update: version scifor lib

arpack:
	@make -C ./arpack serial

blas:
	@make -C ./blas

fftpack:
	@make -C ./fftpack

lapack:
	@make -C ./lapack

minpack:
	@make -C ./minpack

parpack:
	@make -C ./arpack parallel

quadpack:
	@make -C ./quadpack

sf_init:
	@make -C ./SF_INIT

sf_mpi:
	@make -C ./SF_MPI

sf_iotools:
	@make -C ./SF_IOTOOLS

sf_derivate:
	@make -C ./SF_DERIVATE

sf_linalg:
	@make -C ./SF_LINALG

sf_sp_linalg:
	@make -C ./SF_SP_LINALG

sf_optimize:
	@make -C ./SF_OPTIMIZE

sf_integrate:
	@make -C ./SF_INTEGRATE

sf_interpolate:
	@make -C ./SF_INTERPOLATE

sf_special:
	@make -C ./SF_SPECIAL

sf_fft:
	@make -C ./SF_FFT

sf_random:
	@make -C ./SF_RANDOM

sf_parse_input:
	@make -C ./SF_PARSE_INPUT

scifor: SCIFOR_VERSION.o sf_init sf_mpi sf_iotools sf_fft  sf_derivate sf_optimize sf_special sf_integrate sf_interpolate sf_random sf_linalg sf_sp_linalg  sf_parse_input SCIFOR.o


lib:
	@rsync -avP *.mod $(INC_TARGET)/
	ar crv $(LIB_SCIFOR) `ls *.o | sort | uniq`
	ranlib $(LIB_SCIFOR)
	@echo 'SciFortran library is done.........'

version:
	@echo "Actual version is: $(REV)"
	@echo 'character(len=41),parameter,public :: version_scifor = "$(REV)"' > scifor_version.inc


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
	$(FC) -c $(FCOPT) -I./ $(MOPT)./ $< -o ./$@

.f.o:
	$(FC) -c $(FCOPT) $< -o ./$@
