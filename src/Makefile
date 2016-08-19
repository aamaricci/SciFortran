include ../make.inc

REV = $(shell git rev-parse HEAD)

DIR = ./scifor_portable

.SUFFIXES:.f90 .f

.PHONY: arpack blas fftpack lapack minpack parpack quadpack scifor


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

scifor: SCIFOR_VERSION.o SF_COLORS.o SF_PAULI.o SF_FONTS.o SF_CONSTANTS.o SF_TIMER.o IOFILE.o IOPLOT.o IOREAD.o SF_IOTOOLS.o SF_FFT.o SF_MISC.o SF_ARRAYS.o SF_DERIVATE.o optimize_broyden_routines.o OPTIMIZE_ROOT_FINDING.o minimize.o cgfam.o cgsearch.o cgblas.o optimize_cgfit_routines.o OPTIMIZE_ROOT_FINDING.o OPTIMIZE_MINIMIZE.o SF_OPTIMIZE.o SF_SPECIAL.o SF_INTEGRATE.o SF_LINALG.o SF_SP_LINALG.o SF_INTERPOLATE.o SF_RANDOM.o SF_STAT.o LIST_INPUT.o SF_PARSE_INPUT.o SCIFOR.o


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
	$(FC) -c $(FFLAGS) -I./ $(MOPT)./ $< -o ./$@

.f.o:
	$(FC) -c $(FFLAGS) $< -o ./$@
