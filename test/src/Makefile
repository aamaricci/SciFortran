FC=mpif90
FFLAG = -O2 -ffree-line-length-none
FPPMPI=-cpp -D_MPI
FFLAGS= ${FFLAG} ${FPPMPI}

.SUFFIXES: .f90
.PHONY: clean

GLOB_INC:=$(shell pkg-config --cflags scifor)

OBJS=ASSERTING.o

all: $(OBJS) init integrate parse_input iotools derivate fft interpolate optimize

init:
	make -C ./SF_INIT

integrate:
	make -C ./SF_INTEGRATE

parse_input:
	make -C ./SF_PARSE_INPUT

iotools:
	make -C ./SF_IOTOOLS

derivate:
	make -C ./SF_DERIVATE

fft:
	make -C ./SF_FFT

interpolate:
	make -C ./SF_INTERPOLATE

optimize:
	make -C ./SF_OPTIMIZE

.f90.o:	
	$(FC) $(FFLAGS) -c $< ${GLOB_INC}

clean:
	rm -rf *.mod *.o 
