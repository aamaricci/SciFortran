FC=mpif90
FFLAG = -O2 -ffree-line-length-none
FPPMPI=-cpp -D_MPI
FFLAGS= ${FFLAG} ${FPPMPI}

.SUFFIXES: .f90
.PHONY: clean

GLOB_INC:=$(shell pkg-config --cflags scifor)

all: init integrate parse_input

init:
	make -C ./SF_INIT

integrate:
	make -C ./SF_INTEGRATE

parse_input:
	make -C ./SF_PARSE_INPUT

.f90.o:	
	$(FC) $(FFLAGS) -c $< ${GLOB_INC}

clean:
	rm -rf *.mod *.o 
