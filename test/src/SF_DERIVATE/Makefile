FC=mpif90
FFLAG = -O2 -ffree-line-length-none
FPPMPI=-cpp -D_MPI
FFLAGS= ${FFLAG} ${FPPMPI}

TSTBIN=../../bin/SF_DERIVATE

GLOB_INC:=$(shell pkg-config --cflags scifor)
GLOB_LIB:=$(shell pkg-config --libs scifor | sed  "s/;/ /g"  | sed 's/\\/  /g' )


OBJS=../ASSERTING.o
SF_INIT_SOURCES=$(wildcard *.f90)

SF_INIT_EXES=$(SF_INIT_SOURCES:.f90=)

all: $(SF_INIT_EXES)

$(SF_INIT_EXES):$(SF_INIT_SOURCES)
	$(FC) $(FFLAGS) $(OBJS) $@.f90  -o $(TSTBIN)/$@.x -I../ $(GLOB_INC) $(GLOB_LIB)

clean:
	rm -rf *.mod *.o 
