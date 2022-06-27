FC=mpif90
FFLAG = -O2 -ffree-line-length-none
DFLAG = -O2 -p -g -fimplicit-none -Wsurprising  -Waliasing -fwhole-file -fcheck=all -pedantic -fbacktrace -ffree-line-length-none
OFLAG = -O3 -ffast-math -march=native -funroll-loops -ffree-line-length-none
FPPSERIAL =-cpp -D_
FPPMPI=-cpp -D_MPI
FFLAGS= ${FFLAG} ${FPPMPI}

TSTBIN=../bin

GLOB_INC:=$(shell pkg-config --cflags scifor)
GLOB_LIB:=$(shell pkg-config --libs scifor | sed  "s/;/ /g"  | sed 's/\\/  /g' )

##$ Extends the implicit support of the Makefile to .f90 files
.SUFFIXES: .f90

OBJS=asserting.o

all: $(OBJS)
	$(FC) $(FFLAGS) $(OBJS) test_sf_eigh.f90 -o $(TSTBIN)/test_eigh.x $(GLOB_INC) $(GLOB_LIB)
	$(FC) $(FFLAGS) $(OBJS) test_sf_trapz.f90 -o $(TSTBIN)/test_trapz.x $(GLOB_INC) $(GLOB_LIB)
	$(FC) $(FFLAGS) $(OBJS) test_sf_parsing.f90 -o $(TSTBIN)/test_parsing.x $(GLOB_INC) $(GLOB_LIB)


.f90.o:	
	$(FC) $(FFLAGS) -c $< ${GLOB_INC}

clean:
	rm -rf *.mod *.o
