ifeq ($(FC),ifort)
OPT =  -O3 -ftz -assume nobuffered_io -openmp 
STD =  -O2 -assume nobuffered_io
DEB =  -p -O0 -g -debug -fpe0 -traceback -check all
FFLAG += $(STD) -static-intel
DFLAG += $(DEB) -static-intel
MOPT=-module
endif
ifeq ($(FC),gfortran)
OPT=-O3 -funroll-all-loops -fno-f2c
STD=-O2
DEB=-O0 -p -g -Wall
FFLAG += $(STD) -static
DFLAG += $(DEB) -static
MOPT=-J
endif

