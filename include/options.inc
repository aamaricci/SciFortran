OPT =  -O3 -ftz  -openmp 
STD =  -O2 
DEB =  -p -O0 -g -debug -fpe0 -traceback -check all,noarg_temp_created
FFLAG = $(STD) -static-intel
DFLAG = $(DEB) -static-intel
MOPT=-module #leave a space here!

ifeq ($(FC),gfortran)
OPT=-O3 -funroll-all-loops -fno-f2c
STD=-O2
DEB=-O0 -p -g -Wall -fPIC -fmax-errors=1 -g -fcheck=all -fbacktrace
FFLAG = $(STD) -static
DFLAG = $(DEB) -static
MOPT=-J
endif

ifeq ($(FC),xlf90)
OPT =  -O3 -qarch=qp -qtune=qp
STD =  -O2
DEB =  -O0 -C
FFLAG = $(STD)
DFLAG = $(DEB)
MOPT =-qmoddir=
endif



