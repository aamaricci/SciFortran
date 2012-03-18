#ifort options
#OPT =  -O3 -ftz -assume nobuffered_io -openmp #-parallel
#STD =  -O2 -assume nobuffered_io
#DEB =  -p -traceback -O0 -g -debug -fpe0 -traceback  #-static-intel -check all


#gfortran options
OPT = -pg -O3 #find and add other specific to fortran compiler in use
STD = -pg -O1  
DEB = -pg -O0 -g3 -fbounds-check -fbacktrace -Wall -Wextra -Wconversion -pedantic

#########################################################################

FFLAG += $(STD) -static
DFLAG += $(DEB) -static
