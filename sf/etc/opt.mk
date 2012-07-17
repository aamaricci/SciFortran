OPT =  -O3 -ftz -assume nobuffered_io -openmp #-parallel
STD =  -O2 -assume nobuffered_io
DEB =  -p -traceback -O0 -g -debug -fpe0 -traceback  #-static-intel -check all
FFLAG += $(STD) -static
DFLAG += $(DEB) -static
MOPT=-module
