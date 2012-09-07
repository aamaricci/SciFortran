OPT =  -O3 -ftz -assume nobuffered_io -openmp #-parallel
STD =  -O2 -assume nobuffered_io
DEB =  -p -O0 -g -debug -fpe0 -traceback  #-static-intel -check all
FFLAG += $(STD) -static-intel
DFLAG += $(DEB) -static-intel
MOPT=-module
