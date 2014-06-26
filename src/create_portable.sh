#!/bin/bash
shopt -s expand_aliases

if [ "$1" = "--help" -o "$1" = "-h" ];then
   echo "usage $0 <path_to_directory>"
   exit
fi

if [ -z "$1" ];then
    echo "usage $0 <path_to_directory>"
    exit
fi

if [ ! -d "$1" ];then
    mkdir -p $1
fi

DIR_ROOT=$1
DIR_LIB=./lib
mkdir -p $DIR_LIB
FILE_LIST_TO_COPY=`ls *.f90 *.f|paste -sd' '`

rsync -avPhHz $FILE_LIST_TO_COPY $DIR_ROOT/lib/

cat <<EOF > $DIR_ROOT/edit_makefile
#specify compiler here
FC=
#specify where the source to be compiled is
DIR=.
#specify where the executable should go
DIREXE=.
#specify which is the source name without extension [EXE=main not EXE=main.f90]
EXE=

ifeq (\$(FC),ifort)
STD =  -O2 -assume nobuffered_io
DEB =  -p -O0 -g -debug -fpe0 -traceback -check all,noarg_temp_created
MOPT=-module #leave a space here!
endif

ifeq (\$(FC),gfortran)
STD=-O2
DEB=-O0 -p -g -Wall
MOPT=-J
endif

.SUFFIXES: .f90 

#Modify this line if you do not have MKL installed or accessible.
#this variable contains mathematical libraries to which link against
#standard lapack & blas are enough otherwise
ARGS= -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread -lm

#add here your modules if any otherwise leave blank
OBJS=


#SET HERE THE COMPILATION FLAG: DEB is for debugging
FLAG=\$(STD)
#FLAG=\$(DEB)

compile: \$(OBJS)
	@echo " ..................... compile \$(EXE)........................... "
	\$(FC) \$(FLAG) $DIR_LIB/*.o \$(OBJS) \$(DIR)/\$(EXE).f90 -o \$(DIREXE)/\$(EXE) \$(ARGS) -I$DIR_LIB/
	@echo " ......................... done ................................ "
	@echo ""
	@echo ""
	@echo "created" \$(DIREXE)/\$(EXE)



#compile the library modules
libmod: COMVARS PARSE_INPUT TIMER IOTOOLS FFTGF TOOLS ERROR ARRAYS  DERIVATE OPTIMIZE FUNCTIONS GREENFUNX INTEGRATE VECTORS SQUARE_LATTICE LIST MATRIX INTERPOLATE RANDOM STATISTICS 



.f90.o:	
	\$(FC) \$(FLAG) -c \$<


COMVARS: 
	\$(FC) -c \$(FLAG) -I$DIR_LIB/ \$(MOPT)$DIR_LIB/ -o $DIR_LIB/\$@.o  $DIR_LIB/\$@.f90

PARSE_INPUT: 
	\$(FC) -c \$(FLAG) -I$DIR_LIB/ \$(MOPT)$DIR_LIB/ -o $DIR_LIB/PARSE_LIST_INPUT.o $DIR_LIB/PARSE_LIST_INPUT.f90
	\$(FC) -c \$(FLAG) -I$DIR_LIB/ \$(MOPT)$DIR_LIB/ -o $DIR_LIB/\$@.o $DIR_LIB/\$@.f90 

TIMER: 
	\$(FC) -c \$(FLAG) -I$DIR_LIB/ \$(MOPT)$DIR_LIB/ -o $DIR_LIB/\$@.o $DIR_LIB/\$@.f90

IOTOOLS: 
	\$(FC) -c \$(FLAG) -I$DIR_LIB/ \$(MOPT)$DIR_LIB/ -o $DIR_LIB/IOFILE.o $DIR_LIB/IOFILE.f90
	\$(FC) -c \$(FLAG) -I$DIR_LIB/ \$(MOPT)$DIR_LIB/ -o $DIR_LIB/IOPLOT.o $DIR_LIB/IOPLOT.f90  
	\$(FC) -c \$(FLAG) -I$DIR_LIB/ \$(MOPT)$DIR_LIB/ -o $DIR_LIB/IOREAD.o $DIR_LIB/IOREAD.f90 
	\$(FC) -c \$(FLAG) -I$DIR_LIB/ \$(MOPT)$DIR_LIB/ -o $DIR_LIB/\$@.o $DIR_LIB/\$@.f90
  
FFTGF:  COMVARS INTERPOLATE
	\$(FC) -c \$(FLAG) -I$DIR_LIB/ \$(MOPT)$DIR_LIB/ -o $DIR_LIB/\$@.o $DIR_LIB/\$@.f90

TOOLS:  
	\$(FC) -c \$(FLAG) -I$DIR_LIB/ \$(MOPT)$DIR_LIB/ -o $DIR_LIB/\$@.o $DIR_LIB/\$@.f90

ERROR:  
	\$(FC) -c \$(FLAG) -I$DIR_LIB/ \$(MOPT)$DIR_LIB/ -o $DIR_LIB/\$@.o $DIR_LIB/\$@.f90


ARRAYS:  
	\$(FC) -c \$(FLAG) -I$DIR_LIB/ \$(MOPT)$DIR_LIB/ -o $DIR_LIB/\$@.o $DIR_LIB/\$@.f90

DERIVATE: 
	\$(FC) -c \$(FLAG) -I$DIR_LIB/ \$(MOPT)$DIR_LIB/ -o $DIR_LIB/\$@.o $DIR_LIB/\$@.f90


OPTIMIZE_ROOT_FINDING: 
	\$(FC) -c \$(FLAG) -I$DIR_LIB/ \$(MOPT)$DIR_LIB/ -o $DIR_LIB/optimize_broyden_routines.o $DIR_LIB/optimize_broyden_routines.f90    
	\$(FC) -c \$(FLAG) -I$DIR_LIB/ \$(MOPT)$DIR_LIB/ -o $DIR_LIB/minpack.o $DIR_LIB/minpack.f90       
	\$(FC) -c \$(FLAG) -I$DIR_LIB/ \$(MOPT)$DIR_LIB/ -o $DIR_LIB/\$@.o $DIR_LIB/\$@.f90
OPTIMIZE_MINIMIZE: 
	\$(FC) -c \$(FLAG) -I$DIR_LIB/ \$(MOPT)$DIR_LIB/ -o $DIR_LIB/minimize.o $DIR_LIB/minimize.f 
	\$(FC) -c \$(FLAG) -I$DIR_LIB/ \$(MOPT)$DIR_LIB/ -o $DIR_LIB/CGFAM.o  $DIR_LIB/CGFAM.f
	\$(FC) -c \$(FLAG) -I$DIR_LIB/ \$(MOPT)$DIR_LIB/ -o $DIR_LIB/CGSEARCH.o  $DIR_LIB/CGSEARCH.f 
	\$(FC) -c \$(FLAG) -I$DIR_LIB/ \$(MOPT)$DIR_LIB/ -o $DIR_LIB/CGBLAS.o   $DIR_LIB/CGBLAS.f
	\$(FC) -c \$(FLAG) -I$DIR_LIB/ \$(MOPT)$DIR_LIB/ -o $DIR_LIB/optimize_cgfit_routines.o $DIR_LIB/optimize_cgfit_routines.f90
	\$(FC) -c \$(FLAG) -I$DIR_LIB/ \$(MOPT)$DIR_LIB/ -o $DIR_LIB/\$@.o $DIR_LIB/\$@.f90
OPTIMIZE: OPTIMIZE_ROOT_FINDING OPTIMIZE_MINIMIZE 
	\$(FC) -c \$(FLAG) -I$DIR_LIB/ \$(MOPT)$DIR_LIB/ -o $DIR_LIB/\$@.o $DIR_LIB/\$@.f90
 
FUNCTIONS:  
	\$(FC) -c \$(FLAG) -I$DIR_LIB/ \$(MOPT)$DIR_LIB/ -o $DIR_LIB/\$@.o $DIR_LIB/\$@.f90

GREENFUNX: 
	\$(FC) -c \$(FLAG) -I$DIR_LIB/ \$(MOPT)$DIR_LIB/ -o $DIR_LIB/\$@.o $DIR_LIB/\$@.f90

INTEGRATE: 
	\$(FC) -c \$(FLAG) -I$DIR_LIB/ \$(MOPT)$DIR_LIB/ -o $DIR_LIB/integrate_quadpack.o  $DIR_LIB/integrate_quadpack.f90  
	\$(FC) -c \$(FLAG) -I$DIR_LIB/ \$(MOPT)$DIR_LIB/ -o $DIR_LIB/\$@.o $DIR_LIB/\$@.f90


VECTORS: 
	\$(FC) -c \$(FLAG) -I$DIR_LIB/ \$(MOPT)$DIR_LIB/ -o $DIR_LIB/\$@.o $DIR_LIB/\$@.f90

SQUARE_LATTICE: 
	\$(FC) -c \$(FLAG) -I$DIR_LIB/ \$(MOPT)$DIR_LIB/ -o $DIR_LIB/\$@.o $DIR_LIB/\$@.f90

LIST: 
	\$(FC) -c \$(FLAG) -I$DIR_LIB/ \$(MOPT)$DIR_LIB/ -o $DIR_LIB/LIST_D_ORDERED.o   $DIR_LIB/LIST_D_ORDERED.f90      
	\$(FC) -c \$(FLAG) -I$DIR_LIB/ \$(MOPT)$DIR_LIB/ -o $DIR_LIB/LIST_D_UNORDERED.o $DIR_LIB/LIST_D_UNORDERED.f90 
	\$(FC) -c \$(FLAG) -I$DIR_LIB/ \$(MOPT)$DIR_LIB/ -o $DIR_LIB/LIST_Z_UNORDERED.o $DIR_LIB/LIST_Z_UNORDERED.f90    

MATRIX: 
	\$(FC) -c \$(FLAG) -I$DIR_LIB/ \$(MOPT)$DIR_LIB/ -o $DIR_LIB/\$@.o $DIR_LIB/\$@.f90

INTERPOLATE: 
	\$(FC) -c \$(FLAG) -I$DIR_LIB/ \$(MOPT)$DIR_LIB/ -o $DIR_LIB/\$@.o $DIR_LIB/\$@.f90

RANDOM: 
	\$(FC) -c \$(FLAG) -I$DIR_LIB/ \$(MOPT)$DIR_LIB/ -o $DIR_LIB/\$@.o $DIR_LIB/\$@.f90

STATISTICS: 
	\$(FC) -c \$(FLAG) -I$DIR_LIB/ \$(MOPT)$DIR_LIB/ -o $DIR_LIB/\$@.o $DIR_LIB/\$@.f90


clean: 
	@echo "Cleaning:"
	@rm -vf *.mod *.o *~ $DIR_LIB/*.o $DIR_LIB/*.mod

EOF
