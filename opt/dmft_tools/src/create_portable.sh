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
STD =  -O2 
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
FLAG=\$(STD) -I$DIR_LIB/
#FLAG=\$(DEB)

compile: \$(OBJS)
	@echo " ..................... compile \$(EXE)........................... "
	\$(FC) \$(FLAG) $DIR_LIB/*.o \$(OBJS) \$(DIR)/\$(EXE).f90 -o \$(DIREXE)/\$(EXE) \$(ARGS) 
	@echo " ......................... done ................................ "
	@echo ""
	@echo ""
	@echo "created" \$(DIREXE)/\$(EXE)



#compile the library modules
lib:  PARSE_INPUT FFTGF TOOLS ERROR VECTORS SQUARE_LATTICE 

.f90.o:	
	\$(FC) \$(FLAG) -c \$<


PARSE_INPUT: 
	\$(FC) -c \$(FLAG)  \$(MOPT)$DIR_LIB/ -o $DIR_LIB/PARSE_LIST_INPUT.o $DIR_LIB/PARSE_LIST_INPUT.f90
	\$(FC) -c \$(FLAG)  \$(MOPT)$DIR_LIB/ -o $DIR_LIB/\$@.o $DIR_LIB/\$@.f90 

FFTGF: 
	\$(FC) -c \$(FLAG)  \$(MOPT)$DIR_LIB/ -o $DIR_LIB/\$@.o $DIR_LIB/\$@.f90

ERROR:  
	\$(FC) -c \$(FLAG)  \$(MOPT)$DIR_LIB/ -o $DIR_LIB/\$@.o $DIR_LIB/\$@.f90

VECTORS: 
	\$(FC) -c \$(FLAG)  \$(MOPT)$DIR_LIB/ -o $DIR_LIB/\$@.o $DIR_LIB/\$@.f90

SQUARE_LATTICE: 
	\$(FC) -c \$(FLAG)  \$(MOPT)$DIR_LIB/ -o $DIR_LIB/\$@.o $DIR_LIB/\$@.f90


clean: 
	@echo "Cleaning here:"
	@rm -vf *.mod *.o *~ 

clean_lib: 
	@echo "Cleaning lib:"
	@rm -vf $DIR_LIB/*.o $DIR_LIB/*.mod

EOF
