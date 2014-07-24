#!/bin/bash

#SETUP SCIFOR LIBRARY:
if [ -z "${SFDIR}" ] ; then	# se $SFDIR NON e' definita:
    if [ -z "$1" ]; then 	# check if $1 is given: if not print error
	echo "Could not load SciFor library !"
	echo "Usage:"
	echo "$0 <path_to_SciFor_home>"
	return 1
    else
	export SFROOT=$1
    fi
fi

export SFROOT=$SFDIR
export SFLIB=$SFROOT/lib
export SFMODS=$SFROOT/include

add_library_to_system(){    
    LIB=$1
    if [ -z "${LD_LIBRARY_PATH}" ];then
	export LD_LIBRARY_PATH=${LIB}/lib
    else
	export LD_LIBRARY_PATH=${LIB}/lib:$LD_LIBRARY_PATH
    fi

    if [ -z "${LIBRARY_PATH}" ];then
	export LIBRARY_PATH=${LIB}/lib
    else
	export LIBRARY_PATH=${LIB}/lib:$LIBRARY_PATH
    fi

    if [ -z "${FPATH}" ];then
	export FPATH=${LIB}/include
    else
	export FPATH=${LIB}/include:$FPATH
    fi

    if [ -z "${CPATH}" ];then
	export CPATH=${LIB}/include
    else
	export CPATH=${LIB}/include:$CPATH
    fi

    if [ -z "${MANPATH}" ];then
	export MANPATH=${LIB}/man
    else
	export MANPATH=${LIB}/man:$MANPATH
    fi
}

#ADD SCIFOR TO THE SYSTEM ENVIRONMENT
add_library_to_system $SFROOT
source $SFROOT/etc/library.conf


#IF MKL is not available then set standard LAPACK/BLAS as default MATH library:
#Lapack/Blas are compiled with the chosen compiler at the installation.
if [ -z "$MKLROOT" ];then		# standard mkl variable, if defined you are using MKL
    if [ ! -z "$SF_LAPACK_DIR" ];then
	add_library_to_system $SF_LAPACK_DIR
	export LAPACK_LIB=$SF_LAPACK_DIR
    else
	echo "Not using MKL and can not find LAPACK"
    fi
    if [ ! -z "$SF_BLAS_DIR" ];then
	add_library_to_system $SF_BLAS_DIR
	export BLAS_LIB=$SF_BLAS_DIR
    else
	echo "Not using MKL and can not find BLAS"
    fi
fi


for LIB in $SF_LIST_LIB;
do
    if [  -d "$LIB" ];then
	add_library_to_system $$LIB
	lib_name=$(basename $LIB)
    else
	echo "$LIB does not exist or I can not find it: skip installation..."
    fi
done

#IF LIBRARIES ARE AVAILABLE SET PRE-COMPILATION FLAGS ACCORDINGLY
#availability of these libraries is independent of the fact that
#SciFor is actually adding them to the system, you may want to 
#do it your own for some reason.
#This is just a way to communicate SciFor that you have or don't have 
#access to these libraries so to include or not routines that are 
#based on these libraries.
if [ ! -z $SF_FFTPACK ];then 
    export PRECOMP_FFTPACK=FFTPACK
fi
if [ ! -z $SF_MINPACK ];then 
    export PRECOMP_MINPACK=MINPACK
fi
if [ ! -z $SF_QUADPACK ];then 
    export PRECOMP_QUADTPACK=QUADPACK
fi


# #ADD FFTPACK to ENV
# if [ ! -z "$sf_fftpack_dir" ];then
#     add_library_to_system $sf_fftpack_dir
#     export FFTPACK_LIB=$sf_fftpack_dir
#     export PRECOMP_FFTPACK=SC_FFTPACK
# fi

# #ADD MINPACK to ENV
# if [ ! -z "$sf_minpack_dir" ];then
#     add_library_to_system $sf_minpack_dir
#     export MINPACK_LIB=$sf_minpack_dir
#     export PRECOMP_MINPACK=SC_MINPACK
# fi

