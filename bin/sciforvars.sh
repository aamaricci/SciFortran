#!/bin/bash

#SETUP SCIFOR LIBRARY:
if [ -z "${SFDIR}" ] ; then	# se $SFDIR NON e' definita:
    if [ -z "$1" ]; then 	# check if $1 is given: if not print error
	echo "Could not load SciFor library !" >&2
	echo "Usage:" >&2
	echo "$0 <path_to_SciFor_home>" >&2
	return 1
    else
	export SFROOT=$1
    fi
fi

export SFROOT=$SFDIR

#ADD SCIFOR TO THE SYSTEM ENVIRONMENT
if [ -z "${LD_LIBRARY_PATH}" ];then
    export LD_LIBRARY_PATH=${SFROOT}/lib
else
    export LD_LIBRARY_PATH=${SFROOT}/lib:$LD_LIBRARY_PATH
fi

if [ -z "${LIBRARY_PATH}" ];then
    export LIBRARY_PATH=${SFROOT}/lib
else
    export LIBRARY_PATH=${SFROOT}/lib:$LIBRARY_PATH
fi

if [ -z "${FPATH}" ];then
    export FPATH=${SFROOT}/include
else
    export FPATH=${SFROOT}/include:$FPATH
fi

if [ -z "${CPATH}" ];then
    export CPATH=${SFROOT}/include
else
    export CPATH=${SFROOT}/include:$CPATH
fi

source $SFROOT/etc/library.conf






#TO BE REMOVED

# #If MKL is set (by definition of MKLROOT environment variable) set related precompilation flag 
# if [ ! -z "$MKLROOT" ];then
#     export SF_PRECOMP_MKL=MKL
# fi

# #Start adding required libraries to the system. 
# #If not defined a warning message to std.err is given.
# #BLAS
# if [ ! -z "$BLASROOT" ];then
#     add_library_to_system $BLASROOT
#     export PRECOMP_BLAS=BLAS
# else
#     echo "SciFor: can not add BLAS to system." >&2
# fi
# #LAPACK
# if [ ! -z "$LAPACKROOT" ];then
#     add_library_to_system $LAPACKROOT
#     export PRECOMP_LAPACK=LAPACK
# else
#     echo "SciFor: can not add LAPACK to system." >&2
# fi
# #FFTPACK
# if [ ! -z "$FFTPACKROOT" ];then 
#     add_library_to_system $FFTPACKROOT
#     export PRECOMP_FFTPACK=FFTPACK
# else
#     echo "SciFor: can not add FFTPACK to system." >&2
# fi
# #MINPACK
# if [ ! -z "$MINPACKROOT" ];then 
#     add_library_to_system $MINPACKROOT
#     export PRECOMP_MINPACK=MINPACK
# else
#     echo "SciFor: can not add MINPACK to system." >&2
# fi
# #QUADPACK
# if [ ! -z "$QUADPACKROOT" ];then 
#     add_library_to_system $QUADPACKROOT
#     export PRECOMP_QUADPACK=QUADPACK
# else
#     echo "SciFor: can not add QUADPACK to system." >&2
# fi
