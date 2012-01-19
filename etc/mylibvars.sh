#!/bin/bash

source $HOME/lib/etc/library.conf

#ADD MY LIBRARY TO SYSTEM:
if [ -z "${LIBDIR}" ] ; then	# se $LIBDIR NON e' definita:
    if [ -z "$1" ]; then 	# check if $1 is given: if not print error
	echo "Can not load local library!"
	exit
    else
	LIBDIR=$1
    fi
fi

if [ -z "${LD_LIBRARY_PATH}" ];then
    export LD_LIBRARY_PATH="${LIBDIR}/lib"
else
    export LD_LIBRARY_PATH="${LIBDIR}/lib:$LD_LIBRARY_PATH"
fi

if [ -z "${MANPATH}" ];then
    export MANPATH="${LIBDIR}/man:$(manpath)"
else
    export MANPATH="${LIBDIR}/man:$MANPATH"
fi

if [ -z "${LIBRARY_PATH}" ];then
    export LIBRARY_PATH="${LIBDIR}/lib"
else
    export LIBRARY_PATH="${LIBDIR}/lib:$LIBRARY_PATH"
fi


#ADD MPI
if [ ! -z "${MPIDIR}" ];then
    export PATH=$MPIDIR/bin:$PATH
    export MANPATH=$MPIDIR/share/man:$MANPATH
fi


#ADD FortranGSL (local ver. unless stated otherwise):
if [ -z "${FGSLDIR}" ];then
    export FGSLDIR=${LIBOPT}/fgsl
fi
export INCLUDE="${FGSLDIR}/include:${INCLUDE}"
export LD_LIBRARY_PATH="${FGSLDIR}/lib:$LD_LIBRARY_PATH"
export LIBRARY_PATH="${FGSLDIR}/lib:$LIBRARY_PATH"


#ADD GSL (local ver. unless stated otherwise)
if [ -z "$GSLDIR}" ];then
    export GSLDIR=${LIBOPT}/gsl
fi
export INCLUDE="${GSLDIR}/include:${INCLUDE}"
export LD_LIBRARY_PATH="${GSLDIR}/lib:$LD_LIBRARY_PATH"
export LIBRARY_PATH="${GSLDIR}/lib:$LIBRARY_PATH"
export MANPATH=$GSLDIR/share/man:$MANPATH



#ADD FFTW3 (if needed use local ver.):
if [ -z "${FFTW3DIR}" ];then
    export FFTW3DIR=${LIBOPT}/fftw3
fi
export INCLUDE="${FFTW3DIR}/include:${INCLUDE}"
export LD_LIBRARY_PATH="${FFTW3DIR}/lib:$LD_LIBRARY_PATH"
export LIBRARY_PATH="${FFTW3DIR}/lib:$LIBRARY_PATH"
export MANPATH=$FFTW3DIR/share/man:$MANPATH


if [ -z "${MKLDIR}" ];then
    export MANPATH=${LIBOPT}/lapack/share/man:${LIBOPT}/blas/share/man:$MANPATH
fi


#ADD DISLIN 
if [ ! -z "${DISLIN}" ];then
    export PATH=$DISLIN/bin:$PATH
    export LD_LIBRARY_PATH="$LD_LIBRARY_PATH:$DISLIN"
    export LIBRARY_PATH="$LIBRARY_PATH:$DISLIN"
fi


#ADD CUDA (C.Weber)
if [ -z "${CUDADIR}" ]; then
    export CUDADIR=$LIBOPT/cuda
    export PATH=$CUDADIR/bin:$PATH
    export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$CUDADIR/lib64
    export LIBRARY_PATH=$LIBRARY_PATH:$CUDADIR/lib64
fi



#FILE SIZE STORE LIMIT IN Kb (IOTOOLS)
export STORE_SIZE=2048
set_store_size()
{
    SIZE=$1
    export STORE_SIZE=$SIZE
    echo $STORE_SIZE
}


export MYLIB=$LIBDIR/lib
export MYINCLUDE=$LIBDIR/include
export PATH=$LIBDIR/bin:$PATH
