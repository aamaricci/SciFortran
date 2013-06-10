#!/bin/bash

#SETUP SCIFOR LIBRARY:
if [ -z "${SFDIR}" ] ; then	# se $SFDIR NON e' definita:
    if [ -z "$1" ]; then 	# check if $1 is given: if not print error
	echo "Could not load SciFor library !"
	echo "Usage:"
	echo "$0 <path_to_SciFor_home>"
	return 1
    else
	export SFDIR=$1
    fi
fi
SFROOT=$SFDIR

export SFLIB=$SFDIR/lib
export SFMODS=$SFDIR/include
export SFETC=$SFDIR/etc
export SFBIN=$SFDIR/bin

#ADD SCIFOR TO THE SYSTEM ENVIRONMENT
if [ -z "${LD_LIBRARY_PATH}" ];then
    export LD_LIBRARY_PATH=${SFDIR}/lib
else
    export LD_LIBRARY_PATH=${SFDIR}/lib:$LD_LIBRARY_PATH
fi
if [ -z "${MANPATH}" ];then
    export MANPATH=${SFDIR}/man
else
    export MANPATH=${SFDIR}/man:$MANPATH
fi
if [ -z "${LIBRARY_PATH}" ];then
    export LIBRARY_PATH=${SFDIR}/lib
else
    export LIBRARY_PATH=${SFDIR}/lib:$LIBRARY_PATH
fi

source $SFETC/library.conf

#IF MKL is not available then set standard LAPACK/BLAS as default MATH library:
#Lapack/Blas are compiled with the chosen compiler at the installation.
if [ -z "$MKLROOT" ];then		# standard mkl variable, if defined you are using MKL
    if [ ! -z "$sf_lapack_dir" ];then
	export LD_LIBRARY_PATH=${sf_lapack_dir}:$LD_LIBRARY_PATH
	export LIBRARY_PATH=${sf_lapack_dir}:$LIBRARY_PATH
	export MANPATH=${sf_lapack_dir}/man:$MANPATH
	export LAPACK_LIB=$sf_blas_dir
    else
	echo "Not using MKL and can not find LAPACK"
    fi
    if [ ! -z "$sf_blas_dir" ];then
	export LD_LIBRARY_PATH=${sf_blas_dir}:$LD_LIBRARY_PATH
	export LIBRARY_PATH=${sf_blas_dir}:$LIBRARY_PATH
	export MANPATH=${sf_blas_dir}/man:$MANPATH
	export BLAS_LIB=$sf_blas_dir
    else
	echo "Not using MKL and can not find BLAS"
    fi
fi


#ADD FFTW_3 to ENV
if [ ! -z "$sf_fftw_dir" ];then
    export INCLUDE=$sf_fftw_dir/include:${INCLUDE}
    export LD_LIBRARY_PATH=$sf_fftw_dir/lib:$LD_LIBRARY_PATH
    export LIBRARY_PATH=$sf_fftw_dir/lib:$LIBRARY_PATH
    export FFTW_LIB=$sf_fftw_dir
fi

#ADD FFTPACK to ENV
if [ ! -z "$sf_fftpack_dir" ];then
    export INCLUDE=$sf_fftpack_dir/include:${INCLUDE}
    export LD_LIBRARY_PATH=$sf_fftpack_dir/lib:$LD_LIBRARY_PATH
    export LIBRARY_PATH=$sf_fftpack_dir/lib:$LIBRARY_PATH
    export FFTPACK_LIB=$sf_fftpack_dir
fi


#ADD ARPACK to ENV
if [ ! -z "$sf_arpack_dir" ];then
    export INCLUDE=$sf_arpack_dir/include:${INCLUDE}
    export LD_LIBRARY_PATH=$sf_arpack_dir/lib:$LD_LIBRARY_PATH
    export LIBRARY_PATH=$sf_arpack_dir/lib:$LIBRARY_PATH
    export ARPACK_LIB=$sf_arpack_dir
fi


#ADD MINPACK to ENV
if [ ! -z "$sf_minpack_dir" ];then
    export INCLUDE=$sf_minpack_dir/include:${INCLUDE}
    export LD_LIBRARY_PATH=$sf_minpack_dir/lib:$LD_LIBRARY_PATH
    export LIBRARY_PATH=$sf_minpack_dir/lib:$LIBRARY_PATH
    export MINPACK_LIB=$sf_minpack_dir
fi



##################################################################
#REPO: THIS SOFTWARE SHOULD BE MADE AVAILABLE **EXTERNALLY* W/ SCIFOR LIB.
#SCIFOR DOES NOT RELY ON THESE LIBS, THOUGH SOME ROUTINES MAY USE PART OF THEM. 
#E.G. DLPLOT DEPENDS ON DISLIN
# #ADD DISLIN if defined
# if [ ! -z "${DISLIN}" ];then
#     export PATH=$DISLIN/bin:$PATH
#     export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$DISLIN
#     export LIBRARY_PATH=$LIBRARY_PATH:$DISLIN
# fi

# if [ ! -z "${SFMINPACK}" ];then
#     export LD_LIBRARY_PATH=${SFMINPACK}:$LD_LIBRARY_PATH
#     export LIBRARY_PATH=${SFMINPACK}:$LIBRARY_PATH
# fi

# #ADD CUDA (C.Weber) if defined
# if [ ! -z "${CUDADIR}" ]; then
#     export PATH=$CUDADIR/bin:$PATH
#     export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$CUDADIR/lib64
#     export LIBRARY_PATH=$LIBRARY_PATH:$CUDADIR/lib64
# fi

# #ADD GSL if defined
# if [ ! -z "$GSLDIR}" ];then
#     export INCLUDE=${GSLDIR}/include:${INCLUDE}
#     export LD_LIBRARY_PATH=${GSLDIR}/lib:$LD_LIBRARY_PATH
#     export LIBRARY_PATH=${GSLDIR}/lib:$LIBRARY_PATH
#     export MANPATH=$GSLDIR/share/man:$MANPATH
# fi

# #ADD FortranGSL if definted
# if [ ! -z "${FGSLDIR}" ];then
#     export INCLUDE=${FGSLDIR}/include:${INCLUDE}
#     export LD_LIBRARY_PATH=${FGSLDIR}/lib:$LD_LIBRARY_PATH
#     export LIBRARY_PATH=${FGSLDIR}/lib:$LIBRARY_PATH
# fi
