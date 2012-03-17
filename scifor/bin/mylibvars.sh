#!/bin/bash

#SETUP SCIFOR LIBRARY:
if [ -z "${SFDIR}" ] ; then	# se $SFDIR NON e' definita:
    if [ -z "$1" ]; then 	# check if $1 is given: if not print error
	echo "Could not load SciFor library !"
	return 1
    else
	export SFDIR=$1
    fi
fi

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

source $SFDIR/etc/library.conf


#IF MKL is not available then set standard LAPACK/BLAS as default MATH library:
#Lapack/Blas are compiled with the chosen compiler at the installation.
if [ ! -z "${MKLDIR}" ];then
    source $MKLVARSH
else
    export LD_LIBRARY_PATH=${SFDIR}/local/lapack:${SFDIR}/local/blas:$LD_LIBRARY_PATH
    export LIBRARY_PATH=${SFDIR}/local/lapack:${SFDIR}/local/blas:$LIBRARY_PATH
    export MANPATH=${SFDIR}/local/lapack/share/man:${SFDIR}/local/blas/share/man:$MANPATH
fi


#ADD FFTW3 if defined
if [ ! -z "${FFTW3DIR}" ];then
    export INCLUDE=${FFTW3DIR}/include:${INCLUDE}
    export LD_LIBRARY_PATH=${FFTW3DIR}/lib:$LD_LIBRARY_PATH
    export LIBRARY_PATH=${FFTW3DIR}/lib:$LIBRARY_PATH
    export MANPATH=$FFTW3DIR/share/man:$MANPATH
fi

export SFLIB=$SFDIR/lib
export SFINCLUDE=$SFDIR/include
export PATH=$SFDIR/bin:$PATH

source $SFDIR/bin/bash_tools.sh


##################################################################
#REPO:
# #ADD DISLIN if defined
# if [ ! -z "${DISLIN}" ];then
#     export PATH=$DISLIN/bin:$PATH
#     export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$DISLIN
#     export LIBRARY_PATH=$LIBRARY_PATH:$DISLIN
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
