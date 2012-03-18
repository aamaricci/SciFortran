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
export SFLIB=$SFDIR/lib
export SFINCLUDE=$SFDIR/include
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
source $SFBIN/bash_tools.sh
export PATH=$SFBIN:$PATH




#IF MKL is not available then set standard LAPACK/BLAS as default MATH library:
#Lapack/Blas are compiled with the chosen compiler at the installation.
if [ ! -z "$MKLDIR" ];then
    source $MKLVARSH
else
    export LD_LIBRARY_PATH=${SFBLAS}:${SFLAPACK}:$LD_LIBRARY_PATH
    export LIBRARY_PATH=${SFBLAS}:${SFLAPACK}:$LIBRARY_PATH
    export MANPATH=${SFBLAS}/man:${SFLAPACK}/man:$MANPATH
fi


if [ ! -z "$SFFFTW3" ];then
    export INCLUDE=$SFFFTW3/include:${INCLUDE}
    export LD_LIBRARY_PATH=$SFFFTW3/lib:$LD_LIBRARY_PATH
    export LIBRARY_PATH=$SFFFTW3/lib:$LIBRARY_PATH
    export MANPATH=$SFFFTW3/share/man:$MANPATH
fi



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
