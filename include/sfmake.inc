include options.inc

.SUFFIXES: .f90 

#REVISION SOFTWARE GIT:
REV = $(shell git rev-parse HEAD)
VER = 'character(len=41),parameter :: revision = "$(REV)"' > revision.inc

#SciFor INCLUDE modules:
SFINCLUDE= -I${SFMODS}
SFINCLUDE_DEB= -I${SFMODS}/debug


#SETUP MATHLIB: the main mathematical library. If MKL is available use it, if not use just std lapack and blas
ifdef MKLROOT
 ifeq ($(SFOS),Linux)
    MATHLIB     = -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread -lm
 endif
 ifeq ($(SFOS),Darwin)
    MATHLIB	= -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread -lm
 endif
else
MATHLIB	    = -llapack -lblas
endif

ifdef FFTPACK_LIB
MATHLIB      +=-lfftpack
endif

ifdef MINPACK_LIB
MATHLIB      +=-lminpack
endif

ifdef ARPACK_LIB
MATHLIB     +=-lparpack -larpack #-lparpack_debug -larpack_debug
endif

ifdef FFTW_LIB
MATHLIB      +=-lfftw3
SFINCLUDE    +=-I$(FFTW_LIB)/include/
SFINCLUDE_DEB+=-I$(FFTW_LIB)/include/
endif


#----------------------------------------------------
#SciFor LIBRARY:
SFLIBS 	   =  ${SFINCLUDE} -lscifor ${MATHLIB}
SFLIBS_DEB =  ${SFINCLUDE_DEB} -lscifor_deb ${MATHLIB}
#----------------------------------------------------


###################################################################
#REPO:
# ifdef CUDADIR
# GPU     = -L$(SFLIB)/ 	-lgpu
# CUDA    = -L/opt/cuda_2.3/lib64 -lcufft -lcuda -lcudart -lcublas
# CUDA_S  = -L/opt/cuda_sdk_2.3/lib -lcutil
# GPU_LIBS= ${GPU} ${CUDA_STATIC} ${CUDA} 
# endif


# ifdef DISLIN
# DSL         =  -L$(DISLIN)/lib  -ldislin
# X11 	    =  -L/usr/lib -lX11 -lXt -lXext -lxcb -lX11-xcb -lXau -lXdmcp #-lxcb-xlib 
# DLPLOT      =  -ldlplot
# DLPLOT_DEB  =  -ldlplot_deb
# DSL_LIBS    = ${DLPLOT} ${DSL} ${X11}
# DSL_LIBS_DEB= ${DLPLOT_DEB} ${DSL} ${X11}
# DSL_MODS= -I${DISLIN}/ifc
# endif
