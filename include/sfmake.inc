include options.inc

.SUFFIXES: .f90 

#REVISION SOFTWARE GIT:
REV = $(shell git rev-parse HEAD)
VER = 'character(len=41),parameter :: revision = "$(REV)"' > revision.inc

#BASH COMPLETION SCRIPT FOR SciFor
SFCOMPLETION=${SFROOT}/bin/scifor_completion.sh

#SciFor INCLUDE modules:
SFINCLUDE= -I${SFROOT}/include
SFINCLUDE_DEB= -I${SFROOT}/include/debug



#SciFor LIBRARY compilation arguments:
SFARGS 	   =  -lscifor
SFARGS_DEB =  -lscifor_deb
#----------------------------------------------------



#TO BE REMOVED:

# ifdef MKLROOT
#  ifeq ($(SFOS),Linux)
#     MKLARGS     = -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread -lm
#  endif
#  ifeq ($(SFOS),Darwin)
#     MKLARGS	= -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread -lm
#  endif
# endif

# #SETUP MATHLIB: the main mathematical library. If MKL is available use it, if not use just std lapack and blas
# ifdef MKLROOT
# #
#  ifeq ($(SFOS),Linux)
#     MATHLIB     = -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread -lm
#  endif
#  ifeq ($(SFOS),Darwin)
#     MATHLIB	= -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread -lm
#  endif
# #
# else
# #
#  ifdef BLASROOT
#     MATHLIB     += -lblas
#  endif
#  ifdef LAPACKROOT
#     MATHLIB	+= -llapack
#  endif	
# #
# endif
# ifdef FFTPACKROOT
#    MATHLIB      +=-lfftpack
# endif
# ifdef MINPACKROOT
#    MATHLIB      +=-lminpack
# endif
# ifdef QUADPACKROOT
#    MATHLIB      +=-lquadpack
# endif


# #----------------------------------------------------
# #SciFor LIBRARY:
# SFLIBS 	   =  ${SFINCLUDE} -lscifor #${MATHLIB}
# SFLIBS_DEB =  ${SFINCLUDE_DEB} -lscifor_deb #${MATHLIB}
# #----------------------------------------------------



###################################################################
# ifdef ARPACK_LIB
# MATHLIB     +=-lparpack -larpack #-lparpack_debug -larpack_debug
# endif

# ifdef FFTW_LIB
# MATHLIB      +=-lfftw3
# SFINCLUDE    +=-I$(FFTW_LIB)/include/
# SFINCLUDE_DEB+=-I$(FFTW_LIB)/include/
# endif


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
