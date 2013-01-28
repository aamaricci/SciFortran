include options.inc

#SciFor MODULES:
SFMODS 	= -I$(SFINCLUDE)
SFMODS_DEB= -I$(SFINCLUDE)/debug/

ifdef MKLROOT
 ifeq ($(OS),LINUX)
#IFORT 10: -lmkl_em64t -lmkl_core -liomp5  -lpthread
#IFORT 11,12:
    MATHLIB = -Wl,--start-group  $(MKLROOT)/lib/intel64/libmkl_intel_lp64.a $(MKLROOT)/lib/intel64/libmkl_sequential.a $(MKLROOT)/lib/intel64/libmkl_core.a -Wl,--end-group -lpthread -lm
 endif
 ifeq ($(OS),OSX)
    MATHLIB = $(MKLROOT)/lib/libmkl_intel_lp64.a $(MKLROOT)/lib/libmkl_sequential.a $(MKLROOT)/lib/libmkl_core.a -lpthread -lm
 endif
else
MATHLIB	 = -llapack -lblas
endif

#SciFor LIBRARY:
SFLIBS 	   = -lscifor ${MATHLIB}
SFLIBS_DEB = -lscifor_deb ${MATHLIB}


ifdef SFFFTW3
SFLIBS      += -lfftw3
SFLIBDS_DEB += -lfftw3
SFMODS      += -I$(SFFFTW3)/include/
SFMODS_DEB  += -I$(SFFFTW3)/include/
endif


ifdef DISLIN
DSL         =  -L$(DISLIN)/lib  -ldislin
X11 	    =  -L/usr/lib -lX11 -lXt -lXext -lxcb -lX11-xcb -lXau -lXdmcp #-lxcb-xlib 
DLPLOT      =  -ldlplot
DLPLOT_DEB  =  -ldlplot_deb
DSL_LIBS    = ${DLPLOT} ${DSL} ${X11}
DSL_LIBS_DEB= ${DLPLOT_DEB} ${DSL} ${X11}
DSL_MODS= -I${DISLIN}/ifc
endif


#REVISION SOFTWARE GIT:
REV = $(shell git rev-parse HEAD)
VER='character(len=41),parameter :: revision = "$(REV)"' > revision.inc


###################################################################
#REPO:
# ifdef CUDADIR
# GPU     = -L$(SFLIB)/ 	-lgpu
# CUDA    = -L/opt/cuda_2.3/lib64 -lcufft -lcuda -lcudart -lcublas
# CUDA_S  = -L/opt/cuda_sdk_2.3/lib -lcutil
# GPU_LIBS= ${GPU} ${CUDA_STATIC} ${CUDA} 
# endif


# ifdef FGSLDIR
# FGSL	  = -lfgsl_intel -lgsl -lgslcblas
# FGSL_MODS = -I$(FGSLDIR)/include
# SFLIBS += ${FGSL}
# SFMODS += ${FGSL_MODS}
# endif


