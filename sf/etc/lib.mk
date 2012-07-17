include $(SFDIR)/etc/opt.mk

#SciFor MODULES:
SFMODS 	= -I$(SFINCLUDE)
SFMODS_DEB= -I$(SFINCLUDE)/debug/


ifdef MKLDIR
MATHLIB  = -lmkl_gf_ilp64 -lmkl_intel_thread -lmkl_core -liomp5  -lpthread
#-lmkl_em64t -lmkl_core -liomp5  -lpthread
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


