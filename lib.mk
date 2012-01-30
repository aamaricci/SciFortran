include $(HOME)/lib/opt.mk

MODS 	= -I$(MYINCLUDE)
MODS_DEB= -I$(MYINCLUDE)/debug

#MYLIBS:
FFTGF     =  -lfftgf
IOTOOLS   =  -liotools
COMVARS   =  -lcomvars
LIST	  =  -llist
ROOTFIND  =  -lroot
GREENFUNX =  -lgfunx
LATTICE   =  -llattice
TOOLS     =  -ltools
SPLINE    =  -lspline
INTEGRATE =  -lintegrate
RANDOM    =  -lrandom
MATRIX    =  -lmatrix
CHRONO    =  -lchrono
STRINGPACK=  -lstringpack
STATISTICS=  -lstat
#FFTPACK   =  -L$(LIBOPT)/fftpack -lfftpack


ifdef MKLDIR
MATHLIB  =   -lmkl_em64t -lmkl_core -liomp5  -lpthread
else
MATHLIB	 = -L$(LIBOPT)/lapack -L$(LIBOPT)/blas -llapack_intel -lblas_intel
endif


ifdef FGSLDIR
FGSL	  = -lfgsl_intel -L$(GSLDIR)/lib -lgsl -lgslcblas
MODS 	 += -I$(FGSLDIR)/include
MODS_DEB += -I$(FGSLDIR)/include
endif


ifdef FFTW3DIR
FFTW3	 =  -lfftw3
MODS 	 += -I$(FFTW3DIR)/include/
MODS_DEB += -I$(FFTW3DIR)/include/
endif


ifdef DISLIN
DSL         =  -L$(DISLIN)/lib  -ldislin
X11 	    =  -L$(XLIB) -lX11 -lXt -lXext -lxcb -lX11-xcb -lXau -lXdmcp #-lxcb-xlib 
DLPLOT      =  -ldlplot
DLPLOT_DEB  =  -ldlplot_deb
DSL_LIBS    = ${DLPLOT} ${DSL} ${X11}
DSL_LIBS_DEB= ${DLPLOT_DEB} ${DSL} ${X11}
DSL_MODS= -I${DISLIN}/ifc
endif


ifdef CUDADIR
GPU     = -L$(MYLIB)/ 	-lgpu
CUDA    = -L/opt/cuda_2.3/lib64 -lcufft -lcuda -lcudart -lcublas
CUDA_S  = -L/opt/cuda_sdk_2.3/lib -lcutil
GPU_LIBS= ${GPU} ${CUDA_STATIC} ${CUDA} 
endif

##################################################################################################################################################################################

LIBS 	 = ${LIST} ${FFTGF} ${TOOLS} ${ROOTFIND}  ${SPLINE} ${COMVARS} ${INTEGRATE} ${RANDOM} ${GREENFUNX} ${MATRIX} ${IOTOOLS} ${CHRONO} ${STRINGPACK} ${STATISTICS} ${LATTICE} ${FFTW3} ${MATHLIB} ${FGSL}

LIBS_DEB = ${LIST}_deb ${FFTGF}_deb ${TOOLS}_deb ${ROOTFIND}_deb ${SPLINE}_deb ${COMVARS}_deb ${INTEGRATE}_deb ${RANDOM}_deb ${GREENFUNX}_deb ${MATRIX}_deb ${IOTOOLS}_deb ${CHRONO}_deb ${STRINGPACK}_deb ${STATISTICS}_deb ${LATTICE}_deb ${FFTW3} ${MATHLIB} ${FGSL}






