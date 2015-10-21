FC=gfortran
FFLAGS=-O2 -funroll-all-loops -static
PLAT=gnu
OBJ_INSTALL=/opt/scifor/opt/arpack/src/obj_gnu
BLAS_INSTALL=/opt/scifor/opt/arpack/blas/blas_gnu
LAPACK_INSTALL=/opt/scifor/opt/arpack/lapack/lapack_gnu
UTIL_INSTALL=/opt/scifor/opt/arpack/util/util_gnu

home=/opt/scifor/opt/arpack

BLASdir      = $(home)/blas
LAPACKdir    = $(home)/lapack
UTILdir      = $(home)/util
SRCdir       = $(home)/src
DIRS   = $(BLASdir) $(LAPACKdir) $(UTILdir) $(SRCdir)


#  The name of the libraries to be created/linked to
ARPACKLIB  = /Users/amaricci/opt/scifor/gnu/lib/libarpack.a
LAPACKLIB = 
BLASLIB = 
ALIBS =  $(ARPACKLIB) $(LAPACKLIB) $(BLASLIB) 


#  Make our own suffixes list.
#.SUFFIXES:
#.SUFFIXES:.f .o

#
#  Default command.
#
.DEFAULT:
	@$(ECHO) "Unknown target $@, try:  make help"

#
#  Command to build .o files from .f files.
#
.f.o:
	@$(ECHO) Making $@ from $<
	@$(FC) -c $(FFLAGS) $<

#
#  Various compilation programs and flags.
#  You need to make sure these are correct for your system.
LDFLAGS = 
CD	= cd
AR      = ar 
ARFLAGS  = cvq
CHMOD	= chmod
CHFLAGS	= -f
COMPRESS= compress

CP	= cp
CPP	 = /lib/cpp
CPPFLAGS =
ECHO	 = echo
LN	 = ln
LNFLAGS	 = -s
#MAKE	 = /bin/make
MKDIR	 = mkdir
MDFLAGS	 = -p
MV	 = mv
MVFLAGS	 = -f
RM	 = rm
RMFLAGS  = -f
SHELL	 = /bin/sh
TAR	 = tar
RANLIB   = ranlib


help:
	@$(ECHO) "usage: make ?"

