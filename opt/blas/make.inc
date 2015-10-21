SHELL = /bin/sh
PLAT = gnu
OBJ_INSTALL=/opt/scifor/opt/blas/src/obj_gnu
#  
#  Modify the FORTRAN and OPTS definitions to refer to the
#  compiler and desired compiler options for your machine.  NOOPT
#  refers to the compiler options desired when NO OPTIMIZATION is
#  selected.  Define LOADER and LOADOPTS to refer to the loader
#  and desired load options for your machine.
#
FORTRAN = gfortran
OPTS     = -O2 -funroll-all-loops -static 
DRVOPTS  = 
NOOPT    = -O0
LOADER   = gfortran
LOADOPTS = 

# Timer for the SECOND and DSECND routines
# TIMER     = NONE
#
#  The archiver and the flag(s) to use when building archive (library)
#  If you system has no ranlib, set RANLIB = echo.
#
ARCH     = ar
ARCHFLAGS= cr
RANLIB   = ranlib
#
#  The location of BLAS library for linking the testing programs.
#  The target's machine-specific, optimized BLAS library should be
#  used whenever possible.
#
BLASLIB      = /Users/amaricci/opt/scifor/gnu/lib/libblas.a
#
#  Location of the extended-precision BLAS (XBLAS) Fortran library
#  used for building and testing extended-precision routines.  The
#  relevant routines will be compiled and XBLAS will be linked only if
#  USEXBLAS is defined.
#
# USEXBLAS    = Yes
XBLASLIB     =
# XBLASLIB    = -lxblas
#
#  Names of generated libraries.
#
LAPACKLIB    = lapackgnu.a
TMGLIB       = tmglibgnu.a
EIGSRCLIB    = eigsrcgnu.a
LINSRCLIB    = linsrcgnu.a

