!###############################################################
!     PROGRAM  : VECTORS
!     TYPE     : Module
!     PURPOSE  : Setup a class for vector algebra. 
!     TODO     : exploit fortran2003 extensions to include higher 
! and lower dimensionalities with real class properties
!     AUTHORS  : Adriano Amaricci
!###############################################################
module VECTORS
  implicit none
  private

  !2D space and its algebra:
  !=========================================================
  type,public :: vect2D
     real(8) :: x
     real(8) :: y
  end type vect2D
  type,public :: vect3D
     real(8) :: x
     real(8) :: y
     real(8) :: z
  end type vect3D
  interface operator(+)
     module procedure add2D,add3D
  end interface operator(+)
  interface operator(-)
     module procedure subtract2D,subtract3D
  end interface operator(-)
  interface operator(*)
     module procedure prodL2D, prodR2D, prodL3D, prodR3D
  end interface operator(*)
  interface operator(.dot.)
     module procedure dot2D,dot3D
  end interface operator(.dot.)
  interface assignment(=)
     module procedure vector_eq2D, scalar_eq2D,vector_eq3D, scalar_eq3D
  end interface assignment(=)
  interface modulo
     module procedure modulo2d,modulo3d
  end interface modulo

  type(vect2D),parameter,public :: VZero=vect2D(0.d0,0.d0)
  type(vect2D),parameter,public :: VOne=vect2D(1.d0,1.d0)
  type(vect2D),parameter,public :: Xver=vect2D(1.d0,0.d0)
  type(vect2D),parameter,public :: Yver=vect2D(0.d0,1.d0)

  type(vect3D),parameter,public :: VZero3d=vect3D(0.d0,0.d0,0.d0)
  type(vect3D),parameter,public :: VOne3d=vect3D(1.d0,1.d0,1.d0)
  type(vect3D),parameter,public :: Xver3d=vect3D(1.d0,0.d0,0.d0)
  type(vect3D),parameter,public :: Yver3d=vect3D(0.d0,1.d0,0.d0)
  type(vect3D),parameter,public :: Zver3d=vect3D(0.d0,0.d0,1.d0)

  public :: operator(+)
  public :: operator(-)
  public :: operator(*)
  public :: operator(.dot.)
  public :: assignment(=)
  public :: modulo

contains

  !+-----------------------------------------------------------------+
  !PROGRAM  : 2D VECTOR ALGEBRA:
  !PURPOSE  : generate the necessary 2D vectors algebra 
  !+-----------------------------------------------------------------+
  include 'vectorial_algebra_2d.f90'

  !+-----------------------------------------------------------------+
  !PROGRAM  : 3D VECTOR ALGEBRA:
  !PURPOSE  : generate the necessary 2D vectors algebra 
  !+-----------------------------------------------------------------+
  include 'vectorial_algebra_3d.f90'


end module VECTORS
