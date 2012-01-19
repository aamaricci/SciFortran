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
  interface operator(+)
     module procedure add
  end interface operator(+)
  interface operator(-)
     module procedure subtract
  end interface operator(-)
  interface operator(*)
     module procedure prodL, prodR
  end interface operator(*)
  interface operator(.dot.)
     module procedure dot
  end interface operator(.dot.)
  interface assignment(=)
     module procedure vector_eq, scalar_eq
  end interface assignment(=)

  type(vect2D),parameter,public :: VZero=vect2D(0.d0,0.d0)
  type(vect2D),parameter,public :: VOne=vect2D(1.d0,1.d0)
  type(vect2D),parameter,public :: Xver=vect2D(1.d0,0.d0)
  type(vect2D),parameter,public :: Yver=vect2D(0.d0,1.d0)


  public :: operator(+)
  public :: operator(-)
  public :: operator(*)
  public :: operator(.dot.)
  public :: assignment(=)
  public :: modulo,add,subtract,prodL,prodR,dot,scalar_eq,vector_eq

contains
  !+-----------------------------------------------------------------+
  !PROGRAM  : 2D VECTOR ALGEBRA:
  !TYPE     : Subroutine
  !PURPOSE  : instruct the machine with the necessary 2D vectors algebra 
  !+-----------------------------------------------------------------+
  elemental function add(V,W)
    type(vect2D),intent(in)  :: V,W
    type(vect2D) :: add
    add%x = V%x + W%x
    add%y = V%y + W%y
  end function add
  !-------------------------------------

  !+-----------------------------------------------------------------+
  elemental function modulo(V)
    type(vect2D),intent(in)  :: V
    real(8) :: modulo
    modulo = sqrt(V%x**2 + V%y**2)
  end function modulo


  !-------------------------------------
  elemental function subtract(V,W)
    type(vect2D),intent(in)  :: V,W
    type(vect2D) :: subtract
    subtract%x = V%x - W%x
    subtract%y = V%y - W%y
  end function subtract
  !-------------------------------------

  !-------------------------------------
  elemental function prodL(C,V)
    real(8),intent(in) :: C
    type(vect2D),intent(in)  :: V
    type(vect2D) :: prodL
    prodL%x = C*V%x 
    prodL%y = C*V%y
  end function prodL
  !-------------------------------------

  !-------------------------------------
  elemental function prodR(V,C)
    real(8),intent(in) :: C
    type(vect2D),intent(in)  :: V
    type(vect2D) :: prodR
    prodR%x = V%x*C 
    prodR%y = V%y*C
  end function prodR
  !-------------------------------------

  !-------------------------------------
  elemental function dot(V,W)
    type(vect2D),intent(in)  :: V,W
    real(8) :: dot
    dot = V%x*W%x + V%y*W%y
  end function dot
  !-------------------------------------

  !-------------------------------------
  elemental subroutine scalar_eq(V,C)
    real(8),intent(in)  :: C
    type(vect2D),intent(inout) :: V
    V%x = C
    V%y = C
  end subroutine  scalar_eq
  !-------------------------------------

  !-------------------------------------
  elemental subroutine vector_eq(V,W)
    type(vect2D),intent(in)  :: W
    type(vect2D),intent(inout) :: V
    V%x = W%x
    V%y = W%y
  end subroutine  vector_eq
  !******************************************************************
  !******************************************************************
  !******************************************************************
end module VECTORS
