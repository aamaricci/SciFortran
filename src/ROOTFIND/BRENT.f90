MODULE BRENT
  implicit none
  ! private

  ! public :: fzero
  ! public :: zbrent

contains

  function fzero(func,a,b)
    interface
       function func(x)
         real(8),intent(in) :: x
         real(8)            :: func
       end function func
    end interface
    real(8),intent(in) :: a,b
    real(8)            :: fzero
    real(8)            :: tol
    !real(8),external   :: zbrent
    tol=epsilon(a)
    fzero = zbrent(func,a,b,tol)
  end function fzero

  include "zbrent.f90"

END MODULE BRENT
