include "minpack.f90"
MODULE WRAP_MINPACK
  implicit none
  private
  public :: fsolve,ffsolve
contains
  subroutine fsolve(func,x,tol,info)
    real(8),dimension(:)       :: x
    real(8),dimension(size(x)) :: fvec
    integer                    :: n
    real(8),optional           :: tol
    integer,optional           :: info
    real(8)                    :: tol_
    integer                    :: info_
    external func
    tol_ = 1.d-15;if(present(tol))tol_=tol
    n=size(x)
    call hybrd1(func,n,x,fvec,tol_,info_)
    if(present(info))info=info_
  end subroutine fsolve

  function ffsolve(func,x,tol,info)
    real(8),dimension(:)       :: x
    real(8),dimension(size(x)) :: fvec
    real(8),dimension(size(x)) :: ffsolve
    integer                    :: n
    real(8),optional           :: tol
    integer,optional           :: info
    real(8)                    :: tol_
    integer                    :: info_
    external func
    tol_ = 1.d-15;if(present(tol))tol_=tol
    n=size(x)
    call hybrd1(func,n,x,fvec,tol_,info_)
    if(present(info))info=info_
    ffsolve=x
  end function ffsolve
END MODULE WRAP_MINPACK

MODULE ZEROS
  USE BROYDEN
  USE BRENT
  USE WRAP_MINPACK
  implicit none
END MODULE ZEROS
