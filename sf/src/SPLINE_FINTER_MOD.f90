module SPLINE_FINTER_MOD
  USE SPLINE_NR_MOD
  implicit none
  private
  real(8), allocatable,public :: finterX(:)  !vector with frequencies
  real(8), allocatable,public :: finterF(:) !corresponding vector of functional values
  integer,public              :: finterImin,finterImax,finterN

  public :: finter

contains

  function finter(x)
    real(8) :: finter
    real(8) :: x,y,dy
    integer :: itmp,k
    integer :: n
    n=finterN    !order of polynomial interpolation
    finter=0.d0

    itmp=locate(finterX(FinterImin:finterImax),x)
    k=max(itmp-(N-1)/2,1)

    if (k < finterImin)k=finterImin

    if(k+n+1 <= finterImax)then
       call polint(finterX(k:k+n+1),finterF(k:k+n+1),x,y,dy)
    else
       call polint(finterX(k:finterImax),finterF(k:finterImax),x,y,dy)
    endif
    finter=y
    return
  end function finter


end module SPLINE_FINTER_MOD
