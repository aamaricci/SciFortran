  include "MOD_NR.f90"
  !###############################################################
  ! PROGRAM  : MOD_FINTER
  ! TYPE     : Module
  ! PURPOSE  : contains a function definition that can be used
  ! to interpolate any given input function with polinomial of 
  ! order N.
  !to use it is enough to allocate and define the finterXYZ variables.
  !###############################################################
  module MOD_FINTER
    USE MOD_NR
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


    ! function finter(x)
    !   real(8) :: finter
    !   real(8) :: x,y,dy
    !   integer :: itmp,k
    !   integer :: n
    !   n=finterN    !order of polynomial interpolation
    !   finter=0.d0  !Interpolate discrete Aw to continious F 
    !   if (-finterImin.eq.finterImax) then
    !      call locate(finterX(0:finterImax),finterImax+1,abs(x),itmp)
    !   else
    !      call locate(finterX(FinterImin:finterImax),finterImax-finterImin,x,itmp)
    !   endif
    !   itmp=itmp-1
    !   k = (itmp-(n-1)/2)          !,ifs+1-n) min(.. max(itmp-(n-1)/2,0)!
    !   if (k < finterImin)k=finterImin
    !   if (k+n < finterImax) then
    !      if ((x >= 0.d0).or.(finterImin == 0)) then
    !         call polint(finterX(k:k+n),finterF(k:k+n),n,x,y,dy)
    !      else
    !         call polint(finterX(-k-n:-k),finterF(-k-n:-k),n,x,y,dy)
    !      endif
    !      finter=y
    !   else
    !      finter=0.d0
    !   endif
    !   return
    ! end function finter



  end module MOD_FINTER
