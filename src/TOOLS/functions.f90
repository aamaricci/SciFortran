
!+------------------------------------------------------------------+
!PURPOSE  : calculate the Heaviside  function
!+------------------------------------------------------------------+
pure function heaviside(x)
  real(8),intent(in) :: x
  real(8)            :: heaviside
  if(x < 0.d0) then
     heaviside = 0.0d0
  elseif(x==0.d0)then
     heaviside = 0.50d0
  else
     heaviside = 1.0d0
  endif
end function heaviside


!+------------------------------------------------------------------+
!PURPOSE  : calculate step function
!+------------------------------------------------------------------+
pure function step(x,origin)
  real(8),intent(in)          :: x
  logical,optional,intent(in) :: origin
  real(8)                     :: step
  logical                     :: w0
  step=0.d0
  w0=.true.;if(present(origin))w0=origin
  select case(w0)
  case (.true.)
     if(x>=0.d0)step=1.d0
  case (.false.)
     if(x>0.d0)step=1.d0
  end select
end function step



!*******************************************************************
!*******************************************************************
!*******************************************************************



!+-------------------------------------------------------------------+
!PURPOSE  : calculate the Fermi-Dirac distribution
!+-------------------------------------------------------------------+
elemental function fermi(x,beta)
  real(8),intent(in) :: x, beta 
  real(8)            :: fermi
  if(x*beta > 100.d0)then
     fermi=0.d0
     return
  endif
  fermi = 1.d0/(1.d0+exp(beta*x))
end function fermi


!*******************************************************************
!*******************************************************************
!*******************************************************************



!+-------------------------------------------------------------------+
!PURPOSE  : calculate the non-interacting dos for HYPERCUBIC lattice 
!+-------------------------------------------------------------------+
pure function dens_hyperc(x,t1)
  real(8),optional,intent(in) :: t1
  real(8),intent(in)          :: x
  REAL(8):: dens_hyperc,t1_,pi2,sqrt2
  pi2=2.d0*acos(-1.d0)
  sqrt2=sqrt(2.d0)
  t1_=sqrt2 ; if(present(t1))t1_=t1
  dens_hyperc = (1/(t1_*sqrt(pi2)))*exp(-(x**2)/(2.d0*t1_**2))
  return
end function dens_hyperc


!*******************************************************************
!*******************************************************************
!*******************************************************************



!+-------------------------------------------------------------------+
!PURPOSE:  evaluate the sign of a given number (I,R)
!+-------------------------------------------------------------------+
pure function i_sgn(x) result(sgn)
  integer,intent(in) :: x
  integer            :: sgn
  sgn=x/abs(x)
end function i_sgn
pure function d_sgn(x) result(sgn)
  real(8),intent(in) :: x
  real(8)            :: sgn
  sgn=x/abs(x)
end function d_sgn


!*******************************************************************
!*******************************************************************
!*******************************************************************





!+------------------------------------------------------------------+
!PURPOSE  : Evaluate the Complex Error Functions (Faddeeva function)
! w(x)=exp(-x^2)erfc(-ix)
!+------------------------------------------------------------------+
function wfun(z)
  complex(8):: z,wfun
  real(8)   :: x,y,u,v
  logical   :: flag
  x=real(z,8)
  y=aimag(z)
  call wofz(x,y,u,v,flag)
  wfun=cmplx(u,v)
contains
  include "wofz.f90"
end function wfun


!*******************************************************************
!*******************************************************************
!*******************************************************************


!Double precision complex argument Error function
include "zerf.f90"
