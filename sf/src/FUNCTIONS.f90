!###############################################################
! PROGRAM  : FUNCTIONS
! PURPOSE  : give access to standard functions not implicit in F
!###############################################################  
module FUNCTIONS
  USE COMMON_VARS
  implicit none
  private


  !FUNCTIONS:
  public :: heaviside
  public :: step
  public :: fermi
  interface sgn
     module procedure i_sgn,d_sgn
  end interface sgn
  public :: sgn


  !ERROR FUNCS
  public :: wfun         !complex error function (Faddeeva function)
  public :: zerf   

contains


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
    include "functions_wofz.f90"
  end function wfun


  !*******************************************************************
  !*******************************************************************
  !*******************************************************************


  !Double precision complex argument Error function
  include "functions_zerf.f90"


END module FUNCTIONS
