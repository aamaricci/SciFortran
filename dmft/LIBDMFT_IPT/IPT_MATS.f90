!###############################################################
!     PROGRAM  : FUNX_MATS
!     TYPE     : Module
!     PURPOSE  : Contains common routines
!     AUTHORS  : Adriano Amaricci
!###############################################################
module IPT_MATS
  USE IPT_VARS_GLOBAL
  implicit none
  private

  integer            :: i,LL
  integer,save       :: loop=1 
  type(matsubara_gf) :: fg0,sigma
  real(8)            :: n,n0,xmu0
  public             :: solve_mpt_matsubara
  public             :: solve_ipt_matsubara

contains

  !+-------------------------------------------------------------------+
  !PURPOSE  : 
  !COMMENT  : 
  !+-------------------------------------------------------------------+
  function solve_ipt_matsubara(fg0_) result(sigma_)
    complex(8),dimension(:) :: fg0_
    complex(8),dimension(size(fg0_)) :: sigma_
    LL=size(fg0_)
    if(loop==1)then
       if(.not.fg0%status)call allocate_gf(fg0,LL)
       if(.not.sigma%status)call allocate_gf(sigma,LL)
    endif
    fg0%iw=fg0_
    call fftgf_iw2tau(fg0%iw,fg0%tau,beta)
    forall(i=0:LL)sigma%tau(i)=U**2*(fg0%tau(i))**2*fg0%tau(LL-i)
    call fftgf_tau2iw(sigma%tau,sigma%iw,beta)
    sigma_=sigma%iw
    loop=loop+1
  end function solve_ipt_matsubara



  !*******************************************************************
  !*******************************************************************
  !*******************************************************************



  !+-------------------------------------------------------------------+
  !PURPOSE  : 
  !+-------------------------------------------------------------------+
  function solve_mpt_matsubara(fg0_,n_,n0_,xmu0_) result(sigma_)
    complex(8),dimension(:)    :: fg0_
    complex(8),dimension(size(fg0_))    :: sigma_
    real(8)                    :: n_,n0_,xmu0_
    real(8)                    :: A,B
    LL=size(fg0_)
    if(loop==1) then
       if(.not.fg0%status)call allocate_gf(fg0,LL)
       if(.not.sigma%status)call allocate_gf(sigma,LL)
    endif
    fg0%iw=fg0_ ; n=n_ ; n0=n0_ ; xmu0=xmu0_
    call fftgf_iw2tau(fg0%iw,fg0%tau,beta)
    forall(i=0:LL)sigma%tau(i)=U**2*(fg0%tau(i))**2*fg0%tau(LL-i)
    call fftgf_tau2iw(sigma%tau,sigma%iw,beta)
    call get_A ; call get_B
    sigma_ = U*(n-0.5d0) + A*sigma%iw/(1.d0-B*sigma%iw)
    loop=loop+1
  contains
    subroutine get_A
      real(8) :: A1,A2
      A1=n*(1.d0-n)
      A2=n0*(1.d0-n0)
      A=A1/A2
    end subroutine get_A
    subroutine get_B
      real(8) :: B1,B2
      B1 = (xmu0-xmu) + U*(1.d0-2.d0*n)
      B2 = n0*(1.d0-n0)*U**2
      B=B1/B2
    end subroutine get_B
  end function solve_mpt_matsubara



end module IPT_MATS
