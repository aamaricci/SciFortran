!###############################################################
!     PURPOSE  : Contains common routines
!     AUTHORS  : Adriano Amaricci
!###############################################################
module MPT_SC_MATS
  USE IPT_VARS_GLOBAL
  implicit none
  private

  integer                         :: i
  integer,save                    :: loop=1
  type(matsubara_gf),dimension(2) :: fg0,sigma
  real(8)                         :: n,n0,delta,delta0
  type(matsubara_gf)              :: calG11,calG22,calF
  public                          :: solve_mpt_sc_matsubara

contains


  !+-------------------------------------------------------------------+
  !PURPOSE  : 
  !+-------------------------------------------------------------------+
  function solve_mpt_sc_matsubara(fg0_,n_,n0_,delta_,delta0_) result(sigma_)
    complex(8)                :: fg0_(2,L),sigma_(2,L)
    real(8)                   :: n_,n0_,delta_,delta0_
    if(loop==1)then
       if(.not.fg0(1)%status)call allocate_gf(fg0(1),L)
       if(.not.fg0(2)%status)call allocate_gf(fg0(2),L)
       if(.not.sigma(1)%status)call allocate_gf(sigma(1),L)
       if(.not.sigma(2)%status)call allocate_gf(sigma(2),L)
       if(.not.calG11%status)call allocate_gf(calG11,L)
       if(.not.calG22%status)call allocate_gf(calG22,L)
       if(.not.calF%status)call allocate_gf(calF,L)
    endif
    fg0(1)%iw = fg0_(1,:); fg0(2)%iw = fg0_(2,:)
    n     = n_     ; n0    = n0_
    delta = delta_ ; delta0= delta0_
    call simpurity
    sigma_(1,:)=sigma(1)%iw ; sigma_(2,:)=sigma(2)%iw
    loop=loop+1
  end function solve_mpt_sc_matsubara
  !*******************************************************************
  !*******************************************************************
  !*******************************************************************






  !+-------------------------------------------------------------------+
  !PURPOSE  : Evaluate the 2^nd-order perturbation theory self-energy
  !+-------------------------------------------------------------------+
  subroutine simpurity
    real(8)                 :: A,B
    calG11%iw =  fg0(1)%iw
    calG22%iw = -conjg(fg0(1)%iw)
    calF%iw   =  fg0(2)%iw
    call fftgf_iw2tau(calG11%iw,calG11%tau,beta)
    call fftgf_iw2tau(calG22%iw,calG22%tau,beta)
    call fft_iw2tau(calF%iw,calF%tau,beta,L) !; calF%tau =-calF%tau
    forall(i=0:L)sigma(1)%tau(i)  =  U**2*(calG11%tau(i)*calG22%tau(i) -&
         calF%tau(i)**2)*calG22%tau(L-i)
    forall(i=0:L) sigma(2)%tau(i) =  U**2*(calG11%tau(i)*calG22%tau(i) -&
         calF%tau(i)**2)*calF%tau(i)
    call fftgf_tau2iw(sigma(1)%tau,sigma(1)%iw,beta)
    call fftgf_tau2iw(sigma(2)%tau,sigma(2)%iw,beta)

    A=U**2*n*(1.d0-n)-delta**2
    B=U**2*n0*(1.d0-n0)-delta**2
    sigma(1)%iw = u*(n-0.5d0) + sigma(1)%iw*A/B
    sigma(2)%iw = -delta     + sigma(2)%iw
    return
  end subroutine Simpurity
  !******************************************************************
  !******************************************************************
  !******************************************************************


end module MPT_SC_MATS
