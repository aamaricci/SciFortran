!###############################################################
!     PURPOSE  : Contains common routines
!     AUTHORS  : Adriano Amaricci
!###############################################################
module IPT_SC_MATS
  USE IPT_VARS_GLOBAL
  implicit none
  private

  integer                         :: i
  integer,save                    :: loop=1
  type(matsubara_gf),dimension(2) :: fg0,sigma
  type(matsubara_gf)              :: calG11,calG22,calF
  real(8)                         :: n,n0,delta,delta0
  public                          :: solve_ipt_sc_matsubara
  public                          :: solve_mpt_sc_matsubara

contains


  !+-------------------------------------------------------------------+
  !PURPOSE  : 
  !+-------------------------------------------------------------------+
  function solve_ipt_sc_matsubara(fg0_,delta_) result(sigma_)
    complex(8),dimension(2,L)  :: fg0_,sigma_
    real(8)                    :: delta_
    if(loop==1)then
       if(.not.fg0(1)%status)call allocate_gf(fg0(1),L)
       if(.not.fg0(2)%status)call allocate_gf(fg0(2),L)
       if(.not.sigma(1)%status)call allocate_gf(sigma(1),L)
       if(.not.sigma(2)%status)call allocate_gf(sigma(2),L)
       if(.not.calG11%status)call allocate_gf(calG11,L)
       if(.not.calG22%status)call allocate_gf(calG22,L)
       if(.not.calF%status)call allocate_gf(calF,L)
    endif
    fg0(1)%iw = fg0_(1,:); fg0(2)%iw = fg0_(2,:) ; delta = delta_
    call simpurity
    sigma_(1,:)=sigma(1)%iw
    sigma_(2,:)=sigma(2)%iw - delta
    loop=loop+1
  end function solve_ipt_sc_matsubara



  !*******************************************************************
  !*******************************************************************
  !*******************************************************************



  !+-------------------------------------------------------------------+
  !PURPOSE  : 
  !+-------------------------------------------------------------------+
  function solve_mpt_sc_matsubara(fg0_,n_,n0_,delta_,delta0_) result(sigma_)
    complex(8)                :: fg0_(2,L),sigma_(2,L)
    real(8)                   :: n_,n0_,delta_,delta0_,A,B
    if(loop==1)then
       if(.not.fg0(1)%status)call allocate_gf(fg0(1),L)
       if(.not.fg0(2)%status)call allocate_gf(fg0(2),L)
       if(.not.sigma(1)%status)call allocate_gf(sigma(1),L)
       if(.not.sigma(2)%status)call allocate_gf(sigma(2),L)
       if(.not.calG11%status)call allocate_gf(calG11,L)
       if(.not.calG22%status)call allocate_gf(calG22,L)
       if(.not.calF%status)call allocate_gf(calF,L)
    endif
    fg0(1)%iw = fg0_(1,:) ; fg0(2)%iw = fg0_(2,:)
    n     = n_     ; n0    = n0_
    delta = delta_ ; delta0= delta0_
    call simpurity
    A=U**2*n*(1.d0-n)-delta**2
    B=U**2*n0*(1.d0-n0)-delta0**2
    sigma_(1,:) =-u*(n-0.5d0) + sigma(1)%iw*A/B
    sigma_(2,:) =-delta       + sigma(2)%iw*A/B
    loop=loop+1
  end function solve_mpt_sc_matsubara




  !*******************************************************************
  !*******************************************************************
  !*******************************************************************


  !+-------------------------------------------------------------------+
  !PURPOSE  : Evaluate the 2^nd-order perturbation theory self-energy
  !+-------------------------------------------------------------------+
  subroutine simpurity
    calG11%iw =  fg0(1)%iw
    calG22%iw = -conjg(fg0(1)%iw)
    calF%iw   =  fg0(2)%iw
    call fftgf_iw2tau(calG11%iw,calG11%tau,beta)
    call fftgf_iw2tau(calG22%iw,calG22%tau,beta)
    call fftgf_iw2tau(calF%iw,calF%tau,beta,notail=.true.)!;calF%tau =-calF%tau

    forall(i=0:L)
       sigma(1)%tau(i)=  U**2*(calG11%tau(i)*calG22%tau(i) - calF%tau(i)**2)*calG22%tau(L-i)
       sigma(2)%tau(i)= -U**2*(calF%tau(i)**2 - calG11%tau(i)*calG22%tau(i))*calF%tau(i)
    end forall
    call fftgf_tau2iw(sigma(1)%tau,sigma(1)%iw,beta)
    call fftgf_tau2iw(sigma(2)%tau,sigma(2)%iw,beta)
  end subroutine Simpurity


  !******************************************************************
  !******************************************************************
  !******************************************************************




  ! !+-------------------------------------------------------------------+
  ! !PROGRAM  : 
  ! !TYPE     : Subroutine
  ! !PURPOSE  : 
  ! !COMMENT  : 
  ! !+-------------------------------------------------------------------+
  ! subroutine get_gloc_matsubara
  !   complex(8) :: zeta,det
  !   forall(i=1:2)fg(i)%iw=zero
  !   do i=1,Lm
  !      zeta = xi*wm(i) + xmu - sigma(1)%iw(i)
  !      do ik=1,Lk
  !         det = (zeta-epsik(ik))*(conjg(zeta)-epsik(ik)) + (sigma(2)%iw(i))**2
  !         fg(1)%iw(i)=fg(1)%iw(i) + wt(ik)*(conjg(zeta)-epsik(ik))/det
  !         fg(2)%iw(i)=fg(2)%iw(i) + wt(ik)*(sigma(2)%iw(i))/det
  !      enddo
  !   enddo
  !   fg(2)%iw=real(fg(2)%iw,8) + xi*zero
  ! end subroutine get_gloc_matsubara
  ! !*******************************************************************
  ! !*******************************************************************
  ! !*******************************************************************





  ! !+-------------------------------------------------------------------+
  ! !PROGRAM  : 
  ! !TYPE     : Subroutine
  ! !PURPOSE  : 
  ! !COMMENT  : 
  ! !+-------------------------------------------------------------------+
  ! subroutine update_g0_matsubara(loop)
  !   integer                   :: loop
  !   complex(8)                :: det(Lm)
  !   if(loop==1)then
  !      call allocate_gf(calG11,Lm)
  !      call allocate_gf(calG22,Lm)
  !      call allocate_gf(calF,Lm)
  !   endif
  !   !get order parameter: delta=F(0)*U
  !   call fft_iw2tau(fg(2)%iw,fg(2)%tau,beta,Lm)
  !   deltaSC= u*fg(2)%tau(0)!;print*,deltaSC
  !   !calcola calG0^-1, calF0^-1 (WFs)
  !   det         = fg(1)%iw*conjg(fg(1)%iw) + (fg(2)%iw)**2
  !   fg0(1)%iw = conjg(fg(1)%iw)/det + sigma(1)%iw
  !   fg0(2)%iw = fg(2)%iw/det        + sigma(2)%iw - deltaSC
  !   !calG11 =  calG0^-1*/(abs(calG0^-1)^2 + (calF0^-1)^2)
  !   !calG22 = -calG0^-1/(abs(calG0^-1)^2 + (calF0^-1)^2)
  !   !calF   =  calF0^-1/(abs(calG0^-1)^2 + (calF0^-1)^2)
  !   det       = (fg0(1)%iw*conjg(fg0(1)%iw) + fg0(2)%iw**2)
  !   calG11%iw =  conjg(fg0(1)%iw)/det
  !   calG22%iw =  -fg0(1)%iw/det
  !   calF%iw   =  fg0(2)%iw/det
  !   return
  ! end subroutine update_g0_matsubara
  ! !*******************************************************************
  ! !*******************************************************************
  ! !*******************************************************************

end module IPT_SC_MATS
