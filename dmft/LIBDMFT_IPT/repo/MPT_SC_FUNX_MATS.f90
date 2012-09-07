!###############################################################
!     PROGRAM  : FUNX_IPT
!     TYPE     : Module
!     PURPOSE  : Contains common routines
!     AUTHORS  : Adriano Amaricci
!###############################################################
module MPT_SC_FUNX_MATS
  USE IPT_VARS_GLOBAL
  implicit none
  private

  integer,save                    :: loop=1
  type(matsubara_gf),dimension(2) :: fg,fg0,sigma
  real(8)                         :: n,n0,delta,delta0
  type(matsubara_gf)              :: calG11,calG22,calF
  logical                         :: iprint=.true.

  public :: solve_mpt_sc_matsubara

contains
  !+-------------------------------------------------------------------+
  !PROGRAM  : 
  !TYPE     : Subroutine
  !PURPOSE  : 
  !COMMENT  : 
  !+-------------------------------------------------------------------+
  function solve_mpt_sc_matsubara(fg0_,fg_,n_,n0_,delta_,delta0_,label_,iprint_) result(sigma_)
    complex(8)                :: fg0_(2,L),sigma_(2,L),fg_(2,L)
    real(8)                   :: n_,n0_,delta_,delta0_
    character(len=*),optional :: label_
    logical,optional          :: iprint_
    if(present(label_))label=trim(adjustl(trim(label_)))
    if(present(iprint_))iprint=iprint_
    if(loop==1)then
       if(.not.fg(1)%status)call allocate_gf(fg(1),L)
       if(.not.fg0(1)%status)call allocate_gf(fg0(1),L)
       if(.not.sigma(1)%status)call allocate_gf(sigma(1),L)
       if(.not.fg(2)%status)call allocate_gf(fg(2),L)
       if(.not.fg0(2)%status)call allocate_gf(fg0(2),L)
       if(.not.sigma(2)%status)call allocate_gf(sigma(2),L)
       if(.not.calG11%status)call allocate_gf(calG11,L)
       if(.not.calG22%status)call allocate_gf(calG22,L)
       if(.not.calF%status)call allocate_gf(calF,L)
    endif
    fg0(1)%iw = fg0_(1,:); fg0(2)%iw = fg0_(2,:)
    fg(1)%iw  = fg_(1,:) ; fg(2)%iw  = fg_(2,:)

    n     = n_                !-2.0*real(fg(1)%tau(L))
    n0    = n0_               ! -2.0*real(fg0(1)%tau(L))
    delta = delta_            !u*real(fg(2)%tau(0))
    delta0= delta0_           !u*real(fg0(2)%tau(0))

    call simpurity
    sigma_(1,:)=sigma(1)%iw ; sigma_(2,:)=sigma(2)%iw
    if(iprint)call print_out()
    loop=loop+1
  end function solve_mpt_sc_matsubara
  !*******************************************************************
  !*******************************************************************
  !*******************************************************************






  !+-------------------------------------------------------------------+
  !PROGRAM  : SIMPURITY
  !TYPE     : function
  !PURPOSE  : Evaluate the 2^nd-order perturbation theory self-energy
  !+-------------------------------------------------------------------+
  subroutine simpurity
    real(8)                 :: A,B,nn,nn0

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

    nn =n!-0.5d0
    nn0=n0!-0.5d0
    A=U**2*nn*(1.d0-nn)-delta**2
    B=U**2*nn0*(1.d0-nn0)-delta**2
    sigma(1)%iw = u*(n-0.5d0) + sigma(1)%iw*A/B
    sigma(2)%iw = -delta     + sigma(2)%iw
    return
  end subroutine Simpurity
  !******************************************************************
  !******************************************************************
  !******************************************************************








  !+-------------------------------------------------------------------+
  !PROGRAM  : 
  !TYPE     : Subroutine
  !PURPOSE  : 
  !COMMENT  : 
  !+-------------------------------------------------------------------+
  subroutine print_out()
    real(8) :: nimp!,delta
    call fftgf_iw2tau(fg(1)%iw,fg(1)%tau,beta)
    call fft_iw2tau(fg(2)%iw,fg(2)%tau,beta,L)
    nimp = -2.0*real(fg(1)%tau(L))
    ! delta=  u*real(fg(2)%tau(0))

    call splot("nVSiloop.ipt"//trim(adjustl(trim(label))),iloop,nimp,append=TT)
    call splot("deltaVSiloop.ipt"//trim(adjustl(trim(label))),iloop,delta,append=TT)

    call splot("Sigma_iw.ipt"//trim(adjustl(trim(label))),wm,sigma(1)%iw,append=printf)
    call splot("Self_iw.ipt"//trim(adjustl(trim(label))),wm,sigma(2)%iw,append=printf)
    call splot("G_iw.ipt"//trim(adjustl(trim(label))),wm,fg(1)%iw,append=printf)
    call splot("F_iw.ipt"//trim(adjustl(trim(label))),wm,fg(2)%iw,append=printf)

    call splot("Sigma_tau.ipt"//trim(adjustl(trim(label))),tau,sigma(1)%tau,append=printf)
    call splot("Self_tau.ipt"//trim(adjustl(trim(label))),tau,sigma(2)%tau,append=printf)
    call splot("G_tau.ipt"//trim(adjustl(trim(label))),tau,fg(1)%tau,append=printf)
    call splot("F_tau.ipt"//trim(adjustl(trim(label))),tau,fg(2)%tau,append=printf)

    call splot("calG11_iw.ipt"//trim(adjustl(trim(label))),wm,calG11%iw,append=printf)
    call splot("calG22_iw.ipt"//trim(adjustl(trim(label))),wm,calG22%iw,append=printf)
    call splot("calF_iw.ipt"//trim(adjustl(trim(label))),wm,calF%iw,append=printf)
    return
  end subroutine print_out
  !*******************************************************************
  !*******************************************************************
  !*******************************************************************








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

end module MPT_SC_FUNX_MATS
