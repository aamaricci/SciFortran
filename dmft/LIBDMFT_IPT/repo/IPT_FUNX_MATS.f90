!###############################################################
!     PROGRAM  : FUNX_MATS
!     TYPE     : Module
!     PURPOSE  : Contains common routines
!     AUTHORS  : Adriano Amaricci
!###############################################################
module IPT_FUNX_MATS
  USE IPT_VARS_GLOBAL
  implicit none
  private

  integer,save       :: loop=1 
  type(matsubara_gf) :: fg,fg0,sigma
  logical            :: iprint=.true.
  public :: solve_ipt_matsubara

contains
  !+-------------------------------------------------------------------+
  !PROGRAM  : 
  !TYPE     : Subroutine
  !PURPOSE  : 
  !COMMENT  : 
  !+-------------------------------------------------------------------+
  function solve_ipt_matsubara(fg0_,label_,iprint_) result(sigma_)
    complex(8),dimension(-L:L) :: fg0_,sigma_
    character(len=*),optional  :: label_
    logical,optional           :: iprint_
    if(present(label_))label=trim(adjustl(trim(label_)))
    if(present(iprint_))iprint=iprint_
    if(loop==1)then
       if(.not.fg%status)call allocate_gf(fg,L)
       if(.not.fg0%status)call allocate_gf(fg0,L)
       if(.not.sigma%status)call allocate_gf(sigma,L)
    endif
    fg0%iw=fg0_
    call simpurity
    sigma_=sigma%iw
    if(iprint)call print_out()
    loop=loop+1
    return
  end function solve_ipt_matsubara
  !*******************************************************************
  !*******************************************************************
  !*******************************************************************



  !+-------------------------------------------------------------------+
  !PROGRAM  : SIMPURITY
  !TYPE     : function
  !PURPOSE  : Evaluate the 2^nd-order perturbation theory self-energy
  !+-------------------------------------------------------------------+
  subroutine simpurity
    call fftgf_iw2tau(fg0%iw,fg0%tau,beta)
    forall(i=0:L)sigma%tau(i)=U**2*(fg0%tau(i))**2*fg0%tau(L-i)
    call fftgf_tau2iw(sigma%tau,sigma%iw,beta)
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
    real(8) :: nimp,zeta
    fg%iw = one/(one/fg0%iw - sigma%iw)
    call fftgf_iw2tau(fg%iw,fg%tau,beta)
    call fftgf_iw2tau(fg0%iw,fg0%tau,beta)
    nimp= -2.0*real(fg%tau(L))
    zeta=1.d0/(1.d0+abs(aimag(sigma%iw(1))/wm(1)))
    call splot("nVSiloop.ipt"//trim(adjustl(trim(label))),iloop,nimp,append=TT)
    call splot("zetaVSiloop.ipt"//trim(adjustl(trim(label))),iloop,zeta,append=TT)
    call splot("Sigma_iw.ipt"//trim(adjustl(trim(label))),wm,sigma%iw,append=printf)
    call splot("G_iw.ipt"//trim(adjustl(trim(label))),wm,fg%iw,append=printf)
    call splot("G0_iw.ipt"//trim(adjustl(trim(label))),wm,fg0%iw,append=printf)
    call splot("G_tau.ipt"//trim(adjustl(trim(label))),tau,fg%tau(0:L),append=printf)
    call splot("G0_tau.ipt"//trim(adjustl(trim(label))),tau,fg0%tau(0:L),append=printf)
    return
  end subroutine print_out
  !*******************************************************************
  !*******************************************************************
  !*******************************************************************

end module IPT_FUNX_MATS
