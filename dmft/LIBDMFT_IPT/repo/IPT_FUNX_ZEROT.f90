!###############################################################
!     PROGRAM  : FUNX_IPT
!     TYPE     : Module
!     PURPOSE  : Contains common routines
!     AUTHORS  : Adriano Amaricci
!###############################################################
module IPT_FUNX_ZEROT
  USE IPT_VARS_GLOBAL
  implicit none
  private
  integer,save                     :: loop=1 
  type(real_gf)                    :: fg,fg0,sigma
  real(8),allocatable              :: exa(:)
  logical                          :: iprint=.true.

  public :: solve_ipt_zerot

contains

  !+-------------------------------------------------------------------+
  !PROGRAM  : 
  !TYPE     : Subroutine
  !PURPOSE  : 
  !COMMENT  : 
  !+-------------------------------------------------------------------+
  function solve_ipt_zeroT(fg0_,label_,iprint_) result(sigma_)
    complex(8),dimension(-L:L) :: fg0_,sigma_
    real(8)                    :: ex
    character(len=*),optional  :: label_
    logical,optional           :: iprint_
    if(present(label_))label=trim(adjustl(trim(label_)))
    if(present(iprint_))iprint=iprint_
    if(loop==1)then
       if(.not.fg%status)call allocate_gf(fg,L)
       if(.not.fg0%status)call allocate_gf(fg0,L)
       if(.not.sigma%status)call allocate_gf(sigma,L)
       if(.not.allocated(exa))allocate(exa(-L:L))
       ex=-1.d0
       do i=-L,L
          ex=-ex
          exa(i)=ex
       enddo
    endif
    fg0%w=fg0_
    call simpurity
    sigma_=sigma%w
    if(iprint)call print_out()
    loop=loop+1
    return
  end function solve_ipt_zeroT
  !*******************************************************************
  !*******************************************************************
  !*******************************************************************


  !+-------------------------------------------------------------------+
  ! !PROGRAM  : 
  ! !TYPE     : Subroutine
  ! !PURPOSE  : 
  ! !COMMENT  : 
  ! !+-------------------------------------------------------------------+
  ! subroutine get_gloc_zeroT
  !   integer          :: i,ik
  !   complex(8)       :: zeta,sqroot,ifg
  !   real(8)          :: w,sig,sq
  !   fg%w=zero
  !   do i=0,L
  !      w=wr(i) ; zeta = w + xmu - sigma%w(i)
  !      sq=real(zeta) ; sig=1.d0 ; if(sq<0.d0)sig=-1.d0
  !      ifg=zero
  !      zeta=zeta+sig*xi*eps
  !      do ik=1,Lk
  !         ifg=ifg+wt(ik)/(zeta-epsik(ik))
  !      enddo
  !      fg%w(i)=ifg
  !      ! sqroot=sqrt(zeta**2-d**2)
  !      ! fg%w(i)=2.d0/(zeta+sig*sqroot)
  !   enddo
  ! end subroutine get_gloc_zeroT
  ! !*******************************************************************
  ! !*******************************************************************
  ! !*******************************************************************

  ! !+-------------------------------------------------------------------+
  ! !PROGRAM  : 
  ! !TYPE     : Subroutine
  ! !PURPOSE  : 
  ! !COMMENT  : 
  ! !+-------------------------------------------------------------------+
  ! subroutine update_g0_zeroT()
  !   forall(i=0:L)fg0%w(i) = one/(one/fg%w(i) + sigma%w(i))
  !   forall(i=1:L)fg0%w(-i)=-fg0%w(i)
  !   return
  ! end subroutine update_g0_zeroT
  ! !*******************************************************************
  ! !*******************************************************************
  ! !*******************************************************************


  !+-------------------------------------------------------------------+
  !PROGRAM  : SIMPURITY
  !TYPE     : function
  !PURPOSE  : Evaluate the 2^nd-order perturbation theory self-energy
  !+-------------------------------------------------------------------+
  subroutine simpurity
    call fftgf_rw2rt(fg0%w,fg0%t,L) ; fg0%t = fmesh/pi2*fg0%t
    forall(i=-L:L)sigma%t(i)=U**2*(fg0%t(i))**2*fg0%t(-i)
    call fftgf_rt2rw(sigma%t,sigma%w,L) ; sigma%w= dt*sigma%w
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
    real(8) :: nimp
    fg%w = one/(one/fg0%w - sigma%w)
    nimp = -2.d0*sum(aimag(fg%w(-L:L))*fermi(wr(-L:L),beta))*fmesh/pi
    call fftgf_rw2rt(fg%w,fg%t,L) ; fg%t = fmesh/pi2*fg%t
    call splot("nVSiloop.ipt"//trim(adjustl(trim(label))),loop,nimp,append=TT)
    call splot("DOS.ipt"//trim(adjustl(trim(label))),wr,abs(aimag(fg%w)),append=printf)
    call splot("Sigma_realw.ipt"//trim(adjustl(trim(label))),wr,sigma%w,append=printf)
    call splot("G_realw.ipt"//trim(adjustl(trim(label))),wr,fg%w,append=printf)
    call splot("G0_realw.ipt"//trim(adjustl(trim(label))),wr,fg0%w,append=printf)
    call splot("Sigma_t.ipt"//trim(adjustl(trim(label))),t,exa*sigma%t,append=printf)
    call splot("G_t.ipt"//trim(adjustl(trim(label))),t,exa*fg%t,append=printf)
    call splot("G0_t.ipt"//trim(adjustl(trim(label))),t,exa*fg0%t,append=printf)
    return
  end subroutine print_out
  !*******************************************************************
  !*******************************************************************
  !*******************************************************************

end module IPT_FUNX_ZEROT
