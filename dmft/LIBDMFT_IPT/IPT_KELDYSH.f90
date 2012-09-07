!###############################################################
!     PROGRAM  : FUNX_KELDYSH
!     TYPE     : Module
!     PURPOSE  : Contains global variables
!     AUTHORS  : Adriano Amaricci
!###############################################################
module IPT_KELDYSH
  use IPT_VARS_GLOBAL
  private

  integer,save                 :: loop=1 
  type(keldysh_equilibrium_gf) :: fg0,sigma
  real(8),allocatable          :: t(:),wr(:)
  public                       :: solve_ipt_keldysh

contains

  !+-------------------------------------------------------------------+
  !PROGRAM  : 
  !PURPOSE  : 
  !COMMENT  : 
  !+-------------------------------------------------------------------+
  function solve_ipt_keldysh(fg0_,wr_,t_) result(sigma_)
    complex(8),dimension(2*L) :: fg0_,sigma_
    real(8),dimension(2*L)   :: wr_
    real(8),dimension(-L:L)  :: t_
    if(loop==1)then
       if(.not.fg0%status)  call allocate_gf(fg0,L)
       if(.not.sigma%status)call allocate_gf(sigma,L)
    endif
    fg0%ret%w=fg0_ ; t=t_ ; wr=wr_
    dt=abs(t(2)-t(1)) ; fmesh=abs(wr(2)-wr(1))
    call simpurity
    sigma_=sigma%ret%w
    loop=loop+1
  end function solve_ipt_keldysh
  !*******************************************************************
  !*******************************************************************
  !*******************************************************************







  !+-------------------------------------------------------------------+
  !PROGRAM  : 
  !PURPOSE  : 
  !COMMENT  : 
  !+-------------------------------------------------------------------+
  subroutine simpurity
    !FFT to real time:
    !G0^{<,>}(t) = FT[G0^{<,>}(w)]
    call fftgf_rw2rt(fg0%less%w,fg0%less%t,L);   fg0%less%t=xi*fmesh*fg0%less%t
    call fftgf_rw2rt(fg0%gtr%w ,fg0%gtr%t ,L);   fg0%gtr%t =xi*fmesh*fg0%gtr%t
    fg0%ret%t=ret_component_t(fg0%gtr%t,fg0%less%t,t,L)
    do i=-L,L
       sigma%less%t(i)=(U**2)*(fg0%less%t(i)**2)*fg0%gtr%t(-i)
       sigma%gtr%t(i) =(U**2)*(fg0%gtr%t(i)**2)*fg0%less%t(-i)
    enddo
    sigma%ret%t=ret_component_t(sigma%gtr%t,sigma%less%t,t,L)
    call fftgf_rt2rw(sigma%ret%t,sigma%ret%w,L); sigma%ret%w= dt*sigma%ret%w
    call fftgf_rt2rw(fg0%ret%t,fg0%ret%w,L)    ; fg0%ret%w  = dt*fg0%ret%w
  end subroutine simpurity
  !*******************************************************************
  !*******************************************************************
  !*******************************************************************






  ! !+-------------------------------------------------------------------+
  ! !PROGRAM  : 
  ! !TYPE     : Subroutine
  ! !PURPOSE  : 
  ! !COMMENT  : 
  ! !+-------------------------------------------------------------------+
  ! subroutine get_gloc_keldysh
  !   complex(8) :: zeta,sqroot,ifg
  !   real(8)    :: w
  !   integer    :: ik
  !   fgk%ret%w=zero
  !   do i=-L,L
  !      w = wr(i);zeta=cmplx(w,eps)-sigma%ret%w(i)
  !      ifg=zero
  !      do ik=1,Lk
  !         ifg=ifg+wt(ik)/(zeta-epsik(ik))
  !      enddo
  !      fgk%ret%w(i)=ifg
  !   enddo
  ! end subroutine get_gloc_keldysh
  ! !*******************************************************************
  ! !*******************************************************************
  ! !*******************************************************************


  ! !+-------------------------------------------------------------------+
  ! !PROGRAM  : 
  ! !TYPE     : Subroutine
  ! !PURPOSE  : 
  ! !COMMENT  : 
  ! !+-------------------------------------------------------------------+
  ! subroutine update_g0_keldysh(gf_)
  !   complex(8),dimension(-L:L),optional :: gf_
  !   complex(8),dimension(-L:L)          :: gf
  !   !Update Keldysh components of the Weiss Field:
  !   if(present(gf_))then
  !      gf=gf_
  !   else
  !      gf=one/(one/fgk%ret%w + sigma%ret%w)
  !   endif
  !   !Update Keldysh components of the Weiss Field:
  !   fgk0%less%w=less_component_w(gf,wr,beta)
  !   fgk0%gtr%w =gtr_component_w(gf,wr,beta)
  !   return
  ! end subroutine update_g0_keldysh
  ! !*******************************************************************
  ! !*******************************************************************
  ! !*******************************************************************







  ! !+-------------------------------------------------------------------+
  ! !PROGRAM  : 
  ! !TYPE     : Subroutine
  ! !PURPOSE  : 
  ! !COMMENT  : 
  ! !+-------------------------------------------------------------------+
  ! subroutine print_out()
  !   real(8) :: nimp
  !   fg%ret%w = one/(one/fg0%ret%w - sigma%ret%w)
  !   call getGmats(wr,fg%ret%w,gf%iw,beta)
  !   call fftgf_iw2tau(gf%iw,gf%tau,beta)
  !   nimp=-2.d0*gf%tau(L)
  !   call splot("nVSiloop.ipt"//trim(adjustl(trim(label))),iloop,nimp,append=TT)
  !   call splot('Sret_t.ipt'//trim(adjustl(trim(label))),t,exa*sigma%ret%t,append=printf)
  !   call splot('Sgtr_t.ipt'//trim(adjustl(trim(label))),t,exa*sigma%gtr%t,append=printf)
  !   call splot('Sless_t.ipt'//trim(adjustl(trim(label))),t,exa*sigma%less%t,append=printf)
  !   call splot('G0ret_t.ipt'//trim(adjustl(trim(label))),t,exa*fg0%ret%t,append=printf)
  !   call splot('Gret_realw.ipt'//trim(adjustl(trim(label))),wr,fg%ret%w,append=printf)
  !   call splot('Sret_realw.ipt'//trim(adjustl(trim(label))),wr,sigma%ret%w,append=printf)
  !   call splot('G0ret_realw.ipt'//trim(adjustl(trim(label))),wr,fg0%ret%w,append=printf)
  !   call splot('DOS.ipt'//trim(adjustl(trim(label))),wr,-aimag(fg%ret%w)/pi,append=printf)
  !   call getGmats(wr,fg0%ret%w,gf0%iw,beta)
  !   call getGmats(wr,sigma%ret%w,sf%iw,beta)
  !   call splot('GM_iw.ipt'//trim(adjustl(trim(label))),wm,gf%iw,append=printf)
  !   call splot('G0M_iw.ipt'//trim(adjustl(trim(label))),wm,gf0%iw,append=printf)
  !   call splot('SigmaM_iw.ipt'//trim(adjustl(trim(label))),wm,sf%iw,append=printf)
  ! end subroutine print_out
  ! !*******************************************************************
  ! !*******************************************************************
  ! !*******************************************************************




  ! !+-------------------------------------------------------------------+
  ! !PROGRAM  : 
  ! !TYPE     : Subroutine
  ! !PURPOSE  : 
  ! !COMMENT  : 
  ! !+-------------------------------------------------------------------+
  ! subroutine get_energy_keldysh(logic,totE)
  !   real(8),optional   :: totE
  !   logical            :: logic
  !   real(8)            :: nimp,docc,Ekin,Epot,nk(Lk)
  !   integer            :: Lm
  !   complex(8),allocatable :: gf(:),sf(:),gff(:)
  !   real(8),allocatable    :: gtau(:)

  !   Lm=int(L*beta/pi);if(beta<1)Lm=L
  !   allocate(gf(Lm),sf(Lm),gtau(0:L),gff(L))

  !   call getGmats(wr,sigma%ret%w,sf,beta)
  !   call getGmats(wr,fg%ret%w,gf,beta)

  !   call fftgf_iw2tau(gf,gtau,beta)
  !   nimp=-2.d0*gtau(L)

  !   !Energy && k-dispersed quantities
  !   Ekin=0.0
  !   do ik=1,Lk
  !      gff=zero
  !      do i=1,L
  !         w=pi/beta*dble(2*i-1)
  !         gff(i)=one/(xi*w - sorted_epsik(ik) - sf(i))
  !      enddo
  !      call fftgf_iw2tau(gff,gtau,beta)
  !      nk(ik)=-gtau(L) 
  !      Ekin=Ekin + wt(ik)*nk(ik)*epsik(ik)
  !   enddo

  !   Epot=dot_product(conjg(sf(:)),gf(:))/beta/2.0
  !   docc=0.5*nimp  - 0.25d0
  !   if(u/=0.0)docc = Epot/U + 0.5*nimp - 0.25d0
  !   if(present(totE))totE=Ekin+Epot
  !   if(logic)then
  !      call splot("nkVSepsik.ipt",sorted_epsik(1:Lk),nk(1:Lk))
  !      call splot("EtotVS"//trim(extension),xout,Ekin+Epot,append=TT)
  !      call splot("EkinVS"//trim(extension),xout,Ekin,append=TT)
  !      call splot("EpotVS"//trim(extension),xout,Epot,append=TT)
  !      call splot("doccVS"//trim(extension),xout,docc,append=TT)
  !   end if
  !   deallocate(gf,sf,gtau,gff)
  ! end subroutine get_energy_keldysh
  ! !*******************************************************************
  ! !*******************************************************************
  ! !*******************************************************************




end module IPT_KELDYSH

