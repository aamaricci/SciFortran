MODULE FFTGF_FFTPACK
  USE FFT_FFTPACK
  USE CONSTANTS, only: pi,one,xi,zero
  USE ARRAYS, only: linspace,arange
  USE MATRIX, only: solve_linear_system
  USE INTERPOLATE, only: linear_spline,cubic_spline
  USE FUNCTIONS,only:fermi

  implicit none
  private


  !GREEN'S FUNCTION FFT:
  !- - - - - - - - - - - - - - - - - - - - - - - - - - - - -!
  !Real freq.      <==> real time
  public :: fftgf_rw2rt
  public :: fftgf_rt2rw
  !- - - - - - - - - - - - - - - - - - - - - - - - - - - - -!
  !Matsubara freq. <==> imaginary time
  public :: fftgf_iw2tau,f_fftgf_iw2tau
  public :: fftgf_tau2iw,f_fftgf_tau2iw
  !- - - - - - - - - - - - - - - - - - - - - - - - - - - - -!
  public :: fftff_iw2tau , fftff_tau2iw

contains


  !*********************************************************************
  !               FOURIER TRANSFORM OF GREEN'S FUNCTIONS:
  !*********************************************************************  
  !+-------------------------------------------------------------------+
  !PURPOSE  : Evaluate the FFT of a given Green's function from 
  !real frequencies to real time. 
  ! Transform in place if func_out is not given
  !+-------------------------------------------------------------------+
  subroutine fftgf_rw2rt(func_in,func_out)
    complex(8),dimension(:)                         :: func_in
    complex(8),dimension(size(func_in)),optional    :: func_out
    complex(8),dimension(size(func_in))             :: ftmp
    ftmp = func_in
    call fft(ftmp)
    call fftex(ftmp)
    if(present(func_out))then
       func_out = fftshift(ftmp)*size(ftmp)
    else
       func_in  = fftshift(ftmp)*size(ftmp)
    endif
  end subroutine fftgf_rw2rt






  !+-------------------------------------------------------------------+
  !PURPOSE  : Evaluate the FFT of a given Green's function from 
  !real time to real frequency.
  ! Transform in place if func_out is not given
  !+-------------------------------------------------------------------+
  subroutine fftgf_rt2rw(func_in,func_out)
    complex(8),dimension(:)                      :: func_in
    complex(8),dimension(size(func_in)),optional :: func_out
    complex(8),dimension(size(func_in))          :: ftmp
    ftmp = func_in
    call ifft(ftmp)
    call fftex(ftmp)
    if(present(func_out))then
       func_out = ifftshift(ftmp)
    else
       func_in  = ifftshift(ftmp)
    endif
  end subroutine fftgf_rt2rw









  !+-------------------------------------------------------------------+
  !PURPOSE  : Evaluate the FFT of a given function from Matsubara frequencies
  ! to imaginary time. 
  !Output is only for tau\in[0,beta]
  !if g(-tau) is required this has to be implemented in the calling code
  !using the transformation: g(-tau)=-g(beta-tau) for tau>0
  !+-------------------------------------------------------------------+
  subroutine fftgf_iw2tau(giw,gtau,beta,C)
    complex(8),dimension(:)      :: giw
    real(8),dimension(size(giw)) :: gtau
    real(8),dimension(0:3)       :: C
    real(8)                      :: beta
    gtau = f_fftgf_iw2tau(giw,beta,C)
  end subroutine fftgf_iw2tau

  function f_fftgf_iw2tau(giw,beta,C) result(gtau)
    complex(8),dimension(:)             :: giw
    real(8)                             :: beta
    real(8),dimension(0:3)              :: C
    real(8),dimension(size(giw))        :: gtau
    integer                             :: L
    real(8)                             :: dtau
    complex(8),dimension(:),allocatable :: Tiw,Fiw
    real(8),dimension(:),allocatable    :: Ttau,Ftau
    real(8),dimension(:),allocatable    :: tau,wm
    L=size(giw)
    allocate(wm(L),tau(L))
    wm = pi/beta*(2*arange(1,L)-1)
    tau= linspace(0d0,beta,L,mesh=dtau)
    !
    allocate(Tiw(L),Ttau(L))
    Tiw  = C(0) + C(1)/(xi*wm) + C(2)/(xi*wm)**2 + C(3)/(xi*wm)**3
    Ttau = -C(1)/2 + C(2)*(-beta+2*tau)/4       + C(3)*tau*(beta-tau)/4
    !
    allocate(Fiw(2*(L-1)),Ftau(2*(L-1)))
    Fiw=zero
    Fiw(2::2)=Giw(:)-Tiw(:)
    !
    call fft(Fiw)
    Ftau = dreal(Fiw)*2*size(Fiw)/beta
    !
    Gtau = Ftau(1:L) + Ttau
  end function f_fftgf_iw2tau

  function fft_sigma_iw2tau(giw,beta,C) result(gtau)
    complex(8),dimension(:)             :: giw
    real(8)                             :: beta
    real(8),dimension(0:1)              :: C
    real(8),dimension(size(giw))        :: gtau
    integer                             :: L
    real(8)                             :: dtau
    complex(8),dimension(:),allocatable :: Tiw,Fiw
    real(8),dimension(:),allocatable    :: Ttau,Ftau
    real(8),dimension(:),allocatable    :: tau,wm
    L=size(giw)
    allocate(wm(L),tau(L))
    wm = pi/beta*(2*arange(1,L)-1)
    tau= linspace(0d0,beta,L,mesh=dtau)
    !
    allocate(Tiw(L),Ttau(L))
    Tiw = C(0) + C(1)/(xi*wm)
    Ttau = -C(1)/2
    !
    allocate(Fiw(2*(L-1)),Ftau(2*(L-1)))
    Fiw=zero
    Fiw(2::2)=Giw(:)-Tiw(:)
    !
    call fft(Fiw)
    Ftau = dreal(Fiw)*2*size(Fiw)/beta
    !
    Gtau = Ftau(1:L) + Ttau
  end function fft_sigma_iw2tau

  ! subroutine fftgf_iw2tau(gw,gt,beta,notail,nofix)
  !   complex(8),dimension(:)     :: gw
  !   real(8),dimension(size(gw)) :: gt
  !   real(8)                     :: beta
  !   logical,optional            :: notail,nofix
  !   logical                     :: notail_,nofix_
  !   integer                     :: Liw,Ltau
  !   notail_=.false. ; if(present(notail))notail_=notail
  !   nofix_ =.false. ; if(present(nofix))nofix_=nofix
  !   Liw=size(gw)
  !   Ltau=size(gt)
  !   call fftgf_iw2tau_(gw,Liw,gt,Ltau,beta,notail_,nofix_)
  ! end subroutine fftgf_iw2tau
  ! !
  ! subroutine fftgf_iw2tau_(gw,N,gt,L,beta,notail_,nofix_)    
  !   complex(8),dimension(N)       :: gw
  !   integer                       :: N
  !   real(8),dimension(L)          :: gt
  !   integer                       :: L
  !   real(8)                       :: beta
  !   logical                       :: notail_,nofix_
  !   complex(8),dimension(2*(L-1)) :: tmpGw
  !   real(8),dimension(2*(L-1))    :: tmpGt
  !   integer                       :: i
  !   real(8)                       :: wmax,mues,fmom
  !   complex(8)                    :: tail
  !   real(8)                       :: tau,dtau,At,w
  !   if(N<L)stop "fftgf_iw2tau: Liw must be > Ltau"
  !   dtau=beta/dble(L+1)
  !   wmax=pi/beta*dble(2*L-1)
  !   mues=-dreal(gw(L-1))*wmax**2
  !   fmom=-dimag(gw(L-1))*wmax
  !   tmpGw=dcmplx(0.d0,0.d0)
  !   select case(notail_)
  !   case default
  !      do i=1,L-1
  !         w=pi/beta*dble(2*i-1)
  !         tail=-fmom*dcmplx(mues,w)/(mues**2+w**2)
  !         tmpGw(2*i) =(gw(i)-tail)
  !      enddo
  !      call fft(tmpGw)
  !      tmpGt = dreal(tmpGw)*size(tmpGw)*2.d0/beta
  !      do i=1,L-1
  !         tau=dble(i)*dtau
  !         if(mues > 0.d0)then
  !            if((mues*beta) > 30.d0)then
  !               At = -exp(-mues*tau)
  !            else
  !               At = -exp(-mues*tau)/(1.d0 + exp(-beta*mues))
  !            endif
  !         else
  !            if((mues*beta) < -30.d0)then
  !               At = -exp(mues*(beta-tau))
  !            else
  !               At = -exp(-mues*tau)/(1.d0 + exp(-beta*mues))
  !            endif
  !         endif
  !         gt(i)=tmpGt(i)+fmom*At
  !      enddo
  !      if(.not.nofix_)gt(L)=-(gt(1)+1.d0)
  !   case (.true.)
  !      if(N/=L)stop "fftgf_iw2tau: notail requires Liw==Ltau"
  !      forall(i=1:L-1)tmpGw(2*i)=gw(i)
  !      call fft(tmpGw)
  !      tmpGt = dreal(tmpGw)*size(tmpGw)*2.d0/beta
  !      gt(1:L) = tmpGt(1:L)
  !      if(.not.nofix_)gt(L)=-gt(1)
  !   end select
  ! end subroutine fftgf_iw2tau_

  subroutine fftff_iw2tau(gw,gt,beta)
    integer                             :: i,n,L,itau,M
    complex(8),dimension(:)             :: gw
    real(8),dimension(size(gw))         :: wm
    complex(8),dimension(0:)            :: gt
    real(8)                             :: beta,dtau
    n=size(gw) ; L=size(gt)-1
    if(L/=N)then
       print*,"error in fftgf_iw2tau: call w/ notail and L/=N"
       stop
    endif
    dtau=beta/dble(L) 
    gt = cmplx(0.d0,0.d0)
    forall(i=1:n)wm(i)=pi/beta*dble(2*i-1)
    forall(i=0:L)gt(i)=sum(cos(dble(i)*dtau*wm(:))*gw(:))
    gt=gt*2.d0/beta
  end subroutine fftff_iw2tau









  !+-------------------------------------------------------------------+
  !PROGRAM  : FFTGF_TAU2IW
  !+-------------------------------------------------------------------+
  subroutine fftgf_tau2iw(gtau,giw,beta,C,factor,intflag,bcflag)
    real(8),dimension(:)                :: gtau
    complex(8),dimension(size(gtau))    :: giw
    real(8)                             :: beta
    real(8),dimension(0:3)              :: C
    integer,optional                    :: factor
    integer,optional                    :: intflag
    logical,optional                    :: bcflag
    integer                             :: factor_
    integer                             :: intflag_
    logical                             :: bcflag_
    intflag_=1    ;if(present(intflag))intflag_=intflag
    factor_=50    ;if(present(factor))factor_=factor
    bcflag_=.true.;if(present(bcflag)) bcflag_ =bcflag
    giw = f_fftgf_tau2iw(gtau,beta,C,factor_,intflag_,bcflag_)
  end subroutine fftgf_tau2iw

  function f_fftgf_tau2iw(gtau,beta,C,factor,intflag,bcflag) result(giw)
    real(8),dimension(:)                :: gtau
    complex(8),dimension(size(gtau))    :: giw
    real(8)                             :: beta
    integer,optional                    :: factor
    integer,optional                    :: intflag
    logical,optional                    :: bcflag
    integer                             :: factor_
    integer                             :: intflag_
    logical                             :: bcflag_
    real(8),dimension(0:3)              :: C
    integer                             :: L,N,i
    real(8)                             :: dtau
    complex(8),dimension(:),allocatable :: Tiw,Fiw,Iiw
    real(8),dimension(:),allocatable    :: Ttau,Ftau,Itau,TtauP,TtauM
    real(8),dimension(:),allocatable    :: tau,wm,taup
    real(8)                             :: A,B,Zp,Zm
    intflag_=1    ;if(present(intflag))intflag_=intflag
    factor_=50    ;if(present(factor))factor_=factor
    bcflag_=.true.;if(present(bcflag)) bcflag_ =bcflag
    L=size(Gtau)
    allocate(wm(L),tau(L))
    wm = pi/beta*(2*arange(1,L)-1)
    tau= linspace(0d0,beta,L,mesh=dtau)
    !
    allocate(Tiw(L),Ttau(L),TtauP(L),TtauM(L))
    if(C(2)==0d0.AND.C(3)==0d0)then
       Tiw = C(0) + C(1)/(xi*wm)
       Ttau = -C(1)/2 
    else
       A = C(2)
       B = C(3)-C(2)**2
       Tiw = one/(xi*wm - A - B/(xi*wm))
       !
       Zp = 0.5d0*(A + sqrt(A**2 + 4.d0*B))
       Zm = 0.5d0*(A - sqrt(A**2 + 4.d0*B))
       if(Zp>=0.d0)then
          TtauP = Zp*exp(-Zp*tau)*fermi(-Zp,beta)
       else
          TtauP = Zp*exp(-Zp*(tau-beta))*fermi(Zp,beta)
       end if
       if(Zm>=0.d0)then
          TtauM = Zm*exp(-Zm*tau)*fermi(-Zm,beta)
       else
          TtauM = Zm*exp(-Zm*(tau-beta))*fermi(Zm,beta)
       endif
       Ttau = ( TtauM - TtauP )/(Zp-Zm)
    endif
    !
    allocate(Ftau(L),Fiw(L))
    Ftau = Gtau - Ttau
    !
    N=min(factor_*L,2**24)
    allocate(Itau(2*N),Iiw(2*N))
    select case(intflag_)
    case (1)
       call cubic_interp_gtau(Ftau,Itau(N+1:2*N),beta)
    case (2)
       call cubic_spline_gtau(Ftau,Itau(N+1:2*N),beta,bcflag_)
    end select
    forall(i=1:N)Itau(i)=-Itau(2*N-i+1)  !Valid for every fermionic GF
    Iiw = one*Itau
    call ifft(Iiw)
    Fiw = -Iiw(2::2)/size(Iiw)*beta
    !
    Giw = Fiw + Tiw
  end function f_fftgf_tau2iw

  function fft_sigma_tau2iw(gtau,beta,C,factor,intflag,bcflag) result(giw)
    real(8),dimension(:)                :: gtau
    complex(8),dimension(size(gtau))    :: giw
    real(8)                             :: beta
    integer,optional                    :: factor
    integer,optional                    :: intflag
    logical,optional                    :: bcflag
    integer                             :: factor_
    integer                             :: intflag_
    logical                             :: bcflag_
    real(8),dimension(0:1)              :: C
    integer                             :: L,N,i
    real(8)                             :: dtau
    complex(8),dimension(:),allocatable :: Tiw,Fiw,Iiw
    real(8),dimension(:),allocatable    :: Ttau,Ftau,Itau,TtauP,TtauM
    real(8),dimension(:),allocatable    :: tau,wm,taup
    real(8)                             :: A,B,Zp,Zm
    intflag_=1    ;if(present(intflag))intflag_=intflag
    factor_=20    ;if(present(factor))factor_=factor
    bcflag_=.true.;if(present(bcflag)) bcflag_ =bcflag
    L=size(Gtau)
    allocate(wm(L),tau(L))
    wm = pi/beta*(2*arange(1,L)-1)
    tau= linspace(0d0,beta,L,mesh=dtau)
    !
    allocate(Tiw(L),Ttau(L))
    Tiw = C(0) + C(1)/(xi*wm) 
    Ttau = -C(1)/2
    !
    allocate(Ftau(L),Fiw(L))
    Ftau = Gtau - Ttau
    !
    N=min(factor_*L,2**24)
    allocate(Itau(2*N),Iiw(2*N))
    select case(intflag_)
    case (1)
       call cubic_interp_gtau(Ftau,Itau(N+1:2*N),beta)
    case (2)
       call cubic_spline_gtau(Ftau,Itau(N+1:2*N),beta,bcflag_)
    end select
    forall(i=1:N)Itau(i)=-Itau(2*N-i+1)  !Valid for every fermionic GF
    Iiw = one*Itau
    call ifft(Iiw)
    Fiw = -Iiw(2::2)/size(Iiw)*beta
    !
    Giw = Fiw + Tiw
  end function fft_sigma_tau2iw

  ! subroutine fftgf_tau2iw(gt,gw,beta)
  !   real(8)    :: gt(:)
  !   complex(8) :: gw(:)
  !   integer    :: Liw,Ltau
  !   real(8)    :: beta
  !   Liw = size(gw)
  !   Ltau= size(gt)
  !   call fftgf_tau2iw_(gt,Ltau,gw,Liw,beta)
  ! end subroutine fftgf_tau2iw
  ! subroutine fftgf_tau2iw_(gt,L,gw,N,beta)
  !   real(8)                :: gt(L)
  !   integer                :: L
  !   complex(8)             :: gw(N),one=(1.d0,0.d0)
  !   integer                :: N
  !   real(8)                :: beta
  !   integer                :: i,M
  !   real(8),allocatable    :: xgt(:)
  !   complex(8),allocatable :: Igw(:)
  !   real(8),allocatable    :: Igt(:)
  !   M=32*L
  !   allocate(xgt(2*L),Igt(2*M),Igw(2*M))
  !   forall(i=1:L)
  !      xgt(L+i)= gt(i)
  !      xgt(i)  =-gt(L-i+1)                          !Valid for every fermionic GF (bosonic case not here)  
  !   end forall
  !   call interp_gtau(xgt(L+1:2*L),Igt(M+1:2*M),L,M) !interpolate to treat long tails
  !   forall(i=1:M)Igt(i)=-Igt(2*M-i+1)               !Valid for every fermionic GF (bosonic case not here)  
  !   Igw=one*Igt
  !   call ifft(Igw)
  !   Igw = -Igw/size(Igw)*beta
  !   forall(i=1:L)gw(i)=Igw(2*i)
  !   deallocate(Xgt,Igt,Igw)
  ! contains
  !   subroutine interp_gtau(Fin,Fout,Lin,Lout)
  !     integer :: Lin,Lout
  !     real(8) :: Fin(Lin),Xin(Lin)
  !     real(8) :: Fout(Lout),Xout(Lout)
  !     Xin = linspace(0.d0,beta,Lin)
  !     Xout= linspace(0.d0,beta,Lout)
  !     call cubic_spline(Xin,Fin,Xout,Fout)
  !   end subroutine interp_gtau
  ! end subroutine fftgf_tau2iw_

  subroutine fftff_tau2iw(gt,gw,beta)
    complex(8)             :: gt(0:)
    complex(8)             :: gw(:)
    real(8)                :: beta
    integer                :: i,L,n,M
    real(8),allocatable    :: reGt(:),imGt(:)
    complex(8),allocatable :: Igw(:),Igt(:)
    L=size(gt)-1    ; N=size(gw)
    M=16*L
    allocate(Igt(-M:M),Igw(2*M))
    allocate(reGt(0:M),imGt(0:M))
    call interp_gtau(dreal(gt(0:L)),reGt(0:M),L,M)
    call interp_gtau(dimag(gt(0:L)),imGt(0:M),L,M)
    Igt(0:M)=dcmplx(reGt(0:M),imGt(0:M))
    !
    forall(i=1:M)Igt(-i)=-Igt(M-i) !Valid for every fermionic GF (bosonic case not here)
    call fftgf_rt2rw(Igt,Igw)
    Igw=Igw*beta/dble(M)/2.d0
    forall(i=1:n)gw(i)=Igw(2*i)
    deallocate(Igt,Igw)
  contains
    subroutine interp_gtau(Fin,Fout,Lin,Lout)
      integer :: Lin,Lout
      real(8) :: Fin(Lin),Xin(Lin)
      real(8) :: Fout(Lout),Xout(Lout)
      Xin = linspace(0.d0,beta,Lin)
      Xout= linspace(0.d0,beta,Lout)
      call cubic_spline(Xin,Fin,Xout,Fout)
    end subroutine interp_gtau
  end subroutine fftff_tau2iw



END MODULE FFTGF_FFTPACK
