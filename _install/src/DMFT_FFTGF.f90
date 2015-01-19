module DMFT_FFTGF
  USE SF_FFT_FFTPACK
  USE SF_CONSTANTS, only   : pi, one, xi, zero
  USE SF_ARRAYS, only      : linspace, arange
  USE SF_LINALG, only      : solve_linear_system
  USE SF_INTERPOLATE, only : linear_spline, cubic_spline
  USE SF_SPECIAL,only      : fermi

  implicit none
  private

  !REAL TIME OVERLOAD:
  interface fft_sigma_rw2rt
     module procedure fft_gf_rw2rt
  end interface fft_sigma_rw2rt

  interface f_fft_sigma_rw2rt
     module procedure f_fft_gf_rw2rt
  end interface f_fft_sigma_rw2rt

  interface fft_sigma_rt2rw
     module procedure fft_gf_rt2rw
  end interface fft_sigma_rt2rw

  interface f_fft_sigma_rt2rw
     module procedure f_fft_gf_rt2rw
  end interface f_fft_sigma_rt2rw

  interface fft_delta_rw2rt
     module procedure fft_gf_rw2rt
  end interface fft_delta_rw2rt

  interface f_fft_delta_rw2rt
     module procedure f_fft_gf_rw2rt
  end interface f_fft_delta_rw2rt

  interface fft_delta_rt2rw
     module procedure fft_gf_rt2rw
  end interface fft_delta_rt2rw

  interface f_fft_delta_rt2rw
     module procedure f_fft_gf_rt2rw
  end interface f_fft_delta_rt2rw



  !MATSUBARA OVERLOAD:
  interface fft_delta_iw2tau
     module procedure fft_sigma_iw2tau
  end interface fft_delta_iw2tau

  interface f_fft_delta_iw2tau
     module procedure f_fft_sigma_iw2tau
  end interface f_fft_delta_iw2tau

  interface fft_delta_tau2iw
     module procedure fft_sigma_tau2iw
  end interface fft_delta_tau2iw

  interface f_fft_delta_tau2iw
     module procedure f_fft_sigma_tau2iw
  end interface f_fft_delta_tau2iw


  !GREEN'S FUNCTION FFT:
  !- - - - - - - - - - - - - - - - - - - - - - - - - - - - -!
  !Real freq.      <==> real time
  public :: fft_gf_rw2rt
  public :: fft_sigma_rw2rt
  public :: fft_delta_rw2rt
  !
  public :: f_fft_gf_rw2rt
  public :: f_fft_sigma_rw2rt
  public :: f_fft_delta_rw2rt
  !- - - - - - - - - - - - - - - - - - - - - - - - - - - - -!
  public :: fft_gf_rt2rw
  public :: fft_sigma_rt2rw
  public :: fft_delta_rt2rw
  !
  public :: f_fft_gf_rt2rw
  public :: f_fft_sigma_rt2rw
  public :: f_fft_delta_rt2rw
  !- - - - - - - - - - - - - - - - - - - - - - - - - - - - -!
  !Matsubara freq. <==> imaginary time
  public :: fft_gf_iw2tau
  public :: fft_sigma_iw2tau
  public :: fft_delta_iw2tau
  public :: fft_ff_iw2tau
  !
  public :: f_fft_gf_iw2tau
  public :: f_fft_sigma_iw2tau
  public :: f_fft_delta_iw2tau
  !- - - - - - - - - - - - - - - - - - - - - - - - - - - - -!
  public :: fft_gf_tau2iw
  public :: fft_sigma_tau2iw
  public :: fft_delta_tau2iw
  public :: fft_ff_tau2iw
  !
  public :: f_fft_gf_tau2iw
  public :: f_fft_sigma_tau2iw
  public :: f_fft_delta_tau2iw
  !
  public :: fft_gf_extract
contains


  !*********************************************************************
  !               FOURIER TRANSFORM OF GREEN'S FUNCTIONS:
  !*********************************************************************  
  !+-------------------------------------------------------------------+
  !PURPOSE  : Evaluate the FFT of a given Green's function from 
  !real frequencies to real time. 
  ! Transform in place if func_out is not given
  !+-------------------------------------------------------------------+
  subroutine fft_gf_rw2rt(func_in,func_out)
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
  end subroutine fft_gf_rw2rt

  function f_fft_gf_rw2rt(func_in) result(func_out)
    complex(8),dimension(:)             :: func_in
    complex(8),dimension(size(func_in)) :: func_out
    complex(8),dimension(size(func_in)) :: ftmp
    ftmp = func_in
    call fft(ftmp)
    call fftex(ftmp)
    func_out = fftshift(ftmp)*size(ftmp)
  end function f_fft_gf_rw2rt






  !+-------------------------------------------------------------------+
  !PURPOSE  : Evaluate the FFT of a given Green's function from 
  !real time to real frequency.
  ! Transform in place if func_out is not given
  !+-------------------------------------------------------------------+
  subroutine fft_gf_rt2rw(func_in,func_out)
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
  end subroutine fft_gf_rt2rw

  function f_fft_gf_rt2rw(func_in) result(func_out)
    complex(8),dimension(:)             :: func_in
    complex(8),dimension(size(func_in)) :: func_out
    complex(8),dimension(size(func_in)) :: ftmp
    ftmp = func_in
    call ifft(ftmp)
    call fftex(ftmp)
    func_out = ifftshift(ftmp)
  end function f_fft_gf_rt2rw









  !+-------------------------------------------------------------------+
  !PURPOSE  : Evaluate the FFT of a Green's function from Matsubara frequencies
  ! to imaginary time.
  !Output is only for tau\in[0,beta]
  !if g(-tau) is required this has to be implemented in the calling code
  !using the transformation: g(-tau)=-g(beta-tau) for tau>0
  !+-------------------------------------------------------------------+
  subroutine fft_gf_iw2tau(giw,gtau,beta,C)
    complex(8),dimension(:)             :: giw
    real(8),dimension(:)                :: gtau
    real(8)                             :: beta,wmax,fmom,mues
    real(8),dimension(0:3),optional     :: C
    integer                             :: Liw,L
    complex(8),dimension(:),allocatable :: Tiw,Fiw
    real(8),dimension(:),allocatable    :: Ttau,Ftau
    real(8),dimension(:),allocatable    :: tau,wm
    Liw = size(giw)
    L   = size(gtau)
    if(L-1 > Liw)stop "fft_gf_iw2tau: Liw must be > Ltau-1"
    !
    allocate(wm(L),tau(L))
    wm = pi/beta*(2*arange(1,L)-1)
    tau= linspace(0d0,beta,L)
    !
    allocate(Tiw(L),Ttau(L))
    allocate(Fiw(2*(L-1)),Ftau(2*(L-1)))
    !
    if(present(C))then
       Tiw  = C(0) + C(1)/(xi*wm) + C(2)/(xi*wm)**2 + C(3)/(xi*wm)**3
       Ttau = -C(1)/2 + C(2)*(-beta+2*tau)/4       + C(3)*tau*(beta-tau)/4
    else
       wmax=pi/beta*dble(2*L-1)
       mues=-dreal(Giw(L-1))*wmax**2
       fmom=1d0
       !fmom=-dimag(gw(L-1))*wmax
       Tiw =-dcmplx(mues,wm)/(mues**2+wm**2)
       if(mues > 0.d0)then
          if((mues*beta) > 30.d0)then
             Ttau = -exp(-mues*tau)
          else
             Ttau = -exp(-mues*tau)/(1.d0 + exp(-beta*mues))
          endif
       else
          if((mues*beta) < -30.d0)then
             Ttau = -exp(mues*(beta-tau))
          else
             Ttau = -exp(-mues*tau)/(1.d0 + exp(-beta*mues))
          endif
       endif
       Ttau=fmom*Ttau
    endif
    !
    Fiw=zero
    Fiw(2::2)=Giw(:L-1)-Tiw(:)
    !
    call fft(Fiw)
    Ftau = dreal(Fiw)*2*size(Fiw)/beta
    !
    Gtau = Ftau(1:L) + Ttau
  end subroutine fft_gf_iw2tau
  !
  function f_fft_gf_iw2tau(giw,beta,C) result(gtau)
    complex(8),dimension(:)         :: giw
    real(8),dimension(size(giw))    :: gtau
    real(8)                         :: beta
    real(8),dimension(0:3),optional :: C
    if(present(C))then
       call fft_gf_iw2tau(giw,gtau,beta,C)
    else
       call fft_gf_iw2tau(giw,gtau,beta)
    endif
  end function f_fft_gf_iw2tau



  subroutine fft_sigma_iw2tau(giw,gtau,beta,C)
    complex(8),dimension(:)             :: giw
    real(8),dimension(:)                :: gtau
    real(8)                             :: beta,wmax,fmom,mues
    real(8),dimension(0:1),optional     :: C
    integer                             :: Liw,L
    complex(8),dimension(:),allocatable :: Tiw,Fiw
    real(8),dimension(:),allocatable    :: Ttau,Ftau
    real(8),dimension(:),allocatable    :: tau,wm
    Liw = size(giw)
    L   = size(gtau)
    if(L-1 > Liw)stop "fft_sigma_iw2tau: Liw must be > Ltau-1"
    !
    allocate(wm(L),tau(L))
    wm = pi/beta*(2*arange(1,L)-1)
    tau= linspace(0d0,beta,L)
    !
    allocate(Tiw(L),Ttau(L))
    allocate(Fiw(2*(L-1)),Ftau(2*(L-1)))
    !
    if(present(C))then
       Tiw = C(0) + C(1)/(xi*wm)
       Ttau = -C(1)/2
    else
       wmax=pi/beta*dble(2*L-1)
       mues=-dreal(Giw(L-1))*wmax**2
       fmom=1d0
       !fmom=-dimag(Giw(L-1))*wmax
       Tiw =-dcmplx(mues,wm)/(mues**2+wm**2)
       if(mues > 0.d0)then
          if((mues*beta) > 30.d0)then
             Ttau = -exp(-mues*tau)
          else
             Ttau = -exp(-mues*tau)/(1.d0 + exp(-beta*mues))
          endif
       else
          if((mues*beta) < -30.d0)then
             Ttau = -exp(mues*(beta-tau))
          else
             Ttau = -exp(-mues*tau)/(1.d0 + exp(-beta*mues))
          endif
       endif
       Ttau=fmom*Ttau
    endif
    !
    Fiw=zero
    Fiw(2::2)=Giw(:L-1)-Tiw(:)
    !
    call fft(Fiw)
    Ftau = dreal(Fiw)*2*size(Fiw)/beta
    !
    Gtau = Ftau(1:L) + Ttau
  end subroutine fft_sigma_iw2tau
  !
  function f_fft_sigma_iw2tau(giw,beta,C) result(gtau)
    complex(8),dimension(:)         :: giw
    real(8),dimension(size(giw))    :: gtau
    real(8)                         :: beta
    real(8),dimension(0:1),optional :: C
    if(present(C))then
       call fft_sigma_iw2tau(giw,gtau,beta,C)
    else
       call fft_sigma_iw2tau(giw,gtau,beta)
    endif
  end function f_fft_sigma_iw2tau




  !+-------------------------------------------------------------------+
  !PURPOSE  : Evaluate the FFT of an anomalous Green's function from
  ! Matsubara frequencies to complex anomalous imaginary time . 
  !+-------------------------------------------------------------------+
  subroutine fft_ff_iw2tau(gw,gt,beta)
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
  end subroutine fft_ff_iw2tau

















  !+-------------------------------------------------------------------+
  !PROGRAM  : Evaluate the FFT of a Green's function from imaginary time
  ! to Matsubara frequencies.
  !+-------------------------------------------------------------------+
  subroutine fft_gf_tau2iw(giw,gtau,beta,C,factor,intflag,bcflag)
    complex(8),dimension(:)             :: giw
    real(8),dimension(:)                :: gtau
    real(8)                             :: beta
    real(8),dimension(0:3),optional     :: C
    integer,optional                    :: factor
    integer                             :: factor_
    integer,optional                    :: intflag
    integer                             :: intflag_
    logical,optional                    :: bcflag
    logical                             :: bcflag_
    integer                             :: Liw,L,N,i
    real(8)                             :: dtau
    complex(8),dimension(:),allocatable :: Tiw,Fiw,Iiw
    real(8),dimension(:),allocatable    :: Ttau,Ftau,Itau,TtauP,TtauM
    real(8),dimension(:),allocatable    :: tau,wm,taup
    real(8)                             :: A,B,Zp,Zm
    intflag_= 1     ;if(present(intflag))intflag_=intflag
    factor_ = 50    ;if(present(factor))factor_=factor
    bcflag_ = .true.;if(present(bcflag)) bcflag_ =bcflag
    Liw     = size(giw)
    L       = size(gtau)
    if(L-1 > Liw)stop "fft_gf_tau2iw: Liw must be > Ltau-1"
    allocate(wm(Liw),tau(L))
    wm = pi/beta*(2*arange(1,Liw)-1)
    tau= linspace(0d0,beta,L)
    allocate(Tiw(Liw),Ttau(L),TtauP(L),TtauM(L))
    Tiw  = zero
    Ttau = 0d0
    if(present(C))then
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
    endif
    !
    allocate(Fiw(Liw),Ftau(L))
    Ftau = Gtau - Ttau
    !
    N=max(min(factor_*L,2**24),4*Liw)
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
  end subroutine fft_gf_tau2iw
  !
  function f_fft_gf_tau2iw(gtau,beta,C,factor,intflag,bcflag) result(giw)
    real(8),dimension(:)                :: gtau
    complex(8),dimension(size(gtau))    :: giw
    real(8)                             :: beta
    real(8),dimension(0:3),optional     :: C
    integer,optional                    :: factor
    integer,optional                    :: intflag
    logical,optional                    :: bcflag
    integer                             :: factor_
    integer                             :: intflag_
    logical                             :: bcflag_
    intflag_=1    ;if(present(intflag))intflag_=intflag
    factor_=50    ;if(present(factor))factor_=factor
    bcflag_=.true.;if(present(bcflag)) bcflag_ =bcflag
    if(present(C))then
       call fft_gf_tau2iw(giw,gtau,beta,C,factor=factor_,intflag=intflag_,bcflag=bcflag_)
    else
       call fft_gf_tau2iw(giw,gtau,beta,factor=factor_,intflag=intflag_,bcflag=bcflag_)
    endif
  end function f_fft_gf_tau2iw



  subroutine fft_sigma_tau2iw(giw,gtau,beta,C,factor,intflag,bcflag)
    complex(8),dimension(:)             :: giw
    real(8),dimension(:)                :: gtau
    real(8)                             :: beta
    integer,optional                    :: factor
    integer                             :: factor_
    integer,optional                    :: intflag
    integer                             :: intflag_
    logical,optional                    :: bcflag
    logical                             :: bcflag_
    real(8),dimension(0:1),optional     :: C
    integer                             :: Liw,L,N,i
    complex(8),dimension(:),allocatable :: Tiw,Fiw,Iiw
    real(8),dimension(:),allocatable    :: Ttau,Ftau,Itau,TtauP,TtauM
    real(8),dimension(:),allocatable    :: tau,wm,taup
    real(8)                             :: A,B,Zp,Zm
    intflag_= 1     ;if(present(intflag))intflag_=intflag
    factor_ = 20    ;if(present(factor))factor_=factor
    bcflag_ =.true. ;if(present(bcflag)) bcflag_ =bcflag
    Liw     = size(giw)
    L       = size(gtau)
    if(L-1 > Liw)stop "fft_gf_tau2iw: Liw must be > Ltau-1"
    allocate(wm(Liw),tau(L))
    wm = pi/beta*(2*arange(1,Liw)-1)
    tau= linspace(0d0,beta,L)
    !
    allocate(Tiw(Liw),Ttau(L))
    Tiw = zero
    Ttau= 0d0
    if(present(C))then
       Tiw = C(0) + C(1)/(xi*wm) 
       Ttau = -C(1)/2
    endif
    allocate(Fiw(Liw),Ftau(L))
    Ftau = Gtau - Ttau
    !
    N=max(min(factor_*L,2**24),4*Liw)
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
  end subroutine fft_sigma_tau2iw
  !
  function f_fft_sigma_tau2iw(gtau,beta,C,factor,intflag,bcflag) result(giw)
    real(8),dimension(:)                :: gtau
    complex(8),dimension(size(gtau))    :: giw
    real(8)                             :: beta
    real(8),dimension(0:1),optional     :: C
    integer,optional                    :: factor
    integer,optional                    :: intflag
    logical,optional                    :: bcflag
    integer                             :: factor_
    integer                             :: intflag_
    logical                             :: bcflag_
    intflag_=1    ;if(present(intflag))intflag_=intflag
    factor_=50    ;if(present(factor))factor_=factor
    bcflag_=.true.;if(present(bcflag)) bcflag_ =bcflag
    if(present(C))then
       call fft_sigma_tau2iw(giw,gtau,beta,C,factor=factor_,intflag=intflag_,bcflag=bcflag_)
    else
       call fft_sigma_tau2iw(giw,gtau,beta,factor=factor_,intflag=intflag_,bcflag=bcflag_)
    endif
  end function f_fft_sigma_tau2iw





  !+-------------------------------------------------------------------+
  !PROGRAM  : Evaluate the FFT of a complex anomalous Green's function from
  ! imaginary time to Matsubara frequencies.
  !+-------------------------------------------------------------------+
  subroutine fft_ff_tau2iw(gt,gw,beta)
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
    call fft_gf_rt2rw(Igt,Igw)
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
  end subroutine fft_ff_tau2iw







  ! !+-------------------------------------------------------------------+
  ! !PURPOSE  : Sample a given function G(tau) over Nfak < N points.
  ! this has been superseded by the new FFT with different number of points
  ! !+-------------------------------------------------------------------+
  subroutine fft_gf_extract(g0,g00)
    real(8),dimension(0:)  :: g0 !0:L-1
    real(8),dimension(0:) :: g00 !0:Lfak-1
    integer :: N,Nfak
    integer :: i,ip
    real(8) :: p,mismatch
    N=size(g0)-1
    Nfak=size(g00)-1
    g00(0)=g0(0)
    ! if(g0(0) > 0.d0) g00(Nfak)=1.d0-g0(0)
    ! if(g0(0) < 0.d0) g00(Nfak)=-(g0(0)+1.d0)
    mismatch=dble(N)/dble(Nfak)
    do i=1,Nfak-1
       p=dble(i)*mismatch
       ip=int(p)
       g00(i)=g0(ip)
    enddo
    g00(Nfak)=g0(N)
  end subroutine fft_gf_extract
  ! subroutine fft_gf_extract(g0,g00,beta)
  !   real(8),dimension(:)    :: g0
  !   real(8),dimension(:) :: g00
  !   real(8) :: tau0(size(g0)),tau00(size(g00))   
  !   integer                 :: N,Nfak
  !   integer                 :: i,ip,k,k0
  !   real(8)                 :: p,mismatch,beta
  !   N = size(g0)
  !   Nfak = size(g00)
  !   tau0 = linspace(0d0,beta,N)
  !   tau00 = linspace(0d0,beta,Nfak)
  !   g00(1)   = g0(1)
  !   do i=1,Nfak
  !      ip = locate(tau0,tau00(i))
  !      g00(i) = g0(ip)
  !   enddo
  !   g00(Nfak)   = g0(N)
  ! end subroutine fft_gf_extract
  ! subroutine fft_gf_extract(g0,N,g00,Nfak)
  !   real(8),dimension(N)    :: g0
  !   real(8),dimension(Nfak) :: g00
  !   integer                 :: N,Nfak
  !   integer                 :: i,ip
  !   real(8)                 :: p,mismatch
  !   g00(1)   = g0(1)
  !   g00(Nfak)= g0(N)
  !   ! if(g0(0) > 0.d0) g00(Nfak)=1.d0-g0(0)
  !   ! if(g0(0) < 0.d0) g00(Nfak)=-(g0(0)+1.d0)
  !   mismatch=dble(N)/dble(Nfak)
  !   do i=1,Nfak-1
  !      p=dble(i)*mismatch
  !      ip=int(p)
  !      g00(i)=g0(ip)
  !   enddo
  ! end subroutine fft_gf_extract







  !--------------------------------------------------------------
  ! INTERPOLATION ROUTINES:
  !--------------------------------------------------------------
  subroutine cubic_spline_gtau(Gin,Gout,beta,bcflag)
    real(8),dimension(:)          :: Gin
    real(8),dimension(:)          :: Gout
    real(8)                       :: beta
    logical,optional              :: bcflag
    logical                       :: bcflag_
    real(8),dimension(size(Gin))  :: Xin,DGin2
    real(8),dimension(size(Gout)) :: Xout
    integer                       :: Lin,Lout,i
    real(8)                       :: h
    bcflag_=.true.;if(present(bcflag))bcflag_=bcflag
    Lin = size(Gin)
    Lout= size(Gout)
    Xin = linspace(0.d0,beta,Lin,mesh=h)
    Xout= linspace(0.d0,beta,Lout)
    DGin2 = cubic_dy2(Lin,Gin,h,bcflag_)
    do i=1,Lout
       Gout(i) = cubic_splint(Xin,Gin,DGin2,Xout(i))
    enddo
  end subroutine cubic_spline_gtau

  function cubic_dy2(N,y,h,bcflag_) result(dy2)
    integer                    :: N
    real(8)                    :: h
    real(8),dimension(N)       :: y
    logical                    :: bcflag_
    real(8),dimension(N-2,N-2) :: A_
    real(8),dimension(N-2)     :: b_
    real(8),dimension(N)       :: dy2
    integer                    :: i,L
    if(bcflag_)then
       A_ = spline_matrix_A_sumrule(N-2)
       b_ = spline_term_B_sumrule(N,y)*6d0/h**2
    else
       A_ = spline_matrix_A_natural(N-2)
       b_ = spline_term_B_natural(N,y)*6d0/h**2
    endif
    call solve_linear_system(A_,b_)
    if(bcflag_)then
       dy2(1)     = 1d0/4*( 6d0/h**2*(y(N)-y(N-1) + y(2)-y(1)) - b_(2)+b_(N-1))
       dy2(2:N-1) = b_
       dy2(N)     = -dy2(1)
    else
       dy2(1)     = 0.d0
       dy2(2:N-1) = b_
       dy2(N)     = 0.d0
    endif
  end function cubic_dy2

  function spline_matrix_A_sumrule(N) result(A_)
    integer                :: N,i
    real(8),dimension(N,N) :: A_
    A_ = 0d0
    forall(i=1:N-1)
       A_(i,i)   = 4d0
       A_(i,i+1) = 1d0
       A_(i+1,i) = 1d0
    end forall
    A_(1,1) = 15d0/4
    A_(N,N) = 15d0/4
    A_(1,N) = 1d0/4
    A_(N,1) = 1d0/4
  end function spline_matrix_A_sumrule
  !
  function spline_matrix_A_natural(N) result(A_)
    integer                :: N,i
    real(8),dimension(N,N) :: A_
    A_ = 0d0
    forall(i=1:N-1)
       A_(i,i)   = 4d0
       A_(i,i+1) = 1d0
       A_(i+1,i) = 1d0
    end forall
  end function spline_matrix_A_natural

  function spline_term_B_sumrule(N,y) result(b_)
    integer                :: N,i
    real(8),dimension(N)   :: y
    real(8),dimension(N-2) :: b_
    b_(1)   =  5d0/4*y(1) - 9d0/4*y(2) + y(N-2) + 1d0/4*y(N-1) - 1d0/4*y(N)
    forall(i=3:N-2)b_(i-1) = y(i-1) - 2d0*y(i) + y(i+1)
    b_(N-2) = -1d0/4*y(1) + 1d0/4*y(2) + y(N-2) - 9d0/4*y(N-1) + 5d0/4*y(N)
  end function spline_term_B_sumrule
  !
  function spline_term_B_natural(N,y) result(b_)
    integer                :: N,i
    real(8),dimension(N)   :: y
    real(8),dimension(N-2) :: b_
    forall(i=2:N-1)b_(i-1) = y(i-1) - 2d0*y(i) + y(i+1)
  end function spline_term_B_natural

  function cubic_splint(xa,ya,y2a,x) result(splint)
    real(8), dimension(:)        :: xa
    real(8), dimension(size(xa)) :: ya,y2a
    real(8)                      :: x
    real(8)                      :: splint
    integer                      :: khi,klo,n,i
    real(8)                      :: a,b,c,d,h
    n=size(xa)
    klo=max(min(locate(xa,x),n-1),1)
    klo=locate(xa,x)
    if(klo<1)klo=1
    if(klo>n-1)klo=n-1
    khi=klo+1
    h=xa(khi)-xa(klo)
    if(h == 0.d0) stop 'bad xa input in splint'
    a=(xa(khi)-x)/h
    b=(x-xa(klo))/h
    c=1d0/6*(a**3 - a)*h**2
    d=1d0/6*(b**3 - b)*h**2
    splint =a*ya(klo) + b*ya(khi) + c*y2a(klo) + d*y2a(khi)
  end function cubic_splint

  function locate(xx,x)
    real(8), dimension(:), intent(in) :: xx
    real(8), intent(in) :: x
    integer :: locate
    integer :: n,jl,jm,ju
    logical :: ascnd
    n=size(xx)
    ascnd = (xx(n) >= xx(1))
    jl=0
    ju=n+1
    do
       if (ju-jl <= 1) exit
       jm=(ju+jl)/2
       if (ascnd .eqv. (x >= xx(jm))) then
          jl=jm
       else
          ju=jm
       end if
    end do
    if (x == xx(1)) then
       locate=1
    else if (x == xx(n)) then
       locate=n-1
    else
       locate=jl
    end if
  end function locate









  subroutine cubic_interp_gtau(Gin,Gout,beta)
    real(8),dimension(:)          :: Gin
    real(8),dimension(:)          :: Gout
    real(8)                       :: beta
    real(8),dimension(size(Gin))  :: Xin
    real(8),dimension(size(Gout)) :: Xout
    integer                       :: Lin,Lout
    Lin = size(Gin)
    Lout= size(Gout)
    Xin = linspace(0.d0,beta,Lin)
    Xout= linspace(0.d0,beta,Lout)
    call d_cub_interp_v(Xin,Gin,Xout,Gout)
  end subroutine cubic_interp_gtau

  subroutine d_cub_interp_v(Xin,Fin,Xout,Fout)
    integer                       :: i, j
    integer                       :: Lin, Lout
    real(8),dimension(:)          :: Xin
    real(8),dimension(size(Xin))  :: Fin
    real(8),dimension(:)          :: Xout
    real(8),dimension(size(Xout)) :: Fout
    real(8)                       :: xa(1:size(Fin)), ya(4,1:size(Fin))
    real(8)                       :: x, y
    real(8),external              :: ppvalu
    Lin = size(Fin) ; Lout= size(Fout)
    xa(:)  = Xin(:)
    ya(1,:)= Fin(:)
    call cubspl(xa,ya,Lin,0,0)
    do i=1,Lout
       x = Xout(i)
       Fout(i) = PPVALU(xa,ya,Lin-1,4,x,0)
    enddo
    if(Xin(Lin) >= Xout(Lout))then
       Fout(Lout)=Fin(Lin)
    else
       Fout(Lout) = Fout(Lout-2) + &
            (Xout(Lout)-Xout(Lout-2))/(Xout(Lout-1)-Xout(Lout-2))*(Fout(Lout-1)-Fout(Lout-2))
    endif
  end subroutine d_cub_interp_v

  subroutine cubspl ( tau, c, n, ibcbeg, ibcend )
    !***********************************************************************
    !! CUBSPL defines an interpolatory cubic spline.
    !  Discussion:
    !    A tridiagonal linear system for the unknown slopes S(I) of
    !    F at TAU(I), I=1,..., N, is generated and then solved by Gauss
    !    elimination, with S(I) ending up in C(2,I), for all I.

    !  Parameters:
    !
    !    Input, real ( kind = 8 ) TAU(N), the abscissas or X values of
    !    the data points.  The entries of TAU are assumed to be
    !    strictly increasing.
    !
    !    Input, integer ( kind = 4 ) N, the number of data points.  N is
    !    assumed to be at least 2.
    !
    !    Input/output, real ( kind = 8 ) C(4,N).
    !    On input, if IBCBEG or IBCBEG is 1 or 2, then C(2,1)
    !    or C(2,N) should have been set to the desired derivative
    !    values, as described further under IBCBEG and IBCEND.
    !    On output, C contains the polynomial coefficients of
    !    the cubic interpolating spline with interior knots
    !    TAU(2) through TAU(N-1).
    !    In the interval interval (TAU(I), TAU(I+1)), the spline
    !    F is given by
    !      F(X) = 
    !        C(1,I) + 
    !        C(2,I) * H +
    !        C(3,I) * H**2 / 2 + 
    !        C(4,I) * H**3 / 6.
    !    where H=X-TAU(I).  The routine PPVALU may be used to
    !    evaluate F or its derivatives from TAU, C, L=N-1,
    !    and K=4.
    !
    !    Input, integer ( kind = 4 ) IBCBEG, IBCEND, boundary condition indicators
    !    IBCBEG = 0 means no boundary condition at TAU(1) is given.
    !    In this case, the "not-a-knot condition" is used.  That
    !    is, the jump in the third derivative across TAU(2) is
    !    forced to zero.  Thus the first and the second cubic
    !    polynomial pieces are made to coincide.
    !    IBCBEG = 1 means the slope at TAU(1) is to equal the
    !    input value C(2,1).
    !    IBCBEG = 2 means the second derivative at TAU(1) is
    !    to equal C(2,1).
    !    IBCEND = 0, 1, or 2 has analogous meaning concerning the
    !    boundary condition at TAU(N), with the additional
    !    information taken from C(2,N).
    !
    implicit none
    integer ( kind = 4 ) n
    real ( kind = 8 ) c(4,n)
    real ( kind = 8 ) divdf1
    real ( kind = 8 ) divdf3
    real ( kind = 8 ) dtau
    real ( kind = 8 ) g
    integer ( kind = 4 ) i
    integer ( kind = 4 ) ibcbeg
    integer ( kind = 4 ) ibcend
    real ( kind = 8 ) tau(n)
    !
    !  C(3,*) and C(4,*) are used initially for temporary storage.
    !
    !  Store first differences of the TAU sequence in C(3,*).
    !
    !  Store first divided difference of data in C(4,*).
    !
    do i = 2, n
       c(3,i) = tau(i) - tau(i-1)
    end do
    do i = 2, n 
       c(4,i) = ( c(1,i) - c(1,i-1) ) / ( tau(i) - tau(i-1) )
    end do
    !
    !  Construct the first equation from the boundary condition
    !  at the left endpoint, of the form:
    !
    !    C(4,1) * S(1) + C(3,1) * S(2) = C(2,1)
    !
    !  IBCBEG = 0: Not-a-knot
    !
    if ( ibcbeg == 0 ) then
       if ( n <= 2 ) then
          c(4,1) = 1.0D+00
          c(3,1) = 1.0D+00
          c(2,1) = 2.0D+00 * c(4,2)
          go to 120
       end if
       c(4,1) = c(3,3)
       c(3,1) = c(3,2) + c(3,3)
       c(2,1) = ( ( c(3,2) + 2.0D+00 * c(3,1) ) * c(4,2) * c(3,3) &
            + c(3,2)**2 * c(4,3) ) / c(3,1)
       !
       !  IBCBEG = 1: derivative specified.
       !
    else if ( ibcbeg == 1 ) then
       c(4,1) = 1.0D+00
       c(3,1) = 0.0D+00
       if ( n == 2 ) then
          go to 120
       end if
       !
       !  Second derivative prescribed at left end.
       !
    else
       c(4,1) = 2.0D+00
       c(3,1) = 1.0D+00
       c(2,1) = 3.0D+00 * c(4,2) - c(3,2) / 2.0D+00 * c(2,1)
       if ( n == 2 ) then
          go to 120
       end if
    end if
    !
    !  If there are interior knots, generate the corresponding
    !  equations and carry out the forward pass of Gauss elimination,
    !  after which the I-th equation reads:
    !
    !    C(4,I) * S(I) + C(3,I) * S(I+1) = C(2,I).
    !
    do i = 2, n-1
       g = -c(3,i+1) / c(4,i-1)
       c(2,i) = g * c(2,i-1) + 3.0D+00 * ( c(3,i) * c(4,i+1) + c(3,i+1) * c(4,i) )
       c(4,i) = g * c(3,i-1) + 2.0D+00 * ( c(3,i) + c(3,i+1))
    end do
    !
    !  Construct the last equation from the second boundary condition, of
    !  the form
    !
    !    -G * C(4,N-1) * S(N-1) + C(4,N) * S(N) = C(2,N)
    !
    !  If slope is prescribed at right end, one can go directly to
    !  back-substitution, since the C array happens to be set up just
    !  right for it at this point.
    !
    if ( ibcend == 1 ) then
       go to 160
    end if
    if ( 1 < ibcend ) then
       go to 110
    end if
90  continue
    !
    !  Not-a-knot and 3 <= N, and either 3 < N or also not-a-knot
    !  at left end point.
    !
    if ( n /= 3 .or. ibcbeg /= 0 ) then
       g = c(3,n-1) + c(3,n)
       c(2,n) = ( ( c(3,n) + 2.0D+00 * g ) * c(4,n) * c(3,n-1) + c(3,n)**2 &
            * ( c(1,n-1) - c(1,n-2) ) / c(3,n-1) ) / g
       g = - g / c(4,n-1)
       c(4,n) = c(3,n-1)
       c(4,n) = c(4,n) + g * c(3,n-1)
       c(2,n) = ( g * c(2,n-1) + c(2,n) ) / c(4,n)
       go to 160
    end if
    !
    !  N = 3 and not-a-knot also at left.
    !
100 continue
    c(2,n) = 2.0D+00 * c(4,n)
    c(4,n) = 1.0D+00
    g = -1.0D+00 / c(4,n-1)
    c(4,n) = c(4,n) - c(3,n-1) / c(4,n-1)
    c(2,n) = ( g * c(2,n-1) + c(2,n) ) / c(4,n)
    go to 160
    !
    !  IBCEND = 2: Second derivative prescribed at right endpoint.
    !
110 continue
    c(2,n) = 3.0D+00 * c(4,n) + c(3,n) / 2.0D+00 * c(2,n)
    c(4,n) = 2.0D+00
    g = -1.0D+00 / c(4,n-1)
    c(4,n) = c(4,n) - c(3,n-1) / c(4,n-1)
    c(2,n) = ( g * c(2,n-1) + c(2,n) ) / c(4,n)
    go to 160
    !
    !  N = 2.
    !
120 continue
    if ( ibcend == 2  ) then
       c(2,n) = 3.0D+00 * c(4,n) + c(3,n) / 2.0D+00 * c(2,n)
       c(4,n) = 2.0D+00
       g = -1.0D+00 / c(4,n-1)
       c(4,n) = c(4,n) - c(3,n-1) / c(4,n-1)
       c(2,n) = ( g * c(2,n-1) + c(2,n) ) / c(4,n)
    else if ( ibcend == 0 .and. ibcbeg /= 0 ) then
       c(2,n) = 2.0D+00 * c(4,n)
       c(4,n) = 1.0D+00
       g = -1.0D+00 / c(4,n-1)
       c(4,n) = c(4,n) - c(3,n-1) / c(4,n-1)
       c(2,n) = ( g * c(2,n-1) + c(2,n) ) / c(4,n)
    else if ( ibcend == 0 .and. ibcbeg == 0 ) then
       c(2,n) = c(4,n)
    end if
    !
    !  Back solve the upper triangular system 
    !
    !    C(4,I) * S(I) + C(3,I) * S(I+1) = B(I)
    !
    !  for the slopes C(2,I), given that S(N) is already known.
    !
160 continue
    do i = n-1, 1, -1
       c(2,i) = ( c(2,i) - c(3,i) * c(2,i+1) ) / c(4,i)
    end do
    !
    !  Generate cubic coefficients in each interval, that is, the
    !  derivatives at its left endpoint, from value and slope at its
    !  endpoints.
    !
    do i = 2, n
       dtau = c(3,i)
       divdf1 = ( c(1,i) - c(1,i-1) ) / dtau
       divdf3 = c(2,i-1) + c(2,i) - 2.0D+00 * divdf1
       c(3,i-1) = 2.0D+00 * ( divdf1 - c(2,i-1) - divdf3 ) / dtau
       c(4,i-1) = 6.0D+00 * divdf3 / dtau**2
    end do
    return
  end subroutine cubspl



  function ppvalu ( break, coef, l, k, x, jderiv )
    !***********************************************************************
    !
    !! PPVALU evaluates a piecewise polynomial function or its derivative.
    !
    !  Discussion:
    !
    !    PPVALU calculates the value at X of the JDERIV-th derivative of
    !    the piecewise polynomial function F from its piecewise
    !    polynomial representation.
    !
    !    The interval index I, appropriate for X, is found through a
    !    call to INTERV.  The formula for the JDERIV-th derivative
    !    of F is then evaluated by nested multiplication.
    !
    !    The J-th derivative of F is given by:
    !
    !      (d**J) F(X) = 
    !        COEF(J+1,I) + H * (
    !        COEF(J+2,I) + H * (
    !        ...
    !        COEF(K-1,I) + H * (
    !        COEF(K,  I) / (K-J-1) ) / (K-J-2) ... ) / 2 ) / 1
    !
    !    with
    !
    !      H = X - BREAK(I)
    !
    !    and
    !
    !      I = max ( 1, max ( J, BREAK(J) <= X, 1 <= J <= L ) ).
    !
    !  Modified:
    !
    !    16 February 2007
    !
    !  Author:
    !
    !    Carl DeBoor
    !
    !  Reference:
    !
    !    Carl DeBoor,
    !    A Practical Guide to Splines,
    !    Springer, 2001,
    !    ISBN: 0387953663,
    !    LC: QA1.A647.v27.
    !
    !  Parameters:
    !
    !    Input, real ( kind = 8 ) BREAK(L+1), real COEF(*), integer L, the
    !    piecewise polynomial representation of the function F to be evaluated.
    !
    !    Input, integer ( kind = 4 ) K, the order of the polynomial pieces that 
    !    make up the function F.  
    !    The usual value for K is 4, signifying a piecewise 
    !    cubic polynomial.
    !
    !    Input, real ( kind = 8 ) X, the point at which to evaluate F or
    !    of its derivatives.
    !
    !    Input, integer ( kind = 4 ) JDERIV, the order of the derivative to be
    !    evaluated.  If JDERIV is 0, then F itself is evaluated,
    !    which is actually the most common case.  It is assumed
    !    that JDERIV is zero or positive.
    !
    !    Output, real ( kind = 8 ) PPVALU, the value of the JDERIV-th
    !    derivative of F at X.
    !
    implicit none
    integer ( kind = 4 ) k
    integer ( kind = 4 ) l
    real ( kind = 8 ) break(l+1)
    real ( kind = 8 ) coef(k,l)
    real ( kind = 8 ) fmmjdr
    real ( kind = 8 ) h
    integer ( kind = 4 ) i
    integer ( kind = 4 ) jderiv
    integer ( kind = 4 ) m
    integer ( kind = 4 ) ndummy
    real ( kind = 8 ) ppvalu
    real ( kind = 8 ) value
    real ( kind = 8 ) x
    value = 0.0D+00
    fmmjdr = k - jderiv
    !
    !  Derivatives of order K or higher are identically zero.
    !
    if ( k <= jderiv ) then
       return
    end if
    !
    !  Find the index I of the largest breakpoint to the left of X.
    !
    call interv ( break, l+1, x, i, ndummy )
    !
    !  Evaluate the JDERIV-th derivative of the I-th polynomial piece at X.
    !
    h = x - break(i)
    m = k
    do
       value = ( value / fmmjdr ) * h + coef(m,i)
       m = m - 1
       fmmjdr = fmmjdr - 1.0D+00
       if ( fmmjdr <= 0.0D+00 ) then
          exit
       end if
    end do
    ppvalu = value
    return
  end function ppvalu


  subroutine interv ( xt, lxt, x, left, mflag )
    !********************************************************************
    !
    !! INTERV brackets a real value in an ascending vector of values.
    !
    !  Discussion:
    !
    !    The XT array is a set of increasing values.  The goal of the routine
    !    is to determine the largest index I so that XT(I) <= X.
    !
    !    The routine is designed to be efficient in the common situation
    !    that it is called repeatedly, with X taken from an increasing
    !    or decreasing sequence.
    !
    !    This will happen when a piecewise polynomial is to be graphed.
    !    The first guess for LEFT is therefore taken to be the value
    !    returned at the previous call and stored in the local variable ILO.
    !
    !    A first check ascertains that ILO < LXT.  This is necessary
    !    since the present call may have nothing to do with the previous
    !    call.  Then, if 
    !
    !      XT(ILO) <= X < XT(ILO+1), 
    !
    !    we set LEFT = ILO and are done after just three comparisons.
    !
    !    Otherwise, we repeatedly double the difference ISTEP = IHI - ILO
    !    while also moving ILO and IHI in the direction of X, until
    !
    !      XT(ILO) <= X < XT(IHI)
    !
    !    after which we use bisection to get, in addition, ILO + 1 = IHI.
    !    The value LEFT = ILO is then returned.
    !
    !  Modified:
    !
    !    14 February 2007
    !
    !  Author:
    !
    !    Carl DeBoor
    !
    !  Reference:
    !
    !    Carl DeBoor,
    !    A Practical Guide to Splines,
    !    Springer, 2001,
    !    ISBN: 0387953663,
    !    LC: QA1.A647.v27.
    !
    !  Parameters:
    !
    !    Input, real ( kind = 8 ) XT(LXT), a nondecreasing sequence of values.
    !
    !    Input, integer ( kind = 4 ) LXT, the dimension of XT.
    !
    !    Input, real ( kind = 8 ) X, the point whose location with 
    !    respect to the sequence XT is to be determined.
    !
    !    Output, integer ( kind = 4 ) LEFT, the index of the bracketing value:
    !      1     if             X  <  XT(1)
    !      I     if   XT(I)  <= X  < XT(I+1)
    !      LXT   if  XT(LXT) <= X
    !
    !    Output, integer ( kind = 4 ) MFLAG, indicates whether X lies within the
    !    range of the data.
    !    -1:            X  <  XT(1)
    !     0: XT(I)   <= X  < XT(I+1)
    !    +1: XT(LXT) <= X
    !
    implicit none
    integer ( kind = 4 ) lxt
    integer ( kind = 4 ) left
    integer ( kind = 4 ) mflag
    integer ( kind = 4 ) ihi
    integer ( kind = 4 ), save :: ilo = 1
    integer ( kind = 4 ) istep
    integer ( kind = 4 ) middle
    real ( kind = 8 ) x
    real ( kind = 8 ) xt(lxt)
    ihi = ilo + 1
    if ( lxt <= ihi ) then
       if ( xt(lxt) <= x ) then
          go to 110
       end if
       if ( lxt <= 1 ) then
          mflag = -1
          left = 1
          return
       end if
       ilo = lxt - 1
       ihi = lxt
    end if
    if ( xt(ihi) <= x ) then
       go to 20
    end if
    if ( xt(ilo) <= x ) then
       mflag = 0
       left = ilo
       return
    end if
    !
    !  Now X < XT(ILO).  Decrease ILO to capture X.
    !
    istep =  1
10  continue
    ihi = ilo
    ilo = ihi - istep
    if ( 1 < ilo ) then
       if ( xt(ilo) <= x ) then
          go to 50
       end if
       istep = istep * 2
       go to 10
    end if
    ilo = 1
    if ( x < xt(1) ) then
       mflag = -1
       left = 1
       return
    end if
    go to 50
    !
    !  Now XT(IHI) <= X.  Increase IHI to capture X.
    !
20  continue
    istep = 1
30  continue
    ilo = ihi
    ihi = ilo + istep
    if ( ihi < lxt ) then
       if ( x < xt(ihi) ) then
          go to 50
       end if
       istep = istep * 2
       go to 30
    end if
    if ( xt(lxt) <= x ) then
       go to 110
    end if
    !
    !  Now XT(ILO) < = X < XT(IHI).  Narrow the interval.
    !
    ihi = lxt
50  continue
    do
       middle = ( ilo + ihi ) / 2
       if ( middle == ilo ) then
          mflag = 0
          left = ilo
          return
       end if
       !
       !  It is assumed that MIDDLE = ILO in case IHI = ILO+1.
       !
       if ( xt(middle) <= x ) then
          ilo = middle
       else
          ihi = middle
       end if
    end do
    !
    !  Set output and return.
    !
110 continue
    mflag = 1
    if ( x == xt(lxt) ) then
       mflag = 0
    end if
    do left = lxt, 1, -1
       if ( xt(left) < xt(lxt) ) then
          return
       end if
    end do
    return
  end subroutine interv



END MODULE DMFT_FFTGF
