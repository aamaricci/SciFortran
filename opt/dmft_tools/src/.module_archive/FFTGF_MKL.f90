include "mkl_trig_transforms.f90"
module FFTGF_MKL
  use MKL_TRIG_TRANSFORMS
  USE CONSTANTS, only: pi
  USE ARRAYS, only: linspace
  USE INTERPOLATE, only: linear_spline,cubic_spline
  implicit none 
  private

  integer                        :: status
  type(DFTI_DESCRIPTOR), pointer :: Handle


  !GREEN'S FUNCTION FFT:
  !- - - - - - - - - - - - - - - - - - - - - - - - - - - - -!
  !Real freq.      <==> real time
  public :: fftgf_rw2rt
  public :: fftgf_rt2rw
  !- - - - - - - - - - - - - - - - - - - - - - - - - - - - -!
  !Matsubara freq. <==> imaginary time
  public :: fftgf_iw2tau
  public :: fftgf_tau2iw
  !- - - - - - - - - - - - - - - - - - - - - - - - - - - - -!
  public :: fftff_iw2tau , fftff_tau2iw
  ! public :: fftff_iw2tau_


contains




  !*********************************************************************
  !               FOURIER TRANSFORM OF GREEN'S FUNCTIONS:
  !*********************************************************************  
  !+-------------------------------------------------------------------+
  !PURPOSE  : Evaluate the FFT of a given Green's function from 
  !real frequencies to real time.
  !+-------------------------------------------------------------------+
  subroutine fftgf_rw2rt(func_in,func_out)
    complex(8),dimension(:)                         :: func_in
    complex(8),dimension(size(func_in)),optional    :: func_out
    if(present(func_out))then
       call tfft(func_in,func_out)
    else
       call tfft(func_in)
    endif
  end subroutine fftgf_rw2rt




  !+-------------------------------------------------------------------+
  !PURPOSE  : Evaluate the FFT of a given Green's function from 
  !real time to real frequency.
  !+-------------------------------------------------------------------+
  subroutine fftgf_rt2rw(func_in,func_out)
    complex(8),dimension(:)                      :: func_in
    complex(8),dimension(size(func_in)),optional :: func_out
    if(present(func_out))then
       call itfft(func_in,func_out)
    else
       call itfft(func_in)
    endif
  end subroutine fftgf_rt2rw




  !+-------------------------------------------------------------------+
  !PURPOSE  : Evaluate the FFT of a given function from Matsubara frequencies
  ! to imaginary time. 
  !COMMENT  : The routine implements all the necessary manipulations required 
  !by the FFT in this formalism: tail subtraction and reshaping. 
  !Output is only for tau\in[0,beta]
  !if g(-tau) is required this has to be implemented in the calling code
  !using the transformation: g(-tau)=-g(beta-tau) for tau>0
  !+-------------------------------------------------------------------+
  subroutine fftgf_iw2tau(gw,gt,beta,notail,nofix)
    complex(8),dimension(:)     :: gw
    real(8),dimension(size(gw)) :: gt
    real(8)                     :: beta
    logical,optional            :: notail,nofix
    logical                     :: notail_,nofix_
    integer                     :: Liw,Ltau
    notail_=.false. ; if(present(notail))notail_=notail
    nofix_ =.false. ; if(present(nofix))nofix_=nofix
    Liw=size(gw)
    Ltau=size(gt)
    call fftgf_iw2tau_(gw,Liw,gt,Ltau,beta,notail_,nofix_)
  end subroutine fftgf_iw2tau
  !
  subroutine fftgf_iw2tau_(gw,N,gt,L,beta,notail_,nofix_)    
    complex(8),dimension(N)       :: gw
    integer                       :: N
    real(8),dimension(L)          :: gt
    integer                       :: L
    real(8)                       :: beta
    logical                       :: notail_,nofix_
    complex(8),dimension(2*(L-1)) :: tmpGw
    real(8),dimension(2*(L-1))    :: tmpGt
    integer                       :: i
    real(8)                       :: wmax,mues,fmom
    complex(8)                    :: tail
    real(8)                       :: tau,dtau,At,w
    if(N<L)stop "fftgf_iw2tau: Liw must be > Ltau"
    dtau=beta/dble(L+1)
    wmax=pi/beta*dble(2*L-1)
    mues=-dreal(gw(L-1))*wmax**2
    fmom=-dimag(gw(L-1))*wmax
    tmpGw=dcmplx(0.d0,0.d0)
    select case(notail_)
    case default
       do i=1,L-1
          w=pi/beta*dble(2*i-1)
          tail=-fmom*dcmplx(mues,w)/(mues**2+w**2)
          tmpGw(2*i) =(gw(i)-tail)
       enddo
       call fft(tmpGw)
       tmpGt = dreal(tmpGw)*2.d0/beta
       do i=1,L-1
          tau=dble(i)*dtau
          if(mues > 0.d0)then
             if((mues*beta) > 30.d0)then
                At = -exp(-mues*tau)
             else
                At = -exp(-mues*tau)/(1.d0 + exp(-beta*mues))
             endif
          else
             if((mues*beta) < -30.d0)then
                At = -exp(mues*(beta-tau))
             else
                At = -exp(-mues*tau)/(1.d0 + exp(-beta*mues))
             endif
          endif
          gt(i)=tmpGt(i)+fmom*At
       enddo
       if(.not.nofix_)gt(L)=-(gt(1)+1.d0)
    case (.true.)
       if(N/=L)stop "fftgf_iw2tau: notail requires Liw==Ltau"
       forall(i=1:L-1)tmpGw(2*i)=gw(i)
       call fft(tmpGw)
       tmpGt = dreal(tmpGw)*size(tmpGw)*2.d0/beta
       gt(1:L) = tmpGt(1:L)
       if(.not.nofix_)gt(L)=-gt(1)
    end select
  end subroutine fftgf_iw2tau_


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


  subroutine fftff_iw2tau_(gw,gt,beta)
    integer                             :: i,n,L,itau
    complex(8),dimension(:)             :: gw
    real(8),dimension(2*size(gw))       :: reF,imF,tau_long
    complex(8),dimension(0:size(gw))    :: gt
    real(8)                             :: beta,step
    real(8),dimension(0:size(gw))       :: tau_short
    !
    integer                             :: tt_type,ipar(128),ir
    real(8),dimension(3*size(gw)+2)     :: dpar
    type(dfti_descriptor), pointer      :: handle
    N=size(gw) ; L=size(gt)-1
    if(L/=N)then
       print*,"error in fftff_iw2tau: L/=N"
       stop
    endif
    N=2*L
    reF=0.d0 ; imF=0.d0
    forall(i=1:L)
       reF(2*i) = dreal(gw(i))
       imF(2*i) = dimag(gw(i))
    end forall
    tt_type=1
    call D_INIT_TRIG_TRANSFORM(N-1,tt_type,ipar,dpar,ir)
    call D_COMMIT_TRIG_TRANSFORM(reF,handle,ipar,dpar,ir)
    call D_BACKWARD_TRIG_TRANSFORM(reF,handle,ipar,dpar,ir)
    call D_COMMIT_TRIG_TRANSFORM(imF,handle,ipar,dpar,ir)
    call D_BACKWARD_TRIG_TRANSFORM(imF,handle,ipar,dpar,ir)
    call FREE_TRIG_TRANSFORM(handle,ipar,ir)
    step = beta/dble(2*L-1)
    forall(i=1:2*L)tau_long(i)=dble(i-1)*step
    step = beta/dble(L)
    forall(i=1:L+1)tau_short(i-1)=dble(i-1)*step
    call linear_spline(tau_long,dcmplx(reF,imF),tau_short,gt)
    gt(L)=-gt(0)
    gt=gt/beta*2.d0
  end subroutine fftff_iw2tau_



  !+-------------------------------------------------------------------+
  !PURPOSE  :  
  !+-------------------------------------------------------------------+
  subroutine fftgf_tau2iw(gt,gw,beta)
    real(8)    :: gt(:)
    complex(8) :: gw(:)
    integer    :: Liw,Ltau
    real(8)    :: beta
    Liw = size(gw)
    Ltau= size(gt)
    call fftgf_tau2iw_(gt,Ltau,gw,Liw,beta)
  end subroutine fftgf_tau2iw
  subroutine fftgf_tau2iw_(gt,L,gw,N,beta)
    real(8)                :: gt(L)
    integer                :: L
    complex(8)             :: gw(N),one=(1.d0,0.d0)
    integer                :: N
    real(8)                :: beta
    integer                :: i,M
    real(8),allocatable    :: xgt(:)
    complex(8),allocatable :: Igw(:)
    real(8),allocatable    :: Igt(:)
    M=32*L
    allocate(xgt(2*L),Igt(2*M),Igw(2*M))
    forall(i=1:L)
       xgt(L+i)= gt(i)
       xgt(i)  =-gt(L-i+1)                      !Valid for fermionic GF (no bosonic)
    end forall
    call interp_gtau(xgt(L+1:2*L),Igt(M+1:2*M),L,M) !interpolate to treat long tails
    forall(i=1:M)Igt(i)=-Igt(2*M-i+1)           !Valid for fermionic GF (no bosonic)
    Igw=one*Igt
    call ifft(Igw)
    Igw = -Igw*beta
    forall(i=1:L)gw(i)=Igw(2*i)
    deallocate(Xgt,Igt,Igw)
  contains
    subroutine interp_gtau(Fin,Fout,Lin,Lout)
      integer :: Lin,Lout
      real(8) :: Fin(Lin),Xin(Lin)
      real(8) :: Fout(Lout),Xout(Lout)
      Xin = linspace(0.d0,beta,Lin)
      Xout= linspace(0.d0,beta,Lout)
      call cubic_spline(Xin,Fin,Xout,Fout)
    end subroutine interp_gtau
  end subroutine fftgf_tau2iw_



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

end module FFTGF_MKL
