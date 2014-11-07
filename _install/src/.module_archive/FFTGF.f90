!include "mkl_dfti.f90"
!include "mkl_trig_transforms.f90"
module FFTGF
  !use MKL_DFTI
  !use MKL_DFT_TYPE
  use MKL_TRIG_TRANSFORMS
  !
  use FFT
  use INTERPOLATE
  !
  implicit none 
  private
  public :: fftgf_rw2rt  , fftgf_rt2rw
  public :: fftgf_iw2tau , fftgf_tau2iw
  public :: fftff_iw2tau , fftff_tau2iw
  public :: fftff_iw2tau_
  REAL(8),PARAMETER    :: PI    = 3.14159265358979323846264338327950288419716939937510D0


contains



  !+-------------------------------------------------------------------+
  !PURPOSE  : Evaluate the FFT of a given Green's function from 
  !real frequencies to real time.
  !COMMENTS : INPUT: func_in (of the form 2*M)
  !           INPUT: M size of the half-array
  !           OUTPUT: func_out (-M:M)
  !+-------------------------------------------------------------------+
  subroutine fftgf_rw2rt(func_in,func_out,M)
    integer                                :: M
    complex(8),dimension(2*M),intent(in)   :: func_in
    complex(8),dimension(-M:M),intent(out) :: func_out
    complex(8),dimension(2*M)              :: dummy_in
    dummy_in = func_in
    call cfft_1d_forward(dummy_in(1:2*M))
    func_out = cfft_1d_shift(dummy_in,M)
    func_out(M)=func_out(M-1)
  end subroutine fftgf_rw2rt



  !*******************************************************************
  !*******************************************************************
  !*******************************************************************




  !+-------------------------------------------------------------------+
  !PURPOSE  : Evaluate the FFT of a given Green's function from 
  !real time to real frequency.
  !COMMENTS : INPUT: func_in (of the form -M:M)
  !           INPUT: M size of the half-array
  !           OUTPUT: func_out (2*M)
  !+-------------------------------------------------------------------+
  subroutine fftgf_rt2rw(func_in,func_out,M)
    integer                               :: i,M
    real(8)                               :: ex
    complex(8),dimension(-M:M),intent(in) :: func_in
    complex(8),dimension(2*M),intent(out) :: func_out
    complex(8),dimension(2*M)             :: dummy_out
    forall(i=1:2*M)dummy_out(i) = func_in(i-M-1)
    call cfft_1d_backward(dummy_out)
    ex=-1.d0
    do i=1,2*M
       ex=-ex
       func_out(i)=ex*dummy_out(i)
    enddo
  end subroutine fftgf_rt2rw



  !*******************************************************************
  !*******************************************************************
  !*******************************************************************




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
    implicit none
    integer                             :: i,n,L
    logical,optional                    :: notail,nofix
    logical                             :: notail_,nofix_
    complex(8),dimension(:)             :: gw
    real(8),dimension(0:)               :: gt
    complex(8),dimension(:),allocatable :: tmpGw
    real(8),dimension(:),allocatable    :: tmpGt
    complex(8)                          :: tail
    real(8)                             :: wmax,beta,mues,tau,dtau,At,w
    notail_=.false.;if(present(notail))notail_=notail
    nofix_ =.false.;if(present(nofix))nofix_=nofix
    !
    n=size(gw)     ; L=size(gt)-1 ; dtau=beta/real(L,8) 
    !
    allocate(tmpGw(2*L),tmpGt(-L:L))
    !
    wmax = pi/beta*real(2*N-1,8)
    mues =-dreal(gw(N))*wmax**2
    tmpGw= (0.d0,0.d0)
    !
    select case(notail_)
    case default
       do i=1,L
          w=pi/beta*dble(2*i-1)
          tail=-(mues+w*(0.d0,1.d0))/(mues**2+w**2)
          ! tmpGw(2*i)= gw(i)-tail
          if(i<=n)tmpGw(2*i)= gw(i)-tail
          if(i>n)tmpGw(2*i)= tail
       enddo
       call cfft_1d_forward(tmpGw)
       tmpGt = real(cfft_1d_shift(tmpGw,L),8)*2.d0/beta
       do i=0,L-1
          tau=real(i,8)*dtau
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
          gt(i) = tmpGt(i) + At
       enddo
       if(.not.nofix_)gt(L)=-(gt(0)+1.d0)

    case(.true.)
       if(L/=N)then
          print*,"error in fftgf_iw2tau: call w/ notail and L/=N"
          stop
       endif
       forall(i=1:L)tmpGw(2*i)  = gw(i)
       call cfft_1d_forward(tmpGw)
       tmpGt = real(cfft_1d_shift(tmpGw,L),8)*2.d0/beta
       gt(0:L-1) = tmpGt(0:L-1)
       gt(L)=-gt(0)

    end select
    deallocate(tmpGw,tmpGt)
  end subroutine fftgf_iw2tau



  subroutine fftff_iw2tau_(gw,gt,beta)
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
    dtau=beta/real(L,8) 
    gt = (0.d0,0.d0)
    forall(i=1:n)wm(i)=pi/beta*real(2*i-1,8)
    forall(i=0:L)gt(i)=sum(cos(real(i,8)*dtau*wm(:))*gw(:))
    gt=gt*2.d0/beta
  end subroutine fftff_iw2tau_


  subroutine fftff_iw2tau(gw,gt,beta)
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
       reF(2*i) = real(gw(i),8)
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
  end subroutine fftff_iw2tau

  !*******************************************************************
  !*******************************************************************
  !*******************************************************************



  !+-------------------------------------------------------------------+
  !PURPOSE  :  
  !+-------------------------------------------------------------------+
  subroutine fftgf_tau2iw(gt,gw,beta)
    real(8)                :: gt(0:)
    complex(8)             :: gw(:)
    real(8)                :: beta
    integer                :: i,L,n,M
    complex(8),allocatable :: Igw(:),cIgt(:)
    real(8),allocatable    :: Igt(:)
    L=size(gt)-1    ; N=size(gw)
    M=32*N
    allocate(Igt(-M:M),Igw(2*M),cIgt(-M:M))
    call interp(gt(0:L),Igt(0:M),L,M)
    forall(i=1:M)Igt(-i)=-Igt(M-i) !Valid for every fermionic GF (bosonic case not here)      
    cIgt = cmplx(1.d0,0.d0)*Igt
    call fftgf_rt2rw(cIgt,Igw,M)
    Igw=Igw*beta/real(M,8)/2.d0
    forall(i=1:n)gw(i)=Igw(2*i)
    deallocate(Igt,Igw,cIgt)
  contains
    include "fftgf_spline.f90" !This is taken from SPLINE to make this module independent    
  end subroutine fftgf_tau2iw


  subroutine fftff_tau2iw(gt,gw,beta)
    complex(8)             :: gt(0:)
    complex(8)             :: gw(:)
    real(8)                :: beta
    integer                :: i,L,n,M
    real(8),allocatable    :: reGt(:),imGt(:)
    complex(8),allocatable :: Igw(:),Igt(:)
    L=size(gt)-1    ; N=size(gw)
    M=32*L
    allocate(Igt(-M:M),Igw(2*M))
    allocate(reGt(0:M),imGt(0:M))
    call interp(dreal(gt(0:L)),reGt(0:M),L,M)
    call interp(dimag(gt(0:L)),imGt(0:M),L,M)
    Igt(0:M)=dcmplx(reGt(0:M),imGt(0:M))
    !
    forall(i=1:M)Igt(-i)=-Igt(M-i) !Valid for every fermionic GF (bosonic case not here)
    call fftgf_rt2rw(Igt,Igw,M)
    Igw=Igw*beta/dble(M)/2.d0
    forall(i=1:n)gw(i)=Igw(2*i)
    deallocate(Igt,Igw)
  contains
    include "fftgf_spline.f90" !This is taken from SPLINE to make this module independent    
  end subroutine fftff_tau2iw



  !*******************************************************************
  !*******************************************************************
  !*******************************************************************







end module FFTGF
