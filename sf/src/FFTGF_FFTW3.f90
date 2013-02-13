module FFTGF
  implicit none
  include "fftw3.f"
  private
  public :: cfft_1d_forward,cfft_1d_backward,cfft_1d_shift,swap_fftrt2rw,cfft_1d_ex
  public :: fftgf_rw2rt  , fftgf_rt2rw
  public :: fftgf_iw2tau , fftgf_tau2iw
  public :: fftff_iw2tau , fftff_tau2iw

  REAL(8),PARAMETER    :: PI    = 3.14159265358979323846264338327950288419716939937510D0
  integer*8            :: plan
contains

  !+-------------------------------------------------------------------+
  !PURPOSE  : Evaluate forward and backward FFT using MKL routines.
  !+-------------------------------------------------------------------+
  subroutine cfft_1d_forward(func)
    complex(8),dimension(:),intent(inout) :: func
    complex(8),dimension(size(func))      :: dummy
    call dfftw_plan_dft_1d(plan,size(func),func,dummy,FFTW_FORWARD,FFTW_ESTIMATE)
    call dfftw_execute_dft(plan,func,dummy)
    call dfftw_destroy_plan(plan)
    func=dummy
  end subroutine cfft_1d_forward

  subroutine cfft_1d_backward(func)
    complex(8),dimension(:),intent(inout) :: func
    complex(8),dimension(size(func))      :: dummy
    call dfftw_plan_dft_1d(plan,size(func),func,dummy,FFTW_BACKWARD,FFTW_ESTIMATE)
    call dfftw_execute_dft(plan,func,dummy)
    call dfftw_destroy_plan(plan)
    func=dummy
  end subroutine cfft_1d_backward

  function cfft_1d_shift(fin,L) result(fout)
    complex(8),dimension(2*L) :: fin
    complex(8),dimension(-L:L):: fout,dout
    integer                   :: i,L
    forall(i=1:2*L)dout(i-L-1)=fin(i) ![1,2*L]---> [-L,L-1]
    forall(i=-L:-1)fout(i+L)=dout(i)   !g[0,M-1]<--- x[-M,-1]
    forall(i=0:L-1)fout(i-L)=dout(i)   !g[-L,-1]<--- x[0,L-1]
  end function cfft_1d_shift


  subroutine swap_fftrt2rw(func_in)
    integer                             :: i,Nsize,Nhalf
    complex(8),dimension(:)             :: func_in
    complex(8),dimension(size(func_in)) :: dummy
    Nsize=size(func_in) ; Nhalf=Nsize/2
    dummy=func_in
    do i=1,Nhalf
       func_in(i)=dummy(Nhalf+i)
       func_in(Nhalf+i)=dummy(i)
    enddo
  end subroutine swap_fftrt2rw

  subroutine cfft_1d_ex(func)
    complex(8),dimension(:) :: func
    real(8) :: ex
    integer :: i
    ex=-1.d0
    do i=1,size(func)
       ex=-ex
       func(i)=ex*func(i)
    enddo
  end subroutine cfft_1d_ex
  !*******************************************************************
  !*******************************************************************
  !*******************************************************************


  !+-------------------------------------------------------------------+
  !PURPOSE  : Evaluate the FFT of a given Green's function from 
  !real frequencies to real time.
  !COMMENTS : INPUT: func_in (of the form 2*M)
  !           INPUT: M size of the half-array
  !           OUTPUT: func_out (-M:M)
  !+-------------------------------------------------------------------+
  subroutine fftgf_rw2rt(func_in,func_out,M)
    integer                                :: i,M
    complex(8),dimension(2*M),intent(in)   :: func_in
    complex(8),dimension(-M:M),intent(out) :: func_out
    complex(8),dimension(2*M)              :: dummy_in
    !complex(8),dimension(-M:M)             :: dummy_out
    dummy_in = func_in
    call cfft_1d_forward(dummy_in)
    func_out = cfft_1d_shift(dummy_in,M)
    func_out(M)=func_out(M-1)
    !forall(i=-M:M)func_out(i)=dummy_out(-i)
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
  subroutine fftgf_iw2tau(gw,gt,beta,notail)
    implicit none
    integer                             :: i,n,L
    logical,optional                    :: notail
    logical                             :: notail_
    complex(8),dimension(:)             :: gw
    real(8),dimension(0:)               :: gt
    complex(8),dimension(:),allocatable :: tmpGw,tmpGt
    complex(8)                          :: tail
    real(8)                             :: wmax,beta,mues,tau,dtau,At,w

    notail_=.false.;if(present(notail))notail_=notail
    !
    n=size(gw)     ; L=size(gt)-1 ; dtau=beta/real(L,8) 
    !
    allocate(tmpGw(2*L),tmpGt(-L:L))
    !
    wmax=pi/beta*dble(2*L-1)
    mues=-real(gw(L))*wmax**2
    tmpGw=(0.d0,0.d0)
    !
    select case(notail_)
    case default
       do i=1,L
          w=pi/beta*dble(2*i-1)
          tail=-(mues+w*(0.d0,1.d0))/(mues**2+w**2)
          if(i<=n)tmpGw(2*i)= gw(i)-tail
          if(i>n)tmpGw(2*i)=tail
       enddo
       call cfft_1d_forward(tmpGw)
       tmpGt = real(cfft_1d_shift(tmpGw,L),8)*2.d0/beta
       ! call fftgf_rw2rt(tmpGw,tmpGt,L)
       ! tmpGt=2.d0*tmpGt/beta
       do i=0,L-1
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
          gt(i)=tmpGt(i)+At
       enddo
       gt(L)=-(gt(0)+1.d0)
    case(.true.)
       if(L>n)then
          print*,"error in fftgf_iw2tau: call w/ notail and L>n"
          stop
       endif
       forall(i=1:L)tmpGw(2*i)  = gw(i)
       call cfft_1d_forward(tmpGw)
       tmpGt = real(cfft_1d_shift(tmpGw,L),8)*2.d0/beta
       ! call fftgf_rw2rt(tmpGw,tmpGt,L)
       ! tmpGt=2.d0*tmpGt/beta
       gt(0:L-1) = real(tmpGt(0:L-1),8)
       gt(L)=-gt(0)
    end select
  end subroutine fftgf_iw2tau

  subroutine fftff_iw2tau(gw,gt,beta)
    integer                             :: i,n,L,itau,M
    complex(8),dimension(:)             :: gw
    real(8),dimension(size(gw))         :: wm
    complex(8),dimension(0:)            :: gt
    ! real(8),dimension(:),allocatable    :: tmpGr,tmpGi
    ! complex(8),dimension(:),allocatable :: tmpGt,tmpGw,tmpGw1,tmpGw2
    real(8)                             :: wmax,beta,mues,tau,dtau,At,w
    n=size(gw) ; L=size(gt)-1
    if(L/=N)then
       print*,"error in fftgf_iw2tau: call w/ notail and L/=N"
       stop
    endif
    dtau=beta/real(L,8) 
    ! M=2*L
    ! allocate(tmpGw1(2*L),tmpGw2(2*2*L),tmpGt(-L:L))
    ! !Get full domain of the anomalous F(iw) by symmetry F(iw)=F(-iw) <= time-reversal
    ! tmpGw2=(0.d0,0.d0)
    ! tmpGw1(1:L)=gw(1:L)
    ! forall(i=1:L)tmpGw1(2*L-i+1)=gw(i)
    ! forall(i=1:M)tmpGw2(2*i)=tmpGw1(i)
    ! call cfft_1d_forward(tmpGw2)
    ! tmpGt = cfft_1d_shift(tmpGw2,2*L)/beta
    ! do i=-L,L
    !    write(200,*)real(i,8)*dtau,real(tmpGt(i)),aimag(tmpGt(i))
    ! enddo
    gt = (0.d0,0.d0)
    forall(i=1:n)wm(i)=pi/beta*real(2*i-1,8)
    forall(i=0:L)gt(i)=sum(cos(real(i,8)*dtau*wm(:))*gw(:))
    gt=gt*2.d0/beta
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
    M=32*N!L
    allocate(Igt(-M:M),Igw(2*M),cIgt(-M:M))
    call interp(gt(0:L),Igt(0:M),L,M)
    forall(i=1:M)Igt(-i)=-Igt(M-i)
    cIgt = cmplx(1.d0,0.d0)*Igt
    call fftgf_rt2rw(cIgt,Igw,M)
    Igw=Igw*beta/real(M,8)/2.d0
    forall(i=1:n)gw(i)=Igw(2*i)
    deallocate(Igt,Igw,cIgt)
  contains
    include "splinefft.f90" !This is taken from SPLINE to make this module independent    
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
    Igt(0:M)=cmplx(reGt(0:M),imGt(0:M),8)
    !
    forall(i=1:M)Igt(-i)=-Igt(M-i) !Valid for every fermionic GF (bosonic case not here)
    call fftgf_rt2rw(Igt,Igw,M)
    Igw=Igw*beta/real(M,8)/2.d0
    forall(i=1:n)gw(i)=Igw(2*i)
    deallocate(Igt,Igw)
  contains
    include "splinefft.f90" !This is taken from SPLINE to make this module independent    
  end subroutine fftff_tau2iw

  ! subroutine fftgf_tau2iw(gt,gw,beta,notail)
  !   implicit none
  !   integer                :: i,L,n,M
  !   logical,optional       :: notail
  !   logical                :: notail_
  !   real(8)                :: gt(0:)
  !   complex(8)             :: gw(:),one=(1.d0,0.d0)
  !   real(8)                :: beta
  !   complex(8),allocatable :: xgw(:),Igw(:)
  !   real(8),allocatable    :: xgt(:),Igt(:)
  !   notail_=.false. ; if(present(notail))notail_=notail
  !   L=size(gt)-1    ; n=size(gw)
  !   select case(notail_)
  !   case default
  !      allocate(xgt(-L:L)) ; xgt(0:L)=gt(0:L) ; forall(i=1:L)xgt(-i)=-xgt(L-i)
  !      !Fit to get rid of the 2*i problem
  !      M=4*L                        !long enough to get back gw(1:n)
  !      allocate(Igt(-M:M),Igw(2*M)) !allocate the long arrays 
  !      call interp(xgt,Igt,L,M)     !w/ interpolation to treat long tails
  !      call fftgf_rt2rw(one*Igt,Igw,M)
  !      Igw=Igw*beta/dble(M)/2.d0
  !      forall(i=1:n)gw(i)=Igw(2*i)
  !      deallocate(xgt,Igt,Igw)
  !   case(.true.)
  !      call fft_tau2iw(gt(0:L),gw(1:L),beta,L)
  !   end select
  ! contains
  !   !This is taken from SPLINE to make this module independent    
  !   include "splinefft.f90"
  ! end subroutine fftgf_tau2iw



  !*******************************************************************
  !*******************************************************************
  !*******************************************************************






end module FFTGF
