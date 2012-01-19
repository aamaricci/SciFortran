module FFTGF
  use FGSL
  implicit none 
  private
  public :: fftgf_rw2rt  , fftgf_rt2rw
  public :: fftgf_iw2tau , fftgf_tau2iw
  !public :: fft_iw2tau   , fft_tau2iw
  public :: swap_fftrt2rw
  integer                          :: status
  type(fgsl_fft_complex_wavetable) :: wavetable
  type(fgsl_fft_complex_workspace) :: work

contains

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
    complex(8),dimension(-M:M)             :: dummy_out
    dummy_in = func_in
    wavetable = fgsl_fft_complex_wavetable_alloc(2*M*1_fgsl_size_t)
    work      = fgsl_fft_complex_workspace_alloc(2*M*1_fgsl_size_t)
    status    = fgsl_fft_complex_forward(dummy_in,1,2*M*1_fgsl_size_t,wavetable,work)
    call fgsl_fft_complex_workspace_free(work)
    call fgsl_fft_complex_wavetable_free(wavetable) 
    forall(i=1:2*M)dummy_out(i-M-1)=dummy_in(i) ![1,2*M]---> [-M,M-1]
    forall(i=-M:-1)func_out(i+M)=dummy_out(i)   !g[0,M-1]<--- x[-M,-1]
    forall(i=0:M-1)func_out(i-M)=dummy_out(i)   !g[-L,-1]<--- x[0,L-1]
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
    wavetable = fgsl_fft_complex_wavetable_alloc(2*M*1_fgsl_size_t)
    work      = fgsl_fft_complex_workspace_alloc(2*M*1_fgsl_size_t)
    status    = fgsl_fft_complex_backward(dummy_out,1,2*M*1_fgsl_size_t,wavetable,work)
    call fgsl_fft_complex_workspace_free(work)
    call fgsl_fft_complex_wavetable_free(wavetable) 
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
    real(8)                             :: pi,wmax,beta,mues,tau,dtau,At,w
    notail_=.false.;if(present(notail))notail_=notail
    pi=acos(-1.d0)
    n=size(gw)    ; L=size(gt)-1 ; dtau=beta/dble(L) 
    allocate(tmpGw(2*L),tmpGt(-L:L))
    wmax=pi/beta*dble(2*L-1) ; mues=-real(gw(L))*wmax**2
    tmpGw=(0.d0,0.d0)
    select case(notail_)
    case default
       do i=1,L
          w=pi/beta*dble(2*i-1)
          tail=-(mues+w*(0.d0,1.d0))/(mues**2+w**2)
          if(i<=n)tmpGw(2*i)= gw(i)-tail
          if(i>n)tmpGw(2*i)=tail
       enddo
       call fftgf_rw2rt(tmpGw,tmpGt,L) ;tmpGt=2.d0*tmpGt/beta
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
          gt(i)=real(tmpGt(i))+At
       enddo
       gt(L)=-(gt(0)+1.d0)
    case(.true.)
       if(L>n)call abort("error in fftgf_iw2tau: call w/ notail and L>n")
       call fft_iw2tau(gw(1:L),gt(0:L),beta,L)
    end select
  end subroutine fftgf_iw2tau



  !*******************************************************************
  !*******************************************************************
  !*******************************************************************



  !+-------------------------------------------------------------------+
  !PURPOSE  :  
  !+-------------------------------------------------------------------+
  subroutine fftgf_tau2iw(gt,gw,beta,notail)
    implicit none
    integer                :: i,L,n,M
    logical,optional       :: notail
    logical                :: notail_
    real(8)                :: gt(0:)
    complex(8)             :: gw(:),one=(1.d0,0.d0)
    real(8)                :: beta
    complex(8),allocatable :: xgw(:),Igw(:)
    real(8),allocatable    :: xgt(:),Igt(:)
    notail_=.false. ; if(present(notail))notail_=notail
    L=size(gt)-1    ; n=size(gw)
    select case(notail_)
    case default
       allocate(xgt(-L:L)) ; xgt(0:L)=gt(0:L) ; forall(i=1:L)xgt(-i)=-xgt(L-i)
       !Fit to get rid of the 2*i problem
       M=4*L                        !long enough to get back gw(1:n)
       allocate(Igt(-M:M),Igw(2*M)) !allocate the long arrays 
       call interp(xgt,Igt,L,M)     !w/ interpolation to treat long tails
       call fftgf_rt2rw(one*Igt,Igw,M)
       Igw=Igw*beta/dble(M)/2.d0
       forall(i=1:n)gw(i)=Igw(2*i)
       deallocate(xgt,Igt,Igw)
    case(.true.)
       call fft_tau2iw(gt(0:L),gw(1:L),beta,L)
    end select
  contains
    !This is taken from SPLINE to make this module independent    
    include "splinefft.f90"
  end subroutine fftgf_tau2iw


  !*******************************************************************
  !*******************************************************************
  !*******************************************************************




  !+-------------------------------------------------------------------+
  !PURPOSE  : sometime a swap of the two halfs of the output array 
  ! from cfft_rt2rw is needed. THis is implemented in this routine.
  !COMMENTS : func_in should have dimension (2*N), swapped around N
  !+-------------------------------------------------------------------+
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



  !*******************************************************************
  !*******************************************************************
  !*******************************************************************



  include "fftgf_no_correction.f90"




end module FFTGF
