module FFT_FGSL
  use FGSL
  implicit none 
  private
  public :: fftgf_rw2rt  , fftgf_rt2rw
  public :: fftgf_iw2tau , fftgf_tau2iw

  public :: swap_fftrt2rw

contains

  !+-------------------------------------------------------------------+
  !PROGRAM  : CFFT_RW2RT
  !TYPE     : Subroutine
  !PURPOSE  : Evaluate the FFT of a given Green's function from 
  !real frequencies to real time.
  !COMMENTS : INPUT: func_in (of the form -M:M)
  !           INPUT: M size of the half-array
  !           OUTPUT: func_out (-M:M)
  !+-------------------------------------------------------------------+
  subroutine fftgf_rw2rt(func_in,func_out,M)
    integer                                :: i,M,N
    complex(8),dimension(-M:M),intent(in)  :: func_in
    complex(8),dimension(2*M)              :: dummy_in
    complex(8),dimension(-M:M),intent(out) :: func_out
    complex(8),dimension(-M:M)             :: dummy_out
    integer    :: status
    type(fgsl_fft_complex_wavetable) :: wavetable
    type(fgsl_fft_complex_workspace) :: work

    N=2*M
    forall(i=1:2*M)dummy_in(i) = func_in(i-M-1)
    wavetable = fgsl_fft_complex_wavetable_alloc(N*1_fgsl_size_t)
    work      = fgsl_fft_complex_workspace_alloc(N*1_fgsl_size_t)
    status    = fgsl_fft_complex_forward(dummy_in, 1, N*1_fgsl_size_t, wavetable, work)
    call fgsl_fft_complex_workspace_free(work)
    call fgsl_fft_complex_wavetable_free(wavetable) 

    !Manipulations of the out-array:
    forall(i=1:n)dummy_out(i-M-1)=dummy_in(i) ![1,2*M=N]---> [-M,M-1]
    forall(i=-M:-1)func_out(i+M)=dummy_out(i)!g[0,M-1]<--- x[-M,-1]
    forall(i=0:M-1)func_out(i-M)=dummy_out(i)!g[-L,-1]<--- x[0,L-1]
  end subroutine fftgf_rw2rt
  !*******************************************************************
  !*******************************************************************
  !*******************************************************************




  !+-------------------------------------------------------------------+
  !PROGRAM  : CFFT_RT2RW
  !TYPE     : Subroutine
  !PURPOSE  : Evaluate the FFT of a given Green's function from 
  !real time to real frequency.
  !COMMENTS : INPUT: func_in (of the form -M:M)
  !           INPUT: M size of the half-array
  !           OUTPUT: func_out (-M:M)
  !+-------------------------------------------------------------------+
  subroutine fftgf_rt2rw(func_in,func_out,M)
    integer    :: i,M,N
    real(8)    :: ex
    complex(8),dimension(-M:M),intent(in) :: func_in
    complex(8),dimension(2*M),intent(out) :: func_out
    complex(8),dimension(2*M)             :: dummy_out
    integer    :: status
    type(fgsl_fft_complex_wavetable) :: wavetable
    type(fgsl_fft_complex_workspace) :: work

    N=2*M
    forall(i=1:2*M)dummy_out(i) = func_in(i-M-1)

    wavetable = fgsl_fft_complex_wavetable_alloc(N*1_fgsl_size_t)
    work      = fgsl_fft_complex_workspace_alloc(N*1_fgsl_size_t)
    status = fgsl_fft_complex_backward(dummy_out, 1, N*1_fgsl_size_t, wavetable, work)
    call fgsl_fft_complex_workspace_free(work)
    call fgsl_fft_complex_wavetable_free(wavetable) 

    ex=-1.d0
    do i=1,n
       ex=-ex
       func_out(i)=ex*dummy_out(i)
    enddo
  end subroutine fftgf_rt2rw
  !*******************************************************************
  !*******************************************************************
  !*******************************************************************





  !+-------------------------------------------------------------------+
  !PROGRAM  : CFFT_IW2IT
  !TYPE     : Subroutine
  !PURPOSE  : Evaluate the FFT of a given function from Matsubara frequencies
  ! to imaginary time. 
  !"gw"= the GF to FFT. dimension  NMAX (fg(NMAX))
  !"gt"= output function. dimension 0:n (funct(0:N))
  !"beta"= inverse temperature (directly passed)
  !COMMENT  : The routine implements all the necessary manipulations required 
  !by the FFT in this formalism: tail subtraction and reshaping. 
  !Output is only for tau\in[0,beta]
  !if g(-tau) is required this has to be implemented in the calling code
  !using the transformation: g(-tau)=-g(beta-tau) for tau>0
  !+-------------------------------------------------------------------+
  subroutine fftgf_iw2tau(gw,gt,beta)
    implicit none
    complex(8),dimension(:)                  :: gw
    real(8),dimension(0:)                    :: gt
    complex(8),dimension(:),allocatable      :: tmpGw
    complex(8),dimension(:),allocatable      :: tmpGt
    integer :: i,n,L,nlimit
    real(8) :: pi,wmax,beta,mues,tau,dtau,At,w,zmu
    complex(8) :: tail
    n=size(gw)    ; L=size(gt)-1
    allocate(tmpGw(2*L),tmpGt(-L:L))
    pi=acos(-1.d0) ; dtau=beta/dble(L)
    wmax=pi/beta*dble(2*L-1)
    mues=-real(gw(L))*wmax**2
    tmpGw=(0.0,0.0)
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
    gt(L)=-(gt(0)+1.d0) !fix the end point:
  end subroutine fftgf_iw2tau
  !*******************************************************************
  !*******************************************************************
  !*******************************************************************




  !+-------------------------------------------------------------------+
  !PROGRAM  : CFFT_IT2IW
  !TYPE     : Subroutine
  !PURPOSE  : imaginary time. 
  !+-------------------------------------------------------------------+
  subroutine fftgf_tau2iw(gt,gw,beta)
    implicit none
    real(8)                :: gt(0:)
    complex(8)             :: gw(:),one=(1.d0,0.d0)
    real(8)                :: beta
    integer                :: i,L,n,M
    complex(8),allocatable :: xgw(:),Igw(:)
    real(8),allocatable    :: xgt(:),Igt(:)
    L=size(gt)-1  !
    n=size(gw)    !=2*L

    !eXtend G(tau) 0:L --> -L:L 
    allocate(xgt(-L:L)) 
    xgt(0:L)=gt(0:L) ; forall(i=1:L)xgt(-i)=-xgt(L-i)

    !Fit to get rid of the 2*i problem
    M=4*L !; if(M < 4*n)M=4*n    !long enough to get back gw(1:n)
    allocate(Igt(-M:M),Igw(2*M)) !allocate the long arrays 
    call interp(xgt,Igt,L,M)     !w/ interpolation to treat long tails

    !FFT to Matsubara freq.
    call fftgf_rt2rw(one*Igt,Igw,M)
    Igw=Igw*beta/dble(M)/2.d0
    forall(i=1:n)gw(i)=Igw(2*i)
    deallocate(xgt,Igt,Igw)
  contains
    include "splinefft.f90"
  end subroutine fftgf_tau2iw
  !*******************************************************************
  !*******************************************************************
  !*******************************************************************




  !+-------------------------------------------------------------------+
  !PROGRAM  : FOUR1
  !TYPE     : Subroutine
  !PURPOSE  : adapted from Num. Rec. 
  !+-------------------------------------------------------------------+
  subroutine four1(Fct,nn,isign)
    implicit none
    integer    :: nn, isign
    complex(8) :: Fct(nn)
    integer    :: i, istep, j, m, mmax, n
    real(8)    :: Fcttmp(2*nn)
    real(8)    :: tempi, tempr
    real(8)    :: theta, wi, wpi, wpr, wr, wtemp
    do i=1, nn
       Fcttmp(2*i-1) = real(Fct(i))
       Fcttmp(2*i)   = aimag(Fct(i))
    enddo
    n = 2*nn
    j = 1
    do i=1, n, 2
       if(j .gt. i)then
          tempr = Fcttmp(j)
          tempi = Fcttmp(j+1)
          Fcttmp(j)   = Fcttmp(i)
          Fcttmp(j+1) = Fcttmp(i+1)
          Fcttmp(i)   = tempr
          Fcttmp(i+1) = tempi
       endif
       m = n/2
       do while((m .ge. 2) .and. (j .gt. m))
          j = j-m
          m = m/2
       enddo
       j = j+m
    enddo
    mmax = 2
    do while(n .gt. mmax)
       istep = 2*mmax
       theta = 6.28318530717959d0/(isign*mmax)
       wpr = -2.d0*dsin(0.5d0*theta)**2
       wpi = dsin(theta)
       wr = 1.d0
       wi = 0.d0
       do m=1, mmax, 2
          do i=m, n, istep
             j = i+mmax
             tempr = sngl(wr)*Fcttmp(j)-sngl(wi)*Fcttmp(j+1)
             tempi = sngl(wr)*Fcttmp(j+1)+sngl(wi)*Fcttmp(j)
             Fcttmp(j)   = Fcttmp(i)-tempr
             Fcttmp(j+1) = Fcttmp(i+1)-tempi
             Fcttmp(i)   = Fcttmp(i)+tempr
             Fcttmp(i+1) = Fcttmp(i+1)+tempi
          enddo
          wtemp = wr
          wr = wr*wpr-wi*wpi+wr
          wi = wi*wpr+wtemp*wpi+wi
       enddo
       mmax = istep
    enddo
    do i=1, nn
       Fct(i) = cmplx(Fcttmp(2*i-1),Fcttmp(2*i))
    enddo
  end subroutine four1
  !*******************************************************************
  !*******************************************************************
  !*******************************************************************





















  !+-------------------------------------------------------------------+
  !PROGRAM  : SWAP_FFTRW2RT
  !TYPE     : Subroutine
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



end module FFT_FGSL
