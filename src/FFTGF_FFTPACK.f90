MODULE FFTGF_FFTPACK
  implicit none
  private
  public :: cfft_1d_forward,cfft_1d_backward,cfft_1d_shift,swap_fftrt2rw
  public :: fftgf_rw2rt  , fftgf_rt2rw
  public :: fftgf_iw2tau , fftgf_tau2iw

  real(8),parameter    :: pi    = 3.14159265358979323846264338327950288419716939937510d0

contains

  !+-------------------------------------------------------------------+
  !PURPOSE  : Evaluate forward and backward FFT using MKL routines.
  !+-------------------------------------------------------------------+
  subroutine cfft_1d_forward(func)
    complex(8),dimension(:),intent(inout) :: func
    real(8),dimension(5*size(func))       :: wsave
    integer                               :: L
    L=size(func)
    call zffti(L,wsave)
    call zfftf(L,func,wsave)
  end subroutine cfft_1d_forward

  subroutine cfft_1d_backward(func)
    complex(8),dimension(:),intent(inout) :: func
    real(8),dimension(5*size(func))       :: wsave
    integer                               :: L
    L=size(func)
    call zffti(L,wsave)
    call zfftb(L,func,wsave)
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
    integer                                :: M
    complex(8),dimension(2*M),intent(in)   :: func_in
    complex(8),dimension(-M:M),intent(out) :: func_out
    complex(8),dimension(2*M)              :: dummy_in
    dummy_in = func_in
    call cfft_1d_forward(dummy_in)
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
    complex(8),dimension(:)       :: gw
    real(8),dimension(0:size(gw)) :: gt
    real(8) :: beta
    integer :: L
    L=size(gw)
    call fftgf_iw2tau_(gw,gt,beta,L)
  end subroutine fftgf_iw2tau
  subroutine fftgf_iw2tau_(gw,gt,beta,L)
    implicit none
    complex(8),dimension(L)    :: gw
    real(8),dimension(0:L)     :: gt
    complex(8),dimension(2*L)  :: tmpGw
    complex(8),dimension(-L:L) :: tmpGt
    integer                    :: i,n,L
    real(8)                    :: wmax,beta,mues,tau,dtau,At,w
    complex(8)                 :: tail
    dtau=beta/dble(L)
    wmax=pi/beta*dble(2*L-1)
    mues=-dreal(gw(L))*wmax**2
    tmpGw=(0.d0,0.d0)
    do i=1,L
       w=pi/beta*real(2*i-1,8)
       tail=-cmplx(mues,w,8)/(mues**2+w**2)
       tmpGw(2*i)=gw(i)-tail
    enddo
    call cfft_1d_forward(tmpGw)
    tmpGt = cfft_1d_shift(tmpGw,L)
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
       gt(i)=dreal(tmpGt(i))*2.d0/beta+At
    enddo
    gt(L)=-(gt(0)+1.d0) !fix the end point:
  end subroutine fftgf_iw2tau_
  !*******************************************************************
  !*******************************************************************
  !*******************************************************************




  !+-------------------------------------------------------------------+
  !PROGRAM  : CFFT_IT2IW
  !TYPE     : Subroutine
  !PURPOSE  : imaginary time. 
  !+-------------------------------------------------------------------+
  subroutine fftgf_tau2iw(gt,gw,beta)
    real(8)    :: gt(0:)
    complex(8) :: gw(size(gt)-1)
    integer    :: L
    real(8)    :: beta
    L=size(gt)-1
    call fftgf_tau2iw_(gt,gw,beta,L)
  end subroutine fftgf_tau2iw
  subroutine fftgf_tau2iw_(gt,gw,beta,L)
    implicit none
    real(8)                :: gt(0:L)
    complex(8)             :: gw(L)
    real(8)                :: beta
    integer                :: i,L,n,M
    complex(8),allocatable :: Igw(:),cIgt(:)
    real(8),allocatable    :: Igt(:)
    M=32*L
    allocate(Igt(-M:M),Igw(2*M),cIgt(-M:M))
    call interp(gt(0:L),Igt(0:M),L,M)
    forall(i=1:M)Igt(-i)=-Igt(M-i) !Valid for every fermionic GF (bosonic case not here)      
    cIgt = cmplx(1.d0,0.d0)*Igt
    call fftgf_rt2rw(cIgt,Igw,M)
    Igw=Igw*beta/real(M,8)/2.d0
    forall(i=1:L)gw(i)=Igw(2*i)
    deallocate(Igt,Igw,cIgt)
  contains
    include "splinefft.f90"
  end subroutine fftgf_tau2iw_

  !*******************************************************************
  !*******************************************************************
  !*******************************************************************



END MODULE FFTGF_FFTPACK
