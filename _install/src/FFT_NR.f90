module FFT_NR
  implicit none 
  private

  !- - - - - - - - - - - - - - - - - - - - - - - - - - - - -!
  !Fast Fourier Transforms of time/frequency signal
  !- - - - - - - - - - - - - - - - - - - - - - - - - - - - -!
  interface tfft
     module procedure d_tfft,c_tfft
  end interface tfft
  public :: tfft
  public :: d_tfft
  public :: c_tfft

  interface itfft
     module procedure d_itfft,c_itfft
  end interface itfft
  public :: itfft
  public :: d_itfft
  public :: c_itfft


  !- - - - - - - - - - - - - - - - - - - - - - - - - - - - -!
  !Fast Fourier Transforms
  !- - - - - - - - - - - - - - - - - - - - - - - - - - - - -!
  interface fft
     module procedure rfft_1d_forward,cfft_1d_forward
  end interface fft
  public :: fft
  public :: rfft_1d_forward
  public :: cfft_1d_forward

  interface ifft
     module procedure rfft_1d_backward,cfft_1d_backward
  end interface ifft
  public :: ifft
  public :: rfft_1d_backward
  public :: cfft_1d_backward



  !- - - - - - - - - - - - - - - - - - - - - - - - - - - - -!
  !HELPER FUNCTIONS:
  !- - - - - - - - - - - - - - - - - - - - - - - - - - - - -!
  interface fftshift
     module procedure rfft_1d_shift,cfft_1d_shift
  end interface fftshift
  public :: fftshift
  public :: rfft_1d_shift
  public :: cfft_1d_shift

  interface ifftshift
     module procedure rfft_1d_ishift,cfft_1d_ishift
  end interface ifftshift
  public :: ifftshift
  public :: rfft_1d_ishift
  public :: cfft_1d_ishift

  interface fftex
     module procedure rfft_1d_ex,cfft_1d_ex
  end interface fftex
  public :: fftex
  public :: rfft_1d_ex
  public :: cfft_1d_ex
  !- - - - - - - - - - - - - - - - - - - - - - - - - - - - -!




  !MODIFIED NUMERICAL RECIPES  FFT:
  !- - - - - - - - - - - - - - - - - - - - - - - - - - - - -!
  public :: rfour1
  public :: cfour1
  public :: cosft2
  public :: realft



contains



  !*********************************************************************
  !               TIME <==> FREQUENCY DOMAIN FFT:
  !*********************************************************************  
  !+-------------------------------------------------------------------+
  !PURPOSE  : Evaluate the FFT of a function from time to frequency. 
  !+-------------------------------------------------------------------+
  subroutine d_tfft(func_in,func_out)
    real(8),dimension(:)                         :: func_in
    real(8),dimension(size(func_in)),optional    :: func_out
    complex(8),dimension(size(func_in))             :: ftmp
    ftmp = ifftshift(func_in)  
    call fft(ftmp)
    if(present(func_out))then
       func_out = fftshift(ftmp)
    else
       func_in  = fftshift(ftmp)
    endif
  end subroutine d_tfft
  !
  subroutine c_tfft(func_in,func_out)
    complex(8),dimension(:)                         :: func_in
    complex(8),dimension(size(func_in)),optional    :: func_out
    complex(8),dimension(size(func_in))             :: ftmp
    ftmp = ifftshift(func_in)  
    call fft(ftmp)
    if(present(func_out))then
       func_out = fftshift(ftmp)
    else
       func_in  = fftshift(ftmp)
    endif
  end subroutine c_tfft

  !+-------------------------------------------------------------------+
  !PURPOSE  : Evaluate the FFT of a function from frequency to time.
  !+-------------------------------------------------------------------+
  subroutine d_itfft(func_in,func_out)
    real(8),dimension(:)                      :: func_in
    real(8),dimension(size(func_in)),optional :: func_out
    complex(8),dimension(size(func_in))          :: ftmp
    ftmp = func_in
    call ifft(ftmp)
    call fftex(ftmp)
    if(present(func_out))then
       func_out = ifftshift(ftmp)/size(ftmp)
    else
       func_in  = ifftshift(ftmp)/size(ftmp)
    endif
  end subroutine d_itfft
  !
  subroutine c_itfft(func_in,func_out)
    complex(8),dimension(:)                      :: func_in
    complex(8),dimension(size(func_in)),optional :: func_out
    complex(8),dimension(size(func_in))          :: ftmp
    ftmp = func_in
    call ifft(ftmp)
    call fftex(ftmp)
    if(present(func_out))then
       func_out = ifftshift(ftmp)/size(ftmp)
    else
       func_in  = ifftshift(ftmp)/size(ftmp)
    endif
  end subroutine c_itfft









  !*********************************************************************
  !             FAST FOURIER TRANSFORM FUNCTIONS:
  !*********************************************************************  
  !                       FORWARD TRANSFORM
  !+-------------------------------------------------------------------+
  !PURPOSE  : Evaluate forward 1-DIM FFT using MKL routines.
  !out is :
  ! [y(0), y(1), ..., y(N/2),     y(-N/2+1), ...,   y(-1)]   if N is even
  ! [y(0), y(1), ..., y((N-1)/2), y(-(N-1)/2), ..., y(-1)]   if N is odd
  !+-------------------------------------------------------------------+
  subroutine rfft_1d_forward(func)
    real(8),dimension(:),intent(inout) :: func
    complex(8),dimension(size(func)) :: ftmp
    call nr_radix2_test(size(func))
    ftmp=func
    call cfour1(ftmp,size(func),1)
    func=dreal(ftmp)
  end subroutine rfft_1d_forward
  !
  subroutine cfft_1d_forward(func)
    complex(8),dimension(:),intent(inout) :: func
    call nr_radix2_test(size(func))
    call cfour1(func,size(func),1)
  end subroutine cfft_1d_forward


  !*********************************************************************
  !                       BACKWARD TRANSFORM
  !+-------------------------------------------------------------------+
  !PURPOSE  : Evaluate 1-DIM backward FFT using MKL routines.
  ! The returned real/complex array contains y(0), y(1),..., y(n-1) 
  !+-------------------------------------------------------------------+
  subroutine rfft_1d_backward(func)
    real(8),dimension(:),intent(inout) :: func
    complex(8),dimension(size(func))   :: ftmp
    call nr_radix2_test(size(func))
    ftmp=func
    call cfour1(ftmp,size(func),-1)
    func=dreal(ftmp)
  end subroutine rfft_1d_backward
  !
  subroutine cfft_1d_backward(func)
    complex(8),dimension(:),intent(inout) :: func
    call nr_radix2_test(size(func))
    call cfour1(func,size(func),-1)
  end subroutine cfft_1d_backward









  !*********************************************************************
  !                           HELPER FUNCTIONS:
  !*********************************************************************  
  !+-------------------------------------------------------------------+
  !PURPOSE  : test RADIX 2 condition for the NR based FFT
  !+-------------------------------------------------------------------+
  subroutine nr_radix2_test(M)
    integer :: M
    if(IAND(M,M-1)/=0)stop "num. rec. FFT:  radix_2 condition violated"
  end subroutine nr_radix2_test


  !+-------------------------------------------------------------------+
  !PURPOSE  : Shift the zero-frequency to the center of the 1-DIM array
  ! output of the forward-FFT is:
  ! [y(0), y(1), ..., y(N/2),     y(-N/2+1), ...,   y(-1)]   if N is even
  ! [y(0), y(1), ..., y((N-1)/2), y(-(N-1)/2), ..., y(-1)]   if N is odd
  ! using *shift produces:
  ! [y(-N/2+1),   ..., y(-1), y(0), y(1), ..., y(N/2)]       if N is even
  ! [y(-(N-1)/2), ..., y(-1), y(0), y(1), ..., y((N-1)/2)]   if N is even
  !+-------------------------------------------------------------------+
  function rfft_1d_shift(fin) result(fout)
    real(8),dimension(:)         :: fin
    real(8),dimension(size(fin)) :: fout
    integer                      :: L,p2
    L  = size(fin)
    p2 = floor(dble(L+1)/2.d0)
    fout = [fin(p2+1:),fin(1:p2)]
  end function rfft_1d_shift
  !
  function cfft_1d_shift(fin) result(fout)
    complex(8),dimension(:)          :: fin
    complex(8),dimension(size(fin))  :: fout
    integer                          :: L,p2
    L  = size(fin)
    p2 = floor(dble(L+1)/2.d0)
    fout = [fin(p2+1:),fin(1:p2)]
  end function cfft_1d_shift


  !+-------------------------------------------------------------------+
  !PURPOSE  : The inverse of _shift
  !+-------------------------------------------------------------------+
  function rfft_1d_ishift(fin) result(fout)
    real(8),dimension(:)         :: fin
    real(8),dimension(size(fin)) :: fout
    integer                      :: L,p2
    L  = size(fin)
    p2 = L-floor(dble(L+1)/2.d0)
    fout = [fin(p2+1:),fin(1:p2)]
  end function rfft_1d_ishift
  !
  function cfft_1d_ishift(fin) result(fout)
    complex(8),dimension(:)          :: fin
    complex(8),dimension(size(fin))  :: fout
    integer                          :: L,p2
    L  = size(fin)
    p2 = L-floor(dble(L+1)/2.d0)
    fout = [fin(p2+1:),fin(1:p2)]
  end function cfft_1d_ishift


  !+-------------------------------------------------------------------+
  !PURPOSE  : Alternate sign of the output of a FFT:
  !+-------------------------------------------------------------------+
  subroutine rfft_1d_ex(func)
    real(8),dimension(:)    :: func
    real(8)                 :: ex
    integer                 :: i
    ex=-1.d0
    do i=1,size(func)
       ex=-ex
       func(i)=ex*func(i)
    enddo
  end subroutine rfft_1d_ex
  !
  subroutine cfft_1d_ex(func)
    complex(8),dimension(:) :: func
    real(8)                 :: ex
    integer                 :: i
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
  !PROGRAM  : FOUR1, COSFT2, REALFT
  !TYPE     : Subroutines
  !PURPOSE  : adapted from Num. Rec. 
  !+-------------------------------------------------------------------+
  subroutine cosft2(y,n,isign)
    integer :: i,isign,n
    real(8) :: y(n)
    real(8) ::  sum,sum1,y1,y2,ytemp
    real(8) ::  theta,wi,wi1,wpi,wpr,wr,wr1,wtemp
    theta=0.5d0*pi/n
    wr=1.0d0
    wi=0.0d0
    wr1=cos(theta)
    wi1=sin(theta)
    wpr=-2.0d0*wi1**2
    wpi=sin(2.d0*theta)
    if(isign.eq.1)then
       do i=1,n/2
          y1=0.5*(y(i)+y(n-i+1))
          y2=wi1*(y(i)-y(n-i+1))
          y(i)=y1+y2
          y(n-i+1)=y1-y2
          wtemp=wr1
          wr1=wr1*wpr-wi1*wpi+wr1
          wi1=wi1*wpr+wtemp*wpi+wi1
       enddo
       call realft(y,n,1)
       do i=3,n,2
          wtemp=wr
          wr=wr*wpr-wi*wpi+wr
          wi=wi*wpr+wtemp*wpi+wi
          y1=y(i)*wr-y(i+1)*wi
          y2=y(i+1)*wr+y(i)*wi
          y(i)=y1
          y(i+1)=y2
       enddo
       sum=0.5*y(2)
       do i=n,2,-2
          sum1=sum
          sum=sum+y(i)
          y(i)=sum1
       enddo
    else if(isign.eq.-1)then
       ytemp=y(n)
       do i=n,4,-2
          y(i)=y(i-2)-y(i)
       enddo
       y(2)=2.0*ytemp
       do i=3,n,2
          wtemp=wr
          wr=wr*wpr-wi*wpi+wr
          wi=wi*wpr+wtemp*wpi+wi
          y1=y(i)*wr+y(i+1)*wi
          y2=y(i+1)*wr-y(i)*wi
          y(i)=y1
          y(i+1)=y2
       enddo
       call realft(y,n,-1)
       do i=1,n/2
          y1=y(i)+y(n-i+1)
          y2=(0.5/wi1)*(y(i)-y(n-i+1))
          y(i)=0.5*(y1+y2)
          y(n-i+1)=0.5*(y1-y2)
          wtemp=wr1
          wr1=wr1*wpr-wi1*wpi+wr1
          wi1=wi1*wpr+wtemp*wpi+wi1
       enddo
    endif
  end subroutine cosft2

  subroutine realft(data,n,isign)
    integer :: isign,n
    real(8) :: data(n)
    integer ::  i,i1,i2,i3,i4,n2p3
    real(8) ::  c1,c2,h1i,h1r,h2i,h2r,wis,wrs
    real(8) ::  theta,wi,wpi,wpr,wr,wtemp
    theta=3.141592653589793d0/dble(n/2)
    c1=0.5
    if (isign.eq.1) then
       c2=-0.5
       call rfour1(data,n/2,+1)
    else
       c2=0.5
       theta=-theta
    endif
    wpr=-2.0d0*sin(0.5d0*theta)**2
    wpi=sin(theta)
    wr=1.0d0+wpr
    wi=wpi
    n2p3=n+3
    do i=2,n/4
       i1=2*i-1
       i2=i1+1
       i3=n2p3-i2
       i4=i3+1
       wrs=sngl(wr)
       wis=sngl(wi)
       h1r=c1*(data(i1)+data(i3))
       h1i=c1*(data(i2)-data(i4))
       h2r=-c2*(data(i2)+data(i4))
       h2i=c2*(data(i1)-data(i3))
       data(i1)=h1r+wrs*h2r-wis*h2i
       data(i2)=h1i+wrs*h2i+wis*h2r
       data(i3)=h1r-wrs*h2r+wis*h2i
       data(i4)=-h1i+wrs*h2i+wis*h2r
       wtemp=wr
       wr=wr*wpr-wi*wpi+wr
       wi=wi*wpr+wtemp*wpi+wi
    enddo
    if (isign.eq.1) then
       h1r=data(1)
       data(1)=h1r+data(2)
       data(2)=h1r-data(2)
    else
       h1r=data(1)
       data(1)=c1*(h1r+data(2))
       data(2)=c1*(h1r-data(2))
       call rfour1(data,n/2,-1)
    endif
    return
  end subroutine realft

  !This routine is modified version w/ respect to original NumRec!
  !Below is the original one renamed rfour1 (real input twice the size)
  subroutine cfour1(Fct,nn,isign)
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
       Fct(i) = cmplx(Fcttmp(2*i-1),Fcttmp(2*i),8)
    enddo
  end subroutine cfour1

  subroutine rfour1(data,nn,isign)
    integer :: isign,nn
    real(8) :: data(2*nn)
    integer :: i,istep,j,m,mmax,n
    real(8) ::  tempi,tempr
    real(8) ::  theta,wi,wpi,wpr,wr,wtemp
    n=2*nn
    j=1
    do i=1,n,2
       if(j.gt.i)then
          tempr=data(j)
          tempi=data(j+1)
          data(j)=data(i)
          data(j+1)=data(i+1)
          data(i)=tempr
          data(i+1)=tempi
       endif
       m=nn
1      if ((m.ge.2).and.(j.gt.m)) then
          j=j-m
          m=m/2
          goto 1
       endif
       j=j+m
    enddo
    mmax=2
2   if (n.gt.mmax) then
       istep=2*mmax
       theta=6.28318530717959d0/(isign*mmax)
       wpr=-2.d0*sin(0.5d0*theta)**2
       wpi=sin(theta)
       wr=1.d0
       wi=0.d0
       do m=1,mmax,2
          do i=m,n,istep
             j=i+mmax
             tempr=sngl(wr)*data(j)-sngl(wi)*data(j+1)
             tempi=sngl(wr)*data(j+1)+sngl(wi)*data(j)
             data(j)=data(i)-tempr
             data(j+1)=data(i+1)-tempi
             data(i)=data(i)+tempr
             data(i+1)=data(i+1)+tempi
          enddo
          wtemp=wr
          wr=wr*wpr-wi*wpi+wr
          wi=wi*wpr+wtemp*wpi+wi
       enddo
       mmax=istep
       goto 2
    endif
    return
  end subroutine rfour1

  !*******************************************************************
  !*******************************************************************
  !*******************************************************************


end module FFT_NR
