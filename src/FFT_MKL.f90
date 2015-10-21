include "mkl_dfti.f90"
module FFT_MKL
  use MKL_DFTI
  use MKL_DFT_TYPE
  !use MKL_TRIG_TRANSFORMS
  !USE CONSTANTS, only: pi
  !USE ARRAYS, only: linspace
  !USE INTERPOLATE, only: linear_spline,cubic_spline
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


  integer                        :: status
  type(DFTI_DESCRIPTOR), pointer :: Handle


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
       func_out = ifftshift(ftmp)
    else
       func_in  = ifftshift(ftmp)
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
       func_out = ifftshift(ftmp)
    else
       func_in  = ifftshift(ftmp)
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
    Status   = DftiCreateDescriptor(Handle,DFTI_DOUBLE,DFTI_REAL,1,size(func))
    Status   = DftiCommitDescriptor(Handle)
    Status   = DftiComputeForward(Handle,func)
    Status   = DftiFreeDescriptor(Handle)
  end subroutine rfft_1d_forward
  !
  subroutine cfft_1d_forward(func)
    complex(8),dimension(:),intent(inout) :: func
    Status   = DftiCreateDescriptor(Handle,DFTI_DOUBLE,DFTI_COMPLEX,1,size(func))
    Status   = DftiCommitDescriptor(Handle)
    Status   = DftiComputeForward(Handle,func)
    Status   = DftiFreeDescriptor(Handle)
  end subroutine cfft_1d_forward


  !*********************************************************************
  !                       BACKWARD TRANSFORM
  !+-------------------------------------------------------------------+
  !PURPOSE  : Evaluate 1-DIM backward FFT using MKL routines.
  ! The returned real/complex array contains y(0), y(1),..., y(n-1) 
  !+-------------------------------------------------------------------+
  subroutine rfft_1d_backward(func)
    real(8),dimension(:),intent(inout) :: func
    Status = DftiCreateDescriptor(Handle,DFTI_DOUBLE,DFTI_REAL,1,size(func))
    Status = DftiSetValue(Handle, DFTI_BACKWARD_SCALE, 1.d0/dble(size(func)))
    Status = DftiCommitDescriptor(Handle)
    Status = DftiComputeBackward(Handle,func)
    Status = DftiFreeDescriptor(Handle)
  end subroutine rfft_1d_backward
  !
  subroutine cfft_1d_backward(func)
    complex(8),dimension(:),intent(inout) :: func
    Status = DftiCreateDescriptor(Handle,DFTI_DOUBLE,DFTI_COMPLEX,1,size(func))
    Status = DftiSetValue(Handle, DFTI_BACKWARD_SCALE, 1.d0/dble(size(func)))
    Status = DftiCommitDescriptor(Handle)
    Status = DftiComputeBackward(Handle,func)
    Status = DftiFreeDescriptor(Handle)
  end subroutine cfft_1d_backward






  !*********************************************************************
  !                           HELPER FUNCTIONS:
  !*********************************************************************  
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


end module FFT_MKL
