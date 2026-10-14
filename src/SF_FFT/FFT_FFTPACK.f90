MODULE FFT_FFTPACK
!:synopsis: Fast Fourier transforms of real and complex arrays, built on the FFTPACK 5.1 routines
!Wrappers of the FFTPACK 5.1 routines for the fast Fourier transform (FFT) of real and complex
!arrays, with the sine and cosine transforms and a few helper functions. All the transforms
!work in place, and the program stops if the FFTPACK routines report an error.
!
!The module provides, as generic interfaces:
!
!* :f:func_inline:`fft`, :f:func_inline:`ifft`: forward and backward 1-dimensional FFT;
!* :f:func_inline:`fft2`, :f:func_inline:`ifft2`: forward and backward 2-dimensional FFT;
!* :f:func_inline:`fftn`, :f:func_inline:`ifftn`: forward and backward 1-dimensional FFTs of many
!  sequences stored consecutively in a 1-dimensional array;
!* :f:func_inline:`cosft`, :f:func_inline:`icosft`, :f:func_inline:`sinft`, :f:func_inline:`isinft`
!  and their multiple-sequence versions :f:func_inline:`cosftn`, :f:func_inline:`icosftn`,
!  :f:func_inline:`sinftn`, :f:func_inline:`isinftn`: cosine and sine transforms of real arrays;
!* :f:func_inline:`tfft`, :f:func_inline:`itfft`: FFT between the time and the frequency domain,
!  for arrays with the zero of the domain at the center;
!* :f:func_inline:`fftshift`, :f:func_inline:`ifftshift`, :f:func_inline:`fftex`: helper functions.
!
!The FFT interfaces exist for ``real(8)`` and ``complex(8)`` arrays, the cosine and sine ones for
!``real(8)`` arrays only. The 1-dimensional forward and backward FFTs are not normalized, so
!that a forward transform followed by a backward one returns the original array multiplied by
!its length; in the 2-dimensional FFTs and in the cosine and sine transforms the normalization
!is carried by the forward transform. The FFT of a real array is returned in the FFTPACK
!packed form.
!
!The public-facing module, :f:mod:`sf_fft_fftpack` (file SF_FFT.f90), contains the same
!routines together with the utilities for discretized time/frequency signals.
!
  implicit none
  private


  !- - - - - - - - - - - - - - - - - - - - - - - - - - - - -!
  !Fast Fourier Transforms of time/frequency signal
  !- - - - - - - - - - - - - - - - - - - - - - - - - - - - -!
  interface tfft
  !This subroutine evaluates the fast Fourier transform of a real or complex array
  !:f:var:`func_in` of length :code:`N`, arranged with the zero of the domain at the
  !center of the array. It performs the sequence
  !
  !* :f:func_inline:`ifftshift` of the input
  !* :f:func_inline:`fft` (forward transform)
  !* :f:func_inline:`fftshift` of the result
  !* multiplication by :code:`N`
  !
  !If the optional argument :f:var:`func_out` is present it contains the result and 
  !:f:var:`func_in` is left untouched, otherwise :f:var:`func_in` is overwritten. 
  !In the real case the imaginary part of the transform is discarded.
  !
     module procedure d_tfft,c_tfft
  end interface tfft
  public :: tfft
  public :: d_tfft
  public :: c_tfft

  interface itfft
  !This subroutine evaluates the inverse fast Fourier transform of a real or complex 
  !array :f:var:`func_in` of length :code:`N`, and is the counterpart of :f:func_inline:`tfft`. 
  !It performs the sequence
  !
  !* :f:func_inline:`ifft` (backward transform)
  !* :f:func_inline:`fftex` (alternate the sign of the elements, shifting the output by half a period)
  !* :f:func_inline:`ifftshift` of the result
  !* division by :code:`N`
  !
  !If the optional argument :f:var:`func_out` is present it contains the result and 
  !:f:var:`func_in` is left untouched, otherwise :f:var:`func_in` is overwritten. 
  !In the real case the imaginary part of the transform is discarded.
  !
     module procedure d_itfft,c_itfft
  end interface itfft
  public :: itfft
  public :: d_itfft
  public :: c_itfft


  !- - - - - - - - - - - - - - - - - - - - - - - - - - - - -!
  !Fast Fourier Transforms
  !- - - - - - - - - - - - - - - - - - - - - - - - - - - - -!
  interface fft
  !This subroutine evaluates the forward 1-dimensional FFT of an array of lenght :code:`N` 
  !using the FFTPACK 5.1 routines :f:func_inline:`rfft1i` + :f:func_inline:`rfft1f`
  !for the real case and :f:func_inline:`cfft1i` + :f:func_inline:`cfft1f` for the 
  !complex case. These are called with the parameters
  !
  !* :code:`lensav` = :code:`N + int( log(dble(N))/log(2.d0) ) + 4`
  !* :code:`lenwrk` = :code:`N`
  !* :code:`lenr`  = :code:`N`
  !* :code:`inc` = :code:`1`
  ! 
  !The subroutine modifies the input array to contain
  !
  !* :code:`[y(0), y(1), ..., y(N/2),     y(-N/2+1), ...,   y(-1)]`   if `N` is even
  !* :code:`[y(0), y(1), ..., y((N-1)/2), y(-(N-1)/2), ..., y(-1)]`   if `N` is odd
  !
  !where
  !
  !:math:`y(j) = \sum_{k=0}^{N-1} x(k) \cdot e^{-i 2\pi/N \cdot j \cdot k}`  
  !
  !for :math:`j \in [0,N-1]`
  !
     module procedure rfft_1d_forward,cfft_1d_forward
  end interface fft
  public :: fft
  public :: rfft_1d_forward
  public :: cfft_1d_forward

  interface ifft
  !This subroutine evaluates the backward 1-dimensional FFT of an array of lenght :code:`N` 
  !using the FFTPACK 5.1 routines :f:func_inline:`rfft1i` + :f:func_inline:`rfft1b`
  !for the real case and :f:func_inline:`cfft1i` + :f:func_inline:`cfft1b` for the 
  !complex case. These are called with the parameters
  !
  !* :code:`lensav` = :code:`N + int( log(dble(N))/log(2.d0) ) + 4`
  !* :code:`lenwrk` = :code:`N`
  !* :code:`lenr`  = :code:`N`
  !* :code:`inc` = :code:`1`
  ! 
  !The subroutine modifies the input array to contain
  !
  !* :code:`[y(0), y(1), ..., y(N/2),     y(-N/2+1), ...,   y(-1)]`   if `N` is even
  !* :code:`[y(0), y(1), ..., y((N-1)/2), y(-(N-1)/2), ..., y(-1)]`   if `N` is odd
  !
  !where
  !
  !:math:`y(j) = \sum_{k=-N/2}^{N/2-1} x(k) \cdot a^{i 2\pi/N \cdot j \cdot k}`
  !
  !for :math:`j \in [0,N-1]`
  !
     module procedure rfft_1d_backward,cfft_1d_backward
  end interface ifft
  public :: ifft
  public :: rfft_1d_backward
  public :: cfft_1d_backward

  !- - - - - - - - - - - - - - - - - - - - - - - - - - - - -!

  interface fft2
  !This subroutine evaluates the forward 2-dimensional FFT of an array of size :code:`(L,M)`
  !using the FFTPACK 5.1 routines :f:func_inline:`rfft2i` + :f:func_inline:`rfft2f`
  !for the real case and :f:func_inline:`cfft2i` + :f:func_inline:`cfft2f` for the
  !complex case. These are called with the parameters
  !
  !* :code:`ldim` = :code:`L`
  !* :code:`lensav` = :code:`L + 3*M + int( log(dble(N))/log(2.d0) ) + 2*int( log(dble(M))/log(2.d0) ) + 12` (real)
  !* :code:`lensav` = :code:`2*(L+M) + int( log(dble(N))/log(2.d0) ) + int( log(dble(M))/log(2.d0) ) + 8` (complex)
  !* :code:`lenwrk` = :code:`M*(L+1)` (real), :code:`2*L*M` (complex)
  !
  !The subroutine modifies the input array to contain
  !
  !:math:`y(i,j) = \frac{1}{LM} \sum_{l=0}^{L-1} \sum_{m=0}^{M-1} x(l,m) \cdot e^{-i 2\pi (i l/L + j m/M)}`
  !
  !for :math:`i \in [0,L-1]`, :math:`j \in [0,M-1]`. The real subroutine returns the
  !FFTPACK packed form of the transform of a real array.
  !
     module procedure rfft_2d_forward,cfft_2d_forward
  end interface fft2
  public :: fft2
  public :: rfft_2d_forward
  public :: cfft_2d_forward

  interface ifft2
  !This subroutine evaluates the backward 2-dimensional FFT of an array of size :code:`(L,M)`
  !using the FFTPACK 5.1 routines :f:func_inline:`rfft2i` + :f:func_inline:`rfft2b`
  !for the real case and :f:func_inline:`cfft2i` + :f:func_inline:`cfft2b` for the
  !complex case. These are called with the same parameters as :f:func_inline:`fft2`.
  !
  !The subroutine modifies the input array to contain
  !
  !:math:`x(i,j) = \sum_{l=0}^{L-1} \sum_{m=0}^{M-1} y(l,m) \cdot e^{i 2\pi (i l/L + j m/M)}`
  !
  !for :math:`i \in [0,L-1]`, :math:`j \in [0,M-1]`.
  !
     module procedure rfft_2d_backward,cfft_2d_backward
  end interface ifft2
  public :: ifft2
  public :: rfft_2d_backward
  public :: cfft_2d_backward

  !- - - - - - - - - - - - - - - - - - - - - - - - - - - - -!

  interface fftn
  !This subroutine evaluates :code:`lot` forward 1-dimensional FFTs of length :code:`N`
  !on the sequences stored consecutively in the 1-dimensional array :f:var:`func`,
  !using the FFTPACK 5.1 routines :f:func_inline:`rfftmi` + :f:func_inline:`rfftmf`
  !for the real case and :f:func_inline:`cfftmi` + :f:func_inline:`cfftmf` for the
  !complex case. These are called with the parameters
  !
  !* :code:`lenr` = :code:`N*lot`
  !* :code:`lensav` = :code:`N + int( log(dble(N))/log(2.d0) ) + 4` (real), :code:`2*N + int( log(dble(N))/log(2.d0) ) + 4` (complex)
  !* :code:`lenwrk` = :code:`N*lot` (real), :code:`2*N*lot` (complex)
  !* :code:`jump` = :code:`N`
  !* :code:`inc` = :code:`1`
  !
  !For the :math:`l`-th sequence, :math:`l \in [0,lot-1]`, the subroutine modifies the input array to contain
  !
  !:math:`y(l N + j) = \sum_{k=0}^{N-1} x(l N + k) \cdot e^{-i 2\pi/N \cdot j \cdot k}`
  !
  !for :math:`j \in [0,N-1]`. The program stops if :code:`mod(N*lot,size(func))/=0`.
  !
     module procedure rfft_nd_forward,cfft_nd_forward
  end interface fftn
  public :: fftn
  public :: rfft_nd_forward
  public :: cfft_nd_forward

  interface ifftn
  !This subroutine evaluates :code:`lot` backward 1-dimensional FFTs of length :code:`N`
  !on the sequences stored consecutively in the 1-dimensional array :f:var:`func`,
  !using the FFTPACK 5.1 routines :f:func_inline:`rfftmi` + :f:func_inline:`rfftmb`
  !for the real case and :f:func_inline:`cfftmi` + :f:func_inline:`cfftmb` for the
  !complex case. These are called with the same parameters as :f:func_inline:`fftn`.
  !
  !For the :math:`l`-th sequence, :math:`l \in [0,lot-1]`, the subroutine modifies the input array to contain
  !
  !:math:`x(l N + j) = \sum_{k=0}^{N-1} y(l N + k) \cdot e^{i 2\pi/N \cdot j \cdot k}`
  !
  !for :math:`j \in [0,N-1]`.
  !
     module procedure rfft_nd_backward,cfft_nd_backward
  end interface ifftn
  public :: ifftn
  public :: rfft_nd_backward
  public :: cfft_nd_backward

  !- - - - - - - - - - - - - - - - - - - - - - - - - - - - -!

  interface cosft
  !This subroutine evaluates the forward 1-dimensional cosine FFT of a real array of
  !lenght :code:`N` using the FFTPACK 5.1 routines :f:func_inline:`cost1i` + :f:func_inline:`cost1f`.
  !These are called with the parameters
  !
  !* :code:`lensav` = :code:`2*N + int( log(dble(N))/log(2.d0) ) + 4`
  !* :code:`lenwrk` = :code:`N-1`
  !* :code:`lenr`  = :code:`N`
  !* :code:`inc` = :code:`1`
  !
  !The subroutine modifies the input array in place. The transform is the inverse of
  !:f:func_inline:`icosft`, i.e. the normalization is carried by the forward transform.
  !
     module procedure cost_1d_forward
  end interface cosft
  public :: cosft
  public :: cost_1d_forward

  interface icosft
  !This subroutine evaluates the backward 1-dimensional cosine FFT of a real array of
  !lenght :code:`N` using the FFTPACK 5.1 routines :f:func_inline:`cost1i` + :f:func_inline:`cost1b`.
  !These are called with the same parameters as :f:func_inline:`cosft`.
  !
  !The subroutine modifies the input array to contain
  !
  !:math:`x(i) = y(1) + (-1)^{i-1} y(N) + 2 \sum_{k=2}^{N-1} y(k) \cdot \cos\left(\frac{\pi (k-1)(i-1)}{N-1}\right)`
  !
  !for :math:`i \in [1,N]`.
  !
     module procedure cost_1d_backward
  end interface icosft
  public :: icosft
  public :: cost_1d_backward

  !- - - - - - - - - - - - - - - - - - - - - - - - - - - - -!

  interface cosftn
  !This subroutine evaluates :code:`lot` forward cosine FFTs of length :code:`N`
  !on the sequences stored consecutively in the 1-dimensional array :f:var:`func`,
  !using the FFTPACK 5.1 routines :f:func_inline:`costmi` + :f:func_inline:`costmf`.
  !These are called with the parameters
  !
  !* :code:`lenr` = :code:`N*lot`
  !* :code:`lensav` = :code:`2*N + int( log(dble(N))/log(2.d0) ) + 4`
  !* :code:`lenwrk` = :code:`lot*(N+1)`
  !* :code:`jump` = :code:`N`
  !* :code:`inc` = :code:`1`
  !
  !The program stops if :code:`mod(N*lot,size(func))/=0`.
  !
     module procedure cost_Nd_forward
  end interface cosftn
  public :: cosftn
  public :: cost_nd_forward

  interface icosftn
  !This subroutine evaluates :code:`lot` backward cosine FFTs of length :code:`N`
  !on the sequences stored consecutively in the 1-dimensional array :f:var:`func`,
  !using the FFTPACK 5.1 routines :f:func_inline:`costmi` + :f:func_inline:`costmb`.
  !These are called with the same parameters as :f:func_inline:`cosftn`.
  !
     module procedure cost_Nd_backward
  end interface icosftn
  public :: icosftn
  public :: cost_Nd_backward

  !- - - - - - - - - - - - - - - - - - - - - - - - - - - - -!

  interface sinft
  !This subroutine evaluates the forward 1-dimensional sine FFT of a real array of
  !lenght :code:`N` using the FFTPACK 5.1 routines :f:func_inline:`sint1i` + :f:func_inline:`sint1f`.
  !These are called with the parameters
  !
  !* :code:`lensav` = :code:`N/2 + N + int( log(dble(N))/log(2.d0) ) + 4`
  !* :code:`lenwrk` = :code:`2*(N+1)`
  !* :code:`lenr`  = :code:`N`
  !* :code:`inc` = :code:`1`
  !
  !The subroutine modifies the input array in place. The transform is the inverse of
  !:f:func_inline:`isinft`, i.e. the normalization is carried by the forward transform.
  !
     module procedure sint_1d_forward
  end interface sinft
  public :: sinft
  public :: sint_1d_forward

  interface isinft
  !This subroutine evaluates the backward 1-dimensional sine FFT of a real array of
  !lenght :code:`N` using the FFTPACK 5.1 routines :f:func_inline:`sint1i` + :f:func_inline:`sint1b`.
  !These are called with the same parameters as :f:func_inline:`sinft`.
  !
  !The subroutine modifies the input array to contain
  !
  !:math:`x(i) = 2 \sum_{k=1}^{N} y(k) \cdot \sin\left(\frac{\pi k i}{N+1}\right)`
  !
  !for :math:`i \in [1,N]`.
  !
     module procedure sint_1d_backward
  end interface isinft
  public :: isinft
  public :: sint_1d_backward

  !- - - - - - - - - - - - - - - - - - - - - - - - - - - - -!

  interface sinftn
  !This subroutine evaluates :code:`lot` forward sine FFTs of length :code:`N`
  !on the sequences stored consecutively in the 1-dimensional array :f:var:`func`,
  !using the FFTPACK 5.1 routines :f:func_inline:`sintmi` + :f:func_inline:`sintmf`.
  !These are called with the parameters
  !
  !* :code:`lenr` = :code:`N*lot`
  !* :code:`lensav` = :code:`N/2 + N + int( log(dble(N))/log(2.d0) ) + 4`
  !* :code:`lenwrk` = :code:`2*lot*(N+2)`
  !* :code:`jump` = :code:`N`
  !* :code:`inc` = :code:`1`
  !
  !The program stops if :code:`mod(N*lot,size(func))/=0`.
  !
     module procedure sint_Nd_forward
  end interface sinftn
  public :: sinftn
  public :: sint_nd_forward

  interface isinftn
  !This subroutine evaluates :code:`lot` backward sine FFTs of length :code:`N`
  !on the sequences stored consecutively in the 1-dimensional array :f:var:`func`,
  !using the FFTPACK 5.1 routines :f:func_inline:`sintmi` + :f:func_inline:`sintmb`.
  !These are called with the same parameters as :f:func_inline:`sinftn`.
  !
     module procedure sint_Nd_backward
  end interface isinftn
  public :: isinftn
  public :: sint_Nd_backward

  !- - - - - - - - - - - - - - - - - - - - - - - - - - - - -!





  !- - - - - - - - - - - - - - - - - - - - - - - - - - - - -!
  !HELPER FUNCTIONS:
  !- - - - - - - - - - - - - - - - - - - - - - - - - - - - -!
  interface fftshift
  !This function shifts the zero-frequency component of the output of :f:func_inline:`fft`
  !to the center of the array. The input
  !
  !* :code:`[y(0), y(1), ..., y(N/2),     y(-N/2+1), ...,   y(-1)]`   if `N` is even
  !* :code:`[y(0), y(1), ..., y((N-1)/2), y(-(N-1)/2), ..., y(-1)]`   if `N` is odd
  !
  !is returned as
  !
  !* :code:`[y(-N/2), ..., y(-1), y(0), y(1), ..., y(N/2-1)]`          if `N` is even
  !* :code:`[y(-(N-1)/2), ..., y(-1), y(0), y(1), ..., y((N-1)/2)]`    if `N` is odd
  !
     module procedure rfft_1d_shift,cfft_1d_shift
  end interface fftshift
  public :: fftshift
  public :: rfft_1d_shift
  public :: cfft_1d_shift

  interface ifftshift
  !This function is the inverse of :f:func_inline:`fftshift`: it moves the
  !zero-frequency component from the center of the array back to the first position.
  !
     module procedure rfft_1d_ishift,cfft_1d_ishift
  end interface ifftshift
  public :: ifftshift
  public :: rfft_1d_ishift
  public :: cfft_1d_ishift

  interface fftex
  !This subroutine alternates the sign of the elements of an array,
  !:math:`f(i) \to (-1)^{i-1} f(i)`. It is used by :f:func_inline:`itfft` to shift the output of
  !:f:func_inline:`ifft` by half a period.
  !
     module procedure rfft_1d_ex,cfft_1d_ex
  end interface fftex
  public :: fftex
  public :: rfft_1d_ex
  public :: cfft_1d_ex
  !- - - - - - - - - - - - - - - - - - - - - - - - - - - - -!



contains


  !*********************************************************************
  !               TIME <==> FREQUENCY DOMAIN FFT:
  !*********************************************************************  
  subroutine d_tfft(func_in,func_out)
  !Fast Fourier transform from the time to the frequency domain of a real array, see :f:func_inline:`tfft`.
  !
  !Real version: the imaginary part of the transform is discarded.
  !
    real(8),dimension(:)                         :: func_in  !Array of length :code:`N` to transform, with the zero of the time
                                                             ! domain at the center. Overwritten by the transform if
                                                             ! :f:var:`func_out` is absent
    real(8),dimension(size(func_in)),optional    :: func_out  !If present, receives the transform and :f:var:`func_in` is left
                                                              ! untouched
    complex(8),dimension(size(func_in))             :: ftmp
    ftmp = ifftshift(func_in)  
    call fft(ftmp)
    if(present(func_out))then
       func_out = fftshift(ftmp)*size(ftmp)
    else
       func_in  = fftshift(ftmp)*size(ftmp)
    endif
  end subroutine d_tfft
  !
  subroutine c_tfft(func_in,func_out)
  !Fast Fourier transform from the time to the frequency domain of a complex array, see :f:func_inline:`tfft`.
  !
  !Complex version.
  !
    complex(8),dimension(:)                         :: func_in  !Array of length :code:`N` to transform, with the zero of the
                                                                ! time domain at the center. Overwritten by the transform if
                                                                ! :f:var:`func_out` is absent
    complex(8),dimension(size(func_in)),optional    :: func_out  !If present, receives the transform and :f:var:`func_in` is
                                                                 ! left untouched
    complex(8),dimension(size(func_in))             :: ftmp
    ftmp = ifftshift(func_in)  
    call fft(ftmp)
    if(present(func_out))then
       func_out = fftshift(ftmp)*size(ftmp)
    else
       func_in  = fftshift(ftmp)*size(ftmp)
    endif
  end subroutine c_tfft

  subroutine d_itfft(func_in,func_out)
  !Fast Fourier transform from the frequency to the time domain of a real array, see :f:func_inline:`itfft`.
  !
  !Real version: the imaginary part of the transform is discarded.
  !
    real(8),dimension(:)                      :: func_in  !Array of length :code:`N` to transform, with the zero of the
                                                          ! frequency domain at the center. Overwritten by the inverse
                                                          ! transform if :f:var:`func_out` is absent
    real(8),dimension(size(func_in)),optional :: func_out  !If present, receives the inverse transform and :f:var:`func_in` is
                                                           ! left untouched
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
  !Fast Fourier transform from the frequency to the time domain of a complex array, see :f:func_inline:`itfft`.
  !
  !Complex version.
  !
    complex(8),dimension(:)                      :: func_in  !Array of length :code:`N` to transform, with the zero of the
                                                             ! frequency domain at the center. Overwritten by the inverse
                                                             ! transform if :f:var:`func_out` is absent
    complex(8),dimension(size(func_in)),optional :: func_out  !If present, receives the inverse transform and :f:var:`func_in`
                                                              ! is left untouched
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
  !PURPOSE  : Evaluate forward 1-DIM FFT using FFTPACK 5.1 routines.
  !out is :
  ! [y(0), y(1), ..., y(N/2),     y(-N/2+1), ...,   y(-1)]   if N is even
  ! [y(0), y(1), ..., y((N-1)/2), y(-(N-1)/2), ..., y(-1)]   if N is odd
  !where:
  ! y(j) = sum[k=0..N-1] x[k] * exp(-xi*j*k* 2*pi/N), j = 0..N-1
  ! For index J*INC+1 where J=0,...,N-1 (that is, for the Jth
  !          element of the sequence):
  !               N-1
  ! C(J*INC+1) =  SUM C(K*INC+1)*EXP(-XI*J*K*2*PI/N)
  !               K=0
  !+-------------------------------------------------------------------+
  subroutine rfft_1d_forward(func)
  !This subroutine evaluates the forward 1-dimensional FFT of a real array :f:var:`func`
  !of length :code:`N`, in place, using the FFTPACK 5.1 routines :f:func_inline:`rfft1i` + 
  !:f:func_inline:`rfft1f`. These are called with the parameters
  !
  !* :code:`lensav` = :code:`N + int( log(dble(N))/log(2.d0) ) + 4`
  !* :code:`lenwrk` = :code:`N`
  !* :code:`lenr`  = :code:`N`
  !* :code:`inc` = :code:`1`
  !
  !The transform is not normalized. Since the transform of a real sequence is 
  !Hermitian, the subroutine returns the FFTPACK packed form of
  !
  !:math:`y(j) = \sum_{k=0}^{N-1} x(k) \cdot e^{-i 2\pi/N \cdot j \cdot k}`
  !
  !that is, the array contains
  !
  !* :code:`[Re y(0), Re y(1), Im y(1), ..., Re y(N/2-1), Im y(N/2-1), Re y(N/2)]`   if `N` is even
  !* :code:`[Re y(0), Re y(1), Im y(1), ..., Re y((N-1)/2), Im y((N-1)/2)]`           if `N` is odd
  !
  !The program stops if the FFTPACK routines return an error code.
  !
    real(8),dimension(:),intent(inout) :: func      !Input/output array of size :code:`N`
    real(8),dimension(:),allocatable   :: wsave,work
    integer                            :: N,lenwrk,lensav,lenr,inc,ier
    N      = size(func)
    lenwrk = N
    lensav = N + int( log(dble(N))/log(2.d0) ) + 4
    lenr   = N
    inc    = 1
    allocate(wsave(lensav))
    allocate(work(lenwrk))
    call rfft1i(N,wsave,lensav,ier)
    if(ier==2)stop "rfft_1d_forward: LENSAV not big enough"
    call rfft1f(N,inc,func,lenr,wsave,lensav,work,lenwrk,ier)
    deallocate(wsave,work)
    select case(ier)
    case (0)
       return
    case (1)
       stop "rfft_1d_forward: LENR not big enough"
    case (2)
       stop "rfft_1d_forward: LENSAV not big enough"
    case (3)
       stop "rfft_1d_forward: LENWRK not big enough"
    case (20)
       stop "rfft_1d_forward: input error returned by lower level routine"
    end select
  end subroutine rfft_1d_forward
  !
  subroutine cfft_1d_forward(func)
  !This subroutine evaluates the forward 1-dimensional FFT of a complex array :f:var:`func`
  !of length :code:`N`, in place, using the FFTPACK 5.1 routines :f:func_inline:`cfft1i` + 
  !:f:func_inline:`cfft1f`. These are called with the parameters
  !
  !* :code:`lensav` = :code:`2*N + int( log(dble(N))/log(2.d0) ) + 4`
  !* :code:`lenwrk` = :code:`2*N`
  !* :code:`lenc`  = :code:`N`
  !* :code:`inc` = :code:`1`
  !
  !The transform is not normalized. The subroutine modifies the input array to contain
  !
  !:math:`y(j) = \sum_{k=0}^{N-1} x(k) \cdot e^{-i 2\pi/N \cdot j \cdot k}`
  !
  !for :math:`j \in [0,N-1]`, i.e. the array is returned in the aliased order
  !
  !* :code:`[y(0), y(1), ..., y(N/2),     y(-N/2+1), ...,   y(-1)]`   if `N` is even
  !* :code:`[y(0), y(1), ..., y((N-1)/2), y(-(N-1)/2), ..., y(-1)]`   if `N` is odd
  !
  !The program stops if the FFTPACK routines return an error code.
  !
    complex(8),dimension(:),intent(inout) :: func  !Input/output array of size :code:`N`, overwritten by its transform
    real(8),dimension(:),allocatable      :: wsave,work
    integer                               :: N,lenwrk,lensav,lenc,inc,ier
    N      = size(func)
    lenwrk = 2*N
    lensav = 2*N + int( log(dble(N))/log(2.d0) ) + 4
    lenc   = N
    inc    = 1
    allocate(wsave(lensav))
    allocate(work(lenwrk))
    call cfft1i(N,wsave,lensav,ier)
    if(ier==2)stop "cfft_1d_forward: LENSAV not big enough"
    call cfft1f(N,inc,func,lenc,wsave,lensav,work,lenwrk,ier)
    deallocate(wsave,work)
    select case(ier)
    case (0)
       return
    case (1)
       stop "cfft_1d_forward: LENC not big enough"
    case (2)
       stop "cfft_1d_forward: LENSAV not big enough"
    case (3)
       stop "cfft_1d_forward: LENWRK not big enough"
    case (20)
       stop "cfft_1d_forward: input error returned by lower level routine"
    end select
  end subroutine cfft_1d_forward


  !+-------------------------------------------------------------------+
  !PURPOSE  : Evaluate forward 2-DIM FFT using FFTPACK 5.1 routines.
  ! The full complex transform of r(i,j) is given by:
  !                    L-1  M-1
  ! C(I,J) =   1/(L*M)*SUM  SUM  C(L1,M1)*EXP(-XI*2*PI*(I*L1/L + J*M1/M))
  !                    L1=0 M1=0
  !
  !+-------------------------------------------------------------------+
  subroutine rfft_2d_forward(func)
  !This subroutine evaluates the forward 2-dimensional FFT of a real array :f:var:`func`
  !of size :code:`(L,M)`, in place, using the FFTPACK 5.1 routines :f:func_inline:`rfft2i` + 
  !:f:func_inline:`rfft2f`. These are called with the parameters
  !
  !* :code:`ldim` = :code:`L`
  !* :code:`lensav` = :code:`L + 3*M + int( log(dble(L))/log(2.d0) ) + 2*int( log(dble(M))/log(2.d0) ) + 12`
  !* :code:`lenwrk` = :code:`M*(L+1)`
  !
  !The subroutine computes
  !
  !:math:`y(i,j) = \frac{1}{LM} \sum_{l=0}^{L-1} \sum_{m=0}^{M-1} x(l,m) \cdot e^{-i 2\pi (i l/L + j m/M)}`
  !
  !for :math:`i \in [0,L-1]`, :math:`j \in [0,M-1]`. The transform of a real array is 
  !Hermitian, so the array is overwritten with the FFTPACK packed form of :math:`y`.
  !
  !The program stops if the FFTPACK routines return an error code.
  !
    real(8),dimension(:,:),intent(inout) :: func  !Input/output array of size :code:`(L,M)`, overwritten by the FFTPACK packed
                                                  ! form of its transform
    real(8),dimension(:),allocatable   :: wsave,work
    integer                            :: L,M,lenwrk,lensav,ldim,ier
    L      = size(func,1)
    M      = size(func,2)
    lenwrk = M*(L+1)
    lensav = L+3*M  + int(log(dble(L))/log(2.d0))+2*int(log(dble(M))/log(2.d0)) + 12
    ldim   = L
    allocate(wsave(lensav))
    allocate(work(lenwrk))
    call rfft2i(L,M,wsave,lensav,ier)
    if(ier==2)stop "rfft_2d_forward: LENSAV not big enough"
    if(ier==20)stop "rfft_2d_forward: input error returned by lower level routine"
    call rfft2f(ldim,L,M,func,wsave,lensav,work,lenwrk,ier)
    deallocate(wsave,work)
    select case(ier)
    case (0)
       return
    case (6)
       stop "rfft_2d_forward: LDIM is less than 2*INT((L+1)/2)"
    case (2)
       stop "rfft_2d_forward: LENSAV not big enough"
    case (3)
       stop "rfft_2d_forward: LENWRK not big enough"
    case (20)
       stop "rfft_2d_forward: input error returned by lower level routine"
    end select
  end subroutine rfft_2d_forward
  !
  subroutine cfft_2d_forward(func)
  !This subroutine evaluates the forward 2-dimensional FFT of a complex array :f:var:`func`
  !of size :code:`(L,M)`, in place, using the FFTPACK 5.1 routines :f:func_inline:`cfft2i` + 
  !:f:func_inline:`cfft2f`. These are called with the parameters
  !
  !* :code:`ldim` = :code:`L`
  !* :code:`lensav` = :code:`2*(L+M) + int( log(dble(L))/log(2.d0) ) + int( log(dble(M))/log(2.d0) ) + 8`
  !* :code:`lenwrk` = :code:`2*L*M`
  !
  !The subroutine modifies the input array to contain
  !
  !:math:`y(i,j) = \frac{1}{LM} \sum_{l=0}^{L-1} \sum_{m=0}^{M-1} x(l,m) \cdot e^{-i 2\pi (i l/L + j m/M)}`
  !
  !for :math:`i \in [0,L-1]`, :math:`j \in [0,M-1]`.
  !
  !The program stops if the FFTPACK routines return an error code.
  !
    complex(8),dimension(:,:),intent(inout) :: func  !Input/output array of size :code:`(L,M)`, overwritten by its transform
    real(8),dimension(:),allocatable      :: wsave,work
    integer                               :: L,M,lenwrk,lensav,ldim,inc,ier
    L      = size(func,1)
    M      = size(func,2)
    lenwrk = 2*L*M
    lensav = 2*(L+M) + int(log(dble(L))/log(2.d0))+int(log(dble(M))/log(2.d0)) + 8
    ldim   = L
    allocate(wsave(lensav))
    allocate(work(lenwrk))
    call cfft2i(L,M,wsave,lensav,ier)
    if(ier==2)stop "cfft_2d_forward: LENSAV not big enough"
    if(ier==20)stop "cfft_2d_forward: input error returned by lower level routine"
    call cfft2f(ldim,L,M,func,wsave,lensav,work,lenwrk,ier)
    deallocate(wsave,work)
    select case(ier)
    case (0)
       return
    case (5)
       stop "cfft_2d_forward: L > LDIM"
    case (2)
       stop "cfft_2d_forward: LENSAV not big enough"
    case (3)
       stop "cfft_2d_forward: LENWRK not big enough"
    case (20)
       stop "cfft_2d_forward: input error returned by lower level routine"
    end select
  end subroutine cfft_2d_forward


  !+-------------------------------------------------------------------+
  !PURPOSE  : Evaluate forward N-DIM FFT using FFTPACK 5.1 routines.
  ! For index L*JUMP + J*INC +1 where J=0,...,N-1 and 
  ! L=0,...,LOT-1, (that is, for the Jth element of the Lth sequence),
  !                       N-1
  !  C(L*JUMP+J*INC+1) = SUM C(L*JUMP+K*INC+1)*EXP(-XI*J*K*2*PI/N)
  !                       K=0
  !+-------------------------------------------------------------------+
  subroutine rfft_nd_forward(func,n,lot)
  !This subroutine evaluates :code:`lot` forward 1-dimensional FFTs of length :code:`N`
  !on the real sequences stored consecutively in the 1-dimensional array :f:var:`func`,
  !in place, using the FFTPACK 5.1 routines :f:func_inline:`rfftmi` + :f:func_inline:`rfftmf`.
  !These are called with the parameters
  !
  !* :code:`lenr` = :code:`N*lot`
  !* :code:`lensav` = :code:`N + int( log(dble(N))/log(2.d0) ) + 4`
  !* :code:`lenwrk` = :code:`N*lot`
  !* :code:`jump` = :code:`N`
  !* :code:`inc` = :code:`1`
  !
  !The transforms are not normalized. For the :math:`l`-th sequence, :math:`l \in [0,lot-1]`,
  !the subroutine computes
  !
  !:math:`y(l N + j) = \sum_{k=0}^{N-1} x(l N + k) \cdot e^{-i 2\pi/N \cdot j \cdot k}`
  !
  !for :math:`j \in [0,N-1]`. Since the transform of a real sequence is Hermitian, each 
  !sequence is overwritten with the FFTPACK packed form of :math:`y`, as in 
  !:f:func_inline:`rfft_1d_forward`.
  !
  !The program stops if :code:`mod(N*lot,size(func))/=0`, or if the FFTPACK routines 
  !return an error code.
  !
    real(8),dimension(:),intent(inout) :: func  !Input/output array storing :code:`lot` sequences of length :code:`N`
                                                ! consecutively, overwritten by their transforms (FFTPACK packed form)
    integer,intent(in)                 :: N  !Length of each sequence
    integer,intent(in)                 :: lot  !Number of sequences
    real(8),dimension(:),allocatable   :: wsave,work
    integer                            :: L,lenwrk,lensav,lenr,ier,inc,jump
    L=size(func)
    if(mod(N*lot,L)/=0)stop "rfft_Nd_forward: incommensurate values of parameters." 
    !
    lenr   = N*lot
    lenwrk = N*lot
    lensav = N  + int(log(dble(N))/log(2.d0)) + 4
    jump   = N
    inc    = 1
    allocate(wsave(lensav))
    allocate(work(lenwrk))
    call rfftmi(N,wsave,lensav,ier)
    if(ier==2)stop "rfft_Nd_forward: LENSAV not big enough"
    call rfftmf(lot,jump,N,inc,func,lenr,wsave,lensav,work,lenwrk,ier)
    deallocate(wsave,work)
    select case(ier)
    case (0)
       return
    case (1)
       stop "rfft_Nd_forward: LENR   not big enough"
    case (2)
       stop "rfft_Nd_forward: LENSAV not big enough"
    case (3)
       stop "rfft_Nd_forward: LENWRK not big enough"
    case (4)
       stop "rfft_Nd_forward: INC,JUMP,N,LOT are not consistent"
    end select
  end subroutine rfft_nd_forward
  !
  subroutine cfft_nd_forward(func,n,lot)
  !This subroutine evaluates :code:`lot` forward 1-dimensional FFTs of length :code:`N`
  !on the complex sequences stored consecutively in the 1-dimensional array :f:var:`func`,
  !in place, using the FFTPACK 5.1 routines :f:func_inline:`cfftmi` + :f:func_inline:`cfftmf`.
  !These are called with the parameters
  !
  !* :code:`lenc` = :code:`N*lot`
  !* :code:`lensav` = :code:`2*N + int( log(dble(N))/log(2.d0) ) + 4`
  !* :code:`lenwrk` = :code:`2*N*lot`
  !* :code:`jump` = :code:`N`
  !* :code:`inc` = :code:`1`
  !
  !The transforms are not normalized. For the :math:`l`-th sequence, :math:`l \in [0,lot-1]`,
  !the subroutine modifies the input array to contain
  !
  !:math:`y(l N + j) = \sum_{k=0}^{N-1} x(l N + k) \cdot e^{-i 2\pi/N \cdot j \cdot k}`
  !
  !for :math:`j \in [0,N-1]`. Each sequence is returned in the aliased order described 
  !in :f:func_inline:`cfft_1d_forward`.
  !
  !The program stops if :code:`mod(N*lot,size(func))/=0`, or if the FFTPACK routines 
  !return an error code.
  !
    complex(8),dimension(:),intent(inout) :: func  !Input/output array storing :code:`lot` sequences of length :code:`N`
                                                   ! consecutively, overwritten by their transforms
    integer,intent(in)                    :: n  !Length of each sequence
    integer,intent(in)                    :: lot  !Number of sequences
    real(8),dimension(:),allocatable      :: wsave,work
    integer                               :: L,lenwrk,lensav,lenc,ier,inc,jump
    L=size(func)
    if(mod(N*lot,L)/=0)stop "cfft_Nd_forward: incommensurate values of parameters." 
    lenc   = N*lot
    lenwrk = 2*N*lot
    lensav = 2*N + int(log(dble(N))/log(2.d0)) + 4
    jump   = N
    inc    = 1
    allocate(wsave(lensav))
    allocate(work(lenwrk))
    call cfftmi(N,wsave,lensav,ier)
    if(ier==2)stop "cfft_Nd_forward: LENSAV not big enough"
    call cfftmf(lot,jump,N,inc,func,lenc,wsave,lensav,work,lenwrk,ier)
    deallocate(wsave,work)
    select case(ier)
    case (0)
       return
    case (1)
       stop "cfft_Nd_forward: LENR   not big enough"
    case (2)
       stop "cfft_Nd_forward: LENSAV not big enough"
    case (3)
       stop "cfft_Nd_forward: LENWRK not big enough"
    case (4)
       stop "cfft_Nd_forward: INC,JUMP,N,LOT are not consistent"
    end select
  end subroutine cfft_nd_forward


  !+-------------------------------------------------------------------+
  !PURPOSE  : Evaluate forward 1-DIM COS-FFT using FFTPACK 5.1 routines.
  !+-------------------------------------------------------------------+
  subroutine cost_1d_forward(func)
  !This subroutine evaluates the forward 1-dimensional cosine FFT of a real array :f:var:`func`
  !of length :code:`N`, in place, using the FFTPACK 5.1 routines :f:func_inline:`cost1i` + 
  !:f:func_inline:`cost1f`. These are called with the parameters
  !
  !* :code:`lensav` = :code:`2*N + int( log(dble(N))/log(2.d0) ) + 4`
  !* :code:`lenwrk` = :code:`N-1`
  !* :code:`lenr`  = :code:`N`
  !* :code:`inc` = :code:`1`
  !
  !The transform is the inverse of :f:func_inline:`cost_1d_backward`, i.e. the 
  !normalization is carried by the forward transform. The subroutine modifies the 
  !input array to contain the coefficients :math:`y(k)` such that the backward 
  !transform :f:func_inline:`cost_1d_backward` reconstructs
  !
  !:math:`x(i) = y(1) + (-1)^{i-1} y(N) + 2 \sum_{k=2}^{N-1} y(k) \cdot \cos\left(\frac{\pi (k-1)(i-1)}{N-1}\right)`
  !
  !for :math:`i \in [1,N]`.
  !
  !The program stops if the FFTPACK routines return an error code.
  !
    real(8),dimension(:),intent(inout) :: func  !Input/output array of size :code:`N`, overwritten by its cosine transform
    real(8),dimension(:),allocatable   :: wsave,work
    integer                            :: N,lenwrk,lensav,lenr,inc,ier
    N      = size(func)
    lensav = 2*N + int( log(dble(N))/log(2.d0) ) + 4
    lenwrk = N-1
    lenr   = N
    inc    = 1
    allocate(wsave(lensav))
    allocate(work(lenwrk))
    call cost1i(N,wsave,lensav,ier)
    if(ier==2)stop "cost_1d_forward: LENSAV not big enough"
    if(ier==20)stop "cost_1d_forward: error returned by lower level routine"
    call cost1f(N,inc,func,lenr,wsave,lensav,work,lenwrk,ier)
    deallocate(wsave,work)
    select case(ier)
    case (0)
       return
    case (1)
       stop "cost_1d_forward: LENR not big enough"
    case (2)
       stop "cost_1d_forward: LENSAV not big enough"
    case (3)
       stop "cost_1d_forward: LENWRK not big enough"
    case (20)
       stop "cost_1d_forward: input error returned by lower level routine"
    end select
  end subroutine cost_1d_forward

  !+-------------------------------------------------------------------+
  !PURPOSE  : Evaluate forward 1-DIM SIN-FFT using FFTPACK 5.1 routines.
  !+-------------------------------------------------------------------+
  subroutine sint_1d_forward(func)
  !This subroutine evaluates the forward 1-dimensional sine FFT of a real array :f:var:`func`
  !of length :code:`N`, in place, using the FFTPACK 5.1 routines :f:func_inline:`sint1i` + 
  !:f:func_inline:`sint1f`. These are called with the parameters
  !
  !* :code:`lensav` = :code:`N/2 + N + int( log(dble(N))/log(2.d0) ) + 4`
  !* :code:`lenwrk` = :code:`2*(N+1)`
  !* :code:`lenr`  = :code:`N`
  !* :code:`inc` = :code:`1`
  !
  !The transform is the inverse of :f:func_inline:`sint_1d_backward`, i.e. the 
  !normalization is carried by the forward transform. The subroutine modifies the 
  !input array to contain the coefficients :math:`y(k)` such that the backward 
  !transform :f:func_inline:`sint_1d_backward` reconstructs
  !
  !:math:`x(i) = 2 \sum_{k=1}^{N} y(k) \cdot \sin\left(\frac{\pi k i}{N+1}\right)`
  !
  !for :math:`i \in [1,N]`.
  !
  !The program stops if the FFTPACK routines return an error code.
  !
    real(8),dimension(:),intent(inout) :: func  !Input/output array of size :code:`N`, overwritten by its sine transform
    real(8),dimension(:),allocatable   :: wsave,work
    integer                            :: N,lenwrk,lensav,lenr,inc,ier
    N      = size(func)
    lensav = N/2 + N + int( log(dble(N))/log(2.d0) ) + 4
    lenwrk = 2*(N+1)
    lenr   = N
    inc    = 1
    allocate(wsave(lensav))
    allocate(work(lenwrk))
    call sint1i(N,wsave,lensav,ier)
    if(ier==2)stop "sint_1d_forward: LENSAV not big enough"
    if(ier==20)stop "sint_1d_forward: error returned by lower level routine"
    call sint1f(N,inc,func,lenr,wsave,lensav,work,lenwrk,ier)
    deallocate(wsave,work)
    select case(ier)
    case (0)
       return
    case (1)
       stop "sint_1d_forward: LENR not big enough"
    case (2)
       stop "sint_1d_forward: LENSAV not big enough"
    case (3)
       stop "sint_1d_forward: LENWRK not big enough"
    case (20)
       stop "sint_1d_forward: input error returned by lower level routine"
    end select
  end subroutine sint_1d_forward

  !+-------------------------------------------------------------------+
  !PURPOSE  : Evaluate forward M-DIM COS-FFT using FFTPACK 5.1 routines.
  !+-------------------------------------------------------------------+
  subroutine cost_nd_forward(func,n,lot)
  !This subroutine evaluates :code:`lot` forward 1-dimensional cosine FFTs of length :code:`N`
  !on the real sequences stored consecutively in the 1-dimensional array :f:var:`func`,
  !in place, using the FFTPACK 5.1 routines :f:func_inline:`costmi` + :f:func_inline:`costmf`.
  !These are called with the parameters
  !
  !* :code:`lenr` = :code:`N*lot`
  !* :code:`lensav` = :code:`2*N + int( log(dble(N))/log(2.d0) ) + 4`
  !* :code:`lenwrk` = :code:`lot*(N+1)`
  !* :code:`jump` = :code:`N`
  !* :code:`inc` = :code:`1`
  !
  !For the :math:`l`-th sequence, :math:`l \in [0,lot-1]`, the subroutine acts as
  !:f:func_inline:`cost_1d_forward` on the elements :math:`l N + 1, ..., (l+1) N` of 
  !the array, so the transform is the inverse of :f:func_inline:`cost_nd_backward` 
  !and the normalization is carried by the forward transform.
  !
  !The program stops if :code:`mod(N*lot,size(func))/=0`, or if the FFTPACK routines 
  !return an error code.
  !
    real(8),dimension(:),intent(inout) :: func  !Input/output array storing :code:`lot` sequences of length :code:`N`
                                                ! consecutively, overwritten by their cosine transforms
    integer,intent(in)                 :: N  !Length of each sequence
    integer,intent(in)                 :: lot  !Number of sequences
    real(8),dimension(:),allocatable   :: wsave,work
    integer                            :: L,lenwrk,lensav,lenr,inc,ier,jump
    L=size(func)
    if(mod(N*lot,L)/=0)stop "cost_Nd_forward: incommensurate values of parameters." 
    !
    lenr   = N*lot
    lenwrk = lot*(N+1)
    lensav = 2*N  + int(log(dble(N))/log(2.d0)) + 4
    jump   = N
    inc    = 1
    allocate(wsave(lensav),work(lenwrk))
    call costmi(N,wsave,lensav,ier)
    if(ier==2)stop "cost_Nd_forward: LENSAV not big enough"
    if(ier==20)stop "cost_Nd_forward: error returned by lower level routine"
    call costmf(lot,jump,N,inc,func,lenr,wsave,lensav,work,lenwrk,ier)
    deallocate(wsave,work)
    select case(ier)
    case (0)
       return
    case (1)
       stop "cost_Nd_forward: LENR not big enough"
    case (2)
       stop "cost_Nd_forward: LENSAV not big enough"
    case (3)
       stop "cost_Nd_forward: LENWRK not big enough"
    case (4)
       stop "cost_Nd_forward: INC,JUMP,N,LOT are not consistent" 
    case (20)
       stop "cost_Nd_forward: input error returned by lower level routine"
    end select
  end subroutine cost_nd_forward

  !+-------------------------------------------------------------------+
  !PURPOSE  : Evaluate forward M-DIM SIN-FFT using FFTPACK 5.1 routines.
  !+-------------------------------------------------------------------+
  subroutine sint_nd_forward(func,n,lot)
  !This subroutine evaluates :code:`lot` forward 1-dimensional sine FFTs of length :code:`N`
  !on the real sequences stored consecutively in the 1-dimensional array :f:var:`func`,
  !in place, using the FFTPACK 5.1 routines :f:func_inline:`sintmi` + :f:func_inline:`sintmf`.
  !These are called with the parameters
  !
  !* :code:`lenr` = :code:`N*lot`
  !* :code:`lensav` = :code:`N/2 + N + int( log(dble(N))/log(2.d0) ) + 4`
  !* :code:`lenwrk` = :code:`2*lot*(N+2)`
  !* :code:`jump` = :code:`N`
  !* :code:`inc` = :code:`1`
  !
  !For the :math:`l`-th sequence, :math:`l \in [0,lot-1]`, the subroutine acts as
  !:f:func_inline:`sint_1d_forward` on the elements :math:`l N + 1, ..., (l+1) N` of 
  !the array, so the transform is the inverse of :f:func_inline:`sint_nd_backward` 
  !and the normalization is carried by the forward transform.
  !
  !The program stops if :code:`mod(N*lot,size(func))/=0`, or if the FFTPACK routines 
  !return an error code.
  !
    real(8),dimension(:),intent(inout) :: func  !Input/output array storing :code:`lot` sequences of length :code:`N`
                                                ! consecutively, overwritten by their sine transforms
    integer,intent(in)                 :: N  !Length of each sequence
    integer,intent(in)                 :: lot  !Number of sequences
    real(8),dimension(:),allocatable   :: wsave,work
    integer                            :: L,lenwrk,lensav,lenr,inc,ier,jump
    L=size(func)
    if(mod(N*lot,L)/=0)stop "cost_Nd_forward: incommensurate values of parameters." 
    !
    lenr   = N*lot
    lenwrk = lot*2*(N+2)
    lensav = N/2 + N  + int(log(dble(N))/log(2.d0)) + 4
    jump   = N
    inc    = 1
    allocate(wsave(lensav),work(lenwrk))
    call sintmi(N,wsave,lensav,ier)
    if(ier==2)stop "sint_Nd_forward: LENSAV not big enough"
    if(ier==20)stop "sint_Nd_forward: error returned by lower level routine"
    call sintmf(lot,jump,N,inc,func,lenr,wsave,lensav,work,lenwrk,ier)
    deallocate(wsave,work)
    select case(ier)
    case (0)
       return
    case (1)
       stop "sint_Nd_forward: LENR not big enough"
    case (2)
       stop "sint_Nd_forward: LENSAV not big enough"
    case (3)
       stop "sint_Nd_forward: LENWRK not big enough"
    case (4)
       stop "sint_Nd_forward: INC,JUMP,N,LOT are not consistent" 
    case (20)
       stop "sint_Nd_forward: input error returned by lower level routine"
    end select
  end subroutine sint_nd_forward



  !*********************************************************************
  !                       BACKWARD TRANSFORM
  !+-------------------------------------------------------------------+
  !PURPOSE  : Evaluate 1-DIM backward FFT using FFTPACK 5.1 routines.
  ! The returned real/complex array contains y(0), y(1),..., y(n-1) where
  ! y(j) = sum_{k=-N/2,...,N/2-1}(x(k) * exp(2*pi*XI*j*k/N))
  ! For index J*INC+1 where J=0,...,N-1,
  !              N-1
  ! C(J*INC+1) = SUM C(K*INC+1)*EXP(XI*J*K*2*PI/N)
  !              K=0
  !+-------------------------------------------------------------------+
  subroutine rfft_1d_backward(func)
  !This subroutine evaluates the backward 1-dimensional FFT of a real array :f:var:`func`
  !of length :code:`N`, in place, using the FFTPACK 5.1 routines :f:func_inline:`rfft1i` + 
  !:f:func_inline:`rfft1b`. These are called with the parameters
  !
  !* :code:`lensav` = :code:`N + int( log(dble(N))/log(2.d0) ) + 4`
  !* :code:`lenwrk` = :code:`N`
  !* :code:`lenr`  = :code:`N`
  !* :code:`inc` = :code:`1`
  !
  !The transform is not normalized, so applying it after :f:func_inline:`rfft_1d_forward`
  !returns the original array multiplied by :code:`N`. The input array is expected in the 
  !FFTPACK packed form produced by :f:func_inline:`rfft_1d_forward`, that is
  !
  !* :code:`[Re y(0), Re y(1), Im y(1), ..., Re y(N/2-1), Im y(N/2-1), Re y(N/2)]`   if `N` is even
  !* :code:`[Re y(0), Re y(1), Im y(1), ..., Re y((N-1)/2), Im y((N-1)/2)]`           if `N` is odd
  !
  !and the subroutine modifies the input array to contain the real sequence
  !
  !:math:`x(j) = \sum_{k=0}^{N-1} y(k) \cdot e^{i 2\pi/N \cdot j \cdot k}`
  !
  !for :math:`j \in [0,N-1]`, where :math:`y(k)` is the Hermitian sequence encoded by the packed input.
  !
  !The program stops if the FFTPACK routines return an error code.
  !
    real(8),dimension(:),intent(inout) :: func  !Input/output array of size :code:`N`, in the FFTPACK packed form on input,
                                                ! overwritten by the real sequence
    real(8),dimension(:),allocatable   :: wsave,work
    integer                            :: N,lenwrk,lensav,lenr,inc,ier
    N      = size(func)
    lenwrk = N
    lensav = N + int( log(dble(N))/log(2.d0) ) + 4
    lenr   = N
    inc    = 1
    allocate(wsave(lensav))
    call rfft1i(N,wsave,lensav,ier)
    if(ier==2)stop "rfft_1d_backward: LENSAV not big enough"
    allocate(work(lenwrk))
    call rfft1b(N,inc,func,lenr,wsave,lensav,work,lenwrk,ier)
    deallocate(wsave,work)
    select case(ier)
    case (0)
       return
    case (1)
       stop "rfft_1d_backward: LENR not big enough"
    case (2)
       stop "rfft_1d_backward: LENSAV not big enough"
    case (3)
       stop "rfft_1d_backward: LENWRK not big enough"
    case (20)
       stop "rfft_1d_backward: input error returned by lower level routine"
    end select
  end subroutine rfft_1d_backward
  !
  subroutine cfft_1d_backward(func)
  !This subroutine evaluates the backward 1-dimensional FFT of a complex array :f:var:`func`
  !of length :code:`N`, in place, using the FFTPACK 5.1 routines :f:func_inline:`cfft1i` + 
  !:f:func_inline:`cfft1b`. These are called with the parameters
  !
  !* :code:`lensav` = :code:`2*N + int( log(dble(N))/log(2.d0) ) + 4`
  !* :code:`lenwrk` = :code:`2*N`
  !* :code:`lenc`  = :code:`N`
  !* :code:`inc` = :code:`1`
  !
  !The transform is not normalized, so applying it after :f:func_inline:`cfft_1d_forward`
  !returns the original array multiplied by :code:`N`. The input array is expected in the 
  !aliased order produced by :f:func_inline:`cfft_1d_forward`, that is
  !
  !* :code:`[y(0), y(1), ..., y(N/2),     y(-N/2+1), ...,   y(-1)]`   if `N` is even
  !* :code:`[y(0), y(1), ..., y((N-1)/2), y(-(N-1)/2), ..., y(-1)]`   if `N` is odd
  !
  !and the subroutine modifies the input array to contain
  !
  !:math:`x(j) = \sum_{k=0}^{N-1} y(k) \cdot e^{i 2\pi/N \cdot j \cdot k}`
  !
  !for :math:`j \in [0,N-1]`.
  !
  !The program stops if the FFTPACK routines return an error code.
  !
    complex(8),dimension(:),intent(inout) :: func  !Input/output array of size :code:`N`, overwritten by its backward transform
    real(8),dimension(:),allocatable      :: wsave,work
    integer                               :: N,lenwrk,lensav,lenc,inc,ier
    N      = size(func)
    lenwrk = 2*N
    lensav = 2*N + int( log(dble(N))/log(2.d0) ) + 4
    lenc   = N
    inc    = 1
    allocate(wsave(lensav))
    call cfft1i(N,wsave,lensav,ier)
    if(ier==2)stop "cfft_1d_backward: LENSAV not big enough"
    allocate(work(lenwrk))
    call cfft1b(N,inc,func,lenc,wsave,lensav,work,lenwrk,ier)
    deallocate(wsave,work)
    select case(ier)
    case (0)
       return
    case (1)
       stop "cfft_1d_backward: LENC not big enough"
    case (2)
       stop "cfft_1d_backward: LENSAV not big enough"
    case (3)
       stop "cfft_1d_backward: LENWRK not big enough"
    case (20)
       stop "cfft_1d_backward: input error returned by lower level routine"
    end select
  end subroutine cfft_1d_backward


  !+-------------------------------------------------------------------+
  !PURPOSE  : Evaluate backward 2-DIM FFT using FFTPACK 5.1 routines.
  ! For purposes of exposition, assume the index ranges of array C are defined by
  !  C(0:L-1,0:M-1).
  ! For I=0,...,L-1 and J=0,...,M-1, the C(I,J)'s are given in the traditional 
  ! aliased form by
  !           L-1  M-1
  !  C(I,J) = SUM  SUM  C(L1,M1)*EXP(XI*2*PI*(I*L1/L + J*M1/M))
  !           L1=0 M1=0
  !+-------------------------------------------------------------------+
  subroutine rfft_2d_backward(func)
  !This subroutine evaluates the backward 2-dimensional FFT of a real array :f:var:`func`
  !of size :code:`(L,M)`, in place, using the FFTPACK 5.1 routines :f:func_inline:`rfft2i` + 
  !:f:func_inline:`rfft2b`. These are called with the parameters
  !
  !* :code:`ldim` = :code:`L`
  !* :code:`lensav` = :code:`L + 3*M + int( log(dble(L))/log(2.d0) ) + 2*int( log(dble(M))/log(2.d0) ) + 12`
  !* :code:`lenwrk` = :code:`M*(L+1)`
  !
  !The input array is expected in the FFTPACK packed form produced by 
  !:f:func_inline:`rfft_2d_forward`, and the subroutine modifies it to contain the real array
  !
  !:math:`x(i,j) = \sum_{l=0}^{L-1} \sum_{m=0}^{M-1} y(l,m) \cdot e^{i 2\pi (i l/L + j m/M)}`
  !
  !for :math:`i \in [0,L-1]`, :math:`j \in [0,M-1]`, where :math:`y(l,m)` is the Hermitian 
  !array encoded by the packed input. The transform carries no normalization factor, 
  !because the :math:`1/(LM)` factor is carried by :f:func_inline:`rfft_2d_forward`.
  !
  !The program stops if the FFTPACK routines return an error code.
  !
    real(8),dimension(:,:),intent(inout) :: func  !Input/output array of size :code:`(L,M)`, in the FFTPACK packed form on
                                                  ! input, overwritten by the real array
    real(8),dimension(:),allocatable   :: wsave,work
    integer                            :: L,M,lenwrk,lensav,ldim,ier
    L      = size(func,1)
    M      = size(func,2)
    lenwrk = M*(L+1)
    lensav = L+3*M  + int(log(dble(L))/log(2.d0))+2*int(log(dble(M))/log(2.d0)) + 12
    ldim   = L
    allocate(wsave(lensav))
    allocate(work(lenwrk))
    call rfft2i(L,M,wsave,lensav,ier)
    if(ier==2)stop "rfft_2d_backward: LENSAV not big enough"
    if(ier==20)stop "rfft_2d_backward: input error returned by lower level routine"
    call rfft2b(ldim,L,M,func,wsave,lensav,work,lenwrk,ier)
    deallocate(wsave,work)
    select case(ier)
    case (0)
       return
    case (6)
       stop "rfft_2d_backward: LDIM is less than 2*INT((L+1)/2)"
    case (2)
       stop "rfft_2d_backward: LENSAV not big enough"
    case (3)
       stop "rfft_2d_backward: LENWRK not big enough"
    case (20)
       stop "rfft_2d_backward: input error returned by lower level routine"
    end select
  end subroutine rfft_2d_backward
  !
  subroutine cfft_2d_backward(func)
  !This subroutine evaluates the backward 2-dimensional FFT of a complex array :f:var:`func`
  !of size :code:`(L,M)`, in place, using the FFTPACK 5.1 routines :f:func_inline:`cfft2i` + 
  !:f:func_inline:`cfft2b`. These are called with the parameters
  !
  !* :code:`ldim` = :code:`L`
  !* :code:`lensav` = :code:`2*(L+M) + int( log(dble(L))/log(2.d0) ) + int( log(dble(M))/log(2.d0) ) + 8`
  !* :code:`lenwrk` = :code:`2*L*M`
  !
  !The input array is expected in the form produced by :f:func_inline:`cfft_2d_forward`,
  !and the subroutine modifies it to contain
  !
  !:math:`x(i,j) = \sum_{l=0}^{L-1} \sum_{m=0}^{M-1} y(l,m) \cdot e^{i 2\pi (i l/L + j m/M)}`
  !
  !for :math:`i \in [0,L-1]`, :math:`j \in [0,M-1]`. The transform carries no normalization 
  !factor, because the :math:`1/(LM)` factor is carried by :f:func_inline:`cfft_2d_forward`.
  !
  !The program stops if the FFTPACK routines return an error code.
  !
    complex(8),dimension(:,:),intent(inout) :: func  !Input/output array of size :code:`(L,M)`, overwritten by its backward
                                                     ! transform
    real(8),dimension(:),allocatable      :: wsave,work
    integer                               :: L,M,lenwrk,lensav,ldim,inc,ier
    L      = size(func,1)
    M      = size(func,2)
    lenwrk = 2*L*M
    lensav = 2*(L+M) + int(log(dble(L))/log(2.d0))+int(log(dble(M))/log(2.d0)) + 8
    ldim   = L
    allocate(wsave(lensav))
    allocate(work(lenwrk))
    call cfft2i(L,M,wsave,lensav,ier)
    if(ier==2)stop "cfft_2d_backward: LENSAV not big enough"
    if(ier==20)stop "cfft_2d_backward: input error returned by lower level routine"
    call cfft2b(ldim,L,M,func,wsave,lensav,work,lenwrk,ier)
    deallocate(wsave,work)
    select case(ier)
    case (0)
       return
    case (5)
       stop "cfft_2d_backward: L > LDIM"
    case (2)
       stop "cfft_2d_backward: LENSAV not big enough"
    case (3)
       stop "cfft_2d_backward: LENWRK not big enough"
    case (20)
       stop "cfft_2d_backward: input error returned by lower level routine"
    end select
  end subroutine cfft_2d_backward


  !+-------------------------------------------------------------------+
  !PURPOSE  : Evaluate backward N-DIM FFT using FFTPACK 5.1 routines.
  ! For index L*JUMP+J*INC+1 where J=0,...,N-1 and 
  ! L=0,...,LOT-1, (that is, for the Jth element of the Lth sequence),
  !                      N-1
  !  C(L*JUMP+J*INC+1) = SUM C(L*JUMP+K*INC+1)*EXP(XI*J*K*2*PI/N)
  !                      K=0
  !+-------------------------------------------------------------------+
  subroutine rfft_nd_backward(func,n,lot)
  !This subroutine evaluates :code:`lot` backward 1-dimensional FFTs of length :code:`N`
  !on the real sequences stored consecutively in the 1-dimensional array :f:var:`func`,
  !in place, using the FFTPACK 5.1 routines :f:func_inline:`rfftmi` + :f:func_inline:`rfftmb`.
  !These are called with the parameters
  !
  !* :code:`lenr` = :code:`N*lot`
  !* :code:`lensav` = :code:`N + int( log(dble(N))/log(2.d0) ) + 4`
  !* :code:`lenwrk` = :code:`N*lot`
  !* :code:`jump` = :code:`N`
  !* :code:`inc` = :code:`1`
  !
  !Each sequence is expected in the FFTPACK packed form produced by 
  !:f:func_inline:`rfft_nd_forward`. The transforms are not normalized, so applying 
  !them after :f:func_inline:`rfft_nd_forward` returns the original array multiplied 
  !by :code:`N`. For the :math:`l`-th sequence, :math:`l \in [0,lot-1]`, the subroutine 
  !modifies the input array to contain the real sequence
  !
  !:math:`x(l N + j) = \sum_{k=0}^{N-1} y(l N + k) \cdot e^{i 2\pi/N \cdot j \cdot k}`
  !
  !for :math:`j \in [0,N-1]`, where :math:`y(l N + k)` is the Hermitian sequence encoded 
  !by the packed input.
  !
  !The program stops if :code:`mod(N*lot,size(func))/=0`, or if the FFTPACK routines 
  !return an error code.
  !
    real(8),dimension(:),intent(inout) :: func  !Input/output array storing :code:`lot` sequences of length :code:`N`
                                                ! consecutively, in the FFTPACK packed form on input, overwritten by the real
                                                ! sequences
    integer,intent(in)                 :: N  !Length of each sequence
    integer,intent(in)                 :: lot  !Number of sequences
    real(8),dimension(:),allocatable   :: wsave,work
    integer                            :: L,lenwrk,lensav,lenr,ier,inc,jump
    L=size(func)
    if(mod(N*lot,L)/=0)stop "rfft_Nd_backward: incommensurate values of parameters." 
    !
    lenr   = N*lot
    lenwrk = N*lot
    lensav = N  + int(log(dble(N))/log(2.d0)) + 4
    allocate(wsave(lensav))
    allocate(work(lenwrk))
    call rfftmi(N,wsave,lensav,ier)
    if(ier==2)stop "rfft_Nd_backward: LENSAV not big enough"
    jump   = N
    inc    = 1
    call rfftmb(lot,jump,N,inc,func,lenr,wsave,lensav,work,lenwrk,ier)
    deallocate(wsave,work)
    select case(ier)
    case (0)
       return
    case (1)
       stop "rfft_Nd_backward: LENR   not big enough"
    case (2)
       stop "rfft_Nd_backward: LENSAV not big enough"
    case (3)
       stop "rfft_Nd_backward: LENWRK not big enough"
    case (4)
       stop "rfft_Nd_backward: INC,JUMP,N,LOT are not consistent"
    end select
  end subroutine rfft_nd_backward
  !
  subroutine cfft_nd_backward(func,n,lot)
  !This subroutine evaluates :code:`lot` backward 1-dimensional FFTs of length :code:`N`
  !on the complex sequences stored consecutively in the 1-dimensional array :f:var:`func`,
  !in place, using the FFTPACK 5.1 routines :f:func_inline:`cfftmi` + :f:func_inline:`cfftmb`.
  !These are called with the parameters
  !
  !* :code:`lenc` = :code:`N*lot`
  !* :code:`lensav` = :code:`2*N + int( log(dble(N))/log(2.d0) ) + 4`
  !* :code:`lenwrk` = :code:`2*N*lot`
  !* :code:`jump` = :code:`N`
  !* :code:`inc` = :code:`1`
  !
  !Each sequence is expected in the aliased order produced by 
  !:f:func_inline:`cfft_nd_forward`. The transforms are not normalized, so applying 
  !them after :f:func_inline:`cfft_nd_forward` returns the original array multiplied 
  !by :code:`N`. For the :math:`l`-th sequence, :math:`l \in [0,lot-1]`, the subroutine 
  !modifies the input array to contain
  !
  !:math:`x(l N + j) = \sum_{k=0}^{N-1} y(l N + k) \cdot e^{i 2\pi/N \cdot j \cdot k}`
  !
  !for :math:`j \in [0,N-1]`.
  !
  !The program stops if :code:`mod(N*lot,size(func))/=0`, or if the FFTPACK routines 
  !return an error code.
  !
    complex(8),dimension(:),intent(inout) :: func  !Input/output array storing :code:`lot` sequences of length :code:`N`
                                                   ! consecutively, overwritten by their backward transforms
    integer,intent(in)                    :: n  !Length of each sequence
    integer,intent(in)                    :: lot  !Number of sequences
    real(8),dimension(:),allocatable      :: wsave,work
    integer                               :: L,lenwrk,lensav,lenc,ier,inc,jump
    L=size(func)
    if(mod(N*lot,L)/=0)stop "cfft_Nd_backward: incommensurate values of parameters." 
    !
    lenc   = N*lot
    lenwrk = 2*N*lot
    lensav = 2*N + int(log(dble(N))/log(2.d0)) + 4
    allocate(wsave(lensav))
    allocate(work(lenwrk))
    call cfftmi(N,wsave,lensav,ier)
    if(ier==2)stop "cfft_Nd_backward: LENSAV not big enough"
    !
    jump   = N
    inc    = 1
    call cfftmb(lot,jump,N,inc,func,lenc,wsave,lensav,work,lenwrk,ier)
    deallocate(wsave,work)
    select case(ier)
    case (0)
       return
    case (1)
       stop "cfft_Nd_backward: LENR   not big enough"
    case (2)
       stop "cfft_Nd_backward: LENSAV not big enough"
    case (3)
       stop "cfft_Nd_backward: LENWRK not big enough"
    case (4)
       stop "cfft_Nd_backward: INC,JUMP,N,LOT are not consistent"
    end select
  end subroutine cfft_nd_backward


  !+-------------------------------------------------------------------+
  !PURPOSE  : Evaluate backward 1-DIM COS-FFT using FFTPACK 5.1 routines.
  !+-------------------------------------------------------------------+
  subroutine cost_1d_backward(func)
  !This subroutine evaluates the backward 1-dimensional cosine FFT of a real array :f:var:`func`
  !of length :code:`N`, in place, using the FFTPACK 5.1 routines :f:func_inline:`cost1i` + 
  !:f:func_inline:`cost1b`. These are called with the parameters
  !
  !* :code:`lensav` = :code:`2*N + int( log(dble(N))/log(2.d0) ) + 4`
  !* :code:`lenwrk` = :code:`N-1`
  !* :code:`lenr`  = :code:`N`
  !* :code:`inc` = :code:`1`
  !
  !The transform is the inverse of :f:func_inline:`cost_1d_forward`, i.e. the 
  !normalization is carried by the forward transform, so no scaling is applied here. 
  !The input array contains the coefficients :math:`y(k)` produced by 
  !:f:func_inline:`cost_1d_forward`, and the subroutine modifies it to contain
  !
  !:math:`x(i) = y(1) + (-1)^{i-1} y(N) + 2 \sum_{k=2}^{N-1} y(k) \cdot \cos\left(\frac{\pi (k-1)(i-1)}{N-1}\right)`
  !
  !for :math:`i \in [1,N]`.
  !
  !The program stops if the FFTPACK routines return an error code.
  !
    real(8),dimension(:),intent(inout) :: func  !Input/output array of size :code:`N`, overwritten by its backward cosine
                                                ! transform
    real(8),dimension(:),allocatable   :: wsave,work
    integer                            :: N,lenwrk,lensav,lenr,inc,ier
    N      = size(func)
    lensav = 2*N + int( log(dble(N))/log(2.d0) ) + 4
    lenwrk = N-1
    lenr   = N
    inc    = 1
    allocate(wsave(lensav))
    allocate(work(lenwrk))
    call cost1i(N,wsave,lensav,ier)
    if(ier==2)stop "cost_1d_backward: LENSAV not big enough"
    if(ier==20)stop "cost_1d_backward: error returned by lower level routine"
    call cost1b(N,inc,func,lenr,wsave,lensav,work,lenwrk,ier)
    deallocate(wsave,work)
    select case(ier)
    case (0)
       return
    case (1)
       stop "cost_1d_backward: LENR not big enough"
    case (2)
       stop "cost_1d_backward: LENSAV not big enough"
    case (3)
       stop "cost_1d_backward: LENWRK not big enough"
    case (20)
       stop "cost_1d_backward: input error returned by lower level routine"
    end select
  end subroutine cost_1d_backward

  !+-------------------------------------------------------------------+
  !PURPOSE  : Evaluate backward 1-DIM SIN-FFT using FFTPACK 5.1 routines.
  !+-------------------------------------------------------------------+
  subroutine sint_1d_backward(func)
  !This subroutine evaluates the backward 1-dimensional sine FFT of a real array :f:var:`func`
  !of length :code:`N`, in place, using the FFTPACK 5.1 routines :f:func_inline:`sint1i` + 
  !:f:func_inline:`sint1b`. These are called with the parameters
  !
  !* :code:`lensav` = :code:`N/2 + N + int( log(dble(N))/log(2.d0) ) + 4`
  !* :code:`lenwrk` = :code:`2*(N+1)`
  !* :code:`lenr`  = :code:`N`
  !* :code:`inc` = :code:`1`
  !
  !The transform is the inverse of :f:func_inline:`sint_1d_forward`, i.e. the 
  !normalization is carried by the forward transform, so no scaling is applied here. 
  !The input array contains the coefficients :math:`y(k)` produced by 
  !:f:func_inline:`sint_1d_forward`, and the subroutine modifies it to contain
  !
  !:math:`x(i) = 2 \sum_{k=1}^{N} y(k) \cdot \sin\left(\frac{\pi k i}{N+1}\right)`
  !
  !for :math:`i \in [1,N]`.
  !
  !The program stops if the FFTPACK routines return an error code.
  !
    real(8),dimension(:),intent(inout) :: func  !Input/output array of size :code:`N`, overwritten by its backward sine
                                                ! transform
    real(8),dimension(:),allocatable   :: wsave,work
    integer                            :: N,lenwrk,lensav,lenr,inc,ier
    N      = size(func)
    lensav = N/2 + N + int( log(dble(N))/log(2.d0) ) + 4
    lenwrk = 2*(N+1)
    lenr   = N
    inc    = 1
    allocate(wsave(lensav))
    allocate(work(lenwrk))
    call sint1i(N,wsave,lensav,ier)
    if(ier==2)stop "sint_1d_backward: LENSAV not big enough"
    if(ier==20)stop "sint_1d_backward: error returned by lower level routine"
    call sint1b(N,inc,func,lenr,wsave,lensav,work,lenwrk,ier)
    deallocate(wsave,work)
    select case(ier)
    case (0)
       return
    case (1)
       stop "sint_1d_backward: LENR not big enough"
    case (2)
       stop "sint_1d_backward: LENSAV not big enough"
    case (3)
       stop "sint_1d_backward: LENWRK not big enough"
    case (20)
       stop "sint_1d_backward: input error returned by lower level routine"
    end select
  end subroutine sint_1d_backward

  !+-------------------------------------------------------------------+
  !PURPOSE  : Evaluate backward M-DIM COS-FFT using FFTPACK 5.1 routines.
  !+-------------------------------------------------------------------+
  subroutine cost_nd_backward(func,n,lot)
  !This subroutine evaluates :code:`lot` backward 1-dimensional cosine FFTs of length :code:`N`
  !on the real sequences stored consecutively in the 1-dimensional array :f:var:`func`,
  !in place, using the FFTPACK 5.1 routines :f:func_inline:`costmi` + :f:func_inline:`costmb`.
  !These are called with the parameters
  !
  !* :code:`lenr` = :code:`N*lot`
  !* :code:`lensav` = :code:`2*N + int( log(dble(N))/log(2.d0) ) + 4`
  !* :code:`lenwrk` = :code:`lot*(N+1)`
  !* :code:`jump` = :code:`N`
  !* :code:`inc` = :code:`1`
  !
  !For the :math:`l`-th sequence, :math:`l \in [0,lot-1]`, the subroutine acts as
  !:f:func_inline:`cost_1d_backward` on the elements :math:`l N + 1, ..., (l+1) N` of 
  !the array, so the transform is the inverse of :f:func_inline:`cost_nd_forward` 
  !and the normalization is carried by the forward transform, so no scaling is applied here.
  !
  !The program stops if :code:`mod(N*lot,size(func))/=0`, or if the FFTPACK routines 
  !return an error code.
  !
    real(8),dimension(:),intent(inout) :: func  !Input/output array storing :code:`lot` sequences of length :code:`N`
                                                ! consecutively, overwritten by their backward cosine transforms
    integer,intent(in)                 :: N  !Length of each sequence
    integer,intent(in)                 :: lot  !Number of sequences
    real(8),dimension(:),allocatable   :: wsave,work
    integer                            :: L,lenwrk,lensav,lenr,inc,ier,jump
    L=size(func)
    if(mod(N*lot,L)/=0)stop "cost_Nd_backward: incommensurate values of parameters." 
    !
    lenr   = N*lot
    lenwrk = lot*(N+1)
    lensav = 2*N  + int(log(dble(N))/log(2.d0)) + 4
    jump   = N
    inc    = 1
    allocate(wsave(lensav),work(lenwrk))
    call costmi(N,wsave,lensav,ier)
    if(ier==2)stop "cost_Nd_backward: LENSAV not big enough"
    if(ier==20)stop "cost_Nd_backward: error returned by lower level routine"
    call costmb(lot,jump,N,inc,func,lenr,wsave,lensav,work,lenwrk,ier)
    deallocate(wsave,work)
    select case(ier)
    case (0)
       return
    case (1)
       stop "cost_Nd_backward: LENR not big enough"
    case (2)
       stop "cost_Nd_backward: LENSAV not big enough"
    case (3)
       stop "cost_Nd_backward: LENWRK not big enough"
    case (4)
       stop "cost_Nd_backward: INC,JUMP,N,LOT are not consistent" 
    case (20)
       stop "cost_Nd_backward: input error returned by lower level routine"
    end select
  end subroutine cost_nd_backward

  !+-------------------------------------------------------------------+
  !PURPOSE  : Evaluate backward M-DIM SIN-FFT using FFTPACK 5.1 routines.
  !+-------------------------------------------------------------------+
  subroutine sint_nd_backward(func,n,lot)
  !This subroutine evaluates :code:`lot` backward 1-dimensional sine FFTs of length :code:`N`
  !on the real sequences stored consecutively in the 1-dimensional array :f:var:`func`,
  !in place, using the FFTPACK 5.1 routines :f:func_inline:`sintmi` + :f:func_inline:`sintmb`.
  !These are called with the parameters
  !
  !* :code:`lenr` = :code:`N*lot`
  !* :code:`lensav` = :code:`N/2 + N + int( log(dble(N))/log(2.d0) ) + 4`
  !* :code:`lenwrk` = :code:`2*lot*(N+2)`
  !* :code:`jump` = :code:`N`
  !* :code:`inc` = :code:`1`
  !
  !For the :math:`l`-th sequence, :math:`l \in [0,lot-1]`, the subroutine acts as
  !:f:func_inline:`sint_1d_backward` on the elements :math:`l N + 1, ..., (l+1) N` of 
  !the array, so the transform is the inverse of :f:func_inline:`sint_nd_forward` 
  !and the normalization is carried by the forward transform, so no scaling is applied here.
  !
  !The program stops if :code:`mod(N*lot,size(func))/=0`, or if the FFTPACK routines 
  !return an error code.
  !
    real(8),dimension(:),intent(inout) :: func  !Input/output array storing :code:`lot` sequences of length :code:`N`
                                                ! consecutively, overwritten by their backward sine transforms
    integer,intent(in)                 :: N  !Length of each sequence
    integer,intent(in)                 :: lot  !Number of sequences
    real(8),dimension(:),allocatable   :: wsave,work
    integer                            :: L,lenwrk,lensav,lenr,inc,ier,jump
    L=size(func)
    if(mod(N*lot,L)/=0)stop "cost_Nd_backward: incommensurate values of parameters."
    !
    lenr   = N*lot
    lenwrk = 2*lot*(N+2)
    lensav = N/2 + N  + int(log(dble(N))/log(2.d0)) + 4
    jump   = N
    inc    = 1
    allocate(wsave(lensav),work(lenwrk))
    call sintmi(N,wsave,lensav,ier)
    if(ier==2)stop "cost_Nd_backward: LENSAV not big enough"
    if(ier==20)stop "cost_Nd_backward: error returned by lower level routine"
    call sintmb(lot,jump,N,inc,func,lenr,wsave,lensav,work,lenwrk,ier)
    deallocate(wsave,work)
    select case(ier)
    case (0)
       return
    case (1)
       stop "cost_Nd_backward: LENR not big enough"
    case (2)
       stop "cost_Nd_backward: LENSAV not big enough"
    case (3)
       stop "cost_Nd_backward: LENWRK not big enough"
    case (4)
       stop "cost_Nd_backward: INC,JUMP,N,LOT are not consistent" 
    case (20)
       stop "cost_Nd_backward: input error returned by lower level routine"
    end select
  end subroutine sint_nd_backward














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
  !This function shifts the zero-frequency component of a real array :f:var:`fin`
  !of length :code:`N` to the center of the array, and returns the result as an array 
  !of the same size. It is a circular rotation by :code:`floor((N+1)/2)` elements, 
  !so the input
  !
  !* :code:`[y(0), y(1), ..., y(N/2),     y(-N/2+1), ...,   y(-1)]`   if `N` is even
  !* :code:`[y(0), y(1), ..., y((N-1)/2), y(-(N-1)/2), ..., y(-1)]`   if `N` is odd
  !
  !is returned as
  !
  !* :code:`[y(-N/2), ..., y(-1), y(0), y(1), ..., y(N/2-1)]`          if `N` is even
  !* :code:`[y(-(N-1)/2), ..., y(-1), y(0), y(1), ..., y((N-1)/2)]`    if `N` is odd
  !
  !where :code:`y(N/2)` and :code:`y(-N/2)` are the same element for even :code:`N`.
  !
    real(8),dimension(:)         :: fin  !Array to shift, of size :code:`N`, in the aliased order returned by
                                         ! :f:func_inline:`fft`
    real(8),dimension(size(fin)) :: fout  !Shifted array, of the same size as :f:var:`fin`, with the zero-frequency component
                                          ! at the center
    integer                      :: L,p2
    L  = size(fin)
    p2 = floor(dble(L+1)/2.d0)
    fout = [fin(p2+1:),fin(1:p2)]
  end function rfft_1d_shift
  !
  function cfft_1d_shift(fin) result(fout)
  !This function shifts the zero-frequency component of a complex array :f:var:`fin`
  !of length :code:`N` to the center of the array, and returns the result as an array 
  !of the same size. It acts as :f:func_inline:`rfft_1d_shift` on complex data.
  !
    complex(8),dimension(:)          :: fin  !Array to shift, of size :code:`N`, in the aliased order returned by
                                             ! :f:func_inline:`fft`
    complex(8),dimension(size(fin))  :: fout  !Shifted array, of the same size as :f:var:`fin`, with the zero-frequency
                                              ! component at the center
    integer                          :: L,p2
    L  = size(fin)
    p2 = floor(dble(L+1)/2.d0)
    fout = [fin(p2+1:),fin(1:p2)]
  end function cfft_1d_shift


  !+-------------------------------------------------------------------+
  !PURPOSE  : The inverse of _shift
  !+-------------------------------------------------------------------+
  function rfft_1d_ishift(fin) result(fout)
  !This function is the inverse of :f:func_inline:`rfft_1d_shift`. It moves the
  !zero-frequency component of a real array :f:var:`fin` of length :code:`N` from the 
  !center of the array back to the first position, and returns the result as an array 
  !of the same size. It is a circular rotation by :code:`N-floor((N+1)/2)` elements, 
  !so the input
  !
  !* :code:`[y(-N/2), ..., y(-1), y(0), y(1), ..., y(N/2-1)]`          if `N` is even
  !* :code:`[y(-(N-1)/2), ..., y(-1), y(0), y(1), ..., y((N-1)/2)]`    if `N` is odd
  !
  !is returned as
  !
  !* :code:`[y(0), y(1), ..., y(N/2-1), y(-N/2), ..., y(-1)]`          if `N` is even
  !* :code:`[y(0), y(1), ..., y((N-1)/2), y(-(N-1)/2), ..., y(-1)]`    if `N` is odd
  !
  !that is, the aliased order returned by :f:func_inline:`fft`.
  !
    real(8),dimension(:)         :: fin  !Array of size :code:`N` with the zero-frequency component at the center
    real(8),dimension(size(fin)) :: fout  !Array in the aliased order returned by :f:func_inline:`fft`, of the same size as
                                          ! :f:var:`fin`
    integer                      :: L,p2
    L  = size(fin)
    p2 = L-floor(dble(L+1)/2.d0)
    fout = [fin(p2+1:),fin(1:p2)]
  end function rfft_1d_ishift
  !
  function cfft_1d_ishift(fin) result(fout)
  !This function is the inverse of :f:func_inline:`cfft_1d_shift`. It moves the
  !zero-frequency component of a complex array :f:var:`fin` of length :code:`N` from the 
  !center of the array back to the first position, and returns the result as an array 
  !of the same size. It acts as :f:func_inline:`rfft_1d_ishift` on complex data.
  !
    complex(8),dimension(:)          :: fin  !Array of size :code:`N` with the zero-frequency component at the center
    complex(8),dimension(size(fin))  :: fout  !Array in the aliased order returned by :f:func_inline:`fft`, of the same size as
                                              ! :f:var:`fin`
    integer                          :: L,p2
    L  = size(fin)
    p2 = L-floor(dble(L+1)/2.d0)
    fout = [fin(p2+1:),fin(1:p2)]
  end function cfft_1d_ishift


  !+-------------------------------------------------------------------+
  !PURPOSE  : Alternate sign of the output of a FFT:
  !+-------------------------------------------------------------------+
  subroutine rfft_1d_ex(func)
  !This subroutine alternates the sign of the elements of a real array :f:var:`func`
  !in place,
  !
  !:math:`f(i) \to (-1)^{i-1} f(i)`
  !
  !for :math:`i \in [1,N]`, where :code:`N=size(func)`. Applying it twice returns 
  !the original array.
  !
    real(8),dimension(:)    :: func  !Input/output array, overwritten by :math:`(-1)^{i-1} f(i)`
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
  !This subroutine alternates the sign of the elements of a complex array :f:var:`func`
  !in place,
  !
  !:math:`f(i) \to (-1)^{i-1} f(i)`
  !
  !for :math:`i \in [1,N]`, where :code:`N=size(func)`. Applying it twice returns 
  !the original array.
  !
    complex(8),dimension(:) :: func  !Input/output array, overwritten by :math:`(-1)^{i-1} f(i)`
    real(8)                 :: ex
    integer                 :: i
    ex=-1.d0
    do i=1,size(func)
       ex=-ex
       func(i)=ex*func(i)
    enddo
  end subroutine cfft_1d_ex


END MODULE FFT_FFTPACK
