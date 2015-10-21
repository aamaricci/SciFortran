  include "mkl_dfti.f90"
  include "mkl_trig_transforms.f90"
  module FFT
    use MKL_DFTI
    use MKL_DFT_TYPE
    ! use MKL_TRIG_TRANSFORMS
    ! use INTERPOLATE
    implicit none 
    private
    public :: cfft_1d_forward
    public :: cfft_1d_backward
    public :: cfft_1d_shift
    public :: swap_fftrt2rw,cfft_1d_ex
    REAL(8),PARAMETER    :: PI    = 3.14159265358979323846264338327950288419716939937510D0

    integer                        :: status
    type(DFTI_DESCRIPTOR), pointer :: Handle


  contains



    !+-------------------------------------------------------------------+
    !PURPOSE  : Evaluate forward and backward FFT using MKL routines.
    !+-------------------------------------------------------------------+
    subroutine cfft_1d_forward(func)
      complex(8),dimension(:),intent(inout) :: func
      Status   = DftiCreateDescriptor(Handle,DFTI_DOUBLE,DFTI_COMPLEX,1,size(func))
      Status   = DftiCommitDescriptor(Handle)
      Status   = DftiComputeForward(Handle,func)
      Status   = DftiFreeDescriptor(Handle)
    end subroutine cfft_1d_forward

    subroutine cfft_1d_backward(func)
      complex(8),dimension(:),intent(inout) :: func
      Status = DftiCreateDescriptor(Handle,DFTI_DOUBLE,DFTI_COMPLEX,1,size(func))
      Status = DftiCommitDescriptor(Handle)
      Status = DftiComputeBackward(Handle,func)
      Status = DftiFreeDescriptor(Handle)
    end subroutine cfft_1d_backward

    function cfft_1d_shift(fin,L) result(fout)
      integer                   :: i,L
      complex(8),dimension(2*L) :: fin
      complex(8),dimension(-L:L):: fout,dout
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


  end module FFT
  
