module fortran_cuda_routines


 INTERFACE cula__getrf
  MODULE PROCEDURE cula__cgetrf_,cula__zgetrf_,cula__sgetrf_,cula__dgetrf_ 
 END INTERFACE

 INTERFACE cula__getri
  MODULE PROCEDURE cula__cgetri_,cula__zgetri_,cula__sgetri_,cula__dgetri_
 END INTERFACE

 INTERFACE cula_eigenvector_square
  MODULE PROCEDURE cula_eigenvector_square_r,cula_eigenvector_square_d,cula_eigenvector_square_cs,cula_eigenvector_square_c
 END INTERFACE

 INTERFACE cula_svd
  MODULE PROCEDURE cula_svd_cs,cula_svd_c,cula_svd_d,cula_svd_s
 END INTERFACE

 real(8),private,parameter :: shift_diag_positive_def = 100.


contains

!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************

#include 'fortran_cuda_routines_cula.h'

!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************

        !---------------!

      subroutine init_gpu_device
      implicit none
      integer          :: dev,idevices,i,j,status
      character(20)    ::  n
      integer,external :: cula_initialize
      external         :: cula_shutdown
      interface 
       subroutine cuda_init_flags
       end subroutine
      end interface
      write(*,*) '====================================='
      write(*,*) '    --- Existing CUDA devices ---  '
      write(*,*) '====================================='
#ifdef _CULA   
      STATUS = CULA_INITIALIZE()
#endif
      call cuInit(0)
      idevices = 0
      call cuDeviceGetCount(idevices)
      do i=1,idevices
            call cuDeviceGet(dev, i-1)
            call cuDeviceGetName(n, 256, dev)
            write(*,*) '  Device  ' ,i-1,' : ',n
      enddo
      write(*,*) '====================================='
      call cuda_init_flags
      return
      end subroutine

        !---------------!

      subroutine finalize_gpu
      integer          :: status
      integer,external :: cula_shutdown
#ifdef _CULA
       STATUS = CULA_SHUTDOWN()
#endif
      end subroutine

        !---------------!

!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************

end module
