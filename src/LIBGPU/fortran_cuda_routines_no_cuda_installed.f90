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

 contains

      !------------------------------------------!
 subroutine cula__cgetrf_(n,a,ipiv)
 implicit none
 integer    :: n
 complex(4) :: a(n,n)
 integer    :: ipiv(n),status
 end subroutine
       !------------------------------------------!
 subroutine cula__zgetrf_(n,a,ipiv)
 implicit none
 integer(4)       :: n
 complex(8)       :: a(n,n)
 integer(4)       :: ipiv(n),INFO,status
 end subroutine
       !------------------------------------------!
 subroutine cula__dgetrf_(n,a,ipiv)
 implicit none
 integer    :: n,status
 real(8)    :: a(n,n)
 integer(4) :: ipiv(n)
 end subroutine
       !------------------------------------------!
 subroutine cula__sgetrf_(n,a,ipiv)
 implicit none
 integer    :: n,status
 real(4)    :: a(n,n)
 integer(4) :: ipiv(n)
 end subroutine
       !------------------------------------------!
 subroutine cula__cgetri_(n,a,ipiv)
 implicit none
 integer    :: n,status
 complex(4) :: a(n,n)
 integer(4)  :: ipiv(n)
 end subroutine
       !------------------------------------------!
 subroutine cula__zgetri_(n,a,ipiv)
 implicit none
 integer          :: n,status
 complex(8)       :: a(n,n),work(2*n)
 integer(4)       :: ipiv(n)
 end subroutine
       !------------------------------------------!
 subroutine cula__dgetri_(n,a,ipiv)
 implicit none
 integer    :: n,status
 real(8)    :: a(n,n)
 integer    :: ipiv(n)
 end subroutine
       !------------------------------------------!
 subroutine cula__sgetri_(n,a,ipiv)
 implicit none
 integer    :: n,status
 real(4)    :: a(n,n)
 integer    :: ipiv(n)
 end subroutine
       !------------------------------------------!

!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************

  subroutine cula_svd_cs(M,N,A,values,U,V)
  implicit none
  integer          :: M,N
  complex(4)       :: A(M,N),U(M,M),V(N,N)
  real(4)          :: values(N)
  end subroutine

  subroutine cula_svd_c(M,N,A,values,U,V)
  implicit none
  integer          :: M,N
  complex(8)       :: A(M,N),U(M,M),V(N,N)
  real(8)          :: values(N)
  end subroutine

  subroutine cula_svd_d(M,N,A,values,U,V)
  implicit none
  integer          :: M,N
  real(8)          :: A(M,N),U(M,M),V(N,N)
  real(8)          :: values(N)
  end subroutine

  subroutine cula_svd_s(M,N,A,values,U,V)
  implicit none
  integer          :: M,N
  real(4)          :: A(M,N),U(M,M),V(N,N)
  real(4)          :: values(N)
  end subroutine


  subroutine cula_eigenvector_square_cs(lsize,A,values,copyA)
  implicit none
  integer             :: lsize
  complex(4)          :: A(lsize,lsize)
  real(4)             :: values(lsize)
  complex(4)          :: copyA(lsize,lsize)
  end subroutine

  subroutine cula_eigenvector_square_c(lsize,A,values,copyA)
  implicit none
  integer             :: lsize
  complex(8)          :: A(lsize,lsize)
  real(8)             :: values(lsize)
  complex(8)          :: copyA(lsize,lsize)
  end subroutine

  subroutine cula_eigenvector_square_d(lsize,A,values,copyA)
  implicit none
  integer          :: lsize
  real(8)          :: A(lsize,lsize),values(lsize)
  real(8)          :: copyA(lsize,lsize)
  end subroutine

  subroutine cula_eigenvector_square_r(lsize,A,values,copyA)
  implicit none
  integer          :: lsize
  real(4)          :: A(lsize,lsize),values(lsize)
  real(4)          :: copyA(lsize,lsize)
  end subroutine

  subroutine cula_sgerf_test
  end subroutine

  subroutine cula_CHECK_STATUS(STATUS)
  integer :: status
  end subroutine

  subroutine init_gpu_device
  end subroutine

!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************
!**************************************************************************

end module
     
