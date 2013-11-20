!###############################################################
! PROGRAM  : DERIVATE
! PURPOSE  : A set of routines to perform specific integrals
!###############################################################
MODULE DERIVATE
  implicit none
  private

  complex(8),parameter :: one= (1.d0,0.d0)
  complex(8),parameter :: zero=(0.d0,0.d0)
  complex(8),parameter :: xi = (0.d0,1.d0)
  real(8),parameter    :: pi=3.14159265358979323846264338327950288419716939937510D0

  interface djacobian
     module procedure fdjac_nn_func, fdjac_nn_sub, fdjac_mn_func,fdjac_mn_sub
  end interface djacobian

  interface dgradient
     module procedure fdjac_1n_func, fdjac_1n_sub
  end interface dgradient

  interface f_djacobian
     module procedure f_jac_nn_func, f_jac_nn_sub,  f_jac_mn_func, f_jac_mn_sub
  end interface f_djacobian

  interface f_dgradient
     module procedure f_jac_1n_func, f_jac_1n_sub
  end interface f_dgradient


  interface cjacobian
     module procedure c_fdjac_nn_func, c_fdjac_nn_sub, c_fdjac_mn_func, c_fdjac_mn_sub
  end interface cjacobian

  interface cgradient
     module procedure c_fdjac_1n_func,c_fdjac_1n_sub
  end interface cgradient

  interface f_cjacobian
     module procedure c_f_jac_nn_func , c_f_jac_nn_sub , c_f_jac_mn_func , c_f_jac_mn_sub
  end interface f_cjacobian

  interface f_cgradient
     module procedure c_f_jac_1n_func , c_f_jac_1n_sub
  end interface f_cgradient


  public :: djacobian
  public :: dgradient
  public :: f_djacobian
  public :: f_dgradient
  !
  public :: cjacobian
  public :: cgradient
  public :: f_cjacobian
  public :: f_cgradient
  !
  public :: deriv


contains

  function deriv(f,dh) result(df)
    real(8),dimension(:),intent(in) :: f
    real(8),intent(in)              :: dh
    real(8),dimension(size(f))      :: df
    integer                         :: i,L
    L=size(f)
    df(1)= (f(2)-f(1))/dh
    do i=2,L-1
       df(i) = (f(i+1)-f(i-1))/(2.d0*dh)
    enddo
    df(L)= (f(L)-f(L-1))/dh
  end function deriv



  !------------------------------------------------------------------------------
  ! Purpose: estimates an N/M/1 by N jacobian matrix using forward differences.
  !   computes a forward-difference approximation
  !   to the N by N jacobian matrix associated with a specified
  !   problem of N functions in N variables. If the jacobian has
  !   a banded form, then function evaluations are saved by only
  !   approximating the nonzero terms.
  ! Arguments:
  !    -Input FCN: the name of the user-supplied function which
  !    calculates the functions. The routines accept functions/routines.
  !    Check out their explicit form in the interfaces of the subroutines
  !    -Input, real(8) X(:) the point where the jacobian is evaluated.
  !    -Output, real(8) FJAC(N,N) the approximate jacobian matrix.
  !    -Input,optional integer ML, MU, specify the number of subdiagonals and
  !    superdiagonals within the band of the jacobian matrix. Default values N-1
  !    -Input,optional real(8) EPSFCN, is used in determining a suitable step
  !    length for the forward-difference approximation.  This approximation
  !    assumes that the relative errors in the functions are of the order of
  !    EPSFCN.  If EPSFCN is less than the machine precision, it is assumed that
  !    the relative errors in the functions are of the order of the machine
  !    precision.
  !------------------------------------------------------------------------------

  !  DOUBLE PRECISION
  !------------------------
  include 'derivate_fjacobian_d.f90'

  !  DOUBLE COMPLEX
  !------------------------
  include 'derivate_fjacobian_c.f90'


END MODULE DERIVATE
