MODULE SF_PAULI
  implicit none
  private

  complex(8),parameter :: zero=(0.d0,0.d0)
  complex(8),parameter :: xi=(0.d0,1.d0)
  complex(8),parameter :: one=(1.d0,0.d0)

  complex(8),dimension(2,2),parameter,public :: pauli_0=reshape([one,zero,zero,one],[2,2])
  complex(8),dimension(2,2),parameter,public :: pauli_x=reshape([zero,one,one,zero],[2,2])
  complex(8),dimension(2,2),parameter,public :: pauli_y=reshape([zero,xi,-xi,zero],[2,2])
  complex(8),dimension(2,2),parameter,public :: pauli_z=reshape([one,zero,zero,-one],[2,2])
  !
  complex(8),dimension(2,2),parameter,public :: pauli_1=pauli_x
  complex(8),dimension(2,2),parameter,public :: pauli_2=pauli_y
  complex(8),dimension(2,2),parameter,public :: pauli_3=pauli_z
  !
  complex(8),dimension(2,2),parameter,public :: pauli_tau_0=pauli_0
  complex(8),dimension(2,2),parameter,public :: pauli_tau_x=pauli_x
  complex(8),dimension(2,2),parameter,public :: pauli_tau_y=pauli_y
  complex(8),dimension(2,2),parameter,public :: pauli_tau_z=pauli_z
  !
  complex(8),dimension(2,2),parameter,public :: pauli_tau_1=pauli_x
  complex(8),dimension(2,2),parameter,public :: pauli_tau_2=pauli_y
  complex(8),dimension(2,2),parameter,public :: pauli_tau_3=pauli_z
  !
  complex(8),dimension(2,2),parameter,public :: pauli_sigma_0=pauli_0
  complex(8),dimension(2,2),parameter,public :: pauli_sigma_x=pauli_x
  complex(8),dimension(2,2),parameter,public :: pauli_sigma_y=pauli_y
  complex(8),dimension(2,2),parameter,public :: pauli_sigma_z=pauli_z
  !
  complex(8),dimension(2,2),parameter,public :: pauli_sigma_1=pauli_x
  complex(8),dimension(2,2),parameter,public :: pauli_sigma_2=pauli_y
  complex(8),dimension(2,2),parameter,public :: pauli_sigma_3=pauli_z
  !
  !
  complex(8),dimension(2,2),parameter,public :: pauli_sigma_plus =pauli_x+xi*pauli_y
  complex(8),dimension(2,2),parameter,public :: pauli_sigma_minus=pauli_x-xi*pauli_y
  !
  complex(8),dimension(2,2),parameter,public :: pauli_tau_plus =pauli_x+xi*pauli_y
  complex(8),dimension(2,2),parameter,public :: pauli_tau_minus=pauli_x-xi*pauli_y
  

  interface kron_pauli
     module procedure kronecker_product_pauli_matrices
  end interface kron_pauli
  interface kron_sigma
     module procedure kronecker_product_pauli_matrices
  end interface kron_sigma
  public :: kronecker_product_pauli_matrices
  public :: kron_pauli
  public :: kron_sigma

  interface kron_pauli_recursive
     module procedure kronecker_product_pauli_recursive
  end interface kron_pauli_recursive
  interface kron_sigma_recursive
     module procedure kronecker_product_pauli_recursive
  end interface kron_sigma_recursive
  public :: kronecker_product_pauli_recursive
  public :: kron_pauli_recursive
  public :: kron_sigma_recursive



contains


  !---------------------------------------------------------------------
  !PURPOSE: return the Kronecker's product of 2 Pauli's matrices. 
  !---------------------------------------------------------------------
  function kronecker_product_pauli_matrices(sigma1,sigma2) result(gamma)
    complex(8),dimension(2,2),intent(in) :: sigma1,sigma2
    complex(8),dimension(4,4)            :: gamma
    gamma = c_kronecker_product(sigma1,2,2,sigma2,2,2)
  end function kronecker_product_pauli_matrices



  !---------------------------------------------------------------------
  !PURPOSE: return the Kronecker's product of n Pauli's matrices. 
  ! The especification of the order of the matrices is given as input on the 
  ! vector vec_ord_pm, that has dimension npm.
  !---------------------------------------------------------------------
  function kronecker_product_pauli_recursive(vec_ord_pm) result(kron_prod_n_pauli_mat)
    integer, intent(in)     :: vec_ord_pm(:)
    complex(8)              :: kron_prod_n_pauli_mat(2**size(vec_ord_pm),2**size(vec_ord_pm))
    integer                 :: d2
    complex(8)              :: M2(2,2)
    complex(8), allocatable :: M1(:,:), M1_kp_M2(:,:)
    integer                 :: npm,d1, i
    npm=size(vec_ord_pm)
    d2=2
    do i=1,npm-1
       select case(vec_ord_pm(i+1))
       case (0)
          M2 = pauli_sigma_0
       case (1)
          M2 = pauli_sigma_1
       case (2)
          M2 = pauli_sigma_2
       case (3)
          M2 = pauli_sigma_3
       end select
       d1 = 2**i
       if(i==1) then
          allocate(M1(d1,d1))
          select case(vec_ord_pm(i))
          case (0) 
             M1 = pauli_sigma_0
          case (1) 
             M1 = pauli_sigma_1
          case (2) 
             M1 = pauli_sigma_2
          case (3) 
             M1 = pauli_sigma_3
          end select
       endif
       allocate(M1_kp_M2(d1*d2,d1*d2))
       M1_kp_M2 = c_kronecker_product(M1,d1,d1,M2,d2,d2)  
       deallocate(M1)
       allocate(M1(1:d1*d2,1:d1*d2))
       M1 = M1_kp_M2
       deallocate(M1_kp_M2)
    end do
    kron_prod_n_pauli_mat = M1
    deallocate(M1)
  end function kronecker_product_pauli_recursive



  !---------------------------------------------------------------------
  !PURPOSE: Function to compute the tensor product (M1_kp_M2) of 
  ! two complex matrices M1 and M2. nr1(nr2) and nc1(nc2) are 
  ! the number of rows and columns of the Matrix M1 and M2
  !---------------------------------------------------------------------
  function c_kronecker_product(M1, nr1, nc1, M2, nr2, nc2) result(M1_kp_M2)
    integer               :: i, j
    integer,intent(in)    :: nr1,nc1,nr2,nc2
    complex(8),intent(in) :: M1(nr1,nc1), M2(nr2,nc2)
    complex(8)            :: M1_kp_M2(nr1*nr2,nc1*nc2)
    M1_kp_M2 = zero
    forall(i =1:nr1,j=1:nc1)
       M1_kp_M2(nr2*(i-1)+1 : nr2*i , nc2*(j-1)+1 : nc2*j)  =  M1(i,j)*M2
    end forall
  end function c_kronecker_product


END MODULE SF_PAULI
