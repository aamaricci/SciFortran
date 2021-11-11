!-------------------------------------------------------------------------------------------
!PURPOSE: compute the determinant of a matrix using an LU factorization
!-------------------------------------------------------------------------------------------
function ddet(A) result(x)
  ! compute the determinant of a real matrix using an LU factorization
  real(8), intent(in)  :: A(:, :)
  real(8)              :: x
  integer              :: i
  integer              :: info, n
  integer, allocatable :: ipiv(:)
  real(8), allocatable :: At(:,:)
  n = size(A(1,:))
  call assert_shape(A, [n, n], "det", "A")
  allocate(At(n,n), ipiv(n))
  At = A
  call dgetrf(n, n, At, n, ipiv, info)
  if(info /= 0) then
     print *, "dgetrf returned info =", info
     if (info < 0) then
        print *, "the", -info, "-th argument had an illegal value"
     else
        print *, "U(", info, ",", info, ") is exactly zero; The factorization"
        print *, "has been completed, but the factor U is exactly"
        print *, "singular, and division by zero will occur if it is used"
        print *, "to solve a system of equations."
     end if
     stop 'det: dgetrf error'
  end if
  ! At now contains the LU of the factorization A = PLU
  ! as L has unit diagonal entries, the determinant can be computed
  ! from the product of U's diagonal entries. Additional sign changes
  ! stemming from the permutations P have to be taken into account as well.
  x = 1d0
  do i = 1,n
     if(ipiv(i) /= i) then  ! additional sign change
        x = -x*At(i,i)
     else
        x = x*At(i,i)
     endif
  end do
end function ddet
!
function zdet(A) result(x)
  ! compute the determinant of a complex matrix using an LU factorization
  complex(8), intent(in)  :: A(:, :)
  complex(8)              :: x
  integer                 :: i
  integer                 :: info, n
  integer, allocatable    :: ipiv(:)
  complex(8), allocatable :: At(:,:)
  n = size(A(1,:))
  call assert_shape(A, [n, n], "det", "A")
  allocate(At(n,n), ipiv(n))
  At = A
  call zgetrf(n, n, At, n, ipiv, info)
  if(info /= 0) then
     print *, "zgetrf returned info =", info
     if (info < 0) then
        print *, "the", -info, "-th argument had an illegal value"
     else
        print *, "U(", info, ",", info, ") is exactly zero; The factorization"
        print *, "has been completed, but the factor U is exactly"
        print *, "singular, and division by zero will occur if it is used"
        print *, "to solve a system of equations."
     end if
     stop 'zdet error: zgetrf '
  end if
  ! for details on the computation, compare the comment in ddet().
  x = one
  do i = 1,n
     if(ipiv(i) /= i) then  ! additional sign change
        x = -x*At(i,i)
     else
        x = x*At(i,i)
     endif
  end do
end function zdet





!-------------------------------------------------------------------------------------------
!PURPOSE:  construct real matrix from diagonal elements
!-------------------------------------------------------------------------------------------
pure function ddiag(x) result(A)
  real(8), intent(in)  :: x(:)
  real(8), allocatable :: A(:,:)
  integer              :: i, n
  n = size(x)
  allocate(A(n,n))
  A(:,:) = 0d0
  forall(i=1:n) A(i,i) = x(i)
end function ddiag
!-------------------------------------------------------------------------------------------
!PURPOSE:  construct complex matrix from diagonal elements
!-------------------------------------------------------------------------------------------!
pure function zdiag(x) result(A)
  complex(8), intent(in)  :: x(:)
  complex(8), allocatable :: A(:,:)
  integer                 :: i, n
  n = size(x)
  allocate(A(n,n))
  A(:,:) = zero
  forall(i=1:n) A(i,i) = x(i)
end function zdiag





!-------------------------------------------------------------------------------------------
!PURPOSE:  return the diagonal of a matrix [real]
!-------------------------------------------------------------------------------------------
pure function d_diagonal(A) result(dd)
  real(8),intent(in)           :: A(:,:)
  real(8),dimension(size(A,1)) :: dd
  integer                      :: i
  do i = 1,size(A,1)
     dd(i) = A(i,i)
  end do
end function d_diagonal
!-------------------------------------------------------------------------------------------
!PURPOSE:  return the diagonal of a matrix [complex]
!-------------------------------------------------------------------------------------------
pure function z_diagonal(A) result(dd)
  complex(8),intent(in)           :: A(:,:)
  complex(8),dimension(size(A,1)) :: dd
  integer                         :: i
  do i = 1,size(A,1)
     dd(i) = A(i,i)
  end do
end function z_diagonal




!-------------------------------------------------------------------------------------------
!PURPOSE:  return trace along the main diagonal [real]
!-------------------------------------------------------------------------------------------
pure function dtrace(A) result(t)
  real(8), intent(in) :: A(:,:)
  real(8)             :: t
  integer             :: i
  t = 0d0
  do i = 1,minval(shape(A))
     t = t + A(i,i)
  end do
end function dtrace
!-------------------------------------------------------------------------------------------
!PURPOSE:  return trace along the main diagonal [complex]
!-------------------------------------------------------------------------------------------
pure function ztrace(A) result(t)
  complex(8), intent(in) :: A(:,:)
  complex(8)             :: t
  integer                :: i
  t = zero
  do i = 1,minval(shape(A))
     t = t + A(i,i)
  end do
end function ztrace





!-------------------------------------------------------------------------------------------
!PURPOSE:  Returns the identity matrix of size n x n and type real.
!-------------------------------------------------------------------------------------------
pure function deye(n) result(A)
  integer, intent(in) :: n
  real(8)             :: A(n, n)
  integer             :: i
  A = 0d0
  do i = 1, n
     A(i,i) = 1d0
  end do
end function deye
!-------------------------------------------------------------------------------------------
!PURPOSE:  Returns the identity matrix of size n x n and type complex.
!-------------------------------------------------------------------------------------------
pure function zeye(n) result(A)
  integer, intent(in) :: n
  complex(8)          :: A(n, n)
  integer             :: i
  A = zero
  do i = 1, n
     A(i,i) = one
  end do
end function zeye






!-------------------------------------------------------------------------------------------
!PURPOSE:  Returns an array of zeros of specified size from 1 to 7 dimension
!-------------------------------------------------------------------------------------------
pure function zzeros_1(n) result(A)
  integer, intent(in) :: n
  complex(8)          :: A(n)
  A = zero
end function zzeros_1
!
pure function zzeros_2(n1,n2) result(A)
  integer, intent(in) :: n1,n2
  complex(8)          :: A(n1,n2)
  A = zero
end function zzeros_2
!
pure function zzeros_3(n1,n2,n3) result(A)
  integer, intent(in) :: n1,n2,n3
  complex(8)          :: A(n1,n2,n3)
  A = zero
end function zzeros_3
!
pure function zzeros_4(n1,n2,n3,n4) result(A)
  integer, intent(in) :: n1,n2,n3,n4
  complex(8)          :: A(n1,n2,n3,n4)
  A = zero
end function zzeros_4
!
pure function zzeros_5(n1,n2,n3,n4,n5) result(A)
  integer, intent(in) :: n1,n2,n3,n4,n5
  complex(8)          :: A(n1,n2,n3,n4,n5)
  A = zero
end function zzeros_5
!
pure function zzeros_6(n1,n2,n3,n4,n5,n6) result(A)
  integer, intent(in) :: n1,n2,n3,n4,n5,n6
  complex(8)          :: A(n1,n2,n3,n4,n5,n6)
  A = zero
end function zzeros_6
!
pure function zzeros_7(n1,n2,n3,n4,n5,n6,n7) result(A)
  integer, intent(in) :: n1,n2,n3,n4,n5,n6,n7
  complex(8)          :: A(n1,n2,n3,n4,n5,n6,n7)
  A = zero
end function zzeros_7



pure function zones_1(n) result(A)
  integer, intent(in) :: n
  complex(8)          :: A(n)
  A = one
end function zones_1
!
pure function zones_2(n1,n2) result(A)
  integer, intent(in) :: n1,n2
  complex(8)          :: A(n1,n2)
  A = one
end function zones_2
!
pure function zones_3(n1,n2,n3) result(A)
  integer, intent(in) :: n1,n2,n3
  complex(8)          :: A(n1,n2,n3)
  A = one
end function zones_3
!
pure function zones_4(n1,n2,n3,n4) result(A)
  integer, intent(in) :: n1,n2,n3,n4
  complex(8)          :: A(n1,n2,n3,n4)
  A = one
end function zones_4
!
pure function zones_5(n1,n2,n3,n4,n5) result(A)
  integer, intent(in) :: n1,n2,n3,n4,n5
  complex(8)          :: A(n1,n2,n3,n4,n5)
  A = one
end function zones_5
!
pure function zones_6(n1,n2,n3,n4,n5,n6) result(A)
  integer, intent(in) :: n1,n2,n3,n4,n5,n6
  complex(8)          :: A(n1,n2,n3,n4,n5,n6)
  A = one
end function zones_6
!
pure function zones_7(n1,n2,n3,n4,n5,n6,n7) result(A)
  integer, intent(in) :: n1,n2,n3,n4,n5,n6,n7
  complex(8)          :: A(n1,n2,n3,n4,n5,n6,n7)
  A = one
end function zones_7





