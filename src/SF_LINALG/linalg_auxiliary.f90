!-------------------------------------------------------------------------------------------
!PURPOSE: compute the determinant of a real matrix using an LU factorization
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
  ! compute the determinant of a real matrix using an LU factorization
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
function ddiag(x) result(A)
  real(8), intent(in)  :: x(:)
  real(8), allocatable :: A(:,:)
  integer              :: i, n
  n = size(x)
  allocate(A(n,n))
  A(:,:) = 0d0
  forall(i=1:n) A(i,i) = x(i)
end function ddiag
!
function zdiag(x) result(A)
  complex(8), intent(in)  :: x(:)
  complex(8), allocatable :: A(:,:)
  integer                 :: i, n
  n = size(x)
  allocate(A(n,n))
  A(:,:) = zero
  forall(i=1:n) A(i,i) = x(i)
end function zdiag








!-------------------------------------------------------------------------------------------
!PURPOSE:  return trace along the main diagonal
!-------------------------------------------------------------------------------------------
function dtrace(A) result(t)
  real(8), intent(in) :: A(:,:)
  real(8)             :: t
  integer             :: i
  t = 0d0
  do i = 1,minval(shape(A))
     t = t + A(i,i)
  end do
end function dtrace

function ztrace(A) result(t)
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
function deye(n) result(A)
  integer, intent(in) :: n
  real(8)             :: A(n, n)
  integer             :: i
  A = 0d0
  do i = 1, n
     A(i,i) = 1d0
  end do
end function deye

function zeye(n) result(A)
  integer, intent(in) :: n
  complex(8)          :: A(n, n)
  integer             :: i
  A = zero
  do i = 1, n
     A(i,i) = one
  end do
end function zeye
