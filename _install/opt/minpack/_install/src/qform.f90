subroutine qform ( m, n, q, ldq )

!*****************************************************************************80
!
!! QFORM produces the explicit QR factorization of a matrix.
!
!  Discussion:
!
!    The QR factorization of a matrix is usually accumulated in implicit
!    form, that is, as a series of orthogonal transformations of the
!    original matrix.  This routine carries out those transformations,
!    to explicitly exhibit the factorization construced by QRFAC.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 April 2010
!
!  Author:
!
!    Original FORTRAN77 version by Jorge More, Burton Garbow, Kenneth Hillstrom.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Jorge More, Burton Garbow, Kenneth Hillstrom,
!    User Guide for MINPACK-1,
!    Technical Report ANL-80-74,
!    Argonne National Laboratory, 1980.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, is a positive integer input variable set
!    to the number of rows of A and the order of Q.
!
!    Input, integer ( kind = 4 ) N, is a positive integer input variable set
!    to the number of columns of A.
!
!    Input/output, real ( kind = 8 ) Q(LDQ,M).  Q is an M by M array.
!    On input the full lower trapezoid in the first min(M,N) columns of Q
!    contains the factored form.
!    On output, Q has been accumulated into a square matrix.
!
!    Input, integer ( kind = 4 ) LDQ, is a positive integer input variable not less
!    than M which specifies the leading dimension of the array Q.
!
  implicit none

  integer ( kind = 4 ) ldq
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) l
  integer ( kind = 4 ) minmn
  real    ( kind = 8 ) q(ldq,m)
  real    ( kind = 8 ) temp
  real    ( kind = 8 ) wa(m)

  minmn = min ( m, n )

  do j = 2, minmn
    q(1:j-1,j) = 0.0D+00
  end do
!
!  Initialize remaining columns to those of the identity matrix.
!
  q(1:m,n+1:m) = 0.0D+00

  do j = n+1, m
    q(j,j) = 1.0D+00
  end do
!
!  Accumulate Q from its factored form.
!
  do l = 1, minmn

    k = minmn - l + 1

    wa(k:m) = q(k:m,k)

    q(k:m,k) = 0.0D+00
    q(k,k) = 1.0D+00

    if ( wa(k) /= 0.0D+00 ) then

      do j = k, m
        temp = dot_product ( wa(k:m), q(k:m,j) ) / wa(k)
        q(k:m,j) = q(k:m,j) - temp * wa(k:m)
      end do

    end if

  end do

  return
end
