subroutine r1mpyq ( m, n, a, lda, v, w )

!*****************************************************************************80
!
!! R1MPYQ computes A*Q, where Q is the product of Householder transformations.
!
!  Discussion:
!
!    Given an M by N matrix A, this subroutine computes A*Q where
!    Q is the product of 2*(N - 1) transformations
!
!      GV(N-1)*...*GV(1)*GW(1)*...*GW(N-1)
!
!    and GV(I), GW(I) are Givens rotations in the (I,N) plane which
!    eliminate elements in the I-th and N-th planes, respectively.
!    Q itself is not given, rather the information to recover the
!    GV, GW rotations is supplied.
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
!    Input, integer ( kind = 4 ) M, the number of rows of A.
!
!    Input, integer ( kind = 4 ) N, the number of columns of A.
!
!    Input/output, real ( kind = 8 ) A(LDA,N), the M by N array.
!    On input, the matrix A to be postmultiplied by the orthogonal matrix Q.
!    On output, the value of A*Q.
!
!    Input, integer ( kind = 4 ) LDA, the leading dimension of A, which must not
!    be less than M.
!
!    Input, real ( kind = 8 ) V(N), W(N), contain the information necessary
!    to recover the Givens rotations GV and GW.
!
  implicit none

  integer ( kind = 4 ) lda
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real    ( kind = 8 ) a(lda,n)
  real    ( kind = 8 ) c
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real    ( kind = 8 ) s
  real    ( kind = 8 ) temp
  real    ( kind = 8 ) v(n)
  real    ( kind = 8 ) w(n)
!
!  Apply the first set of Givens rotations to A.
!
  do j = n-1, 1, -1

     if ( 1.0D+00 < abs ( v(j) ) ) then
       c = 1.0D+00 / v(j)
       s = sqrt ( 1.0D+00 - c**2 )
     else
       s = v(j)
       c = sqrt ( 1.0D+00 - s**2 )
     end if

     do i = 1, m
        temp =   c * a(i,j) - s * a(i,n)
        a(i,n) = s * a(i,j) + c * a(i,n)
        a(i,j) = temp
     end do

  end do
!
!  Apply the second set of Givens rotations to A.
!
  do j = 1, n-1

     if ( abs ( w(j) ) > 1.0D+00 ) then
       c = 1.0D+00 / w(j)
       s = sqrt ( 1.0D+00 - c**2 )
     else
       s = w(j)
       c = sqrt ( 1.0D+00 - s**2 )
     end if

     do i = 1, m
        temp =     c * a(i,j) + s * a(i,n)
        a(i,n) = - s * a(i,j) + c * a(i,n)
        a(i,j) = temp
     end do

  end do

  return
end
