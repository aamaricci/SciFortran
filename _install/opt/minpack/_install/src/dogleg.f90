subroutine dogleg ( n, r, lr, diag, qtb, delta, x )

!*****************************************************************************80
!
!! DOGLEG finds the minimizing combination of Gauss-Newton and gradient steps.
!
!  Discussion:
!
!    Given an M by N matrix A, an N by N nonsingular diagonal
!    matrix D, an M-vector B, and a positive number DELTA, the
!    problem is to determine the convex combination X of the
!    Gauss-Newton and scaled gradient directions that minimizes
!    (A*X - B) in the least squares sense, subject to the
!    restriction that the euclidean norm of D*X be at most DELTA.
!
!    This subroutine completes the solution of the problem
!    if it is provided with the necessary information from the
!    QR factorization of A.  That is, if A = Q*R, where Q has
!    orthogonal columns and R is an upper triangular matrix,
!    then DOGLEG expects the full upper triangle of R and
!    the first N components of Q'*B.
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
!    Input, integer ( kind = 4 ) N, the order of the matrix R.
!
!    Input, real ( kind = 8 ) R(LR), the upper triangular matrix R stored
!    by rows.
!
!    Input, integer ( kind = 4 ) LR, the size of the R array, which must be no less
!    than (N*(N+1))/2.
!
!    Input, real ( kind = 8 ) DIAG(N), the diagonal elements of the matrix D.
!
!    Input, real ( kind = 8 ) QTB(N), the first N elements of the vector Q'* B.
!
!    Input, real ( kind = 8 ) DELTA, is a positive upper bound on the
!    euclidean norm of D*X(1:N).
!
!    Output, real ( kind = 8 ) X(N), the desired convex combination of the
!    Gauss-Newton direction and the scaled gradient direction.
!
  implicit none

  integer ( kind = 4 ) lr
  integer ( kind = 4 ) n

  real    ( kind = 8 ) alpha
  real    ( kind = 8 ) bnorm
  real    ( kind = 8 ) delta
  real    ( kind = 8 ) diag(n)
  real    ( kind = 8 ) enorm
  real    ( kind = 8 ) epsmch
  real    ( kind = 8 ) gnorm
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jj
  integer ( kind = 4 ) k
  integer ( kind = 4 ) l
  real    ( kind = 8 ) qnorm
  real    ( kind = 8 ) qtb(n)
  real    ( kind = 8 ) r(lr)
  real    ( kind = 8 ) sgnorm
  real    ( kind = 8 ) sum2
  real    ( kind = 8 ) temp
  real    ( kind = 8 ) wa1(n)
  real    ( kind = 8 ) wa2(n)
  real    ( kind = 8 ) x(n)

  epsmch = epsilon ( epsmch )
!
!  Calculate the Gauss-Newton direction.
!
  jj = ( n * ( n + 1 ) ) / 2 + 1

  do k = 1, n

     j = n - k + 1
     jj = jj - k
     l = jj + 1
     sum2 = 0.0D+00

     do i = j + 1, n
       sum2 = sum2 + r(l) * x(i)
       l = l + 1
     end do

     temp = r(jj)

     if ( temp == 0.0D+00 ) then

       l = j
       do i = 1, j
         temp = max ( temp, abs ( r(l)) )
         l = l + n - i
       end do

       if ( temp == 0.0D+00 ) then
         temp = epsmch
       else
         temp = epsmch * temp
       end if

     end if

     x(j) = ( qtb(j) - sum2 ) / temp

  end do
!
!  Test whether the Gauss-Newton direction is acceptable.
!
  wa1(1:n) = 0.0D+00
  wa2(1:n) = diag(1:n) * x(1:n)
  qnorm = enorm ( n, wa2 )

  if ( qnorm <= delta ) then
    return
  end if
!
!  The Gauss-Newton direction is not acceptable.
!  Calculate the scaled gradient direction.
!
  l = 1
  do j = 1, n
     temp = qtb(j)
     do i = j, n
       wa1(i) = wa1(i) + r(l) * temp
       l = l + 1
     end do
     wa1(j) = wa1(j) / diag(j)
  end do
!
!  Calculate the norm of the scaled gradient.
!  Test for the special case in which the scaled gradient is zero.
!
  gnorm = enorm ( n, wa1 )
  sgnorm = 0.0D+00
  alpha = delta / qnorm

  if ( gnorm /= 0.0D+00 ) then
!
!  Calculate the point along the scaled gradient which minimizes the quadratic.
!
    wa1(1:n) = ( wa1(1:n) / gnorm ) / diag(1:n)

    l = 1
    do j = 1, n
      sum2 = 0.0D+00
      do i = j, n
        sum2 = sum2 + r(l) * wa1(i)
        l = l + 1
      end do
      wa2(j) = sum2
    end do

    temp = enorm ( n, wa2 )
    sgnorm = ( gnorm / temp ) / temp
!
!  Test whether the scaled gradient direction is acceptable.
!
    alpha = 0.0D+00
!
!  The scaled gradient direction is not acceptable.
!  Calculate the point along the dogleg at which the quadratic is minimized.
!
    if ( sgnorm < delta ) then

      bnorm = enorm ( n, qtb )
      temp = ( bnorm / gnorm ) * ( bnorm / qnorm ) * ( sgnorm / delta )
      temp = temp - ( delta / qnorm ) * ( sgnorm / delta)**2 &
        + sqrt ( ( temp - ( delta / qnorm ) )**2 &
        + ( 1.0D+00 - ( delta / qnorm )**2 ) &
        * ( 1.0D+00 - ( sgnorm / delta )**2 ) )

      alpha = ( ( delta / qnorm ) * ( 1.0D+00 - ( sgnorm / delta )**2 ) ) / temp

    end if

  end if
!
!  Form appropriate convex combination of the Gauss-Newton
!  direction and the scaled gradient direction.
!
  temp = ( 1.0D+00 - alpha ) * min ( sgnorm, delta )

  x(1:n) = temp * wa1(1:n) + alpha * x(1:n)

  return
end
