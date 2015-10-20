subroutine lmpar ( n, r, ldr, ipvt, diag, qtb, delta, par, x, sdiag )

!*****************************************************************************80
!
!! LMPAR computes a parameter for the Levenberg-Marquardt method.
!
!  Discussion:
!
!    Given an M by N matrix A, an N by N nonsingular diagonal
!    matrix D, an M-vector B, and a positive number DELTA,
!    the problem is to determine a value for the parameter
!    PAR such that if X solves the system
!
!      A*X = B,
!      sqrt ( PAR ) * D * X = 0,
!
!    in the least squares sense, and DXNORM is the euclidean
!    norm of D*X, then either PAR is zero and
!
!      ( DXNORM - DELTA ) <= 0.1 * DELTA,
!
!    or PAR is positive and
!
!      abs ( DXNORM - DELTA) <= 0.1 * DELTA.
!
!    This subroutine completes the solution of the problem
!    if it is provided with the necessary information from the
!    QR factorization, with column pivoting, of A.  That is, if
!    A*P = Q*R, where P is a permutation matrix, Q has orthogonal
!    columns, and R is an upper triangular matrix with diagonal
!    elements of nonincreasing magnitude, then LMPAR expects
!    the full upper triangle of R, the permutation matrix P,
!    and the first N components of Q'*B.  On output
!    LMPAR also provides an upper triangular matrix S such that
!
!      P' * ( A' * A + PAR * D * D ) * P = S'* S.
!
!    S is employed within LMPAR and may be of separate interest.
!
!    Only a few iterations are generally needed for convergence
!    of the algorithm.  If, however, the limit of 10 iterations
!    is reached, then the output PAR will contain the best
!    value obtained so far.
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
!    Input, integer ( kind = 4 ) N, the order of R.
!
!    Input/output, real ( kind = 8 ) R(LDR,N),the N by N matrix.  The full
!    upper triangle must contain the full upper triangle of the matrix R.
!    On output the full upper triangle is unaltered, and the strict lower
!    triangle contains the strict upper triangle (transposed) of the upper
!    triangular matrix S.
!
!    Input, integer ( kind = 4 ) LDR, the leading dimension of R.  LDR must be
!    no less than N.
!
!    Input, integer ( kind = 4 ) IPVT(N), defines the permutation matrix P such that
!    A*P = Q*R.  Column J of P is column IPVT(J) of the identity matrix.
!
!    Input, real ( kind = 8 ) DIAG(N), the diagonal elements of the matrix D.
!
!    Input, real ( kind = 8 ) QTB(N), the first N elements of the vector Q'*B.
!
!    Input, real ( kind = 8 ) DELTA, an upper bound on the euclidean norm of D*X.
!    DELTA should be positive.
!
!    Input/output, real ( kind = 8 ) PAR.  On input an initial estimate of the
!    Levenberg-Marquardt parameter.  On output the final estimate.
!    PAR should be nonnegative.
!
!    Output, real ( kind = 8 ) X(N), the least squares solution of the system
!    A*X = B, sqrt(PAR)*D*X = 0, for the output value of PAR.
!
!    Output, real ( kind = 8 ) SDIAG(N), the diagonal elements of the upper
!    triangular matrix S.
!
  implicit none

  integer ( kind = 4 ) ldr
  integer ( kind = 4 ) n

  real    ( kind = 8 ) delta
  real    ( kind = 8 ) diag(n)
  real    ( kind = 8 ) dwarf
  real    ( kind = 8 ) dxnorm
  real    ( kind = 8 ) enorm
  real    ( kind = 8 ) gnorm
  real    ( kind = 8 ) fp
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ipvt(n)
  integer ( kind = 4 ) iter
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) l
  integer ( kind = 4 ) nsing
  real    ( kind = 8 ) par
  real    ( kind = 8 ) parc
  real    ( kind = 8 ) parl
  real    ( kind = 8 ) paru
  real    ( kind = 8 ) qnorm
  real    ( kind = 8 ) qtb(n)
  real    ( kind = 8 ) r(ldr,n)
  real    ( kind = 8 ) sdiag(n)
  real    ( kind = 8 ) sum2
  real    ( kind = 8 ) temp
  real    ( kind = 8 ) wa1(n)
  real    ( kind = 8 ) wa2(n)
  real    ( kind = 8 ) x(n)
!
!  DWARF is the smallest positive magnitude.
!
  dwarf = tiny ( dwarf )
!
!  Compute and store in X the Gauss-Newton direction.
!
!  If the jacobian is rank-deficient, obtain a least squares solution.
!
  nsing = n

  do j = 1, n
    wa1(j) = qtb(j)
    if ( r(j,j) == 0.0D+00 .and. nsing == n ) then
      nsing = j - 1
    end if
    if ( nsing < n ) then
      wa1(j) = 0.0D+00
    end if
  end do

  do k = 1, nsing
    j = nsing - k + 1
    wa1(j) = wa1(j) / r(j,j)
    temp = wa1(j)
    wa1(1:j-1) = wa1(1:j-1) - r(1:j-1,j) * temp
  end do

  do j = 1, n
    l = ipvt(j)
    x(l) = wa1(j)
  end do
!
!  Initialize the iteration counter.
!  Evaluate the function at the origin, and test
!  for acceptance of the Gauss-Newton direction.
!
  iter = 0
  wa2(1:n) = diag(1:n) * x(1:n)
  dxnorm = enorm ( n, wa2 )
  fp = dxnorm - delta

  if ( fp <= 0.1D+00 * delta ) then
    go to 220
  end if
!
!  If the jacobian is not rank deficient, the Newton
!  step provides a lower bound, PARL, for the zero of
!  the function.
!
!  Otherwise set this bound to zero.
!
  parl = 0.0D+00

  if ( nsing >= n ) then

    do j = 1, n
      l = ipvt(j)
      wa1(j) = diag(l) * ( wa2(l) / dxnorm )
    end do

    do j = 1, n
      sum2 = dot_product ( wa1(1:j-1), r(1:j-1,j) )
      wa1(j) = ( wa1(j) - sum2 ) / r(j,j)
    end do

    temp = enorm ( n, wa1 )
    parl = ( ( fp / delta ) / temp ) / temp

  end if
!
!  Calculate an upper bound, PARU, for the zero of the function.
!
  do j = 1, n
    sum2 = dot_product ( qtb(1:j), r(1:j,j) )
    l = ipvt(j)
    wa1(j) = sum2 / diag(l)
  end do

  gnorm = enorm ( n, wa1 )
  paru = gnorm / delta
  if ( paru == 0.0D+00 ) then
    paru = dwarf / min ( delta, 0.1D+00 )
  end if
!
!  If the input PAR lies outside of the interval (PARL, PARU),
!  set PAR to the closer endpoint.
!
  par = max ( par, parl )
  par = min ( par, paru )
  if ( par == 0.0D+00 ) then
    par = gnorm / dxnorm
  end if
!
!  Beginning of an iteration.
!
  150 continue

     iter = iter + 1
!
!  Evaluate the function at the current value of PAR.
!
     if ( par == 0.0D+00 ) then
       par = max ( dwarf, 0.001D+00 * paru )
     end if

     wa1(1:n) = sqrt ( par ) * diag(1:n)

     call qrsolv ( n, r, ldr, ipvt, wa1, qtb, x, sdiag )

     wa2(1:n) = diag(1:n) * x(1:n)
     dxnorm = enorm ( n, wa2 )
     temp = fp
     fp = dxnorm - delta
!
!  If the function is small enough, accept the current value of PAR.
!
    if ( abs ( fp ) <= 0.1D+00 * delta ) then
      go to 220
    end if
!
!  Test for the exceptional cases where PARL
!  is zero or the number of iterations has reached 10.
!
    if ( parl == 0.0D+00 .and. fp <= temp .and. temp < 0.0D+00 ) then
      go to 220
    else if ( iter == 10 ) then
      go to 220
    end if
!
!  Compute the Newton correction.
!
     do j = 1, n
       l = ipvt(j)
       wa1(j) = diag(l) * ( wa2(l) / dxnorm )
     end do

     do j = 1, n
       wa1(j) = wa1(j) / sdiag(j)
       temp = wa1(j)
       wa1(j+1:n) = wa1(j+1:n) - r(j+1:n,j) * temp
     end do

     temp = enorm ( n, wa1 )
     parc = ( ( fp / delta ) / temp ) / temp
!
!  Depending on the sign of the function, update PARL or PARU.
!
     if ( 0.0D+00 < fp ) then
       parl = max ( parl, par )
     else if ( fp < 0.0D+00 ) then
       paru = min ( paru, par )
     end if
!
!  Compute an improved estimate for PAR.
!
     par = max ( parl, par + parc )
!
!  End of an iteration.
!
     go to 150

220  continue
!
!  Termination.
!
  if ( iter == 0 ) then
    par = 0.0D+00
  end if

  return
end
