subroutine lmstr1 ( fcn, m, n, x, fvec, fjac, ldfjac, tol, info )

!*****************************************************************************80
!
!! LMSTR1 minimizes M functions in N variables using the Levenberg-Marquardt method.
!
!  Discussion:
!
!    LMSTR1 minimizes the sum of the squares of M nonlinear functions in
!    N variables by a modification of the Levenberg-Marquardt algorithm
!    which uses minimal storage.
!
!    This is done by using the more general least-squares solver
!    LMSTR.  The user must provide a subroutine which calculates
!    the functions and the rows of the jacobian.
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
!    Input, external FCN, the name of the user-supplied subroutine which
!    calculates the functions and the rows of the jacobian.
!    FCN should have the form:
!
!      subroutine fcn ( m, n, x, fvec, fjrow, iflag )
!      integer ( kind = 4 ) m,n,iflag
!      integer ( kind = 4 ) n
!      real fjrow(n)
!      real fvec(m)
!      integer ( kind = 4 ) iflag
!      real x(n)
!
!    If the input value of IFLAG is 1, calculate the functions at X and
!    return this vector in FVEC.
!    If the input value of IFLAG is I > 1, calculate the (I-1)-st row of
!    the jacobian at X, and return this vector in FJROW.
!    To terminate the algorithm, set the output value of IFLAG negative.
!
!    Input, integer ( kind = 4 ) M, the number of functions.
!
!    Input, integer ( kind = 4 ) N, the number of variables.  N must not exceed M.
!
!    Input/output, real ( kind = 8 ) X(N).  On input, X must contain an initial
!    estimate of the solution vector.  On output X contains the final
!    estimate of the solution vector.
!
!    Output, real ( kind = 8 ) FVEC(M), the functions evaluated at the output X.
!
!    Output, real ( kind = 8 ) FJAC(LDFJAC,N), an N by N array.  The upper
!    triangle contains an upper triangular matrix R such that
!
!      P' * ( JAC' * JAC ) * P = R' * R,
!
!    where P is a permutation matrix and JAC is the final calculated
!    jacobian.  Column J of P is column IPVT(J) of the identity matrix.
!    The lower triangular part of FJAC contains information generated
!    during the computation of R.
!
!    Input, integer ( kind = 4 ) LDFJAC, the leading dimension of the array FJAC.
!    LDFJAC must be at least N.
!
!    Input, real ( kind = 8 ) TOL. Termination occurs when the algorithm estimates
!    either that the relative error in the sum of squares is at most TOL
!    or that the relative error between X and the solution is at most TOL.
!    TOL should be nonnegative.
!
!    Output, integer ( kind = 4 ) INFO, error flag.  If the user has terminated
!    execution, INFO is set to the (negative) value of IFLAG. See the description
!    of FCN.  Otherwise, INFO is set as follows:
!    0, improper input parameters.
!    1, algorithm estimates that the relative error in the sum of squares
!       is at most TOL.
!    2, algorithm estimates that the relative error between X and the
!       solution is at most TOL.
!    3, conditions for INFO = 1 and INFO = 2 both hold.
!    4, FVEC is orthogonal to the columns of the jacobian to machine precision.
!    5, number of calls to FCN with IFLAG = 1 has reached 100*(N+1).
!    6, TOL is too small.  No further reduction in the sum of squares
!       is possible.
!    7, TOL is too small.  No further improvement in the approximate
!       solution X is possible.
!
  implicit none

  integer ( kind = 4 ) ldfjac
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real    ( kind = 8 ) diag(n)
  real    ( kind = 8 ) factor
  external             fcn
  real    ( kind = 8 ) fjac(ldfjac,n)
  real    ( kind = 8 ) ftol
  real    ( kind = 8 ) fvec(m)
  real    ( kind = 8 ) gtol
  integer ( kind = 4 ) info
  integer ( kind = 4 ) ipvt(n)
  integer ( kind = 4 ) maxfev
  integer ( kind = 4 ) mode
  integer ( kind = 4 ) nfev
  integer ( kind = 4 ) njev
  integer ( kind = 4 ) nprint
  real    ( kind = 8 ) qtf(n)
  real    ( kind = 8 ) tol
  real    ( kind = 8 ) x(n)
  real    ( kind = 8 ) xtol

  info = 0

  if ( n <= 0 ) then
    return
  else if ( m < n ) then
    return
  else if ( ldfjac < n ) then
    return
  else if ( tol < 0.0D+00 ) then
    return
  end if

  factor = 100.0D+00
  maxfev = 100 * ( n + 1 )
  ftol = tol
  xtol = tol
  gtol = 0.0D+00
  mode = 1
  nprint = 0

  call lmstr ( fcn, m, n, x, fvec, fjac, ldfjac, ftol, xtol, gtol, maxfev, &
    diag, mode, factor, nprint, info, nfev, njev, ipvt, qtf )

  if ( info == 8 ) then
    info = 4
  end if

  return
end
