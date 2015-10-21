subroutine hybrj1 ( fcn, n, x, fvec, fjac, ldfjac, tol, info )

!*****************************************************************************80
!
!! HYBRJ1 seeks a zero of N nonlinear equations in N variables by Powell's method.
!
!  Discussion:
!
!    HYBRJ1 finds a zero of a system of N nonlinear functions in N variables
!    by a modification of the Powell hybrid method.  This is done by using the
!    more general nonlinear equation solver HYBRJ.  The user
!    must provide a subroutine which calculates the functions
!    and the jacobian.
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
!    calculates the functions and the jacobian.  FCN should have the form:
!
!      subroutine fcn ( n, x, fvec, fjac, ldfjac, iflag )
!      integer ( kind = 4 ) ldfjac
!      integer ( kind = 4 ) n
!      real fjac(ldfjac,n)
!      real fvec(n)
!      integer ( kind = 4 ) iflag
!      real x(n)
!
!    If IFLAG = 1 on intput, FCN should calculate the functions at X and
!    return this vector in FVEC.
!    If IFLAG = 2 on input, FCN should calculate the jacobian at X and
!    return this matrix in FJAC.
!    To terminate the algorithm, the user may set IFLAG negative.
!
!    Input, integer ( kind = 4 ) N, the number of functions and variables.
!
!    Input/output, real ( kind = 8 ) X(N).  On input, X must contain an initial
!    estimate of the solution vector.  On output X contains the final
!    estimate of the solution vector.
!
!    Output, real ( kind = 8 ) FVEC(N), the functions evaluated at the output X.
!
!    Output, real ( kind = 8 ) FJAC(LDFJAC,N), an N by N array which contains
!    the orthogonal matrix Q produced by the QR factorization of the final
!    approximate jacobian.
!
!    Input, integer ( kind = 4 ) LDFJAC, the leading dimension of the array FJAC.
!    LDFJAC must be at least N.
!
!    Input, real ( kind = 8 ) TOL.  Termination occurs when the algorithm
!    estimates that the relative error between X and the solution is at most
!    TOL.  TOL should be nonnegative.
!
!    Output, integer ( kind = 4 ) INFO, error flag.  If the user has terminated
!    execution, INFO is set to the (negative) value of IFLAG. See the description
!    of FCN.  Otherwise, INFO is set as follows:
!    0, improper input parameters.
!    1, algorithm estimates that the relative error between X and the
!       solution is at most TOL.
!    2, number of calls to FCN with IFLAG = 1 has reached 100*(N+1).
!    3, TOL is too small.  No further improvement in the approximate
!       solution X is possible.
!    4, iteration is not making good progress.
!
  implicit none

  integer ( kind = 4 ) ldfjac
  integer ( kind = 4 ) n

  real    ( kind = 8 ) diag(n)
  real    ( kind = 8 ) factor
  external             fcn
  real    ( kind = 8 ) fjac(ldfjac,n)
  real    ( kind = 8 ) fvec(n)
  integer ( kind = 4 ) info
  integer ( kind = 4 ) j
  integer ( kind = 4 ) lr
  integer ( kind = 4 ) maxfev
  integer ( kind = 4 ) mode
  integer ( kind = 4 ) nfev
  integer ( kind = 4 ) njev
  integer ( kind = 4 ) nprint
  real    ( kind = 8 ) qtf(n)
  real    ( kind = 8 ) r((n*(n+1))/2)
  real    ( kind = 8 ) tol
  real    ( kind = 8 ) x(n)
  real    ( kind = 8 ) xtol

  info = 0

  if ( n <= 0 ) then
    return
  else if ( ldfjac < n ) then
    return
  else if ( tol < 0.0D+00 ) then
    return
  end if

  maxfev = 100 * ( n + 1 )
  xtol = tol
  mode = 2
  diag(1:n) = 1.0D+00
  factor = 100.0D+00
  nprint = 0
  lr = ( n * ( n + 1 ) ) / 2

  call hybrj ( fcn, n, x, fvec, fjac, ldfjac, xtol, maxfev, diag, mode, &
    factor, nprint, info, nfev, njev, r, lr, qtf )

  if ( info == 5 ) then
    info = 4
  end if

  return
end
