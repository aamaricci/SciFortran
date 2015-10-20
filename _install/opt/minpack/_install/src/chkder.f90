subroutine chkder ( m, n, x, fvec, fjac, ldfjac, xp, fvecp, mode, err )

!*****************************************************************************80
!
!! CHKDER checks the gradients of M functions of N variables.
!
!  Discussion:
!
!    CHKDER checks the gradients of M nonlinear functions in N variables,
!    evaluated at a point X, for consistency with the functions themselves.
!
!    The user calls CHKDER twice, first with MODE = 1 and then with MODE = 2.
!
!    MODE = 1.
!      On input,
!        X contains the point of evaluation.
!      On output,
!        XP is set to a neighboring point.
!
!    Now the user must evaluate the function and gradients at X, and the
!    function at XP.  Then the subroutine is called again:
!
!    MODE = 2.
!      On input,
!        FVEC contains the function values at X,
!        FJAC contains the function gradients at X.
!        FVECP contains the functions evaluated at XP.
!      On output,
!        ERR contains measures of correctness of the respective gradients.
!
!    The subroutine does not perform reliably if cancellation or
!    rounding errors cause a severe loss of significance in the
!    evaluation of a function.  Therefore, none of the components
!    of X should be unusually small (in particular, zero) or any
!    other value which may cause loss of significance.
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
!    Input, integer ( kind = 4 ) M, is the number of functions.
!
!    Input, integer ( kind = 4 ) N, is the number of variables.
!
!    Input, real ( kind = 8 ) X(N), the point at which the jacobian is to be
!    evaluated.
!
!    Input, real ( kind = 8 ) FVEC(M), is used only when MODE = 2.
!    In that case, it should contain the function values at X.
!
!    Input, real ( kind = 8 ) FJAC(LDFJAC,N), an M by N array.  When MODE = 2,
!    FJAC(I,J) should contain the value of dF(I)/dX(J).
!
!    Input, integer ( kind = 4 ) LDFJAC, the leading dimension of the array FJAC.
!    LDFJAC must be at least M.
!
!    Output, real ( kind = 8 ) XP(N), on output with MODE = 1, is a neighboring
!    point of X, at which the function is to be evaluated.
!
!    Input, real ( kind = 8 ) FVECP(M), on input with MODE = 2, is the function
!    value at XP.
!
!    Input, integer ( kind = 4 ) MODE, should be set to 1 on the first call, and
!    2 on the second.
!
!    Output, real ( kind = 8 ) ERR(M).  On output when MODE = 2, ERR contains
!    measures of correctness of the respective gradients.  If there is no
!    severe loss of significance, then if ERR(I):
!      = 1.0D+00, the I-th gradient is correct,
!      = 0.0D+00, the I-th gradient is incorrect.
!      > 0.5D+00, the I-th gradient is probably correct.
!      < 0.5D+00, the I-th gradient is probably incorrect.
!
  implicit none

  integer ( kind = 4 ) ldfjac
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real    ( kind = 8 ) eps
  real    ( kind = 8 ) epsf
  real    ( kind = 8 ) epslog
  real    ( kind = 8 ) epsmch
  real    ( kind = 8 ) err(m)
  real    ( kind = 8 ) fjac(ldfjac,n)
  real    ( kind = 8 ) fvec(m)
  real    ( kind = 8 ) fvecp(m)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) mode
  real    ( kind = 8 ) temp
  real    ( kind = 8 ) x(n)
  real    ( kind = 8 ) xp(n)

  epsmch = epsilon ( epsmch )
  eps = sqrt ( epsmch )
!
!  MODE = 1.
!
  if ( mode == 1 ) then

     do j = 1, n
       temp = eps * abs ( x(j) )
       if ( temp == 0.0D+00 ) then
         temp = eps
       end if
       xp(j) = x(j) + temp
     end do
!
!  MODE = 2.
!
  else if ( mode == 2 ) then

     epsf = 100.0D+00 * epsmch
     epslog = log10 ( eps )

     err = 0.0D+00

     do j = 1, n
       temp = abs ( x(j) )
       if ( temp == 0.0D+00 ) then
         temp = 1.0D+00
       end if
       err(1:m) = err(1:m) + temp * fjac(1:m,j)
     end do

     do i = 1, m

       temp = 1.0D+00

       if ( fvec(i) /= 0.0D+00 .and. fvecp(i) /= 0.0D+00 .and. &
         abs ( fvecp(i)-fvec(i)) >= epsf * abs ( fvec(i) ) ) then
         temp = eps * abs ( (fvecp(i)-fvec(i)) / eps - err(i) ) &
           / ( abs ( fvec(i) ) + abs ( fvecp(i) ) )
       end if

       err(i) = 1.0D+00

       if ( epsmch < temp .and. temp < eps ) then
         err(i) = ( log10 ( temp ) - epslog ) / epslog
       end if

       if ( eps <= temp ) then
         err(i) = 0.0D+00
       end if

     end do

  end if

  return
end
