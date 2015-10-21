subroutine hybrd ( fcn, n, x, fvec, xtol, maxfev, ml, mu, epsfcn, diag, mode, &
  factor, nprint, info, nfev, fjac, ldfjac, r, lr, qtf )

!*****************************************************************************80
!
!! HYBRD seeks a zero of N nonlinear equations in N variables.
!
!  Discussion:
!
!    HYBRD finds a zero of a system of N nonlinear functions in N variables
!    by a modification of the Powell hybrid method.  The user must provide a
!    subroutine which calculates the functions.  The jacobian is
!    then calculated by a forward-difference approximation.
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
!    calculates the functions.  The routine should have the form:
!
!      subroutine fcn ( n, x, fvec, iflag )
!      integer ( kind = 4 ) n
!      real fvec(n)
!      integer ( kind = 4 ) iflag
!      real x(n)
!
!    The value of IFLAG should not be changed by FCN unless
!    the user wants to terminate execution of the routine.
!    In this case set IFLAG to a negative integer.
!
!    Input, integer ( kind = 4 ) N, the number of functions and variables.
!
!    Input/output, real ( kind = 8 ) X(N).  On input, X must contain an initial
!    estimate of the solution vector.  On output X contains the final
!    estimate of the solution vector.
!
!    Output, real ( kind = 8 ) FVEC(N), the functions evaluated at the output X.
!
!    Input, real ( kind = 8 ) XTOL.  Termination occurs when the relative error
!    between two consecutive iterates is at most XTOL.  XTOL should be
!    nonnegative.
!
!    Input, integer ( kind = 4 ) MAXFEV.  Termination occurs when the number of
!    calls to FCN is at least MAXFEV by the end of an iteration.
!
!    Input, integer ( kind = 4 ) ML, MU, specify the number of subdiagonals and
!    superdiagonals within the band of the jacobian matrix.  If the jacobian
!    is not banded, set ML and MU to at least n - 1.
!
!    Input, real ( kind = 8 ) EPSFCN, is used in determining a suitable step
!    length for the forward-difference approximation.  This approximation
!    assumes that the relative errors in the functions are of the order of
!    EPSFCN.  If EPSFCN is less than the machine precision, it is assumed that
!    the relative errors in the functions are of the order of the machine
!    precision.
!
!    Input/output, real ( kind = 8 ) DIAG(N).  If MODE = 1, then DIAG is set
!    internally.  If MODE = 2, then DIAG must contain positive entries that
!    serve as multiplicative scale factors for the variables.
!
!    Input, integer ( kind = 4 ) MODE, scaling option.
!    1, variables will be scaled internally.
!    2, scaling is specified by the input DIAG vector.
!
!    Input, real ( kind = 8 ) FACTOR, determines the initial step bound.  This
!    bound is set to the product of FACTOR and the euclidean norm of DIAG*X if
!    nonzero, or else to FACTOR itself.  In most cases, FACTOR should lie
!    in the interval (0.1, 100) with 100 the recommended value.
!
!    Input, integer ( kind = 4 ) NPRINT, enables controlled printing of
!    iterates if it is positive.  In this case, FCN is called with IFLAG = 0 at
!    the beginning of the first iteration and every NPRINT iterations thereafter
!    and immediately prior to return, with X and FVEC available
!    for printing.  If NPRINT is not positive, no special calls
!    of FCN with IFLAG = 0 are made.
!
!    Output, integer ( kind = 4 ) INFO, error flag.  If the user has terminated
!    execution, INFO is set to the (negative) value of IFLAG. See the description
!    of FCN.
!    Otherwise, INFO is set as follows:
!    0, improper input parameters.
!    1, relative error between two consecutive iterates is at most XTOL.
!    2, number of calls to FCN has reached or exceeded MAXFEV.
!    3, XTOL is too small.  No further improvement in the approximate
!       solution X is possible.
!    4, iteration is not making good progress, as measured by the improvement
!       from the last five jacobian evaluations.
!    5, iteration is not making good progress, as measured by the improvement
!       from the last ten iterations.
!
!    Output, integer ( kind = 4 ) NFEV, the number of calls to FCN.
!
!    Output, real ( kind = 8 ) FJAC(LDFJAC,N), an N by N array which contains
!    the orthogonal matrix Q produced by the QR factorization of the final
!    approximate jacobian.
!
!    Input, integer ( kind = 4 ) LDFJAC, the leading dimension of the array FJAC.
!    LDFJAC must be at least N.
!
!    Output, real ( kind = 8 ) R(LR), the upper triangular matrix produced by
!    the QR factorization of the final approximate jacobian, stored rowwise.
!
!    Input, integer ( kind = 4 ) LR, the size of the R array, which must be no
!    less than (N*(N+1))/2.
!
!    Output, real ( kind = 8 ) QTF(N), contains the vector Q'*FVEC.
!
  implicit none

  integer ( kind = 4 ) ldfjac
  integer ( kind = 4 ) lr
  integer ( kind = 4 ) n

  real    ( kind = 8 ) actred
  real    ( kind = 8 ) delta
  real    ( kind = 8 ) diag(n)
  real    ( kind = 8 ) enorm
  real    ( kind = 8 ) epsfcn
  real    ( kind = 8 ) epsmch
  real    ( kind = 8 ) factor
  external             fcn
  real    ( kind = 8 ) fjac(ldfjac,n)
  real    ( kind = 8 ) fnorm
  real    ( kind = 8 ) fnorm1
  real    ( kind = 8 ) fvec(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) iflag
  integer ( kind = 4 ) info
  integer ( kind = 4 ) iter
  integer ( kind = 4 ) iwa(1)
  integer ( kind = 4 ) j
  logical              jeval
  integer ( kind = 4 ) l
  integer ( kind = 4 ) maxfev
  integer ( kind = 4 ) ml
  integer ( kind = 4 ) mode
  integer ( kind = 4 ) msum
  integer ( kind = 4 ) mu
  integer ( kind = 4 ) ncfail
  integer ( kind = 4 ) nslow1
  integer ( kind = 4 ) nslow2
  integer ( kind = 4 ) ncsuc
  integer ( kind = 4 ) nfev
  integer ( kind = 4 ) nprint
  real    ( kind = 8 ) pnorm
  real    ( kind = 8 ) prered
  real    ( kind = 8 ) qtf(n)
  real    ( kind = 8 ) r(lr)
  real    ( kind = 8 ) ratio
  logical              sing
  real    ( kind = 8 ) sum2
  real    ( kind = 8 ) temp
  real    ( kind = 8 ) wa1(n)
  real    ( kind = 8 ) wa2(n)
  real    ( kind = 8 ) wa3(n)
  real    ( kind = 8 ) wa4(n)
  real    ( kind = 8 ) x(n)
  real    ( kind = 8 ) xnorm
  real    ( kind = 8 ) xtol

  epsmch = epsilon ( epsmch )

  info = 0
  iflag = 0
  nfev = 0
!
!  Check the input parameters for errors.
!
  if ( n <= 0 ) then
    return
  else if ( xtol < 0.0D+00 ) then
    return
  else if ( maxfev <= 0 ) then
    return
  else if ( ml < 0 ) then
    return
  else if ( mu < 0 ) then
    return
  else if ( factor <= 0.0D+00 ) then
    return
  else if ( ldfjac < n ) then
    return
  else if ( lr < ( n * (n + 1) ) / 2 ) then
    return
  end if

  if ( mode == 2 ) then

    do j = 1, n
      if ( diag(j) <= 0.0D+00 ) then
        go to 300
      end if
    end do

  end if
!
!  Evaluate the function at the starting point
!  and calculate its norm.
!
  iflag = 1
  call fcn ( n, x, fvec, iflag )
  nfev = 1

  if ( iflag < 0 ) then
    go to 300
  end if

  fnorm = enorm ( n, fvec )
!
!  Determine the number of calls to FCN needed to compute the jacobian matrix.
!
  msum = min ( ml+mu+1, n )
!
!  Initialize iteration counter and monitors.
!
  iter = 1
  ncsuc = 0
  ncfail = 0
  nslow1 = 0
  nslow2 = 0
!
!  Beginning of the outer loop.
!
30 continue

     jeval = .true.
!
!  Calculate the jacobian matrix.
!
     iflag = 2
     call fdjac1 ( fcn, n, x, fvec, fjac, ldfjac, iflag, ml, mu, epsfcn )

     nfev = nfev + msum
     if ( iflag < 0 ) go to 300
!
!  Compute the QR factorization of the jacobian.
!
     call qrfac ( n, n, fjac, ldfjac, .false., iwa, 1, wa1, wa2 )
!
!  On the first iteration, if MODE is 1, scale according
!  to the norms of the columns of the initial jacobian.
!
    if ( iter == 1 ) then

      if ( mode /= 2 ) then

        diag(1:n) = wa2(1:n)
        do j = 1, n
          if ( wa2(j) == 0.0D+00 ) then
            diag(j) = 1.0D+00
          end if
        end do

      end if
!
!  On the first iteration, calculate the norm of the scaled X
!  and initialize the step bound DELTA.
!
      wa3(1:n) = diag(1:n) * x(1:n)
      xnorm = enorm ( n, wa3 )
      delta = factor * xnorm
      if ( delta == 0.0D+00 ) then
        delta = factor
      end if

    end if
!
!  Form Q' * FVEC and store in QTF.
!
     qtf(1:n) = fvec(1:n)

     do j = 1, n

       if ( fjac(j,j) /= 0.0D+00 ) then
         temp = - dot_product ( qtf(j:n), fjac(j:n,j) ) / fjac(j,j)
         qtf(j:n) = qtf(j:n) + fjac(j:n,j) * temp
       end if

     end do
!
!  Copy the triangular factor of the QR factorization into R.
!
     sing = .false.

     do j = 1, n
        l = j
        do i = 1, j-1
          r(l) = fjac(i,j)
          l = l + n - i
        end do
        r(l) = wa1(j)
        if ( wa1(j) == 0.0D+00 ) sing = .true.
     end do
!
!  Accumulate the orthogonal factor in FJAC.
!
     call qform ( n, n, fjac, ldfjac )
!
!  Rescale if necessary.
!
     if ( mode /= 2 ) then
       do j = 1, n
         diag(j) = max ( diag(j), wa2(j) )
       end do
     end if
!
!  Beginning of the inner loop.
!
180    continue
!
!  If requested, call FCN to enable printing of iterates.
!
        if ( 0 < nprint ) then
          iflag = 0
          if ( mod ( iter-1, nprint ) == 0) then
            call fcn ( n, x, fvec, iflag )
          end if
          if ( iflag < 0 ) go to 300
        end if
!
!  Determine the direction P.
!
        call dogleg ( n, r, lr, diag, qtf, delta, wa1 )
!
!  Store the direction P and X + P.
!  Calculate the norm of P.
!
        wa1(1:n) = - wa1(1:n)
        wa2(1:n) = x(1:n) + wa1(1:n)
        wa3(1:n) = diag(1:n) * wa1(1:n)

        pnorm = enorm ( n, wa3 )
!
!  On the first iteration, adjust the initial step bound.
!
        if ( iter == 1 ) then
          delta = min ( delta, pnorm )
        end if
!
!  Evaluate the function at X + P and calculate its norm.
!
        iflag = 1
        call fcn ( n, wa2, wa4, iflag )
        nfev = nfev + 1

        if ( iflag < 0 ) then
          go to 300
        end if

        fnorm1 = enorm ( n, wa4 )
!
!  Compute the scaled actual reduction.
!
        actred = -1.0D+00
        if ( fnorm1 < fnorm ) then
          actred = 1.0D+00 - ( fnorm1 / fnorm )**2
        endif
!
!  Compute the scaled predicted reduction.
!
        l = 1
        do i = 1, n
          sum2 = 0.0D+00
          do j = i, n
            sum2 = sum2 + r(l) * wa1(j)
            l = l + 1
          end do
          wa3(i) = qtf(i) + sum2
        end do

        temp = enorm ( n, wa3 )
        prered = 0.0D+00
        if ( temp < fnorm ) then
          prered = 1.0D+00 - ( temp / fnorm )**2
        end if
!
!  Compute the ratio of the actual to the predicted reduction.
!
        ratio = 0.0D+00
        if ( 0.0D+00 < prered ) then
          ratio = actred / prered
        end if
!
!  Update the step bound.
!
        if ( ratio < 0.1D+00 ) then

          ncsuc = 0
          ncfail = ncfail + 1
          delta = 0.5D+00 * delta

        else

          ncfail = 0
          ncsuc = ncsuc + 1

          if ( 0.5D+00 <= ratio .or. 1 < ncsuc ) then
            delta = max ( delta, pnorm / 0.5D+00 )
          end if

          if ( abs ( ratio - 1.0D+00 ) <= 0.1D+00 ) then
            delta = pnorm / 0.5D+00
          end if

        end if
!
!  Test for successful iteration.
!
!  Successful iteration.
!  Update X, FVEC, and their norms.
!
        if ( 0.0001D+00 <= ratio ) then
          x(1:n) = wa2(1:n)
          wa2(1:n) = diag(1:n) * x(1:n)
          fvec(1:n) = wa4(1:n)
          xnorm = enorm ( n, wa2 )
          fnorm = fnorm1
          iter = iter + 1
        end if
!
!  Determine the progress of the iteration.
!
        nslow1 = nslow1 + 1
        if ( 0.001D+00 <= actred ) then
          nslow1 = 0
        end if

        if ( jeval ) then
          nslow2 = nslow2 + 1
        end if

        if ( 0.1D+00 <= actred ) then
          nslow2 = 0
        end if
!
!  Test for convergence.
!
        if ( delta <= xtol * xnorm .or. fnorm == 0.0D+00 ) then
          info = 1
        end if

        if ( info /= 0 ) then
          go to 300
        end if
!
!  Tests for termination and stringent tolerances.
!
        if ( maxfev <= nfev ) then
          info = 2
        end if

        if ( 0.1D+00 * max ( 0.1D+00 * delta, pnorm ) <= epsmch * xnorm ) then
          info = 3
        end if

        if ( nslow2 == 5 ) then
          info = 4
        end if

        if ( nslow1 == 10 ) info = 5
        if ( info /= 0 ) go to 300
!
!  Criterion for recalculating jacobian approximation
!  by forward differences.
!
        if ( ncfail == 2 ) go to 290
!
!  Calculate the rank one modification to the jacobian
!  and update QTF if necessary.
!
        do j = 1, n
          sum2 = dot_product ( wa4(1:n), fjac(1:n,j) )
          wa2(j) = ( sum2 - wa3(j) ) / pnorm
          wa1(j) = diag(j) * ( ( diag(j) * wa1(j) ) / pnorm )
          if ( 0.0001D+00 <= ratio ) then
            qtf(j) = sum2
          end if
        end do
!
!  Compute the QR factorization of the updated jacobian.
!
        call r1updt ( n, n, r, lr, wa1, wa2, wa3, sing )
        call r1mpyq ( n, n, fjac, ldfjac, wa2, wa3 )
        call r1mpyq ( 1, n, qtf, 1, wa2, wa3 )
!
!  End of the inner loop.
!
        jeval = .false.
        go to 180
  290   continue
!
!  End of the outer loop.
!
     go to 30
  300 continue
!
!  Termination, either normal or user imposed.
!
  if ( iflag < 0 ) then
    info = iflag
  end if

  iflag = 0

  if ( 0 < nprint ) then
    call fcn ( n, x, fvec, iflag )
  end if

  return
end
