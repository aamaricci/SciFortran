!===============================================================================
!***** Quadrature Adaptive Generic [Extended] *****!
!===============================================================================
subroutine qag ( f, a, b, epsabs, epsrel, key, result, abserr, neval, ier )
  !*****************************************************************************80
  !! QAG approximates an integral over a finite interval.
  !  Discussion:
  !    The routine calculates an approximation RESULT to a definite integral   
  !      I = integral of F over (A,B),
  !    hopefully satisfying
  !      || I - RESULT || <= max ( EPSABS, EPSREL * ||I|| ).
  !
  !    QAG is a simple globally adaptive integrator using the strategy of 
  !    Aind (Piessens, 1973).  It is possible to choose between 6 pairs of
  !    Gauss-Kronrod quadrature formulae for the rule evaluation component. 
  !    The pairs of high degree of precision are suitable for handling
  !    integration difficulties due to a strongly oscillating integrand.
  !
  !  Parameters:
  !    Input, external real(8) :: F, the name of the function routine, of the form
  !      function f ( x )
  !      real(8) :: f
  !      real(8) :: x
  !    which evaluates the integrand function.
  !
  !    Input, real(8) :: A, B, the limits of integration.
  !
  !    Input, real(8) :: EPSABS, EPSREL, the absolute and relative accuracy requested.
  !
  !    Input, integer KEY, chooses the order of the local integration rule:
  !    1,  7 Gauss points, 15 Gauss-Kronrod points,
  !    2, 10 Gauss points, 21 Gauss-Kronrod points,
  !    3, 15 Gauss points, 31 Gauss-Kronrod points,
  !    4, 20 Gauss points, 41 Gauss-Kronrod points,
  !    5, 25 Gauss points, 51 Gauss-Kronrod points,
  !    6, 30 Gauss points, 61 Gauss-Kronrod points.
  !
  !    Output, real(8) :: RESULT, the estimated value of the integral.
  !
  !    Output, real(8) :: ABSERR, an estimate of || I - RESULT ||.
  !
  !    Output, integer NEVAL, the number of times the integral was evaluated.
  !
  !    Output, integer IER, return code.
  !    0, normal and reliable termination of the routine.  It is assumed that the 
  !      requested accuracy has been achieved.
  !    1, maximum number of subdivisions allowed has been achieved.  One can 
  !      allow more subdivisions by increasing the value of LIMIT in QAG. 
  !      However, if this yields no improvement it is advised to analyze the
  !      integrand to determine the integration difficulties.  If the position
  !      of a local difficulty can be determined, such as a singularity or
  !      discontinuity within the interval) one will probably gain from 
  !      splitting up the interval at this point and calling the integrator 
  !      on the subranges.  If possible, an appropriate special-purpose 
  !      integrator should be used which is designed for handling the type 
  !      of difficulty involved.
  !    2, the occurrence of roundoff error is detected, which prevents the
  !      requested tolerance from being achieved.
  !    3, extremely bad integrand behavior occurs at some points of the
  !      integration interval.
  !    6, the input is invalid, because EPSABS < 0 and EPSREL < 0.
  !
  !  Local parameters:
  !
  !    LIMIT is the maximum number of subintervals allowed in
  !    the subdivision process of QAGE.
  implicit none
  integer, parameter :: limit = 500
  real(8)            :: a
  real(8)            :: abserr
  real(8)            :: alist(limit)
  real(8)            :: b
  real(8)            :: blist(limit)
  real(8)            :: elist(limit)
  real(8)            :: epsabs
  real(8)            :: epsrel
  real(8), external  :: f
  integer ier
  integer iord(limit)
  integer key
  integer last
  integer neval
  real(8)            :: result
  real(8)            :: rlist(limit)
  call qage ( f, a, b, epsabs, epsrel, key, limit, result, abserr, neval, &
       ier, alist, blist, rlist, elist, iord, last )
  return
end subroutine qag

!
!
!

subroutine qage ( f, a, b, epsabs, epsrel, key, limit, result, abserr, neval, &
     ier, alist, blist, rlist, elist, iord, last )
  !*****************************************************************************80
  !! QAGE estimates a definite integral.
  !  Discussion:
  !
  !    The routine calculates an approximation RESULT to a definite integral   
  !      I = integral of F over (A,B),
  !    hopefully satisfying
  !      || I - RESULT || <= max ( EPSABS, EPSREL * ||I|| ).
  !
  !  Parameters:
  !
  !    Input, external real(8) :: F, the name of the function routine, of the form
  !      function f ( x )
  !      real(8) :: f
  !      real(8) :: x
  !    which evaluates the integrand function.
  !
  !    Input, real(8) :: A, B, the limits of integration.
  !
  !    Input, real(8) :: EPSABS, EPSREL, the absolute and relative accuracy requested.
  !
  !    Input, integer KEY, chooses the order of the local integration rule:
  !    1,  7 Gauss points, 15 Gauss-Kronrod points,
  !    2, 10 Gauss points, 21 Gauss-Kronrod points,
  !    3, 15 Gauss points, 31 Gauss-Kronrod points,
  !    4, 20 Gauss points, 41 Gauss-Kronrod points,
  !    5, 25 Gauss points, 51 Gauss-Kronrod points,
  !    6, 30 Gauss points, 61 Gauss-Kronrod points.
  !
  !    Input, integer LIMIT, the maximum number of subintervals that
  !    can be used.
  !
  !    Output, real(8) :: RESULT, the estimated value of the integral.
  !
  !    Output, real(8) :: ABSERR, an estimate of || I - RESULT ||.
  !
  !    Output, integer NEVAL, the number of times the integral was evaluated.
  !
  !    Output, integer IER, return code.
  !    0, normal and reliable termination of the routine.  It is assumed that the 
  !      requested accuracy has been achieved.
  !    1, maximum number of subdivisions allowed has been achieved.  One can 
  !      allow more subdivisions by increasing the value of LIMIT in QAG. 
  !      However, if this yields no improvement it is advised to analyze the
  !      integrand to determine the integration difficulties.  If the position
  !      of a local difficulty can be determined, such as a singularity or
  !      discontinuity within the interval) one will probably gain from 
  !      splitting up the interval at this point and calling the integrator 
  !      on the subranges.  If possible, an appropriate special-purpose 
  !      integrator should be used which is designed for handling the type 
  !      of difficulty involved.
  !    2, the occurrence of roundoff error is detected, which prevents the
  !      requested tolerance from being achieved.
  !    3, extremely bad integrand behavior occurs at some points of the
  !      integration interval.
  !    6, the input is invalid, because EPSABS < 0 and EPSREL < 0.
  !
  !    Workspace, real(8) :: ALIST(LIMIT), BLIST(LIMIT), contains in entries 1 
  !    through LAST the left and right ends of the partition subintervals.
  !
  !    Workspace, real(8) :: RLIST(LIMIT), contains in entries 1 through LAST
  !    the integral approximations on the subintervals.
  !
  !    Workspace, real(8) :: ELIST(LIMIT), contains in entries 1 through LAST
  !    the absolute error estimates on the subintervals.
  !
  !    Output, integer IORD(LIMIT), the first K elements of which are pointers 
  !    to the error estimates over the subintervals, such that
  !    elist(iord(1)), ..., elist(iord(k)) form a decreasing sequence, with
  !    k = last if last <= (limit/2+2), and k = limit+1-last otherwise.
  !
  !    Output, integer LAST, the number of subintervals actually produced 
  !    in the subdivision process.
  !
  !  Local parameters:
  !    alist     - list of left end points of all subintervals
  !                       considered up to now
  !    blist     - list of right end points of all subintervals
  !                       considered up to now
  !    elist(i)  - error estimate applying to rlist(i)
  !    maxerr    - pointer to the interval with largest error estimate
  !    errmax    - elist(maxerr)
  !    area      - sum of the integrals over the subintervals
  !    errsum    - sum of the errors over the subintervals
  !    errbnd    - requested accuracy max(epsabs,epsrel*abs(result))
  !    *****1    - variable for the left subinterval
  !    *****2    - variable for the right subinterval
  !    last      - index for subdivision
  implicit none
  integer limit
  real(8) :: a
  real(8) :: abserr
  real(8) :: alist(limit)
  real(8) :: area
  real(8) :: area1
  real(8) :: area12
  real(8) :: area2
  real(8) :: a1
  real(8) :: a2
  real(8) :: b
  real(8) :: blist(limit)
  real(8) :: b1
  real(8) :: b2
  real(8) :: c
  real(8) :: defabs
  real(8) :: defab1
  real(8) :: defab2
  real(8) :: elist(limit)
  real(8) :: epsabs
  real(8) :: epsrel
  real(8) :: errbnd
  real(8) :: errmax
  real(8) :: error1
  real(8) :: error2
  real(8) :: erro12
  real(8) :: errsum
  real(8), external :: f
  integer ier
  integer iord(limit)
  integer iroff1
  integer iroff2
  integer key
  integer keyf
  integer last
  integer maxerr
  integer neval
  integer nrmax
  real(8) :: resabs
  real(8) :: result
  real(8) :: rlist(limit)
  !
  !  Test on validity of parameters.
  !
  ier = 0
  neval = 0
  last = 0
  result = 0.0e+00
  abserr = 0.0e+00
  alist(1) = a
  blist(1) = b
  rlist(1) = 0.0e+00
  elist(1) = 0.0e+00
  iord(1) = 0
  if ( epsabs < 0.0e+00 .and. epsrel < 0.0e+00 ) then
     ier = 6
     return
  end if
  !
  !  First approximation to the integral.
  !
  keyf = key
  keyf = max ( keyf, 1 )
  keyf = min ( keyf, 6 )
  c = keyf
  neval = 0
  if ( keyf == 1 ) then
     call qk15 ( f, a, b, result, abserr, defabs, resabs )
  else if ( keyf == 2 ) then
     call qk21 ( f, a, b, result, abserr, defabs, resabs )
  else if ( keyf == 3 ) then
     call qk31 ( f, a, b, result, abserr, defabs, resabs )
  else if ( keyf == 4 ) then
     call qk41 ( f, a, b, result, abserr, defabs, resabs )
  else if ( keyf == 5 ) then
     call qk51 ( f, a, b, result, abserr, defabs, resabs )
  else if ( keyf == 6 ) then
     call qk61 ( f, a, b, result, abserr, defabs, resabs )
  end if
  last = 1
  rlist(1) = result
  elist(1) = abserr
  iord(1) = 1
  !
  !  Test on accuracy.
  !
  errbnd = max ( epsabs, epsrel * abs ( result ) )
  if ( abserr <= 5.0e+01 * epsilon ( defabs ) * defabs .and. &
       errbnd < abserr ) then
     ier = 2
  end if
  if ( limit == 1 ) then
     ier = 1
  end if
  if ( ier /= 0 .or. &
       ( abserr <= errbnd .and. abserr /= resabs ) .or. &
       abserr == 0.0e+00 ) then
     if ( keyf /= 1 ) then
        neval = (10*keyf+1) * (2*neval+1)
     else
        neval = 30 * neval + 15
     end if
     return
  end if
  !
  !  Initialization.
  !
  errmax = abserr
  maxerr = 1
  area = result
  errsum = abserr
  nrmax = 1
  iroff1 = 0
  iroff2 = 0
  do last = 2, limit
     !
     !  Bisect the subinterval with the largest error estimate.
     !
     a1 = alist(maxerr)
     b1 = 0.5E+00 * ( alist(maxerr) + blist(maxerr) )
     a2 = b1
     b2 = blist(maxerr)
     if ( keyf == 1 ) then
        call qk15 ( f, a1, b1, area1, error1, resabs, defab1 )
     else if ( keyf == 2 ) then
        call qk21 ( f, a1, b1, area1, error1, resabs, defab1 )
     else if ( keyf == 3 ) then
        call qk31 ( f, a1, b1, area1, error1, resabs, defab1 )
     else if ( keyf == 4 ) then
        call qk41 ( f, a1, b1, area1, error1, resabs, defab1)
     else if ( keyf == 5 ) then
        call qk51 ( f, a1, b1, area1, error1, resabs, defab1 )
     else if ( keyf == 6 ) then
        call qk61 ( f, a1, b1, area1, error1, resabs, defab1 )
     end if
     if ( keyf == 1 ) then
        call qk15 ( f, a2, b2, area2, error2, resabs, defab2 )
     else if ( keyf == 2 ) then
        call qk21 ( f, a2, b2, area2, error2, resabs, defab2 )
     else if ( keyf == 3 ) then
        call qk31 ( f, a2, b2, area2, error2, resabs, defab2 )
     else if ( keyf == 4 ) then
        call qk41 ( f, a2, b2, area2, error2, resabs, defab2 )
     else if ( keyf == 5 ) then
        call qk51 ( f, a2, b2, area2, error2, resabs, defab2 )
     else if ( keyf == 6 ) then
        call qk61 ( f, a2, b2, area2, error2, resabs, defab2 )
     end if
     !
     !  Improve previous approximations to integral and error and
     !  test for accuracy.
     !
     neval = neval + 1
     area12 = area1 + area2
     erro12 = error1 + error2
     errsum = errsum + erro12 - errmax
     area = area + area12 - rlist(maxerr)
     if ( defab1 /= error1 .and. defab2 /= error2 ) then
        if ( abs ( rlist(maxerr) - area12 ) <= 1.0e-05 * abs ( area12 ) &
             .and. 9.9e-01 * errmax <= erro12 ) then
           iroff1 = iroff1 + 1
        end if
        if ( 10 < last .and. errmax < erro12 ) then
           iroff2 = iroff2 + 1
        end if
     end if
     rlist(maxerr) = area1
     rlist(last) = area2
     errbnd = max ( epsabs, epsrel * abs ( area ) )
     !
     !  Test for roundoff error and eventually set error flag.
     !
     if ( errbnd < errsum ) then
        if ( 6 <= iroff1 .or. 20 <= iroff2 ) then
           ier = 2
        end if
        !
        !  Set error flag in the case that the number of subintervals
        !  equals limit.
        !
        if ( last == limit ) then
           ier = 1
        end if
        !
        !  Set error flag in the case of bad integrand behavior
        !  at a point of the integration range.
        !
        if ( max ( abs ( a1 ), abs ( b2 ) ) <= ( 1.0e+00 + c * 1.0e+03 * &
             epsilon ( a1 ) ) * ( abs ( a2 ) + 1.0e+04 * tiny ( a2 ) ) ) then
           ier = 3
        end if
     end if
     !
     !  Append the newly-created intervals to the list.
     !
     if ( error2 <= error1 ) then
        alist(last) = a2
        blist(maxerr) = b1
        blist(last) = b2
        elist(maxerr) = error1
        elist(last) = error2
     else
        alist(maxerr) = a2
        alist(last) = a1
        blist(last) = b1
        rlist(maxerr) = area2
        rlist(last) = area1
        elist(maxerr) = error2
        elist(last) = error1
     end if
     !
     !  Call QSORT to maintain the descending ordering
     !  in the list of error estimates and select the subinterval
     !  with the largest error estimate (to be bisected next).
     !
     call qsort ( limit, last, maxerr, errmax, elist, iord, nrmax )
     if ( ier /= 0 .or. errsum <= errbnd ) then
        exit
     end if
  end do
  !
  !  Compute final result.
  !
  result = sum ( rlist(1:last) )
  abserr = errsum
  if ( keyf /= 1 ) then
     neval = ( 10 * keyf + 1 ) * ( 2 * neval + 1 )
  else
     neval = 30 * neval + 15
  end if
  return
end subroutine qage
