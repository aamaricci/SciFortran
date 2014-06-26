subroutine qagi ( f, bound, inf, epsabs, epsrel, result, abserr, neval, ier )
  !*****************************************************************************80
  !! QAGI estimates an integral over a semi-infinite or infinite interval.
  !
  !  Discussion:
  !
  !    The routine calculates an approximation RESULT to a definite integral   
  !      I = integral of F over (A, +Infinity), 
  !    or 
  !      I = integral of F over (-Infinity,A)
  !    or 
  !      I = integral of F over (-Infinity,+Infinity),
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
  !    Input, real(8) :: BOUND, the value of the finite endpoint of the integration
  !    range, if any, that is, if INF is 1 or -1.
  !
  !    Input, integer INF, indicates the type of integration range.
  !    1:  (  BOUND,    +Infinity),
  !    -1: ( -Infinity,  BOUND),
  !    2:  ( -Infinity, +Infinity).
  !
  !    Input, real(8) :: EPSABS, EPSREL, the absolute and relative accuracy requested.
  !
  !    Output, real(8) :: RESULT, the estimated value of the integral.
  !
  !    Output, real(8) :: ABSERR, an estimate of || I - RESULT ||.
  !
  !    Output, integer NEVAL, the number of times the integral was evaluated.
  !
  !    Output, integer IER, error indicator.
  !    0, normal and reliable termination of the routine.  It is assumed that 
  !      the requested accuracy has been achieved.
  !    > 0,  abnormal termination of the routine.  The estimates for result
  !      and error are less reliable.  It is assumed that the requested
  !      accuracy has not been achieved.
  !    1, maximum number of subdivisions allowed has been achieved.  One can 
  !      allow more subdivisions by increasing the data value of LIMIT in QAGI
  !      (and taking the according dimension adjustments into account).
  !      However, if this yields no improvement it is advised to analyze the
  !      integrand in order to determine the integration difficulties.  If the
  !      position of a local difficulty can be determined (e.g. singularity,
  !      discontinuity within the interval) one will probably gain from
  !      splitting up the interval at this point and calling the integrator 
  !      on the subranges.  If possible, an appropriate special-purpose 
  !      integrator should be used, which is designed for handling the type
  !      of difficulty involved.
  !    2, the occurrence of roundoff error is detected, which prevents the
  !      requested tolerance from being achieved.  The error may be
  !      under-estimated.
  !    3, extremely bad integrand behavior occurs at some points of the
  !      integration interval.
  !    4, the algorithm does not converge.  Roundoff error is detected in the
  !      extrapolation table.  It is assumed that the requested tolerance
  !      cannot be achieved, and that the returned result is the best which 
  !      can be obtained.
  !    5, the integral is probably divergent, or slowly convergent.  It must 
  !      be noted that divergence can occur with any other value of IER.
  !    6, the input is invalid, because INF /= 1 and INF /= -1 and INF /= 2, or
  !      epsabs < 0 and epsrel < 0.  result, abserr, neval are set to zero.
  !
  !  Local parameters:
  !
  !            the dimension of rlist2 is determined by the value of
  !            limexp in QEXTR.
  !
  !           alist     - list of left end points of all subintervals
  !                       considered up to now
  !           blist     - list of right end points of all subintervals
  !                       considered up to now
  !           rlist(i)  - approximation to the integral over
  !                       (alist(i),blist(i))
  !           rlist2    - array of dimension at least (limexp+2),
  !                       containing the part of the epsilon table
  !                       which is still needed for further computations
  !           elist(i)  - error estimate applying to rlist(i)
  !           maxerr    - pointer to the interval with largest error
  !                       estimate
  !           errmax    - elist(maxerr)
  !           erlast    - error on the interval currently subdivided
  !                       (before that subdivision has taken place)
  !           area      - sum of the integrals over the subintervals
  !           errsum    - sum of the errors over the subintervals
  !           errbnd    - requested accuracy max(epsabs,epsrel*
  !                       abs(result))
  !           *****1    - variable for the left subinterval
  !           *****2    - variable for the right subinterval
  !           last      - index for subdivision
  !           nres      - number of calls to the extrapolation routine
  !           numrl2    - number of elements currently in rlist2. if an
  !                       appropriate approximation to the compounded
  !                       integral has been obtained, it is put in
  !                       rlist2(numrl2) after numrl2 has been increased
  !                       by one.
  !           small     - length of the smallest interval considered up
  !                       to now, multiplied by 1.5
  !           erlarg    - sum of the errors over the intervals larger
  !                       than the smallest interval considered up to now
  !           extrap    - logical variable denoting that the routine
  !                       is attempting to perform extrapolation. i.e.
  !                       before subdividing the smallest interval we
  !                       try to decrease the value of erlarg.
  !           noext     - logical variable denoting that extrapolation
  !                       is no longer allowed (true-value)
  !
  implicit none
  integer, parameter :: limit = 500
  real(8)            :: abseps
  real(8)            :: abserr
  real(8)            :: alist(limit)
  real(8)            :: area
  real(8)            :: area1
  real(8)            :: area12
  real(8)            :: area2
  real(8)            :: a1
  real(8)            :: a2
  real(8)            :: blist(limit)
  real(8)            :: boun
  real(8)            :: bound
  real(8)            :: b1
  real(8)            :: b2
  real(8)            :: correc
  real(8)            :: defabs
  real(8)            :: defab1
  real(8)            :: defab2
  real(8)            :: dres
  real(8)            :: elist(limit)
  real(8)            :: epsabs
  real(8)            :: epsrel
  real(8)            :: erlarg
  real(8)            :: erlast
  real(8)            :: errbnd
  real(8)            :: errmax
  real(8)            :: error1
  real(8)            :: error2
  real(8)            :: erro12
  real(8)            :: errsum
  real(8)            :: ertest
  logical extrap
  real(8), external  :: f
  integer id
  integer ier
  integer ierro
  integer inf
  integer iord(limit)
  integer iroff1
  integer iroff2
  integer iroff3
  integer jupbnd
  integer k
  integer ksgn
  integer ktmin
  integer last
  integer maxerr
  integer neval
  logical noext
  integer nres
  integer nrmax
  integer numrl2
  real(8)            :: resabs
  real(8)            :: reseps
  real(8)            :: result
  real(8)            :: res3la(3)
  real(8)            :: rlist(limit)
  real(8)            :: rlist2(52)
  real(8)            :: small
  !
  !  Test on validity of parameters.
  !
  ier = 0
  neval = 0
  last = 0
  result = 0.0e+00
  abserr = 0.0e+00
  alist(1) = 0.0e+00
  blist(1) = 1.0e+00
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
  !  Determine the interval to be mapped onto (0,1).
  !  If INF = 2 the integral is computed as i = i1+i2, where
  !  i1 = integral of f over (-infinity,0),
  !  i2 = integral of f over (0,+infinity).
  !
  if ( inf == 2 ) then
     boun = 0.0e+00
  else
     boun = bound
  end if
  call qk15i ( f, boun, inf, 0.0d+00, 1.0d+00, result, abserr, defabs, resabs )
  !
  !  Test on accuracy.
  !
  last = 1
  rlist(1) = result
  elist(1) = abserr
  iord(1) = 1
  dres = abs ( result )
  errbnd = max ( epsabs, epsrel * dres )
  if ( abserr <= 100.0E+00 * epsilon ( defabs ) * defabs .and. &
       errbnd < abserr ) then
     ier = 2
  end if
  if ( limit == 1 ) then
     ier = 1
  end if
  if ( ier /= 0 .or. (abserr <= errbnd .and. abserr /= resabs ) .or. &
       abserr == 0.0e+00 ) go to 130
  !
  !  Initialization.
  !
  rlist2(1) = result
  errmax = abserr
  maxerr = 1
  area = result
  errsum = abserr
  abserr = huge ( abserr )
  nrmax = 1
  nres = 0
  ktmin = 0
  numrl2 = 2
  extrap = .false.
  noext = .false.
  ierro = 0
  iroff1 = 0
  iroff2 = 0
  iroff3 = 0
  if ( ( 1.0e+00 - 5.0e+01 * epsilon ( defabs ) ) * defabs <= dres ) then
     ksgn = 1
  else
     ksgn = -1
  end if
  do last = 2, limit
     !
     !  Bisect the subinterval with nrmax-th largest error estimate.
     !
     a1 = alist(maxerr)
     b1 = 5.0e-01 * ( alist(maxerr) + blist(maxerr) )
     a2 = b1
     b2 = blist(maxerr)
     erlast = errmax
     call qk15i ( f, boun, inf, a1, b1, area1, error1, resabs, defab1 )
     call qk15i ( f, boun, inf, a2, b2, area2, error2, resabs, defab2 )
     !
     !  Improve previous approximations to integral and error
     !  and test for accuracy.
     !
     area12 = area1 + area2
     erro12 = error1 + error2
     errsum = errsum + erro12 - errmax
     area = area + area12 - rlist(maxerr)
     if ( defab1 /= error1 .and. defab2 /= error2 ) then
        if ( abs ( rlist(maxerr) - area12 ) <= 1.0e-05 * abs ( area12 ) &
             .and. 9.9e-01 * errmax <= erro12 ) then
           if ( extrap ) then
              iroff2 = iroff2 + 1
           end if
           if ( .not. extrap ) then
              iroff1 = iroff1 + 1
           end if
        end if
        if ( 10 < last .and. errmax < erro12 ) then
           iroff3 = iroff3 + 1
        end if
     end if
     rlist(maxerr) = area1
     rlist(last) = area2
     errbnd = max ( epsabs, epsrel * abs ( area ) )
     !
     !  Test for roundoff error and eventually set error flag.
     !
     if ( 10 <= iroff1 + iroff2 .or. 20 <= iroff3 ) then
        ier = 2
     end if
     if ( 5 <= iroff2 ) then
        ierro = 3
     end if
     !
     !  Set error flag in the case that the number of subintervals equals LIMIT.
     !
     if ( last == limit ) then
        ier = 1
     end if
     !
     !  Set error flag in the case of bad integrand behavior
     !  at some points of the integration range.
     !
     if ( max ( abs(a1), abs(b2) ) <= (1.0e+00 + 1.0e+03 * epsilon ( a1 ) ) * &
          ( abs(a2) + 1.0e+03 * tiny ( a2 ) )) then
        ier = 4
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
     !  with NRMAX-th largest error estimate (to be bisected next).
     !
     call qsort ( limit, last, maxerr, errmax, elist, iord, nrmax )
     if ( errsum <= errbnd ) go to 115
     if ( ier /= 0 ) then
        exit
     end if
     if ( last == 2 ) then
        small = 3.75e-01
        erlarg = errsum
        ertest = errbnd
        rlist2(2) = area
        cycle
     end if
     if ( noext ) then
        cycle
     end if
     erlarg = erlarg - erlast
     if ( small < abs ( b1 - a1 ) ) then
        erlarg = erlarg + erro12
     end if
     !
     !  Test whether the interval to be bisected next is the
     !  smallest interval.
     !
     if ( .not. extrap ) then
        if ( small < abs ( blist(maxerr) - alist(maxerr) ) ) then
           cycle
        end if
        extrap = .true.
        nrmax = 2
     end if
     if ( ierro == 3 .or. erlarg <= ertest ) then
        go to 60
     end if
     !
     !  The smallest interval has the largest error.
     !  before bisecting decrease the sum of the errors over the
     !  larger intervals (erlarg) and perform extrapolation.
     !
     id = nrmax
     jupbnd = last
     if ( (2+limit/2) < last ) then
        jupbnd = limit + 3 - last
     end if
     do k = id, jupbnd
        maxerr = iord(nrmax)
        errmax = elist(maxerr)
        if ( small < abs ( blist(maxerr) - alist(maxerr) ) ) then
           go to 90
        end if
        nrmax = nrmax + 1
     end do
     !
     !  Extrapolate.
     !
60   continue
     numrl2 = numrl2 + 1
     rlist2(numrl2) = area
     call qextr ( numrl2, rlist2, reseps, abseps, res3la, nres ) 
     ktmin = ktmin+1
     if ( 5 < ktmin .and. abserr < 1.0e-03 * errsum ) then
        ier = 5
     end if
     if ( abseps < abserr ) then
        ktmin = 0
        abserr = abseps
        result = reseps
        correc = erlarg
        ertest = max ( epsabs, epsrel * abs(reseps) )
        if ( abserr <= ertest ) then
           exit
        end if
     end if
     !
     !  Prepare bisection of the smallest interval.
     !
     if ( numrl2 == 1 ) then
        noext = .true.
     end if
     if ( ier == 5 ) then
        exit
     end if
     maxerr = iord(1)
     errmax = elist(maxerr)
     nrmax = 1
     extrap = .false.
     small = small * 5.0e-01
     erlarg = errsum
90   continue
  end do
  !
  !  Set final result and error estimate.
  !
  if ( abserr == huge ( abserr ) ) then
     go to 115
  end if
  if ( ( ier + ierro ) == 0 ) then
     go to 110
  end if
  if ( ierro == 3 ) then
     abserr = abserr + correc
  end if
  if ( ier == 0 ) then
     ier = 3
  end if
  if ( result /= 0.0e+00 .and. area /= 0.0e+00) then
     go to 105
  end if
  if ( errsum < abserr ) then
     go to 115
  end if
  if ( area == 0.0e+00 ) then
     go to 130
  end if
  go to 110
105 continue
  if ( errsum / abs ( area ) < abserr / abs ( result )  ) then
     go to 115
  end if
  !
  !  Test on divergence
  !
110 continue
  if ( ksgn == (-1) .and. &
       max ( abs(result), abs(area) ) <=  defabs * 1.0e-02) go to 130
  if ( 1.0e-02 > (result/area) .or. &
       (result/area) > 1.0e+02 .or. &
       errsum > abs(area)) then
     ier = 6
  end if
  go to 130
  !
  !  Compute global integral sum.
  !
115 continue
  result = sum ( rlist(1:last) )
  abserr = errsum
130 continue
  neval = 30 * last - 15
  if ( inf == 2 ) then
     neval = 2 * neval
  end if
  if ( 2 < ier ) then
     ier = ier - 1
  end if
  return
end subroutine qagi
