subroutine qagp ( f, a, b, npts2, points, epsabs, epsrel, result, abserr, &
     neval, ier )
  !*****************************************************************************80
  !! QAGP computes a definite integral.
  !  Discussion:
  !
  !    The routine calculates an approximation RESULT to a definite integral   
  !      I = integral of F over (A,B),
  !    hopefully satisfying
  !      || I - RESULT || <= max ( EPSABS, EPSREL * ||I|| ).
  !
  !    Interior break points of the integration interval,
  !    where local difficulties of the integrand may occur, such as
  !    singularities or discontinuities, are provided by the user.
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
  !    Input, integer NPTS2, the number of user-supplied break points within 
  !    the integration range, plus 2.  NPTS2 must be at least 2.
  !
  !    Input/output, real(8) :: POINTS(NPTS2), contains the user provided interior
  !    breakpoints in entries 1 through NPTS2-2.  If these points are not
  !    in ascending order on input, they will be sorted.
  !
  !    Input, real(8) :: EPSABS, EPSREL, the absolute and relative accuracy requested.
  !
  !    Output, real(8) :: RESULT, the estimated value of the integral.
  !
  !    Output, real(8) :: ABSERR, an estimate of || I - RESULT ||.
  !
  !    Output, integer NEVAL, the number of times the integral was evaluated.
  !
  !    Output, integer IER, return flag.
  !                     ier = 0 normal and reliable termination of the
  !                             routine. it is assumed that the requested
  !                             accuracy has been achieved.
  !                     ier > 0 abnormal termination of the routine.
  !                             the estimates for integral and error are
  !                             less reliable. it is assumed that the
  !                             requested accuracy has not been achieved.
  !                     ier = 1 maximum number of subdivisions allowed
  !                             has been achieved. one can allow more
  !                             subdivisions by increasing the data value
  !                             of limit in qagp(and taking the according
  !                             dimension adjustments into account).
  !                             however, if this yields no improvement
  !                             it is advised to analyze the integrand
  !                             in order to determine the integration
  !                             difficulties. if the position of a local
  !                             difficulty can be determined (i.e.
  !                             singularity, discontinuity within the
  !                             interval), it should be supplied to the
  !                             routine as an element of the vector
  !                             points. if necessary, an appropriate
  !                             special-purpose integrator must be used,
  !                             which is designed for handling the type
  !                             of difficulty involved.
  !                         = 2 the occurrence of roundoff error is
  !                             detected, which prevents the requested
  !                             tolerance from being achieved.
  !                             the error may be under-estimated.
  !                         = 3 extremely bad integrand behavior occurs
  !                             at some points of the integration
  !                             interval.
  !                         = 4 the algorithm does not converge. roundoff
  !                             error is detected in the extrapolation
  !                             table. it is presumed that the requested
  !                             tolerance cannot be achieved, and that
  !                             the returned result is the best which
  !                             can be obtained.
  !                         = 5 the integral is probably divergent, or
  !                             slowly convergent. it must be noted that
  !                             divergence can occur with any other value
  !                             of ier > 0.
  !                         = 6 the input is invalid because
  !                             npts2 < 2 or
  !                             break points are specified outside
  !                             the integration range or
  !                             epsabs < 0 and epsrel < 0,
  !                             or limit < npts2.
  !                             result, abserr, neval are set to zero.
  !
  !  Local parameters:
  !
  !            the dimension of rlist2 is determined by the value of
  !            limexp in QEXTR (rlist2 should be of dimension
  !            (limexp+2) at least).
  !
  !           alist     - list of left end points of all subintervals
  !                       considered up to now
  !           blist     - list of right end points of all subintervals
  !                       considered up to now
  !           rlist(i)  - approximation to the integral over
  !                       (alist(i),blist(i))
  !           rlist2    - array of dimension at least limexp+2
  !                       containing the part of the epsilon table which
  !                       is still needed for further computations
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
  !           numrl2    - number of elements in rlist2. if an appropriate
  !                       approximation to the compounded integral has
  !                       obtained, it is put in rlist2(numrl2) after
  !                       numrl2 has been increased by one.
  !           erlarg    - sum of the errors over the intervals larger
  !                       than the smallest interval considered up to now
  !           extrap    - logical variable denoting that the routine
  !                       is attempting to perform extrapolation. i.e.
  !                       before subdividing the smallest interval we
  !                       try to decrease the value of erlarg.
  !           noext     - logical variable denoting that extrapolation is
  !                       no longer allowed (true-value)
  !
  implicit none
  integer, parameter :: limit = 500
  real(8)            :: a
  real(8)            :: abseps
  real(8)            :: abserr
  real(8)            :: alist(limit)
  real(8)            :: area
  real(8)            :: area1
  real(8)            :: area12
  real(8)            :: area2
  real(8)            :: a1
  real(8)            :: a2
  real(8)            :: b
  real(8)            :: blist(limit)
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
  real(8)            :: erro12
  real(8)            :: error2
  real(8)            :: errsum
  real(8)            :: ertest
  logical extrap
  real(8), external  :: f
  integer i
  integer id
  integer ier
  integer ierro
  integer ind1
  integer ind2
  integer iord(limit)
  integer iroff1
  integer iroff2
  integer iroff3
  integer j
  integer jlow
  integer jupbnd
  integer k
  integer ksgn
  integer ktmin
  integer last
  integer levcur
  integer level(limit)
  integer levmax
  integer maxerr
  integer ndin(40)
  integer neval
  integer nint
  logical noext
  integer npts
  integer npts2
  integer nres
  integer nrmax
  integer numrl2
  real(8)            :: points(npts2)
  real(8)            :: pts(npts2)
  real(8)            :: resa
  real(8)            :: resabs
  real(8)            :: reseps
  real(8)            :: result
  real(8)            :: res3la(3)
  real(8)            :: rlist(limit)
  real(8)            :: rlist2(52)
  real(8)            :: sign
  real(8)            :: temp
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
  level(1) = 0
  npts = npts2 - 2
  if ( npts2 < 2 ) then
     ier = 6
     return
  else if ( limit <= npts .or. ( epsabs < 0.0e+00 .and. &
       epsrel < 0.0e+00) ) then
     ier = 6
     return
  end if
  !
  !  If any break points are provided, sort them into an
  !  ascending sequence.
  !
  if ( b < a ) then
     sign = -1.0e+00
  else
     sign = +1.0E+00
  end if
  pts(1) = min ( a, b )
  do i = 1, npts
     pts(i+1) = points(i)
  end do
  pts(npts+2) = max ( a, b )
  nint = npts+1
  a1 = pts(1)
  if ( npts /= 0 ) then
     do i = 1, nint
        do j = i+1, nint+1
           if ( pts(j) < pts(i) ) then
              temp   = pts(i)
              pts(i) = pts(j)
              pts(j) = temp
           end if
        end do
     end do
     if ( pts(1) /= min ( a, b ) .or. pts(nint+1) /= max ( a, b ) ) then
        ier = 6
        return
     end if
  end if
  !
  !  Compute first integral and error approximations.
  !
  resabs = 0.0e+00
  do i = 1, nint
     b1 = pts(i+1)
     call qk21 ( f, a1, b1, area1, error1, defabs, resa )
     abserr = abserr + error1
     result = result + area1
     ndin(i) = 0
     if ( error1 == resa .and. error1 /= 0.0e+00 ) then
        ndin(i) = 1
     end if
     resabs = resabs + defabs
     level(i) = 0
     elist(i) = error1
     alist(i) = a1
     blist(i) = b1
     rlist(i) = area1
     iord(i) = i
     a1 = b1
  end do
  errsum = 0.0e+00
  do i = 1, nint
     if ( ndin(i) == 1 ) then
        elist(i) = abserr
     end if
     errsum = errsum + elist(i)
  end do
  !
  !  Test on accuracy.
  !
  last = nint
  neval = 21 * nint
  dres = abs ( result )
  errbnd = max ( epsabs, epsrel * dres )
  if ( abserr <= 1.0e+02 * epsilon ( resabs ) * resabs .and. &
       abserr > errbnd ) then
     ier = 2
  end if
  if ( nint /= 1 ) then
     do i = 1, npts
        jlow = i+1
        ind1 = iord(i)
        do j = jlow, nint
           ind2 = iord(j)
           if ( elist(ind1) <= elist(ind2) ) then
              ind1 = ind2
              k = j
           end if
        end do
        if ( ind1 /= iord(i) ) then
           iord(k) = iord(i)
           iord(i) = ind1
        end if
     end do
     if ( limit < npts2 ) then
        ier = 1
     end if
  end if
  if ( ier /= 0 .or. abserr <= errbnd ) then
     return
  end if
  !
  !  Initialization
  !
  rlist2(1) = result
  maxerr = iord(1)
  errmax = elist(maxerr)
  area = result
  nrmax = 1
  nres = 0
  numrl2 = 1
  ktmin = 0
  extrap = .false.
  noext = .false.
  erlarg = errsum
  ertest = errbnd
  levmax = 1
  iroff1 = 0
  iroff2 = 0
  iroff3 = 0
  ierro = 0
  abserr = huge ( abserr )
  if ( dres >= ( 1.0e+00 - 0.5E+00 * epsilon ( resabs ) ) * resabs ) then
     ksgn = 1
  else
     ksgn = -1
  end if
  do last = npts2, limit
     !
     !  Bisect the subinterval with the NRMAX-th largest error estimate.
     !
     levcur = level(maxerr) + 1
     a1 = alist(maxerr)
     b1 = 0.5E+00 * ( alist(maxerr) + blist(maxerr) )
     a2 = b1
     b2 = blist(maxerr)
     erlast = errmax
     call qk21 ( f, a1, b1, area1, error1, resa, defab1 )
     call qk21 ( f, a2, b2, area2, error2, resa, defab2 )
     !
     !  Improve previous approximations to integral and error
     !  and test for accuracy.
     !
     neval = neval + 42
     area12 = area1 + area2
     erro12 = error1 + error2
     errsum = errsum + erro12 -errmax
     area = area + area12 - rlist(maxerr)
     if ( defab1 /= error1 .and. defab2 /= error2 ) then
        if ( abs ( rlist ( maxerr ) - area12 ) <= 1.0e-05 * abs(area12) .and. &
             erro12 >= 9.9e-01 * errmax ) then
           if ( extrap ) then
              iroff2 = iroff2+1
           else
              iroff1 = iroff1+1
           end if
        end if
        if ( last > 10 .and. erro12 > errmax ) then
           iroff3 = iroff3 + 1
        end if
     end if
     level(maxerr) = levcur
     level(last) = levcur
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
     !  Set error flag in the case that the number of subintervals
     !  equals limit.
     !
     if ( last == limit ) then
        ier = 1
     end if
     !
     !  Set error flag in the case of bad integrand behavior
     !  at a point of the integration range
     !
     if ( max ( abs(a1), abs(b2)) <= ( 1.0e+00 + 1.0e+03 * epsilon ( a1 ) )* &
          ( abs ( a2 ) + 1.0e+03 * tiny ( a2 ) ) ) then
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
     !  with nrmax-th largest error estimate (to be bisected next).
     !
     call qsort ( limit, last, maxerr, errmax, elist, iord, nrmax )
     if ( errsum <= errbnd ) then
        go to 190
     end if
     if ( ier /= 0 ) then
        exit
     end if
     if ( noext ) then
        cycle
     end if
     erlarg = erlarg - erlast
     if ( levcur+1 <= levmax ) then
        erlarg = erlarg + erro12
     end if
     !
     !  Test whether the interval to be bisected next is the
     !  smallest interval.
     !
     if ( .not. extrap ) then
        if ( level(maxerr)+1 <= levmax ) then
           cycle
        end if
        extrap = .true.
        nrmax = 2
     end if
     !
     !  The smallest interval has the largest error.
     !  Before bisecting decrease the sum of the errors over the
     !  larger intervals (erlarg) and perform extrapolation.
     !
     if ( ierro /= 3 .and. erlarg > ertest ) then
        id = nrmax
        jupbnd = last
        if ( last > (2+limit/2) ) then
           jupbnd = limit+3-last
        end if
        do k = id, jupbnd
           maxerr = iord(nrmax)
           errmax = elist(maxerr)
           if ( level(maxerr)+1 <= levmax ) then
              go to 160
           end if
           nrmax = nrmax + 1
        end do
     end if
     !
     !  Perform extrapolation.
     !
     numrl2 = numrl2 + 1
     rlist2(numrl2) = area
     if ( numrl2 <= 2 ) then
        go to 155
     end if
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
        if ( abserr < ertest ) then
           exit
        end if
     end if
     !
     !  Prepare bisection of the smallest interval.
     !
     if ( numrl2 == 1 ) then
        noext = .true.
     end if
     if ( 5 <= ier ) then
        exit
     end if
155  continue
     maxerr = iord(1)
     errmax = elist(maxerr)
     nrmax = 1
     extrap = .false.
     levmax = levmax + 1
     erlarg = errsum
160  continue
  end do
  !
  !  Set the final result.
  !
  if ( abserr == huge ( abserr ) ) then
     go to 190
  end if
  if ( ( ier + ierro ) == 0 ) then
     go to 180
  end if
  if ( ierro == 3 ) then
     abserr = abserr + correc
  end if
  if ( ier == 0 ) then
     ier = 3
  end if
  if ( result /= 0.0e+00 .and. area /= 0.0e+00 ) then
     go to 175
  end if
  if ( errsum < abserr ) then
     go to 190
  end if
  if ( area == 0.0e+00 ) then
     go to 210
  end if
  go to 180
175 continue
  if ( abserr / abs(result) > errsum / abs(area) ) then
     go to 190
  end if
  !
  !  Test on divergence.
  !
180 continue
  if ( ksgn == (-1) .and. max ( abs(result),abs(area)) <=  &
       resabs*1.0e-02 ) go to 210
  if ( 1.0e-02 > (result/area) .or. (result/area) > 1.0e+02 .or. &
       errsum > abs(area) ) then
     ier = 6
  end if
  go to 210
  !
  !  Compute global integral sum.
  !
190 continue
  result = sum ( rlist(1:last) )
  abserr = errsum
210 continue
  if ( 2 < ier ) then
     ier = ier - 1
  end if
  result = result * sign
  return
end subroutine qagp
