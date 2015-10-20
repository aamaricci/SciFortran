subroutine qags ( f, a, b, epsabs, epsrel, result, abserr, neval, ier )
  !*****************************************************************************80
  !! QAGS estimates the integral of a function.
  !  Discussion:
  !
  !    The routine calculates an approximation RESULT to a definite integral   
  !      I = integral of F over (A,B),
  !    hopefully satisfying
  !      || I - RESULT || <= max ( EPSABS, EPSREL * ||I|| ).
  !
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
  !    Output, real(8) :: RESULT, the estimated value of the integral.
  !
  !    Output, real(8) :: ABSERR, an estimate of || I - RESULT ||.
  !
  !    Output, integer NEVAL, the number of times the integral was evaluated.
  !
  !    Output, integer IER, error flag.
  !                     ier = 0 normal and reliable termination of the
  !                             routine. it is assumed that the requested
  !                             accuracy has been achieved.
  !                     ier > 0 abnormal termination of the routine
  !                             the estimates for integral and error are
  !                             less reliable. it is assumed that the
  !                             requested accuracy has not been achieved.
  !                         = 1 maximum number of subdivisions allowed
  !                             has been achieved. one can allow more sub-
  !                             divisions by increasing the data value of
  !                             limit in qags (and taking the according
  !                             dimension adjustments into account).
  !                             however, if this yields no improvement
  !                             it is advised to analyze the integrand
  !                             in order to determine the integration
  !                             difficulties. if the position of a
  !                             local difficulty can be determined (e.g.
  !                             singularity, discontinuity within the
  !                             interval) one will probably gain from
  !                             splitting up the interval at this point
  !                             and calling the integrator on the sub-
  !                             ranges. if possible, an appropriate
  !                             special-purpose integrator should be used,
  !                             which is designed for handling the type
  !                             of difficulty involved.
  !                         = 2 the occurrence of roundoff error is detec-
  !                             ted, which prevents the requested
  !                             tolerance from being achieved.
  !                             the error may be under-estimated.
  !                         = 3 extremely bad integrand behavior occurs
  !                             at some  points of the integration
  !                             interval.
  !                         = 4 the algorithm does not converge. roundoff
  !                             error is detected in the extrapolation
  !                             table. it is presumed that the requested
  !                             tolerance cannot be achieved, and that the
  !                             returned result is the best which can be
  !                             obtained.
  !                         = 5 the integral is probably divergent, or
  !                             slowly convergent. it must be noted that
  !                             divergence can occur with any other value
  !                             of ier.
  !                         = 6 the input is invalid, because
  !                             epsabs < 0 and epsrel < 0,
  !                             result, abserr and neval are set to zero.
  !
  !  Local Parameters:
  !
  !           alist     - list of left end points of all subintervals
  !                       considered up to now
  !           blist     - list of right end points of all subintervals
  !                       considered up to now
  !           rlist(i)  - approximation to the integral over
  !                       (alist(i),blist(i))
  !           rlist2    - array of dimension at least limexp+2 containing
  !                       the part of the epsilon table which is still
  !                       needed for further computations
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
  !           *****1    - variable for the left interval
  !           *****2    - variable for the right interval
  !           last      - index for subdivision
  !           nres      - number of calls to the extrapolation routine
  !           numrl2    - number of elements currently in rlist2. if an
  !                       appropriate approximation to the compounded
  !                       integral has been obtained it is put in
  !                       rlist2(numrl2) after numrl2 has been increased
  !                       by one.
  !           small     - length of the smallest interval considered
  !                       up to now, multiplied by 1.5
  !           erlarg    - sum of the errors over the intervals larger
  !                       than the smallest interval considered up to now
  !           extrap    - logical variable denoting that the routine is
  !                       attempting to perform extrapolation i.e. before
  !                       subdividing the smallest interval we try to
  !                       decrease the value of erlarg.
  !           noext     - logical variable denoting that extrapolation
  !                       is no longer allowed (true value)
  !
  implicit none
  integer, parameter :: limit = 500
  real(8) :: a
  real(8) :: abseps
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
  real(8) :: correc
  real(8) :: defabs
  real(8) :: defab1
  real(8) :: defab2
  real(8) :: dres
  real(8) :: elist(limit)
  real(8) :: epsabs
  real(8) :: epsrel
  real(8) :: erlarg
  real(8) :: erlast
  real(8) :: errbnd
  real(8) :: errmax
  real(8) :: error1
  real(8) :: error2
  real(8) :: erro12
  real(8) :: errsum
  real(8) :: ertest
  logical extrap
  real(8), external :: f
  integer id
  integer ier
  integer ierro
  integer iord(limit)
  integer iroff1
  integer iroff2
  integer iroff3
  integer jupbnd
  integer k
  integer ksgn
  integer ktmin
  integer last
  logical noext
  integer maxerr
  integer neval
  integer nres
  integer nrmax
  integer numrl2
  real(8) :: resabs
  real(8) :: reseps
  real(8) :: result
  real(8) :: res3la(3)
  real(8) :: rlist(limit)
  real(8) :: rlist2(52)
  real(8) :: small
  !
  !  The dimension of rlist2 is determined by the value of
  !  limexp in QEXTR (rlist2 should be of dimension
  !  (limexp+2) at least).
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
  if ( epsabs < 0.0e+00 .and. epsrel < 0.0e+00 ) then
     ier = 6
     return
  end if
  !
  !  First approximation to the integral.
  !
  ierro = 0
  call qk21 ( f, a, b, result, abserr, defabs, resabs )
  !
  !  Test on accuracy.
  !
  dres = abs ( result )
  errbnd = max ( epsabs, epsrel * dres )
  last = 1
  rlist(1) = result
  elist(1) = abserr
  iord(1) = 1
  if ( abserr <= 1.0e+02 * epsilon ( defabs ) * defabs .and. &
       abserr > errbnd ) then
     ier = 2
  end if
  if ( limit == 1 ) then
     ier = 1
  end if
  if ( ier /= 0 .or. (abserr <= errbnd .and. abserr /= resabs ) .or. &
       abserr == 0.0e+00 ) go to 140
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
  numrl2 = 2
  ktmin = 0
  extrap = .false.
  noext = .false.
  iroff1 = 0
  iroff2 = 0
  iroff3 = 0
  if ( dres >= (1.0e+00-5.0e+01* epsilon ( defabs ) ) * defabs ) then
     ksgn = 1
  else
     ksgn = -1
  end if
  do last = 2, limit
     !
     !  Bisect the subinterval with the nrmax-th largest error estimate.
     !
     a1 = alist(maxerr)
     b1 = 5.0e-01 * ( alist(maxerr) + blist(maxerr) )
     a2 = b1
     b2 = blist(maxerr)
     erlast = errmax
     call qk21 ( f, a1, b1, area1, error1, resabs, defab1 )
     call qk21 ( f, a2, b2, area2, error2, resabs, defab2 )
     !
     !  Improve previous approximations to integral and error
     !  and test for accuracy.
     !
     area12 = area1+area2
     erro12 = error1+error2
     errsum = errsum+erro12-errmax
     area = area+area12-rlist(maxerr)
     if ( defab1 == error1 .or. defab2 == error2 ) go to 15
     if ( abs ( rlist(maxerr) - area12) > 1.0e-05 * abs(area12) &
          .or. erro12 < 9.9e-01 * errmax ) go to 10
     if ( extrap ) then
        iroff2 = iroff2+1
     else
        iroff1 = iroff1+1
     end if
10   continue
     if ( last > 10 .and. erro12 > errmax ) then
        iroff3 = iroff3+1
     end if
15   continue
     rlist(maxerr) = area1
     rlist(last) = area2
     errbnd = max ( epsabs, epsrel*abs(area) )
     !
     !  Test for roundoff error and eventually set error flag.
     !
     if ( iroff1+iroff2 >= 10 .or. iroff3 >= 20 ) then
        ier = 2
     end if
     if ( iroff2 >= 5 ) then
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
     !  at a point of the integration range.
     !
     if ( max ( abs(a1),abs(b2)) <= (1.0e+00+1.0e+03* epsilon ( a1 ) )* &
          (abs(a2)+1.0e+03* tiny ( a2 ) ) ) then
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
     if ( errsum <= errbnd ) go to 115
     if ( ier /= 0 ) then
        exit
     end if
     if ( last == 2 ) go to 80
     if ( noext ) go to 90
     erlarg = erlarg-erlast
     if ( abs(b1-a1) > small ) then
        erlarg = erlarg+erro12
     end if
     !
     !  Test whether the interval to be bisected next is the
     !  smallest interval.
     !
     if ( .not. extrap ) then
        if ( abs(blist(maxerr)-alist(maxerr)) > small ) go to 90
        extrap = .true.
        nrmax = 2
     end if
     !40  continue
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
           if ( abs(blist(maxerr)-alist(maxerr)) > small ) then
              go to 90
           end if
           nrmax = nrmax+1
        end do
     end if
     !
     !  Perform extrapolation.
     !
     !60  continue
     numrl2 = numrl2+1
     rlist2(numrl2) = area
     call qextr ( numrl2, rlist2, reseps, abseps, res3la, nres )
     ktmin = ktmin+1
     if ( ktmin > 5 .and. abserr < 1.0e-03 * errsum ) then
        ier = 5
     end if
     if ( abseps < abserr ) then
        ktmin = 0
        abserr = abseps
        result = reseps
        correc = erlarg
        ertest = max ( epsabs,epsrel*abs(reseps))
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
     go to 90
80   continue
     small = abs ( b - a ) * 3.75e-01
     erlarg = errsum
     ertest = errbnd
     rlist2(2) = area
90   continue
  end do
  !
  !  Set final result and error estimate.
  !
  if ( abserr == huge ( abserr ) ) then
     go to 115
  end if
  if ( ier + ierro == 0 ) then
     go to 110
  end if
  if ( ierro == 3 ) then
     abserr = abserr + correc
  end if
  if ( ier == 0 ) then
     ier = 3
  end if
  if ( result /= 0.0e+00.and.area /= 0.0e+00 ) then
     go to 105
  end if
  if ( abserr > errsum ) go to 115
  if ( area == 0.0e+00 ) go to 130
  go to 110
105 continue
  if ( abserr/abs(result) > errsum/abs(area) ) go to 115
  !
  !  Test on divergence.
  !
110 continue
  if ( ksgn == (-1).and.max ( abs(result),abs(area)) <=  &
       defabs*1.0e-02 ) go to 130
  if ( 1.0e-02 > (result/area) .or. (result/area) > 1.0e+02 &
       .or. errsum > abs(area) ) then
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
  if ( 2 < ier ) then
     ier = ier - 1
  end if
140 continue
  neval = 42*last-21
  return
end subroutine qags
