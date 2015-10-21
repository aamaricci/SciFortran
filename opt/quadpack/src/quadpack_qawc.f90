subroutine qawc ( f, a, b, c, epsabs, epsrel, result, abserr, neval, ier )
  !*****************************************************************************80
  !! QAWC computes a Cauchy principal value.
  !  Discussion:
  !
  !    The routine calculates an approximation RESULT to a Cauchy principal
  !    value 
  !      I = integral of F*W over (A,B),
  !    with
  !      W(X) = 1 / (X-C),
  !    with C distinct from A and B, hopefully satisfying
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
  !    Input, real(8) :: C, a parameter in the weight function, which must
  !    not be equal to A or B.
  !
  !    Input, real(8) :: EPSABS, EPSREL, the absolute and relative accuracy requested.
  !
  !    Output, real(8) :: RESULT, the estimated value of the integral.
  !
  !    Output, real(8) :: ABSERR, an estimate of || I - RESULT ||.
  !
  !    Output, integer NEVAL, the number of times the integral was evaluated.
  !
  !            ier    - integer
  !                     ier = 0 normal and reliable termination of the
  !                             routine. it is assumed that the requested
  !                             accuracy has been achieved.
  !                     ier > 0 abnormal termination of the routine
  !                             the estimates for integral and error are
  !                             less reliable. it is assumed that the
  !                             requested accuracy has not been achieved.
  !                     ier = 1 maximum number of subdivisions allowed
  !                             has been achieved. one can allow more sub-
  !                             divisions by increasing the data value of
  !                             limit in qawc (and taking the according
  !                             dimension adjustments into account).
  !                             however, if this yields no improvement it
  !                             is advised to analyze the integrand in
  !                             order to determine the integration
  !                             difficulties. if the position of a local
  !                             difficulty can be determined (e.g.
  !                             singularity, discontinuity within the
  !                             interval one will probably gain from
  !                             splitting up the interval at this point
  !                             and calling appropriate integrators on the
  !                             subranges.
  !                         = 2 the occurrence of roundoff error is detec-
  !                             ted, which prevents the requested
  !                             tolerance from being achieved.
  !                         = 3 extremely bad integrand behavior occurs
  !                             at some points of the integration
  !                             interval.
  !                         = 6 the input is invalid, because
  !                             c = a or c = b or
  !                             epsabs < 0 and epsrel < 0,
  !                             result, abserr, neval are set to zero.
  !
  !  Local parameters:
  !
  !    LIMIT is the maximum number of subintervals allowed in the
  !    subdivision process of qawce. take care that limit >= 1.
  !
  implicit none
  integer, parameter :: limit = 500
  real(8) :: a
  real(8) :: abserr
  real(8) :: alist(limit)
  real(8) :: b
  real(8) :: blist(limit)
  real(8) :: elist(limit)
  real(8) :: c
  real(8) :: epsabs
  real(8) :: epsrel
  real(8), external :: f
  integer ier
  integer iord(limit)
  integer last
  integer neval
  real(8) :: result
  real(8) :: rlist(limit)
  call qawce ( f, a, b, c, epsabs, epsrel, limit, result, abserr, neval, ier, &
       alist, blist, rlist, elist, iord, last )
  return
end subroutine qawc



subroutine qawce ( f, a, b, c, epsabs, epsrel, limit, result, abserr, neval, &
     ier, alist, blist, rlist, elist, iord, last )
  !*****************************************************************************80
  !! QAWCE computes a Cauchy principal value.
  !  Discussion:
  !
  !    The routine calculates an approximation RESULT to a Cauchy principal
  !    value   
  !      I = integral of F*W over (A,B),
  !    with
  !      W(X) = 1 / ( X - C ),
  !    with C distinct from A and B, hopefully satisfying
  !      | I - RESULT | <= max ( EPSABS, EPSREL * |I| ).
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
  !    Input, real(8) :: C, a parameter in the weight function, which cannot be
  !    equal to A or B.
  !
  !    Input, real(8) :: EPSABS, EPSREL, the absolute and relative accuracy requested.
  !
  !    Input, integer LIMIT, the upper bound on the number of subintervals that
  !    will be used in the partition of [A,B].  LIMIT is typically 500.
  !
  !    Output, real(8) :: RESULT, the estimated value of the integral.
  !
  !    Output, real(8) :: ABSERR, an estimate of || I - RESULT ||.
  !
  !    Output, integer NEVAL, the number of times the integral was evaluated.
  !
  !            ier    - integer
  !                     ier = 0 normal and reliable termination of the
  !                             routine. it is assumed that the requested
  !                             accuracy has been achieved.
  !                     ier > 0 abnormal termination of the routine
  !                             the estimates for integral and error are
  !                             less reliable. it is assumed that the
  !                             requested accuracy has not been achieved.
  !                     ier = 1 maximum number of subdivisions allowed
  !                             has been achieved. one can allow more sub-
  !                             divisions by increasing the value of
  !                             limit. however, if this yields no
  !                             improvement it is advised to analyze the
  !                             integrand, in order to determine the
  !                             integration difficulties.  if the position
  !                             of a local difficulty can be determined
  !                             (e.g. singularity, discontinuity within
  !                             the interval) one will probably gain
  !                             from splitting up the interval at this
  !                             point and calling appropriate integrators
  !                             on the subranges.
  !                         = 2 the occurrence of roundoff error is detec-
  !                             ted, which prevents the requested
  !                             tolerance from being achieved.
  !                         = 3 extremely bad integrand behavior occurs
  !                             at some interior points of the integration
  !                             interval.
  !                         = 6 the input is invalid, because
  !                             c = a or c = b or
  !                             epsabs < 0 and epsrel < 0,
  !                             or limit < 1.
  !                             result, abserr, neval, rlist(1), elist(1),
  !                             iord(1) and last are set to zero.
  !                             alist(1) and blist(1) are set to a and b
  !                             respectively.
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
  !            iord    - integer
  !                      vector of dimension at least limit, the first k
  !                      elements of which are pointers to the error
  !                      estimates over the subintervals, so that
  !                      elist(iord(1)), ...,  elist(iord(k)) with
  !                      k = last if last <= (limit/2+2), and
  !                      k = limit+1-last otherwise, form a decreasing
  !                      sequence.
  !
  !            last    - integer
  !                      number of subintervals actually produced in
  !                      the subdivision process
  !
  !  Local parameters:
  !
  !           alist     - list of left end points of all subintervals
  !                       considered up to now
  !           blist     - list of right end points of all subintervals
  !                       considered up to now
  !           rlist(i)  - approximation to the integral over
  !                       (alist(i),blist(i))
  !           elist(i)  - error estimate applying to rlist(i)
  !           maxerr    - pointer to the interval with largest error
  !                       estimate
  !           errmax    - elist(maxerr)
  !           area      - sum of the integrals over the subintervals
  !           errsum    - sum of the errors over the subintervals
  !           errbnd    - requested accuracy max(epsabs,epsrel*
  !                       abs(result))
  !           *****1    - variable for the left subinterval
  !           *****2    - variable for the right subinterval
  !           last      - index for subdivision
  !
  implicit none
  integer limit
  real(8) :: a
  real(8) :: aa
  real(8) :: abserr
  real(8) :: alist(limit)
  real(8) :: area
  real(8) :: area1
  real(8) :: area12
  real(8) :: area2
  real(8) :: a1
  real(8) :: a2
  real(8) :: b
  real(8) :: bb
  real(8) :: blist(limit)
  real(8) :: b1
  real(8) :: b2
  real(8) :: c
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
  integer krule
  integer last
  integer maxerr
  integer nev
  integer neval
  integer nrmax
  real(8) :: result
  real(8) :: rlist(limit)
  !
  !  Test on validity of parameters.
  !
  ier = 0
  neval = 0
  last = 0
  alist(1) = a
  blist(1) = b
  rlist(1) = 0.0e+00
  elist(1) = 0.0e+00
  iord(1) = 0
  result = 0.0e+00
  abserr = 0.0e+00
  if ( c == a  ) then
     ier = 6
     return
  else if ( c == b ) then
     ier = 6
     return
  else if ( epsabs < 0.0e+00 .and. epsrel < 0.0e+00 ) then
     ier = 6
     return
  end if
  !
  !  First approximation to the integral.
  !
  if ( a <= b ) then
     aa = a
     bb = b
  else
     aa = b
     bb = a
  end if
  krule = 1
  call qc25c ( f, aa, bb, c, result, abserr, krule, neval )
  last = 1
  rlist(1) = result
  elist(1) = abserr
  iord(1) = 1
  alist(1) = a
  blist(1) = b
  !
  !  Test on accuracy.
  !
  errbnd = max ( epsabs, epsrel * abs(result) )
  if ( limit == 1 ) then
     ier = 1
     go to 70
  end if
  if ( abserr < min ( 1.0e-02 * abs(result), errbnd)  ) then
     go to 70
  end if
  !
  !  Initialization
  !
  alist(1) = aa
  blist(1) = bb
  rlist(1) = result
  errmax = abserr
  maxerr = 1
  area = result
  errsum = abserr
  nrmax = 1
  iroff1 = 0
  iroff2 = 0
  do last = 2, limit
     !
     !  Bisect the subinterval with nrmax-th largest error estimate.
     !
     a1 = alist(maxerr)
     b1 = 5.0e-01*(alist(maxerr)+blist(maxerr))
     b2 = blist(maxerr)
     if ( c <= b1 .and. a1 < c ) then
        b1 = 5.0e-01*(c+b2)
     end if
     if ( b1 < c .and. c < b2 ) then
        b1 = 5.0e-01 * ( a1 + c )
     end if
     a2 = b1
     krule = 2
     call qc25c ( f, a1, b1, c, area1, error1, krule, nev )
     neval = neval+nev
     call qc25c ( f, a2, b2, c, area2, error2, krule, nev )
     neval = neval+nev
     !
     !  Improve previous approximations to integral and error
     !  and test for accuracy.
     !
     area12 = area1 + area2
     erro12 = error1 + error2
     errsum = errsum + erro12 - errmax
     area = area + area12 - rlist(maxerr)
     if ( abs ( rlist(maxerr)-area12) < 1.0e-05 * abs(area12) &
          .and. erro12 >= 9.9e-01 * errmax .and. krule == 0 ) &
          iroff1 = iroff1+1
     if ( last > 10.and.erro12 > errmax .and. krule == 0 ) then
        iroff2 = iroff2+1
     end if
     rlist(maxerr) = area1
     rlist(last) = area2
     errbnd = max ( epsabs, epsrel * abs(area) )
     if ( errsum > errbnd ) then
        !
        !  Test for roundoff error and eventually set error flag.
        !
        if ( iroff1 >= 6 .and. iroff2 > 20 ) then
           ier = 2
        end if
        !
        !  Set error flag in the case that number of interval
        !  bisections exceeds limit.
        !
        if ( last == limit ) then
           ier = 1
        end if
        !
        !  Set error flag in the case of bad integrand behavior at
        !  a point of the integration range.
        !
        if ( max ( abs(a1), abs(b2) ) <= ( 1.0e+00 + 1.0e+03 * epsilon ( a1 ) ) &
             *( abs(a2)+1.0e+03* tiny ( a2 ) )) then
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
     !  with NRMAX-th largest error estimate (to be bisected next).
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
70 continue 
  if ( aa == b ) then
     result = - result
  end if
  return
end subroutine qawce
