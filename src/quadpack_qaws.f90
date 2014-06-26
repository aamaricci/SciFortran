subroutine qaws ( f, a, b, alfa, beta, integr, epsabs, epsrel, result, &
     abserr, neval, ier )
  !*****************************************************************************80
  !! QAWS estimates integrals with algebraico-logarithmic endpoint singularities.
  !  Discussion:
  !
  !    This routine calculates an approximation RESULT to a given
  !    definite integral   
  !      I = integral of f*w over (a,b) 
  !    where w shows a singular behavior at the end points, see parameter
  !    integr, hopefully satisfying following claim for accuracy
  !      abs(i-result) <= max(epsabs,epsrel*abs(i)).
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
  !    Input, real(8) :: ALFA, BETA, parameters used in the weight function.
  !    ALFA and BETA should be greater than -1.
  !
  !    Input, integer INTEGR, indicates which weight function is to be used
  !    = 1  (x-a)**alfa*(b-x)**beta
  !    = 2  (x-a)**alfa*(b-x)**beta*log(x-a)
  !    = 3  (x-a)**alfa*(b-x)**beta*log(b-x)
  !    = 4  (x-a)**alfa*(b-x)**beta*log(x-a)*log(b-x)
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
  !                             the estimates for the integral and error
  !                             are less reliable. it is assumed that the
  !                             requested accuracy has not been achieved.
  !                     ier = 1 maximum number of subdivisions allowed
  !                             has been achieved. one can allow more
  !                             subdivisions by increasing the data value
  !                             of limit in qaws (and taking the according
  !                             dimension adjustments into account).
  !                             however, if this yields no improvement it
  !                             is advised to analyze the integrand, in
  !                             order to determine the integration
  !                             difficulties which prevent the requested
  !                             tolerance from being achieved. in case of
  !                             a jump discontinuity or a local
  !                             singularity of algebraico-logarithmic type
  !                             at one or more interior points of the
  !                             integration range, one should proceed by
  !                             splitting up the interval at these points
  !                             and calling the integrator on the
  !                             subranges.
  !                         = 2 the occurrence of roundoff error is
  !                             detected, which prevents the requested
  !                             tolerance from being achieved.
  !                         = 3 extremely bad integrand behavior occurs
  !                             at some points of the integration
  !                             interval.
  !                         = 6 the input is invalid, because
  !                             b <= a or alfa <= (-1) or beta <= (-1) or
  !                             integr < 1 or integr > 4 or
  !                             epsabs < 0 and epsrel < 0,
  !                             result, abserr, neval are set to zero.
  !
  !  Local parameters:
  !
  !    LIMIT is the maximum number of subintervals allowed in the
  !    subdivision process of qawse. take care that limit >= 2.
  !
  implicit none
  integer, parameter :: limit = 500
  real(8) :: a
  real(8) :: abserr
  real(8) :: alfa
  real(8) :: alist(limit)
  real(8) :: b
  real(8) :: blist(limit)
  real(8) :: beta
  real(8) :: elist(limit)
  real(8) :: epsabs
  real(8) :: epsrel
  real(8), external :: f
  integer ier
  integer integr
  integer iord(limit)
  integer last
  integer neval
  real(8) :: result
  real(8) :: rlist(limit)
  call qawse ( f, a, b, alfa, beta, integr, epsabs, epsrel, limit, result, &
       abserr, neval, ier, alist, blist, rlist, elist, iord, last )
  return
end subroutine qaws



subroutine qawse ( f, a, b, alfa, beta, integr, epsabs, epsrel, limit, &
     result, abserr, neval, ier, alist, blist, rlist, elist, iord, last )
  !*****************************************************************************80
  !! QAWSE estimates integrals with algebraico-logarithmic endpoint singularities.
  !  Discussion:
  !    This routine calculates an approximation RESULT to an integral
  !      I = integral of F(X) * W(X) over (a,b), 
  !    where W(X) shows a singular behavior at the endpoints, hopefully 
  !    satisfying:
  !      | I - RESULT | <= max ( epsabs, epsrel * |I| ).
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
  !    Input, real(8) :: ALFA, BETA, parameters used in the weight function.
  !    ALFA and BETA should be greater than -1.
  !
  !    Input, integer INTEGR, indicates which weight function is used:
  !    = 1  (x-a)**alfa*(b-x)**beta
  !    = 2  (x-a)**alfa*(b-x)**beta*log(x-a)
  !    = 3  (x-a)**alfa*(b-x)**beta*log(b-x)
  !    = 4  (x-a)**alfa*(b-x)**beta*log(x-a)*log(b-x)
  !
  !    Input, real(8) :: EPSABS, EPSREL, the absolute and relative accuracy requested.
  !
  !    Input, integer LIMIT, an upper bound on the number of subintervals
  !    in the partition of (A,B), LIMIT >= 2.  If LIMIT < 2, the routine 
  !     will end with IER = 6.
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
  !                             the estimates for the integral and error
  !                             are less reliable. it is assumed that the
  !                             requested accuracy has not been achieved.
  !                         = 1 maximum number of subdivisions allowed
  !                             has been achieved. one can allow more
  !                             subdivisions by increasing the value of
  !                             limit. however, if this yields no
  !                             improvement it is advised to analyze the
  !                             integrand, in order to determine the
  !                             integration difficulties which prevent
  !                             the requested tolerance from being
  !                             achieved. in case of a jump discontinuity
  !                             or a local singularity of algebraico-
  !                             logarithmic type at one or more interior
  !                             points of the integration range, one
  !                             should proceed by splitting up the
  !                             interval at these points and calling the
  !                             integrator on the subranges.
  !                         = 2 the occurrence of roundoff error is
  !                             detected, which prevents the requested
  !                             tolerance from being achieved.
  !                         = 3 extremely bad integrand behavior occurs
  !                             at some points of the integration
  !                             interval.
  !                         = 6 the input is invalid, because
  !                             b <= a or alfa <= (-1) or beta <= (-1) or
  !                             integr < 1 or integr > 4, or
  !                             epsabs < 0 and epsrel < 0,
  !                             or limit < 2.
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
  !            iord   - integer
  !                     vector of dimension at least limit, the first k
  !                     elements of which are pointers to the error
  !                     estimates over the subintervals, so that
  !                     elist(iord(1)), ..., elist(iord(k)) with k = last
  !                     if last <= (limit/2+2), and k = limit+1-last
  !                     otherwise, form a decreasing sequence.
  !
  !    Output, integer LAST, the number of subintervals actually produced in 
  !    the subdivision process.
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
  real(8) :: abserr
  real(8) :: alfa
  real(8) :: alist(limit)
  real(8) :: area
  real(8) :: area1
  real(8) :: area12
  real(8) :: area2
  real(8) :: a1
  real(8) :: a2
  real(8) :: b
  real(8) :: beta
  real(8) :: blist(limit)
  real(8) :: b1
  real(8) :: b2
  real(8) :: centre
  real(8) :: elist(limit)
  real(8) :: epsabs
  real(8) :: epsrel
  real(8) :: errbnd
  real(8) :: errmax
  real(8) :: error1
  real(8) :: erro12
  real(8) :: error2
  real(8) :: errsum
  real(8), external :: f
  integer ier
  integer integr
  integer iord(limit)
  integer iroff1
  integer iroff2
  integer last
  integer maxerr
  integer nev
  integer neval
  integer nrmax
  real(8) :: resas1
  real(8) :: resas2
  real(8) :: result
  real(8) :: rg(25)
  real(8) :: rh(25)
  real(8) :: ri(25)
  real(8) :: rj(25)
  real(8) :: rlist(limit)
  !
  !  Test on validity of parameters.
  !
  ier = 0
  neval = 0
  last = 0
  rlist(1) = 0.0e+00
  elist(1) = 0.0e+00
  iord(1) = 0
  result = 0.0e+00
  abserr = 0.0e+00
  if ( b <= a .or. &
       (epsabs < 0.0e+00 .and. epsrel < 0.0e+00) .or. &
       alfa <= (-1.0e+00) .or. &
       beta <= (-1.0e+00) .or. &
       integr < 1 .or. &
       integr > 4 .or. &
       limit < 2 ) then
     ier = 6
     return
  end if
  !
  !  Compute the modified Chebyshev moments.
  !
  call qmomo ( alfa, beta, ri, rj, rg, rh, integr )
  !
  !  Integrate over the intervals (a,(a+b)/2) and ((a+b)/2,b).
  !
  centre = 5.0e-01 * ( b + a )
  call qc25s ( f, a, b, a, centre, alfa, beta, ri, rj, rg, rh, area1, &
       error1, resas1, integr, nev )
  neval = nev
  call qc25s ( f, a, b, centre, b, alfa, beta, ri, rj, rg, rh, area2, &
       error2, resas2, integr, nev )
  last = 2
  neval = neval+nev
  result = area1+area2
  abserr = error1+error2
  !
  !  Test on accuracy.
  !
  errbnd = max ( epsabs, epsrel * abs ( result ) )
  !
  !  Initialization.
  !
  if ( error2 <= error1 ) then
     alist(1) = a
     alist(2) = centre
     blist(1) = centre
     blist(2) = b
     rlist(1) = area1
     rlist(2) = area2
     elist(1) = error1
     elist(2) = error2
  else
     alist(1) = centre
     alist(2) = a
     blist(1) = b
     blist(2) = centre
     rlist(1) = area2
     rlist(2) = area1
     elist(1) = error2
     elist(2) = error1
  end if
  iord(1) = 1
  iord(2) = 2
  if ( limit == 2 ) then
     ier = 1
     return
  end if
  if ( abserr <= errbnd ) then
     return
  end if
  errmax = elist(1)
  maxerr = 1
  nrmax = 1
  area = result
  errsum = abserr
  iroff1 = 0
  iroff2 = 0
  do last = 3, limit
     !
     !  Bisect the subinterval with largest error estimate.
     !
     a1 = alist(maxerr)
     b1 = 5.0e-01 * ( alist(maxerr) + blist(maxerr) )
     a2 = b1
     b2 = blist(maxerr)
     call qc25s ( f, a, b, a1, b1, alfa, beta, ri, rj, rg, rh, area1, &
          error1, resas1, integr, nev )
     neval = neval + nev
     call qc25s ( f, a, b, a2, b2, alfa, beta, ri, rj, rg, rh, area2, &
          error2, resas2, integr, nev )
     neval = neval + nev
     !
     !  Improve previous approximations integral and error and
     !  test for accuracy.
     !
     area12 = area1+area2
     erro12 = error1+error2
     errsum = errsum+erro12-errmax
     area = area+area12-rlist(maxerr)
     !
     !  Test for roundoff error.
     !
     if ( a /= a1 .and. b /= b2 ) then
        if ( resas1 /= error1 .and. resas2 /= error2 ) then
           if ( abs ( rlist(maxerr) - area12 ) < 1.0e-05 * abs ( area12 ) &
                .and.erro12 >= 9.9e-01*errmax ) then
              iroff1 = iroff1 + 1
           end if
           if ( last > 10 .and. erro12 > errmax ) then
              iroff2 = iroff2 + 1
           end if
        end if
     end if
     rlist(maxerr) = area1
     rlist(last) = area2
     !
     !  Test on accuracy.
     !
     errbnd = max ( epsabs, epsrel * abs ( area ) )
     if ( errsum > errbnd ) then
        !
        !  Set error flag in the case that the number of interval
        !  bisections exceeds limit.
        !
        if ( last == limit ) then
           ier = 1
        end if
        !
        !  Set error flag in the case of roundoff error.
        !
        if ( iroff1 >= 6 .or. iroff2 >= 20 ) then
           ier = 2
        end if
        !
        !  Set error flag in the case of bad integrand behavior
        !  at interior points of integration range.
        !
        if ( max ( abs(a1),abs(b2)) <= (1.0e+00+1.0e+03* epsilon ( a1 ) )* &
             ( abs(a2) + 1.0e+03* tiny ( a2) )) then
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
     !  with largest error estimate (to be bisected next).
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
  return
end subroutine qawse
