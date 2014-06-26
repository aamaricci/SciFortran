subroutine qawf ( f, a, omega, integr, epsabs, result, abserr, neval, ier )
  !*****************************************************************************80
  !! QAWF computes Fourier integrals over the interval [ A, +Infinity ).
  !  Discussion:
  !
  !    The routine calculates an approximation RESULT to a definite integral  
  ! 
  !      I = integral of F*COS(OMEGA*X) 
  !    or 
  !      I = integral of F*SIN(OMEGA*X) 
  !
  !    over the interval [A,+Infinity), hopefully satisfying
  !
  !      || I - RESULT || <= EPSABS.
  !
  !    If OMEGA = 0 and INTEGR = 1, the integral is calculated by means 
  !    of QAGI, and IER has the meaning as described in the comments of QAGI.
  !
  !  Parameters:
  !
  !    Input, external real(8) :: F, the name of the function routine, of the form
  !      function f ( x )
  !      real(8) :: f
  !      real(8) :: x
  !    which evaluates the integrand function.
  !
  !    Input, real(8) :: A, the lower limit of integration.
  !
  !    Input, real(8) :: OMEGA, the parameter in the weight function.
  !
  !    Input, integer INTEGR, indicates which weight functions is used
  !    = 1, w(x) = cos(omega*x)
  !    = 2, w(x) = sin(omega*x)
  !
  !    Input, real(8) :: EPSABS, the absolute accuracy requested.
  !
  !    Output, real(8) :: RESULT, the estimated value of the integral.
  !
  !    Output, real(8) :: ABSERR, an estimate of || I - RESULT ||.
  !
  !    Output, integer NEVAL, the number of times the integral was evaluated.
  !
  !            ier    - integer
  !                     ier = 0 normal and reliable termination of the
  !                             routine. it is assumed that the
  !                             requested accuracy has been achieved.
  !                     ier > 0 abnormal termination of the routine.
  !                             the estimates for integral and error are
  !                             less reliable. it is assumed that the
  !                             requested accuracy has not been achieved.
  !                    if omega /= 0
  !                     ier = 6 the input is invalid because
  !                             (integr /= 1 and integr /= 2) or
  !                              epsabs <= 0
  !                              result, abserr, neval, lst are set to
  !                              zero.
  !                         = 7 abnormal termination of the computation
  !                             of one or more subintegrals
  !                         = 8 maximum number of cycles allowed
  !                             has been achieved, i.e. of subintervals
  !                             (a+(k-1)c,a+kc) where
  !                             c = (2*int(abs(omega))+1)*pi/abs(omega),
  !                             for k = 1, 2, ...
  !                         = 9 the extrapolation table constructed for
  !                             convergence acceleration of the series
  !                             formed by the integral contributions
  !                             over the cycles, does not converge to
  !                             within the requested accuracy.
  !
  !  Local parameters:
  !
  !    Integer LIMLST, gives an upper bound on the number of cycles, LIMLST >= 3.
  !    if limlst < 3, the routine will end with ier = 6.
  !
  !    Integer MAXP1, an upper bound on the number of Chebyshev moments which 
  !    can be stored, i.e. for the intervals of lengths abs(b-a)*2**(-l), 
  !    l = 0,1, ..., maxp1-2, maxp1 >= 1.  if maxp1 < 1, the routine will end
  !    with ier = 6.
  !
  implicit none
  integer, parameter :: limit = 500
  integer, parameter :: limlst = 50
  integer, parameter :: maxp1 = 21
  real(8) :: a
  real(8) :: abserr
  real(8) :: alist(limit)
  real(8) :: blist(limit)
  real(8) :: chebmo(maxp1,25)
  real(8) :: elist(limit)
  real(8) :: epsabs
  real(8) :: erlst(limlst)
  real(8), external :: f
  integer ier
  integer integr
  integer iord(limit)
  integer ierlst(limlst)
  integer lst
  integer neval
  integer nnlog(limit)
  real(8) :: omega
  real(8) :: result
  real(8) :: rlist(limit)
  real(8) :: rslst(limlst)
  ier = 6
  neval = 0
  result = 0.0e+00
  abserr = 0.0e+00
  if ( limlst < 3 .or. maxp1 < 1 ) then
     return
  end if
  call qawfe ( f, a, omega, integr, epsabs, limlst, limit, maxp1, &
       result, abserr, neval, ier, rslst, erlst, ierlst, lst, alist, blist, &
       rlist, elist, iord, nnlog, chebmo )
  return
end subroutine qawf




subroutine qawfe ( f, a, omega, integr, epsabs, limlst, limit, maxp1, &
     result, abserr, neval, ier, rslst, erlst, ierlst, lst, alist, blist, &
     rlist, elist, iord, nnlog, chebmo )
  !*****************************************************************************80
  !! QAWFE computes Fourier integrals.
  !  Discussion:
  !
  !    The routine calculates an approximation RESULT to a definite integral   
  !      I = integral of F*COS(OMEGA*X) or F*SIN(OMEGA*X) over (A,+Infinity),
  !    hopefully satisfying
  !      || I - RESULT || <= EPSABS.
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
  !    Input, real(8) :: A, the lower limit of integration.
  !
  !    Input, real(8) :: OMEGA, the parameter in the weight function.
  !
  !    Input, integer INTEGR, indicates which weight function is used
  !    = 1      w(x) = cos(omega*x)
  !    = 2      w(x) = sin(omega*x)
  !
  !    Input, real(8) :: EPSABS, the absolute accuracy requested.
  !
  !    Input, integer LIMLST, an upper bound on the number of cycles.
  !    LIMLST must be at least 1.  In fact, if LIMLST < 3, the routine 
  !    will end with IER= 6.
  !
  !    Input, integer LIMIT, an upper bound on the number of subintervals 
  !    allowed in the partition of each cycle, limit >= 1.
  !
  !            maxp1  - integer
  !                     gives an upper bound on the number of
  !                     Chebyshev moments which can be stored, i.e.
  !                     for the intervals of lengths abs(b-a)*2**(-l),
  !                     l=0,1, ..., maxp1-2, maxp1 >= 1
  !
  !    Output, real(8) :: RESULT, the estimated value of the integral.
  !
  !    Output, real(8) :: ABSERR, an estimate of || I - RESULT ||.
  !
  !    Output, integer NEVAL, the number of times the integral was evaluated.
  !
  !            ier    - ier = 0 normal and reliable termination of
  !                             the routine. it is assumed that the
  !                             requested accuracy has been achieved.
  !                     ier > 0 abnormal termination of the routine
  !                             the estimates for integral and error
  !                             are less reliable. it is assumed that
  !                             the requested accuracy has not been
  !                             achieved.
  !                    if omega /= 0
  !                     ier = 6 the input is invalid because
  !                             (integr /= 1 and integr /= 2) or
  !                              epsabs <= 0 or limlst < 3.
  !                              result, abserr, neval, lst are set
  !                              to zero.
  !                         = 7 bad integrand behavior occurs within
  !                             one or more of the cycles. location
  !                             and type of the difficulty involved
  !                             can be determined from the vector ierlst.
  !                             here lst is the number of cycles actually
  !                             needed (see below).
  !                             ierlst(k) = 1 the maximum number of
  !                                           subdivisions (= limit)
  !                                           has been achieved on the
  !                                           k th cycle.
  !                                       = 2 occurence of roundoff
  !                                           error is detected and
  !                                           prevents the tolerance
  !                                           imposed on the k th cycle
  !                                           from being acheived.
  !                                       = 3 extremely bad integrand
  !                                           behavior occurs at some
  !                                           points of the k th cycle.
  !                                       = 4 the integration procedure
  !                                           over the k th cycle does
  !                                           not converge (to within the
  !                                           required accuracy) due to
  !                                           roundoff in the
  !                                           extrapolation procedure
  !                                           invoked on this cycle. it
  !                                           is assumed that the result
  !                                           on this interval is the
  !                                           best which can be obtained.
  !                                       = 5 the integral over the k th
  !                                           cycle is probably divergent
  !                                           or slowly convergent. it
  !                                           must be noted that
  !                                           divergence can occur with
  !                                           any other value of
  !                                           ierlst(k).
  !                         = 8 maximum number of  cycles  allowed
  !                             has been achieved, i.e. of subintervals
  !                             (a+(k-1)c,a+kc) where
  !                             c = (2*int(abs(omega))+1)*pi/abs(omega),
  !                             for k = 1, 2, ..., lst.
  !                             one can allow more cycles by increasing
  !                             the value of limlst (and taking the
  !                             according dimension adjustments into
  !                             account).
  !                             examine the array iwork which contains
  !                             the error flags over the cycles, in order
  !                             to eventual look for local integration
  !                             difficulties.
  !                             if the position of a local difficulty can
  !                             be determined (e.g. singularity,
  !                             discontinuity within the interval)
  !                             one will probably gain from splitting
  !                             up the interval at this point and
  !                             calling appopriate integrators on the
  !                             subranges.
  !                         = 9 the extrapolation table constructed for
  !                             convergence acceleration of the series
  !                             formed by the integral contributions
  !                             over the cycles, does not converge to
  !                             within the required accuracy.
  !                             as in the case of ier = 8, it is advised
  !                             to examine the array iwork which contains
  !                             the error flags on the cycles.
  !                    if omega = 0 and integr = 1,
  !                    the integral is calculated by means of qagi
  !                    and ier = ierlst(1) (with meaning as described
  !                    for ierlst(k), k = 1).
  !
  !            rslst  - real
  !                     vector of dimension at least limlst
  !                     rslst(k) contains the integral contribution
  !                     over the interval (a+(k-1)c,a+kc) where
  !                     c = (2*int(abs(omega))+1)*pi/abs(omega),
  !                     k = 1, 2, ..., lst.
  !                     note that, if omega = 0, rslst(1) contains
  !                     the value of the integral over (a,infinity).
  !
  !            erlst  - real
  !                     vector of dimension at least limlst
  !                     erlst(k) contains the error estimate
  !                     corresponding with rslst(k).
  !
  !            ierlst - integer
  !                     vector of dimension at least limlst
  !                     ierlst(k) contains the error flag corresponding
  !                     with rslst(k). for the meaning of the local error
  !                     flags see description of output parameter ier.
  !
  !            lst    - integer
  !                     number of subintervals needed for the integration
  !                     if omega = 0 then lst is set to 1.
  !
  !            alist, blist, rlist, elist - real
  !                     vector of dimension at least limit,
  !
  !            iord, nnlog - integer
  !                     vector of dimension at least limit, providing
  !                     space for the quantities needed in the
  !                     subdivision process of each cycle
  !
  !            chebmo - real
  !                     array of dimension at least (maxp1,25),
  !                     providing space for the Chebyshev moments
  !                     needed within the cycles
  !
  !  Local parameters:
  !
  !           c1, c2    - end points of subinterval (of length
  !                       cycle)
  !           cycle     - (2*int(abs(omega))+1)*pi/abs(omega)
  !           psum      - vector of dimension at least (limexp+2)
  !                       (see routine qextr)
  !                       psum contains the part of the epsilon table
  !                       which is still needed for further computations.
  !                       each element of psum is a partial sum of
  !                       the series which should sum to the value of
  !                       the integral.
  !           errsum    - sum of error estimates over the
  !                       subintervals, calculated cumulatively
  !           epsa      - absolute tolerance requested over current
  !                       subinterval
  !           chebmo    - array containing the modified Chebyshev
  !                       moments (see also routine qc25o)
  !
  implicit none
  integer limit
  integer limlst
  integer maxp1
  real(8) :: a
  real(8) :: abseps
  real(8) :: abserr
  real(8) :: alist(limit)
  real(8) :: blist(limit)
  real(8) :: chebmo(maxp1,25)
  real(8) :: correc
  real(8) :: cycle
  real(8) :: c1
  real(8) :: c2
  real(8) :: dl
  ! real(8) :: dla
  real(8) :: drl
  real(8) :: elist(limit)
  real(8) :: ep
  real(8) :: eps
  real(8) :: epsa
  real(8) :: epsabs
  real(8) :: erlst(limlst)
  real(8) :: errsum
  real(8), external :: f
  real(8) :: fact
  integer ier
  integer ierlst(limlst)
  integer integr
  integer iord(limit)
  integer ktmin
  integer l
  integer ll
  integer lst
  integer momcom
  integer nev
  integer neval
  integer nnlog(limit)
  integer nres
  integer numrl2
  real(8) :: omega
  real(8), parameter :: p = 0.9E+00
  real(8), parameter :: pi = 3.1415926535897932E+00
  real(8) :: p1
  real(8) :: psum(52)
  real(8) :: reseps
  real(8) :: result
  real(8) :: res3la(3)
  real(8) :: rlist(limit)
  real(8) :: rslst(limlst)
  !
  !  The dimension of  psum  is determined by the value of
  !  limexp in QEXTR (psum must be
  !  of dimension (limexp+2) at least).
  !
  !  Test on validity of parameters.
  !
  result = 0.0e+00
  abserr = 0.0e+00
  neval = 0
  lst = 0
  ier = 0
  if ( (integr /= 1 .and. integr /= 2 ) .or. &
       epsabs <= 0.0e+00 .or. &
       limlst < 3 ) then
     ier = 6
     return
  end if
  if ( omega == 0.0e+00 ) then
     if ( integr == 1 ) then
        call qagi ( f, 0.d+00, 1, epsabs, 0.d+00, result, abserr, neval, ier )
     else
        result = 0.0E+00
        abserr = 0.0E+00
        neval = 0
        ier = 0
     end if
     rslst(1) = result
     erlst(1) = abserr
     ierlst(1) = ier
     lst = 1
     return
  end if
  !
  !  Initializations.
  !
  l = int ( abs ( omega ) )
  dl = 2 * l + 1
  cycle = dl * pi / abs ( omega )
  ier = 0
  ktmin = 0
  neval = 0
  numrl2 = 0
  nres = 0
  c1 = a
  c2 = cycle+a
  p1 = 1.0e+00-p
  eps = epsabs
  if ( epsabs > tiny ( epsabs ) / p1 ) then
     eps = epsabs * p1
  end if
  ep = eps
  fact = 1.0e+00
  correc = 0.0e+00
  abserr = 0.0e+00
  errsum = 0.0e+00
  do lst = 1, limlst
     !
     !  Integrate over current subinterval.
     !
     !   dla = lst
     epsa = eps * fact
     call qfour ( f, c1, c2, omega, integr, epsa, 0.0d+00, limit, lst, maxp1, &
          rslst(lst), erlst(lst), nev, ierlst(lst), alist, blist, rlist, elist, &
          iord, nnlog, momcom, chebmo )
     neval = neval + nev
     fact = fact * p
     errsum = errsum + erlst(lst)
     drl = 5.0e+01 * abs(rslst(lst))
     !
     !  Test on accuracy with partial sum.
     !
     if ((errsum+drl) <= epsabs.and.lst >= 6) then
        go to 80
     end if
     correc = max ( correc,erlst(lst))
     if ( ierlst(lst) /= 0 ) then
        eps = max ( ep,correc*p1)
        ier = 7
     end if
     if ( ier == 7 .and. (errsum+drl) <= correc*1.0e+01.and. lst > 5) go to 80
     numrl2 = numrl2+1
     if ( lst <= 1 ) then
        psum(1) = rslst(1)
        go to 40
     end if
     psum(numrl2) = psum(ll) + rslst(lst)
     if ( lst == 2 ) then
        go to 40
     end if
     !
     !  Test on maximum number of subintervals
     !
     if ( lst == limlst ) then
        ier = 8
     end if
     !
     !  Perform new extrapolation
     !
     call qextr ( numrl2, psum, reseps, abseps, res3la, nres )
     !
     !  Test whether extrapolated result is influenced by roundoff
     !
     ktmin = ktmin + 1
     if ( ktmin >= 15 .and. abserr <= 1.0e-03 * (errsum+drl) ) then
        ier = 9
     end  if
     if ( abseps <= abserr .or. lst == 3 ) then
        abserr = abseps
        result = reseps
        ktmin = 0
        !
        !  If IER is not 0, check whether direct result (partial
        !  sum) or extrapolated result yields the best integral
        !  approximation
        !
        if ( ( abserr + 1.0e+01 * correc ) <= epsabs ) then
           exit
        end if
        if ( abserr <= epsabs .and. 1.0e+01 * correc >= epsabs ) then
           exit
        end if
     end if
     if ( ier /= 0 .and. ier /= 7 ) then
        exit
     end if
40   continue
     ll = numrl2
     c1 = c2
     c2 = c2+cycle
  end do
  !
  !  Set final result and error estimate.
  !
  !60 continue
  abserr = abserr + 1.0e+01 * correc
  if ( ier == 0 ) then
     return
  end if
  if ( result /= 0.0e+00 .and. psum(numrl2) /= 0.0e+00) go to 70
  if ( abserr > errsum ) then
     go to 80
  end if
  if ( psum(numrl2) == 0.0e+00 ) then
     return
  end if
70 continue
  if ( abserr / abs(result) <= (errsum+drl)/abs(psum(numrl2)) ) then
     if ( ier >= 1 .and. ier /= 7 ) then
        abserr = abserr + drl
     end if
     return
  end if
80 continue
  result = psum(numrl2)
  abserr = errsum + drl
  return
end subroutine qawfe
