
subroutine qawo ( f, a, b, omega, integr, epsabs, epsrel, result, abserr, &
     neval, ier )
  !*****************************************************************************80
  !! QAWO computes the integrals of oscillatory integrands.
  !  Discussion:
  !
  !    The routine calculates an approximation RESULT to a given
  !    definite integral
  !      I = Integral ( A <= X <= B ) F(X) * cos ( OMEGA * X ) dx
  !    or 
  !      I = Integral ( A <= X <= B ) F(X) * sin ( OMEGA * X ) dx
  !    hopefully satisfying following claim for accuracy
  !      | I - RESULT | <= max ( epsabs, epsrel * |I| ).
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
  !    Input, real(8) :: OMEGA, the parameter in the weight function.
  !
  !    Input, integer INTEGR, specifies the weight function:
  !    1, W(X) = cos ( OMEGA * X )
  !    2, W(X) = sin ( OMEGA * X )
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
  !                             routine. it is assumed that the
  !                             requested accuracy has been achieved.
  !                   - ier > 0 abnormal termination of the routine.
  !                             the estimates for integral and error are
  !                             less reliable. it is assumed that the
  !                             requested accuracy has not been achieved.
  !                     ier = 1 maximum number of subdivisions allowed
  !                             (= leniw/2) has been achieved. one can
  !                             allow more subdivisions by increasing the
  !                             value of leniw (and taking the according
  !                             dimension adjustments into account).
  !                             however, if this yields no improvement it
  !                             is advised to analyze the integrand in
  !                             order to determine the integration
  !                             difficulties. if the position of a local
  !                             difficulty can be determined (e.g.
  !                             singularity, discontinuity within the
  !                             interval) one will probably gain from
  !                             splitting up the interval at this point
  !                             and calling the integrator on the
  !                             subranges. if possible, an appropriate
  !                             special-purpose integrator should
  !                             be used which is designed for handling
  !                             the type of difficulty involved.
  !                         = 2 the occurrence of roundoff error is
  !                             detected, which prevents the requested
  !                             tolerance from being achieved.
  !                             the error may be under-estimated.
  !                         = 3 extremely bad integrand behavior occurs
  !                             at some interior points of the integration
  !                             interval.
  !                         = 4 the algorithm does not converge. roundoff
  !                             error is detected in the extrapolation
  !                             table. it is presumed that the requested
  !                             tolerance cannot be achieved due to
  !                             roundoff in the extrapolation table,
  !                             and that the returned result is the best
  !                             which can be obtained.
  !                         = 5 the integral is probably divergent, or
  !                             slowly convergent. it must be noted that
  !                             divergence can occur with any other value
  !                             of ier.
  !                         = 6 the input is invalid, because
  !                             epsabs < 0 and epsrel < 0,
  !                             result, abserr, neval are set to zero.
  !
  !  Local parameters:
  !
  !    limit is the maximum number of subintervals allowed in the
  !    subdivision process of QFOUR. take care that limit >= 1.
  !
  !    maxp1 gives an upper bound on the number of Chebyshev moments
  !    which can be stored, i.e. for the intervals of lengths
  !    abs(b-a)*2**(-l), l = 0, 1, ... , maxp1-2. take care that
  !    maxp1 >= 1.
  implicit none
  integer, parameter :: limit = 500
  integer, parameter :: maxp1 = 21
  real(8) :: a
  real(8) :: abserr
  real(8) :: alist(limit)
  real(8) :: b
  real(8) :: blist(limit)
  real(8) :: chebmo(maxp1,25)
  real(8) :: elist(limit)
  real(8) :: epsabs
  real(8) :: epsrel
  real(8), external :: f
  integer ier
  integer integr
  integer iord(limit)
  integer momcom
  integer neval
  integer nnlog(limit)
  real(8) :: omega
  real(8) :: result
  real(8) :: rlist(limit)
  call qfour ( f, a, b, omega, integr, epsabs, epsrel, limit, 1, maxp1, &
       result, abserr, neval, ier, alist, blist, rlist, elist, iord, nnlog, &
       momcom, chebmo )
  return
end subroutine qawo
