  subroutine qc25c ( f, a, b, c, result, abserr, krul, neval )

    !*****************************************************************************80
    !
    !! QC25C returns integration rules for Cauchy Principal Value integrals.
    !
    !  Discussion:
    !
    !    This routine estimates 
    !      I = integral of F(X) * W(X) over (a,b) 
    !    with error estimate, where 
    !      w(x) = 1/(x-c)
    !
    !  Author:
    !
    !    Robert Piessens, Elise de Doncker-Kapenger, 
    !    Christian Ueberhuber, David Kahaner
    !
    !  Reference:
    !
    !    Robert Piessens, Elise de Doncker-Kapenger, 
    !    Christian Ueberhuber, David Kahaner,
    !    QUADPACK, a Subroutine Package for Automatic Integration,
    !    Springer Verlag, 1983
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
    !    Input, real(8) :: C, the parameter in the weight function.
    !
    !    Output, real(8) :: RESULT, the estimated value of the integral.
    !    RESULT is computed by using a generalized Clenshaw-Curtis method if
    !    C lies within ten percent of the integration interval.  In the 
    !    other case the 15-point Kronrod rule obtained by optimal addition
    !    of abscissae to the 7-point Gauss rule, is applied.
    !
    !    Output, real(8) :: ABSERR, an estimate of || I - RESULT ||.
    !
    !           krul   - integer
    !                    key which is decreased by 1 if the 15-point
    !                    Gauss-Kronrod scheme has been used
    !
    !    Output, integer NEVAL, the number of times the integral was evaluated.
    !
    !  Local parameters:
    !
    !           fval   - value of the function f at the points
    !                    cos(k*pi/24),  k = 0, ..., 24
    !           cheb12 - Chebyshev series expansion coefficients, for the
    !                    function f, of degree 12
    !           cheb24 - Chebyshev series expansion coefficients, for the
    !                    function f, of degree 24
    !           res12  - approximation to the integral corresponding to the
    !                    use of cheb12
    !           res24  - approximation to the integral corresponding to the
    !                    use of cheb24
    !           qwgtc  - external function subprogram defining the weight
    !                    function
    !           hlgth  - half-length of the interval
    !           centr  - mid point of the interval
    !
    implicit none

    real(8) :: a
    real(8) :: abserr
    real(8) :: ak22
    real(8) :: amom0
    real(8) :: amom1
    real(8) :: amom2
    real(8) :: b
    real(8) :: c
    real(8) :: cc
    real(8) :: centr
    real(8) :: cheb12(13)
    real(8) :: cheb24(25)
    real(8), external :: f
    real(8) :: fval(25)
    real(8) :: hlgth
    integer i
    integer isym
    integer k 
    integer kp
    integer krul
    integer neval
    real(8) :: p2
    real(8) :: p3
    real(8) :: p4
    real(8), external :: qwgtc
    real(8) :: resabs
    real(8) :: resasc
    real(8) :: result
    real(8) :: res12
    real(8) :: res24
    real(8) :: u
    real(8), parameter, dimension ( 11 ) :: x = (/ &
         9.914448613738104e-01, 9.659258262890683e-01, &
         9.238795325112868e-01, 8.660254037844386e-01, &
         7.933533402912352e-01, 7.071067811865475e-01, &
         6.087614290087206e-01, 5.000000000000000e-01, &
         3.826834323650898e-01, 2.588190451025208e-01, &
         1.305261922200516e-01 /)
    !
    !  Check the position of C.
    !
    cc = ( 2.0e+00 * c - b - a ) / ( b - a )
    !
    !  Apply the 15-point Gauss-Kronrod scheme.
    !
    if ( abs ( cc ) >= 1.1e+00 ) then
       krul = krul - 1
       call qk15w ( f, qwgtc, c, p2, p3, p4, kp, a, b, result, abserr, &
            resabs, resasc )
       neval = 15
       if ( resasc == abserr ) then
          krul = krul+1
       end if
       return
    end if
    !
    !  Use the generalized Clenshaw-Curtis method.
    !
    hlgth = 5.0e-01 * ( b - a )
    centr = 5.0e-01 * ( b + a )
    neval = 25
    fval(1) = 5.0e-01 * f(hlgth+centr)
    fval(13) = f(centr)
    fval(25) = 5.0e-01 * f(centr-hlgth)

    do i = 2, 12
       u = hlgth * x(i-1)
       isym = 26 - i
       fval(i) = f(u+centr)
       fval(isym) = f(centr-u)
    end do
    !
    !  Compute the Chebyshev series expansion.
    !
    call qcheb ( x, fval, cheb12, cheb24 )
    !
    !  The modified Chebyshev moments are computed by forward
    !  recursion, using AMOM0 and AMOM1 as starting values.
    !
    amom0 = log ( abs ( ( 1.0e+00 - cc ) / ( 1.0e+00 + cc ) ) )
    amom1 = 2.0e+00 + cc * amom0
    res12 = cheb12(1) * amom0 + cheb12(2) * amom1
    res24 = cheb24(1) * amom0 + cheb24(2) * amom1

    do k = 3, 13
       amom2 = 2.0e+00 * cc * amom1 - amom0
       ak22 = ( k - 2 ) * ( k - 2 )
       if ( ( k / 2 ) * 2 == k ) then
          amom2 = amom2 - 4.0e+00 / ( ak22 - 1.0e+00 )
       end if
       res12 = res12 + cheb12(k) * amom2
       res24 = res24 + cheb24(k) * amom2
       amom0 = amom1
       amom1 = amom2
    end do

    do k = 14, 25
       amom2 = 2.0e+00 * cc * amom1 - amom0
       ak22 = ( k - 2 ) * ( k - 2 )
       if ( ( k / 2 ) * 2 == k ) then
          amom2 = amom2 - 4.0e+00 / ( ak22 - 1.0e+00 )
       end if
       res24 = res24 + cheb24(k) * amom2
       amom0 = amom1
       amom1 = amom2
    end do

    result = res24
    abserr = abs ( res24 - res12 )

    return
  end subroutine qc25c
  subroutine qc25o ( f, a, b, omega, integr, nrmom, maxp1, ksave, result, &
       abserr, neval, resabs, resasc, momcom, chebmo )

    !*****************************************************************************80
    !
    !! QC25O returns integration rules for integrands with a COS or SIN factor.
    !
    !  Discussion:
    !
    !    This routine estimates the integral
    !      I = integral of f(x) * w(x) over (a,b)
    !    where
    !      w(x) = cos(omega*x)
    !    or 
    !      w(x) = sin(omega*x),
    !    and estimates
    !      J = integral ( A <= X <= B ) |F(X)| dx.
    !
    !    For small values of OMEGA or small intervals (a,b) the 15-point
    !    Gauss-Kronrod rule is used.  In all other cases a generalized
    !    Clenshaw-Curtis method is used, that is, a truncated Chebyshev 
    !    expansion of the function F is computed on (a,b), so that the 
    !    integrand can be written as a sum of terms of the form W(X)*T(K,X), 
    !    where T(K,X) is the Chebyshev polynomial of degree K.  The Chebyshev
    !    moments are computed with use of a linear recurrence relation.
    !
    !  Author:
    !
    !    Robert Piessens, Elise de Doncker-Kapenger, 
    !    Christian Ueberhuber, David Kahaner
    !
    !  Reference:
    !
    !    Robert Piessens, Elise de Doncker-Kapenger, 
    !    Christian Ueberhuber, David Kahaner,
    !    QUADPACK, a Subroutine Package for Automatic Integration,
    !    Springer Verlag, 1983
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
    !    Input, integer INTEGR, indicates which weight function is to be used
    !    = 1, w(x) = cos(omega*x)
    !    = 2, w(x) = sin(omega*x)
    !
    !    ?, integer NRMOM, the length of interval (a,b) is equal to the length
    !    of the original integration interval divided by
    !    2**nrmom (we suppose that the routine is used in an
    !    adaptive integration process, otherwise set
    !    nrmom = 0).  nrmom must be zero at the first call.
    !
    !           maxp1  - integer
    !                    gives an upper bound on the number of Chebyshev
    !                    moments which can be stored, i.e. for the intervals
    !                    of lengths abs(bb-aa)*2**(-l), l = 0,1,2, ...,
    !                    maxp1-2.
    !
    !           ksave  - integer
    !                    key which is one when the moments for the
    !                    current interval have been computed
    !
    !    Output, real(8) :: RESULT, the estimated value of the integral.
    !
    !           abserr - real
    !                    estimate of the modulus of the absolute
    !                    error, which should equal or exceed abs(i-result)
    !
    !    Output, integer NEVAL, the number of times the integral was evaluated.
    !
    !    Output, real(8) :: RESABS, approximation to the integral J.
    !
    !    Output, real(8) :: RESASC, approximation to the integral of abs(F-I/(B-A)).
    !
    !         on entry and return
    !           momcom - integer
    !                    for each interval length we need to compute
    !                    the Chebyshev moments. momcom counts the number
    !                    of intervals for which these moments have already
    !                    been computed. if nrmom < momcom or ksave = 1,
    !                    the Chebyshev moments for the interval (a,b)
    !                    have already been computed and stored, otherwise
    !                    we compute them and we increase momcom.
    !
    !           chebmo - real
    !                    array of dimension at least (maxp1,25) containing
    !                    the modified Chebyshev moments for the first momcom
    !                    interval lengths
    !
    !  Local parameters:
    !
    !    maxp1 gives an upper bound
    !           on the number of Chebyshev moments which can be
    !           computed, i.e. for the interval (bb-aa), ...,
    !           (bb-aa)/2**(maxp1-2).
    !           should this number be altered, the first dimension of
    !           chebmo needs to be adapted.
    !
    !    x contains the values cos(k*pi/24)
    !           k = 1, ...,11, to be used for the Chebyshev expansion of f
    !
    !           centr  - mid point of the integration interval
    !           hlgth  - half length of the integration interval
    !           fval   - value of the function f at the points
    !                    (b-a)*0.5*cos(k*pi/12) + (b+a)*0.5
    !                    k = 0, ...,24
    !           cheb12 - coefficients of the Chebyshev series expansion
    !                    of degree 12, for the function f, in the
    !                    interval (a,b)
    !           cheb24 - coefficients of the Chebyshev series expansion
    !                    of degree 24, for the function f, in the
    !                    interval (a,b)
    !           resc12 - approximation to the integral of
    !                    cos(0.5*(b-a)*omega*x)*f(0.5*(b-a)*x+0.5*(b+a))
    !                    over (-1,+1), using the Chebyshev series
    !                    expansion of degree 12
    !           resc24 - approximation to the same integral, using the
    !                    Chebyshev series expansion of degree 24
    !           ress12 - the analogue of resc12 for the sine
    !           ress24 - the analogue of resc24 for the sine
    !
    implicit none

    integer maxp1

    real(8) :: a
    real(8) :: abserr
    real(8) :: ac
    real(8) :: an
    real(8) :: an2
    real(8) :: as
    real(8) :: asap
    real(8) :: ass
    real(8) :: b
    real(8) :: centr
    real(8) :: chebmo(maxp1,25)
    real(8) :: cheb12(13)
    real(8) :: cheb24(25)
    real(8) :: conc
    real(8) :: cons
    real(8) :: cospar
    real(8) :: d(28)
    real(8) :: d1(28)
    real(8) :: d2(28)
    real(8) :: d3(28)
    real(8) :: estc
    real(8) :: ests
    real(8), external :: f
    real(8) :: fval(25)
    real(8) :: hlgth
    integer i
    integer integr
    integer isym
    integer j
    integer k
    integer ksave
    integer m
    integer momcom
    integer neval
    integer, parameter :: nmac = 28
    integer noeq1
    integer noequ
    integer nrmom
    real(8) :: omega
    real(8) :: parint
    real(8) :: par2
    real(8) :: par22
    real(8) :: p2
    real(8) :: p3
    real(8) :: p4
    real(8), external :: qwgto
    real(8) :: resabs
    real(8) :: resasc
    real(8) :: resc12
    real(8) :: resc24
    real(8) :: ress12
    real(8) :: ress24
    real(8) :: result
    real(8) :: sinpar
    real(8) :: v(28)
    real(8), dimension ( 11 ) :: x = (/ &
         9.914448613738104e-01,     9.659258262890683e-01, &
         9.238795325112868e-01,     8.660254037844386e-01, &
         7.933533402912352e-01,     7.071067811865475e-01, &
         6.087614290087206e-01,     5.000000000000000e-01, &
         3.826834323650898e-01,     2.588190451025208e-01, &
         1.305261922200516e-01 /)

    centr = 5.0e-01 * ( b + a )
    hlgth = 5.0e-01 * ( b - a )
    parint = omega * hlgth
    !
    !  Compute the integral using the 15-point Gauss-Kronrod
    !  formula if the value of the parameter in the integrand
    !  is small or if the length of the integration interval
    !  is less than (bb-aa)/2**(maxp1-2), where (aa,bb) is the
    !  original integration interval.
    !
    if ( abs ( parint ) <= 2.0e+00 ) then

       call qk15w ( f, qwgto, omega, p2, p3, p4, integr, a, b, result, &
            abserr, resabs, resasc )

       neval = 15
       return

    end if
    !
    !  Compute the integral using the generalized clenshaw-curtis method.
    !
    conc = hlgth * cos(centr*omega)
    cons = hlgth * sin(centr*omega)
    resasc = huge ( resasc )
    neval = 25
    !
    !  Check whether the Chebyshev moments for this interval
    !  have already been computed.
    !
    if ( nrmom < momcom .or. ksave == 1 ) then
       go to 140
    end if
    !
    !  Compute a new set of Chebyshev moments.
    !
    m = momcom + 1
    par2 = parint * parint
    par22 = par2 + 2.0e+00
    sinpar = sin(parint)
    cospar = cos(parint)
    !
    !  Compute the Chebyshev moments with respect to cosine.
    !
    v(1) = 2.0e+00 * sinpar / parint
    v(2) = (8.0e+00*cospar+(par2+par2-8.0e+00)*sinpar/ parint)/par2
    v(3) = (3.2e+01*(par2-1.2e+01)*cospar+(2.0e+00* &
         ((par2-8.0e+01)*par2+1.92e+02)*sinpar)/ &
         parint)/(par2*par2)
    ac = 8.0e+00*cospar
    as = 2.4e+01*parint*sinpar

    if ( abs ( parint ) > 2.4e+01 ) then
       go to 70
    end if
    !
    !  Compute the Chebyshev moments as the solutions of a boundary value 
    !  problem with one initial value (v(3)) and one end value computed
    !  using an asymptotic formula.
    !
    noequ = nmac-3
    noeq1 = noequ-1
    an = 6.0e+00

    do k = 1, noeq1
       an2 = an*an
       d(k) = -2.0e+00*(an2-4.0e+00) * (par22-an2-an2)
       d2(k) = (an-1.0e+00)*(an-2.0e+00) * par2
       d1(k) = (an+3.0e+00)*(an+4.0e+00) * par2
       v(k+3) = as-(an2-4.0e+00) * ac
       an = an+2.0e+00
    end do

    an2 = an*an
    d(noequ) = -2.0e+00*(an2-4.0e+00) * (par22-an2-an2)
    v(noequ+3) = as - ( an2 - 4.0e+00 ) * ac
    v(4) = v(4) - 5.6e+01 * par2 * v(3)
    ass = parint * sinpar
    asap = (((((2.10e+02*par2-1.0e+00)*cospar-(1.05e+02*par2 &
         -6.3e+01)*ass)/an2-(1.0e+00-1.5e+01*par2)*cospar &
         +1.5e+01*ass)/an2-cospar+3.0e+00*ass)/an2-cospar)/an2
    v(noequ+3) = v(noequ+3)-2.0e+00*asap*par2*(an-1.0e+00)* &
         (an-2.0e+00)
    !
    !  Solve the tridiagonal system by means of Gaussian
    !  elimination with partial pivoting.
    !
    d3(1:noequ) = 0.0e+00

    d2(noequ) = 0.0e+00

    do i = 1, noeq1

       if ( abs(d1(i)) > abs(d(i)) ) then
          an = d1(i)
          d1(i) = d(i)
          d(i) = an
          an = d2(i)
          d2(i) = d(i+1)
          d(i+1) = an
          d3(i) = d2(i+1)
          d2(i+1) = 0.0e+00
          an = v(i+4)
          v(i+4) = v(i+3)
          v(i+3) = an
       end if

       d(i+1) = d(i+1)-d2(i)*d1(i)/d(i)
       d2(i+1) = d2(i+1)-d3(i)*d1(i)/d(i)
       v(i+4) = v(i+4)-v(i+3)*d1(i)/d(i)

    end do

    v(noequ+3) = v(noequ+3) / d(noequ)
    v(noequ+2) = (v(noequ+2)-d2(noeq1)*v(noequ+3))/d(noeq1)

    do i = 2, noeq1
       k = noequ-i
       v(k+3) = (v(k+3)-d3(k)*v(k+5)-d2(k)*v(k+4))/d(k)
    end do

    go to 90
    !
    !  Compute the Chebyshev moments by means of forward recursion
    !
70  continue

    an = 4.0e+00

    do i = 4, 13
       an2 = an*an
       v(i) = ((an2-4.0e+00)*(2.0e+00*(par22-an2-an2)*v(i-1)-ac) &
            +as-par2*(an+1.0e+00)*(an+2.0e+00)*v(i-2))/ &
            (par2*(an-1.0e+00)*(an-2.0e+00))
       an = an+2.0e+00
    end do

90  continue

    do j = 1, 13
       chebmo(m,2*j-1) = v(j)
    end do
    !
    !  Compute the Chebyshev moments with respect to sine.
    !
    v(1) = 2.0e+00*(sinpar-parint*cospar)/par2
    v(2) = (1.8e+01-4.8e+01/par2)*sinpar/par2 &
         +(-2.0e+00+4.8e+01/par2)*cospar/parint
    ac = -2.4e+01*parint*cospar
    as = -8.0e+00*sinpar
    chebmo(m,2) = v(1)
    chebmo(m,4) = v(2)

    if ( abs(parint) <= 2.4e+01 ) then

       do k = 3, 12
          an = k
          chebmo(m,2*k) = -sinpar/(an*(2.0e+00*an-2.0e+00)) &
               -2.5e-01*parint*(v(k+1)/an-v(k)/(an-1.0e+00))
       end do
       !
       !  Compute the Chebyshev moments by means of forward recursion.
       !
    else

       an = 3.0e+00

       do i = 3, 12
          an2 = an*an
          v(i) = ((an2-4.0e+00)*(2.0e+00*(par22-an2-an2)*v(i-1)+as) &
               +ac-par2*(an+1.0e+00)*(an+2.0e+00)*v(i-2)) &
               /(par2*(an-1.0e+00)*(an-2.0e+00))
          an = an+2.0e+00
          chebmo(m,2*i) = v(i)
       end do

    end if

140 continue

    if ( nrmom < momcom ) then
       m = nrmom + 1
    end if

    if ( momcom < maxp1 - 1 .and. nrmom >= momcom ) then
       momcom = momcom + 1
    end if
    !
    !  Compute the coefficients of the Chebyshev expansions
    !  of degrees 12 and 24 of the function F.
    !
    fval(1) = 5.0e-01 * f(centr+hlgth)
    fval(13) = f(centr)
    fval(25) = 5.0e-01 * f(centr-hlgth)

    do i = 2, 12
       isym = 26-i
       fval(i) = f(hlgth*x(i-1)+centr)
       fval(isym) = f(centr-hlgth*x(i-1))
    end do

    call qcheb ( x, fval, cheb12, cheb24 )
    !
    !  Compute the integral and error estimates.
    !
    resc12 = cheb12(13) * chebmo(m,13)
    ress12 = 0.0e+00
    estc = abs ( cheb24(25)*chebmo(m,25))+abs((cheb12(13)- &
         cheb24(13))*chebmo(m,13) )
    ests = 0.0e+00
    k = 11

    do j = 1, 6
       resc12 = resc12+cheb12(k)*chebmo(m,k)
       ress12 = ress12+cheb12(k+1)*chebmo(m,k+1)
       estc = estc+abs((cheb12(k)-cheb24(k))*chebmo(m,k))
       ests = ests+abs((cheb12(k+1)-cheb24(k+1))*chebmo(m,k+1))
       k = k-2
    end do

    resc24 = cheb24(25)*chebmo(m,25)
    ress24 = 0.0e+00
    resabs = abs(cheb24(25))
    k = 23

    do j = 1, 12

       resc24 = resc24+cheb24(k)*chebmo(m,k)
       ress24 = ress24+cheb24(k+1)*chebmo(m,k+1)
       resabs = resabs+abs(cheb24(k))+abs(cheb24(k+1))

       if ( j <= 5 ) then
          estc = estc+abs(cheb24(k)*chebmo(m,k))
          ests = ests+abs(cheb24(k+1)*chebmo(m,k+1))
       end if

       k = k-2

    end do

    resabs = resabs * abs ( hlgth )

    if ( integr == 1 ) then
       result = conc * resc24-cons*ress24
       abserr = abs ( conc * estc ) + abs ( cons * ests )
    else
       result = conc*ress24+cons*resc24
       abserr = abs(conc*ests)+abs(cons*estc)
    end if

    return
  end subroutine qc25o
  subroutine qc25s ( f, a, b, bl, br, alfa, beta, ri, rj, rg, rh, result, &
       abserr, resasc, integr, neval )

    !*****************************************************************************80
    !
    !! QC25S returns rules for algebraico-logarithmic end point singularities.
    !
    !  Discussion:
    !
    !    This routine computes 
    !      i = integral of F(X) * W(X) over (bl,br), 
    !    with error estimate, where the weight function W(X) has a singular
    !    behavior of algebraico-logarithmic type at the points
    !    a and/or b. 
    !
    !    The interval (bl,br) is a subinterval of (a,b).
    !
    !  Author:
    !
    !    Robert Piessens, Elise de Doncker-Kapenger, 
    !    Christian Ueberhuber, David Kahaner
    !
    !  Reference:
    !
    !    Robert Piessens, Elise de Doncker-Kapenger, 
    !    Christian Ueberhuber, David Kahaner,
    !    QUADPACK, a Subroutine Package for Automatic Integration,
    !    Springer Verlag, 1983
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
    !    Input, real(8) :: BL, BR, the lower and upper limits of integration.
    !    A <= BL < BR <= B.
    !
    !    Input, real(8) :: ALFA, BETA, parameters in the weight function.
    !
    !    Input, real(8) :: RI(25), RJ(25), RG(25), RH(25), modified Chebyshev moments 
    !    for the application of the generalized Clenshaw-Curtis method,
    !    computed in QMOMO.
    !
    !    Output, real(8) :: RESULT, the estimated value of the integral, computed by 
    !    using a generalized clenshaw-curtis method if b1 = a or br = b.
    !    In all other cases the 15-point Kronrod rule is applied, obtained by
    !    optimal addition of abscissae to the 7-point Gauss rule.
    !
    !    Output, real(8) :: ABSERR, an estimate of || I - RESULT ||.
    !
    !    Output, real(8) :: RESASC, approximation to the integral of abs(F*W-I/(B-A)).
    !
    !    Input, integer INTEGR,  determines the weight function
    !    1, w(x) = (x-a)**alfa*(b-x)**beta
    !    2, w(x) = (x-a)**alfa*(b-x)**beta*log(x-a)
    !    3, w(x) = (x-a)**alfa*(b-x)**beta*log(b-x)
    !    4, w(x) = (x-a)**alfa*(b-x)**beta*log(x-a)*log(b-x)
    !
    !    Output, integer NEVAL, the number of times the integral was evaluated.
    !
    !  Local Parameters:
    !
    !           fval   - value of the function f at the points
    !                    (br-bl)*0.5*cos(k*pi/24)+(br+bl)*0.5
    !                    k = 0, ..., 24
    !           cheb12 - coefficients of the Chebyshev series expansion
    !                    of degree 12, for the function f, in the interval
    !                    (bl,br)
    !           cheb24 - coefficients of the Chebyshev series expansion
    !                    of degree 24, for the function f, in the interval
    !                    (bl,br)
    !           res12  - approximation to the integral obtained from cheb12
    !           res24  - approximation to the integral obtained from cheb24
    !           qwgts  - external function subprogram defining the four
    !                    possible weight functions
    !           hlgth  - half-length of the interval (bl,br)
    !           centr  - mid point of the interval (bl,br)
    !
    !           the vector x contains the values cos(k*pi/24)
    !           k = 1, ..., 11, to be used for the computation of the
    !           Chebyshev series expansion of f.
    !
    implicit none

    real(8) :: a
    real(8) :: abserr
    real(8) :: alfa
    real(8) :: b
    real(8) :: beta
    real(8) :: bl
    real(8) :: br
    real(8) :: centr
    real(8) :: cheb12(13)
    real(8) :: cheb24(25)
    real(8) :: dc
    real(8), external :: f
    real(8) :: factor
    real(8) :: fix
    real(8) :: fval(25)
    real(8) :: hlgth
    integer i
    integer integr
    integer isym
    integer neval
    real(8), external :: qwgts
    real(8) :: resabs
    real(8) :: resasc
    real(8) :: result
    real(8) :: res12
    real(8) :: res24
    real(8) :: rg(25)
    real(8) :: rh(25)
    real(8) :: ri(25)
    real(8) :: rj(25)
    real(8) :: u
    real(8), dimension ( 11 ) :: x = (/ &
         9.914448613738104e-01,     9.659258262890683e-01, &
         9.238795325112868e-01,     8.660254037844386e-01, &
         7.933533402912352e-01,     7.071067811865475e-01, &
         6.087614290087206e-01,     5.000000000000000e-01, &
         3.826834323650898e-01,     2.588190451025208e-01, &
         1.305261922200516e-01 /)

    neval = 25

    if ( bl == a .and. (alfa /= 0.0e+00 .or. integr == 2 .or. integr == 4)) then
       go to 10
    end if

    if ( br == b .and. (beta /= 0.0e+00 .or. integr == 3 .or. integr == 4)) &
         go to 140
    !
    !  If a > bl and b < br, apply the 15-point Gauss-Kronrod scheme.
    !
    call qk15w ( f, qwgts, a, b, alfa, beta, integr, bl, br, result, abserr, &
         resabs, resasc )

    neval = 15
    return
    !
    !  This part of the program is executed only if a = bl.
    !
    !  Compute the Chebyshev series expansion of the function
    !  f1 = (0.5*(b+b-br-a)-0.5*(br-a)*x)**beta*f(0.5*(br-a)*x+0.5*(br+a))
    !
10  continue

    hlgth = 5.0e-01*(br-bl)
    centr = 5.0e-01*(br+bl)
    fix = b-centr
    fval(1) = 5.0e-01*f(hlgth+centr)*(fix-hlgth)**beta
    fval(13) = f(centr)*(fix**beta)
    fval(25) = 5.0e-01*f(centr-hlgth)*(fix+hlgth)**beta

    do i = 2, 12
       u = hlgth*x(i-1)
       isym = 26-i
       fval(i) = f(u+centr)*(fix-u)**beta
       fval(isym) = f(centr-u)*(fix+u)**beta
    end do

    factor = hlgth**(alfa+1.0e+00)
    result = 0.0e+00
    abserr = 0.0e+00
    res12 = 0.0e+00
    res24 = 0.0e+00

    if ( integr > 2 ) go to 70

    call qcheb ( x, fval, cheb12, cheb24 )
    !
    !  integr = 1  (or 2)
    !
    do i = 1, 13
       res12 = res12+cheb12(i)*ri(i)
       res24 = res24+cheb24(i)*ri(i)
    end do

    do i = 14, 25
       res24 = res24 + cheb24(i) * ri(i)
    end do

    if ( integr == 1 ) go to 130
    !
    !  integr = 2
    !
    dc = log ( br - bl )
    result = res24 * dc
    abserr = abs((res24-res12)*dc)
    res12 = 0.0e+00
    res24 = 0.0e+00

    do i = 1, 13
       res12 = res12+cheb12(i)*rg(i)
       res24 = res24+cheb24(i)*rg(i)
    end do

    do i = 14, 25
       res24 = res24+cheb24(i)*rg(i)
    end do

    go to 130
    !
    !  Compute the Chebyshev series expansion of the function
    !  F4 = f1*log(0.5*(b+b-br-a)-0.5*(br-a)*x)
    !
70  continue

    fval(1) = fval(1) * log ( fix - hlgth )
    fval(13) = fval(13) * log ( fix )
    fval(25) = fval(25) * log ( fix + hlgth )

    do i = 2, 12
       u = hlgth*x(i-1)
       isym = 26-i
       fval(i) = fval(i) * log ( fix - u )
       fval(isym) = fval(isym) * log ( fix + u )
    end do

    call qcheb ( x, fval, cheb12, cheb24 )
    !
    !  integr = 3  (or 4)
    !
    do i = 1, 13
       res12 = res12+cheb12(i)*ri(i)
       res24 = res24+cheb24(i)*ri(i)
    end do

    do i = 14, 25
       res24 = res24+cheb24(i)*ri(i)
    end do

    if ( integr == 3 ) then
       go to 130
    end if
    !
    !  integr = 4
    !
    dc = log ( br - bl )
    result = res24*dc
    abserr = abs((res24-res12)*dc)
    res12 = 0.0e+00
    res24 = 0.0e+00

    do i = 1, 13
       res12 = res12+cheb12(i)*rg(i)
       res24 = res24+cheb24(i)*rg(i)
    end do

    do i = 14, 25
       res24 = res24+cheb24(i)*rg(i)
    end do

130 continue

    result = (result+res24)*factor
    abserr = (abserr+abs(res24-res12))*factor
    go to 270
    !
    !  This part of the program is executed only if b = br.
    !
    !  Compute the Chebyshev series expansion of the function
    !  f2 = (0.5*(b+bl-a-a)+0.5*(b-bl)*x)**alfa*f(0.5*(b-bl)*x+0.5*(b+bl))
    !
140 continue

    hlgth = 5.0e-01*(br-bl)
    centr = 5.0e-01*(br+bl)
    fix = centr-a
    fval(1) = 5.0e-01*f(hlgth+centr)*(fix+hlgth)**alfa
    fval(13) = f(centr)*(fix**alfa)
    fval(25) = 5.0e-01*f(centr-hlgth)*(fix-hlgth)**alfa

    do i = 2, 12
       u = hlgth*x(i-1)
       isym = 26-i
       fval(i) = f(u+centr)*(fix+u)**alfa
       fval(isym) = f(centr-u)*(fix-u)**alfa
    end do

    factor = hlgth**(beta+1.0e+00)
    result = 0.0e+00
    abserr = 0.0e+00
    res12 = 0.0e+00
    res24 = 0.0e+00

    if ( integr == 2 .or. integr == 4 ) then
       go to 200
    end if
    !
    !  integr = 1  (or 3)
    !
    call qcheb ( x, fval, cheb12, cheb24 )

    do i = 1, 13
       res12 = res12+cheb12(i)*rj(i)
       res24 = res24+cheb24(i)*rj(i)
    end do

    do i = 14, 25
       res24 = res24+cheb24(i)*rj(i)
    end do

    if ( integr == 1 ) go to 260
    !
    !  integr = 3
    !
    dc = log ( br - bl )
    result = res24*dc
    abserr = abs((res24-res12)*dc)
    res12 = 0.0e+00
    res24 = 0.0e+00

    do i = 1, 13
       res12 = res12+cheb12(i)*rh(i)
       res24 = res24+cheb24(i)*rh(i)
    end do

    do i = 14, 25
       res24 = res24+cheb24(i)*rh(i)
    end do

    go to 260
    !
    !  Compute the Chebyshev series expansion of the function
    !  f3 = f2*log(0.5*(b-bl)*x+0.5*(b+bl-a-a))
    !
200 continue

    fval(1) = fval(1) * log ( hlgth + fix )
    fval(13) = fval(13) * log ( fix )
    fval(25) = fval(25) * log ( fix - hlgth )

    do i = 2, 12
       u = hlgth*x(i-1)
       isym = 26-i
       fval(i) = fval(i) * log(u+fix)
       fval(isym) = fval(isym) * log(fix-u)
    end do

    call qcheb ( x, fval, cheb12, cheb24 )
    !
    !  integr = 2  (or 4)
    !
    do i = 1, 13
       res12 = res12+cheb12(i)*rj(i)
       res24 = res24+cheb24(i)*rj(i)
    end do

    do i = 14, 25
       res24 = res24+cheb24(i)*rj(i)
    end do

    if ( integr == 2 ) go to 260

    dc = log(br-bl)
    result = res24*dc
    abserr = abs((res24-res12)*dc)
    res12 = 0.0e+00
    res24 = 0.0e+00
    !
    !  integr = 4
    !
    do i = 1, 13
       res12 = res12+cheb12(i)*rh(i)
       res24 = res24+cheb24(i)*rh(i)
    end do

    do i = 14, 25
       res24 = res24+cheb24(i)*rh(i)
    end do

260 continue

    result = (result+res24)*factor
    abserr = (abserr+abs(res24-res12))*factor

270 continue

    return
  end subroutine qc25s
  subroutine qcheb ( x, fval, cheb12, cheb24 )

    !*****************************************************************************80
    !
    !! QCHEB computes the Chebyshev series expansion.
    !
    !  Discussion:
    !
    !    This routine computes the Chebyshev series expansion
    !    of degrees 12 and 24 of a function using a fast Fourier transform method
    !
    !      f(x) = sum(k=1, ...,13) (cheb12(k)*t(k-1,x)),
    !      f(x) = sum(k=1, ...,25) (cheb24(k)*t(k-1,x)),
    !
    !    where T(K,X) is the Chebyshev polynomial of degree K.
    !
    !  Author:
    !
    !    Robert Piessens, Elise de Doncker-Kapenger, 
    !    Christian Ueberhuber, David Kahaner
    !
    !  Reference:
    !
    !    Robert Piessens, Elise de Doncker-Kapenger, 
    !    Christian Ueberhuber, David Kahaner,
    !    QUADPACK, a Subroutine Package for Automatic Integration,
    !    Springer Verlag, 1983
    !
    !  Parameters:
    !
    !    Input, real(8) :: X(11), contains the values of COS(K*PI/24), for K = 1 to 11.
    !
    !    Input/output, real(8) :: FVAL(25), the function values at the points
    !    (b+a+(b-a)*cos(k*pi/24))/2, k = 0, ...,24, where (a,b) is the 
    !    approximation interval.  FVAL(1) and FVAL(25) are divided by two
    !    These values are destroyed at output.
    !
    !    Output, real(8) :: CHEB12(13), the Chebyshev coefficients for degree 12.
    !
    !    Output, real(8) :: CHEB24(25), the Chebyshev coefficients for degree 24.
    !
    implicit none

    real(8) :: alam
    real(8) :: alam1
    real(8) :: alam2
    real(8) :: cheb12(13)
    real(8) :: cheb24(25)
    real(8) :: fval(25)
    integer i
    integer j
    real(8) :: part1
    real(8) :: part2
    real(8) :: part3
    real(8) :: v(12)
    real(8) :: x(11)

    do i = 1, 12
       j = 26-i
       v(i) = fval(i)-fval(j)
       fval(i) = fval(i)+fval(j)
    end do

    alam1 = v(1)-v(9)
    alam2 = x(6)*(v(3)-v(7)-v(11))
    cheb12(4) = alam1+alam2
    cheb12(10) = alam1-alam2
    alam1 = v(2)-v(8)-v(10)
    alam2 = v(4)-v(6)-v(12)
    alam = x(3)*alam1+x(9)*alam2
    cheb24(4) = cheb12(4)+alam
    cheb24(22) = cheb12(4)-alam
    alam = x(9)*alam1-x(3)*alam2
    cheb24(10) = cheb12(10)+alam
    cheb24(16) = cheb12(10)-alam
    part1 = x(4)*v(5)
    part2 = x(8)*v(9)
    part3 = x(6)*v(7)
    alam1 = v(1)+part1+part2
    alam2 = x(2)*v(3)+part3+x(10)*v(11)
    cheb12(2) = alam1+alam2
    cheb12(12) = alam1-alam2
    alam = x(1)*v(2)+x(3)*v(4)+x(5)*v(6)+x(7)*v(8) &
         +x(9)*v(10)+x(11)*v(12)
    cheb24(2) = cheb12(2)+alam
    cheb24(24) = cheb12(2)-alam
    alam = x(11)*v(2)-x(9)*v(4)+x(7)*v(6)-x(5)*v(8) &
         +x(3)*v(10)-x(1)*v(12)
    cheb24(12) = cheb12(12)+alam
    cheb24(14) = cheb12(12)-alam
    alam1 = v(1)-part1+part2
    alam2 = x(10)*v(3)-part3+x(2)*v(11)
    cheb12(6) = alam1+alam2
    cheb12(8) = alam1-alam2
    alam = x(5)*v(2)-x(9)*v(4)-x(1)*v(6) &
         -x(11)*v(8)+x(3)*v(10)+x(7)*v(12)
    cheb24(6) = cheb12(6)+alam
    cheb24(20) = cheb12(6)-alam
    alam = x(7)*v(2)-x(3)*v(4)-x(11)*v(6)+x(1)*v(8) &
         -x(9)*v(10)-x(5)*v(12)
    cheb24(8) = cheb12(8)+alam
    cheb24(18) = cheb12(8)-alam

    do i = 1, 6
       j = 14-i
       v(i) = fval(i)-fval(j)
       fval(i) = fval(i)+fval(j)
    end do

    alam1 = v(1)+x(8)*v(5)
    alam2 = x(4)*v(3)
    cheb12(3) = alam1+alam2
    cheb12(11) = alam1-alam2
    cheb12(7) = v(1)-v(5)
    alam = x(2)*v(2)+x(6)*v(4)+x(10)*v(6)
    cheb24(3) = cheb12(3)+alam
    cheb24(23) = cheb12(3)-alam
    alam = x(6)*(v(2)-v(4)-v(6))
    cheb24(7) = cheb12(7)+alam
    cheb24(19) = cheb12(7)-alam
    alam = x(10)*v(2)-x(6)*v(4)+x(2)*v(6)
    cheb24(11) = cheb12(11)+alam
    cheb24(15) = cheb12(11)-alam

    do i = 1, 3
       j = 8-i
       v(i) = fval(i)-fval(j)
       fval(i) = fval(i)+fval(j)
    end do

    cheb12(5) = v(1)+x(8)*v(3)
    cheb12(9) = fval(1)-x(8)*fval(3)
    alam = x(4)*v(2)
    cheb24(5) = cheb12(5)+alam
    cheb24(21) = cheb12(5)-alam
    alam = x(8)*fval(2)-fval(4)
    cheb24(9) = cheb12(9)+alam
    cheb24(17) = cheb12(9)-alam
    cheb12(1) = fval(1)+fval(3)
    alam = fval(2)+fval(4)
    cheb24(1) = cheb12(1)+alam
    cheb24(25) = cheb12(1)-alam
    cheb12(13) = v(1)-v(3)
    cheb24(13) = cheb12(13)
    alam = 1.0e+00/6.0e+00

    do i = 2, 12
       cheb12(i) = cheb12(i)*alam
    end do

    alam = 5.0e-01*alam
    cheb12(1) = cheb12(1)*alam
    cheb12(13) = cheb12(13)*alam

    do i = 2, 24
       cheb24(i) = cheb24(i)*alam
    end do

    cheb24(1) = 0.5E+00 * alam*cheb24(1)
    cheb24(25) = 0.5E+00 * alam*cheb24(25)

    return
  end subroutine qcheb
  subroutine qextr ( n, epstab, result, abserr, res3la, nres )

    !*****************************************************************************80
    !
    !! QEXTR carries out the Epsilon extrapolation algorithm.
    !
    !  Discussion:
    !
    !    The routine determines the limit of a given sequence of approximations, 
    !    by means of the epsilon algorithm of P. Wynn.  An estimate of the 
    !    absolute error is also given.  The condensed epsilon table is computed.
    !    Only those elements needed for the computation of the next diagonal
    !    are preserved.
    !
    !  Author:
    !
    !    Robert Piessens, Elise de Doncker-Kapenger, 
    !    Christian Ueberhuber, David Kahaner
    !
    !  Reference:
    !
    !    Robert Piessens, Elise de Doncker-Kapenger, 
    !    Christian Ueberhuber, David Kahaner,
    !    QUADPACK, a Subroutine Package for Automatic Integration,
    !    Springer Verlag, 1983
    !
    !  Parameters:
    !
    !    Input, integer N, indicates the entry of EPSTAB which contains
    !    the new element in the first column of the epsilon table.
    !
    !    Input/output, real(8) :: EPSTAB(52), the two lower diagonals of the triangular
    !    epsilon table.  The elements are numbered starting at the right-hand 
    !    corner of the triangle.
    !
    !    Output, real(8) :: RESULT, the estimated value of the integral.
    !
    !    Output, real(8) :: ABSERR, estimate of the absolute error computed from
    !    RESULT and the 3 previous results.
    !
    !    ?, real(8) :: RES3LA(3), the last 3 results.
    !
    !    Input/output, integer NRES, the number of calls to the routine.  This
    !    should be zero on the first call, and is automatically updated
    !    before return.
    !
    !  Local Parameters:
    !
    !           e0     - the 4 elements on which the
    !           e1       computation of a new element in
    !           e2       the epsilon table is based
    !           e3                 e0
    !                        e3    e1    new
    !                              e2
    !           newelm - number of elements to be computed in the new
    !                    diagonal
    !           error  - error = abs(e1-e0)+abs(e2-e1)+abs(new-e2)
    !           result - the element in the new diagonal with least value
    !                    of error
    !           limexp is the maximum number of elements the epsilon table
    !           can contain. if this number is reached, the upper diagonal
    !           of the epsilon table is deleted.
    !
    implicit none

    real(8) :: abserr
    real(8) :: delta1
    real(8) :: delta2
    real(8) :: delta3
    real(8) :: epsinf
    real(8) :: epstab(52)
    real(8) :: error
    real(8) :: err1
    real(8) :: err2
    real(8) :: err3
    real(8) :: e0
    real(8) :: e1
    real(8) :: e1abs
    real(8) :: e2
    real(8) :: e3
    integer i
    integer ib
    integer ib2
    integer ie
    integer indx
    integer k1
    integer k2
    integer k3
    integer limexp
    integer n
    integer newelm
    integer nres
    integer num
    real(8) :: res
    real(8) :: result
    real(8) :: res3la(3)
    real(8) :: ss
    real(8) :: tol1
    real(8) :: tol2
    real(8) :: tol3

    nres = nres+1
    abserr = huge ( abserr )
    result = epstab(n)

    if ( n < 3 ) then
       abserr = max ( abserr,0.5e+00* epsilon ( result ) *abs(result))
       return
    end if

    limexp = 50
    epstab(n+2) = epstab(n)
    newelm = (n-1)/2
    epstab(n) = huge ( epstab(n) )
    num = n
    k1 = n

    do i = 1, newelm

       k2 = k1-1
       k3 = k1-2
       res = epstab(k1+2)
       e0 = epstab(k3)
       e1 = epstab(k2)
       e2 = res
       e1abs = abs(e1)
       delta2 = e2-e1
       err2 = abs(delta2)
       tol2 = max ( abs(e2),e1abs)* epsilon ( e2 )
       delta3 = e1-e0
       err3 = abs(delta3)
       tol3 = max ( e1abs,abs(e0))* epsilon ( e0 )
       !
       !  If e0, e1 and e2 are equal to within machine accuracy, convergence 
       !  is assumed.
       !
       if ( err2 <= tol2 .and. err3 <= tol3 ) then
          result = res
          abserr = err2+err3
          abserr = max ( abserr,0.5e+00* epsilon ( result ) *abs(result))
          return
       end if

       e3 = epstab(k1)
       epstab(k1) = e1
       delta1 = e1-e3
       err1 = abs(delta1)
       tol1 = max ( e1abs,abs(e3))* epsilon ( e3 )
       !
       !  If two elements are very close to each other, omit a part
       !  of the table by adjusting the value of N.
       !
       if ( err1 <= tol1 .or. err2 <= tol2 .or. err3 <= tol3 ) go to 20

       ss = 1.0e+00/delta1+1.0e+00/delta2-1.0e+00/delta3
       epsinf = abs ( ss*e1 )
       !
       !  Test to detect irregular behavior in the table, and
       !  eventually omit a part of the table adjusting the value of N.
       !
       if ( epsinf > 1.0e-04 ) go to 30

20     continue

       n = i+i-1
       exit
       !
       !  Compute a new element and eventually adjust the value of RESULT.
       !
30     continue

       res = e1+1.0e+00/ss
       epstab(k1) = res
       k1 = k1-2
       error = err2+abs(res-e2)+err3

       if ( error <= abserr ) then
          abserr = error
          result = res
       end if

    end do
    !
    !  Shift the table.
    !
    if ( n == limexp ) then
       n = 2*(limexp/2)-1
    end if

    if ( (num/2)*2 == num ) then
       ib = 2
    else
       ib = 1
    end if

    ie = newelm+1

    do i = 1, ie
       ib2 = ib+2
       epstab(ib) = epstab(ib2)
       ib = ib2
    end do

    if ( num /= n ) then

       indx = num-n+1

       do i = 1, n
          epstab(i)= epstab(indx)
          indx = indx+1
       end do

    end if

    if ( nres < 4 ) then
       res3la(nres) = result
       abserr = huge ( abserr )
    else
       abserr = abs(result-res3la(3))+abs(result-res3la(2)) &
            +abs(result-res3la(1))
       res3la(1) = res3la(2)
       res3la(2) = res3la(3)
       res3la(3) = result
    end if

    abserr = max ( abserr,0.5e+00* epsilon ( result ) *abs(result))

    return
  end subroutine qextr
  subroutine qfour ( f, a, b, omega, integr, epsabs, epsrel, limit, icall, &
       maxp1, result, abserr, neval, ier, alist, blist, rlist, elist, iord, &
       nnlog, momcom, chebmo )

    !*****************************************************************************80
    !
    !! QFOUR estimates the integrals of oscillatory functions.
    !
    !  Discussion:
    !
    !    This routine calculates an approximation RESULT to a definite integral
    !      I = integral of F(X) * COS(OMEGA*X) 
    !    or
    !      I = integral of F(X) * SIN(OMEGA*X) 
    !    over (A,B), hopefully satisfying:
    !      | I - RESULT | <= max ( epsabs, epsrel * |I| ) ).
    !
    !    QFOUR is called by QAWO and QAWF.  It can also be called directly in 
    !    a user-written program.  In the latter case it is possible for the 
    !    user to determine the first dimension of array CHEBMO(MAXP1,25).
    !    See also parameter description of MAXP1.  Additionally see
    !    parameter description of ICALL for eventually re-using
    !    Chebyshev moments computed during former call on subinterval
    !    of equal length abs(B-A).
    !
    !  Author:
    !
    !    Robert Piessens, Elise de Doncker-Kapenger, 
    !    Christian Ueberhuber, David Kahaner
    !
    !  Reference:
    !
    !    Robert Piessens, Elise de Doncker-Kapenger, 
    !    Christian Ueberhuber, David Kahaner,
    !    QUADPACK, a Subroutine Package for Automatic Integration,
    !    Springer Verlag, 1983
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
    !    Input, real(8) :: OMEGA, the multiplier of X in the weight function.
    !
    !    Input, integer INTEGR, indicates the weight functions to be used.
    !    = 1, w(x) = cos(omega*x)
    !    = 2, w(x) = sin(omega*x)
    !
    !    Input, real(8) :: EPSABS, EPSREL, the absolute and relative accuracy requested.
    !
    !    Input, integer LIMIT, the maximum number of subintervals of [A,B]
    !    that can be generated.
    !
    !    icall  - integer
    !                     if qfour is to be used only once, ICALL must
    !                     be set to 1.  assume that during this call, the
    !                     Chebyshev moments (for clenshaw-curtis integration
    !                     of degree 24) have been computed for intervals of
    !                     lenghts (abs(b-a))*2**(-l), l=0,1,2,...momcom-1.
    !                     the Chebyshev moments already computed can be
    !                     re-used in subsequent calls, if qfour must be
    !                     called twice or more times on intervals of the
    !                     same length abs(b-a). from the second call on, one
    !                     has to put then ICALL > 1.
    !                     if ICALL < 1, the routine will end with ier = 6.
    !
    !            maxp1  - integer
    !                     gives an upper bound on the number of
    !                     Chebyshev moments which can be stored, i.e.
    !                     for the intervals of lenghts abs(b-a)*2**(-l),
    !                     l=0,1, ..., maxp1-2, maxp1 >= 1.
    !                     if maxp1 < 1, the routine will end with ier = 6.
    !                     increasing (decreasing) the value of maxp1
    !                     decreases (increases) the computational time but
    !                     increases (decreases) the required memory space.
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
    !                             has been achieved. one can allow more
    !                             subdivisions by increasing the value of
    !                             limit (and taking according dimension
    !                             adjustments into account). however, if
    !                             this yields no improvement it is advised
    !                             to analyze the integrand, in order to
    !                             determine the integration difficulties.
    !                             if the position of a local difficulty can
    !                             be determined (e.g. singularity,
    !                             discontinuity within the interval) one
    !                             will probably gain from splitting up the
    !                             interval at this point and calling the
    !                             integrator on the subranges. if possible,
    !                             an appropriate special-purpose integrator
    !                             should be used which is designed for
    !                             handling the type of difficulty involved.
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
    !                             tolerance cannot be achieved due to
    !                             roundoff in the extrapolation table, and
    !                             that the returned result is the best which
    !                             can be obtained.
    !                         = 5 the integral is probably divergent, or
    !                             slowly convergent. it must be noted that
    !                             divergence can occur with any other value
    !                             of ier > 0.
    !                         = 6 the input is invalid, because
    !                             epsabs < 0 and epsrel < 0,
    !                             or (integr /= 1 and integr /= 2) or
    !                             ICALL < 1 or maxp1 < 1.
    !                             result, abserr, neval, last, rlist(1),
    !                             elist(1), iord(1) and nnlog(1) are set to
    !                             zero. alist(1) and blist(1) are set to a
    !                             and b respectively.
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
    !                     estimates over the subintervals, such that
    !                     elist(iord(1)), ..., elist(iord(k)), form
    !                     a decreasing sequence, with k = last
    !                     if last <= (limit/2+2), and
    !                     k = limit+1-last otherwise.
    !
    !            nnlog  - integer
    !                     vector of dimension at least limit, indicating the
    !                     subdivision levels of the subintervals, i.e.
    !                     iwork(i) = l means that the subinterval numbered
    !                     i is of length abs(b-a)*2**(1-l)
    !
    !         on entry and return
    !            momcom - integer
    !                     indicating that the Chebyshev moments have been
    !                     computed for intervals of lengths
    !                     (abs(b-a))*2**(-l), l=0,1,2, ..., momcom-1,
    !                     momcom < maxp1
    !
    !            chebmo - real
    !                     array of dimension (maxp1,25) containing the
    !                     Chebyshev moments
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
    !                       been obtained it is put in rlist2(numrl2) after
    !                       numrl2 has been increased by one
    !           small     - length of the smallest interval considered
    !                       up to now, multiplied by 1.5
    !           erlarg    - sum of the errors over the intervals larger
    !                       than the smallest interval considered up to now
    !           extrap    - logical variable denoting that the routine is
    !                       attempting to perform extrapolation, i.e. before
    !                       subdividing the smallest interval we try to
    !                       decrease the value of erlarg
    !           noext     - logical variable denoting that extrapolation
    !                       is no longer allowed (true value)
    !
    implicit none

    integer limit
    integer maxp1

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
    real(8) :: chebmo(maxp1,25)
    real(8) :: correc
    real(8) :: defab1
    real(8) :: defab2
    real(8) :: defabs
    real(8) :: domega
    real(8) :: dres
    real(8) :: elist(limit)
    real(8) :: epsabs
    real(8) :: epsrel
    real(8) :: erlarg
    real(8) :: erlast
    real(8) :: errbnd
    real(8) :: errmax
    real(8) :: error1
    real(8) :: erro12
    real(8) :: error2
    real(8) :: errsum
    real(8) :: ertest
    logical extall
    logical extrap
    real(8), external :: f
    integer icall
    integer id
    integer ier
    integer ierro
    integer integr
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
    integer momcom
    integer nev
    integer neval
    integer nnlog(limit)
    logical noext
    integer nres
    integer nrmax
    integer nrmom
    integer numrl2
    real(8) :: omega
    real(8) :: resabs
    real(8) :: reseps
    real(8) :: result
    real(8) :: res3la(3)
    real(8) :: rlist(limit)
    real(8) :: rlist2(52)
    real(8) :: small
    real(8) :: width
    !
    !  the dimension of rlist2 is determined by  the value of
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
    iord(1) = 0
    nnlog(1) = 0

    if ( (integr /= 1.and.integr /= 2) .or. (epsabs < 0.0e+00.and. &
         epsrel < 0.0e+00) .or. icall < 1 .or. maxp1 < 1 ) then
       ier = 6
       return
    end if
    !
    !  First approximation to the integral.
    !
    domega = abs ( omega )
    nrmom = 0

    if ( icall <= 1 ) then
       momcom = 0
    end if

    call qc25o ( f, a, b, domega, integr, nrmom, maxp1, 0, result, abserr, &
         neval, defabs, resabs, momcom, chebmo )
    !
    !  Test on accuracy.
    !
    dres = abs(result)
    errbnd = max ( epsabs,epsrel*dres)
    rlist(1) = result
    elist(1) = abserr
    iord(1) = 1
    if ( abserr <= 1.0e+02* epsilon ( defabs ) *defabs .and. &
         abserr > errbnd ) ier = 2

    if ( limit == 1 ) then
       ier = 1
    end if

    if ( ier /= 0 .or. abserr <= errbnd ) then
       go to 200
    end if
    !
    !  Initializations
    !
    errmax = abserr
    maxerr = 1
    area = result
    errsum = abserr
    abserr = huge ( abserr )
    nrmax = 1
    extrap = .false.
    noext = .false.
    ierro = 0
    iroff1 = 0
    iroff2 = 0
    iroff3 = 0
    ktmin = 0
    small = abs(b-a)*7.5e-01
    nres = 0
    numrl2 = 0
    extall = .false.

    if ( 5.0e-01*abs(b-a)*domega <= 2.0e+00) then
       numrl2 = 1
       extall = .true.
       rlist2(1) = result
    end if

    if ( 2.5e-01 * abs(b-a) * domega <= 2.0e+00 ) then
       extall = .true.
    end if

    if ( dres >= (1.0e+00-5.0e+01* epsilon ( defabs ) )*defabs ) then
       ksgn = 1
    else
       ksgn = -1
    end if
    !
    !  main do-loop
    !
    do last = 2, limit
       !
       !  Bisect the subinterval with the nrmax-th largest error estimate.
       !
       nrmom = nnlog(maxerr)+1
       a1 = alist(maxerr)
       b1 = 5.0e-01*(alist(maxerr)+blist(maxerr))
       a2 = b1
       b2 = blist(maxerr)
       erlast = errmax

       call qc25o ( f, a1, b1, domega, integr, nrmom, maxp1, 0, area1, &
            error1, nev, resabs, defab1, momcom, chebmo )

       neval = neval+nev

       call qc25o ( f, a2, b2, domega, integr, nrmom, maxp1, 1, area2, &
            error2, nev, resabs, defab2, momcom, chebmo )

       neval = neval+nev
       !
       !  Improve previous approximations to integral and error and
       !  test for accuracy.
       !
       area12 = area1+area2
       erro12 = error1+error2
       errsum = errsum+erro12-errmax
       area = area+area12-rlist(maxerr)
       if ( defab1 == error1 .or. defab2 == error2 ) go to 25
       if ( abs(rlist(maxerr)-area12) > 1.0e-05*abs(area12) &
            .or. erro12 < 9.9e-01*errmax ) go to 20
       if ( extrap ) iroff2 = iroff2+1

       if ( .not.extrap ) then
          iroff1 = iroff1+1
       end if

20     continue

       if ( last > 10.and.erro12 > errmax ) iroff3 = iroff3+1

25     continue

       rlist(maxerr) = area1
       rlist(last) = area2
       nnlog(maxerr) = nrmom
       nnlog(last) = nrmom
       errbnd = max ( epsabs,epsrel*abs(area))
       !
       !  Test for roundoff error and eventually set error flag
       !
       if ( iroff1+iroff2 >= 10 .or. iroff3 >= 20 ) ier = 2

       if ( iroff2 >= 5) ierro = 3
       !
       !  Set error flag in the case that the number of subintervals
       !  equals limit.
       !
       if ( last == limit ) then
          ier = 1
       end if
       !
       !  Set error flag in the case of bad integrand behavior at
       !  a point of the integration range.
       !
       if ( max ( abs(a1),abs(b2)) <= (1.0e+00+1.0e+03* epsilon ( a1 ) ) &
            *(abs(a2)+1.0e+03* tiny ( a2 ) )) then
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
          go to 170
       end if

       if ( ier /= 0 ) then
          exit
       end if

       if ( last == 2 .and. extall ) go to 120

       if ( noext ) then
          cycle
       end if

       if ( .not. extall ) go to 50
       erlarg = erlarg-erlast
       if ( abs(b1-a1) > small ) erlarg = erlarg+erro12
       if ( extrap ) go to 70
       !
       !  Test whether the interval to be bisected next is the
       !  smallest interval.
       !
50     continue

       width = abs(blist(maxerr)-alist(maxerr))

       if ( width > small ) then
          cycle
       end if

       if ( extall ) go to 60
       !
       !  Test whether we can start with the extrapolation procedure
       !  (we do this if we integrate over the next interval with
       !  use of a Gauss-Kronrod rule - see QC25O).
       !
       small = small*5.0e-01

       if ( 2.5e-01*width*domega > 2.0e+00 ) then
          cycle
       end if

       extall = .true.
       go to 130

60     continue

       extrap = .true.
       nrmax = 2

70     continue

       if ( ierro == 3 .or. erlarg <= ertest ) go to 90
       !
       !  The smallest interval has the largest error.
       !  Before bisecting decrease the sum of the errors over the
       !  larger intervals (ERLARG) and perform extrapolation.
       !
       jupbnd = last

       if ( last > (limit/2+2) ) then
          jupbnd = limit+3-last
       end if

       id = nrmax

       do k = id, jupbnd
          maxerr = iord(nrmax)
          errmax = elist(maxerr)
          if ( abs(blist(maxerr)-alist(maxerr)) > small ) go to 140
          nrmax = nrmax+1
       end do
       !
       !  Perform extrapolation.
       !
90     continue

       numrl2 = numrl2+1
       rlist2(numrl2) = area

       if ( numrl2 < 3 ) go to 110

       call qextr ( numrl2, rlist2, reseps, abseps, res3la, nres )
       ktmin = ktmin+1

       if ( ktmin > 5.and.abserr < 1.0e-03*errsum ) then
          ier = 5
       end if

       if ( abseps >= abserr ) go to 100

       ktmin = 0
       abserr = abseps
       result = reseps
       correc = erlarg
       ertest = max ( epsabs, epsrel*abs(reseps))

       if ( abserr <= ertest ) then
          exit
       end if
       !
       !  Prepare bisection of the smallest interval.
       !
100    continue

       if ( numrl2 == 1 ) then
          noext = .true.
       end if

       if ( ier == 5 ) then
          exit
       end if

110    continue

       maxerr = iord(1)
       errmax = elist(maxerr)
       nrmax = 1
       extrap = .false.
       small = small*5.0e-01
       erlarg = errsum
       cycle

120    continue

       small = small * 5.0e-01
       numrl2 = numrl2 + 1
       rlist2(numrl2) = area

130    continue

       ertest = errbnd
       erlarg = errsum

140    continue

    end do
    !
    !  set the final result.
    !
    if ( abserr == huge ( abserr ) .or. nres == 0 ) then
       go to 170
    end if

    if ( ier+ierro == 0 ) go to 165
    if ( ierro == 3 ) abserr = abserr+correc
    if ( ier == 0 ) ier = 3
    if ( result /= 0.0e+00.and.area /= 0.0e+00 ) go to 160
    if ( abserr > errsum ) go to 170
    if ( area == 0.0e+00 ) go to 190
    go to 165

160 continue

    if ( abserr/abs(result) > errsum/abs(area) ) go to 170
    !
    !  Test on divergence.
    !
165 continue

    if ( ksgn == (-1) .and. max ( abs(result),abs(area)) <=  &
         defabs*1.0e-02 ) go to 190

    if ( 1.0e-02 > (result/area) .or. (result/area) > 1.0e+02 &
         .or. errsum >= abs(area) ) ier = 6

    go to 190
    !
    !  Compute global integral sum.
    !
170 continue

    result = sum ( rlist(1:last) )

    abserr = errsum

190 continue

    if (ier > 2) ier=ier-1

200 continue

    if ( integr == 2 .and. omega < 0.0e+00 ) then
       result = -result
    end if

    return
  end subroutine qfour
  subroutine qk15 ( f, a, b, result, abserr, resabs, resasc )

    !*****************************************************************************80
    !
    !! QK15 carries out a 15 point Gauss-Kronrod quadrature rule.
    !
    !  Discussion:
    !
    !    This routine approximates
    !      I = integral ( A <= X <= B ) F(X) dx
    !    with an error estimate, and
    !      J = integral ( A <= X <= B ) | F(X) | dx
    !
    !  Author:
    !
    !    Robert Piessens, Elise de Doncker-Kapenger, 
    !    Christian Ueberhuber, David Kahaner
    !
    !  Reference:
    !
    !    Robert Piessens, Elise de Doncker-Kapenger, 
    !    Christian Ueberhuber, David Kahaner,
    !    QUADPACK, a Subroutine Package for Automatic Integration,
    !    Springer Verlag, 1983
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
    !    Output, real(8) :: RESULT, the estimated value of the integral.
    !    RESULT is computed by applying the 15-point Kronrod rule (RESK) 
    !    obtained by optimal addition of abscissae to the 7-point Gauss rule 
    !    (RESG).
    !
    !    Output, real(8) :: ABSERR, an estimate of | I - RESULT |.
    !
    !    Output, real(8) :: RESABS, approximation to the integral of the absolute
    !    value of F.
    !
    !    Output, real(8) :: RESASC, approximation to the integral | F-I/(B-A) | 
    !    over [A,B].
    !
    !  Local Parameters:
    !
    !           the abscissae and weights are given for the interval (-1,1).
    !           because of symmetry only the positive abscissae and their
    !           corresponding weights are given.
    !
    !           xgk    - abscissae of the 15-point Kronrod rule
    !                    xgk(2), xgk(4), ...  abscissae of the 7-point
    !                    Gauss rule
    !                    xgk(1), xgk(3), ...  abscissae which are optimally
    !                    added to the 7-point Gauss rule
    !
    !           wgk    - weights of the 15-point Kronrod rule
    !
    !           wg     - weights of the 7-point Gauss rule
    !
    !           centr  - mid point of the interval
    !           hlgth  - half-length of the interval
    !           absc   - abscissa
    !           fval*  - function value
    !           resg   - result of the 7-point Gauss formula
    !           resk   - result of the 15-point Kronrod formula
    !           reskh  - approximation to the mean value of f over (a,b),
    !                    i.e. to i/(b-a)
    !
    implicit none

    real(8) :: a
    real(8) :: absc
    real(8) :: abserr
    real(8) :: b
    real(8) :: centr
    real(8) :: dhlgth
    real(8), external :: f
    real(8) :: fc
    real(8) :: fsum
    real(8) :: fval1
    real(8) :: fval2
    real(8) :: fv1(7)
    real(8) :: fv2(7)
    real(8) :: hlgth
    integer j
    integer jtw
    integer jtwm1
    real(8) :: resabs
    real(8) :: resasc
    real(8) :: resg
    real(8) :: resk
    real(8) :: reskh
    real(8) :: result
    real(8) :: wg(4)
    real(8) :: wgk(8)
    real(8) :: xgk(8)

    data xgk(1),xgk(2),xgk(3),xgk(4),xgk(5),xgk(6),xgk(7),xgk(8)/ &
         9.914553711208126e-01,   9.491079123427585e-01, &
         8.648644233597691e-01,   7.415311855993944e-01, &
         5.860872354676911e-01,   4.058451513773972e-01, &
         2.077849550078985e-01,   0.0e+00              /
    data wgk(1),wgk(2),wgk(3),wgk(4),wgk(5),wgk(6),wgk(7),wgk(8)/ &
         2.293532201052922e-02,   6.309209262997855e-02, &
         1.047900103222502e-01,   1.406532597155259e-01, &
         1.690047266392679e-01,   1.903505780647854e-01, &
         2.044329400752989e-01,   2.094821410847278e-01/
    data wg(1),wg(2),wg(3),wg(4)/ &
         1.294849661688697e-01,   2.797053914892767e-01, &
         3.818300505051189e-01,   4.179591836734694e-01/
    !
    centr = 5.0e-01*(a+b)
    hlgth = 5.0e-01*(b-a)
    dhlgth = abs(hlgth)
    !
    !  Compute the 15-point Kronrod approximation to the integral,
    !  and estimate the absolute error.
    !
    fc = f(centr)
    resg = fc*wg(4)
    resk = fc*wgk(8)
    resabs = abs(resk)

    do j = 1, 3
       jtw = j*2
       absc = hlgth*xgk(jtw)
       fval1 = f(centr-absc)
       fval2 = f(centr+absc)
       fv1(jtw) = fval1
       fv2(jtw) = fval2
       fsum = fval1+fval2
       resg = resg+wg(j)*fsum
       resk = resk+wgk(jtw)*fsum
       resabs = resabs+wgk(jtw)*(abs(fval1)+abs(fval2))
    end do

    do j = 1, 4
       jtwm1 = j*2-1
       absc = hlgth*xgk(jtwm1)
       fval1 = f(centr-absc)
       fval2 = f(centr+absc)
       fv1(jtwm1) = fval1
       fv2(jtwm1) = fval2
       fsum = fval1+fval2
       resk = resk+wgk(jtwm1)*fsum
       resabs = resabs+wgk(jtwm1)*(abs(fval1)+abs(fval2))
    end do

    reskh = resk * 5.0e-01
    resasc = wgk(8)*abs(fc-reskh)

    do j = 1, 7
       resasc = resasc+wgk(j)*(abs(fv1(j)-reskh)+abs(fv2(j)-reskh))
    end do

    result = resk*hlgth
    resabs = resabs*dhlgth
    resasc = resasc*dhlgth
    abserr = abs((resk-resg)*hlgth)

    if ( resasc /= 0.0e+00.and.abserr /= 0.0e+00 ) then
       abserr = resasc*min ( 1.0e+00,(2.0e+02*abserr/resasc)**1.5e+00)
    end if

    if ( resabs > tiny ( resabs ) / (5.0e+01* epsilon ( resabs ) ) ) then
       abserr = max (( epsilon ( resabs ) *5.0e+01)*resabs,abserr)
    end if

    return
  end subroutine qk15
  subroutine qk15i ( f, boun, inf, a, b, result, abserr, resabs, resasc )

    !*****************************************************************************80
    !
    !! QK15I applies a 15 point Gauss-Kronrod quadrature on an infinite interval.
    !
    !  Discussion:
    !
    !    The original infinite integration range is mapped onto the interval 
    !    (0,1) and (a,b) is a part of (0,1).  The routine then computes:
    !
    !    i = integral of transformed integrand over (a,b),
    !    j = integral of abs(transformed integrand) over (a,b).
    !
    !  Author:
    !
    !    Robert Piessens, Elise de Doncker-Kapenger, 
    !    Christian Ueberhuber, David Kahaner
    !
    !  Reference:
    !
    !    Robert Piessens, Elise de Doncker-Kapenger, 
    !    Christian Ueberhuber, David Kahaner,
    !    QUADPACK, a Subroutine Package for Automatic Integration,
    !    Springer Verlag, 1983
    !
    !  Parameters:
    !
    !    Input, external real(8) :: F, the name of the function routine, of the form
    !      function f ( x )
    !      real(8) :: f
    !      real(8) :: x
    !    which evaluates the integrand function.
    !
    !    Input, real(8) :: BOUN, the finite bound of the original integration range,
    !    or zero if INF is 2.
    !
    !    Input, integer INF, indicates the type of the interval.
    !    -1: the original interval is (-infinity,BOUN),
    !    +1, the original interval is (BOUN,+infinity),
    !    +2, the original interval is (-infinity,+infinity) and
    !    the integral is computed as the sum of two integrals, one 
    !    over (-infinity,0) and one over (0,+infinity).
    !
    !    Input, real(8) :: A, B, the limits of integration, over a subrange of [0,1].
    !
    !    Output, real(8) :: RESULT, the estimated value of the integral.
    !    RESULT is computed by applying the 15-point Kronrod rule (RESK) obtained 
    !    by optimal addition of abscissae to the 7-point Gauss rule (RESG).
    !
    !    Output, real(8) :: ABSERR, an estimate of | I - RESULT |.
    !
    !    Output, real(8) :: RESABS, approximation to the integral of the absolute
    !    value of F.
    !
    !    Output, real(8) :: RESASC, approximation to the integral of the
    !    transformated integrand | F-I/(B-A) | over [A,B].
    !
    !  Local Parameters:
    !
    !           centr  - mid point of the interval
    !           hlgth  - half-length of the interval
    !           absc*  - abscissa
    !           tabsc* - transformed abscissa
    !           fval*  - function value
    !           resg   - result of the 7-point Gauss formula
    !           resk   - result of the 15-point Kronrod formula
    !           reskh  - approximation to the mean value of the transformed
    !                    integrand over (a,b), i.e. to i/(b-a)
    !
    implicit none

    real(8) :: a
    real(8) :: absc
    real(8) :: absc1
    real(8) :: absc2
    real(8) :: abserr
    real(8) :: b
    real(8) :: boun
    real(8) :: centr
    real(8) :: dinf
    real(8), external :: f
    real(8) :: fc
    real(8) :: fsum
    real(8) :: fval1
    real(8) :: fval2
    real(8) :: fv1(7)
    real(8) :: fv2(7)
    real(8) :: hlgth
    integer inf
    integer j
    real(8) :: resabs
    real(8) :: resasc
    real(8) :: resg
    real(8) :: resk
    real(8) :: reskh
    real(8) :: result
    real(8) :: tabsc1
    real(8) :: tabsc2
    real(8) :: wg(8)
    real(8) :: wgk(8)
    real(8) :: xgk(8)
    !
    !  the abscissae and weights are supplied for the interval
    !  (-1,1).  because of symmetry only the positive abscissae and
    !  their corresponding weights are given.
    !
    !           xgk    - abscissae of the 15-point Kronrod rule
    !                    xgk(2), xgk(4), ... abscissae of the 7-point Gauss
    !                    rule
    !                    xgk(1), xgk(3), ...  abscissae which are optimally
    !                    added to the 7-point Gauss rule
    !
    !           wgk    - weights of the 15-point Kronrod rule
    !
    !           wg     - weights of the 7-point Gauss rule, corresponding
    !                    to the abscissae xgk(2), xgk(4), ...
    !                    wg(1), wg(3), ... are set to zero.
    !
    data xgk(1),xgk(2),xgk(3),xgk(4),xgk(5),xgk(6),xgk(7),xgk(8)/ &
         9.914553711208126e-01,     9.491079123427585e-01, &
         8.648644233597691e-01,     7.415311855993944e-01, &
         5.860872354676911e-01,     4.058451513773972e-01, &
         2.077849550078985e-01,     0.0000000000000000e+00/

    data wgk(1),wgk(2),wgk(3),wgk(4),wgk(5),wgk(6),wgk(7),wgk(8)/ &
         2.293532201052922e-02,     6.309209262997855e-02, &
         1.047900103222502e-01,     1.406532597155259e-01, &
         1.690047266392679e-01,     1.903505780647854e-01, &
         2.044329400752989e-01,     2.094821410847278e-01/

    data wg(1),wg(2),wg(3),wg(4),wg(5),wg(6),wg(7),wg(8)/ &
         0.0000000000000000e+00,     1.294849661688697e-01, &
         0.0000000000000000e+00,     2.797053914892767e-01, &
         0.0000000000000000e+00,     3.818300505051189e-01, &
         0.0000000000000000e+00,     4.179591836734694e-01/

    dinf = min ( 1, inf )

    centr = 5.0e-01*(a+b)
    hlgth = 5.0e-01*(b-a)
    tabsc1 = boun+dinf*(1.0e+00-centr)/centr
    fval1 = f(tabsc1)
    if ( inf == 2 ) fval1 = fval1+f(-tabsc1)
    fc = (fval1/centr)/centr
    !
    !  Compute the 15-point Kronrod approximation to the integral,
    !  and estimate the error.
    !
    resg = wg(8)*fc
    resk = wgk(8)*fc
    resabs = abs(resk)

    do j = 1, 7

       absc = hlgth*xgk(j)
       absc1 = centr-absc
       absc2 = centr+absc
       tabsc1 = boun+dinf*(1.0e+00-absc1)/absc1
       tabsc2 = boun+dinf*(1.0e+00-absc2)/absc2
       fval1 = f(tabsc1)
       fval2 = f(tabsc2)

       if ( inf == 2 ) then
          fval1 = fval1+f(-tabsc1)
          fval2 = fval2+f(-tabsc2)
       end if

       fval1 = (fval1/absc1)/absc1
       fval2 = (fval2/absc2)/absc2
       fv1(j) = fval1
       fv2(j) = fval2
       fsum = fval1+fval2
       resg = resg+wg(j)*fsum
       resk = resk+wgk(j)*fsum
       resabs = resabs+wgk(j)*(abs(fval1)+abs(fval2))
    end do

    reskh = resk * 5.0e-01
    resasc = wgk(8) * abs(fc-reskh)

    do j = 1, 7
       resasc = resasc + wgk(j)*(abs(fv1(j)-reskh)+abs(fv2(j)-reskh))
    end do

    result = resk * hlgth
    resasc = resasc * hlgth
    resabs = resabs * hlgth
    abserr = abs ( ( resk - resg ) * hlgth )

    if ( resasc /= 0.0e+00.and.abserr /= 0.0e+00) then
       abserr = resasc* min ( 1.0e+00,(2.0e+02*abserr/resasc)**1.5e+00)
    end if

    if ( resabs > tiny ( resabs ) / ( 5.0e+01 * epsilon ( resabs ) ) ) then
       abserr = max (( epsilon ( resabs ) *5.0e+01)*resabs,abserr)
    end if

    return
  end subroutine qk15i
  subroutine qk15w ( f, w, p1, p2, p3, p4, kp, a, b, result, abserr, resabs, &
       resasc )

    !*****************************************************************************80
    !
    !! QK15W applies a 15 point Gauss-Kronrod rule for a weighted integrand.
    !
    !  Discussion:
    !
    !    This routine approximates 
    !      i = integral of f*w over (a,b), 
    !    with error estimate, and
    !      j = integral of abs(f*w) over (a,b)
    !
    !  Author:
    !
    !    Robert Piessens, Elise de Doncker-Kapenger, 
    !    Christian Ueberhuber, David Kahaner
    !
    !  Reference:
    !
    !    Robert Piessens, Elise de Doncker-Kapenger, 
    !    Christian Ueberhuber, David Kahaner,
    !    QUADPACK, a Subroutine Package for Automatic Integration,
    !    Springer Verlag, 1983
    !
    !  Parameters:
    !
    !    Input, external real(8) :: F, the name of the function routine, of the form
    !      function f ( x )
    !      real(8) :: f
    !      real(8) :: x
    !    which evaluates the integrand function.
    !
    !              w      - real
    !                       function subprogram defining the integrand
    !                       weight function w(x). the actual name for w
    !                       needs to be declared e x t e r n a l in the
    !                       calling program.
    !
    !    ?, real(8) :: P1, P2, P3, P4, parameters in the weight function
    !
    !    Input, integer KP, key for indicating the type of weight function
    !
    !    Input, real(8) :: A, B, the limits of integration.
    !
    !    Output, real(8) :: RESULT, the estimated value of the integral.
    !    RESULT is computed by applying the 15-point Kronrod rule (RESK) obtained by
    !    optimal addition of abscissae to the 7-point Gauss rule (RESG).
    !
    !    Output, real(8) :: ABSERR, an estimate of | I - RESULT |.
    !
    !    Output, real(8) :: RESABS, approximation to the integral of the absolute
    !    value of F.
    !
    !    Output, real(8) :: RESASC, approximation to the integral | F-I/(B-A) | 
    !    over [A,B].
    !
    !  Local Parameters:
    !
    !           centr  - mid point of the interval
    !           hlgth  - half-length of the interval
    !           absc*  - abscissa
    !           fval*  - function value
    !           resg   - result of the 7-point Gauss formula
    !           resk   - result of the 15-point Kronrod formula
    !           reskh  - approximation to the mean value of f*w over (a,b),
    !                    i.e. to i/(b-a)
    !
    implicit none

    real(8) :: a
    real(8) :: absc
    real(8) :: absc1
    real(8) :: absc2
    real(8) :: abserr
    real(8) :: b
    real(8) :: centr
    real(8) :: dhlgth
    real(8), external :: f
    real(8) :: fc
    real(8) :: fsum
    real(8) :: fval1
    real(8) :: fval2
    real(8) :: fv1(7)
    real(8) :: fv2(7)
    real(8) :: hlgth
    integer j
    integer jtw
    integer jtwm1
    integer kp
    real(8) :: p1
    real(8) :: p2
    real(8) :: p3
    real(8) :: p4
    real(8) :: resabs
    real(8) :: resasc
    real(8) :: resg
    real(8) :: resk
    real(8) :: reskh
    real(8) :: result
    real(8), external :: w
    real(8), dimension ( 4 ) :: wg = (/ &
         1.294849661688697e-01,     2.797053914892767e-01, &
         3.818300505051889e-01,     4.179591836734694e-01 /)
    real(8) :: wgk(8)
    real(8) :: xgk(8)
    !
    !  the abscissae and weights are given for the interval (-1,1).
    !  because of symmetry only the positive abscissae and their
    !  corresponding weights are given.
    !
    !           xgk    - abscissae of the 15-point Gauss-Kronrod rule
    !                    xgk(2), xgk(4), ... abscissae of the 7-point Gauss
    !                    rule
    !                    xgk(1), xgk(3), ... abscissae which are optimally
    !                    added to the 7-point Gauss rule
    !
    !           wgk    - weights of the 15-point Gauss-Kronrod rule
    !
    !           wg     - weights of the 7-point Gauss rule
    !
    data xgk(1),xgk(2),xgk(3),xgk(4),xgk(5),xgk(6),xgk(7),xgk(8)/ &
         9.914553711208126e-01,     9.491079123427585e-01, &
         8.648644233597691e-01,     7.415311855993944e-01, &
         5.860872354676911e-01,     4.058451513773972e-01, &
         2.077849550789850e-01,     0.000000000000000e+00/

    data wgk(1),wgk(2),wgk(3),wgk(4),wgk(5),wgk(6),wgk(7),wgk(8)/ &
         2.293532201052922e-02,     6.309209262997855e-02, &
         1.047900103222502e-01,     1.406532597155259e-01, &
         1.690047266392679e-01,     1.903505780647854e-01, &
         2.044329400752989e-01,     2.094821410847278e-01/
    !
    centr = 5.0e-01*(a+b)
    hlgth = 5.0e-01*(b-a)
    dhlgth = abs(hlgth)
    !
    !  Compute the 15-point Kronrod approximation to the integral,
    !  and estimate the error.
    !
    fc = f(centr)*w(centr,p1,p2,p3,p4,kp)
    resg = wg(4)*fc
    resk = wgk(8)*fc
    resabs = abs(resk)

    do j = 1, 3
       jtw = j*2
       absc = hlgth*xgk(jtw)
       absc1 = centr-absc
       absc2 = centr+absc
       fval1 = f(absc1)*w(absc1,p1,p2,p3,p4,kp)
       fval2 = f(absc2)*w(absc2,p1,p2,p3,p4,kp)
       fv1(jtw) = fval1
       fv2(jtw) = fval2
       fsum = fval1+fval2
       resg = resg+wg(j)*fsum
       resk = resk+wgk(jtw)*fsum
       resabs = resabs+wgk(jtw)*(abs(fval1)+abs(fval2))
    end do

    do j = 1, 4
       jtwm1 = j*2-1
       absc = hlgth*xgk(jtwm1)
       absc1 = centr-absc
       absc2 = centr+absc
       fval1 = f(absc1)*w(absc1,p1,p2,p3,p4,kp)
       fval2 = f(absc2)*w(absc2,p1,p2,p3,p4,kp)
       fv1(jtwm1) = fval1
       fv2(jtwm1) = fval2
       fsum = fval1+fval2
       resk = resk+wgk(jtwm1)*fsum
       resabs = resabs+wgk(jtwm1)*(abs(fval1)+abs(fval2))
    end do

    reskh = resk*5.0e-01
    resasc = wgk(8)*abs(fc-reskh)

    do j = 1, 7
       resasc = resasc+wgk(j)*(abs(fv1(j)-reskh)+abs(fv2(j)-reskh))
    end do

    result = resk*hlgth
    resabs = resabs*dhlgth
    resasc = resasc*dhlgth
    abserr = abs((resk-resg)*hlgth)

    if ( resasc /= 0.0e+00.and.abserr /= 0.0e+00) then
       abserr = resasc*min ( 1.0e+00,(2.0e+02*abserr/resasc)**1.5e+00)
    end if

    if ( resabs > tiny ( resabs ) /(5.0e+01* epsilon ( resabs ) ) ) then
       abserr = max ( ( epsilon ( resabs ) * 5.0e+01)*resabs,abserr)
    end if

    return
  end subroutine qk15w
  subroutine qk21 ( f, a, b, result, abserr, resabs, resasc )

    !*****************************************************************************80
    !
    !! QK21 carries out a 21 point Gauss-Kronrod quadrature rule.
    !
    !  Discussion:
    !
    !    This routine approximates
    !      I = integral ( A <= X <= B ) F(X) dx
    !    with an error estimate, and
    !      J = integral ( A <= X <= B ) | F(X) | dx
    !
    !  Author:
    !
    !    Robert Piessens, Elise de Doncker-Kapenger, 
    !    Christian Ueberhuber, David Kahaner
    !
    !  Reference:
    !
    !    Robert Piessens, Elise de Doncker-Kapenger, 
    !    Christian Ueberhuber, David Kahaner,
    !    QUADPACK, a Subroutine Package for Automatic Integration,
    !    Springer Verlag, 1983
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
    !    Output, real(8) :: RESULT, the estimated value of the integral.
    !    RESULT is computed by applying the 21-point Kronrod rule (resk) 
    !    obtained by optimal addition of abscissae to the 10-point Gauss 
    !    rule (resg).
    !
    !    Output, real(8) :: ABSERR, an estimate of | I - RESULT |.
    !
    !    Output, real(8) :: RESABS, approximation to the integral of the absolute
    !    value of F.
    !
    !    Output, real(8) :: RESASC, approximation to the integral | F-I/(B-A) | 
    !    over [A,B].
    !
    implicit none

    real(8) :: a
    real(8) :: absc
    real(8) :: abserr
    real(8) :: b
    real(8) :: centr
    real(8) :: dhlgth
    real(8), external :: f
    real(8) :: fc
    real(8) :: fsum
    real(8) :: fval1
    real(8) :: fval2
    real(8) :: fv1(10)
    real(8) :: fv2(10)
    real(8) :: hlgth
    integer j
    integer jtw
    integer jtwm1
    real(8) :: resabs
    real(8) :: resasc
    real(8) :: resg
    real(8) :: resk
    real(8) :: reskh
    real(8) :: result
    real(8) :: wg(5)
    real(8) :: wgk(11)
    real(8) :: xgk(11)
    !
    !           the abscissae and weights are given for the interval (-1,1).
    !           because of symmetry only the positive abscissae and their
    !           corresponding weights are given.
    !
    !           xgk    - abscissae of the 21-point Kronrod rule
    !                    xgk(2), xgk(4), ...  abscissae of the 10-point
    !                    Gauss rule
    !                    xgk(1), xgk(3), ...  abscissae which are optimally
    !                    added to the 10-point Gauss rule
    !
    !           wgk    - weights of the 21-point Kronrod rule
    !
    !           wg     - weights of the 10-point Gauss rule
    !
    data xgk(1),xgk(2),xgk(3),xgk(4),xgk(5),xgk(6),xgk(7),xgk(8), &
         xgk(9),xgk(10),xgk(11)/ &
         9.956571630258081e-01,     9.739065285171717e-01, &
         9.301574913557082e-01,     8.650633666889845e-01, &
         7.808177265864169e-01,     6.794095682990244e-01, &
         5.627571346686047e-01,     4.333953941292472e-01, &
         2.943928627014602e-01,     1.488743389816312e-01, &
         0.000000000000000e+00/
    !
    data wgk(1),wgk(2),wgk(3),wgk(4),wgk(5),wgk(6),wgk(7),wgk(8), &
         wgk(9),wgk(10),wgk(11)/ &
         1.169463886737187e-02,     3.255816230796473e-02, &
         5.475589657435200e-02,     7.503967481091995e-02, &
         9.312545458369761e-02,     1.093871588022976e-01, &
         1.234919762620659e-01,     1.347092173114733e-01, &
         1.427759385770601e-01,     1.477391049013385e-01, &
         1.494455540029169e-01/
    !
    data wg(1),wg(2),wg(3),wg(4),wg(5)/ &
         6.667134430868814e-02,     1.494513491505806e-01, &
         2.190863625159820e-01,     2.692667193099964e-01, &
         2.955242247147529e-01/
    !
    !
    !           list of major variables
    !
    !           centr  - mid point of the interval
    !           hlgth  - half-length of the interval
    !           absc   - abscissa
    !           fval*  - function value
    !           resg   - result of the 10-point Gauss formula
    !           resk   - result of the 21-point Kronrod formula
    !           reskh  - approximation to the mean value of f over (a,b),
    !                    i.e. to i/(b-a)
    !
    centr = 5.0e-01*(a+b)
    hlgth = 5.0e-01*(b-a)
    dhlgth = abs(hlgth)
    !
    !  Compute the 21-point Kronrod approximation to the
    !  integral, and estimate the absolute error.
    !
    resg = 0.0e+00
    fc = f(centr)
    resk = wgk(11)*fc
    resabs = abs(resk)

    do j = 1, 5
       jtw = 2*j
       absc = hlgth*xgk(jtw)
       fval1 = f(centr-absc)
       fval2 = f(centr+absc)
       fv1(jtw) = fval1
       fv2(jtw) = fval2
       fsum = fval1+fval2
       resg = resg+wg(j)*fsum
       resk = resk+wgk(jtw)*fsum
       resabs = resabs+wgk(jtw)*(abs(fval1)+abs(fval2))
    end do

    do j = 1, 5
       jtwm1 = 2*j-1
       absc = hlgth*xgk(jtwm1)
       fval1 = f(centr-absc)
       fval2 = f(centr+absc)
       fv1(jtwm1) = fval1
       fv2(jtwm1) = fval2
       fsum = fval1+fval2
       resk = resk+wgk(jtwm1)*fsum
       resabs = resabs+wgk(jtwm1)*(abs(fval1)+abs(fval2))
    end do

    reskh = resk*5.0e-01
    resasc = wgk(11)*abs(fc-reskh)

    do j = 1, 10
       resasc = resasc+wgk(j)*(abs(fv1(j)-reskh)+abs(fv2(j)-reskh))
    end do

    result = resk*hlgth
    resabs = resabs*dhlgth
    resasc = resasc*dhlgth
    abserr = abs((resk-resg)*hlgth)

    if ( resasc /= 0.0e+00.and.abserr /= 0.0e+00) then
       abserr = resasc*min ( 1.0e+00,(2.0e+02*abserr/resasc)**1.5e+00)
    end if

    if ( resabs > tiny ( resabs ) /(5.0e+01* epsilon ( resabs ) )) then
       abserr = max (( epsilon ( resabs ) *5.0e+01)*resabs,abserr)
    end if

    return
  end subroutine qk21
  subroutine qk31 ( f, a, b, result, abserr, resabs, resasc )

    !*****************************************************************************80
    !
    !! QK31 carries out a 31 point Gauss-Kronrod quadrature rule.
    !
    !  Discussion:
    !
    !    This routine approximates
    !      I = integral ( A <= X <= B ) F(X) dx
    !    with an error estimate, and
    !      J = integral ( A <= X <= B ) | F(X) | dx
    !
    !  Author:
    !
    !    Robert Piessens, Elise de Doncker-Kapenger, 
    !    Christian Ueberhuber, David Kahaner
    !
    !  Reference:
    !
    !    Robert Piessens, Elise de Doncker-Kapenger, 
    !    Christian Ueberhuber, David Kahaner,
    !    QUADPACK, a Subroutine Package for Automatic Integration,
    !    Springer Verlag, 1983
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
    !    Output, real(8) :: RESULT, the estimated value of the integral.
    !                       result is computed by applying the 31-point
    !                       Gauss-Kronrod rule (resk), obtained by optimal
    !                       addition of abscissae to the 15-point Gauss
    !                       rule (resg).
    !
    !    Output, real(8) :: ABSERR, an estimate of | I - RESULT |.
    !
    !    Output, real(8) :: RESABS, approximation to the integral of the absolute
    !    value of F.
    !
    !    Output, real(8) :: RESASC, approximation to the integral | F-I/(B-A) | 
    !    over [A,B].
    !
    implicit none

    real(8) :: a
    real(8) :: absc
    real(8) :: abserr
    real(8) :: b
    real(8) :: centr
    real(8) :: dhlgth
    real(8), external :: f
    real(8) :: fc
    real(8) :: fsum
    real(8) :: fval1
    real(8) :: fval2
    real(8) :: fv1(15)
    real(8) :: fv2(15)
    real(8) :: hlgth
    integer j
    integer jtw
    integer jtwm1
    real(8) :: resabs
    real(8) :: resasc
    real(8) :: resg
    real(8) :: resk
    real(8) :: reskh
    real(8) :: result
    real(8) :: wg(8)
    real(8) :: wgk(16)
    real(8) :: xgk(16)
    !
    !           the abscissae and weights are given for the interval (-1,1).
    !           because of symmetry only the positive abscissae and their
    !           corresponding weights are given.
    !
    !           xgk    - abscissae of the 31-point Kronrod rule
    !                    xgk(2), xgk(4), ...  abscissae of the 15-point
    !                    Gauss rule
    !                    xgk(1), xgk(3), ...  abscissae which are optimally
    !                    added to the 15-point Gauss rule
    !
    !           wgk    - weights of the 31-point Kronrod rule
    !
    !           wg     - weights of the 15-point Gauss rule
    !
    data xgk(1),xgk(2),xgk(3),xgk(4),xgk(5),xgk(6),xgk(7),xgk(8), &
         xgk(9),xgk(10),xgk(11),xgk(12),xgk(13),xgk(14),xgk(15),xgk(16)/ &
         9.980022986933971e-01,   9.879925180204854e-01, &
         9.677390756791391e-01,   9.372733924007059e-01, &
         8.972645323440819e-01,   8.482065834104272e-01, &
         7.904185014424659e-01,   7.244177313601700e-01, &
         6.509967412974170e-01,   5.709721726085388e-01, &
         4.850818636402397e-01,   3.941513470775634e-01, &
         2.991800071531688e-01,   2.011940939974345e-01, &
         1.011420669187175e-01,   0.0e+00               /
    data wgk(1),wgk(2),wgk(3),wgk(4),wgk(5),wgk(6),wgk(7),wgk(8), &
         wgk(9),wgk(10),wgk(11),wgk(12),wgk(13),wgk(14),wgk(15),wgk(16)/ &
         5.377479872923349e-03,   1.500794732931612e-02, &
         2.546084732671532e-02,   3.534636079137585e-02, &
         4.458975132476488e-02,   5.348152469092809e-02, &
         6.200956780067064e-02,   6.985412131872826e-02, &
         7.684968075772038e-02,   8.308050282313302e-02, &
         8.856444305621177e-02,   9.312659817082532e-02, &
         9.664272698362368e-02,   9.917359872179196e-02, &
         1.007698455238756e-01,   1.013300070147915e-01/
    data wg(1),wg(2),wg(3),wg(4),wg(5),wg(6),wg(7),wg(8)/ &
         3.075324199611727e-02,   7.036604748810812e-02, &
         1.071592204671719e-01,   1.395706779261543e-01, &
         1.662692058169939e-01,   1.861610000155622e-01, &
         1.984314853271116e-01,   2.025782419255613e-01/
    !
    !
    !           list of major variables
    !
    !           centr  - mid point of the interval
    !           hlgth  - half-length of the interval
    !           absc   - abscissa
    !           fval*  - function value
    !           resg   - result of the 15-point Gauss formula
    !           resk   - result of the 31-point Kronrod formula
    !           reskh  - approximation to the mean value of f over (a,b),
    !                    i.e. to i/(b-a)
    !
    centr = 5.0e-01*(a+b)
    hlgth = 5.0e-01*(b-a)
    dhlgth = abs(hlgth)
    !
    !  Compute the 31-point Kronrod approximation to the integral,
    !  and estimate the absolute error.
    !
    fc = f(centr)
    resg = wg(8)*fc
    resk = wgk(16)*fc
    resabs = abs(resk)

    do j = 1, 7
       jtw = j*2
       absc = hlgth*xgk(jtw)
       fval1 = f(centr-absc)
       fval2 = f(centr+absc)
       fv1(jtw) = fval1
       fv2(jtw) = fval2
       fsum = fval1+fval2
       resg = resg+wg(j)*fsum
       resk = resk+wgk(jtw)*fsum
       resabs = resabs+wgk(jtw)*(abs(fval1)+abs(fval2))
    end do

    do j = 1, 8
       jtwm1 = j*2-1
       absc = hlgth*xgk(jtwm1)
       fval1 = f(centr-absc)
       fval2 = f(centr+absc)
       fv1(jtwm1) = fval1
       fv2(jtwm1) = fval2
       fsum = fval1+fval2
       resk = resk+wgk(jtwm1)*fsum
       resabs = resabs+wgk(jtwm1)*(abs(fval1)+abs(fval2))
    end do

    reskh = resk*5.0e-01
    resasc = wgk(16)*abs(fc-reskh)

    do j = 1, 15
       resasc = resasc+wgk(j)*(abs(fv1(j)-reskh)+abs(fv2(j)-reskh))
    end do

    result = resk*hlgth
    resabs = resabs*dhlgth
    resasc = resasc*dhlgth
    abserr = abs((resk-resg)*hlgth)

    if ( resasc /= 0.0e+00.and.abserr /= 0.0e+00) &
         abserr = resasc*min ( 1.0e+00,(2.0e+02*abserr/resasc)**1.5e+00)

    if ( resabs > tiny ( resabs ) /(5.0e+01* epsilon ( resabs ) )) then
       abserr = max (( epsilon ( resabs ) *5.0e+01)*resabs,abserr)
    end if

    return
  end subroutine qk31
  subroutine qk41 ( f, a, b, result, abserr, resabs, resasc )

    !*****************************************************************************80
    !
    !! QK41 carries out a 41 point Gauss-Kronrod quadrature rule.
    !
    !  Discussion:
    !
    !    This routine approximates
    !      I = integral ( A <= X <= B ) F(X) dx
    !    with an error estimate, and
    !      J = integral ( A <= X <= B ) | F(X) | dx
    !
    !  Author:
    !
    !    Robert Piessens, Elise de Doncker-Kapenger, 
    !    Christian Ueberhuber, David Kahaner
    !
    !  Reference:
    !
    !    Robert Piessens, Elise de Doncker-Kapenger, 
    !    Christian Ueberhuber, David Kahaner,
    !    QUADPACK, a Subroutine Package for Automatic Integration,
    !    Springer Verlag, 1983
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
    !    Output, real(8) :: RESULT, the estimated value of the integral.
    !                       result is computed by applying the 41-point
    !                       Gauss-Kronrod rule (resk) obtained by optimal
    !                       addition of abscissae to the 20-point Gauss
    !                       rule (resg).
    !
    !    Output, real(8) :: ABSERR, an estimate of | I - RESULT |.
    !
    !    Output, real(8) :: RESABS, approximation to the integral of the absolute
    !    value of F.
    !
    !    Output, real(8) :: RESASC, approximation to the integral | F-I/(B-A) | 
    !    over [A,B].
    !
    !  Local Parameters:
    !
    !           centr  - mid point of the interval
    !           hlgth  - half-length of the interval
    !           absc   - abscissa
    !           fval*  - function value
    !           resg   - result of the 20-point Gauss formula
    !           resk   - result of the 41-point Kronrod formula
    !           reskh  - approximation to mean value of f over (a,b), i.e.
    !                    to i/(b-a)
    !
    implicit none

    real(8) :: a
    real(8) :: absc
    real(8) :: abserr
    real(8) :: b
    real(8) :: centr
    real(8) :: dhlgth
    real(8), external :: f
    real(8) :: fc
    real(8) :: fsum
    real(8) :: fval1
    real(8) :: fval2
    real(8) :: fv1(20)
    real(8) :: fv2(20)
    real(8) :: hlgth
    integer j
    integer jtw
    integer jtwm1
    real(8) :: resabs
    real(8) :: resasc
    real(8) :: resg
    real(8) :: resk
    real(8) :: reskh
    real(8) :: result
    real(8) :: wg(10)
    real(8) :: wgk(21)
    real(8) :: xgk(21)
    !
    !           the abscissae and weights are given for the interval (-1,1).
    !           because of symmetry only the positive abscissae and their
    !           corresponding weights are given.
    !
    !           xgk    - abscissae of the 41-point Gauss-Kronrod rule
    !                    xgk(2), xgk(4), ...  abscissae of the 20-point
    !                    Gauss rule
    !                    xgk(1), xgk(3), ...  abscissae which are optimally
    !                    added to the 20-point Gauss rule
    !
    !           wgk    - weights of the 41-point Gauss-Kronrod rule
    !
    !           wg     - weights of the 20-point Gauss rule
    !
    data xgk(1),xgk(2),xgk(3),xgk(4),xgk(5),xgk(6),xgk(7),xgk(8), &
         xgk(9),xgk(10),xgk(11),xgk(12),xgk(13),xgk(14),xgk(15),xgk(16), &
         xgk(17),xgk(18),xgk(19),xgk(20),xgk(21)/ &
         9.988590315882777e-01,   9.931285991850949e-01, &
         9.815078774502503e-01,   9.639719272779138e-01, &
         9.408226338317548e-01,   9.122344282513259e-01, &
         8.782768112522820e-01,   8.391169718222188e-01, &
         7.950414288375512e-01,   7.463319064601508e-01, &
         6.932376563347514e-01,   6.360536807265150e-01, &
         5.751404468197103e-01,   5.108670019508271e-01, &
         4.435931752387251e-01,   3.737060887154196e-01, &
         3.016278681149130e-01,   2.277858511416451e-01, &
         1.526054652409227e-01,   7.652652113349733e-02, &
         0.0e+00               /
    data wgk(1),wgk(2),wgk(3),wgk(4),wgk(5),wgk(6),wgk(7),wgk(8), &
         wgk(9),wgk(10),wgk(11),wgk(12),wgk(13),wgk(14),wgk(15),wgk(16), &
         wgk(17),wgk(18),wgk(19),wgk(20),wgk(21)/ &
         3.073583718520532e-03,   8.600269855642942e-03, &
         1.462616925697125e-02,   2.038837346126652e-02, &
         2.588213360495116e-02,   3.128730677703280e-02, &
         3.660016975820080e-02,   4.166887332797369e-02, &
         4.643482186749767e-02,   5.094457392372869e-02, &
         5.519510534828599e-02,   5.911140088063957e-02, &
         6.265323755478117e-02,   6.583459713361842e-02, &
         6.864867292852162e-02,   7.105442355344407e-02, &
         7.303069033278667e-02,   7.458287540049919e-02, &
         7.570449768455667e-02,   7.637786767208074e-02, &
         7.660071191799966e-02/
    data wg(1),wg(2),wg(3),wg(4),wg(5),wg(6),wg(7),wg(8),wg(9),wg(10)/ &
         1.761400713915212e-02,   4.060142980038694e-02, &
         6.267204833410906e-02,   8.327674157670475e-02, &
         1.019301198172404e-01,   1.181945319615184e-01, &
         1.316886384491766e-01,   1.420961093183821e-01, &
         1.491729864726037e-01,   1.527533871307259e-01/
    !
    centr = 5.0e-01*(a+b)
    hlgth = 5.0e-01*(b-a)
    dhlgth = abs(hlgth)
    !
    !  Compute 41-point Gauss-Kronrod approximation to the
    !  the integral, and estimate the absolute error.
    !
    resg = 0.0e+00
    fc = f(centr)
    resk = wgk(21)*fc
    resabs = abs(resk)

    do j = 1, 10
       jtw = j*2
       absc = hlgth*xgk(jtw)
       fval1 = f(centr-absc)
       fval2 = f(centr+absc)
       fv1(jtw) = fval1
       fv2(jtw) = fval2
       fsum = fval1+fval2
       resg = resg+wg(j)*fsum
       resk = resk+wgk(jtw)*fsum
       resabs = resabs+wgk(jtw)*(abs(fval1)+abs(fval2))
    end do

    do j = 1, 10
       jtwm1 = j*2-1
       absc = hlgth*xgk(jtwm1)
       fval1 = f(centr-absc)
       fval2 = f(centr+absc)
       fv1(jtwm1) = fval1
       fv2(jtwm1) = fval2
       fsum = fval1+fval2
       resk = resk+wgk(jtwm1)*fsum
       resabs = resabs+wgk(jtwm1)*(abs(fval1)+abs(fval2))
    end do

    reskh = resk*5.0e-01
    resasc = wgk(21)*abs(fc-reskh)

    do j = 1, 20
       resasc = resasc+wgk(j)*(abs(fv1(j)-reskh)+abs(fv2(j)-reskh))
    end do

    result = resk*hlgth
    resabs = resabs*dhlgth
    resasc = resasc*dhlgth
    abserr = abs((resk-resg)*hlgth)

    if ( resasc /= 0.0e+00.and.abserr /= 0.0e+00) &
         abserr = resasc*min ( 1.0e+00,(2.0e+02*abserr/resasc)**1.5e+00)

    if ( resabs > tiny ( resabs ) /(5.0e+01* epsilon ( resabs ) )) then
       abserr = max (( epsilon ( resabs ) *5.0e+01)*resabs,abserr)
    end if

    return
  end subroutine qk41
  subroutine qk51 ( f, a, b, result, abserr, resabs, resasc )

    !*****************************************************************************80
    !
    !! QK51 carries out a 51 point Gauss-Kronrod quadrature rule.
    !
    !  Discussion:
    !
    !    This routine approximates
    !      I = integral ( A <= X <= B ) F(X) dx
    !    with an error estimate, and
    !      J = integral ( A <= X <= B ) | F(X) | dx
    !
    !  Author:
    !
    !    Robert Piessens, Elise de Doncker-Kapenger, 
    !    Christian Ueberhuber, David Kahaner
    !
    !  Reference:
    !
    !    Robert Piessens, Elise de Doncker-Kapenger, 
    !    Christian Ueberhuber, David Kahaner,
    !    QUADPACK, a Subroutine Package for Automatic Integration,
    !    Springer Verlag, 1983
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
    !    Output, real(8) :: RESULT, the estimated value of the integral.
    !                       result is computed by applying the 51-point
    !                       Kronrod rule (resk) obtained by optimal addition
    !                       of abscissae to the 25-point Gauss rule (resg).
    !
    !    Output, real(8) :: ABSERR, an estimate of | I - RESULT |.
    !
    !    Output, real(8) :: RESABS, approximation to the integral of the absolute
    !    value of F.
    !
    !    Output, real(8) :: RESASC, approximation to the integral | F-I/(B-A) | 
    !    over [A,B].
    !
    !  Local Parameters:
    !
    !           centr  - mid point of the interval
    !           hlgth  - half-length of the interval
    !           absc   - abscissa
    !           fval*  - function value
    !           resg   - result of the 25-point Gauss formula
    !           resk   - result of the 51-point Kronrod formula
    !           reskh  - approximation to the mean value of f over (a,b),
    !                    i.e. to i/(b-a)
    !
    implicit none

    real(8) :: a
    real(8) :: absc
    real(8) :: abserr
    real(8) :: b
    real(8) :: centr
    real(8) :: dhlgth
    real(8), external :: f
    real(8) :: fc
    real(8) :: fsum
    real(8) :: fval1
    real(8) :: fval2
    real(8) :: fv1(25)
    real(8) :: fv2(25)
    real(8) :: hlgth
    integer j
    integer jtw
    integer jtwm1
    real(8) :: resabs
    real(8) :: resasc
    real(8) :: resg
    real(8) :: resk
    real(8) :: reskh
    real(8) :: result
    real(8) :: wg(13)
    real(8) :: wgk(26)
    real(8) :: xgk(26)
    !
    !           the abscissae and weights are given for the interval (-1,1).
    !           because of symmetry only the positive abscissae and their
    !           corresponding weights are given.
    !
    !           xgk    - abscissae of the 51-point Kronrod rule
    !                    xgk(2), xgk(4), ...  abscissae of the 25-point
    !                    Gauss rule
    !                    xgk(1), xgk(3), ...  abscissae which are optimally
    !                    added to the 25-point Gauss rule
    !
    !           wgk    - weights of the 51-point Kronrod rule
    !
    !           wg     - weights of the 25-point Gauss rule
    !
    data xgk(1),xgk(2),xgk(3),xgk(4),xgk(5),xgk(6),xgk(7),xgk(8), &
         xgk(9),xgk(10),xgk(11),xgk(12),xgk(13),xgk(14)/ &
         9.992621049926098e-01,   9.955569697904981e-01, &
         9.880357945340772e-01,   9.766639214595175e-01, &
         9.616149864258425e-01,   9.429745712289743e-01, &
         9.207471152817016e-01,   8.949919978782754e-01, &
         8.658470652932756e-01,   8.334426287608340e-01, &
         7.978737979985001e-01,   7.592592630373576e-01, &
         7.177664068130844e-01,   6.735663684734684e-01/
    data xgk(15),xgk(16),xgk(17),xgk(18),xgk(19),xgk(20),xgk(21), &
         xgk(22),xgk(23),xgk(24),xgk(25),xgk(26)/ &
         6.268100990103174e-01,   5.776629302412230e-01, &
         5.263252843347192e-01,   4.730027314457150e-01, &
         4.178853821930377e-01,   3.611723058093878e-01, &
         3.030895389311078e-01,   2.438668837209884e-01, &
         1.837189394210489e-01,   1.228646926107104e-01, &
         6.154448300568508e-02,   0.0e+00               /
    data wgk(1),wgk(2),wgk(3),wgk(4),wgk(5),wgk(6),wgk(7),wgk(8), &
         wgk(9),wgk(10),wgk(11),wgk(12),wgk(13),wgk(14)/ &
         1.987383892330316e-03,   5.561932135356714e-03, &
         9.473973386174152e-03,   1.323622919557167e-02, &
         1.684781770912830e-02,   2.043537114588284e-02, &
         2.400994560695322e-02,   2.747531758785174e-02, &
         3.079230016738749e-02,   3.400213027432934e-02, &
         3.711627148341554e-02,   4.008382550403238e-02, &
         4.287284502017005e-02,   4.550291304992179e-02/
    data wgk(15),wgk(16),wgk(17),wgk(18),wgk(19),wgk(20),wgk(21), &
         wgk(22),wgk(23),wgk(24),wgk(25),wgk(26)/ &
         4.798253713883671e-02,   5.027767908071567e-02, &
         5.236288580640748e-02,   5.425112988854549e-02, &
         5.595081122041232e-02,   5.743711636156783e-02, &
         5.868968002239421e-02,   5.972034032417406e-02, &
         6.053945537604586e-02,   6.112850971705305e-02, &
         6.147118987142532e-02,   6.158081806783294e-02/
    data wg(1),wg(2),wg(3),wg(4),wg(5),wg(6),wg(7),wg(8),wg(9),wg(10), &
         wg(11),wg(12),wg(13)/ &
         1.139379850102629e-02,   2.635498661503214e-02, &
         4.093915670130631e-02,   5.490469597583519e-02, &
         6.803833381235692e-02,   8.014070033500102e-02, &
         9.102826198296365e-02,   1.005359490670506e-01, &
         1.085196244742637e-01,   1.148582591457116e-01, &
         1.194557635357848e-01,   1.222424429903100e-01, &
         1.231760537267155e-01/
    !
    centr = 5.0e-01*(a+b)
    hlgth = 5.0e-01*(b-a)
    dhlgth = abs(hlgth)
    !
    !  Compute the 51-point Kronrod approximation to the integral,
    !  and estimate the absolute error.
    !
    fc = f(centr)
    resg = wg(13)*fc
    resk = wgk(26)*fc
    resabs = abs(resk)

    do j = 1, 12
       jtw = j*2
       absc = hlgth*xgk(jtw)
       fval1 = f(centr-absc)
       fval2 = f(centr+absc)
       fv1(jtw) = fval1
       fv2(jtw) = fval2
       fsum = fval1+fval2
       resg = resg+wg(j)*fsum
       resk = resk+wgk(jtw)*fsum
       resabs = resabs+wgk(jtw)*(abs(fval1)+abs(fval2))
    end do

    do j = 1, 13
       jtwm1 = j*2-1
       absc = hlgth*xgk(jtwm1)
       fval1 = f(centr-absc)
       fval2 = f(centr+absc)
       fv1(jtwm1) = fval1
       fv2(jtwm1) = fval2
       fsum = fval1+fval2
       resk = resk+wgk(jtwm1)*fsum
       resabs = resabs+wgk(jtwm1)*(abs(fval1)+abs(fval2))
    end do

    reskh = resk*5.0e-01
    resasc = wgk(26)*abs(fc-reskh)

    do j = 1, 25
       resasc = resasc+wgk(j)*(abs(fv1(j)-reskh)+abs(fv2(j)-reskh))
    end do

    result = resk*hlgth
    resabs = resabs*dhlgth
    resasc = resasc*dhlgth
    abserr = abs((resk-resg)*hlgth)

    if ( resasc /= 0.0e+00.and.abserr /= 0.0e+00) then
       abserr = resasc*min ( 1.0e+00,(2.0e+02*abserr/resasc)**1.5e+00)
    end if

    if ( resabs > tiny ( resabs ) / (5.0e+01* epsilon ( resabs ) ) ) then
       abserr = max (( epsilon ( resabs ) *5.0e+01)*resabs,abserr)
    end if

    return
  end subroutine qk51
  subroutine qk61 ( f, a, b, result, abserr, resabs, resasc ) 

    !*****************************************************************************80
    !
    !! QK61 carries out a 61 point Gauss-Kronrod quadrature rule.
    !
    !  Discussion:
    !
    !    This routine approximates
    !      I = integral ( A <= X <= B ) F(X) dx
    !    with an error estimate, and
    !      J = integral ( A <= X <= B ) | F(X) | dx
    !
    !  Author:
    !
    !    Robert Piessens, Elise de Doncker-Kapenger, 
    !    Christian Ueberhuber, David Kahaner
    !
    !  Reference:
    !
    !    Robert Piessens, Elise de Doncker-Kapenger, 
    !    Christian Ueberhuber, David Kahaner,
    !    QUADPACK, a Subroutine Package for Automatic Integration,
    !    Springer Verlag, 1983
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
    !    Output, real(8) :: RESULT, the estimated value of the integral.
    !                    result is computed by applying the 61-point
    !                    Kronrod rule (resk) obtained by optimal addition of
    !                    abscissae to the 30-point Gauss rule (resg).
    !
    !    Output, real(8) :: ABSERR, an estimate of | I - RESULT |.
    !
    !    Output, real(8) :: RESABS, approximation to the integral of the absolute
    !    value of F.
    !
    !    Output, real(8) :: RESASC, approximation to the integral | F-I/(B-A) | 
    !    over [A,B].
    !
    !  Local Parameters:
    !
    !           centr  - mid point of the interval
    !           hlgth  - half-length of the interval
    !           absc   - abscissa
    !           fval*  - function value
    !           resg   - result of the 30-point Gauss rule
    !           resk   - result of the 61-point Kronrod rule
    !           reskh  - approximation to the mean value of f
    !                    over (a,b), i.e. to i/(b-a)
    !
    implicit none

    real(8) :: a
    real(8) :: absc
    real(8) :: abserr
    real(8) :: b
    real(8) :: centr
    real(8) :: dhlgth
    real(8), external :: f
    real(8) :: fc
    real(8) :: fsum
    real(8) :: fval1
    real(8) :: fval2
    real(8) :: fv1(30)
    real(8) :: fv2(30)
    real(8) :: hlgth
    integer j
    integer jtw
    integer jtwm1
    real(8) :: resabs
    real(8) :: resasc
    real(8) :: resg
    real(8) :: resk
    real(8) :: reskh
    real(8) :: result
    real(8) :: wg(15)
    real(8) :: wgk(31)
    real(8) :: xgk(31)
    !
    !           the abscissae and weights are given for the
    !           interval (-1,1). because of symmetry only the positive
    !           abscissae and their corresponding weights are given.
    !
    !           xgk   - abscissae of the 61-point Kronrod rule
    !                   xgk(2), xgk(4)  ... abscissae of the 30-point
    !                   Gauss rule
    !                   xgk(1), xgk(3)  ... optimally added abscissae
    !                   to the 30-point Gauss rule
    !
    !           wgk   - weights of the 61-point Kronrod rule
    !
    !           wg    - weigths of the 30-point Gauss rule
    !
    data xgk(1),xgk(2),xgk(3),xgk(4),xgk(5),xgk(6),xgk(7),xgk(8), &
         xgk(9),xgk(10)/ &
         9.994844100504906e-01,     9.968934840746495e-01, &
         9.916309968704046e-01,     9.836681232797472e-01, &
         9.731163225011263e-01,     9.600218649683075e-01, &
         9.443744447485600e-01,     9.262000474292743e-01, &
         9.055733076999078e-01,     8.825605357920527e-01/
    data xgk(11),xgk(12),xgk(13),xgk(14),xgk(15),xgk(16),xgk(17), &
         xgk(18),xgk(19),xgk(20)/ &
         8.572052335460611e-01,     8.295657623827684e-01, &
         7.997278358218391e-01,     7.677774321048262e-01, &
         7.337900624532268e-01,     6.978504947933158e-01, &
         6.600610641266270e-01,     6.205261829892429e-01, &
         5.793452358263617e-01,     5.366241481420199e-01/
    data xgk(21),xgk(22),xgk(23),xgk(24),xgk(25),xgk(26),xgk(27), &
         xgk(28),xgk(29),xgk(30),xgk(31)/ &
         4.924804678617786e-01,     4.470337695380892e-01, &
         4.004012548303944e-01,     3.527047255308781e-01, &
         3.040732022736251e-01,     2.546369261678898e-01, &
         2.045251166823099e-01,     1.538699136085835e-01, &
         1.028069379667370e-01,     5.147184255531770e-02, &
         0.0e+00                   /
    data wgk(1),wgk(2),wgk(3),wgk(4),wgk(5),wgk(6),wgk(7),wgk(8), &
         wgk(9),wgk(10)/ &
         1.389013698677008e-03,     3.890461127099884e-03, &
         6.630703915931292e-03,     9.273279659517763e-03, &
         1.182301525349634e-02,     1.436972950704580e-02, &
         1.692088918905327e-02,     1.941414119394238e-02, &
         2.182803582160919e-02,     2.419116207808060e-02/
    data wgk(11),wgk(12),wgk(13),wgk(14),wgk(15),wgk(16),wgk(17), &
         wgk(18),wgk(19),wgk(20)/ &
         2.650995488233310e-02,     2.875404876504129e-02, &
         3.090725756238776e-02,     3.298144705748373e-02, &
         3.497933802806002e-02,     3.688236465182123e-02, &
         3.867894562472759e-02,     4.037453895153596e-02, &
         4.196981021516425e-02,     4.345253970135607e-02/
    data wgk(21),wgk(22),wgk(23),wgk(24),wgk(25),wgk(26),wgk(27), &
         wgk(28),wgk(29),wgk(30),wgk(31)/ &
         4.481480013316266e-02,     4.605923827100699e-02, &
         4.718554656929915e-02,     4.818586175708713e-02, &
         4.905543455502978e-02,     4.979568342707421e-02, &
         5.040592140278235e-02,     5.088179589874961e-02, &
         5.122154784925877e-02,     5.142612853745903e-02, &
         5.149472942945157e-02/
    data wg(1),wg(2),wg(3),wg(4),wg(5),wg(6),wg(7),wg(8)/ &
         7.968192496166606e-03,     1.846646831109096e-02, &
         2.878470788332337e-02,     3.879919256962705e-02, &
         4.840267283059405e-02,     5.749315621761907e-02, &
         6.597422988218050e-02,     7.375597473770521e-02/
    data wg(9),wg(10),wg(11),wg(12),wg(13),wg(14),wg(15)/ &
         8.075589522942022e-02,     8.689978720108298e-02, &
         9.212252223778613e-02,     9.636873717464426e-02, &
         9.959342058679527e-02,     1.017623897484055e-01, &
         1.028526528935588e-01/

    centr = 5.0e-01*(b+a)
    hlgth = 5.0e-01*(b-a)
    dhlgth = abs(hlgth)
    !
    !  Compute the 61-point Kronrod approximation to the integral,
    !  and estimate the absolute error.
    !
    resg = 0.0e+00
    fc = f(centr)
    resk = wgk(31)*fc
    resabs = abs(resk)

    do j = 1, 15
       jtw = j*2
       absc = hlgth*xgk(jtw)
       fval1 = f(centr-absc)
       fval2 = f(centr+absc)
       fv1(jtw) = fval1
       fv2(jtw) = fval2
       fsum = fval1+fval2
       resg = resg+wg(j)*fsum
       resk = resk+wgk(jtw)*fsum
       resabs = resabs+wgk(jtw)*(abs(fval1)+abs(fval2))
    end do

    do j = 1, 15
       jtwm1 = j*2-1
       absc = hlgth*xgk(jtwm1)
       fval1 = f(centr-absc)
       fval2 = f(centr+absc)
       fv1(jtwm1) = fval1
       fv2(jtwm1) = fval2
       fsum = fval1+fval2
       resk = resk+wgk(jtwm1)*fsum
       resabs = resabs+wgk(jtwm1)*(abs(fval1)+abs(fval2))
    end do

    reskh = resk * 5.0e-01
    resasc = wgk(31)*abs(fc-reskh)

    do j = 1, 30
       resasc = resasc+wgk(j)*(abs(fv1(j)-reskh)+abs(fv2(j)-reskh))
    end do

    result = resk*hlgth
    resabs = resabs*dhlgth
    resasc = resasc*dhlgth
    abserr = abs((resk-resg)*hlgth)

    if ( resasc /= 0.0e+00 .and. abserr /= 0.0e+00) then
       abserr = resasc*min ( 1.0e+00,(2.0e+02*abserr/resasc)**1.5e+00)
    end if

    if ( resabs > tiny ( resabs ) / (5.0e+01* epsilon ( resabs ) )) then
       abserr = max ( ( epsilon ( resabs ) *5.0e+01)*resabs, abserr )
    end if

    return
  end subroutine qk61
  subroutine qmomo ( alfa, beta, ri, rj, rg, rh, integr )

    !*****************************************************************************80
    !
    !! QMOMO computes modified Chebyshev moments.
    !
    !  Discussion:
    !
    !    This routine computes modified Chebyshev moments.
    !    The K-th modified Chebyshev moment is defined as the
    !    integral over (-1,1) of W(X)*T(K,X), where T(K,X) is the
    !    Chebyshev polynomial of degree K.
    !
    !  Author:
    !
    !    Robert Piessens, Elise de Doncker-Kapenger, 
    !    Christian Ueberhuber, David Kahaner
    !
    !  Reference:
    !
    !    Robert Piessens, Elise de Doncker-Kapenger, 
    !    Christian Ueberhuber, David Kahaner,
    !    QUADPACK, a Subroutine Package for Automatic Integration,
    !    Springer Verlag, 1983
    !
    !  Parameters:
    !
    !    Input, real(8) :: ALFA, a parameter in the weight function w(x), ALFA > -1.
    !
    !    Input, real(8) :: BETA, a parameter in the weight function w(x), BETA > -1.
    !
    !           ri     - real
    !                    vector of dimension 25
    !                    ri(k) is the integral over (-1,1) of
    !                    (1+x)**alfa*t(k-1,x), k = 1, ..., 25.
    !
    !           rj     - real
    !                    vector of dimension 25
    !                    rj(k) is the integral over (-1,1) of
    !                    (1-x)**beta*t(k-1,x), k = 1, ..., 25.
    !
    !           rg     - real
    !                    vector of dimension 25
    !                    rg(k) is the integral over (-1,1) of
    !                    (1+x)**alfa*log((1+x)/2)*t(k-1,x), k = 1, ...,25.
    !
    !           rh     - real
    !                    vector of dimension 25
    !                    rh(k) is the integral over (-1,1) of
    !                    (1-x)**beta*log((1-x)/2)*t(k-1,x), k = 1, ..., 25.
    !
    !           integr - integer
    !                    input parameter indicating the modified moments
    !                    to be computed
    !                    integr = 1 compute ri, rj
    !                           = 2 compute ri, rj, rg
    !                           = 3 compute ri, rj, rh
    !                           = 4 compute ri, rj, rg, rh
    !
    implicit none

    real(8) :: alfa
    real(8) :: alfp1
    real(8) :: alfp2
    real(8) :: an
    real(8) :: anm1
    real(8) :: beta
    real(8) :: betp1
    real(8) :: betp2
    integer i
    integer im1
    integer integr
    real(8) :: ralf
    real(8) :: rbet
    real(8) :: rg(25)
    real(8) :: rh(25)
    real(8) :: ri(25)
    real(8) :: rj(25)
    !
    alfp1 = alfa+1.0e+00
    betp1 = beta+1.0e+00
    alfp2 = alfa+2.0e+00
    betp2 = beta+2.0e+00
    ralf = 2.0e+00**alfp1
    rbet = 2.0e+00**betp1
    !
    !  Compute RI, RJ using a forward recurrence relation.
    !
    ri(1) = ralf/alfp1
    rj(1) = rbet/betp1
    ri(2) = ri(1)*alfa/alfp2
    rj(2) = rj(1)*beta/betp2
    an = 2.0e+00
    anm1 = 1.0e+00

    do i = 3, 25
       ri(i) = -(ralf+an*(an-alfp2)*ri(i-1))/(anm1*(an+alfp1))
       rj(i) = -(rbet+an*(an-betp2)*rj(i-1))/(anm1*(an+betp1))
       anm1 = an
       an = an+1.0e+00
    end do

    if ( integr == 1 ) go to 70
    if ( integr == 3 ) go to 40
    !
    !  Compute RG using a forward recurrence relation.
    !
    rg(1) = -ri(1)/alfp1
    rg(2) = -(ralf+ralf)/(alfp2*alfp2)-rg(1)
    an = 2.0e+00
    anm1 = 1.0e+00
    im1 = 2

    do i = 3, 25
       rg(i) = -(an*(an-alfp2)*rg(im1)-an*ri(im1)+anm1*ri(i))/ &
            (anm1*(an+alfp1))
       anm1 = an
       an = an+1.0e+00
       im1 = i
    end do

    if ( integr == 2 ) go to 70
    !
    !  Compute RH using a forward recurrence relation.
    !
40  continue

    rh(1) = -rj(1) / betp1
    rh(2) = -(rbet+rbet)/(betp2*betp2)-rh(1)
    an = 2.0e+00
    anm1 = 1.0e+00
    im1 = 2

    do i = 3, 25
       rh(i) = -(an*(an-betp2)*rh(im1)-an*rj(im1)+ &
            anm1*rj(i))/(anm1*(an+betp1))
       anm1 = an
       an = an+1.0e+00
       im1 = i
    end do

    do i = 2, 25, 2
       rh(i) = -rh(i)
    end do

70  continue

    do i = 2, 25, 2
       rj(i) = -rj(i)
    end do

    !  90 continue

    return
  end subroutine qmomo
  subroutine qsort ( limit, last, maxerr, ermax, elist, iord, nrmax )

    !*****************************************************************************80
    !
    !! QSORT maintains the order of a list of local error estimates.
    !
    !  Discussion:
    !
    !    This routine maintains the descending ordering in the list of the 
    !    local error estimates resulting from the interval subdivision process. 
    !    At each call two error estimates are inserted using the sequential 
    !    search top-down for the largest error estimate and bottom-up for the
    !    smallest error estimate.
    !
    !  Author:
    !
    !    Robert Piessens, Elise de Doncker-Kapenger, 
    !    Christian Ueberhuber, David Kahaner
    !
    !  Reference:
    !
    !    Robert Piessens, Elise de Doncker-Kapenger, 
    !    Christian Ueberhuber, David Kahaner,
    !    QUADPACK, a Subroutine Package for Automatic Integration,
    !    Springer Verlag, 1983
    !
    !  Parameters:
    !
    !    Input, integer LIMIT, the maximum number of error estimates the list can
    !    contain.
    !
    !    Input, integer LAST, the current number of error estimates.
    !
    !    Input/output, integer MAXERR, the index in the list of the NRMAX-th 
    !    largest error.
    !
    !    Output, real(8) :: ERMAX, the NRMAX-th largest error = ELIST(MAXERR).
    !
    !    Input, real(8) :: ELIST(LIMIT), contains the error estimates.
    !
    !    Input/output, integer IORD(LAST).  The first K elements contain 
    !    pointers to the error estimates such that ELIST(IORD(1)) through
    !    ELIST(IORD(K)) form a decreasing sequence, with
    !      K = LAST 
    !    if 
    !      LAST <= (LIMIT/2+2), 
    !    and otherwise
    !      K = LIMIT+1-LAST.
    !
    !    Input/output, integer NRMAX.
    !
    implicit none

    integer last

    real(8) :: elist(last)
    real(8) :: ermax
    real(8) :: errmax
    real(8) :: errmin
    integer i
    integer ibeg
    integer iord(last)
    integer isucc
    integer j
    integer jbnd
    integer jupbn
    integer k
    integer limit
    integer maxerr
    integer nrmax
    !
    !  Check whether the list contains more than two error estimates.
    !
    if ( last <= 2 ) then
       iord(1) = 1
       iord(2) = 2
       go to 90
    end if
    !
    !  This part of the routine is only executed if, due to a
    !  difficult integrand, subdivision increased the error
    !  estimate. in the normal case the insert procedure should
    !  start after the nrmax-th largest error estimate.
    !
    errmax = elist(maxerr)

    do i = 1, nrmax-1

       isucc = iord(nrmax-1)

       if ( errmax <= elist(isucc) ) then
          exit
       end if

       iord(nrmax) = isucc
       nrmax = nrmax-1

    end do
    !
    !  Compute the number of elements in the list to be maintained
    !  in descending order.  This number depends on the number of
    !  subdivisions still allowed.
    !
    jupbn = last

    if ( (limit/2+2) < last ) then
       jupbn = limit+3-last
    end if

    errmin = elist(last)
    !
    !  Insert errmax by traversing the list top-down, starting
    !  comparison from the element elist(iord(nrmax+1)).
    !
    jbnd = jupbn-1
    ibeg = nrmax+1

    do i = ibeg, jbnd
       isucc = iord(i)
       if ( elist(isucc) <= errmax ) then
          go to 60
       end if
       iord(i-1) = isucc
    end do

    iord(jbnd) = maxerr
    iord(jupbn) = last
    go to 90
    !
    !  Insert errmin by traversing the list bottom-up.
    !
60  continue

    iord(i-1) = maxerr
    k = jbnd

    do j = i, jbnd
       isucc = iord(k)
       if ( errmin < elist(isucc) ) then
          go to 80
       end if
       iord(k+1) = isucc
       k = k-1
    end do

    iord(i) = last
    go to 90

80  continue

    iord(k+1) = last
    !
    !  Set maxerr and ermax.
    !
90  continue

    maxerr = iord(nrmax)
    ermax = elist(maxerr)

    return
  end subroutine qsort
  function qwgtc ( x, c, p2, p3, p4, kp )

    !*****************************************************************************80
    !
    !! QWGTC defines the weight function used by QC25C.
    !
    !  Discussion:
    !
    !    The weight function has the form 1 / ( X - C ).
    !
    !  Author:
    !
    !    Robert Piessens, Elise de Doncker-Kapenger, 
    !    Christian Ueberhuber, David Kahaner
    !
    !  Reference:
    !
    !    Robert Piessens, Elise de Doncker-Kapenger, 
    !    Christian Ueberhuber, David Kahaner,
    !    QUADPACK, a Subroutine Package for Automatic Integration,
    !    Springer Verlag, 1983
    !
    !  Parameters:
    !
    !    Input, real(8) :: X, the point at which the weight function is evaluated.
    !
    !    Input, real(8) :: C, the location of the singularity.
    !
    !    Input, real(8) :: P2, P3, P4, parameters that are not used.
    !
    !    Input, integer KP, a parameter that is not used.
    !
    !    Output, real(8) :: QWGTC, the value of the weight function at X.
    !
    implicit none

    real(8) :: c
    integer kp
    real(8) :: p2
    real(8) :: p3
    real(8) :: p4
    real(8) :: qwgtc
    real(8) :: x

    qwgtc = 1.0E+00 / ( x - c )

    return
  end function qwgtc
  function qwgto ( x, omega, p2, p3, p4, integr )

    !*****************************************************************************80
    !
    !! QWGTO defines the weight functions used by QC25O.
    !
    !  Author:
    !
    !    Robert Piessens, Elise de Doncker-Kapenger, 
    !    Christian Ueberhuber, David Kahaner
    !
    !  Reference:
    !
    !    Robert Piessens, Elise de Doncker-Kapenger, 
    !    Christian Ueberhuber, David Kahaner,
    !    QUADPACK, a Subroutine Package for Automatic Integration,
    !    Springer Verlag, 1983
    !
    !  Parameters:
    !
    !    Input, real(8) :: X, the point at which the weight function is evaluated.
    !
    !    Input, real(8) :: OMEGA, the factor multiplying X.
    !
    !    Input, real(8) :: P2, P3, P4, parameters that are not used.
    !
    !    Input, integer INTEGR, specifies which weight function is used:
    !    1. W(X) = cos ( OMEGA * X )
    !    2, W(X) = sin ( OMEGA * X )
    !
    !    Output, real(8) :: QWGTO, the value of the weight function at X.
    !
    implicit none

    integer integr
    real(8) :: omega
    real(8) :: p2
    real(8) :: p3
    real(8) :: p4
    real(8) :: qwgto
    real(8) :: x

    if ( integr == 1 ) then
       qwgto = cos ( omega * x )
    else if ( integr == 2 ) then
       qwgto = sin ( omega * x )
    end if

    return
  end function qwgto
  function qwgts ( x, a, b, alfa, beta, integr )

    !*****************************************************************************80
    !
    !! QWGTS defines the weight functions used by QC25S.
    !
    !  Author:
    !
    !    Robert Piessens, Elise de Doncker-Kapenger, 
    !    Christian Ueberhuber, David Kahaner
    !
    !  Reference:
    !
    !    Robert Piessens, Elise de Doncker-Kapenger, 
    !    Christian Ueberhuber, David Kahaner,
    !    QUADPACK, a Subroutine Package for Automatic Integration,
    !    Springer Verlag, 1983
    !
    !  Parameters:
    !
    !    Input, real(8) :: X, the point at which the weight function is evaluated.
    !
    !    Input, real(8) :: A, B, the endpoints of the integration interval.
    !
    !    Input, real(8) :: ALFA, BETA, exponents that occur in the weight function.
    !
    !    Input, integer INTEGR, specifies which weight function is used:
    !    1. W(X) = (X-A)**ALFA * (B-X)**BETA
    !    2, W(X) = (X-A)**ALFA * (B-X)**BETA * log (X-A)
    !    3, W(X) = (X-A)**ALFA * (B-X)**BETA * log (B-X)
    !    4, W(X) = (X-A)**ALFA * (B-X)**BETA * log (X-A) * log(B-X)
    !
    !    Output, real(8) :: QWGTS, the value of the weight function at X.
    !
    implicit none

    real(8) :: a
    real(8) :: alfa
    real(8) :: b
    real(8) :: beta
    integer integr
    real(8) :: qwgts
    real(8) :: x

    if ( integr == 1 ) then
       qwgts = ( x - a )**alfa * ( b - x )**beta
    else if ( integr == 2 ) then
       qwgts = ( x - a )**alfa * ( b - x )**beta * log ( x - a )
    else if ( integr == 3 ) then
       qwgts = ( x - a )**alfa * ( b - x )**beta * log ( b - x )
    else if ( integr == 4 ) then
       qwgts = ( x - a )**alfa * ( b - x )**beta * log ( x - a ) * log ( b - x )
    end if

    return
  end function qwgts
