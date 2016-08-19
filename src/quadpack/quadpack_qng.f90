subroutine qng ( f, a, b, epsabs, epsrel, result, abserr, neval, ier )
  !*****************************************************************************80
  !! QNG estimates an integral, using non-adaptive integration.
  !  Discussion:
  !
  !    The routine calculates an approximation RESULT to a definite integral   
  !      I = integral of F over (A,B),
  !    hopefully satisfying
  !      || I - RESULT || <= max ( EPSABS, EPSREL * ||I|| ).
  !
  !    The routine is a simple non-adaptive automatic integrator, based on
  !    a sequence of rules with increasing degree of algebraic
  !    precision (Patterson, 1968).
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
  !    RESULT is obtained by applying the 21-point Gauss-Kronrod rule (RES21)
  !    obtained  by optimal addition of abscissae to the 10-point Gauss rule
  !    (RES10), or by applying the 43-point rule (RES43) obtained by optimal
  !    addition of abscissae to the 21-point Gauss-Kronrod rule, or by 
  !    applying the 87-point rule (RES87) obtained by optimal addition of
  !    abscissae to the 43-point rule.
  !
  !    Output, real(8) :: ABSERR, an estimate of || I - RESULT ||.
  !
  !    Output, integer NEVAL, the number of times the integral was evaluated.
  !
  !           ier    - ier = 0 normal and reliable termination of the
  !                            routine. it is assumed that the requested
  !                            accuracy has been achieved.
  !                    ier > 0 abnormal termination of the routine. it is
  !                            assumed that the requested accuracy has
  !                            not been achieved.
  !                    ier = 1 the maximum number of steps has been
  !                            executed. the integral is probably too
  !                            difficult to be calculated by qng.
  !                        = 6 the input is invalid, because
  !                            epsabs < 0 and epsrel < 0,
  !                            result, abserr and neval are set to zero.
  !
  !  Local Parameters:
  !
  !           centr  - mid point of the integration interval
  !           hlgth  - half-length of the integration interval
  !           fcentr - function value at mid point
  !           absc   - abscissa
  !           fval   - function value
  !           savfun - array of function values which have already
  !                    been computed
  !           res10  - 10-point Gauss result
  !           res21  - 21-point Kronrod result
  !           res43  - 43-point result
  !           res87  - 87-point result
  !           resabs - approximation to the integral of abs(f)
  !           resasc - approximation to the integral of abs(f-i/(b-a))
  !
  implicit none
  real(8) :: a
  real(8) :: absc
  real(8) :: abserr
  real(8) :: b
  real(8) :: centr
  real(8) :: dhlgth
  real(8) :: epsabs
  real(8) :: epsrel
  real, external :: f
  real(8) :: fcentr
  real(8) :: fval
  real(8) :: fval1
  real(8) :: fval2
  real(8) :: fv1(5)
  real(8) :: fv2(5)
  real(8) :: fv3(5)
  real(8) :: fv4(5)
  real(8) :: hlgth
  integer ier
  integer ipx
  integer k
  integer l
  integer neval
  real(8) :: result
  real(8) :: res10
  real(8) :: res21
  real(8) :: res43
  real(8) :: res87
  real(8) :: resabs
  real(8) :: resasc
  real(8) :: reskh
  real(8) :: savfun(21)
  real(8) :: w10(5)
  real(8) :: w21a(5)
  real(8) :: w21b(6)
  real(8) :: w43a(10)
  real(8) :: w43b(12)
  real(8) :: w87a(21)
  real(8) :: w87b(23)
  real(8) :: x1(5)
  real(8) :: x2(5)
  real(8) :: x3(11)
  real(8) :: x4(22)
  !
  !           the following data statements contain the abscissae
  !           and weights of the integration rules used.
  !
  !           x1      abscissae common to the 10-, 21-, 43- and 87-point
  !                   rule
  !           x2      abscissae common to the 21-, 43- and 87-point rule
  !           x3      abscissae common to the 43- and 87-point rule
  !           x4      abscissae of the 87-point rule
  !           w10     weights of the 10-point formula
  !           w21a    weights of the 21-point formula for abscissae x1
  !           w21b    weights of the 21-point formula for abscissae x2
  !           w43a    weights of the 43-point formula for absissae x1, x3
  !           w43b    weights of the 43-point formula for abscissae x3
  !           w87a    weights of the 87-point formula for abscissae x1,
  !                   x2 and x3
  !           w87b    weights of the 87-point formula for abscissae x4
  !
  data x1(1),x1(2),x1(3),x1(4),x1(5)/ &
       9.739065285171717e-01,     8.650633666889845e-01, &
       6.794095682990244e-01,     4.333953941292472e-01, &
       1.488743389816312e-01/
  data x2(1),x2(2),x2(3),x2(4),x2(5)/ &
       9.956571630258081e-01,     9.301574913557082e-01, &
       7.808177265864169e-01,     5.627571346686047e-01, &
       2.943928627014602e-01/
  data x3(1),x3(2),x3(3),x3(4),x3(5),x3(6),x3(7),x3(8),x3(9),x3(10), &
       x3(11)/ &
       9.993333609019321e-01,     9.874334029080889e-01, &
       9.548079348142663e-01,     9.001486957483283e-01, &
       8.251983149831142e-01,     7.321483889893050e-01, &
       6.228479705377252e-01,     4.994795740710565e-01, &
       3.649016613465808e-01,     2.222549197766013e-01, &
       7.465061746138332e-02/
  data x4(1),x4(2),x4(3),x4(4),x4(5),x4(6),x4(7),x4(8),x4(9),x4(10), &
       x4(11),x4(12),x4(13),x4(14),x4(15),x4(16),x4(17),x4(18),x4(19), &
       x4(20),x4(21),x4(22)/         9.999029772627292e-01, &
       9.979898959866787e-01,     9.921754978606872e-01, &
       9.813581635727128e-01,     9.650576238583846e-01, &
       9.431676131336706e-01,     9.158064146855072e-01, &
       8.832216577713165e-01,     8.457107484624157e-01, &
       8.035576580352310e-01,     7.570057306854956e-01, &
       7.062732097873218e-01,     6.515894665011779e-01, &
       5.932233740579611e-01,     5.314936059708319e-01, &
       4.667636230420228e-01,     3.994248478592188e-01, &
       3.298748771061883e-01,     2.585035592021616e-01, &
       1.856953965683467e-01,     1.118422131799075e-01, &
       3.735212339461987e-02/
  data w10(1),w10(2),w10(3),w10(4),w10(5)/ &
       6.667134430868814e-02,     1.494513491505806e-01, &
       2.190863625159820e-01,     2.692667193099964e-01, &
       2.955242247147529e-01/
  data w21a(1),w21a(2),w21a(3),w21a(4),w21a(5)/ &
       3.255816230796473e-02,     7.503967481091995e-02, &
       1.093871588022976e-01,     1.347092173114733e-01, &
       1.477391049013385e-01/
  data w21b(1),w21b(2),w21b(3),w21b(4),w21b(5),w21b(6)/ &
       1.169463886737187e-02,     5.475589657435200e-02, &
       9.312545458369761e-02,     1.234919762620659e-01, &
       1.427759385770601e-01,     1.494455540029169e-01/
  data w43a(1),w43a(2),w43a(3),w43a(4),w43a(5),w43a(6),w43a(7), &
       w43a(8),w43a(9),w43a(10)/     1.629673428966656e-02, &
       3.752287612086950e-02,     5.469490205825544e-02, &
       6.735541460947809e-02,     7.387019963239395e-02, &
       5.768556059769796e-03,     2.737189059324884e-02, &
       4.656082691042883e-02,     6.174499520144256e-02, &
       7.138726726869340e-02/
  data w43b(1),w43b(2),w43b(3),w43b(4),w43b(5),w43b(6),w43b(7), &
       w43b(8),w43b(9),w43b(10),w43b(11),w43b(12)/ &
       1.844477640212414e-03,     1.079868958589165e-02, &
       2.189536386779543e-02,     3.259746397534569e-02, &
       4.216313793519181e-02,     5.074193960018458e-02, &
       5.837939554261925e-02,     6.474640495144589e-02, &
       6.956619791235648e-02,     7.282444147183321e-02, &
       7.450775101417512e-02,     7.472214751740301e-02/
  data w87a(1),w87a(2),w87a(3),w87a(4),w87a(5),w87a(6),w87a(7), &
       w87a(8),w87a(9),w87a(10),w87a(11),w87a(12),w87a(13),w87a(14), &
       w87a(15),w87a(16),w87a(17),w87a(18),w87a(19),w87a(20),w87a(21)/ &
       8.148377384149173e-03,     1.876143820156282e-02, &
       2.734745105005229e-02,     3.367770731163793e-02, &
       3.693509982042791e-02,     2.884872430211531e-03, &
       1.368594602271270e-02,     2.328041350288831e-02, &
       3.087249761171336e-02,     3.569363363941877e-02, &
       9.152833452022414e-04,     5.399280219300471e-03, &
       1.094767960111893e-02,     1.629873169678734e-02, &
       2.108156888920384e-02,     2.537096976925383e-02, &
       2.918969775647575e-02,     3.237320246720279e-02, &
       3.478309895036514e-02,     3.641222073135179e-02, &
       3.725387550304771e-02/
  data w87b(1),w87b(2),w87b(3),w87b(4),w87b(5),w87b(6),w87b(7), &
       w87b(8),w87b(9),w87b(10),w87b(11),w87b(12),w87b(13),w87b(14), &
       w87b(15),w87b(16),w87b(17),w87b(18),w87b(19),w87b(20),w87b(21), &
       w87b(22),w87b(23)/            2.741455637620724e-04, &
       1.807124155057943e-03,     4.096869282759165e-03, &
       6.758290051847379e-03,     9.549957672201647e-03, &
       1.232944765224485e-02,     1.501044734638895e-02, &
       1.754896798624319e-02,     1.993803778644089e-02, &
       2.219493596101229e-02,     2.433914712600081e-02, &
       2.637450541483921e-02,     2.828691078877120e-02, &
       3.005258112809270e-02,     3.164675137143993e-02, &
       3.305041341997850e-02,     3.425509970422606e-02, &
       3.526241266015668e-02,     3.607698962288870e-02, &
       3.669860449845609e-02,     3.712054926983258e-02, &
       3.733422875193504e-02,     3.736107376267902e-02/
  !
  !  Test on validity of parameters.
  !
  result = 0.0e+00
  abserr = 0.0e+00
  neval = 0
  if ( epsabs < 0.0e+00 .and. epsrel < 0.0e+00 ) then
     ier = 6
     return
  end if
  hlgth = 5.0e-01 * ( b - a )
  dhlgth = abs ( hlgth )
  centr = 5.0e-01 * ( b + a )
  fcentr = f(centr)
  neval = 21
  ier = 1
  !
  !  Compute the integral using the 10- and 21-point formula.
  !
  do l = 1, 3
     if ( l == 1 ) then
        res10 = 0.0e+00
        res21 = w21b(6) * fcentr
        resabs = w21b(6) * abs(fcentr)
        do k = 1, 5
           absc = hlgth * x1(k)
           fval1 = f(centr+absc)
           fval2 = f(centr-absc)
           fval = fval1 + fval2
           res10 = res10 + w10(k)*fval
           res21 = res21 + w21a(k)*fval
           resabs = resabs + w21a(k)*(abs(fval1)+abs(fval2))
           savfun(k) = fval
           fv1(k) = fval1
           fv2(k) = fval2
        end do
        ipx = 5
        do k = 1, 5
           ipx = ipx + 1
           absc = hlgth * x2(k)
           fval1 = f(centr+absc)
           fval2 = f(centr-absc)
           fval = fval1 + fval2
           res21 = res21 + w21b(k) * fval
           resabs = resabs + w21b(k) * ( abs ( fval1 ) + abs ( fval2 ) )
           savfun(ipx) = fval
           fv3(k) = fval1
           fv4(k) = fval2
        end do
        !
        !  Test for convergence.
        !
        result = res21 * hlgth
        resabs = resabs * dhlgth
        reskh = 5.0e-01 * res21
        resasc = w21b(6) * abs ( fcentr - reskh )
        do k = 1, 5
           resasc = resasc+w21a(k)*(abs(fv1(k)-reskh)+abs(fv2(k)-reskh)) &
                +w21b(k)*(abs(fv3(k)-reskh)+abs(fv4(k)-reskh))
        end do
        abserr = abs ( ( res21 - res10 ) * hlgth )
        resasc = resasc * dhlgth
        !
        !  Compute the integral using the 43-point formula.
        !
     else if ( l == 2 ) then
        res43 = w43b(12)*fcentr
        neval = 43
        do k = 1, 10
           res43 = res43 + savfun(k) * w43a(k)
        end do
        do k = 1, 11
           ipx = ipx + 1
           absc = hlgth * x3(k)
           fval = f(absc+centr) + f(centr-absc)
           res43 = res43 + fval * w43b(k)
           savfun(ipx) = fval
        end do
        !
        !  Test for convergence.
        !
        result = res43 * hlgth
        abserr = abs((res43-res21)*hlgth)
        !
        !  Compute the integral using the 87-point formula.
        !
     else if ( l == 3 ) then
        res87 = w87b(23) * fcentr
        neval = 87
        do k = 1, 21
           res87 = res87 + savfun(k) * w87a(k)
        end do
        do k = 1, 22
           absc = hlgth * x4(k)
           res87 = res87 + w87b(k) * ( f(absc+centr) + f(centr-absc) )
        end do
        result = res87 * hlgth
        abserr = abs ( ( res87 - res43) * hlgth )
     end if
     if ( resasc /= 0.0e+00.and.abserr /= 0.0e+00 ) then
        abserr = resasc * min ( 1.0e+00,(2.0e+02*abserr/resasc)**1.5e+00)
     end if
     if ( resabs > tiny ( resabs ) / ( 5.0e+01 * epsilon ( resabs ) ) ) then
        abserr = max (( epsilon ( resabs ) *5.0e+01) * resabs, abserr )
     end if
     if ( abserr <= max ( epsabs, epsrel*abs(result))) then
        ier = 0
     end if
     if ( ier == 0 ) then
        exit
     end if
  end do
  return
end subroutine qng
