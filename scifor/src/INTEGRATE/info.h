    !*************************************************************************
    !
    !! AAAA is a dummy subroutine with QUADPACK documentation in its comments.
    !
    ! 1. introduction
    !
    !    quadpack is a fortran subroutine package for the numerical
    !    computation of definite one-dimensional integrals. it originated
    !    from a joint project of r. piessens and e. de doncker (appl.
    !    math. and progr. div.- k.u.leuven, belgium), c. ueberhuber (inst.
    !    fuer math.- techn.u.wien, austria), and d. kahaner (nation. bur.
    !    of standards- washington d.c., u.s.a.).
    !
    ! 2. survey
    !
    !    - qags : is an integrator based on globally adaptive interval
    !             subdivision in connection with extrapolation (de doncker,
    !             1978) by the epsilon algorithm (wynn, 1956).
    !
    !    - qagp : serves the same purposes as qags, but also allows
    !             for eventual user-supplied information, i.e. the
    !             abscissae of internal singularities, discontinuities
    !             and other difficulties of the integrand function.
    !             the algorithm is a modification of that in qags.
    !
    !    - qagi : handles integration over infinite intervals. the
    !             infinite range is mapped onto a finite interval and
    !             then the same strategy as in qags is applied.
    !
    !    - qawo : is a routine for the integration of cos(omega*x)*f(x)
    !             or sin(omega*x)*f(x) over a finite interval (a,b).
    !             omega is is specified by the user
    !             the rule evaluation component is based on the
    !             modified clenshaw-curtis technique.
    !             an adaptive subdivision scheme is used connected with
    !             an extrapolation procedure, which is a modification
    !             of that in qags and provides the possibility to deal
    !             even with singularities in f.
    !
    !    - qawf : calculates the fourier cosine or fourier sine
    !             transform of f(x), for user-supplied interval (a,
    !             infinity), omega, and f. the procedure of qawo is
    !             used on successive finite intervals, and convergence
    !             acceleration by means of the epsilon algorithm (wynn,
    !             1956) is applied to the series of the integral
    !             contributions.
    !
    !    - qaws : integrates w(x)*f(x) over (a,b) with a < b finite,
    !             and   w(x) = ((x-a)**alfa)*((b-x)**beta)*v(x)
    !             where v(x) = 1 or log(x-a) or log(b-x)
    !                            or log(x-a)*log(b-x)
    !             and   -1 < alfa, -1 < beta.
    !             the user specifies a, b, alfa, beta and the type of
    !             the function v.
    !             a globally adaptive subdivision strategy is applied,
    !             with modified clenshaw-curtis integration on the
    !             subintervals which contain a or b.
    !
    !    - qawc : computes the cauchy principal value of f(x)/(x-c)
    !             over a finite interval (a,b) and for
    !             user-determined c.
    !             the strategy is globally adaptive, and modified
    !             clenshaw-curtis integration is used on the subranges
    !             which contain the point x = c.
    !
    !  each of the routines above also has a "more detailed" version
    !    with a name ending in e, as qage.  these provide more
    !    information and control than the easier versions.
    !
    !
    !   the preceeding routines are all automatic.  that is, the user
    !      inputs his problem and an error tolerance.  the routine
    !      attempts to perform the integration to within the requested
    !      absolute or relative error.
    !   there are, in addition, a number of non-automatic integrators.
    !      these are most useful when the problem is such that the
    !      user knows that a fixed rule will provide the accuracy
    !      required.  typically they return an error estimate but make
    !      no attempt to satisfy any particular input error request.
    !
    !      qk15
    !      qk21
    !      qk31
    !      qk41
    !      qk51
    !      qk61
    !           estimate the integral on [a,b] using 15, 21,..., 61
    !           point rule and return an error estimate.
    !      qk15i 15 point rule for (semi)infinite interval.
    !      qk15w 15 point rule for special singular weight functions.
    !      qc25c 25 point rule for cauchy principal values
    !      qc25o 25 point rule for sin/cos integrand.
    !      qmomo integrates k-th degree chebychev polynomial times
    !            function with various explicit singularities.
    !
    ! 3. guidelines for the use of quadpack
    !
    !    here it is not our purpose to investigate the question when
    !    automatic quadrature should be used. we shall rather attempt
    !    to help the user who already made the decision to use quadpack,
    !    with selecting an appropriate routine or a combination of
    !    several routines for handling his problem.
    !
    !    for both quadrature over finite and over infinite intervals,
    !    one of the first questions to be answered by the user is
    !    related to the amount of computer time he wants to spend,
    !    versus his -own- time which would be needed, for example, for
    !    manual subdivision of the interval or other analytic
    !    manipulations.
    !
    !    (1) the user may not care about computer time, or not be
    !        willing to do any analysis of the problem. especially when
    !        only one or a few integrals must be calculated, this attitude
    !        can be perfectly reasonable. in this case it is clear that
    !        either the most sophisticated of the routines for finite
    !        intervals, qags, must be used, or its analogue for infinite
    !        intervals, qagi. these routines are able to cope with
    !        rather difficult, even with improper integrals.
    !        this way of proceeding may be expensive. but the integrator
    !        is supposed to give you an answer in return, with additional
    !        information in the case of a failure, through its error
    !        estimate and flag. yet it must be stressed that the programs
    !        cannot be totally reliable.
    !
    !    (2) the user may want to examine the integrand function.
    !        if bad local difficulties occur, such as a discontinuity, a
    !        singularity, derivative singularity or high peak at one or
    !        more points within the interval, the first advice is to
    !        split up the interval at these points. the integrand must
    !        then be examinated over each of the subintervals separately,
    !        so that a suitable integrator can be selected for each of
    !        them. if this yields problems involving relative accuracies
    !        to be imposed on -finite- subintervals, one can make use of
    !        qagp, which must be provided with the positions of the local
    !        difficulties. however, if strong singularities are present
    !        and a high accuracy is requested, application of qags on the
    !        subintervals may yield a better result.
    !
    !        for quadrature over finite intervals we thus dispose of qags
    !        and
    !        - qng for well-behaved integrands,
    !        - qag for functions with an oscillating behavior of a non
    !          specific type,
    !        - qawo for functions, eventually singular, containing a
    !          factor cos(omega*x) or sin(omega*x) where omega is known,
    !        - qaws for integrands with algebraico-logarithmic end point
    !          singularities of known type,
    !        - qawc for cauchy principal values.
    !
    !        remark
    !
    !        on return, the work arrays in the argument lists of the
    !        adaptive integrators contain information about the interval
    !        subdivision process and hence about the integrand behavior:
    !        the end points of the subintervals, the local integral
    !        contributions and error estimates, and eventually other
    !        characteristics. for this reason, and because of its simple
    !        globally adaptive nature, the routine qag in particular is
    !        well-suited for integrand examination. difficult spots can
    !        be located by investigating the error estimates on the
    !        subintervals.
    !
    !        for infinite intervals we provide only one general-purpose
    !        routine, qagi. it is based on the qags algorithm applied
    !        after a transformation of the original interval into (0,1).
    !        yet it may eventuate that another type of transformation is
    !        more appropriate, or one might prefer to break up the
    !        original interval and use qagi only on the infinite part
    !        and so on. these kinds of actions suggest a combined use of
    !        different quadpack integrators. note that, when the only
    !        difficulty is an integrand singularity at the finite
    !        integration limit, it will in general not be necessary to
    !        break up the interval, as qagi deals with several types of
    !        singularity at the boundary point of the integration range.
    !        it also handles slowly convergent improper integrals, on
    !        the condition that the integrand does not oscillate over
    !        the entire infinite interval. if it does we would advise
    !        to sum succeeding positive and negative contributions to
    !        the integral -e.g. integrate between the zeros- with one
    !        or more of the finite-range integrators, and apply
    !        convergence acceleration eventually by means of quadpack
    !        subroutine qelg which implements the epsilon algorithm.
    !        such quadrature problems include the fourier transform as
    !        a special case. yet for the latter we have an automatic
    !        integrator available, qawf.
    !
