subroutine cubspl ( tau, c, n, ibcbeg, ibcend )
  !***********************************************************************
  !! CUBSPL defines an interpolatory cubic spline.
  !  Discussion:
  !    A tridiagonal linear system for the unknown slopes S(I) of
  !    F at TAU(I), I=1,..., N, is generated and then solved by Gauss
  !    elimination, with S(I) ending up in C(2,I), for all I.

  !  Parameters:
  !
  !    Input, real ( kind = 8 ) TAU(N), the abscissas or X values of
  !    the data points.  The entries of TAU are assumed to be
  !    strictly increasing.
  !
  !    Input, integer ( kind = 4 ) N, the number of data points.  N is
  !    assumed to be at least 2.
  !
  !    Input/output, real ( kind = 8 ) C(4,N).
  !    On input, if IBCBEG or IBCBEG is 1 or 2, then C(2,1)
  !    or C(2,N) should have been set to the desired derivative
  !    values, as described further under IBCBEG and IBCEND.
  !    On output, C contains the polynomial coefficients of
  !    the cubic interpolating spline with interior knots
  !    TAU(2) through TAU(N-1).
  !    In the interval interval (TAU(I), TAU(I+1)), the spline
  !    F is given by
  !      F(X) = 
  !        C(1,I) + 
  !        C(2,I) * H +
  !        C(3,I) * H**2 / 2 + 
  !        C(4,I) * H**3 / 6.
  !    where H=X-TAU(I).  The routine PPVALU may be used to
  !    evaluate F or its derivatives from TAU, C, L=N-1,
  !    and K=4.
  !
  !    Input, integer ( kind = 4 ) IBCBEG, IBCEND, boundary condition indicators
  !    IBCBEG = 0 means no boundary condition at TAU(1) is given.
  !    In this case, the "not-a-knot condition" is used.  That
  !    is, the jump in the third derivative across TAU(2) is
  !    forced to zero.  Thus the first and the second cubic
  !    polynomial pieces are made to coincide.
  !    IBCBEG = 1 means the slope at TAU(1) is to equal the
  !    input value C(2,1).
  !    IBCBEG = 2 means the second derivative at TAU(1) is
  !    to equal C(2,1).
  !    IBCEND = 0, 1, or 2 has analogous meaning concerning the
  !    boundary condition at TAU(N), with the additional
  !    information taken from C(2,N).
  !
  implicit none
  integer ( kind = 4 ) n
  real ( kind = 8 ) c(4,n)
  real ( kind = 8 ) divdf1
  real ( kind = 8 ) divdf3
  real ( kind = 8 ) dtau
  real ( kind = 8 ) g
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ibcbeg
  integer ( kind = 4 ) ibcend
  real ( kind = 8 ) tau(n)
  !
  !  C(3,*) and C(4,*) are used initially for temporary storage.
  !
  !  Store first differences of the TAU sequence in C(3,*).
  !
  !  Store first divided difference of data in C(4,*).
  !
  do i = 2, n
     c(3,i) = tau(i) - tau(i-1)
  end do
  do i = 2, n 
     c(4,i) = ( c(1,i) - c(1,i-1) ) / ( tau(i) - tau(i-1) )
  end do
  !
  !  Construct the first equation from the boundary condition
  !  at the left endpoint, of the form:
  !
  !    C(4,1) * S(1) + C(3,1) * S(2) = C(2,1)
  !
  !  IBCBEG = 0: Not-a-knot
  !
  if ( ibcbeg == 0 ) then
     if ( n <= 2 ) then
        c(4,1) = 1.0D+00
        c(3,1) = 1.0D+00
        c(2,1) = 2.0D+00 * c(4,2)
        go to 120
     end if
     c(4,1) = c(3,3)
     c(3,1) = c(3,2) + c(3,3)
     c(2,1) = ( ( c(3,2) + 2.0D+00 * c(3,1) ) * c(4,2) * c(3,3) &
          + c(3,2)**2 * c(4,3) ) / c(3,1)
     !
     !  IBCBEG = 1: derivative specified.
     !
  else if ( ibcbeg == 1 ) then
     c(4,1) = 1.0D+00
     c(3,1) = 0.0D+00
     if ( n == 2 ) then
        go to 120
     end if
     !
     !  Second derivative prescribed at left end.
     !
  else
     c(4,1) = 2.0D+00
     c(3,1) = 1.0D+00
     c(2,1) = 3.0D+00 * c(4,2) - c(3,2) / 2.0D+00 * c(2,1)
     if ( n == 2 ) then
        go to 120
     end if
  end if
  !
  !  If there are interior knots, generate the corresponding
  !  equations and carry out the forward pass of Gauss elimination,
  !  after which the I-th equation reads:
  !
  !    C(4,I) * S(I) + C(3,I) * S(I+1) = C(2,I).
  !
  do i = 2, n-1
     g = -c(3,i+1) / c(4,i-1)
     c(2,i) = g * c(2,i-1) + 3.0D+00 * ( c(3,i) * c(4,i+1) + c(3,i+1) * c(4,i) )
     c(4,i) = g * c(3,i-1) + 2.0D+00 * ( c(3,i) + c(3,i+1))
  end do
  !
  !  Construct the last equation from the second boundary condition, of
  !  the form
  !
  !    -G * C(4,N-1) * S(N-1) + C(4,N) * S(N) = C(2,N)
  !
  !  If slope is prescribed at right end, one can go directly to
  !  back-substitution, since the C array happens to be set up just
  !  right for it at this point.
  !
  if ( ibcend == 1 ) then
     go to 160
  end if
  if ( 1 < ibcend ) then
     go to 110
  end if
90 continue
  !
  !  Not-a-knot and 3 <= N, and either 3 < N or also not-a-knot
  !  at left end point.
  !
  if ( n /= 3 .or. ibcbeg /= 0 ) then
     g = c(3,n-1) + c(3,n)
     c(2,n) = ( ( c(3,n) + 2.0D+00 * g ) * c(4,n) * c(3,n-1) + c(3,n)**2 &
          * ( c(1,n-1) - c(1,n-2) ) / c(3,n-1) ) / g
     g = - g / c(4,n-1)
     c(4,n) = c(3,n-1)
     c(4,n) = c(4,n) + g * c(3,n-1)
     c(2,n) = ( g * c(2,n-1) + c(2,n) ) / c(4,n)
     go to 160
  end if
  !
  !  N = 3 and not-a-knot also at left.
  !
100 continue
  c(2,n) = 2.0D+00 * c(4,n)
  c(4,n) = 1.0D+00
  g = -1.0D+00 / c(4,n-1)
  c(4,n) = c(4,n) - c(3,n-1) / c(4,n-1)
  c(2,n) = ( g * c(2,n-1) + c(2,n) ) / c(4,n)
  go to 160
  !
  !  IBCEND = 2: Second derivative prescribed at right endpoint.
  !
110 continue
  c(2,n) = 3.0D+00 * c(4,n) + c(3,n) / 2.0D+00 * c(2,n)
  c(4,n) = 2.0D+00
  g = -1.0D+00 / c(4,n-1)
  c(4,n) = c(4,n) - c(3,n-1) / c(4,n-1)
  c(2,n) = ( g * c(2,n-1) + c(2,n) ) / c(4,n)
  go to 160
  !
  !  N = 2.
  !
120 continue
  if ( ibcend == 2  ) then
     c(2,n) = 3.0D+00 * c(4,n) + c(3,n) / 2.0D+00 * c(2,n)
     c(4,n) = 2.0D+00
     g = -1.0D+00 / c(4,n-1)
     c(4,n) = c(4,n) - c(3,n-1) / c(4,n-1)
     c(2,n) = ( g * c(2,n-1) + c(2,n) ) / c(4,n)
  else if ( ibcend == 0 .and. ibcbeg /= 0 ) then
     c(2,n) = 2.0D+00 * c(4,n)
     c(4,n) = 1.0D+00
     g = -1.0D+00 / c(4,n-1)
     c(4,n) = c(4,n) - c(3,n-1) / c(4,n-1)
     c(2,n) = ( g * c(2,n-1) + c(2,n) ) / c(4,n)
  else if ( ibcend == 0 .and. ibcbeg == 0 ) then
     c(2,n) = c(4,n)
  end if
  !
  !  Back solve the upper triangular system 
  !
  !    C(4,I) * S(I) + C(3,I) * S(I+1) = B(I)
  !
  !  for the slopes C(2,I), given that S(N) is already known.
  !
160 continue
  do i = n-1, 1, -1
     c(2,i) = ( c(2,i) - c(3,i) * c(2,i+1) ) / c(4,i)
  end do
  !
  !  Generate cubic coefficients in each interval, that is, the
  !  derivatives at its left endpoint, from value and slope at its
  !  endpoints.
  !
  do i = 2, n
     dtau = c(3,i)
     divdf1 = ( c(1,i) - c(1,i-1) ) / dtau
     divdf3 = c(2,i-1) + c(2,i) - 2.0D+00 * divdf1
     c(3,i-1) = 2.0D+00 * ( divdf1 - c(2,i-1) - divdf3 ) / dtau
     c(4,i-1) = 6.0D+00 * divdf3 / dtau**2
  end do
  return
end subroutine cubspl



function ppvalu ( break, coef, l, k, x, jderiv )
  !***********************************************************************
  !
  !! PPVALU evaluates a piecewise polynomial function or its derivative.
  !
  !  Discussion:
  !
  !    PPVALU calculates the value at X of the JDERIV-th derivative of
  !    the piecewise polynomial function F from its piecewise
  !    polynomial representation.
  !
  !    The interval index I, appropriate for X, is found through a
  !    call to INTERV.  The formula for the JDERIV-th derivative
  !    of F is then evaluated by nested multiplication.
  !
  !    The J-th derivative of F is given by:
  !
  !      (d**J) F(X) = 
  !        COEF(J+1,I) + H * (
  !        COEF(J+2,I) + H * (
  !        ...
  !        COEF(K-1,I) + H * (
  !        COEF(K,  I) / (K-J-1) ) / (K-J-2) ... ) / 2 ) / 1
  !
  !    with
  !
  !      H = X - BREAK(I)
  !
  !    and
  !
  !      I = max ( 1, max ( J, BREAK(J) <= X, 1 <= J <= L ) ).
  !
  !  Modified:
  !
  !    16 February 2007
  !
  !  Author:
  !
  !    Carl DeBoor
  !
  !  Reference:
  !
  !    Carl DeBoor,
  !    A Practical Guide to Splines,
  !    Springer, 2001,
  !    ISBN: 0387953663,
  !    LC: QA1.A647.v27.
  !
  !  Parameters:
  !
  !    Input, real ( kind = 8 ) BREAK(L+1), real COEF(*), integer L, the
  !    piecewise polynomial representation of the function F to be evaluated.
  !
  !    Input, integer ( kind = 4 ) K, the order of the polynomial pieces that 
  !    make up the function F.  
  !    The usual value for K is 4, signifying a piecewise 
  !    cubic polynomial.
  !
  !    Input, real ( kind = 8 ) X, the point at which to evaluate F or
  !    of its derivatives.
  !
  !    Input, integer ( kind = 4 ) JDERIV, the order of the derivative to be
  !    evaluated.  If JDERIV is 0, then F itself is evaluated,
  !    which is actually the most common case.  It is assumed
  !    that JDERIV is zero or positive.
  !
  !    Output, real ( kind = 8 ) PPVALU, the value of the JDERIV-th
  !    derivative of F at X.
  !
  implicit none
  integer ( kind = 4 ) k
  integer ( kind = 4 ) l
  real ( kind = 8 ) break(l+1)
  real ( kind = 8 ) coef(k,l)
  real ( kind = 8 ) fmmjdr
  real ( kind = 8 ) h
  integer ( kind = 4 ) i
  integer ( kind = 4 ) jderiv
  integer ( kind = 4 ) m
  integer ( kind = 4 ) ndummy
  real ( kind = 8 ) ppvalu
  real ( kind = 8 ) value
  real ( kind = 8 ) x
  value = 0.0D+00
  fmmjdr = k - jderiv
  !
  !  Derivatives of order K or higher are identically zero.
  !
  if ( k <= jderiv ) then
     return
  end if
  !
  !  Find the index I of the largest breakpoint to the left of X.
  !
  call interv ( break, l+1, x, i, ndummy )
  !
  !  Evaluate the JDERIV-th derivative of the I-th polynomial piece at X.
  !
  h = x - break(i)
  m = k
  do
     value = ( value / fmmjdr ) * h + coef(m,i)
     m = m - 1
     fmmjdr = fmmjdr - 1.0D+00
     if ( fmmjdr <= 0.0D+00 ) then
        exit
     end if
  end do
  ppvalu = value
  return
end function ppvalu


subroutine interv ( xt, lxt, x, left, mflag )
  !********************************************************************
  !
  !! INTERV brackets a real value in an ascending vector of values.
  !
  !  Discussion:
  !
  !    The XT array is a set of increasing values.  The goal of the routine
  !    is to determine the largest index I so that XT(I) <= X.
  !
  !    The routine is designed to be efficient in the common situation
  !    that it is called repeatedly, with X taken from an increasing
  !    or decreasing sequence.
  !
  !    This will happen when a piecewise polynomial is to be graphed.
  !    The first guess for LEFT is therefore taken to be the value
  !    returned at the previous call and stored in the local variable ILO.
  !
  !    A first check ascertains that ILO < LXT.  This is necessary
  !    since the present call may have nothing to do with the previous
  !    call.  Then, if 
  !
  !      XT(ILO) <= X < XT(ILO+1), 
  !
  !    we set LEFT = ILO and are done after just three comparisons.
  !
  !    Otherwise, we repeatedly double the difference ISTEP = IHI - ILO
  !    while also moving ILO and IHI in the direction of X, until
  !
  !      XT(ILO) <= X < XT(IHI)
  !
  !    after which we use bisection to get, in addition, ILO + 1 = IHI.
  !    The value LEFT = ILO is then returned.
  !
  !  Modified:
  !
  !    14 February 2007
  !
  !  Author:
  !
  !    Carl DeBoor
  !
  !  Reference:
  !
  !    Carl DeBoor,
  !    A Practical Guide to Splines,
  !    Springer, 2001,
  !    ISBN: 0387953663,
  !    LC: QA1.A647.v27.
  !
  !  Parameters:
  !
  !    Input, real ( kind = 8 ) XT(LXT), a nondecreasing sequence of values.
  !
  !    Input, integer ( kind = 4 ) LXT, the dimension of XT.
  !
  !    Input, real ( kind = 8 ) X, the point whose location with 
  !    respect to the sequence XT is to be determined.
  !
  !    Output, integer ( kind = 4 ) LEFT, the index of the bracketing value:
  !      1     if             X  <  XT(1)
  !      I     if   XT(I)  <= X  < XT(I+1)
  !      LXT   if  XT(LXT) <= X
  !
  !    Output, integer ( kind = 4 ) MFLAG, indicates whether X lies within the
  !    range of the data.
  !    -1:            X  <  XT(1)
  !     0: XT(I)   <= X  < XT(I+1)
  !    +1: XT(LXT) <= X
  !
  implicit none
  integer ( kind = 4 ) lxt
  integer ( kind = 4 ) left
  integer ( kind = 4 ) mflag
  integer ( kind = 4 ) ihi
  integer ( kind = 4 ), save :: ilo = 1
  integer ( kind = 4 ) istep
  integer ( kind = 4 ) middle
  real ( kind = 8 ) x
  real ( kind = 8 ) xt(lxt)
  ihi = ilo + 1
  if ( lxt <= ihi ) then
     if ( xt(lxt) <= x ) then
        go to 110
     end if
     if ( lxt <= 1 ) then
        mflag = -1
        left = 1
        return
     end if
     ilo = lxt - 1
     ihi = lxt
  end if
  if ( xt(ihi) <= x ) then
     go to 20
  end if
  if ( xt(ilo) <= x ) then
     mflag = 0
     left = ilo
     return
  end if
  !
  !  Now X < XT(ILO).  Decrease ILO to capture X.
  !
  istep =  1
10 continue
  ihi = ilo
  ilo = ihi - istep
  if ( 1 < ilo ) then
     if ( xt(ilo) <= x ) then
        go to 50
     end if
     istep = istep * 2
     go to 10
  end if
  ilo = 1
  if ( x < xt(1) ) then
     mflag = -1
     left = 1
     return
  end if
  go to 50
  !
  !  Now XT(IHI) <= X.  Increase IHI to capture X.
  !
20 continue
  istep = 1
30 continue
  ilo = ihi
  ihi = ilo + istep
  if ( ihi < lxt ) then
     if ( x < xt(ihi) ) then
        go to 50
     end if
     istep = istep * 2
     go to 30
  end if
  if ( xt(lxt) <= x ) then
     go to 110
  end if
  !
  !  Now XT(ILO) < = X < XT(IHI).  Narrow the interval.
  !
  ihi = lxt
50 continue
  do
     middle = ( ilo + ihi ) / 2
     if ( middle == ilo ) then
        mflag = 0
        left = ilo
        return
     end if
     !
     !  It is assumed that MIDDLE = ILO in case IHI = ILO+1.
     !
     if ( xt(middle) <= x ) then
        ilo = middle
     else
        ihi = middle
     end if
  end do
  !
  !  Set output and return.
  !
110 continue
  mflag = 1
  if ( x == xt(lxt) ) then
     mflag = 0
  end if
  do left = lxt, 1, -1
     if ( xt(left) < xt(lxt) ) then
        return
     end if
  end do
  return
end subroutine interv
