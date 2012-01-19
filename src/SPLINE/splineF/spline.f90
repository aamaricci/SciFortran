subroutine basis_function_b_val ( tdata, tval, yval )

!*****************************************************************************80
!
!! BASIS_FUNCTION_B_VAL evaluates the B spline basis function.
!
!  Discussion:
!
!    The B spline basis function is a piecewise cubic which
!    has the properties that:
!
!    * it equals 2/3 at TDATA(3), 1/6 at TDATA(2) and TDATA(4);
!    * it is 0 for TVAL <= TDATA(1) or TDATA(5) <= TVAL;
!    * it is strictly increasing from TDATA(1) to TDATA(3),
!      and strictly decreasing from TDATA(3) to TDATA(5);
!    * the function and its first two derivatives are continuous
!      at each node TDATA(I).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    06 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Alan Davies, Philip Samuels,
!    An Introduction to Computational Geometry for Curves and Surfaces,
!    Clarendon Press, 1996,
!    ISBN: 0-19-851478-6,
!    LC: QA448.D38.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) TDATA(5), the nodes associated with the 
!    basis function.  The entries of TDATA are assumed to be distinct 
!    and increasing.
!
!    Input, real ( kind = 8 ) TVAL, a point at which the B spline basis 
!    function is to be evaluated.
!
!    Output, real ( kind = 8 ) YVAL, the value of the function at TVAL.
!
  implicit none

  integer ( kind = 4 ), parameter :: ndata = 5

  integer ( kind = 4 ) left
  integer ( kind = 4 ) right
  real    ( kind = 8 ) tdata(ndata)
  real    ( kind = 8 ) tval
  real    ( kind = 8 ) u
  real    ( kind = 8 ) yval

  if ( tval <= tdata(1) .or. tdata(ndata) <= tval ) then
    yval = 0.0D+00
    return
  end if
!
!  Find the interval [ TDATA(LEFT), TDATA(RIGHT) ] containing TVAL.
!
  call r8vec_bracket ( ndata, tdata, tval, left, right )
!
!  U is the normalized coordinate of TVAL in this interval.
!
  u = ( tval - tdata(left) ) / ( tdata(right) - tdata(left) )
!
!  Now evaluate the function.
!
  if ( tval < tdata(2) ) then
    yval = u**3 / 6.0D+00
  else if ( tval < tdata(3) ) then
    yval = ( ( (    - 3.0D+00   &
                * u + 3.0D+00 ) &
                * u + 3.0D+00 ) &
                * u + 1.0D+00 ) / 6.0D+00
  else if ( tval < tdata(4) ) then
    yval = ( ( (    + 3.0D+00   &
                * u - 6.0D+00 ) &
                * u + 0.0D+00 ) &
                * u + 4.0D+00 ) / 6.0D+00
  else if ( tval < tdata(5) ) then
    yval = ( 1.0D+00 - u )**3 / 6.0D+00
  end if

  return
end
subroutine basis_function_beta_val ( beta1, beta2, tdata, tval, yval )

!*****************************************************************************80
!
!! BASIS_FUNCTION_BETA_VAL evaluates the beta spline basis function.
!
!  Discussion:
!
!    With BETA1 = 1 and BETA2 = 0, the beta spline basis function 
!    equals the B spline basis function.
!
!    With BETA1 large, and BETA2 = 0, the beta spline basis function
!    skews to the right, that is, its maximum increases, and occurs
!    to the right of the center.
!
!    With BETA1 = 1 and BETA2 large, the beta spline becomes more like
!    a linear basis function; that is, its value in the outer two intervals
!    goes to zero, and its behavior in the inner two intervals approaches
!    a piecewise linear function that is 0 at the second node, 1 at the
!    third, and 0 at the fourth.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    09 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Alan Davies, Philip Samuels,
!    An Introduction to Computational Geometry for Curves and Surfaces,
!    Clarendon Press, 1996,
!    ISBN: 0-19-851478-6,
!    LC: QA448.D38.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) BETA1, the skew or bias parameter.
!    BETA1 = 1 for no skew or bias.
!
!    Input, real ( kind = 8 ) BETA2, the tension parameter.
!    BETA2 = 0 for no tension.
!
!    Input, real ( kind = 8 ) TDATA(5), the nodes associated with the 
!    basis function.  The entries of TDATA are assumed to be distinct 
!    and increasing.
!
!    Input, real ( kind = 8 ) TVAL, a point at which the B spline 
!    basis function is to be evaluated.
!
!    Output, real ( kind = 8 ) YVAL, the value of the function at TVAL.
!
  implicit none

  integer ( kind = 4 ), parameter :: ndata = 5

  real    ( kind = 8 ) a
  real    ( kind = 8 ) b
  real    ( kind = 8 ) beta1
  real    ( kind = 8 ) beta2
  real    ( kind = 8 ) c
  real    ( kind = 8 ) d
  integer ( kind = 4 ) left
  integer ( kind = 4 ) right
  real    ( kind = 8 ) tdata(ndata)
  real    ( kind = 8 ) tval
  real    ( kind = 8 ) u
  real    ( kind = 8 ) yval

  if ( tval <= tdata(1) .or. tdata(ndata) <= tval ) then
    yval = 0.0D+00
    return
  end if
!
!  Find the interval [ TDATA(LEFT), TDATA(RIGHT) ] containing TVAL.
!
  call r8vec_bracket ( ndata, tdata, tval, left, right )
!
!  U is the normalized coordinate of TVAL in this interval.
!
  u = ( tval - tdata(left) ) / ( tdata(right) - tdata(left) )
!
!  Now evaluate the function.
!
  if ( tval < tdata(2) ) then

    yval = 2.0D+00 * u**3 

  else if ( tval < tdata(3) ) then

    a = beta2 + 4.0D+00 * beta1 + 4.0D+00 * beta1 * beta1 &
      + 6.0D+00 * ( 1.0D+00 - beta1 * beta1 ) &
      - 3.0D+00 * ( 2.0D+00 + beta2 + 2.0D+00 * beta1 ) &
      + 2.0D+00 * ( 1.0D+00 + beta2 + beta1 + beta1 * beta1 )

    b = - 6.0D+00 * ( 1.0D+00 - beta1 * beta1 ) &
        + 6.0D+00 * ( 2.0D+00 + beta2 + 2.0D+00 * beta1 ) &
        - 6.0D+00 * ( 1.0D+00 + beta2 + beta1 + beta1 * beta1 )

    c = - 3.0D+00 * ( 2.0D+00 + beta2 + 2.0D+00 * beta1 ) &
        + 6.0D+00 * ( 1.0D+00 + beta2 + beta1 + beta1 * beta1 )

    d = - 2.0D+00 * ( 1.0D+00 + beta2 + beta1 + beta1 * beta1 )

    yval = ( ( d * u + c ) * u + b ) * u + a

  else if ( tval < tdata(4) ) then

    a = beta2 + 4.0D+00 * beta1 + 4.0D+00 * beta1 * beta1

    b = - 6.0D+00 * beta1 * ( 1.0D+00 - beta1 * beta1 )

    c = - 3.0D+00 * ( beta2 + 2.0D+00 * beta1**2 + 2.0D+00 * beta1**3 )

    d = 2.0D+00 * ( beta2 + beta1 + beta1**2 + beta1**3 )

    yval = ( ( d * u + c ) * u + b ) * u + a

  else if ( tval < tdata(5) ) then

    yval = 2.0D+00 * beta1**3 * ( 1.0D+00 - u )**3

  end if

  yval = yval / ( 2.0D+00 + beta2 + 4.0D+00 * beta1 + 4.0D+00 * beta1**2 &
    + 2.0D+00 * beta1**3 )

  return
end
subroutine basis_matrix_b_uni ( mbasis )

!*****************************************************************************80
!
!! BASIS_MATRIX_B_UNI sets up the uniform B spline basis matrix.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    17 June 2006
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    James Foley, Andries vanDam, Steven Feiner, John Hughes,
!    Computer Graphics, Principles and Practice,
!    Second Edition,
!    Addison Wesley, 1995,
!    ISBN: 0201848406,
!    LC: T385.C5735.
!
!  Parameters:
!
!    Output, real ( kind = 8 ) MBASIS(4,4), the basis matrix.
!
  implicit none

  real    ( kind = 8 ) mbasis(4,4)
!
!  In the following statement, the matrix appears as though it
!  has been transposed.
!
  mbasis(1:4,1:4) = real ( reshape ( &
    (/ -1,  3, -3, 1,                &
        3, -6,  0, 4,                &
       -3,  3,  3, 1,                &
        1,  0,  0, 0 /),             &
    (/ 4, 4 /) ), kind = 8 ) / 6.0D+00

  return
end
subroutine basis_matrix_beta_uni ( beta1, beta2, mbasis )

!*****************************************************************************80
!
!! BASIS_MATRIX_BETA_UNI sets up the uniform beta spline basis matrix.
!
!  Discussion:
!
!    If BETA1 = 1 and BETA2 = 0, then the beta spline reduces to 
!    the B spline.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    12 February 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    James Foley, Andries vanDam, Steven Feiner, John Hughes,
!    Computer Graphics, Principles and Practice,
!    Second Edition,
!    Addison Wesley, 1995,
!    ISBN: 0201848406,
!    LC: T385.C5735.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) BETA1, the skew or bias parameter.
!    BETA1 = 1 for no skew or bias.
!
!    Input, real ( kind = 8 ) BETA2, the tension parameter.
!    BETA2 = 0 for no tension.
!
!    Output, real ( kind = 8 ) MBASIS(4,4), the basis matrix.
!
  implicit none

  real    ( kind = 8 ) beta1
  real    ( kind = 8 ) beta2
  real    ( kind = 8 ) delta
  real    ( kind = 8 ) mbasis(4,4)

  mbasis(1,1) = - 2.0D+00 * beta1 * beta1 * beta1
  mbasis(1,2) =   2.0D+00 * beta2 &
    + 2.0 * beta1 * ( beta1 * beta1 + beta1 + 1.0D+00 )
  mbasis(1,3) = - 2.0D+00 * ( beta2 + beta1 * beta1 + beta1 + 1.0D+00 )
  mbasis(1,4) =   2.0D+00

  mbasis(2,1) =   6.0D+00 * beta1 * beta1 * beta1
  mbasis(2,2) = - 3.0D+00 * beta2 &
    - 6.0D+00 * beta1 * beta1 * ( beta1 + 1.0D+00 )
  mbasis(2,3) =   3.0D+00 * beta2 + 6.0D+00 * beta1 * beta1
  mbasis(2,4) =   0.0D+00

  mbasis(3,1) = - 6.0D+00 * beta1 * beta1 * beta1
  mbasis(3,2) =   6.0D+00 * beta1 * ( beta1 - 1.0D+00 ) * ( beta1 + 1.0D+00 )
  mbasis(3,3) =   6.0D+00 * beta1
  mbasis(3,4) =   0.0D+00

  mbasis(4,1) =   2.0D+00 * beta1 * beta1 * beta1
  mbasis(4,2) =   4.0D+00 * beta1 * ( beta1 + 1.0D+00 ) + beta2
  mbasis(4,3) =   2.0D+00
  mbasis(4,4) =   0.0D+00

  delta = ( ( 2.0D+00   &
    * beta1 + 4.0D+00 ) &
    * beta1 + 4.0D+00 ) &
    * beta1 + 2.0D+00 + beta2

  mbasis(1:4,1:4) = mbasis(1:4,1:4) / delta

  return
end
subroutine basis_matrix_bezier ( mbasis )

!*****************************************************************************80
!
!! BASIS_MATRIX_BEZIER sets up the cubic Bezier spline basis matrix.
!
!  Discussion:
!
!    This basis matrix assumes that the data points are stored as
!    ( P1, P2, P3, P4 ).  P1 is the function value at T = 0, while
!    P2 is used to approximate the derivative at T = 0 by
!    dP/dt = 3 * ( P2 - P1 ).  Similarly, P4 is the function value
!    at T = 1, and P3 is used to approximate the derivative at T = 1
!    by dP/dT = 3 * ( P4 - P3 ).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    05 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    James Foley, Andries vanDam, Steven Feiner, John Hughes,
!    Computer Graphics, Principles and Practice,
!    Second Edition,
!    Addison Wesley, 1995,
!    ISBN: 0201848406,
!    LC: T385.C5735.
!
!  Parameters:
!
!    Output, real ( kind = 8 ) MBASIS(4,4), the basis matrix.
!
  implicit none

  real    ( kind = 8 ) mbasis(4,4)

  mbasis(1,1) = -1.0D+00
  mbasis(1,2) =  3.0D+00
  mbasis(1,3) = -3.0D+00
  mbasis(1,4) =  1.0D+00

  mbasis(2,1) =  3.0D+00
  mbasis(2,2) = -6.0D+00
  mbasis(2,3) =  3.0D+00
  mbasis(2,4) =  0.0D+00

  mbasis(3,1) = -3.0D+00
  mbasis(3,2) =  3.0D+00
  mbasis(3,3) =  0.0D+00
  mbasis(3,4) =  0.0D+00

  mbasis(4,1) =  1.0D+00
  mbasis(4,2) =  0.0D+00
  mbasis(4,3) =  0.0D+00
  mbasis(4,4) =  0.0D+00

  return
end
subroutine basis_matrix_hermite ( mbasis )

!*****************************************************************************80
!
!! BASIS_MATRIX_HERMITE sets up the Hermite spline basis matrix.
!
!  Discussion:
!
!    This basis matrix assumes that the data points are stored as
!    ( P1, P2, P1', P2' ), with P1 and P1' being the data value and 
!    the derivative dP/dT at T = 0, while P2 and P2' apply at T = 1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    06 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    James Foley, Andries vanDam, Steven Feiner, John Hughes,
!    Computer Graphics, Principles and Practice,
!    Second Edition,
!    Addison Wesley, 1995,
!    ISBN: 0201848406,
!    LC: T385.C5735.
!
!  Parameters:
!
!    Output, real ( kind = 8 ) MBASIS(4,4), the basis matrix.
!
  implicit none

  real    ( kind = 8 ) mbasis(4,4)

  mbasis(1,1) =  2.0D+00
  mbasis(1,2) = -2.0D+00
  mbasis(1,3) =  1.0D+00
  mbasis(1,4) =  1.0D+00

  mbasis(2,1) = -3.0D+00
  mbasis(2,2) =  3.0D+00
  mbasis(2,3) = -2.0D+00
  mbasis(2,4) = -1.0D+00

  mbasis(3,1) =  0.0D+00
  mbasis(3,2) =  0.0D+00
  mbasis(3,3) =  1.0D+00
  mbasis(3,4) =  0.0D+00

  mbasis(4,1) =  1.0D+00
  mbasis(4,2) =  0.0D+00
  mbasis(4,3) =  0.0D+00
  mbasis(4,4) =  0.0D+00

  return
end
subroutine basis_matrix_overhauser_nonuni ( alpha, beta, mbasis )

!*****************************************************************************80
!
!! BASIS_MATRIX_OVERHAUSER_NONUNI: nonuniform Overhauser spline basis matrix.
!
!  Discussion:
!
!    This basis matrix assumes that the data points P1, P2, P3 and
!    P4 are not uniformly spaced in T, and that P2 corresponds to T = 0,
!    and P3 to T = 1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    05 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) ALPHA, BETA.
!    ALPHA = || P2 - P1 || / ( || P3 - P2 || + || P2 - P1 || )
!    BETA  = || P3 - P2 || / ( || P4 - P3 || + || P3 - P2 || ).
!
!    Output, real ( kind = 8 ) MBASIS(4,4), the basis matrix.
!
  implicit none

  real    ( kind = 8 ) alpha
  real    ( kind = 8 ) beta
  real    ( kind = 8 ) mbasis(4,4)

  mbasis(1,1) = - ( 1.0D+00 - alpha ) * ( 1.0D+00 - alpha ) / alpha
  mbasis(1,2) =   beta + ( 1.0D+00 - alpha ) / alpha
  mbasis(1,3) =   alpha - 1.0D+00 / ( 1.0D+00 - beta )
  mbasis(1,4) =   beta * beta / ( 1.0D+00 - beta )

  mbasis(2,1) =   2.0D+00 * ( 1.0D+00 - alpha ) * ( 1.0D+00 - alpha ) / alpha
  mbasis(2,2) = ( - 2.0D+00 * ( 1.0D+00 - alpha ) - alpha * beta ) / alpha
  mbasis(2,3) = ( 2.0D+00 * ( 1.0D+00 - alpha ) &
    - beta * ( 1.0D+00 - 2.0D+00 * alpha ) ) / ( 1.0D+00 - beta )
  mbasis(2,4) = - beta * beta / ( 1.0D+00 - beta )

  mbasis(3,1) = - ( 1.0D+00 - alpha ) * ( 1.0D+00 - alpha ) / alpha
  mbasis(3,2) =   ( 1.0D+00 - 2.0D+00 * alpha ) / alpha
  mbasis(3,3) =   alpha
  mbasis(3,4) =   0.0D+00

  mbasis(4,1) =   0.0D+00
  mbasis(4,2) =   1.0D+00
  mbasis(4,3) =   0.0D+00
  mbasis(4,4) =   0.0D+00

  return
end
subroutine basis_matrix_overhauser_nul ( alpha, mbasis )

!*****************************************************************************80
!
!! BASIS_MATRIX_OVERHAUSER_NUL sets the nonuniform left Overhauser basis matrix.
!
!  Discussion:
!
!    This basis matrix assumes that the data points P1, P2, and
!    P3 are not uniformly spaced in T, and that P1 corresponds to T = 0,
!    and P2 to T = 1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    27 August 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) ALPHA.
!    ALPHA = || P2 - P1 || / ( || P3 - P2 || + || P2 - P1 || )
!
!    Output, real ( kind = 8 ) MBASIS(3,3), the basis matrix.
!
  implicit none

  real    ( kind = 8 ) alpha
  real    ( kind = 8 ) mbasis(3,3)

  mbasis(1,1) =   1.0D+00 / alpha
  mbasis(1,2) = - 1.0D+00 / ( alpha * ( 1.0D+00 - alpha ) )
  mbasis(1,3) =   1.0D+00 / ( 1.0D+00 - alpha )

  mbasis(2,1) = - ( 1.0D+00 + alpha ) / alpha
  mbasis(2,2) =   1.0D+00 / ( alpha * ( 1.0D+00 - alpha ) )
  mbasis(2,3) = - alpha / ( 1.0D+00 - alpha )

  mbasis(3,1) =   1.0D+00
  mbasis(3,2) =   0.0D+00
  mbasis(3,3) =   0.0D+00

  return
end
subroutine basis_matrix_overhauser_nur ( beta, mbasis )

!*****************************************************************************80
!
!! BASIS_MATRIX_OVERHAUSER_NUR: the nonuniform right Overhauser basis matrix.
!
!  Discussion:
!
!    This basis matrix assumes that the data points PN-2, PN-1, and
!    PN are not uniformly spaced in T, and that PN-1 corresponds to T = 0,
!    and PN to T = 1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    27 August 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) BETA.
!    BETA = || P(N) - P(N-1) || / ( || P(N) - P(N-1) || + || P(N-1) - P(N-2) || )
!
!    Output, real ( kind = 8 ) MBASIS(3,3), the basis matrix.
!
  implicit none

  real    ( kind = 8 ) beta
  real    ( kind = 8 ) mbasis(3,3)

  mbasis(1,1) =   1.0D+00 / beta
  mbasis(1,2) = - 1.0D+00 / ( beta * ( 1.0D+00 - beta ) )
  mbasis(1,3) =   1.0D+00 / ( 1.0D+00 - beta )

  mbasis(2,1) = - ( 1.0D+00 + beta ) / beta
  mbasis(2,2) =   1.0D+00 / ( beta * ( 1.0D+00 - beta ) )
  mbasis(2,3) = - beta / ( 1.0D+00 - beta )

  mbasis(3,1) =   1.0D+00
  mbasis(3,2) =   0.0D+00
  mbasis(3,3) =   0.0D+00

  return
end
subroutine basis_matrix_overhauser_uni ( mbasis )

!*****************************************************************************80
!
!! BASIS_MATRIX_OVERHAUSER_UNI sets the uniform Overhauser spline basis matrix.
!
!  Discussion:
!
!    This basis matrix assumes that the data points P1, P2, P3 and
!    P4 are uniformly spaced in T, and that P2 corresponds to T = 0,
!    and P3 to T = 1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    05 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    James Foley, Andries vanDam, Steven Feiner, John Hughes,
!    Computer Graphics, Principles and Practice,
!    Second Edition,
!    Addison Wesley, 1995,
!    ISBN: 0201848406,
!    LC: T385.C5735.
!
!  Parameters:
!
!    Output, real ( kind = 8 ) MBASIS(4,4), the basis matrix.
!
  implicit none

  real    ( kind = 8 ) mbasis(4,4)

  mbasis(1,1) = - 1.0D+00 / 2.0D+00
  mbasis(1,2) =   3.0D+00 / 2.0D+00
  mbasis(1,3) = - 3.0D+00 / 2.0D+00
  mbasis(1,4) =   1.0D+00 / 2.0D+00

  mbasis(2,1) =   2.0D+00 / 2.0D+00
  mbasis(2,2) = - 5.0D+00 / 2.0D+00
  mbasis(2,3) =   4.0D+00 / 2.0D+00
  mbasis(2,4) = - 1.0D+00 / 2.0D+00

  mbasis(3,1) = - 1.0D+00 / 2.0D+00
  mbasis(3,2) =   0.0D+00
  mbasis(3,3) =   1.0D+00 / 2.0D+00
  mbasis(3,4) =   0.0D+00

  mbasis(4,1) =   0.0D+00
  mbasis(4,2) =   2.0D+00 / 2.0D+00
  mbasis(4,3) =   0.0D+00
  mbasis(4,4) =   0.0D+00

  return
end
subroutine basis_matrix_overhauser_uni_l ( mbasis )

!*****************************************************************************80
!
!! BASIS_MATRIX_OVERHAUSER_UNI_L sets the left uniform Overhauser basis matrix.
!
!  Discussion:
!
!    This basis matrix assumes that the data points P1, P2, and P3
!    are not uniformly spaced in T, and that P1 corresponds to T = 0,
!    and P2 to T = 1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    05 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) MBASIS(3,3), the basis matrix.
!
  implicit none

  real    ( kind = 8 ) mbasis(3,3)

  mbasis(1,1) =   2.0D+00
  mbasis(1,2) = - 4.0D+00
  mbasis(1,3) =   2.0D+00

  mbasis(2,1) = - 3.0D+00
  mbasis(2,2) =   4.0D+00
  mbasis(2,3) = - 1.0D+00

  mbasis(3,1) =   1.0D+00
  mbasis(3,2) =   0.0D+00
  mbasis(3,3) =   0.0D+00

  return
end
subroutine basis_matrix_overhauser_uni_r ( mbasis )

!*****************************************************************************80
!
!! BASIS_MATRIX_OVERHAUSER_UNI_R sets the right uniform Overhauser basis matrix.
!
!  Discussion:
!
!    This basis matrix assumes that the data points P(N-2), P(N-1), 
!    and P(N) are uniformly spaced in T, and that P(N-1) corresponds to 
!    T = 0, and P(N) to T = 1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    05 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) MBASIS(3,3), the basis matrix.
!
  implicit none

  real    ( kind = 8 ) mbasis(3,3)

  mbasis(1,1) =   2.0D+00
  mbasis(1,2) = - 4.0D+00
  mbasis(1,3) =   2.0D+00

  mbasis(2,1) = - 3.0D+00
  mbasis(2,2) =   4.0D+00
  mbasis(2,3) = - 1.0D+00

  mbasis(3,1) =   1.0D+00
  mbasis(3,2) =   0.0D+00
  mbasis(3,3) =   0.0D+00

  return
end
subroutine basis_matrix_tmp ( left, n, mbasis, ndata, tdata, ydata, tval, yval )

!*****************************************************************************80
!
!! BASIS_MATRIX_TMP computes Q = T * MBASIS * P
!
!  Discussion:
!
!    YDATA is a vector of data values, most frequently the values of some
!    function sampled at uniformly spaced points.  MBASIS is the basis
!    matrix for a particular kind of spline.  T is a vector of the
!    powers of the normalized difference between TVAL and the left
!    endpoint of the interval.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 February 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) LEFT, indicats that TVAL is in the interval
!    [ TDATA(LEFT), TDATA(LEFT+1) ], or that this is the "nearest"
!    interval to TVAL.
!    For TVAL < TDATA(1), use LEFT = 1.
!    For TDATA(NDATA) < TVAL, use LEFT = NDATA - 1.
!
!    Input, integer ( kind = 4 ) N, the order of the basis matrix.
!
!    Input, real ( kind = 8 ) MBASIS(N,N), the basis matrix.
!
!    Input, integer ( kind = 4 ) NDATA, the dimension of the vectors TDATA 
!    and YDATA.
!
!    Input, real ( kind = 8 ) TDATA(NDATA), the abscissa values.  This routine
!    assumes that the TDATA values are uniformly spaced, with an
!    increment of 1.0.
!
!    Input, real ( kind = 8 ) YDATA(NDATA), the data values to be 
!    interpolated or approximated.
!
!    Input, real ( kind = 8 ) TVAL, the value of T at which the spline is to be
!    evaluated.
!
!    Output, real ( kind = 8 ) YVAL, the value of the spline at TVAL.
!
  implicit none

  integer ( kind = 4 ), parameter :: maxn = 4
  integer ( kind = 4 ) n
  integer ( kind = 4 ) ndata

  real    ( kind = 8 ) arg
  integer ( kind = 4 ) first
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) left
  real    ( kind = 8 ) mbasis(n,n)
  real    ( kind = 8 ) tdata(ndata)
  real    ( kind = 8 ) tval
  real    ( kind = 8 ) tvec(maxn)
  real    ( kind = 8 ) ydata(ndata)
  real    ( kind = 8 ) yval

  if ( left == 1 ) then
    arg = 0.5D+00 * ( tval - tdata(left) )
    first = left
  else if ( left < ndata - 1 ) then
    arg = tval - tdata(left)
    first = left - 1
  else if ( left == ndata - 1 ) then
    arg = 0.5D+00 * ( 1.0D+00 + tval - tdata(left) )
    first = left - 1
  end if
!
!  TVEC(I) = ARG**(N-I).
!
  tvec(n) = 1.0D+00
  do i = n-1, 1, -1
    tvec(i) = arg * tvec(i+1)
  end do

  yval = 0.0D+00
  do j = 1, n
    yval = yval + dot_product ( tvec(1:n), mbasis(1:n,j) ) &
      * ydata(first - 1 + j)
  end do

  return
end
subroutine bc_val ( n, t, xcon, ycon, xval, yval )

!*****************************************************************************80
!
!! BC_VAL evaluates a parameterized N-th degree Bezier curve in 2D.
!
!  Discussion:
!
!    BC_VAL(T) is the value of a vector function of the form
!
!      BC_VAL(T) = ( X(T), Y(T) )
!
!    where
!
!      X(T) = sum ( 0 <= I <= N ) XCON(I) * BERN(I,N)(T)
!      Y(T) = sum ( 0 <= I <= N ) YCON(I) * BERN(I,N)(T)
!
!    BERN(I,N)(T) is the I-th Bernstein polynomial of order N
!    defined on the interval [0,1],
!
!    XCON(0:N) and YCON(0:N) are the coordinates of N+1 "control points".
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    12 February 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    David Kahaner, Cleve Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1989,
!    ISBN: 0-13-627258-4,
!    LC: TA345.K34.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the degree of the Bezier curve.
!    N must be at least 0.
!
!    Input, real ( kind = 8 ) T, the point at which the Bezier curve should
!    be evaluated.  The best results are obtained within the interval
!    [0,1] but T may be anywhere.
!
!    Input, real ( kind = 8 ) XCON(0:N), YCON(0:N), the X and Y coordinates
!    of the control points.  The Bezier curve will pass through
!    the points ( XCON(0), YCON(0) ) and ( XCON(N), YCON(N) ), but
!    generally NOT through the other control points.
!
!    Output, real ( kind = 8 ) XVAL, YVAL, the X and Y coordinates of the point
!    on the Bezier curve corresponding to the given T value.
!
  implicit none

  integer ( kind = 4 ) n

  real    ( kind = 8 ) bval(0:n)
  real    ( kind = 8 ) t
  real    ( kind = 8 ) xcon(0:n)
  real    ( kind = 8 ) xval
  real    ( kind = 8 ) ycon(0:n)
  real    ( kind = 8 ) yval

  call bp01 ( n, t, bval )

  xval = dot_product ( xcon(0:n), bval(0:n) )
  yval = dot_product ( ycon(0:n), bval(0:n) )

  return
end
function bez_val ( n, x, a, b, y )

!*****************************************************************************80
!
!! BEZ_VAL evaluates an N-th degree Bezier function at a point.
!
!  Discussion:
!
!    The Bezier function has the form:
!
!      BEZ(X) = sum ( 0 <= I <= N ) Y(I) * BERN(N,I)( (X-A)/(B-A) )
!
!    BERN(N,I)(X) is the I-th Bernstein polynomial of order N
!    defined on the interval [0,1], 
!
!    Y(0:N) is a set of coefficients,
!
!    and if, for I = 0 to N, we define the N+1 points
!
!      X(I) = ( (N-I)*A + I*B) / N, 
!
!    equally spaced in [A,B], the pairs ( X(I), Y(I) ) can be regarded as
!    "control points".
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    12 February 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    David Kahaner, Cleve Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1989,
!    ISBN: 0-13-627258-4,
!    LC: TA345.K34.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the degree of the Bezier function.
!    N must be at least 0.  
!
!    Input, real ( kind = 8 ) X, the point at which the Bezier function should
!    be evaluated.  The best results are obtained within the interval
!    [A,B] but X may be anywhere.
!
!    Input, real ( kind = 8 ) A, B, the interval over which the Bezier function
!    has been defined.  This is the interval in which the control
!    points have been set up.  Note BEZ(A) = Y(0) and BEZ(B) = Y(N),
!    although BEZ will not, in general pass through the other
!    control points.  A and B must not be equal.
!
!    Input, real ( kind = 8 ) Y(0:N), a set of data defining the Y coordinates
!    of the control points.
!
!    Output, real ( kind = 8 ) BEZ_VAL, the value of the Bezier function at X.
!
  implicit none

  integer ( kind = 4 ) n

  real    ( kind = 8 ) a
  real    ( kind = 8 ) b
  real    ( kind = 8 ) bez_val
  real    ( kind = 8 ) bval(0:n)
  real    ( kind = 8 ) x
  real    ( kind = 8 ) x01
  real    ( kind = 8 ) y(0:n)

  if ( b - a == 0.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'BEZ_VAL - Fatal error!'
    write ( *, '(a,g14.6)' ) '  Null interval, A = B = ', a
    stop
  end if
!
!  X01 lies in [0,1], in the same relative position as X in [A,B].
!
  x01 = ( x - a ) / ( b - a )

  call bp01 ( n, x01, bval )

  bez_val = dot_product ( y(0:n), bval(0:n) )

  return
end
subroutine bp_approx ( n, a, b, ydata, xval, yval )

!*****************************************************************************80
!
!! BP_APPROX evaluates the Bernstein polynomial approximant to F(X) on [A,B].
!
!  Formula:
!
!    BERN(F)(X) = sum ( 0 <= I <= N ) F(X(I)) * B_BASE(I,X)
!
!    where
!
!      X(I) = ( ( N - I ) * A + I * B ) / N
!      B_BASE(I,X) is the value of the I-th Bernstein basis polynomial at X.
!
!  Discussion:
!
!    The Bernstein polynomial BERN(F) for F(X) is an approximant, not an
!    interpolant; in other words, its value is not guaranteed to equal
!    that of F at any particular point.  However, for a fixed interval
!    [A,B], if we let N increase, the Bernstein polynomial converges
!    uniformly to F everywhere in [A,B], provided only that F is continuous.
!    Even if F is not continuous, but is bounded, the polynomial converges
!    pointwise to F(X) at all points of continuity.  On the other hand,
!    the convergence is quite slow compared to other interpolation
!    and approximation schemes.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    12 February 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    David Kahaner, Cleve Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1989,
!    ISBN: 0-13-627258-4,
!    LC: TA345.K34.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the degree of the Bernstein polynomial 
!    to be used.  N must be at least 0.
!
!    Input, real ( kind = 8 ) A, B, the endpoints of the interval on which the
!    approximant is based.  A and B should not be equal.
!
!    Input, real ( kind = 8 ) YDATA(0:N), the data values at N+1 equally 
!    spaced points in [A,B].  If N = 0, then the evaluation point should 
!    be 0.5 * ( A + B).  Otherwise, evaluation point I should be 
!    ( (N-I)*A + I*B ) / N ).
!
!    Input, real ( kind = 8 ) XVAL, the point at which the Bernstein polynomial 
!    approximant is to be evaluated.  XVAL does not have to lie in the 
!    interval [A,B].
!
!    Output, real ( kind = 8 ) YVAL, the value of the Bernstein polynomial
!    approximant for F, based in [A,B], evaluated at XVAL.
!
  implicit none

  integer ( kind = 4 ) n

  real    ( kind = 8 ) a
  real    ( kind = 8 ) b
  real    ( kind = 8 ) bvec(0:n)
  real    ( kind = 8 ) xval
  real    ( kind = 8 ) ydata(0:n)
  real    ( kind = 8 ) yval
!
!  Evaluate the Bernstein basis polynomials at XVAL.
!
  call bpab ( n, a, b, xval, bvec )
!
!  Now compute the sum of YDATA(I) * BVEC(I).
!
  yval = dot_product ( ydata(0:n), bvec(0:n) )

  return
end
subroutine bp01 ( n, x, bern )

!*****************************************************************************80
!
!! BP01 evaluates the Bernstein basis polynomials for [0,1] at a point.
!
!  Discussion:
!
!    For any N greater than or equal to 0, there is a set of N+1 Bernstein
!    basis polynomials, each of degree N, which form a basis for
!    all polynomials of degree N on [0,1].
!
!    BERN(N,I,X) = [N!/(I!*(N-I)!)] * (1-X)**(N-I) * X**I
!
!    N is the degree;
!
!    0 <= I <= N indicates which of the N+1 basis polynomials 
!    of degree N to choose;
!
!    X is a point in [0,1] at which to evaluate the basis polynomial.
!
!  First values:
!
!    B(0,0,X) = 1
!
!    B(1,0,X) =      1-X
!    B(1,1,X) =                X
!
!    B(2,0,X) =     (1-X)**2
!    B(2,1,X) = 2 * (1-X)    * X
!    B(2,2,X) =                X**2
!
!    B(3,0,X) =     (1-X)**3
!    B(3,1,X) = 3 * (1-X)**2 * X
!    B(3,2,X) = 3 * (1-X)    * X**2
!    B(3,3,X) =                X**3
!
!    B(4,0,X) =     (1-X)**4
!    B(4,1,X) = 4 * (1-X)**3 * X
!    B(4,2,X) = 6 * (1-X)**2 * X**2
!    B(4,3,X) = 4 * (1-X)    * X**3
!    B(4,4,X) =                X**4
!
!  Special values:
!
!    B(N,I,1/2) = C(N,K) / 2**N
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    12 February 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    David Kahaner, Cleve Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1989,
!    ISBN: 0-13-627258-4,
!    LC: TA345.K34.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the degree of the Bernstein basis polynomials.
!    N must be at least 0.
!
!    Input, real ( kind = 8 ) X, the evaluation point.
!
!    Output, real ( kind = 8 ) BERN(0:N), the values of the N+1 Bernstein basis
!    polynomials at X.
!
  implicit none

  integer ( kind = 4 ) n

  real    ( kind = 8 ) bern(0:n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real    ( kind = 8 ) x

  if ( n == 0 ) then

    bern(0) = 1.0D+00

  else if ( 0 < n ) then

    bern(0) = 1.0D+00 - x
    bern(1) = x

    do i = 2, n
      bern(i) = x * bern(i-1)
      do j = i-1, 1, -1
        bern(j) = x * bern(j-1) + ( 1.0D+00 - x ) * bern(j)
      end do
      bern(0) = ( 1.0D+00 - x ) * bern(0)
    end do

  end if

  return
end
subroutine bpab ( n, a, b, x, bern )

!*****************************************************************************80
!
!! BPAB evaluates the Bernstein basis polynomials for [A,B] at a point.
!
!  Discussion:
!
!    BERN(N,I,X) = [N!/(I!*(N-I)!)] * (B-X)**(N-I) * (X-A)**I / (B-A)**N
!
!    B(0,0,X) =   1
!
!    B(1,0,X) = (      B-X                ) / (B-A)
!    B(1,1,X) = (                 X-A     ) / (B-A)
!
!    B(2,0,X) = (     (B-X)**2            ) / (B-A)**2
!    B(2,1,X) = ( 2 * (B-X)    * (X-A)    ) / (B-A)**2
!    B(2,2,X) = (                (X-A)**2 ) / (B-A)**2
!
!    B(3,0,X) = (     (B-X)**3            ) / (B-A)**3
!    B(3,1,X) = ( 3 * (B-X)**2 * (X-A)    ) / (B-A)**3
!    B(3,2,X) = ( 3 * (B-X)    * (X-A)**2 ) / (B-A)**3
!    B(3,3,X) = (                (X-A)**3 ) / (B-A)**3
!
!    B(4,0,X) = (     (B-X)**4            ) / (B-A)**4
!    B(4,1,X) = ( 4 * (B-X)**3 * (X-A)    ) / (B-A)**4
!    B(4,2,X) = ( 6 * (B-X)**2 * (X-A)**2 ) / (B-A)**4
!    B(4,3,X) = ( 4 * (B-X)    * (X-A)**3 ) / (B-A)**4
!    B(4,4,X) = (                (X-A)**4 ) / (B-A)**4
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    12 February 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    David Kahaner, Cleve Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1989,
!    ISBN: 0-13-627258-4,
!    LC: TA345.K34.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the degree of the Bernstein basis polynomials.
!    There is a set of N+1 Bernstein basis polynomials, each of degree N, 
!    which form a basis for polynomials of degree N on [A,B].  N must
!    be at least 0.
!
!    Input, real ( kind = 8 ) A, B, the endpoints of the interval on which the
!    polynomials are to be based.  A and B should not be equal.
!
!    Input, real ( kind = 8 ) X, the point at which the polynomials are to be
!    evaluated.  X need not lie in the interval [A,B].
!
!    Output, real ( kind = 8 ) BERN(0:N), the values of the N+1 Bernstein basis
!    polynomials at X.
!
  implicit none

  integer ( kind = 4 ) n

  real    ( kind = 8 ) a
  real    ( kind = 8 ) b
  real    ( kind = 8 ) bern(0:n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real    ( kind = 8 ) x

  if ( b == a ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'BPAB - Fatal error!'
    write ( *, '(a,g14.6)' ) '  A = B = ', a
    stop
  end if

  if ( n == 0 ) then

    bern(0) = 1.0D+00

  else if ( 0 < n ) then

    bern(0) = ( b - x ) / ( b - a )
    bern(1) = ( x - a ) / ( b - a )

    do i = 2, n
      bern(i) = ( x - a ) * bern(i-1) / ( b - a )
      do j = i-1, 1, -1
        bern(j) = ( ( b - x ) * bern(j) + ( x - a ) * bern(j-1) ) / ( b - a )
      end do
      bern(0) = ( b - x ) * bern(0) / ( b - a )
    end do

  end if

  return
end
subroutine chfev ( x1, x2, f1, f2, d1, d2, ne, xe, fe, next, ierr )

!*****************************************************************************80
!
!! CHFEV evaluates a cubic polynomial given in Hermite form.
!
!  Discussion:
!
!    This routine evaluates a cubic polynomial given in Hermite form at an
!    array of points.  While designed for use by SPLINE_PCHIP_VAL, it may
!    be useful directly as an evaluator for a piecewise cubic
!    Hermite function in applications, such as graphing, where
!    the interval is known in advance.
!
!    The cubic polynomial is determined by function values
!    F1, F2 and derivatives D1, D2 on the interval [X1,X2].
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    25 June 2008
!
!  Author:  
!
!    Original FORTRAN77 version by Fred Fritsch.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Fred Fritsch, Ralph Carlson, 
!    Monotone Piecewise Cubic Interpolation,
!    SIAM Journal on Numerical Analysis,
!    Volume 17, Number 2, April 1980, pages 238-246.
!
!    David Kahaner, Cleve Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1989,
!    ISBN: 0-13-627258-4,
!    LC: TA345.K34.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X1, X2, the endpoints of the interval of
!    definition of the cubic.  X1 and X2 must be distinct.
!
!    Input, real ( kind = 8 ) F1, F2, the values of the function at X1 and
!    X2, respectively.
!
!    Input, real ( kind = 8 ) D1, D2, the derivative values at X1 and
!    X2, respectively.
!
!    Input, integer ( kind = 4 ) NE, the number of evaluation points.
!
!    Input, real ( kind = 8 ) XE(NE), the points at which the function is to
!    be evaluated.  If any of the XE are outside the interval
!    [X1,X2], a warning error is returned in NEXT.
!
!    Output, real ( kind = 8 ) FE(NE), the value of the cubic function
!    at the points XE.
!
!    Output, integer ( kind = 4 ) NEXT(2), indicates the number of 
!    extrapolation points:
!    NEXT(1) = number of evaluation points to the left of interval.
!    NEXT(2) = number of evaluation points to the right of interval.
!
!    Output, integer ( kind = 4 ) IERR, error flag.
!    0, no errors.
!    -1, NE < 1.
!    -2, X1 == X2.
!
  implicit none

  integer ( kind = 4 ) ne

  real    ( kind = 8 ) c2
  real    ( kind = 8 ) c3
  real    ( kind = 8 ) d1
  real    ( kind = 8 ) d2
  real    ( kind = 8 ) del1
  real    ( kind = 8 ) del2
  real    ( kind = 8 ) delta
  real    ( kind = 8 ) f1
  real    ( kind = 8 ) f2
  real    ( kind = 8 ) fe(ne)
  real    ( kind = 8 ) h
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierr
  integer ( kind = 4 ) next(2)
  real    ( kind = 8 ) x
  real    ( kind = 8 ) x1
  real    ( kind = 8 ) x2
  real    ( kind = 8 ) xe(ne)
  real    ( kind = 8 ) xma
  real    ( kind = 8 ) xmi

  if ( ne < 1 ) then
    ierr = -1
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'CHFEV - Fatal error!'
    write ( *, '(a)' ) '  Number of evaluation points is less than 1.'
    write ( *, '(a,i8)' ) '  NE = ', ne
    stop
  end if

  h = x2 - x1

  if ( h == 0.0D+00 ) then
    ierr = -2
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'CHFEV - Fatal error!'
    write ( *, '(a)' ) '  The interval [X1,X2] is of zero length.'
    stop
  end if
!
!  Initialize.
!
  ierr = 0
  next(1) = 0
  next(2) = 0
  xmi = min ( 0.0D+00, h )
  xma = max ( 0.0D+00, h )
!
!  Compute cubic coefficients expanded about X1.
!
  delta = ( f2 - f1 ) / h
  del1 = ( d1 - delta ) / h
  del2 = ( d2 - delta ) / h
  c2 = -( del1 + del1 + del2 )
  c3 = ( del1 + del2 ) / h
!
!  Evaluation loop.
!
  do i = 1, ne

    x = xe(i) - x1
    fe(i) = f1 + x * ( d1 + x * ( c2 + x * c3 ) )
!
!  Count the extrapolation points.
!
    if ( x < xmi ) then
      next(1) = next(1) + 1
    end if

    if ( xma < x ) then
      next(2) = next(2) + 1
    end if

  end do

  return
end
subroutine data_to_dif ( ntab, xtab, ytab, diftab )

!*****************************************************************************80
!
!! DATA_TO_DIF sets up a divided difference table from raw data.
!
!  Discussion:
!
!    Space can be saved by using a single array for both the DIFTAB and
!    YTAB dummy parameters.  In that case, the divided difference table will
!    overwrite the Y data without interfering with the computation.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    04 September 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Carl deBoor,
!    A Practical Guide to Splines,
!    Springer, 2001,
!    ISBN: 0387953663.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NTAB, the number of pairs of points
!    (XTAB(I),YTAB(I)) which are to be used as data.  The
!    number of entries to be used in DIFTAB, XTAB and YTAB.
!
!    Input, real ( kind = 8 ) XTAB(NTAB), the X values at which data was taken.
!    These values must be distinct.
!
!    Input, real ( kind = 8 ) YTAB(NTAB), the corresponding Y values.
!
!    Output, real ( kind = 8 ) DIFTAB(NTAB), the divided difference coefficients
!    corresponding to the input (XTAB,YTAB).
!
  implicit none

  integer ( kind = 4 ) ntab

  real    ( kind = 8 ) diftab(ntab)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  logical r8vec_distinct
  real    ( kind = 8 ) xtab(ntab)
  real    ( kind = 8 ) ytab(ntab)

  if ( .not. r8vec_distinct ( ntab, xtab ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'DATA_TO_DIF - Fatal error!'
    write ( *, '(a)' ) '  Two entries of XTAB are equal!'
    return
  end if
!
!  Copy the data values into DIFTAB.
!
  diftab(1:ntab) = ytab(1:ntab)
!
!  Compute the divided differences.
!
  do i = 2, ntab
    do j = ntab, i, -1

      diftab(j) = ( diftab(j) - diftab(j-1) ) / ( xtab(j) - xtab(j+1-i) )

    end do
  end do

  return
end
subroutine dif_val ( ntab, xtab, diftab, xval, yval )

!*****************************************************************************80
!
!! DIF_VAL evaluates a divided difference polynomial at a point.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    05 September 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Carl deBoor,
!    A Practical Guide to Splines,
!    Springer, 2001,
!    ISBN: 0387953663.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NTAB, the number of divided difference
!    coefficients in DIFTAB, and the number of points XTAB.
!
!    Input, real ( kind = 8 ) XTAB(NTAB), the X values upon which the
!    divided difference polynomial is based.
!
!    Input, real ( kind = 8 ) DIFTAB(NTAB), the divided difference 
!    polynomial coefficients.
!
!    Input, real ( kind = 8 ) XVAL, the value where the polynomial 
!    is to be evaluated.
!
!    Output, real ( kind = 8 ) YVAL, the value of the polynomial at XVAL.
!
  implicit none

  integer ( kind = 4 ) ntab

  real    ( kind = 8 ) diftab(ntab)
  integer ( kind = 4 ) i
  real    ( kind = 8 ) xtab(ntab)
  real    ( kind = 8 ) xval
  real    ( kind = 8 ) yval

  yval = diftab(ntab)
  do i = 1, ntab-1
    yval = diftab(ntab-i) + ( xval - xtab(ntab-i) ) * yval
  end do

  return
end
subroutine least_set_old ( ntab, xtab, ytab, ndeg, ptab, b, c, d, eps, ierror )

!*****************************************************************************80
!
!! LEAST_SET_OLD constructs the least squares polynomial approximation to data.
!
!  Discussion:
!
!    The least squares polynomial is not returned directly as a simple
!    polynomial.  Instead, it is represented in terms of a set of
!    orthogonal polynomials appopriate for the given data.  This makes
!    the computation more accurate, but means that the user can not
!    easily evaluate the computed polynomial.  Instead, the routine 
!    LEAST_EVAL should be used to evaluate the least squares polynomial
!    at any point.  (However, the value of the least squares polynomial
!    at each of the data points is returned as part of this computation.)
!
!
!    A discrete unweighted inner product is used, so that
!
!      ( F(X), G(X) ) = sum ( 1 <= I <= NTAB ) F(XTAB(I)) * G(XTAB(I)).
!
!    The least squares polynomial is determined using a set of
!    orthogonal polynomials PHI.  These polynomials can be defined
!    recursively by:
!
!      PHI(0)(X) = 1
!      PHI(1)(X) = X - B(1)
!      PHI(I)(X) = ( X - B(I) ) * PHI(I-1)(X) - D(I) * PHI(I-2)(X)
!
!    The array B(1:NDEG) contains the values
!
!      B(I) = ( X*PHI(I-1), PHI(I-1) ) / ( PHI(I-1), PHI(I-1) )
!
!    The array D(2:NDEG) contains the values
!
!      D(I) = ( PHI(I-1), PHI(I-1) ) / ( PHI(I-2), PHI(I-2) )
!
!    Using this basis, the least squares polynomial can be represented as
!
!      P(X)(I) = sum ( 0 <= I <= NDEG ) C(I) * PHI(I)(X)
!
!    The array C(0:NDEG) contains the values
!
!      C(I) = ( YTAB(I), PHI(I) ) / ( PHI(I), PHI(I) )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    15 May 2004
!
!  Author:
!
!    FORTRAN90 version by John Burkardt
!
!  Reference:
!
!    Gisela Engeln-Muellges, Frank Uhlig,
!    Numerical Algorithms with C,
!    Springer, 1996,
!    ISBN: 3-540-60530-4.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NTAB, the number of data points.
!
!    Input, real ( kind = 8 ) XTAB(NTAB), the X data.  The values in XTAB
!    should be distinct, and in increasing order.
!
!    Input, real ( kind = 8 ) YTAB(NTAB), the Y data values corresponding
!    to the X data in XTAB.
!
!    Input, integer ( kind = 4 ) NDEG, the degree of the polynomial which the
!    program is to use.  NDEG must be at least 0, and less than or 
!    equal to NTAB-1.
!
!    Output, real ( kind = 8 ) PTAB(NTAB), the value of the least 
!    squares polynomial at the points XTAB(1:NTAB).
!
!    Output, real ( kind = 8 ) B(1:NDEG), C(0:NDEG), D(2:NDEG), arrays 
!    needed to evaluate the polynomial.
!
!    Output, real ( kind = 8 ) EPS, the root-mean-square discrepancy of the
!    polynomial fit.
!
!    Output, integer ( kind = 4 ) IERROR, error flag.
!    zero, no error occurred;
!    nonzero, an error occurred, and the polynomial could not be computed.
!
  implicit none

  integer ( kind = 4 ) ndeg
  integer ( kind = 4 ) ntab

  real    ( kind = 8 ) b(1:ndeg)
  real    ( kind = 8 ) c(0:ndeg)
  real    ( kind = 8 ) d(2:ndeg)
  real    ( kind = 8 ) eps
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i0l1
  integer ( kind = 4 ) i1l1
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) it
  integer ( kind = 4 ) k
  integer ( kind = 4 ) mdeg
  real    ( kind = 8 ) ptab(ntab)
  real    ( kind = 8 ) rn0
  real    ( kind = 8 ) rn1
  real    ( kind = 8 ) s
  real    ( kind = 8 ) sum2
  real    ( kind = 8 ) xtab(ntab)
  real    ( kind = 8 ) y_sum
  real    ( kind = 8 ) ytab(ntab)
  real    ( kind = 8 ) ztab(2*ntab)

  ierror = 0
!
!  Check NDEG.
!
  if ( ndeg < 0 ) then
    ierror = 1
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'LEAST_SET_OLD - Fatal error!'
    write ( *, '(a)' ) '  NDEG < 0.'
    stop
  end if

  if ( ntab <= ndeg ) then
    ierror = 1
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'LEAST_SET_OLD - Fatal error!'
    write ( *, '(a)' ) '  NTAB <= NDEG.'
    stop
  end if
!
!  Check that the abscissas are strictly increasing.
!
  do i = 1, ntab-1
    if ( xtab(i+1) <= xtab(i) ) then
      ierror = 1
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'LEAST_SET_OLD - Fatal error!'
      write ( *, '(a)' ) '  XTAB must be strictly increasing, but'
      write ( *, '(a,i8,a,g14.6)' ) '  XTAB(', i, ') = ', xtab(i)
      write ( *, '(a,i8,a,g14.6)' ) '  XTAB(', i+1, ') = ', xtab(i+1)
      stop
    end if
  end do

  i0l1 = 0
  i1l1 = ntab
!
!  The polynomial is of degree at least 0.
!
  y_sum = sum ( ytab(1:ntab) )
  rn0 = ntab
  c(0) = y_sum / real ( ntab, kind = 8 )

  ptab(1:ntab) = y_sum / real ( ntab, kind = 8 )

  if ( ndeg == 0 ) then
    eps = sum ( ptab(1:ntab) - ytab(1:ntab) )**2
    eps = sqrt ( eps / real ( ntab, kind = 8 ) )
    return
  end if
!
!  The polynomial is of degree at least 1.
!
  b(1) = sum ( xtab(1:ntab) ) / real ( ntab, kind = 8 )

  s = 0.0D+00
  sum2 = 0.0D+00
  do i = 1, ntab
    ztab(i1l1+i) = xtab(i) - b(1)
    s = s + ztab(i1l1+i)**2
    sum2 = sum2 + ztab(i1l1+i) * ( ytab(i) - ptab(i) )
  end do

  rn1 = s
  c(1) = sum2 / s

  do i = 1, ntab
    ptab(i) = ptab(i) + c(1) * ztab(i1l1+i)
  end do

  if ( ndeg == 1 ) then
    eps = sum ( ( ptab(1:ntab) - ytab(1:ntab) )**2 )
    eps = sqrt ( eps / real ( ntab, kind = 8 ) )
    return
  end if

  ztab(1:ntab) = 1.0D+00

  mdeg = 2
  k = 2

  do

    d(k) = rn1 / rn0

    sum2 = 0.0D+00
    do i = 1, ntab
      sum2 = sum2 + xtab(i) * ztab(i1l1+i) * ztab(i1l1+i)
    end do

    b(k) = sum2 / rn1

    s = 0.0D+00
    sum2 = 0.0D+00
    do i = 1, ntab
      ztab(i0l1+i) = ( xtab(i) - b(k) ) * ztab(i1l1+i) &
        - d(k) * ztab(i0l1+i)
      s = s + ztab(i0l1+i) * ztab(i0l1+i)
      sum2 = sum2 + ztab(i0l1+i) * ( ytab(i) - ptab(i) )
    end do

    rn0 = rn1
    rn1 = s

    c(k) = sum2 / rn1

    it = i0l1
    i0l1 = i1l1
    i1l1 = it

    do i = 1, ntab
      ptab(i) = ptab(i) + c(k) * ztab(i1l1+i)
    end do

    if ( ndeg <= mdeg ) then
      exit
    end if

    mdeg = mdeg + 1
    k = k + 1

  end do
!
!  Compute the RMS error.
!
  eps = sum ( ( ptab(1:ntab) - ytab(1:ntab) )**2 )
  eps = sqrt ( eps / real ( ntab, kind = 8 ) )

  return
end
subroutine least_val_old ( x, ndeg, b, c, d, value )

!*****************************************************************************80
!
!! LEAST_VAL_OLD evaluates a least squares polynomial defined by LEAST_SET_OLD.
!
!  Discussion:
!
!    This is an "old" version of the routine.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    15 May 2004
!
!  Author:
!
!    FORTRAN90 version by John Burkardt
!
!  Reference:
!
!    Gisela Engeln-Muellges, Frank Uhlig,
!    Numerical Algorithms with C,
!    Springer, 1996,
!    ISBN: 3-540-60530-4.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the point at which the polynomial is 
!    to be evaluated.
!
!    Input, integer ( kind = 4 ) NDEG, the degree of the least squares 
!    polynomial.
!
!    Input, real ( kind = 8 ) B(1:NDEG), C(0:NDEG), D(2:NDEG), arrays 
!    defined by LEAST_SET_OLD, and needed to evaluate the polynomial.
!
!    Output, real ( kind = 8 ) VALUE, the value of the polynomial at X.
!
  implicit none

  integer ( kind = 4 ) ndeg

  real    ( kind = 8 ) b(1:ndeg)
  real    ( kind = 8 ) c(0:ndeg)
  real    ( kind = 8 ) d(2:ndeg)
  integer ( kind = 4 ) k
  real    ( kind = 8 ) sk
  real    ( kind = 8 ) skp1
  real    ( kind = 8 ) skp2
  real    ( kind = 8 ) value
  real    ( kind = 8 ) x

  if ( ndeg <= 0 ) then

    value = c(0)

  else if ( ndeg == 1 ) then

    value = c(0) + c(1) * ( x - b(1) )

  else

    skp2 = c(ndeg)
    skp1 = c(ndeg-1) + ( x - b(ndeg) ) * skp2

    do k = ndeg-2, 0, -1
      sk = c(k) + ( x - b(k+1) ) * skp1 - d(k+2) * skp2
      skp2 = skp1
      skp1 = sk
    end do

    value = sk

  end if

  return
end
subroutine least_set ( point_num, x, f, w, nterms, b, c, d )

!*****************************************************************************80
!
!! LEAST_SET defines a least squares polynomial for given data.
!
!  Discussion:
!
!    This routine is based on ORTPOL by Conte and deBoor.
!
!    The polynomial may be evaluated at any point X by calling LEAST_VAL.
!
!    Thanks to Andrew Telford for pointing out a mistake in the form of
!    the check that there are enough unique data points, 25 June 2008.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    25 June 2008
!
!  Author:
!
!    Original FORTRAN77 version by Samuel Conte, Carl deBoor.
!    FORTRAN90 version by John Burkardt
!
!  Reference:
!
!    Samuel Conte, Carl deBoor,
!    Elementary Numerical Analysis,
!    Second Edition,
!    McGraw Hill, 1972,
!    ISBN: 07-012446-4,
!    LC: QA297.C65.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) POINT_NUM, the number of data values.
!
!    Input, real ( kind = 8 ) X(POINT_NUM), the abscissas of the data points.
!    At least NTERMS of the values in X must be distinct.
!
!    Input, real ( kind = 8 ) F(POINT_NUM), the data values at the points X(*).
!
!    Input, real ( kind = 8 ) W(POINT_NUM), the weights associated with
!    the data points.  Each entry of W should be positive.
!
!    Input, integer ( kind = 4 ) NTERMS, the number of terms to use in the
!    approximating polynomial.  NTERMS must be at least 1.
!    The degree of the polynomial is NTERMS-1.
!
!    Output, real ( kind = 8 ) B(NTERMS), C(NTERMS), D(NTERMS), are quantities
!    defining the least squares polynomial for the input data,
!    which will be needed to evaluate the polynomial.
!
  implicit none

  integer ( kind = 4 ) point_num
  integer ( kind = 4 ) nterms

  real    ( kind = 8 ) b(nterms)
  real    ( kind = 8 ) c(nterms)
  real    ( kind = 8 ) d(nterms)
  real    ( kind = 8 ) f(point_num)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real    ( kind = 8 ) p
  real    ( kind = 8 ) pj(point_num)
  real    ( kind = 8 ) pjm1(point_num)
  real    ( kind = 8 ) s(nterms)
  real    ( kind = 8 ), parameter :: tol = 0.0D+00
  integer ( kind = 4 ) unique_num
  real    ( kind = 8 ) w(point_num)
  real    ( kind = 8 ) x(point_num)
!
!  Make sure at least NTERMS X values are unique.
!
  call r8vec_unique_count ( point_num, x, tol, unique_num )

  if ( unique_num < nterms ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'LEAST_SET - Fatal error!'
    write ( *, '(a)' ) '  The number of distinct X values must be'
    write ( *, '(a,i8)') '  at least NTERMS = ', nterms
    write ( *, '(a,i8)' ) '  but the input data has only ', unique_num
    write ( *, '(a)' ) '  distinct entries.'
    return
  end if
!
!  Make sure all W entries are positive.
!
  do i = 1, point_num
    if ( w(i) <= 0.0D+00 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'LEAST_SET - Fatal error!'
      write ( *, '(a)' ) '  All weights W must be positive,'
      write ( *, '(a,i8)' ) '  but weight ', i
      write ( *, '(a,g14.6)' ) '  is ', w(i)
      return
    end if
  end do
!
!  Start inner product summations at zero.
!
  b(1:nterms) = 0.0D+00
  c(1:nterms) = 0.0D+00
  d(1:nterms) = 0.0D+00
  s(1:nterms) = 0.0D+00
!
!  Set the values of P(-1,X) and P(0,X) at all data points.
!
  pjm1(1:point_num) = 0.0D+00
  pj(1:point_num) = 1.0D+00
!
!  Now compute the value of P(J,X(I)) as
!
!    P(J,X(I)) = ( X(I) - B(J) ) * P(J-1,X(I)) - C(J) * P(J-2,X(I))
!
!  where
!
!    S(J) = < P(J,X), P(J,X) >
!    B(J) = < x*P(J,X), P(J,X) > / < P(J,X), P(J,X) >
!    C(J) = S(J) / S(J-1)
!
!  and the least squares coefficients are
!
!    D(J) = < F(X), P(J,X) > / < P(J,X), P(J,X) >
!
  do j = 1, nterms

    d(j) = d(j) + sum ( w(1:point_num) * f(1:point_num) * pj(1:point_num) )
    b(j) = b(j) + sum ( w(1:point_num) * x(1:point_num) * pj(1:point_num)**2 )
    s(j) = s(j) + sum ( w(1:point_num) * pj(1:point_num)**2 )

    d(j) = d(j) / s(j)

    if ( j == nterms ) then
      c(j) = 0.0D+00
      return
    end if

    b(j) = b(j) / s(j)

    if ( j == 1 ) then
      c(j) = 0.0D+00
    else
      c(j) = s(j) / s(j-1)
    end if

    do i = 1, point_num
      p = pj(i)
      pj(i) = ( x(i) - b(j) ) * pj(i) - c(j) * pjm1(i)
      pjm1(i) = p
    end do

  end do

  return
end
subroutine least_val ( nterms, b, c, d, x, px )

!*****************************************************************************80
!
!! LEAST_VAL evaluates a least squares polynomial defined by LEAST_SET.
!
!  Discussion:
!
!    The least squares polynomial is assumed to be defined as a sum
!
!      P(X) = sum ( 1 <= I <= NTERMS ) D(I) * P(I-1,X)
!
!    where the orthogonal basis polynomials P(I,X) satisfy the following
!    three term recurrence:
!
!      P(-1,X) = 0
!      P(0,X) = 1
!      P(I,X) = ( X - B(I-1) ) * P(I-1,X) - C(I-1) * P(I-2,X)
!
!    Therefore, the least squares polynomial can be evaluated as follows:
!
!    If NTERMS is 1, then the value of P(X) is D(1) * P(0,X) = D(1).
!
!    Otherwise, P(X) is defined as the sum of NTERMS > 1 terms.  We can
!    reduce the number of terms by 1, because the polynomial P(NTERMS,X)
!    can be rewritten as a sum of polynomials;  Therefore, P(NTERMS,X)
!    can be eliminated from the sum, and its coefficient merged in with
!    those of other polynomials.  Repeat this process for P(NTERMS-1,X)
!    and so on until a single term remains.
!    P(NTERMS,X) of P(NTERMS-1,X)
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    24 May 2005
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Samuel Conte, Carl deBoor,
!    Elementary Numerical Analysis,
!    Second Edition,
!    McGraw Hill, 1972,
!    ISBN: 07-012446-4,
!    LC: QA297.C65.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NTERMS, the number of terms in the least 
!    squares polynomial.  NTERMS must be at least 1.  The input value of NTERMS
!    may be reduced from the value given to LEAST_SET.  This will
!    evaluate the least squares polynomial of the lower degree specified.
!
!    Input, real ( kind = 8 ) B(NTERMS), C(NTERMS), D(NTERMS), the information
!    computed by LEAST_SET.
!
!    Input, real ( kind = 8 ) X, the point at which the least squares polynomial
!    is to be evaluated.
!
!    Output, real ( kind = 8 ) PX, the value of the least squares 
!    polynomial at X.
!
  implicit none

  integer ( kind = 4 ) nterms

  real    ( kind = 8 ) b(nterms)
  real    ( kind = 8 ) c(nterms)
  real    ( kind = 8 ) d(nterms)
  integer ( kind = 4 ) i
  real    ( kind = 8 ) prev
  real    ( kind = 8 ) prev2
  real    ( kind = 8 ) px
  real    ( kind = 8 ) x

  px = d(nterms)
  prev = 0.0D+00

  do i = nterms - 1, 1, -1

    prev2 = prev
    prev = px

    if ( i == nterms-1 ) then
      px = d(i) + ( x - b(i) ) * prev
    else
      px = d(i) + ( x - b(i) ) * prev - c(i+1) * prev2
    end if

  end do

  return
end
subroutine least_val2 ( nterms, b, c, d, x, px, pxp )

!*****************************************************************************80
!
!! LEAST_VAL2 evaluates a least squares polynomial defined by LEAST_SET.
!
!  Discussion:
!
!    This routine also computes the derivative of the polynomial.
!
!    The least squares polynomial is assumed to be defined as a sum
!
!      P(X) = sum ( 1 <= I <= NTERMS ) D(I) * P(I-1,X)
!
!    where the orthogonal basis polynomials P(I,X) satisfy the following
!    three term recurrence:
!
!      P(-1,X) = 0
!      P(0,X) = 1
!      P(I,X) = ( X - B(I-1) ) * P(I-1,X) - C(I-1) * P(I-2,X)
!
!    Therefore, the least squares polynomial can be evaluated as follows:
!
!    If NTERMS is 1, then the value of P(X) is D(1) * P(0,X) = D(1).
!
!    Otherwise, P(X) is defined as the sum of NTERMS > 1 terms.  We can
!    reduce the number of terms by 1, because the polynomial P(NTERMS,X)
!    can be rewritten as a sum of polynomials;  Therefore, P(NTERMS,X)
!    can be eliminated from the sum, and its coefficient merged in with
!    those of other polynomials.  Repeat this process for P(NTERMS-1,X)
!    and so on until a single term remains.
!    P(NTERMS,X) of P(NTERMS-1,X)
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    24 May 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NTERMS, the number of terms in the least squares
!    polynomial.  NTERMS must be at least 1.  The value of NTERMS
!    may be reduced from the value given to LEAST_SET.
!    This will cause LEAST_VAL to evaluate the least squares polynomial
!    of the lower degree specified.
!
!    Input, real ( kind = 8 ) B(NTERMS), C(NTERMS), D(NTERMS), the information
!    computed by LEAST_SET.
!
!    Input, real ( kind = 8 ) X, the point at which the least squares polynomial
!    is to be evaluated.
!
!    Output, real ( kind = 8 ) PX, PXP, the value and derivative of the least
!    squares polynomial at X.
!
  implicit none

  integer ( kind = 4 ) nterms

  real    ( kind = 8 ) b(nterms)
  real    ( kind = 8 ) c(nterms)
  real    ( kind = 8 ) d(nterms)
  integer ( kind = 4 ) i
  real    ( kind = 8 ) px
  real    ( kind = 8 ) pxm1
  real    ( kind = 8 ) pxm2
  real    ( kind = 8 ) pxp
  real    ( kind = 8 ) pxpm1
  real    ( kind = 8 ) pxpm2
  real    ( kind = 8 ) x

  px = d(nterms)
  pxp = 0.0D+00
  pxm1 = 0.0D+00
  pxpm1 = 0.0D+00

  do i = nterms - 1, 1, -1

    pxm2 = pxm1
    pxpm2 = pxpm1
    pxm1 = px
    pxpm1 = pxp

    if ( i == nterms - 1 ) then
      px = d(i) + ( x - b(i) ) * pxm1
      pxp = pxm1 + ( x - b(i) ) * pxpm1
    else
      px = d(i) + ( x - b(i) ) * pxm1 - c(i+1) * pxm2
      pxp = pxm1 + ( x - b(i) ) * pxpm1 - c(i+1) * pxpm2
    end if

  end do

  return
end
subroutine parabola_val2 ( dim_num, ndata, tdata, ydata, left, tval, yval )

!*****************************************************************************80
!
!! PARABOLA_VAL2 evaluates a parabolic interpolant through tabular data.
!
!  Discussion:
!
!    This routine is a utility routine used by OVERHAUSER_SPLINE_VAL.
!    It constructs the parabolic interpolant through the data in
!    3 consecutive entries of a table and evaluates this interpolant
!    at a given abscissa value.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    26 January 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the dimension of a single data point.
!    DIM_NUM must be at least 1.
!
!    Input, integer ( kind = 4 ) NDATA, the number of data points.
!    NDATA must be at least 3.
!
!    Input, real ( kind = 8 ) TDATA(NDATA), the abscissas of the data 
!    points.  The values in TDATA must be in strictly ascending order.
!
!    Input, real ( kind = 8 ) YDATA(DIM_NUM,NDATA), the data points 
!    corresponding to the abscissas.
!
!    Input, integer ( kind = 4 ) LEFT, the location of the first of the three
!    consecutive data points through which the parabolic interpolant
!    must pass.  1 <= LEFT <= NDATA - 2.
!
!    Input, real ( kind = 8 ) TVAL, the value of T at which the parabolic
!    interpolant is to be evaluated.  Normally, TDATA(1) <= TVAL <= T(NDATA),
!    and the data will be interpolated.  For TVAL outside this range,
!    extrapolation will be used.
!
!    Output, real ( kind = 8 ) YVAL(DIM_NUM), the value of the parabolic
!    interpolant at TVAL.
!
  implicit none

  integer ( kind = 4 ) ndata
  integer ( kind = 4 ) dim_num

  real    ( kind = 8 ) dif1
  real    ( kind = 8 ) dif2
  integer ( kind = 4 ) i
  integer ( kind = 4 ) left
  real    ( kind = 8 ) t1
  real    ( kind = 8 ) t2
  real    ( kind = 8 ) t3
  real    ( kind = 8 ) tval
  real    ( kind = 8 ) tdata(ndata)
  real    ( kind = 8 ) ydata(dim_num,ndata)
  real    ( kind = 8 ) y1
  real    ( kind = 8 ) y2
  real    ( kind = 8 ) y3
  real    ( kind = 8 ) yval(dim_num)
!
!  Check.
!
  if ( left < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PARABOLA_VAL2 - Fatal error!'
    write ( *, '(a)' ) '  LEFT < 1.'
    write ( *, '(a,i8)' ) '  LEFT = ', left
    stop
  end if

  if ( ndata-2 < left ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PARABOLA_VAL2 - Fatal error!'
    write ( *, '(a)' ) '  NDATA-2 < LEFT.'
    write ( *, '(a,i8)' ) '  NDATA = ', ndata
    write ( *, '(a,i8)' ) '  LEFT =  ', left
    stop
  end if

  if ( dim_num < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PARABOLA_VAL2 - Fatal error!'
    write ( *, '(a)' ) '  DIM_NUM < 1.'
    stop
  end if
!
!  Copy out the three abscissas.
!
  t1 = tdata(left)
  t2 = tdata(left+1)
  t3 = tdata(left+2)

  if ( t2 <= t1 .or. t3 <= t2 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PARABOLA_VAL2 - Fatal error!'
    write ( *, '(a)' ) '  T2 <= T1 or T3 <= T2.'
    stop
  end if
!
!  Construct and evaluate a parabolic interpolant for the data
!  in each dimension.
!
  do i = 1, dim_num

    y1 = ydata(i,left)
    y2 = ydata(i,left+1)
    y3 = ydata(i,left+2)

    dif1 = ( y2 - y1 ) / ( t2 - t1 )
    dif2 = ( ( y3 - y1 ) / ( t3 - t1 ) &
         - ( y2 - y1 ) / ( t2 - t1 ) ) / ( t3 - t2 )

    yval(i) = y1 + ( tval - t1 ) * ( dif1 + ( tval - t2 ) * dif2 )

  end do

  return
end
function pchst ( arg1, arg2 )

!*****************************************************************************80
!
!! PCHST: PCHIP sign-testing routine.
!
!  Discussion:
!
!    This routine essentially computes the sign of ARG1 * ARG2.
!
!    The object is to do this without multiplying ARG1 * ARG2, to avoid
!    possible over/underflow problems.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    25 June 2008
!
!  Author:
!
!    Original FORTRAN77 version by Fred Fritsch.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Fred Fritsch, Ralph Carlson, 
!    Monotone Piecewise Cubic Interpolation,
!    SIAM Journal on Numerical Analysis,
!    Volume 17, Number 2, April 1980, pages 238-246.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) ARG1, ARG2, two values to check.
!
!    Output, real ( kind = 8 ) PCHST,
!    -1.0, if ARG1 and ARG2 are of opposite sign.
!     0.0, if either argument is zero.
!    +1.0, if ARG1 and ARG2 are of the same sign.
!
  implicit none

  real    ( kind = 8 ) arg1
  real    ( kind = 8 ) arg2
  real    ( kind = 8 ) pchst

  pchst = sign ( 1.0D+00, arg1 ) * sign ( 1.0D+00, arg2 )

  if ( arg1 == 0.0D+00 .or. arg2 == 0.0D+00 ) then
    pchst = 0.0D+00
  end if

  return
end
subroutine r8_swap ( x, y )

!*****************************************************************************80
!
!! R8_SWAP swaps two real values.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    01 May 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, real ( kind = 8 ) X, Y.  On output, the values of X and
!    Y have been interchanged.
!
  implicit none

  real    ( kind = 8 ) x
  real    ( kind = 8 ) y
  real    ( kind = 8 ) z

  z = x
  x = y
  y = z

  return
end
function r8_uniform_01 ( seed )

!*****************************************************************************80
!
!! R8_UNIFORM_01 is a portable pseudorandom number generator.
!
!  Discussion:
!
!    This routine implements the recursion
!
!      seed = 16807 * seed mod ( 2**31 - 1 )
!      unif = seed / ( 2**31 - 1 )
!
!    The integer arithmetic never requires more than 32 bits,
!    including a sign bit.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    11 August 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Paul Bratley, Bennett Fox, Linus Schrage,
!    A Guide to Simulation,
!    Springer Verlag, pages 201-202, 1983.
!
!    Bennett Fox,
!    Algorithm 647:
!    Implementation and Relative Efficiency of Quasirandom
!    Sequence Generators,
!    ACM Transactions on Mathematical Software,
!    Volume 12, Number 4, pages 362-376, 1986.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) SEED, the "seed" value, which should
!    NOT be 0. (Otherwise, the output values of SEED and UNIFORM will be zero.)
!    On output, SEED has been updated.
!
!    Output, real ( kind = 8 ) R8_UNIFORM_01, a new pseudorandom variate,
!    strictly between 0 and 1.
!
  implicit none

  integer ( kind = 4 ) k
  integer ( kind = 4 ) seed
  real    ( kind = 8 ) r8_uniform_01

  k = seed / 127773

  seed = 16807 * ( seed - k * 127773 ) - k * 2836

  if ( seed < 0 ) then
    seed = seed + 2147483647
  end if
!
!  Although SEED can be represented exactly as a 32 bit integer,
!  it generally cannot be represented exactly as a 32 bit real number!
!
  r8_uniform_01 = real ( seed, kind = 8 ) * 4.656612875E-10

  return
end
subroutine r83_mxv ( n, a, x, b )

!*****************************************************************************80
!
!! R83_MXV multiplies an R83 matrix times a vector.
!
!  Discussion:
!
!    The R83 storage format is used for a tridiagonal matrix.
!    The superdiagonal is stored in entries (1,2:N), the diagonal in
!    entries (2,1:N), and the subdiagonal in (3,1:N-1).  Thus, the
!    original matrix is "collapsed" vertically into the array.
!
!  Example:
!
!    Here is how an R83 matrix of order 5 would be stored:
!
!       *  A12 A23 A34 A45
!      A11 A22 A33 A44 A55
!      A21 A32 A43 A54  *
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    02 November 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the linear system.
!
!    Input, real ( kind = 8 ) A(3,N), the R83 matrix.
!
!    Input, real ( kind = 8 ) X(N), the vector to be multiplied by A.
!
!    Output, real ( kind = 8 ) B(N), the product A * x.
!
  implicit none

  integer ( kind = 4 ) n

  real    ( kind = 8 ) a(3,n)
  real    ( kind = 8 ) b(n)
  real    ( kind = 8 ) x(n)

  b(1:n)   =            a(2,1:n)   * x(1:n)
  b(1:n-1) = b(1:n-1) + a(1,2:n)   * x(2:n)
  b(2:n)   = b(2:n)   + a(3,1:n-1) * x(1:n-1)

  return
end
subroutine r83_np_fs ( n, a, b, x )

!*****************************************************************************80
!
!! R83_NP_FS factors and solves an R83 system.
!
!  Discussion:
!
!    The R83 storage format is used for a tridiagonal matrix.
!    The superdiagonal is stored in entries (1,2:N), the diagonal in
!    entries (2,1:N), and the subdiagonal in (3,1:N-1).  Thus, the
!    original matrix is "collapsed" vertically into the array.
!
!    This algorithm requires that each diagonal entry be nonzero.
!    It does not use pivoting, and so can fail on systems that
!    are actually nonsingular.
!
!    Here is how an R83 matrix of order 5 would be stored:
!
!       *  A12 A23 A34 A45
!      A11 A22 A33 A44 A55
!      A21 A32 A43 A54  *
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    02 November 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the linear system.
!
!    Input/output, real ( kind = 8 ) A(3,N).
!    On input, the tridiagonal matrix.
!    On output, the data in these vectors has been overwritten
!    by factorization information.
!
!    Input, real ( kind = 8 ) B(N), the right hand side of the linear system.
!
!    Output, real ( kind = 8 ) X(N), the solution of the linear system.
!
  implicit none

  integer ( kind = 4 ) n

  real    ( kind = 8 ) a(3,n)
  real    ( kind = 8 ) b(n)
  integer ( kind = 4 ) i
  real    ( kind = 8 ) x(n)
  real    ( kind = 8 ) xmult
!
!  The diagonal entries can't be zero.
!
  do i = 1, n
    if ( a(2,i) == 0.0D+00 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'R83_NP_FS - Fatal error!'
      write ( *, '(a,i8,a)' ) '  A(2,', i, ') = 0.'
      return
    end if
  end do

  x(1:n) = b(1:n)

  do i = 2, n
    xmult = a(3,i-1) / a(2,i-1)
    a(2,i) = a(2,i) - xmult * a(1,i)
    x(i)   = x(i)   - xmult * x(i-1)
  end do

  x(n) = x(n) / a(2,n)
  do i = n-1, 1, -1
    x(i) = ( x(i) - a(1,i+1) * x(i+1) ) / a(2,i)
  end do

  return
end
subroutine r83_uniform ( n, seed, a )

!*****************************************************************************80
!
!! R83_UNIFORM randomizes an R83 matrix.
!
!  Discussion:
!
!    The R83 storage format is used for a tridiagonal matrix.
!    The superdiagonal is stored in entries (1,2:N), the diagonal in
!    entries (2,1:N), and the subdiagonal in (3,1:N-1).  Thus, the
!    original matrix is "collapsed" vertically into the array.
!
!    Here is how an R83 matrix of order 5 would be stored:
!
!       *  A12 A23 A34 A45
!      A11 A22 A33 A44 A55
!      A21 A32 A43 A54  *
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    19 September 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the linear system.
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random number 
!    generator.
!
!    Output, real ( kind = 8 ) A(3,N), the R83 matrix.
!
  implicit none

  integer ( kind = 4 ) n

  real    ( kind = 8 ) a(3,n)
  integer ( kind = 4 ) seed

  a(1,1) = 0.0D+00
  call r8vec_uniform_01 ( n-1, seed, a(1,2:n) )

  call r8vec_uniform_01 ( n,   seed, a(2,1:n) )

  call r8vec_uniform_01 ( n-1, seed, a(3,1:n-1) )
  a(3,n) = 0.0D+00

  return
end
subroutine r8vec_bracket ( n, x, xval, left, right )

!*****************************************************************************80
!
!! R8VEC_BRACKET searches a sorted R8VEC for successive brackets of a value.
!
!  Discussion:
!
!    An R8VEC is an array of double precision real values.
!
!    If the values in the vector are thought of as defining intervals
!    on the real line, then this routine searches for the interval
!    nearest to or containing the given value.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    06 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, length of input array.
!
!    Input, real ( kind = 8 ) X(N), an array sorted into ascending order.
!
!    Input, real ( kind = 8 ) XVAL, a value to be bracketed.
!
!    Output, integer ( kind = 4 ) LEFT, RIGHT, the results of the search.
!    Either:
!      XVAL < X(1), when LEFT = 1, RIGHT = 2;
!      X(N) < XVAL, when LEFT = N-1, RIGHT = N;
!    or
!      X(LEFT) <= XVAL <= X(RIGHT).
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  integer ( kind = 4 ) left
  integer ( kind = 4 ) right
  real    ( kind = 8 ) x(n)
  real    ( kind = 8 ) xval

  do i = 2, n - 1

    if ( xval < x(i) ) then
      left = i - 1
      right = i
      return
    end if

   end do

  left = n - 1
  right = n

  return
end
subroutine r8vec_bracket3 ( n, t, tval, left )

!*****************************************************************************80
!
!! R8VEC_BRACKET3 finds the interval containing or nearest a given value.
!
!  Discussion:
!
!    An R8VEC is an array of double precision real values.
!
!    The routine always returns the index LEFT of the sorted array
!    T with the property that either
!    *  T is contained in the interval [ T(LEFT), T(LEFT+1) ], or
!    *  T < T(LEFT) = T(1), or
!    *  T > T(LEFT+1) = T(N).
!
!    The routine is useful for interpolation problems, where
!    the abscissa must be located within an interval of data
!    abscissas for interpolation, or the "nearest" interval
!    to the (extreme) abscissa must be found so that extrapolation
!    can be carried out.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    05 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, length of the input array.
!
!    Input, real ( kind = 8 ) T(N), an array sorted into ascending order.
!
!    Input, real ( kind = 8 ) TVAL, a value to be bracketed by entries of T.
!
!    Input/output, integer ( kind = 4 ) LEFT.
!
!    On input, if 1 <= LEFT <= N-1, LEFT is taken as a suggestion for the
!    interval [ T(LEFT), T(LEFT+1) ] in which TVAL lies.  This interval
!    is searched first, followed by the appropriate interval to the left
!    or right.  After that, a binary search is used.
!
!    On output, LEFT is set so that the interval [ T(LEFT), T(LEFT+1) ]
!    is the closest to TVAL; it either contains TVAL, or else TVAL
!    lies outside the interval [ T(1), T(N) ].
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) high
  integer ( kind = 4 ) left
  integer ( kind = 4 ) low
  integer ( kind = 4 ) mid
  real    ( kind = 8 ) t(n)
  real    ( kind = 8 ) tval
!
!  Check the input data.
!
  if ( n < 2 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8VEC_BRACKET3 - Fatal error!'
    write ( *, '(a)' ) '  N must be at least 2.'
    stop
  end if
!
!  If LEFT is not between 1 and N-1, set it to the middle value.
!
  if ( left < 1 .or. n - 1 < left ) then
    left = ( n + 1 ) / 2
  end if
!
!  CASE 1: TVAL < T(LEFT):
!  Search for TVAL in [T(I), T(I+1)] for intervals I = 1 to LEFT-1.
!
  if ( tval < t(left) ) then

    if ( left == 1 ) then
      return
    else if ( left == 2 ) then
      left = 1
      return
    else if ( t(left-1) <= tval ) then
      left = left - 1
      return
    else if ( tval <= t(2) ) then
      left = 1
      return
    end if
!
!  ...Binary search for TVAL in [T(I), T(I+1)] for intervals I = 2 to LEFT-2.
!
    low = 2
    high = left - 2

    do

      if ( low == high ) then
        left = low
        return
      end if

      mid = ( low + high + 1 ) / 2

      if ( t(mid) <= tval ) then
        low = mid
      else
        high = mid - 1
      end if

    end do
!
!  CASE2: T(LEFT+1) < TVAL:
!  Search for TVAL in {T(I),T(I+1)] for intervals I = LEFT+1 to N-1.
!
  else if ( t(left+1) < tval ) then

    if ( left == n - 1 ) then
      return
    else if ( left == n - 2 ) then
      left = left + 1
      return
    else if ( tval <= t(left+2) ) then
      left = left + 1
      return
    else if ( t(n-1) <= tval ) then
      left = n - 1
      return
    end if
!
!  ...Binary search for TVAL in [T(I), T(I+1)] for intervals I = LEFT+2 to N-2.
!
    low = left + 2
    high = n - 2

    do

      if ( low == high ) then
        left = low
        return
      end if

      mid = ( low + high + 1 ) / 2

      if ( t(mid) <= tval ) then
        low = mid
      else
        high = mid - 1
      end if

    end do
!
!  CASE3: T(LEFT) <= TVAL <= T(LEFT+1):
!  T is in [T(LEFT), T(LEFT+1)], as the user said it might be.
!
  else

  end if

  return
end
function r8vec_distinct ( n, x )

!*****************************************************************************80
!
!! R8VEC_DISTINCT is true if the entries in an R8VEC are distinct.
!
!  Discussion:
!
!    An R8VEC is an array of double precision real values.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    16 September 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in the vector.
!
!    Input, real ( kind = 8 ) X(N), the vector to be checked.
!
!    Output, logical R8VEC_DISTINCT is TRUE if all N elements of X 
!    are distinct.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  logical r8vec_distinct
  real    ( kind = 8 ) x(n)

  r8vec_distinct = .false.

  do i = 2, n
    do j = 1, i - 1 
      if ( x(i) == x(j) ) then
        return
      end if
    end do
  end do

  r8vec_distinct = .true.

  return
end
subroutine r8vec_even ( n, alo, ahi, a )

!*****************************************************************************80
!
!! R8VEC_EVEN returns N real values, evenly spaced between ALO and AHI.
!
!  Discussion:
!
!    An R8VEC is an array of double precision real values.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    17 February 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of values.
!
!    Input, real ( kind = 8 ) ALO, AHI, the low and high values.
!
!    Output, real ( kind = 8 ) A(N), N evenly spaced values.
!    Normally, A(1) = ALO and A(N) = AHI.
!    However, if N = 1, then A(1) = 0.5*(ALO+AHI).
!
  implicit none

  integer ( kind = 4 ) n

  real    ( kind = 8 ) a(n)
  real    ( kind = 8 ) ahi
  real    ( kind = 8 ) alo
  integer ( kind = 4 ) i

  if ( n == 1 ) then

    a(1) = 0.5D+00 * ( alo + ahi )

  else

    do i = 1, n
      a(i) = ( real ( n - i,     kind = 8 ) * alo   &
             + real (     i - 1, kind = 8 ) * ahi ) &
             / real ( n     - 1, kind = 8 )
    end do

  end if

  return
end
subroutine r8vec_indicator ( n, a )

!*****************************************************************************80
!
!! R8VEC_INDICATOR sets an R8VEC to the indicator vector.
!
!  Discussion:
!
!    An R8VEC is an array of double precision real values.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    09 February 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of elements of A.
!
!    Output, real ( kind = 8 ) A(N), the array to be initialized.
!
  implicit none

  integer ( kind = 4 ) n

  real    ( kind = 8 ) a(n)
  integer ( kind = 4 ) i

  do i = 1, n
    a(i) = real ( i, kind = 8 )
  end do

  return
end
subroutine r8vec_order_type ( n, a, order )

!*****************************************************************************80
!
!! R8VEC_ORDER_TYPE determines the order type of an R8VEC.
!
!  Discussion:
!
!    An R8VEC is an array of double precision real values.
!
!    We assume the array is increasing or decreasing, and we want to
!    verify that.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    20 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries of the array.
!
!    Input, real ( kind = 8 ) A(N), the array to be checked.
!
!    Output, integer ( kind = 4 ) ORDER, order indicator:
!    -1, no discernable order;
!    0, all entries are equal;
!    1, ascending order;
!    2, strictly ascending order;
!    3, descending order;
!    4, strictly descending order.
!
  implicit none

  integer ( kind = 4 ) n

  real    ( kind = 8 ) a(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) order
!
!  Search for the first value not equal to A(1).
!
  i = 1

  do

    i = i + 1

    if ( n < i ) then
      order = 0
      return
    end if

    if ( a(1) < a(i) ) then

      if ( i == 2 ) then
        order = 2
      else
        order = 1
      end if

      exit

    else if ( a(i) < a(1) ) then

      if ( i == 2 ) then
        order = 4
      else
        order = 3
      end if

      exit

    end if

  end do
!
!  Now we have a "direction".  Examine subsequent entries.
!
  do

    i = i + 1
    if ( n < i ) then
      exit
    end if

    if ( order == 1 ) then

      if ( a(i) < a(i-1) ) then
        order = -1
        exit
      end if

    else if ( order == 2 ) then

      if ( a(i) < a(i-1) ) then
        order = -1
        exit
      else if ( a(i) == a(i-1) ) then
        order = 1
      end if

    else if ( order == 3 ) then

      if ( a(i-1) < a(i) ) then
        order = -1
        exit
      end if

    else if ( order == 4 ) then

      if ( a(i-1) < a(i) ) then
        order = -1
        exit
      else if ( a(i) == a(i-1) ) then
        order = 3
      end if

    end if

  end do
 
  return
end
subroutine r8vec_print ( n, a, title )

!*****************************************************************************80
!
!! R8VEC_PRINT prints an R8VEC.
!
!  Discussion:
!
!    An R8VEC is an array of double precision real values.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    16 December 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of components of the vector.
!
!    Input, real ( kind = 8 ) A(N), the vector to be printed.
!
!    Input, character ( len = * ) TITLE, a title to be printed first.
!    TITLE may be blank.
!
  implicit none

  integer ( kind = 4 ) n

  real    ( kind = 8 ) a(n)
  integer ( kind = 4 ) i
  character ( len = * ) title

  if ( title /= ' ' ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) trim ( title )
  end if

  write ( *, '(a)' ) ' '
  do i = 1, n
    write ( *, '(i8,g14.6)' ) i, a(i)
  end do

  return
end
subroutine r8vec_sort_bubble_a ( n, a )

!*****************************************************************************80
!
!! R8VEC_SORT_BUBBLE_A ascending sorts an R8VEC using bubble sort.
!
!  Discussion:
!
!    An R8VEC is an array of double precision real values.
!
!    Bubble sort is simple to program, but inefficient.  It should not
!    be used for large arrays.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    01 February 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in the array.
!
!    Input/output, real ( kind = 8 ) A(N).
!    On input, an unsorted array.
!    On output, the array has been sorted.
!
  implicit none

  integer ( kind = 4 ) n

  real    ( kind = 8 ) a(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j

  do i = 1, n-1
    do j = i+1, n
      if ( a(j) < a(i) ) then
        call r8_swap ( a(i), a(j) )
      end if
    end do
  end do

  return
end
subroutine r8vec_uniform_01 ( n, seed, r )

!*****************************************************************************80
!
!! R8VEC_UNIFORM_01 returns a unit pseudorandom R8VEC.
!
!  Discussion:
!
!    An R8VEC is an array of double precision real values.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    19 August 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Paul Bratley, Bennett Fox, Linus Schrage,
!    A Guide to Simulation,
!    Springer Verlag, pages 201-202, 1983.
!
!    Bennett Fox,
!    Algorithm 647:
!    Implementation and Relative Efficiency of Quasirandom
!    Sequence Generators,
!    ACM Transactions on Mathematical Software,
!    Volume 12, Number 4, pages 362-376, 1986.
!
!    Peter Lewis, Allen Goodman, James Miller,
!    A Pseudo-Random Number Generator for the System/360,
!    IBM Systems Journal,
!    Volume 8, pages 136-143, 1969.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of entries in the vector.
!
!    Input/output, integer ( kind = 4 ) SEED, the "seed" value, which should
!    NOT be 0.  On output, SEED has been updated.
!
!    Output, real ( kind = 8 ) R(N), the vector of pseudorandom values.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  integer ( kind = 4 ) k
  integer ( kind = 4 ) seed
  real    ( kind = 8 ) r(n)

  do i = 1, n

    k = seed / 127773

    seed = 16807 * ( seed - k * 127773 ) - k * 2836

    if ( seed < 0 ) then
      seed = seed + 2147483647
    end if

    r(i) = real ( seed, kind = 8 ) * 4.656612875D-10

  end do

  return
end
subroutine r8vec_unique_count ( n, a, tol, unique_num )

!*****************************************************************************80
!
!! R8VEC_UNIQUE_COUNT counts the unique elements in an unsorted R8VEC.
!
!  Discussion:
!
!    An R8VEC is an array of double precision real values.
!
!    Because the array is unsorted, this algorithm is O(N^2).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    08 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of elements of A.
!
!    Input, real ( kind = 8 ) A(N), the unsorted array to examine.
!
!    Input, real ( kind = 8 ) TOL, a nonnegative tolerance for equality. 
!    Set it to 0.0 for the strictest test.
!
!    Output, integer ( kind = 4 ) UNIQUE_NUM, the number of unique elements of A.
!
  implicit none

  integer ( kind = 4 ) n

  real    ( kind = 8 ) a(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) unique_num
  real    ( kind = 8 ) tol

  unique_num = 0

  do i = 1, n

    unique_num = unique_num + 1

    do j = 1, i - 1

      if ( abs ( a(i) - a(j) ) <= tol ) then
        unique_num = unique_num - 1
        exit
      end if

    end do

  end do

  return
end
subroutine spline_b_val ( ndata, tdata, ydata, tval, yval )

!*****************************************************************************80
!
!! SPLINE_B_VAL evaluates a cubic B spline approximant.
!
!  Discussion:
!
!    The cubic B spline will approximate the data, but is not 
!    designed to interpolate it.
!
!    In effect, two "phantom" data values are appended to the data,
!    so that the spline will interpolate the first and last data values.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    11 February 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Carl deBoor,
!    A Practical Guide to Splines,
!    Springer, 2001,
!    ISBN: 0387953663.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NDATA, the number of data values.
!
!    Input, real ( kind = 8 ) TDATA(NDATA), the abscissas of the data.
!
!    Input, real ( kind = 8 ) YDATA(NDATA), the data values.
!
!    Input, real ( kind = 8 ) TVAL, a point at which the spline is 
!    to be evaluated.
!
!    Output, real ( kind = 8 ) YVAL, the value of the function at TVAL.
!
  implicit none

  integer ( kind = 4 ) ndata

  real    ( kind = 8 ) bval
  integer ( kind = 4 ) left
  integer ( kind = 4 ) right
  real    ( kind = 8 ) tdata(ndata)
  real    ( kind = 8 ) tval
  real    ( kind = 8 ) u
  real    ( kind = 8 ) ydata(ndata)
  real    ( kind = 8 ) yval
!
!  Find the nearest interval [ TDATA(LEFT), TDATA(RIGHT) ] to TVAL.
!
  call r8vec_bracket ( ndata, tdata, tval, left, right )
!
!  Evaluate the 5 nonzero B spline basis functions in the interval,
!  weighted by their corresponding data values.
!
  u = ( tval - tdata(left) ) / ( tdata(right) - tdata(left) )
  yval = 0.0D+00
!
!  B function associated with node LEFT - 1, (or "phantom node"),
!  evaluated in its 4th interval.
!
  bval = ( ( (     - 1.0D+00   &
               * u + 3.0D+00 ) &
               * u - 3.0D+00 ) &
               * u + 1.0D+00 ) / 6.0D+00

  if ( 0 < left-1 ) then
    yval = yval + ydata(left-1) * bval
  else
    yval = yval + ( 2.0D+00 * ydata(1) - ydata(2) ) * bval
  end if
!
!  B function associated with node LEFT,
!  evaluated in its third interval.
!
  bval = ( ( (       3.0D+00   &
               * u - 6.0D+00 ) &
               * u + 0.0D+00 ) &
               * u + 4.0D+00 ) / 6.0D+00

  yval = yval + ydata(left) * bval
!
!  B function associated with node RIGHT,
!  evaluated in its second interval.
!
  bval = ( ( (     - 3.0D+00   &
               * u + 3.0D+00 ) &
               * u + 3.0D+00 ) &
               * u + 1.0D+00 ) / 6.0D+00

  yval = yval + ydata(right) * bval
!
!  B function associated with node RIGHT+1, (or "phantom node"),
!  evaluated in its first interval.
!
  bval = u**3 / 6.0D+00

  if ( right+1 <= ndata ) then
    yval = yval + ydata(right+1) * bval
  else
    yval = yval + ( 2.0D+00 * ydata(ndata) - ydata(ndata-1) ) * bval
  end if

  return
end
subroutine spline_beta_val ( beta1, beta2, ndata, tdata, ydata, tval, yval )

!*****************************************************************************80
!
!! SPLINE_BETA_VAL evaluates a cubic beta spline approximant.
!
!  Discussion:
!
!    The cubic beta spline will approximate the data, but is not 
!    designed to interpolate it.
!
!    If BETA1 = 1 and BETA2 = 0, the cubic beta spline will be the
!    same as the cubic B spline approximant.
!
!    With BETA1 = 1 and BETA2 large, the beta spline becomes more like
!    a linear spline.
!
!    In effect, two "phantom" data values are appended to the data,
!    so that the spline will interpolate the first and last data values.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    11 February 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) BETA1, the skew or bias parameter.
!    BETA1 = 1 for no skew or bias.
!
!    Input, real ( kind = 8 ) BETA2, the tension parameter.
!    BETA2 = 0 for no tension.
!
!    Input, integer ( kind = 4 ) NDATA, the number of data values.
!
!    Input, real ( kind = 8 ) TDATA(NDATA), the abscissas of the data.
!
!    Input, real ( kind = 8 ) YDATA(NDATA), the data values.
!
!    Input, real ( kind = 8 ) TVAL, a point at which the spline is
!    to be evaluated.
!
!    Output, real ( kind = 8 ) YVAL, the value of the function at TVAL.
!
  implicit none

  integer ( kind = 4 ) ndata

  real    ( kind = 8 ) a
  real    ( kind = 8 ) b
  real    ( kind = 8 ) beta1
  real    ( kind = 8 ) beta2
  real    ( kind = 8 ) bval
  real    ( kind = 8 ) c
  real    ( kind = 8 ) d
  real    ( kind = 8 ) delta
  integer ( kind = 4 ) left
  integer ( kind = 4 ) right
  real    ( kind = 8 ) tdata(ndata)
  real    ( kind = 8 ) tval
  real    ( kind = 8 ) u
  real    ( kind = 8 ) ydata(ndata)
  real    ( kind = 8 ) yval
!
!  Find the nearest interval [ TDATA(LEFT), TDATA(RIGHT) ] to TVAL.
!
  call r8vec_bracket ( ndata, tdata, tval, left, right )
!
!  Evaluate the 5 nonzero beta spline basis functions in the interval,
!  weighted by their corresponding data values.
!
  u = ( tval - tdata(left) ) / ( tdata(right) - tdata(left) )

  delta = ( ( 2.0D+00   &
    * beta1 + 4.0D+00 ) &
    * beta1 + 4.0D+00 ) &
    * beta1 + 2.0D+00 + beta2

  yval = 0.0D+00
!
!  Beta function associated with node LEFT - 1, (or "phantom node"),
!  evaluated in its 4th interval.
!
  bval = 2.0D+00 * ( beta1 * ( 1.0D+00 - u ) )**3  / delta

  if ( 0 < left-1 ) then
    yval = yval + ydata(left-1) * bval
  else
    yval = yval + ( 2.0D+00 * ydata(1) - ydata(2) ) * bval
  end if
!
!  Beta function associated with node LEFT,
!  evaluated in its third interval.
!
  a = beta2 + ( 4.0D+00 + 4.0D+00 * beta1 ) * beta1

  b = - 6.0D+00 * beta1 * ( 1.0D+00 - beta1 ) * ( 1.0D+00 + beta1 )

  c = ( (     - 6.0D+00   &
      * beta1 - 6.0D+00 ) &
      * beta1 + 0.0D+00 ) &
      * beta1 - 3.0D+00 * beta2

  d = ( (     + 2.0D+00   &
      * beta1 + 2.0D+00 ) &
      * beta1 + 2.0D+00 ) &
      * beta1 + 2.0D+00 * beta2

  bval = ( a + u * ( b + u * ( c + u * d ) ) ) / delta

  yval = yval + ydata(left) * bval
!
!  Beta function associated with node RIGHT,
!  evaluated in its second interval.
!
  a = 2.0D+00

  b = 6.0D+00 * beta1

  c = 3.0D+00 * beta2 + 6.0D+00 * beta1 * beta1

  d = - 2.0D+00 * ( 1.0D+00 + beta2 + beta1 + beta1 * beta1 )

  bval = ( a + u * ( b + u * ( c + u * d ) ) ) / delta

  yval = yval + ydata(right) * bval
!
!  Beta function associated with node RIGHT+1, (or "phantom node"),
!  evaluated in its first interval.
!
  bval = 2.0D+00 * u**3 / delta

  if ( right + 1 <= ndata ) then
    yval = yval + ydata(right+1) * bval
  else
    yval = yval + ( 2.0D+00 * ydata(ndata) - ydata(ndata-1) ) * bval
  end if

  return
end
subroutine spline_bezier_val ( dim_num, interval_num, data_val, point_num, &
  point_t, point_val )

!*****************************************************************************80
!
!! SPLINE_BEZIER_VAL evaluates a cubic Bezier spline.
!
!  Discussion:
!
!    The cubic Bezier spline of N parts is defined by choosing 
!    3*N+1 equally spaced T-abscissa values in the interval [0,N].
!    This defines N subintervals, each of length 1, and each containing
!    4 successives T abscissa values.
!
!    At each abscissa value, a DIM_NUM-dimensional Bezier control
!    value is assigned.  Over each interval, a Bezier cubic function
!    is used to define the value of the Bezier spline.  To the left of
!    the first interval, or to the right of the last interval,
!    extrapolation may be used to extend the spline definition to
!    the entire real line.
!
!    Note that the Bezier spline will pass through the 1st, 4th,
!    and in general 3*I+1 control values exactly.  The other control
!    values are not interpolating points.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    17 June 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the spatial dimension.
!
!    Input, integer ( kind = 4 ) INTERVAL_NUM, the number of intervals.
!
!    Input, real ( kind = 8 ) DATA_VAL(DIM_NUM,3*INTERVAL_NUM+1), the control
!    values.
!
!    Input, integer ( kind = 4 ) POINT_NUM, the number of sample points at 
!    which the Bezier cubic spline is to be evaluated.
!
!    Input, real ( kind = 8 ) POINT_T(POINT_NUM), the "T" values associated
!    with the points.  A value of T between 0 and 1, for instance,
!    is associated with the first interval, and a value of T between
!    INTERVAL_NUM-1 and INTERVAL_NUM is in the last interval. 
!
!    Output, real ( kind = 8 ) POINT_VAL(DIM_NUM,POINT_NUM), the value
!    of the Bezier cubic spline at the sample points.
!
  implicit none

  integer ( kind = 4 ), parameter :: cubic = 3
  integer ( kind = 4 ) interval_num
  integer ( kind = 4 ) dim_num
  integer ( kind = 4 ) point_num

  real    ( kind = 8 ) bernstein_val(0:cubic)
  real    ( kind = 8 ) data_val(dim_num,cubic*interval_num+1)
  integer ( kind = 4 ) dim
  integer ( kind = 4 ) interval
  integer ( kind = 4 ) offset
  integer ( kind = 4 ) point
  real    ( kind = 8 ) point_t(point_num)
  real    ( kind = 8 ) point_val(dim_num,point_num)
  real    ( kind = 8 ) t
  real    ( kind = 8 ) t_01

  do point = 1, point_num

    t = point_t(point)

    interval = int ( t + 1 )

    interval = max ( interval, 1 )
    interval = min ( interval, interval_num )

    offset = 1 + ( interval - 1 ) * cubic

    t_01 = t - real ( interval - 1, kind = 8 )

    call bp01 ( cubic, t_01, bernstein_val )

    do dim = 1, dim_num
      point_val(dim,point) = dot_product ( &
        data_val(dim,offset:offset+cubic), bernstein_val(0:cubic) )
    end do

  end do

  return
end
subroutine spline_constant_val ( ndata, tdata, ydata, tval, yval )

!*****************************************************************************80
!
!! SPLINE_CONSTANT_VAL evaluates a piecewise constant spline at a point.
!
!  Discussion:
!
!    NDATA-1 points TDATA define NDATA intervals, with the first
!    and last being semi-infinite.
!
!    The value of the spline anywhere in interval I is YDATA(I).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    16 November 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NDATA, the number of data points defining 
!    the spline.  NDATA must be at least 1.
!
!    Input, real ( kind = 8 ) TDATA(NDATA-1), the breakpoints.  The values 
!    of TDATA should be distinct and increasing.
!
!    Input, real ( kind = 8 ) YDATA(NDATA), the values of the spline in
!    the intervals defined by the breakpoints.
!
!    Input, real ( kind = 8 ) TVAL, the point at which the spline is 
!    to be evaluated.
!
!    Output, real ( kind = 8 ) YVAL, the value of the spline at TVAL.  
!
  implicit none

  integer ( kind = 4 ) ndata

  integer ( kind = 4 ) i
  real    ( kind = 8 ) tdata(ndata-1)
  real    ( kind = 8 ) tval
  real    ( kind = 8 ) ydata(ndata)
  real    ( kind = 8 ) yval
!
!  Check NDATA.
!
  if ( ndata < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'SPLINE_CONSTANT_VAL - Fatal error!'
    write ( *, '(a)' ) '  NDATA < 1.'
    stop
  end if

  do i = 1, ndata-1
    if ( tval <= tdata(i) ) then
      yval = ydata(i)
      return
    end if
  end do

  yval = ydata(ndata)

  return
end
subroutine spline_cubic_set ( n, t, y, ibcbeg, ybcbeg, ibcend, ybcend, ypp )

!*****************************************************************************80
!
!! SPLINE_CUBIC_SET computes the second derivatives of a piecewise cubic spline.
!
!  Discussion:
!
!    For data interpolation, the user must call SPLINE_CUBIC_SET to 
!    determine the second derivative data, passing in the data to be 
!    interpolated, and the desired boundary conditions.
!
!    The data to be interpolated, plus the SPLINE_CUBIC_SET output, 
!    defines the spline.  The user may then call SPLINE_CUBIC_VAL to 
!    evaluate the spline at any point.
!
!    The cubic spline is a piecewise cubic polynomial.  The intervals
!    are determined by the "knots" or abscissas of the data to be
!    interpolated.  The cubic spline has continous first and second
!    derivatives over the entire interval of interpolation.  
!
!    For any point T in the interval T(IVAL), T(IVAL+1), the form of
!    the spline is
!
!      SPL(T) = A(IVAL)
!             + B(IVAL) * ( T - T(IVAL) ) 
!             + C(IVAL) * ( T - T(IVAL) )**2
!             + D(IVAL) * ( T - T(IVAL) )**3
!
!    If we assume that we know the values Y(*) and YPP(*), which represent
!    the values and second derivatives of the spline at each knot, then
!    the coefficients can be computed as:
!
!      A(IVAL) = Y(IVAL)
!      B(IVAL) = ( Y(IVAL+1) - Y(IVAL) ) / ( T(IVAL+1) - T(IVAL) )
!        - ( YPP(IVAL+1) + 2 * YPP(IVAL) ) * ( T(IVAL+1) - T(IVAL) ) / 6
!      C(IVAL) = YPP(IVAL) / 2
!      D(IVAL) = ( YPP(IVAL+1) - YPP(IVAL) ) / ( 6 * ( T(IVAL+1) - T(IVAL) ) )
!
!    Since the first derivative of the spline is
!
!      SPL'(T) =     B(IVAL)
!              + 2 * C(IVAL) * ( T - T(IVAL) )
!              + 3 * D(IVAL) * ( T - T(IVAL) )**2,
!
!    the requirement that the first derivative be continuous at interior
!    knot I results in a total of N-2 equations, of the form:
!
!      B(IVAL-1) + 2 C(IVAL-1) * (T(IVAL)-T(IVAL-1)) 
!      + 3 * D(IVAL-1) * (T(IVAL) - T(IVAL-1))**2 = B(IVAL)
!
!    or, setting H(IVAL) = T(IVAL+1) - T(IVAL)
!
!      ( Y(IVAL) - Y(IVAL-1) ) / H(IVAL-1)
!      - ( YPP(IVAL) + 2 * YPP(IVAL-1) ) * H(IVAL-1) / 6
!      + YPP(IVAL-1) * H(IVAL-1)
!      + ( YPP(IVAL) - YPP(IVAL-1) ) * H(IVAL-1) / 2
!      = 
!      ( Y(IVAL+1) - Y(IVAL) ) / H(IVAL)
!      - ( YPP(IVAL+1) + 2 * YPP(IVAL) ) * H(IVAL) / 6
!
!    or
!
!      YPP(IVAL-1) * H(IVAL-1) + 2 * YPP(IVAL) * ( H(IVAL-1) + H(IVAL) )
!      + YPP(IVAL) * H(IVAL) 
!      =
!      6 * ( Y(IVAL+1) - Y(IVAL) ) / H(IVAL)
!      - 6 * ( Y(IVAL) - Y(IVAL-1) ) / H(IVAL-1)    
!
!    Boundary conditions must be applied at the first and last knots.  
!    The resulting tridiagonal system can be solved for the YPP values.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    27 January 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Carl deBoor,
!    A Practical Guide to Splines,
!    Springer, 2001,
!    ISBN: 0387953663.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of data points; N must be 
!    at least 2. 
!
!    Input, real ( kind = 8 ) T(N), the points where data is specified.  
!    The values should be distinct, and increasing.
!
!    Input, real ( kind = 8 ) Y(N), the data values to be interpolated.
!
!    Input, integer ( kind = 4 ) IBCBEG, the left boundary condition flag:
!
!      0: the spline should be a quadratic over the first interval;
!      1: the first derivative at the left endpoint should be YBCBEG;
!      2: the second derivative at the left endpoint should be YBCBEG.
!
!    Input, real ( kind = 8 ) YBCBEG, the left boundary value, if needed.
!
!    Input, integer ( kind = 4 ) IBCEND, the right boundary condition flag:
!
!      0: the spline should be a quadratic over the last interval;
!      1: the first derivative at the right endpoint should be YBCEND;
!      2: the second derivative at the right endpoint should be YBCEND.
!
!    Input, real ( kind = 8 ) YBCEND, the right boundary value, if needed.
!
!    Output, real ( kind = 8 ) YPP(N), the second derivatives of 
!    the cubic spline.
!
  implicit none

  integer ( kind = 4 ) n

  real    ( kind = 8 ) a(3,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ibcbeg
  integer ( kind = 4 ) ibcend
  real    ( kind = 8 ) t(n)
  real    ( kind = 8 ) y(n)
  real    ( kind = 8 ) ybcbeg
  real    ( kind = 8 ) ybcend
  real    ( kind = 8 ) ypp(n)
!
!  Check.
!
  if ( n <= 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'SPLINE_CUBIC_SET - Fatal error!'
    write ( *, '(a)' ) '  The number of knots must be at least 2.'
    write ( *, '(a,i8)' ) '  The input value of N = ', n
    stop
  end if

  do i = 1, n-1
    if ( t(i+1) <= t(i) ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'SPLINE_CUBIC_SET - Fatal error!'
      write ( *, '(a)' ) '  The knots must be strictly increasing, but'
      write ( *, '(a,i8,a,g14.6)' ) '  T(',  i,') = ', t(i)
      write ( *, '(a,i8,a,g14.6)' ) '  T(',i+1,') = ', t(i+1)
      stop
    end if
  end do
!
!  Set the first equation.
!
  if ( ibcbeg == 0 ) then
    ypp(1) = 0.0D+00
    a(2,1) = 1.0D+00
    a(1,2) = -1.0D+00
  else if ( ibcbeg == 1 ) then
    ypp(1) = ( y(2) - y(1) ) / ( t(2) - t(1) ) - ybcbeg
    a(2,1) = ( t(2) - t(1) ) / 3.0D+00 
    a(1,2) = ( t(2) - t(1) ) / 6.0D+00
  else if ( ibcbeg == 2 ) then
    ypp(1) = ybcbeg
    a(2,1) = 1.0D+00
    a(1,2) = 0.0D+00
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'SPLINE_CUBIC_SET - Fatal error!'
    write ( *, '(a)' ) '  The boundary flag IBCBEG must be 0, 1 or 2.'
    write ( *, '(a,i8)' ) '  The input value is IBCBEG = ', ibcbeg
    stop
  end if
!
!  Set the intermediate equations.
!
  do i = 2, n-1
    ypp(i) = ( y(i+1) - y(i) ) / ( t(i+1) - t(i) ) &
           - ( y(i) - y(i-1) ) / ( t(i) - t(i-1) )
    a(3,i-1) = ( t(i) - t(i-1) ) / 6.0D+00
    a(2,i) = ( t(i+1) - t(i-1) ) / 3.0D+00
    a(1,i+1) = ( t(i+1) - t(i) ) / 6.0D+00
  end do
!
!  Set the last equation.
!
  if ( ibcend == 0 ) then
    ypp(n) = 0.0D+00
    a(3,n-1) = -1.0D+00
    a(2,n) = 1.0D+00
  else if ( ibcend == 1 ) then
    ypp(n) = ybcend - ( y(n) - y(n-1) ) / ( t(n) - t(n-1) )
    a(3,n-1) = ( t(n) - t(n-1) ) / 6.0D+00
    a(2,n) = ( t(n) - t(n-1) ) / 3.0D+00
  else if ( ibcend == 2 ) then
    ypp(n) = ybcend
    a(3,n-1) = 0.0D+00
    a(2,n) = 1.0D+00
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'SPLINE_CUBIC_SET - Fatal error!'
    write ( *, '(a)' ) '  The boundary flag IBCEND must be 0, 1 or 2.'
    write ( *, '(a,i8)' ) '  The input value is IBCEND = ', ibcend
    stop
  end if
!
!  Special case:
!    N = 2, IBCBEG = IBCEND = 0.
!
  if ( n == 2 .and. ibcbeg == 0 .and. ibcend == 0 ) then

    ypp(1) = 0.0D+00
    ypp(2) = 0.0D+00
!
!  Solve the linear system.
!
  else

    call r83_np_fs ( n, a, ypp, ypp )

  end if

  return
end
subroutine spline_cubic_val ( n, t, y, ypp, tval, yval, ypval, yppval )

!*****************************************************************************80
!
!! SPLINE_CUBIC_VAL evaluates a piecewise cubic spline at a point.
!
!  Discussion:
!
!    SPLINE_CUBIC_SET must have already been called to define the 
!    values of YPP.
!
!    For any point T in the interval T(IVAL), T(IVAL+1), the form of
!    the spline is
!
!      SPL(T) = A 
!             + B * ( T - T(IVAL) ) 
!             + C * ( T - T(IVAL) )**2
!             + D * ( T - T(IVAL) )**3
!
!    Here:
!      A = Y(IVAL)
!      B = ( Y(IVAL+1) - Y(IVAL) ) / ( T(IVAL+1) - T(IVAL) )
!        - ( YPP(IVAL+1) + 2 * YPP(IVAL) ) * ( T(IVAL+1) - T(IVAL) ) / 6
!      C = YPP(IVAL) / 2
!      D = ( YPP(IVAL+1) - YPP(IVAL) ) / ( 6 * ( T(IVAL+1) - T(IVAL) ) )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    20 November 2000
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Carl deBoor,
!    A Practical Guide to Splines,
!    Springer, 2001,
!    ISBN: 0387953663.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of data values.
!
!    Input, real ( kind = 8 ) T(N), the knot values.
!
!    Input, real ( kind = 8 ) Y(N), the data values at the knots.
!
!    Input, real ( kind = 8 ) YPP(N), the second derivatives of the 
!    spline at the knots.
!
!    Input, real ( kind = 8 ) TVAL, a point, typically between T(1) and 
!    T(N), at which the spline is to be evalulated.  If TVAL lies outside 
!    this range, extrapolation is used.
!
!    Output, real ( kind = 8 ) YVAL, YPVAL, YPPVAL, the value of the spline, and
!    its first two derivatives at TVAL.
!
  implicit none

  integer ( kind = 4 ) n

  real    ( kind = 8 ) dt
  real    ( kind = 8 ) h
  integer ( kind = 4 ) left
  integer ( kind = 4 ) right
  real    ( kind = 8 ) t(n)
  real    ( kind = 8 ) tval
  real    ( kind = 8 ) y(n)
  real    ( kind = 8 ) ypp(n)
  real    ( kind = 8 ) yppval
  real    ( kind = 8 ) ypval
  real    ( kind = 8 ) yval
!
!  Determine the interval [T(LEFT), T(RIGHT)] that contains TVAL.
!  Values below T(1) or above T(N) use extrapolation.
!
  call r8vec_bracket ( n, t, tval, left, right )
!
!  Evaluate the polynomial.
!
  dt = tval - t(left)
  h = t(right) - t(left)

  yval = y(left) &
       + dt * ( ( y(right) - y(left) ) / h &
              - ( ypp(right) / 6.0D+00 + ypp(left) / 3.0D+00 ) * h &
       + dt * ( 0.5D+00 * ypp(left) &
       + dt * ( ( ypp(right) - ypp(left) ) / ( 6.0D+00 * h ) ) ) )

  ypval = ( y(right) - y(left) ) / h &
       - ( ypp(right) / 6.0D+00 + ypp(left) / 3.0D+00 ) * h &
       + dt * ( ypp(left) &
       + dt * ( 0.5D+00 * ( ypp(right) - ypp(left) ) / h ) )

  yppval = ypp(left) + dt * ( ypp(right) - ypp(left) ) / h 

  return
end
subroutine spline_cubic_val2 ( n, t, y, ypp, left, tval, yval, ypval, yppval )

!*****************************************************************************80
!
!! SPLINE_CUBIC_VAL2 evaluates a piecewise cubic spline at a point.
!
!  Discussion:
!
!    This routine is a modification of SPLINE_CUBIC_VAL; it allows the
!    user to speed up the code by suggesting the appropriate T interval
!    to search first.
!
!    SPLINE_CUBIC_SET must have already been called to define the
!    values of YPP.
!
!    In the LEFT interval, let RIGHT = LEFT+1.  The form of the spline is
!
!      SPL(T) = 
!          A
!        + B * ( T - T(LEFT) )
!        + C * ( T - T(LEFT) )**2
!        + D * ( T - T(LEFT) )**3
!
!    Here:
!      A = Y(LEFT)
!      B = ( Y(RIGHT) - Y(LEFT) ) / ( T(RIGHT) - T(LEFT) )
!        - ( YPP(RIGHT) + 2 * YPP(LEFT) ) * ( T(RIGHT) - T(LEFT) ) / 6
!      C = YPP(LEFT) / 2
!      D = ( YPP(RIGHT) - YPP(LEFT) ) / ( 6 * ( T(RIGHT) - T(LEFT) ) )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    24 February 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Carl deBoor,
!    A Practical Guide to Splines,
!    Springer, 2001,
!    ISBN: 0387953663.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of knots.
!
!    Input, real ( kind = 8 ) T(N), the knot values.
!
!    Input, real ( kind = 8 ) Y(N), the data values at the knots.
!
!    Input, real ( kind = 8 ) YPP(N), the second derivatives of the spline at
!    the knots.
!
!    Input/output, integer ( kind = 4 ) LEFT, the suggested T interval to search.
!    LEFT should be between 1 and N-1.  If LEFT is not in this range,
!    then its value will be ignored.  On output, LEFT is set to the
!    actual interval in which TVAL lies.
!
!    Input, real ( kind = 8 ) TVAL, a point, typically between T(1) and T(N), at
!    which the spline is to be evalulated.  If TVAL lies outside
!    this range, extrapolation is used.
!
!    Output, real ( kind = 8 ) YVAL, YPVAL, YPPVAL, the value of the spline, and
!    its first two derivatives at TVAL.
!
  implicit none

  integer ( kind = 4 ) n

  real    ( kind = 8 ) dt
  real    ( kind = 8 ) h
  integer ( kind = 4 ) left
  integer ( kind = 4 ) right
  real    ( kind = 8 ) t(n)
  real    ( kind = 8 ) tval
  real    ( kind = 8 ) y(n)
  real    ( kind = 8 ) ypp(n)
  real    ( kind = 8 ) yppval
  real    ( kind = 8 ) ypval
  real    ( kind = 8 ) yval
!
!  Determine the interval [T(LEFT), T(RIGHT)] that contains TVAL.
!  
!  What you want from R8VEC_BRACKET3 is that TVAL is to be computed
!  by the data in interval {T(LEFT), T(RIGHT)].
!
  call r8vec_bracket3 ( n, t, tval, left )
  right = left + 1
!
!  In the interval LEFT, the polynomial is in terms of a normalized
!  coordinate  ( DT / H ) between 0 and 1.
!
  dt = tval - t(left)
  h = t(right) - t(left)

  yval = y(left) + dt * ( ( y(right) - y(left) ) / h &
              - ( ypp(right) / 6.0D+00 + ypp(left) / 3.0D+00 ) * h &
       + dt * ( 0.5D+00 * ypp(left) &
       + dt * ( ( ypp(right) - ypp(left) ) / ( 6.0D+00 * h ) ) ) )

  ypval = ( y(right) - y(left) ) / h &
      - ( ypp(right) / 6.0D+00 + ypp(left) / 3.0D+00 ) * h &
      + dt * ( ypp(left) &
      + dt * ( 0.5D+00 * ( ypp(right) - ypp(left) ) / h ) )

  yppval = ypp(left) + dt * ( ypp(right) - ypp(left) ) / h

  return
end
subroutine spline_hermite_set ( ndata, tdata, ydata, ypdata, c )

!*****************************************************************************80
!
!! SPLINE_HERMITE_SET sets up a piecewise cubic Hermite interpolant.
!
!  Discussion:
!
!    Once the array C is computed, then in the interval 
!    (TDATA(I), TDATA(I+1)), the interpolating Hermite polynomial 
!    is given by
!
!      SVAL(TVAL) =                 C(1,I)
!         + ( TVAL - TDATA(I) ) * ( C(2,I)
!         + ( TVAL - TDATA(I) ) * ( C(3,I)
!         + ( TVAL - TDATA(I) ) *   C(4,I) ) )
!
!    This is algorithm CALCCF of Conte and deBoor.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    11 February 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Samuel Conte, Carl deBoor,
!    Elementary Numerical Analysis,
!    Second Edition,
!    McGraw Hill, 1972,
!    ISBN: 07-012446-4,
!    LC: QA297.C65.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NDATA, the number of data points.
!    NDATA must be at least 2.
!
!    Input, real ( kind = 8 ) TDATA(NDATA), the abscissas of the data points.
!    The entries of TDATA are assumed to be strictly increasing.
!
!    Input, real ( kind = 8 ) Y(NDATA), YP(NDATA), the value of the
!    function and its derivative at TDATA(1:NDATA).
!
!    Output, real ( kind = 8 ) C(4,NDATA), the coefficients of the 
!    Hermite polynomial.
!    C(1,1:NDATA) = Y(1:NDATA) and C(2,1:NDATA) = YP(1:NDATA).
!    C(3,1:NDATA-1) and C(4,1:NDATA-1) are the quadratic and cubic 
!    coefficients.
!
  implicit none

  integer ( kind = 4 ) ndata

  real    ( kind = 8 ) c(4,ndata)
  real    ( kind = 8 ) divdif1
  real    ( kind = 8 ) divdif3
  real    ( kind = 8 ) dt
  integer ( kind = 4 ) i
  real    ( kind = 8 ) tdata(ndata)
  real    ( kind = 8 ) ydata(ndata)
  real    ( kind = 8 ) ypdata(ndata)
!
!  Check NDATA.
!
  if ( ndata < 2 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'SPLINE_HERMITE_SET - Fatal error!'
    write ( *, '(a)' ) '  NDATA < 2.'
    stop
  end if

  c(1,1:ndata) = ydata(1:ndata)
  c(2,1:ndata) = ypdata(1:ndata)

  do i = 1, ndata-1
    dt = tdata(i+1) - tdata(i)
    divdif1 = ( c(1,i+1) - c(1,i) ) / dt
    divdif3 = c(2,i) + c(2,i+1) - 2.0D+00 * divdif1
    c(3,i) = ( divdif1 - c(2,i) - divdif3 ) / dt
    c(4,i) = divdif3 / ( dt * dt )
  end do

  c(3,ndata) = 0.0D+00
  c(4,ndata) = 0.0D+00

  return
end
subroutine spline_hermite_val ( ndata, tdata, c, tval, sval, spval )

!*****************************************************************************80
!
!! SPLINE_HERMITE_VAL evaluates a piecewise cubic Hermite interpolant.
!
!  Discussion:
!
!    SPLINE_HERMITE_SET must be called first, to set up the
!    spline data from the raw function and derivative data.
!
!    In the interval (TDATA(I), TDATA(I+1)), the interpolating 
!    Hermite polynomial is given by
!
!      SVAL(TVAL) =                 C(1,I)
!         + ( TVAL - TDATA(I) ) * ( C(2,I)
!         + ( TVAL - TDATA(I) ) * ( C(3,I)
!         + ( TVAL - TDATA(I) ) *   C(4,I) ) )
!
!    and
!
!      SVAL'(TVAL) =                    C(2,I)
!         + ( TVAL - TDATA(I) ) * ( 2 * C(3,I)
!         + ( TVAL - TDATA(I) ) *   3 * C(4,I) )
!
!    This is algorithm PCUBIC of Conte and deBoor.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    11 February 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Samuel Conte, Carl deBoor,
!    Elementary Numerical Analysis,
!    Second Edition,
!    McGraw Hill, 1972,
!    ISBN: 07-012446-4,
!    LC: QA297.C65.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NDATA, the number of data points.
!    NDATA must be at least 2.
!
!    Input, real ( kind = 8 ) TDATA(NDATA), the abscissas of the data points.
!    The entries of TDATA are assumed to be strictly increasing.
!
!    Input, real ( kind = 8 ) C(4,NDATA), the coefficient data computed by
!    SPLINE_HERMITE_SET.
!
!    Input, real ( kind = 8 ) TVAL, the point where the interpolant is to
!    be evaluated.
!
!    Output, real ( kind = 8 ) SVAL, SPVAL, the value of the interpolant 
!    and its derivative at TVAL.
!
  implicit none

  integer ( kind = 4 ) ndata

  real    ( kind = 8 ) c(4,ndata)
  real    ( kind = 8 ) dt
  integer ( kind = 4 ) left
  integer ( kind = 4 ) right
  real    ( kind = 8 ) spval
  real    ( kind = 8 ) sval
  real    ( kind = 8 ) tdata(ndata)
  real    ( kind = 8 ) tval
!
!  Check NDATA.
!
  if ( ndata < 2 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'SPLINE_HERMITE_VAL - Fatal error!'
    write ( *, '(a)' ) '  NDATA < 2.'
    stop
  end if
!
!  Find the interval [ TDATA(LEFT), TDATA(RIGHT) ] that contains
!  or is nearest to TVAL.
!
  call r8vec_bracket ( ndata, tdata, tval, left, right )
!
!  Evaluate the cubic polynomial.
!
  dt = tval - tdata(left)

  sval = c(1,left) + dt * ( c(2,left) + dt * ( c(3,left) + dt * c(4,left) ) )

  spval = c(2,left) + dt * ( 2.0D+00 * c(3,left) + dt * 3.0D+00 * c(4,left) )

  return
end
subroutine spline_linear_int ( ndata, tdata, ydata, a, b, int_val )

!*****************************************************************************80
!
!! SPLINE_LINEAR_INT evaluates the integral of a piecewise linear spline.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    01 November 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NDATA, the number of data points defining 
!    the spline.  NDATA must be at least 2.
!
!    Input, real ( kind = 8 ) TDATA(NDATA), YDATA(NDATA), the values of 
!    the independent and dependent variables at the data points.  The
!    values of TDATA should be distinct and increasing.
!
!    Input, real ( kind = 8 ) A, B, the interval over which the 
!    integral is desired.
!
!    Output, real ( kind = 8 ) INT_VAL, the value of the integral.
!
  implicit none

  integer ( kind = 4 ) ndata

  real    ( kind = 8 ) a
  real    ( kind = 8 ) a_copy
  integer ( kind = 4 ) a_left
  integer ( kind = 4 ) a_right
  real    ( kind = 8 ) b
  real    ( kind = 8 ) b_copy
  integer ( kind = 4 ) b_left
  integer ( kind = 4 ) b_right
  integer ( kind = 4 ) i_left
  real    ( kind = 8 ) int_val
  real    ( kind = 8 ) tdata(ndata)
  real    ( kind = 8 ) tval
  real    ( kind = 8 ) ydata(ndata)
  real    ( kind = 8 ) yp
  real    ( kind = 8 ) yval
!
!  Check NDATA.
!
  if ( ndata < 2 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'SPLINE_LINEAR_INT - Fatal error!'
    write ( *, '(a)' ) '  NDATA < 2.'
    stop
  end if

  int_val = 0.0D+00 

  if ( a == b ) then
    return
  end if

  a_copy = min ( a, b )
  b_copy = max ( a, b )
!
!  Find the interval [ TDATA(A_LEFT), TDATA(A_RIGHT) ] that contains, or is
!  nearest to, A.
!
  call r8vec_bracket ( ndata, tdata, a_copy, a_left, a_right )
!
!  Find the interval [ TDATA(B_LEFT), TDATA(B_RIGHT) ] that contains, or is
!  nearest to, B.
!
  call r8vec_bracket ( ndata, tdata, b_copy, b_left, b_right )
!
!  If A and B are in the same interval...
!
  if ( a_left == b_left ) then

    tval = ( a_copy + b_copy ) / 2.0D+00

    yp = ( ydata(a_right) - ydata(a_left) ) / &
         ( tdata(a_right) - tdata(a_left) )

    yval = ydata(a_left) + ( tval - tdata(a_left) ) * yp

    int_val = yval * ( b_copy - a_copy )

    return
  end if
!
!  Otherwise, integrate from:
!
!  A               to TDATA(A_RIGHT),
!  TDATA(A_RIGHT)  to TDATA(A_RIGHT+1),...
!  TDATA(B_LEFT-1) to TDATA(B_LEFT),
!  TDATA(B_LEFT)   to B.
!
!  Use the fact that the integral of a linear function is the
!  value of the function at the midpoint times the width of the interval.
!
  tval = ( a_copy + tdata(a_right) ) / 2.0D+00

  yp = ( ydata(a_right) - ydata(a_left) ) / &
       ( tdata(a_right) - tdata(a_left) )

  yval = ydata(a_left) + ( tval - tdata(a_left) ) * yp

  int_val = int_val + yval * ( tdata(a_right) - a_copy )

  do i_left = a_right, b_left - 1

    tval = ( tdata(i_left+1) + tdata(i_left) ) / 2.0D+00

    yp = ( ydata(i_left+1) - ydata(i_left) ) / &
         ( tdata(i_left+1) - tdata(i_left) )

    yval = ydata(i_left) + ( tval - tdata(i_left) ) * yp

    int_val = int_val + yval * ( tdata(i_left + 1) - tdata(i_left) )

  end do

  tval = ( tdata(b_left) + b_copy ) / 2.0D+00

  yp = ( ydata(b_right) - ydata(b_left) ) / &
       ( tdata(b_right) - tdata(b_left) )

  yval = ydata(b_left) + ( tval - tdata(b_left) ) * yp

  int_val = int_val + yval * ( b_copy - tdata(b_left) )

  if ( b < a ) then
    int_val = - int_val
  end if

  return
end
subroutine spline_linear_intset ( n, int_x, int_v, data_x, data_y )

!*****************************************************************************80
!
!! SPLINE_LINEAR_INTSET sets a piecewise linear spline with given integral properties.
!
!  Discussion:
!
!    The user has in mind an interval, divided by INT_N+1 points into
!    INT_N intervals.  A linear spline is to be constructed,
!    with breakpoints at the centers of each interval, and extending
!    continuously to the left of the first and right of the last
!    breakpoints.  The constraint on the linear spline is that it is
!    required that it have a given integral value over each interval.
!
!    A tridiagonal linear system of equations is solved for the
!    values of the spline at the breakpoints.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    27 January 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of intervals.
!
!    Input, real ( kind = 8 ) INT_X(N+1), the points that define the intervals.
!    Interval I lies between INT_X(I) and INT_X(I+1).
!
!    Input, real ( kind = 8 ) INT_V(N), the desired value of the integral of the
!    linear spline over each interval.
!
!    Output, real ( kind = 8 ) DATA_X(N), DATA_Y(N), the values of the
!    independent and dependent variables at the data points.  The values 
!    of DATA_X are the interval midpoints.  The values of DATA_Y are 
!    determined in such a way that the exact integral of the linear 
!    spline over interval I is equal to INT_V(I).
!
  implicit none

  integer ( kind = 4 ) n

  real    ( kind = 8 ) a(3,n)
  real    ( kind = 8 ) data_x(n)
  real    ( kind = 8 ) data_y(n)
  real    ( kind = 8 ) int_v(n)
  real    ( kind = 8 ) int_x(n+1)
!
!  Set up the easy stuff.
!
  data_x(1:n) = 0.5D+00 * ( int_x(1:n) + int_x(2:n+1) )
!
!  Set up C, D, E, the coefficients of the linear system.
!
  a(3,1:n-2) = 1.0D+00 &
    - ( 0.5D+00 * ( data_x(2:n-1) + int_x(2:n-1) ) - data_x(1:n-2) ) &
    / ( data_x(2:n-1) - data_x(1:n-2) )
  a(3,n-1) = 0.0D+00
  a(3,n) = 0.0D+00

  a(2,1) = int_x(2) - int_x(1)

  a(2,2:n-1) = 1.0D+00 &
    + ( 0.5D+00 * ( data_x(2:n-1) + int_x(2:n-1) ) &
    - data_x(1:n-2) ) &
    / ( data_x(2:n-1) - data_x(1:n-2) ) &
    - ( 0.5D+00 * ( data_x(2:n-1) + int_x(3:n) ) - data_x(2:n-1) ) &
    / ( data_x(3:n) - data_x(2:n-1) )

  a(2,n) = int_x(n+1) - int_x(n)

  a(1,1) = 0.0D+00
  a(1,2) = 0.0D+00

  a(1,3:n) = ( 0.5D+00 * ( data_x(2:n-1) + int_x(3:n) ) &
    - data_x(2:n-1) ) / ( data_x(3:n) - data_x(2:n-1) )
!
!  Set up DATA_Y, which begins as the right hand side of the linear system.
!
  data_y(1) = int_v(1)
  data_y(2:n-1) = 2.0D+00 * int_v(2:n-1) / ( int_x(3:n) - int_x(2:n-1) )
  data_y(n) = int_v(n)
!
!  Solve the linear system.
!
  call r83_np_fs ( n, a, data_y, data_y )

  return
end
subroutine spline_linear_val ( ndata, tdata, ydata, tval, yval, ypval )

!*****************************************************************************80
!
!! SPLINE_LINEAR_VAL evaluates a piecewise linear spline at a point.
!
!  Discussion:
!
!    Because of the extremely simple form of the linear spline,
!    the raw data points ( TDATA(I), YDATA(I)) can be used directly to
!    evaluate the spline at any point.  No processing of the data
!    is required.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    06 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NDATA, the number of data points defining 
!    the spline.  NDATA must be at least 2.
!
!    Input, real ( kind = 8 ) TDATA(NDATA), YDATA(NDATA), the values of 
!    the independent and dependent variables at the data points.  The 
!    values of TDATA should be distinct and increasing.
!
!    Input, real ( kind = 8 ) TVAL, the point at which the spline is
!    to be evaluated.
!
!    Output, real ( kind = 8 ) YVAL, YPVAL, the value of the spline and 
!    its first derivative dYdT at TVAL.  YPVAL is not reliable if TVAL 
!    is exactly equal to TDATA(I) for some I.
!
  implicit none

  integer ( kind = 4 ) ndata

  integer ( kind = 4 ) left
  integer ( kind = 4 ) right
  real    ( kind = 8 ) tdata(ndata)
  real    ( kind = 8 ) tval
  real    ( kind = 8 ) ydata(ndata)
  real    ( kind = 8 ) ypval
  real    ( kind = 8 ) yval
!
!  Check NDATA.
!
  if ( ndata < 2 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'SPLINE_LINEAR_VAL - Fatal error!'
    write ( *, '(a)' ) '  NDATA < 2.'
    stop
  end if
!
!  Find the interval [ TDATA(LEFT), TDATA(RIGHT) ] that contains, or is
!  nearest to, TVAL.
!
  call r8vec_bracket ( ndata, tdata, tval, left, right )
!
!  Now evaluate the piecewise linear function.
!
  ypval = ( ydata(right) - ydata(left) ) / ( tdata(right) - tdata(left) )

  yval = ydata(left) +  ( tval - tdata(left) ) * ypval

  return
end
subroutine spline_overhauser_nonuni_val ( ndata, tdata, ydata, tval, yval )

!*****************************************************************************80
!
!! SPLINE_OVERHAUSER_NONUNI_VAL evaluates the nonuniform Overhauser spline.
!
!  Discussion:
!
!    The nonuniformity refers to the fact that the abscissa values
!    need not be uniformly spaced.
!
!    Thanks to Doug Fortune for pointing out that the point distances
!    used to define ALPHA and BETA should be the Euclidean distances
!    between the points.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    08 May 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NDATA, the number of data points.
!    3 <= NDATA is required.
!
!    Input, real ( kind = 8 ) TDATA(NDATA), the abscissas of the data points.
!    The values of TDATA are assumed to be distinct and increasing.
!
!    Input, real ( kind = 8 ) YDATA(NDATA), the data values.
!
!    Input, real ( kind = 8 ) TVAL, the value where the spline is to
!    be evaluated.
!
!    Output, real ( kind = 8 ) YVAL, the value of the spline at TVAL.
!
  implicit none

  integer ( kind = 4 ) ndata

  real    ( kind = 8 ) alpha
  real    ( kind = 8 ) beta
  real    ( kind = 8 ) d21
  real    ( kind = 8 ) d32
  real    ( kind = 8 ) d43
  integer ( kind = 4 ) left
  real    ( kind = 8 ) mbasis(4,4)
  real    ( kind = 8 ) mbasis_l(3,3)
  real    ( kind = 8 ) mbasis_r(3,3)
  integer ( kind = 4 ) right
  real    ( kind = 8 ) tdata(ndata)
  real    ( kind = 8 ) tval
  real    ( kind = 8 ) ydata(ndata)
  real    ( kind = 8 ) yval
!
!  Check NDATA.
!
  if ( ndata < 3 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'SPLINE_OVERHAUSER_NONUNI_VAL - Fatal error!'
    write ( *, '(a)' ) '  NDATA < 3.'
    stop
  end if
!
!  Find the nearest interval [ TDATA(LEFT), TDATA(RIGHT) ] to TVAL.
!
  call r8vec_bracket ( ndata, tdata, tval, left, right )
!
!  Evaluate the spline in the given interval.
!
  if ( left == 1 ) then

    d21 = sqrt ( ( tdata(2) - tdata(1) )**2 &
               + ( ydata(2) - ydata(1) )**2 )

    d32 = sqrt ( ( tdata(3) - tdata(2) )**2 &
               + ( ydata(3) - ydata(2) )**2 )

    alpha = d21 / ( d32 + d21 )

    call basis_matrix_overhauser_nul ( alpha, mbasis_l )

    call basis_matrix_tmp ( left, 3, mbasis_l, ndata, tdata, ydata, tval, yval )

  else if ( left < ndata-1 ) then

    d21 = sqrt ( ( tdata(left) - tdata(left-1) )**2 &
               + ( ydata(left) - ydata(left-1) )**2 )

    d32 = sqrt ( ( tdata(left+1) - tdata(left) )**2 &
               + ( ydata(left+1) - ydata(left) )**2 )

    d43 = sqrt ( ( tdata(left+2) - tdata(left+1) )**2 &
               + ( ydata(left+2) - ydata(left+1) )**2 )

    alpha = d21 / ( d32 + d21 )
    beta  = d32 / ( d43 + d32 )

    call basis_matrix_overhauser_nonuni ( alpha, beta, mbasis )

    call basis_matrix_tmp ( left, 4, mbasis, ndata, tdata, ydata, tval, yval )

  else if ( left == ndata-1 ) then

    d32 = sqrt ( ( tdata(ndata-1) - tdata(ndata-2) )**2 &
               + ( ydata(ndata-1) - ydata(ndata-2) )**2 )

    d43 = sqrt ( ( tdata(ndata) - tdata(ndata-1) )**2 &
               + ( ydata(ndata) - ydata(ndata-1) )**2 )

    beta  = d32 / ( d43 + d32 )

    call basis_matrix_overhauser_nur ( beta, mbasis_r )

    call basis_matrix_tmp ( left, 3, mbasis_r, ndata, tdata, ydata, tval, yval )

  else

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'SPLINE_OVERHAUSER_NONUNI_VAL - Fatal error!'
    write ( *, '(a,i8)' ) '  Nonsensical value of LEFT = ', left
    write ( *, '(a,i8)' ) '  but 0 < LEFT < NDATA = ', ndata
    write ( *, '(a)' ) '  is required.'
    stop

  end if

  return
end
subroutine spline_overhauser_uni_val ( ndata, tdata, ydata, tval, yval )

!*****************************************************************************80
!
!! SPLINE_OVERHAUSER_UNI_VAL evaluates the uniform Overhauser spline.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    06 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NDATA, the number of data points.
!    NDATA must be at least 3.
!
!    Input, real ( kind = 8 ) TDATA(NDATA), the abscissas of the data points.
!    The values of TDATA are assumed to be distinct and increasing.
!    This routine also assumes that the values of TDATA are uniformly
!    spaced; for instance, TDATA(1) = 10, TDATA(2) = 11, TDATA(3) = 12...
!
!    Input, real ( kind = 8 ) YDATA(NDATA), the data values.
!
!    Input, real ( kind = 8 ) TVAL, the value where the spline is to
!    be evaluated.
!
!    Output, real ( kind = 8 ) YVAL, the value of the spline at TVAL.
!
  implicit none

  integer ( kind = 4 ) ndata

  integer ( kind = 4 ) left
  real    ( kind = 8 ) mbasis(4,4)
  real    ( kind = 8 ) mbasis_l(3,3)
  real    ( kind = 8 ) mbasis_r(3,3)
  integer ( kind = 4 ) right
  real    ( kind = 8 ) tdata(ndata)
  real    ( kind = 8 ) tval
  real    ( kind = 8 ) ydata(ndata)
  real    ( kind = 8 ) yval
!
!  Check NDATA.
!
  if ( ndata < 3 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'SPLINE_OVERHAUSER_UNI_VAL - Fatal error!'
    write ( *, '(a)' ) '  NDATA < 3.'
    stop
  end if
!
!  Find the nearest interval [ TDATA(LEFT), TDATA(RIGHT) ] to TVAL.
!
  call r8vec_bracket ( ndata, tdata, tval, left, right )
!
!  Evaluate the spline in the given interval.
!
  if ( left == 1 ) then

    call basis_matrix_overhauser_uni_l ( mbasis_l )

    call basis_matrix_tmp ( left, 3, mbasis_l, ndata, tdata, ydata, tval, yval )

  else if ( left < ndata-1 ) then

    call basis_matrix_overhauser_uni ( mbasis )

    call basis_matrix_tmp ( left, 4, mbasis, ndata, tdata, ydata, tval, yval )

  else if ( left == ndata-1 ) then

    call basis_matrix_overhauser_uni_r ( mbasis_r )

    call basis_matrix_tmp ( left, 3, mbasis_r, ndata, tdata, ydata, tval, yval )

  end if

  return
end
subroutine spline_overhauser_val ( dim_num, ndata, tdata, ydata, tval, yval )

!*****************************************************************************80
!
!! SPLINE_OVERHAUSER_VAL evaluates an Overhauser spline.
!
!  Discussion:
!
!    Over the first and last intervals, the Overhauser spline is a 
!    quadratic.  In the intermediate intervals, it is a piecewise cubic.
!    The Overhauser spline is also known as the Catmull-Rom spline.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!   08 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    JA Brewer, DC Anderson,
!    Visual Interaction with Overhauser Curves and Surfaces,
!    SIGGRAPH 77,
!    in Proceedings of the 4th Annual Conference on Computer Graphics
!    and Interactive Techniques,
!    ASME, July 1977, pages 132-137.
!
!    Edwin Catmull, Raphael Rom,
!    A Class of Local Interpolating Splines,
!    in Computer Aided Geometric Design,
!    edited by Robert Barnhill, Richard Reisenfeld,
!    Academic Press, 1974, pages 317-326,
!    ISBN: 0120790505.
!
!    David Rogers, Alan Adams,
!    Mathematical Elements of Computer Graphics,
!    Second Edition,
!    McGraw Hill, 1989,
!    ISBN: 0070535299.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the dimension of a single data point.
!    DIM_NUM must be at least 1.
!
!    Input, integer ( kind = 4 ) NDATA, the number of data points.
!    NDATA must be at least 3.
!
!    Input, real ( kind = 8 ) TDATA(NDATA), the abscissas of the data 
!    points.  The values in TDATA must be in strictly ascending order.
!
!    Input, real ( kind = 8 ) YDATA(DIM_NUM,NDATA), the data points 
!    corresponding to the abscissas.
!
!    Input, real ( kind = 8 ) TVAL, the abscissa value at which the spline
!    is to be evaluated.  Normally, TDATA(1) <= TVAL <= T(NDATA), and 
!    the data will be interpolated.  For TVAL outside this range, 
!    extrapolation will be used.
!
!    Output, real ( kind = 8 ) YVAL(DIM_NUM), the value of the spline at TVAL.
!
  implicit none

  integer ( kind = 4 ) ndata
  integer ( kind = 4 ) dim_num

  integer ( kind = 4 ) left
  integer ( kind = 4 ) order
  integer ( kind = 4 ) right
  real    ( kind = 8 ) tdata(ndata)
  real    ( kind = 8 ) tval
  real    ( kind = 8 ) ydata(dim_num,ndata)
  real    ( kind = 8 ) yl(dim_num)
  real    ( kind = 8 ) yr(dim_num)
  real    ( kind = 8 ) yval(dim_num)
!
!  Check.
!
  call r8vec_order_type ( ndata, tdata, order )

  if ( order /= 2 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'SPLINE_OVERHAUSER_VAL - Fatal error!'
    write ( *, '(a)' ) '  The data abscissas are not strictly ascending.'
    stop
  end if

  if ( ndata < 3 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'SPLINE_OVERHAUSER_VAL - Fatal error!'
    write ( *, '(a)' ) '  NDATA < 3.'
    stop
  end if

  if ( dim_num < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'SPLINE_OVERHAUSER_VAL - Fatal error!'
    write ( *, '(a)' ) '  DIM_NUM < 1.'
    stop
  end if
!
!  Locate the abscissa interval T(LEFT), T(LEFT+1) nearest to or 
!  containing TVAL.
!
  call r8vec_bracket ( ndata, tdata, tval, left, right )
!
!  Evaluate the "left hand" quadratic defined at T(LEFT-1), T(LEFT), T(RIGHT).
!
  if ( 0 < left-1 ) then
    call parabola_val2 ( dim_num, ndata, tdata, ydata, left-1, tval, yl )
  end if
!
!  Evaluate the "right hand" quadratic defined at T(LEFT), T(RIGHT), T(RIGHT+1).
!
  if ( right+1 <= ndata ) then
    call parabola_val2 ( dim_num, ndata, tdata, ydata, left, tval, yr )
  end if
!
!  Average the quadratics.
!
  if ( left == 1 ) then

    yval(1:dim_num) = yr(1:dim_num)

  else if ( right < ndata ) then

    yval(1:dim_num) =  &
      ( ( tdata(right) - tval               ) * yl(1:dim_num)   &
      + (                tval - tdata(left) ) * yr(1:dim_num) ) &
      / ( tdata(right)        - tdata(left) )

  else

    yval(1:dim_num) = yl(1:dim_num)

  end if

  return
end
subroutine spline_pchip_set ( n, x, f, d )

!*****************************************************************************80
!
!! SPLINE_PCHIP_SET sets derivatives for a piecewise cubic Hermite interpolant.
!
!  Discussion:
!
!    This routine computes what would normally be called a Hermite 
!    interpolant.  However, the user is only required to supply function
!    values, not derivative values as well.  This routine computes
!    "suitable" derivative values, so that the resulting Hermite interpolant
!    has desirable shape and monotonicity properties.
!
!    The interpolant will have an extremum at each point where
!    monotonicity switches direction.
!
!    The resulting piecewise cubic Hermite function may be evaluated
!    by SPLINE_PCHIP_VAL.
!
!    This routine was originally named "PCHIM".
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    15 December 2008
!
!  Author:
!
!    Original FORTRAN77 version by Fred Fritsch.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Fred Fritsch, Ralph Carlson,
!    Monotone Piecewise Cubic Interpolation,
!    SIAM Journal on Numerical Analysis,
!    Volume 17, Number 2, April 1980, pages 238-246.
!
!    Fred Fritsch, Judy Butland,
!    A Method for Constructing Local Monotone Piecewise Cubic Interpolants,
!    SIAM Journal on Scientific and Statistical Computing,
!    Volume 5, Number 2, 1984, pages 300-304.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of data points.  N must be 
!    at least 2.
!
!    Input, real ( kind = 8 ) X(N), the strictly increasing independent
!    variable values.
!
!    Input, real ( kind = 8 ) F(N), dependent variable values to be
!    interpolated.  F(I) is the value corresponding to X(I).
!    This routine is designed for monotonic data, but it will work for any
!    F array.  It will force extrema at points where monotonicity switches
!    direction.
!
!    Output, real ( kind = 8 ) D(N), the derivative values at the
!    data points.  If the data are monotonic, these values will determine
!    a monotone cubic Hermite function.
!
  implicit none

  integer ( kind = 4 ) n

  real    ( kind = 8 ) d(n)
  real    ( kind = 8 ) del1
  real    ( kind = 8 ) del2
  real    ( kind = 8 ) dmax
  real    ( kind = 8 ) dmin
  real    ( kind = 8 ) drat1
  real    ( kind = 8 ) drat2
  real    ( kind = 8 ) dsave
  real    ( kind = 8 ) f(n)
  real    ( kind = 8 ) h1
  real    ( kind = 8 ) h2
  real    ( kind = 8 ) hsum
  real    ( kind = 8 ) hsumt3
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierr
  integer ( kind = 4 ) nless1
  real    ( kind = 8 ) pchst
  real    ( kind = 8 ) temp
  real    ( kind = 8 ) w1
  real    ( kind = 8 ) w2
  real    ( kind = 8 ) x(n)
!
!  Check the arguments.
!
  if ( n < 2 ) then
    ierr = -1
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'SPLINE_PCHIP_SET - Fatal error!'
    write ( *, '(a)' ) '  Number of data points less than 2.'
    stop
  end if

  do i = 2, n
    if ( x(i) <= x(i-1) ) then
      ierr = -3
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'SPLINE_PCHIP_SET - Fatal error!'
      write ( *, '(a)' ) '  X array not strictly increasing.'
      stop
    end if
  end do

  ierr = 0
  nless1 = n - 1
  h1 = x(2) - x(1)
  del1 = ( f(2) - f(1) ) / h1
  dsave = del1
!
!  Special case N=2, use linear interpolation.
!
  if ( n == 2 ) then
    d(1) = del1
    d(n) = del1
    return
  end if
!
!  Normal case, 3 <= N.
!
  h2 = x(3) - x(2)
  del2 = ( f(3) - f(2) ) / h2
!
!  Set D(1) via non-centered three point formula, adjusted to be
!  shape preserving.
!
  hsum = h1 + h2
  w1 = ( h1 + hsum ) / hsum
  w2 = -h1 / hsum
  d(1) = w1 * del1 + w2 * del2

  if ( pchst ( d(1), del1 ) <= 0.0D+00 ) then

    d(1) = 0.0D+00
!
!  Need do this check only if monotonicity switches.
!
  else if ( pchst ( del1, del2 ) < 0.0D+00 ) then

     dmax = 3.0D+00 * del1

     if ( abs ( dmax ) < abs ( d(1) ) ) then
       d(1) = dmax
     end if

  end if
!
!  Loop through interior points.
!
  do i = 2, nless1

    if ( 2 < i ) then
      h1 = h2
      h2 = x(i+1) - x(i)
      hsum = h1 + h2
      del1 = del2
      del2 = ( f(i+1) - f(i) ) / h2
    end if
!
!  Set D(I)=0 unless data are strictly monotonic.
!
    d(i) = 0.0D+00

    temp = pchst ( del1, del2 )

    if ( temp < 0.0D+00 ) then

      ierr = ierr + 1
      dsave = del2
!
!  Count number of changes in direction of monotonicity.
!
    else if ( temp == 0.0D+00 ) then

      if ( del2 /= 0.0D+00 ) then
        if ( pchst ( dsave, del2 ) < 0.0D+00 ) then
          ierr = ierr + 1
        end if
        dsave = del2
      end if
!
!  Use Brodlie modification of Butland formula.
!
    else

      hsumt3 = 3.0D+00 * hsum
      w1 = ( hsum + h1 ) / hsumt3
      w2 = ( hsum + h2 ) / hsumt3
      dmax = max ( abs ( del1 ), abs ( del2 ) )
      dmin = min ( abs ( del1 ), abs ( del2 ) )
      drat1 = del1 / dmax
      drat2 = del2 / dmax
      d(i) = dmin / ( w1 * drat1 + w2 * drat2 )

    end if

  end do
!
!  Set D(N) via non-centered three point formula, adjusted to be
!  shape preserving.
!
  w1 = -h2 / hsum
  w2 = ( h2 + hsum ) / hsum
  d(n) = w1 * del1 + w2 * del2

  if ( pchst ( d(n), del2 ) <= 0.0D+00 ) then
    d(n) = 0.0D+00
  else if ( pchst ( del1, del2 ) < 0.0D+00 ) then
!
!  Need do this check only if monotonicity switches.
!
    dmax = 3.0D+00 * del2

    if ( abs ( dmax ) < abs ( d(n) ) ) then
      d(n) = dmax
    end if

  end if

  return
end
subroutine spline_pchip_val ( n, x, f, d, ne, xe, fe )

!*****************************************************************************80
!
!! SPLINE_PCHIP_VAL evaluates a piecewise cubic Hermite function.
!
!  Description:
!
!    This routine may be used by itself for Hermite interpolation, or as an
!    evaluator for SPLINE_PCHIP_SET.
!
!    This routine evaluates the cubic Hermite function at the points XE.
!
!    Most of the coding between the call to CHFEV and the end of
!    the IR loop could be eliminated if it were permissible to
!    assume that XE is ordered relative to X.
!
!    CHFEV does not assume that X1 is less than X2.  Thus, it would
!    be possible to write a version of SPLINE_PCHIP_VAL that assumes a strictly
!    decreasing X array by simply running the IR loop backwards
!    and reversing the order of appropriate tests.
!
!    The present code has a minor bug, which I have decided is not
!    worth the effort that would be required to fix it.
!    If XE contains points in [X(N-1),X(N)], followed by points less than
!    X(N-1), followed by points greater than X(N), the extrapolation points
!    will be counted (at least) twice in the total returned in IERR.
!
!    The evaluation will be most efficient if the elements of XE are
!    increasing relative to X; that is, for all J <= K,
!      X(I) <= XE(J)
!    implies
!      X(I) <= XE(K).
!
!    If any of the XE are outside the interval [X(1),X(N)],
!    values are extrapolated from the nearest extreme cubic,
!    and a warning error is returned.
!
!    This routine was originally called "PCHFE".
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    14 August 2005
!
!  Author:
!
!    Original FORTRAN77 version by Fred Fritsch.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Fred Fritsch, Ralph Carlson, 
!    Monotone Piecewise Cubic Interpolation,
!    SIAM Journal on Numerical Analysis,
!    Volume 17, Number 2, April 1980, pages 238-246.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of data points.  N must be 
!    at least 2.
!
!    Input, real ( kind = 8 ) X(N), the strictly increasing independent
!    variable values.
!
!    Input, real ( kind = 8 ) F(N), the function values.
!
!    Input, real ( kind = 8 ) D(N), the derivative values.
!
!    Input, integer ( kind = 4 ) NE, the number of evaluation points.
!
!    Input, real ( kind = 8 ) XE(NE), points at which the function is to
!    be evaluated.
!
!    Output, real ( kind = 8 ) FE(NE), the values of the cubic Hermite
!    function at XE.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) ne

  real    ( kind = 8 ) d(n)
  real    ( kind = 8 ) f(n)
  real    ( kind = 8 ) fe(ne)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierc
  integer ( kind = 4 ) ierr
  integer ( kind = 4 ) ir
  integer ( kind = 4 ) j
  integer ( kind = 4 ) j_first
  integer ( kind = 4 ) j_new
  integer ( kind = 4 ) j_save
  integer ( kind = 4 ) next(2)
  integer ( kind = 4 ) nj
  real    ( kind = 8 ) x(n)
  real    ( kind = 8 ) xe(ne)
!
!  Check arguments.
!
  if ( n < 2 ) then
    ierr = -1
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'SPLINE_PCHIP_VAL - Fatal error!'
    write ( *, '(a)' ) '  Number of data points less than 2.'
    stop
  end if

  do i = 2, n
    if ( x(i) <= x(i-1) ) then
      ierr = -3
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'SPLINE_PCHIP_VAL - Fatal error!'
      write ( *, '(a)' ) '  X array not strictly increasing.'
      stop
    end if
  end do

  if ( ne < 1 ) then
    ierr = -4
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'SPLINE_PCHIP_VAL - Fatal error!'
    write ( *, '(a)' ) '  Number of evaluation points less than 1.'
    return
  end if

  ierr = 0
!
!  Loop over intervals.
!  The interval index is IL = IR-1.
!  The interval is X(IL) <= X < X(IR).
!
  j_first = 1
  ir = 2

  do
!
!  Skip out of the loop if have processed all evaluation points.
!
    if ( ne < j_first ) then
      exit
    end if
!
!  Locate all points in the interval.
!
    j_save = ne + 1

    do j = j_first, ne
      if ( x(ir) <= xe(j) ) then
        j_save = j
        if ( ir == n ) then
          j_save = ne + 1
        end if
        exit
      end if
    end do
!
!  Have located first point beyond interval.
!
    j = j_save

    nj = j - j_first
!
!  Skip evaluation if no points in interval.
!
    if ( nj /= 0 ) then
!
!  Evaluate cubic at XE(J_FIRST:J-1).
!
      call chfev ( x(ir-1), x(ir), f(ir-1), f(ir), d(ir-1), d(ir), &
        nj, xe(j_first:j-1), fe(j_first:j-1), next, ierc )

      if ( ierc < 0 ) then
        ierr = -5
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'SPLINE_PCHIP_VAL - Fatal error!'
        write ( *, '(a)' ) '  Error return from CHFEV.'
        stop
      end if
!
!  In the current set of XE points, there are NEXT(2) to the right of X(IR).
!
      if ( next(2) /= 0 ) then

        if ( ir < n ) then
          ierr = -5
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'SPLINE_PCHIP_VAL - Fatal error!'
          write ( *, '(a)' ) '  IR < N.'
          stop
        end if
!
!  These are actually extrapolation points.
!
        ierr = ierr + next(2)

      end if
!
!  In the current set of XE points, there are NEXT(1) to the left of X(IR-1).
!
      if ( next(1) /= 0 ) then
!
!  These are actually extrapolation points.
!
        if ( ir <= 2 ) then
          ierr = ierr + next(1)
        else

          j_new = -1

          do i = j_first, j-1
            if ( xe(i) < x(ir-1) ) then
              j_new = i
              exit
            end if
          end do

          if ( j_new == -1 ) then
            ierr = -5
            write ( *, '(a)' ) ' '
            write ( *, '(a)' ) 'SPLINE_PCHIP_VAL - Fatal error!'
            write ( *, '(a)' ) '  Could not bracket the data point.'
            stop
          end if
!
!  Reset J.  This will be the new J_FIRST.
!
          j = j_new
!
!  Now find out how far to back up in the X array.
!
          do i = 1, ir-1
            if ( xe(j) < x(i) ) then
              exit
            end if
          end do
!
!  At this point, either XE(J) < X(1) or X(i-1) <= XE(J) < X(I) .
!
!  Reset IR, recognizing that it will be incremented before cycling.
!
          ir = max ( 1, i-1 )

        end if

      end if

      j_first = j

    end if

    ir = ir + 1

    if ( n < ir ) then
      exit
    end if

  end do

  return
end
subroutine spline_quadratic_val ( ndata, tdata, ydata, tval, yval, ypval )

!*****************************************************************************80
!
!! SPLINE_QUADRATIC_VAL evaluates a piecewise quadratic spline at a point.
!
!  Discussion:
!
!    Because of the simple form of a piecewise quadratic spline,
!    the raw data points ( TDATA(I), YDATA(I)) can be used directly to
!    evaluate the spline at any point.  No processing of the data
!    is required.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    24 October 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NDATA, the number of data points defining 
!    the spline.  NDATA should be odd and at least 3.
!
!    Input, real ( kind = 8 ) TDATA(NDATA), YDATA(NDATA), the values of 
!    the independent and dependent variables at the data points.  The 
!    values of TDATA should be distinct and increasing.
!
!    Input, real ( kind = 8 ) TVAL, the point at which the spline is to 
!    be evaluated.
!
!    Output, real ( kind = 8 ) YVAL, YPVAL, the value of the spline and 
!    its first derivative dYdT at TVAL.  YPVAL is not reliable if TVAL 
!    is exactly equal to TDATA(I) for some I.
!
  implicit none

  integer ( kind = 4 ) ndata

  real    ( kind = 8 ) dif1
  real    ( kind = 8 ) dif2
  integer ( kind = 4 ), parameter :: i4_2 = 2
  integer ( kind = 4 ) left
  integer ( kind = 4 ) right
  real    ( kind = 8 ) t1
  real    ( kind = 8 ) t2
  real    ( kind = 8 ) t3
  real    ( kind = 8 ) tdata(ndata)
  real    ( kind = 8 ) tval
  real    ( kind = 8 ) y1
  real    ( kind = 8 ) y2
  real    ( kind = 8 ) y3
  real    ( kind = 8 ) ydata(ndata)
  real    ( kind = 8 ) ypval
  real    ( kind = 8 ) yval

  if ( ndata < 3 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'SPLINE_QUADRATIC_VAL - Fatal error!'
    write ( *, '(a)' ) '  NDATA < 3.'
    stop
  end if

  if ( mod ( ndata, i4_2 ) == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'SPLINE_QUADRATIC_VAL - Fatal error!'
    write ( *, '(a)' ) '  NDATA must be odd.'
    stop
  end if
!
!  Find the interval [ TDATA(LEFT), TDATA(RIGHT) ] that contains, or is
!  nearest to, TVAL.
!
  call r8vec_bracket ( ndata, tdata, tval, left, right )
!
!  Force LEFT to be odd.
!
  if ( mod ( left, i4_2 ) == 0 ) then
    left = left - 1
  end if
!
!  Copy out the three abscissas.
!
  t1 = tdata(left)
  t2 = tdata(left+1)
  t3 = tdata(left+2)

  if ( t2 <= t1 .or. t3 <= t2 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'SPLINE_QUADRATIC_VAL - Fatal error!'
    write ( *, '(a)' ) '  T2 <= T1 or T3 <= T2.'
    stop
  end if
!
!  Construct and evaluate a parabolic interpolant for the data
!  in each dimension.
!
  y1 = ydata(left)
  y2 = ydata(left+1)
  y3 = ydata(left+2)

  dif1 = ( y2 - y1 ) / ( t2 - t1 )

  dif2 = ( ( y3 - y1 ) / ( t3 - t1 ) &
       - ( y2 - y1 ) / ( t2 - t1 ) ) / ( t3 - t2 )

  yval = y1 + ( tval - t1 ) * ( dif1 + ( tval - t2 ) * dif2 )
  ypval = dif1 + dif2 * ( 2.0D+00 * tval - t1 - t2 )

  return
end
subroutine timestamp ( )

!*****************************************************************************80
!
!! TIMESTAMP prints the current YMDHMS date as a time stamp.
!
!  Example:
!
!    31 May 2001   9:45:54.872 AM
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    06 August 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    None
!
  implicit none

  character ( len = 8 ) ampm
  integer ( kind = 4 ) d
  integer ( kind = 4 ) h
  integer ( kind = 4 ) m
  integer ( kind = 4 ) mm
  character ( len = 9 ), parameter, dimension(12) :: month = (/ &
    'January  ', 'February ', 'March    ', 'April    ', &
    'May      ', 'June     ', 'July     ', 'August   ', &
    'September', 'October  ', 'November ', 'December ' /)
  integer ( kind = 4 ) n
  integer ( kind = 4 ) s
  integer ( kind = 4 ) values(8)
  integer ( kind = 4 ) y

  call date_and_time ( values = values )

  y = values(1)
  m = values(2)
  d = values(3)
  h = values(5)
  n = values(6)
  s = values(7)
  mm = values(8)

  if ( h < 12 ) then
    ampm = 'AM'
  else if ( h == 12 ) then
    if ( n == 0 .and. s == 0 ) then
      ampm = 'Noon'
    else
      ampm = 'PM'
    end if
  else
    h = h - 12
    if ( h < 12 ) then
      ampm = 'PM'
    else if ( h == 12 ) then
      if ( n == 0 .and. s == 0 ) then
        ampm = 'Midnight'
      else
        ampm = 'AM'
      end if
    end if
  end if

  write ( *, '(i2,1x,a,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3,1x,a)' ) &
    d, trim ( month(m) ), y, h, ':', n, ':', s, '.', mm, trim ( ampm )

  return
end
