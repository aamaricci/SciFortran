program main

!*****************************************************************************80
!
!! MINPACK_PRB runs the MINPACK tests.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    30 December 2004
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'MINPACK_PRB'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test the MINPACK library.'

  call test01 ( )
  call test02 ( )
  call test03 ( )
  call test04 ( )
  call test05 ( )
  call test06 ( )
  call test07 ( )
  call test08 ( )
  call test09 ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'MINPACK_PRB'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine test01 ( )

!*****************************************************************************80
!
!! TEST01 tests CHKDER.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    30 December 2004
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 5
  integer ( kind = 4 ), parameter :: m = n
  integer ( kind = 4 ), parameter :: ldfjac = n

  real ( kind = 8 ) err(m)
  real ( kind = 8 ) fjac(ldfjac,n)
  real ( kind = 8 ) fvec(m)
  real ( kind = 8 ) fvecp(m)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ido
  integer ( kind = 4 ) iflag
  integer ( kind = 4 ) j
  integer ( kind = 4 ) mode
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) xp(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST01'
  write ( *, '(a)' ) '  CHKDER compares a user supplied jacobian'
  write ( *, '(a)' ) '  and a finite difference approximation to it'
  write ( *, '(a)' ) '  and judges whether the jacobian is correct.'

  do ido = 1, 2

    if ( ido == 1 ) then

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  On the first test, use a correct jacobian.'

    else if ( ido == 2 ) then

       write ( *, '(a)' ) ' '
       write ( *, '(a)' ) '  Repeat the test, but use a "bad" jacobian'
       write ( *, '(a)' ) '  and see if the routine notices!'

     end if
!
!  Set the point at which the test is to be made:
!
    x(1:n) = 0.5D+00

    call r8vec_print ( n, x, '  Evaluation point X:' )

    mode = 1
    call chkder ( m, n, x, fvec, fjac, ldfjac, xp, fvecp, mode, err )

    iflag = 1

    call f01 ( n, x, fvec, fjac, ldfjac, iflag )
    call f01 ( n, xp, fvecp, fjac, ldfjac, iflag )

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Sampled function values F(X) and F(XP)'
    write ( *, '(a)' ) ' '
    do i = 1, m
      write ( *, '(i3,2g14.6)' ) i, fvec(i), fvecp(i)
    end do

    iflag = 2
    call f01 ( n, x, fvec, fjac, ldfjac, iflag )
!
!  Here's where we put a mistake into the jacobian, on purpose.
!
    if ( ido == 2 ) then
      fjac(1,1) = 1.01D+00 * fjac(1,1)
      fjac(2,3) = - fjac(2,3)
    end if

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Computed jacobian'
    write ( *, '(a)' ) ' '
    do i = 1, m
      write ( *, '(5g14.6)' ) fjac(i,1:n)
    end do

    mode = 2
    call chkder ( m, n, x, fvec, fjac, ldfjac, xp, fvecp, mode, err )

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  CHKDER gradient component error estimates:'
    write ( *, '(a)' ) '     > 0.5, the component is probably correct.'
    write ( *, '(a)' ) '     < 0.5, the component is probably incorrect.'
    write ( *, '(a)' ) ' '
    do i = 1, m
      write ( *, '(i6,g14.6)' ) i, err(i)
    end do

  end do

  return
end
subroutine f01 ( n, x, fvec, fjac, ldfjac, iflag )

!*****************************************************************************80
!
!! F01 is a function/jacobian routine.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 May 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of variables.
!
!    Input, real ( kind = 8 ) X(N), the variable values.
!
!    Output, real ( kind = 8 ) FVEC(N), the function values at X,
!    if IFLAG = 1.
!
!    Output, real ( kind = 8 ) FJAC(LDFJAC,N), the N by N jacobian at X,
!    if IFLAG = 2.
!
!    Input, integer ( kind = 4 ) LDFJAC, the leading dimension of FJAC,
!    which must be at least N.
!
!    Input, integer ( kind = 4 ) IFLAG:
!    1, please compute F(I) (X).
!    2, please compute FJAC(I,J) (X).
!
  implicit none

  integer ( kind = 4 ) ldfjac
  integer ( kind = 4 ) n

  real ( kind = 8 ) fjac(ldfjac,n)
  real ( kind = 8 ) fvec(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) iflag
  integer ( kind = 4 ) j
  real ( kind = 8 ) prod
  real ( kind = 8 ) x(n)
!
!  If IFLAG is 1, we are supposed to evaluate F(X).
!
  if ( iflag == 1 ) then

    do i = 1, n - 1
      fvec(i) = x(i) - real ( n + 1, kind = 8 ) + sum ( x(1:n) )
    end do

    fvec(n) = product ( x(1:n) ) - 1.0D+00
!
!  If IFLAG is 2, we are supposed to evaluate FJAC(I,J) = d F(I)/d X(J)
!
  else if ( iflag == 2 ) then

    fjac(1:n-1,1:n) = 1.0D+00

    do i = 1, n - 1
      fjac(i,i) = 2.0D+00
    end do

    prod = product ( x(1:n) )

    do j = 1, n
      fjac(n,j) = prod / x(j)
    end do

  end if

  return
end
subroutine test02 ( )

!*****************************************************************************80
!
!! TEST02 tests HYBRD1.
!
!  Discussion:
!
!    This is an example of what your main program would look
!    like if you wanted to use MINPACK to solve N nonlinear equations
!    in N unknowns.  In this version, we avoid computing the jacobian
!    matrix, and request that MINPACK approximate it for us.
!
!    The set of nonlinear equations is:
!
!      x1**2 - 10 * x1 + x2**2 + 8 = 0
!      x1 * x2**2 + x1 - 10 * x2 + 8 = 0
!
!    with solution x1 = x2 = 1
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    30 December 2004
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 2

  external f02
  real ( kind = 8 ) fvec(n)
  integer ( kind = 4 ) iflag
  integer ( kind = 4 ) info
  real ( kind = 8 ) tol
  real ( kind = 8 ) x(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST02'
  write ( *, '(a)' ) '  HYBRD1 solves a nonlinear system of equations.'

  x(1:2) = (/ 3.0D+00, 0.0D+00 /)
  call r8vec_print ( n, x, '  Initial X:' )
  iflag = 1
  call f02 ( n, x, fvec, iflag )
  call r8vec_print ( n, fvec, '  F(X):' )

  tol = 0.00001D+00

  call hybrd1 ( f02, n, x, fvec, tol, info )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i6)' ) '  Returned value of INFO = ', info
  call r8vec_print ( n, x, '  X:' )
  call r8vec_print ( n, fvec, '  F(X):' )

  return
end
subroutine f02 ( n, x, fvec, iflag )

!*****************************************************************************80
!
!! F02 is a function routine.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    30 December 2004
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) fvec(n)
  integer ( kind = 4 ) iflag
  real ( kind = 8 ) x(n)

  fvec(1) = x(1) * x(1) - 10.0D+00 * x(1) + x(2) * x(2) + 8.0D+00
  fvec(2) = x(1) * x(2) * x(2) + x(1) - 10.0D+00 * x(2) + 8.0D+00

  return
end
subroutine test03 ( )

!*****************************************************************************80
!
!! TEST03 tests HYBRJ1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    30 December 2004
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 2
  integer ( kind = 4 ), parameter :: ldfjac = n

  external f03
  real ( kind = 8 ) fjac(ldfjac,n)
  real ( kind = 8 ) fvec(n)
  integer ( kind = 4 ) iflag
  integer ( kind = 4 ) info
  real ( kind = 8 ) tol
  real ( kind = 8 ) x(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST03'
  write ( *, '(a)' ) '  HYBRJ1 solves a nonlinear system of equations.'

  x(1:2) = (/ 3.0D+00, 0.0D+00 /)
  call r8vec_print ( n, x, '  Initial X:' )
  iflag = 1
  call f02 ( n, x, fvec, iflag )
  call r8vec_print ( n, fvec, '  F(X):' )

  tol = 0.00001D+00

  call hybrj1 ( f03, n, x, fvec, fjac, ldfjac, tol, info )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i6)' ) '  Returned value of INFO = ', info
  call r8vec_print ( n, x, '  X:' )
  call r8vec_print ( n, fvec, '  F(X):' )

  return
end
subroutine f03 ( n, x, fvec, fjac, ldfjac, iflag )

!*****************************************************************************80
!
!! F03 is a function/jacobian routine.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    30 December 2004
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) ldfjac
  integer ( kind = 4 ) n

  real ( kind = 8 ) fjac(ldfjac,n)
  real ( kind = 8 ) fvec(n)
  integer ( kind = 4 ) iflag
  real ( kind = 8 ) x(n)

  if ( iflag == 1 ) then

    fvec(1) = x(1) * x(1) - 10.0D+00 * x(1) + x(2) * x(2) + 8.0D+00
    fvec(2) = x(1) * x(2) * x(2) + x(1) - 10.0D+00 * x(2) + 8.0D+00

  else if ( iflag == 2 ) then

    fjac(1,1) = 2.0D+00 * x(1) - 10.0D+00
    fjac(1,2) = 2.0D+00 * x(2)
    fjac(2,1) = x(2) * x(2) + 1.0D+00
    fjac(2,2) = 2.0D+00 * x(1) * x(2) - 10.0D+00

  end if

  return
end
subroutine test04 ( )

!*****************************************************************************80
!
!! TEST04 tests LMDER1.
!
!  Discussion:
!
!    LMDER1 solves M nonlinear equations in N unknowns, with M >= N.
!
!    LMDER1 seeks a solution X minimizing the euclidean norm of the residual.
!
!    In this example, the set of equations is actually linear, but
!    normally they are nonlinear.
!
!    In this problem, we have a set of pairs of data points, and we
!    seek a functional relationship between them.  We assume the
!    relationship is of the form
!
!      y = a * x + b
!
!    and we want to know the values of a and b.  Therefore, we would like
!    to find numbers a and b which satisfy a set of equations.
!
!    The data points are (2,2), (4,11), (6,28) and (8,40).
!
!    Therefore, the equations we want to satisfy are:
!
!      a * 2 + b -  2 = 0
!      a * 4 + b - 11 = 0
!      a * 6 + b - 28 = 0
!      a * 8 + b - 40 = 0
!
!    The least squares solution of this system is a=6.55, b=-12.5,
!    In other words, the line y=6.55*x-12.5 is the line which "best"
!    models the data in the least squares sense.
!
!    Problems with more variables, or higher degree polynomials, would
!    be solved similarly.  For example, suppose we have (x,y,z) data,
!    and we wish to find a relationship of the form f(x,y,z).  We assume
!    that x and y occur linearly, and z quadratically.  Then the equation
!    we seek has the form:
!
!      a*x+b*y+c*z + d*z*z + e = 0
!
!    and, supposing that our first two points were (1,2,3), (1,3,8), our set of
!    equations would begin:
!
!      a*1+b*2+c*3 + d*9  + e = 0
!      a*1+b*3+c*8 + d*64 + e = 0
!
!    and so on.
!
!    M is the number of equations, which in this case is the number of
!    (x,y) data values.
!
!    N is the number of variables, which in this case is the number of
!    'free' coefficients in the relationship we are trying to determine.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    31 October 2005
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: m = 4
  integer ( kind = 4 ), parameter :: n = 2
  integer ( kind = 4 ), parameter :: ldfjac = m

  external f04
  real ( kind = 8 ) fjac(ldfjac,n)
  real ( kind = 8 ) fvec(m)
  integer ( kind = 4 ) iflag
  integer ( kind = 4 ) info
  real ( kind = 8 ) tol
  real ( kind = 8 ) x(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST04'
  write ( *, '(a)' ) '  LMDER1 minimizes M functions in N variables.'

  x(1:2) = (/ 0.0D+00, 5.0D+00 /)
  call r8vec_print ( n, x, '  Initial X:' )
  iflag = 1
  call f04 ( m, n, x, fvec, fjac, ldfjac, iflag )
  call r8vec_print ( m, fvec, '  F(X):' )

  tol = 0.00001D+00

  call lmder1 ( f04, m, n, x, fvec, fjac, ldfjac, tol, info )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i6)' ) '  Returned value of INFO = ', info
  call r8vec_print ( n, x, '  X:' )
  call r8vec_print ( m, fvec, '  F(X):' )

  return
end
subroutine f04 ( m, n, x, fvec, fjac, ldfjac, iflag )

!*****************************************************************************80
!
!! F04 is a function/jacobian routine.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    30 December 2004
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) ldfjac
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) fjac(ldfjac,n)
  real ( kind = 8 ) fvec(m)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) iflag
  real ( kind = 8 ) x(n)
  real ( kind = 8 ), parameter, dimension ( 4 ) :: xdat = (/ &
    2.0D+00,  4.0D+00,  6.0D+00,  8.0D+00 /)
  real ( kind = 8 ), parameter, dimension ( 4 ) :: ydat = (/ &
    2.0D+00, 11.0D+00, 28.0D+00, 40.0D+00 /)

  if ( iflag == 1 ) then

    fvec(1:m) = x(1) * xdat(1:m) + x(2) - ydat(1:m)

  else if ( iflag == 2 ) then

    fjac(1:m,1) = xdat(1:m)
    fjac(1:m,2) = 1.0D+00

  else

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'F04 - Fatal error!'
    write ( *, '(a,i6)' ) '  Called with unexpected value of IFLAG = ', iflag
    stop

  end if

  return
end
subroutine test05 ( )

!*****************************************************************************80
!
!! TEST05 tests LMDER1.
!
!  Discussion:
!
!    LMDER1 solves M nonlinear equations in n unknowns, where M is greater
!    than N.  The functional fit is nonlinear this time, of the form
!
!      y=a+b*x**c,
!
!    with x and y data, and a, b and c unknown.
!
!    This problem is set up so that the data is exactly fit by by
!    a=1, b=3, c=2.  Normally, the data would only be approximately
!    fit by the best possible solution.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    30 December 2004
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: m = 10
  integer ( kind = 4 ), parameter :: n = 3
  integer ( kind = 4 ), parameter :: ldfjac = m

  external f05
  real ( kind = 8 ) fjac(ldfjac,n)
  real ( kind = 8 ) fvec(m)
  integer ( kind = 4 ) iflag
  integer ( kind = 4 ) info
  real ( kind = 8 ) tol
  real ( kind = 8 ) x(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST05'
  write ( *, '(a)' ) '  LMDER1 minimizes M functions in N variables.'

  x(1:3) = (/ 0.0D+00, 5.0D+00, 1.3D+00 /)
  call r8vec_print ( n, x, '  Initial X:' )
  iflag = 1
  call f05 ( m, n, x, fvec, fjac, ldfjac, iflag )
  call r8vec_print ( m, fvec, '  F(X):' )

  tol = 0.00001D+00

  call lmder1 ( f05, m, n, x, fvec, fjac, ldfjac, tol, info )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i6)' ) '  Returned value of INFO = ', info
  call r8vec_print ( n, x, '  X:' )
  call r8vec_print ( m, fvec, '  F(X):' )

  return
end
subroutine f05 ( m, n, x, fvec, fjac, ldfjac, iflag )

!*****************************************************************************80
!
!! F05 is a function/jacobian routine.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    30 December 2004
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) ldfjac
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) fjac(ldfjac,n)
  real ( kind = 8 ) fvec(m)
  integer ( kind = 4 ) iflag
  real ( kind = 8 ) x(n)
  real ( kind = 8 ), dimension ( 10 ) :: xdat = (/ &
    1.0D+00, 2.0D+00, 3.0D+00, 4.0D+00, 5.0D+00, &
    6.0D+00, 7.0D+00, 8.0D+00, 9.0D+00, 10.0D+00 /)
  real ( kind = 8 ), dimension ( 10 ) :: ydat = (/ &
    4.0D+00, 13.0D+00, 28.0D+00, 49.0D+00, 76.0D+00, &
    109.0D+00, 148.0D+00, 193.0D+00, 244.0D+00, 301.0D+00 /)

  if ( iflag == 1 ) then

    fvec(1:m) = x(1) + x(2) * xdat(1:m)**x(3) - ydat(1:m)

  else if ( iflag == 2 ) then

    fjac(1:m,1) = 1.0D+00
    fjac(1:m,2) = xdat(1:m)**x(3)
    fjac(1:m,3) = x(2) * log ( xdat(1:m) ) * xdat(1:m)**x(3)

  end if

  return
end
subroutine test06 ( )

!*****************************************************************************80
!
!! TEST06 tests LMDIF1.
!
!  Discussion:
!
!    LMDIF1 solves M nonlinear equations in N unknowns, where M is greater
!    than N.  Generally, you cannot get a solution vector x which will satisfy
!    all the equations.  That is, the vector equation f(x)=0 cannot
!    be solved exactly.  Instead, minpack seeks a solution x so that
!    the euclidean norm transpose(f(x))*f(x) is minimized.  The size
!    of the euclidean norm is a measure of how good the solution is.
!
!    In this example, the set of equations is actually linear, but
!    normally they are nonlinear.
!
!    In this problem, we have a set of pairs of data points, and we
!    seek a functional relationship between them.  We assume the
!    relationship is of the form
!
!      y=a*x+b
!
!    and we want to know the values of a and b.  Therefore, we would like
!    to find numbers a and b which satisfy a set of equations.
!
!    The data points are (2,2), (4,11), (6,28) and (8,40).
!
!    Therefore, the equations we want to satisfy are:
!
!      a * 2 + b -  2 = 0
!      a * 4 + b - 11 = 0
!      a * 6 + b - 28 = 0
!      a * 8 + b - 40 = 0
!
!    The least squares solution of this system is a=6.55, b=-12.5,
!    In other words, the line y=6.55*x-12.5 is the line which "best"
!    models the data in the least squares sense.
!
!    Problems with more variables, or higher degree polynomials, would
!    be solved similarly.  For example, suppose we have (x,y,z) data,
!    and we wish to find a relationship of the form f(x,y,z).  We assume
!    that x and y occur linearly, and z quadratically.  Then the equation
!    we seek has the form:
!
!      a*x+b*y+c*z + d*z*z + e = 0
!
!    and, supposing that our first two points were (1,2,3), (1,3,8), our set of
!    equations would begin:
!
!      a*1+b*2+c*3 + d*9  + e = 0
!      a*1+b*3+c*8 + d*64 + e = 0
!
!    and so on.
!
!    M is the number of equations, which in this case is the number of
!    (x,y) data values.
!
!    N is the number of variables, which in this case is the number of
!    'free' coefficients in the relationship we are trying to determine.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    30 December 2004
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: m = 4
  integer ( kind = 4 ), parameter :: n = 2

  external f06
  real ( kind = 8 ) fvec(m)
  integer ( kind = 4 ) iflag
  integer ( kind = 4 ) info
  real ( kind = 8 ) tol
  real ( kind = 8 ) x(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST06'
  write ( *, '(a)' ) '  LMDIF1 minimizes M functions in N variables.'

  x(1:2) = (/ 0.0D+00, 5.0D+00 /)
  call r8vec_print ( n, x, '  Initial X:' )
  iflag = 1
  call f06 ( m, n, x, fvec, iflag )
  call r8vec_print ( m, fvec, '  F(X):' )

  tol = 0.00001D+00

  call lmdif1 ( f06, m, n, x, fvec, tol, info )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i6)' ) '  Returned value of INFO = ', info
  call r8vec_print ( n, x, '  X:' )
  call r8vec_print ( m, fvec, '  F(X):' )

  return
end
subroutine f06 ( m, n, x, fvec, iflag )

!*****************************************************************************80
!
!! F06 is a function routine.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    30 December 2004
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) fvec(m)
  integer ( kind = 4 ) iflag
  real ( kind = 8 ) x(n)
  real ( kind = 8 ), dimension ( 4 ) :: xdat = (/ &
    2.0D+00,  4.0D+00,  6.0D+00,  8.0D+00 /)
  real ( kind = 8 ), dimension ( 4 ) :: ydat = (/ &
    2.0D+00, 11.0D+00, 28.0D+00, 40.0D+00 /)

  fvec(1:m) = x(1) * xdat(1:m) + x(2) - ydat(1:m)

  return
end
subroutine test07 ( )

!*****************************************************************************80
!
!! TEST07 tests LMDIF1.
!
!  Discussion:
!
!    LMDIF1 solves M nonlinear equations in N unknowns, where M is greater
!    than N.  It is similar to test02, except that the functional fit is
!    nonlinear this time, of the form
!
!      y = a + b * x**c,
!
!    with x and y data, and a, b and c unknown.
!
!    This problem is set up so that the data is exactly fit by by
!    a=1, b=3, c=2.  Normally, the data would only be approximately
!    fit by the best possible solution.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    30 December 2004
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: m = 10
  integer ( kind = 4 ), parameter :: n = 3

  external f07
  real ( kind = 8 ) fvec(m)
  integer ( kind = 4 ) iflag
  integer ( kind = 4 ) info
  real ( kind = 8 ) tol
  real ( kind = 8 ) x(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST07'
  write ( *, '(a)' ) '  LMDIF1 minimizes M functions in N variables.'

  x(1:3) = (/ 0.0D+00, 5.0D+00, 1.3D+00 /)
  call r8vec_print ( n, x, '  X:' )
  iflag = 1
  call f07 ( m, n, x, fvec, iflag )
  call r8vec_print ( m, fvec, '  F(X):' )

  tol = 0.00001D+00

  call lmdif1 ( f07, m, n, x, fvec, tol, info )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i6)' ) '  Returned value of INFO = ', info
  call r8vec_print ( n, x, '  X:' )
  call r8vec_print ( m, fvec, '  F(X):' )

  return
end
subroutine f07 ( m, n, x, fvec, iflag )

!*****************************************************************************80
!
!! F07 is a function routine.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    30 December 2004
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) fvec(m)
  integer ( kind = 4 ) iflag
  real ( kind = 8 ) x(n)
  real ( kind = 8 ), dimension ( 10 ) :: xdat = (/ &
    1.0D+00, 2.0D+00, 3.0D+00, 4.0D+00, 5.0D+00, &
    6.0D+00, 7.0D+00, 8.0D+00, 9.0D+00, 10.0D+00 /)
  real ( kind = 8 ), dimension ( 10 ) :: ydat = (/ &
    4.0D+00, 13.0D+00, 28.0D+00, 49.0D+00, 76.0D+00, &
    109.0D+00, 148.0D+00, 193.0D+00, 244.0D+00, 301.0D+00 /)

  fvec(1:m) = x(1) + x(2) * xdat(1:m)**x(3) - ydat(1:m)

  return
end
subroutine test08 ( )

!*****************************************************************************80
!
!! TEST08 tests LMSTR1.
!
!  Discussion:
!
!    LMSTR1 solves M nonlinear equations in N unknowns, where M is greater
!    than N.  Generally, you cannot get a solution vector x which will satisfy
!    all the equations.  That is, the vector equation f(x)=0 cannot
!    be solved exactly.  Instead, minpack seeks a solution x so that
!    the euclidean norm transpose(f(x))*f(x) is minimized.  The size
!    of the euclidean norm is a measure of how good the solution is.
!
!    In this example, the set of equations is actually linear, but
!    normally they are nonlinear.
!
!    In this problem, we have a set of pairs of data points, and we
!    seek a functional relationship between them.  We assume the
!    relationship is of the form
!
!      y=a*x+b
!
!    and we want to know the values of a and b.  Therefore, we would like
!    to find numbers a and b which satisfy a set of equations.
!
!    The data points are (2,2), (4,11), (6,28) and (8,40).
!
!    Therefore, the equations we want to satisfy are:
!
!      a * 2 + b -  2 = 0
!      a * 4 + b - 11 = 0
!      a * 6 + b - 28 = 0
!      a * 8 + b - 40 = 0
!
!    The least squares solution of this system is a=6.55, b=-12.5,
!    In other words, the line y=6.55*x-12.5 is the line which "best"
!    models the data in the least squares sense.
!
!    Problems with more variables, or higher degree polynomials, would
!    be solved similarly.  For example, suppose we have (x,y,z) data,
!    and we wish to find a relationship of the form f(x,y,z).  We assume
!    that x and y occur linearly, and z quadratically.  Then the equation
!    we seek has the form:
!
!      a*x + b*y + c*z + d*z*z + e = 0
!
!    and, supposing that our first two points were (1,2,3), (1,3,8), our set of
!    equations would begin:
!
!      a*1 + b*2 + c*3 + d*9  + e = 0
!      a*1 + b*3 + c*8 + d*64 + e = 0
!
!    and so on.
!
!    M is the number of equations, which in this case is the number of
!    (x,y) data values.
!
!    N is the number of variables, which in this case is the number of
!    'free' coefficients in the relationship we are trying to determine.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    30 December 2004
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: m = 4
  integer ( kind = 4 ), parameter :: n = 2
  integer ( kind = 4 ), parameter :: ldfjac = m

  external f08
  real ( kind = 8 ) fjac(ldfjac,n)
  integer ( kind = 4 ) fjrow
  real ( kind = 8 ) fvec(m)
  integer ( kind = 4 ) iflag
  integer ( kind = 4 ) info
  real ( kind = 8 ) tol
  real ( kind = 8 ) x(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST08'
  write ( *, '(a)' ) '  LMSTR1 minimizes M functions in N variables.'

  x(1:2) = (/ 0.0D+00, 5.0D+00 /)
  call r8vec_print ( n, x, '  Initial X:' )
  iflag = 1
  call f08 ( m, n, x, fvec, fjrow, iflag )
  call r8vec_print ( m, fvec, '  F(X):' )

  tol = 0.00001D+00

  call lmstr1 ( f08, m, n, x, fvec, fjac, ldfjac, tol, info )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i6)' ) '  Returned value of INFO = ', info
  call r8vec_print ( n, x, '  X:' )
  call r8vec_print ( m, fvec, '  F(X):' )

  return
end
subroutine f08 ( m, n, x, fvec, fjrow, iflag )

!*****************************************************************************80
!
!! F08 is a function/jacobian routine.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    30 December 2004
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) fjrow(n)
  real ( kind = 8 ) fvec(m)
  integer ( kind = 4 ) iflag
  real ( kind = 8 ) x(n)
  real ( kind = 8 ), dimension ( 4 ) :: xdat = (/ &
    2.0D+00,  4.0D+00,  6.0D+00,  8.0D+00 /)
  real ( kind = 8 ), dimension ( 4 ) :: ydat = (/ &
    2.0D+00, 11.0D+00, 28.0D+00, 40.0D+00 /)

  if ( iflag == 1 ) then

    fvec(1:m) = x(1) * xdat(1:m) + x(2) - ydat(1:m)

  else

    fjrow(1) = xdat(iflag-1)
    fjrow(2) = 1.0D+00

  end if

  return
end
subroutine test09 ( )

!*****************************************************************************80
!
!! TEST09 tests LMSTR1.
!
!  Discussion:
!
!    LMSTR1 solves M nonlinear equations in N unknowns, where M is greater
!    than N.  This test is similar to test02, except that the functional fit
!    is nonlinear this time, of the form
!
!      y = a + b * x**c,
!
!    with x and y data, and a, b and c unknown.
!
!    This problem is set up so that the data is exactly fit by by
!    a=1, b=3, c=2.  Normally, the data would only be approximately
!    fit by the best possible solution.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    30 December 2004
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: m = 10
  integer ( kind = 4 ), parameter :: n = 3
  integer ( kind = 4 ), parameter :: ldfjac = m

  external f09
  real ( kind = 8 ) fjac(ldfjac,n)
  real ( kind = 8 ) fjrow(n)
  real ( kind = 8 ) fvec(m)
  integer ( kind = 4 ) iflag
  integer ( kind = 4 ) info
  real ( kind = 8 ) tol
  real ( kind = 8 ) x(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST09'
  write ( *, '(a)' ) '  LMSTR1 minimizes M functions in N variables.'

  x(1:3) = (/ 0.0D+00, 5.0D+00, 1.3D+00 /)
  call r8vec_print ( n, x, '  Initial X:' )
  iflag = 1
  call f09 ( m, n, x, fvec, fjrow, iflag )
  call r8vec_print ( m, fvec, '  F(X):' )

  tol = 0.00001D+00

  call lmstr1 ( f09, m, n, x, fvec, fjac, ldfjac, tol, info )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i6)' ) '  Returned value of INFO = ', info
  call r8vec_print ( n, x, '  X:' )
  call r8vec_print ( m, fvec, '  F(X):' )

  return
end
subroutine f09 ( m, n, x, fvec, fjrow, iflag )

!*****************************************************************************80
!
!! F09 is a function/jacobian routine.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    30 December 2004
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) fjrow(n)
  real ( kind = 8 ) fvec(m)
  integer ( kind = 4 ) iflag
  real ( kind = 8 ) x(n)
  real ( kind = 8 ), dimension ( 10 ) :: xdat = (/ &
    1.0D+00, 2.0D+00, 3.0D+00, 4.0D+00, 5.0D+00, &
    6.0D+00, 7.0D+00, 8.0D+00, 9.0D+00, 10.0D+00 /)
  real ( kind = 8 ), dimension ( 10 ) :: ydat = (/ &
    4.0D+00, 13.0D+00, 28.0D+00, 49.0D+00, 76.0D+00, &
    109.0D+00, 148.0D+00, 193.0D+00, 244.0D+00, 301.0D+00 /)

  if ( iflag == 1 ) then

    fvec(1:m) = x(1) + x(2) * xdat(1:m)**x(3) - ydat(1:m)

  else

    fjrow(1) = 1.0D+00
    fjrow(2) = xdat(iflag-1)**x(3)
    fjrow(3) = x(2) * log ( xdat(iflag-1) ) * xdat(iflag-1)**x(3)

  end if

  return
end
