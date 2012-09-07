program main

!*****************************************************************************80
!
!! MAIN is the main program for PPPACK_PRB.
!
!  Discussion:
!
!    PPPACK_PRB calls the PPPACK tests.
!
!  Modified:
!
!    14 February 2007
!
  implicit none

  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'PPPACK_PRB'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test the routines in the PPPACK library.'
 
  call test01 ( )
  call test02 ( )
  call test03 ( )
  call test04 ( )
  call test05 ( )
  call test06 ( )
  call test07 ( )
  call test08 ( )
  call test09 ( )

  call test10 ( )
  call test11 ( )
  call test12 ( )
  call test13 ( )
  call test14 ( )
  call test15 ( )
  call test16 ( )
  call test17 ( )
  call test18 ( )
  call test19 ( )

  call test20 ( )
  call test21 ( )
  call test22 ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'PPPACK_PRB'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine test01 ( )

!*****************************************************************************80
!
!! TEST01 demonstrates polynomial interpolation on the Runge function.
!
!  Discussion:
!
!    This example comes from chapter II of the reference.
!
!  Modified:
!
!    13 February 2007
!
!  Author:
!
!    Carl deBoor
!
!  Reference:
!
!    Carl DeBoor,
!    A Practical Guide to Splines,
!    Springer, 2001,
!    ISBN: 0387953663,
!    LC: QA1.A647.v27.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 20

  real ( kind = 8 ) algerp
  real ( kind = 8 ) aloger
  real ( kind = 8 ) d(n_max)
  real ( kind = 8 ) decay
  real ( kind = 8 ) dx
  real ( kind = 8 ) errmax
  integer ( kind = 4 ) i
  integer ( kind = 4 ) istep
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) n
  real ( kind = 8 ) pnatx
  real ( kind = 8 ), external :: runge
  real ( kind = 8 ) tau(n_max)
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST01'
  write ( *, '(a)' ) '  Demonstrate the failure of polynomial interpolation'
  write ( *, '(a)' ) '  when applied to the Runge function, using equally'
  write ( *, '(a)' ) '  spaced nodes.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Use polynomial interpolation of order N.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '   N   Max error   Decay exponent'
  write ( *, '(a)' ) ' '

  algerp = 0.0D+00
  decay = 0.0D+00
  istep = 20

  do n = 2, n_max, 2
!
!  Choose interpolation points TAU equally spaced in (-1,1).
!
    do i = 1, n
      tau(i) =  ( real ( n - i,     kind = 8 ) * ( -1.0D+00 )   &
                + real (     i - 1, kind = 8 ) * (  1.0D+00 ) ) &
                / real ( n     - 1, kind = 8 )
      d(i) = runge ( tau(i) )
    end do
!
!  Calculate the interpolating polynomial, using divided differences.
!
    do k = 1, n-1
      do i = 1, n-k
        d(i) = ( d(i+1) - d(i) ) / ( tau(i+k) - tau(i) )
      end do
    end do
!
!  Estimate the maximum interpolation error on (-1,1).
!
    errmax = 0.0D+00
!
!  Consider the subinterval ( TAU(I-1), TAU(I) ).
!
    do i = 2, n
!
!  Sample ISTEP points in the subinterval.
!
      do j = 1, istep

        x = ( real ( istep - j, kind = 8 ) * tau(i-1)   &
            + real (         j, kind = 8 ) * tau(i)   ) &
            / real ( istep,     kind = 8 )

        pnatx = d(1)
        do k = 2, n
          pnatx = d(k) + ( x - tau(k) ) * pnatx
        end do

        errmax = max ( errmax, abs ( runge ( x ) - pnatx) )

      end do
    end do

    aloger = log ( errmax )

    if ( 2 < n ) then
      decay = ( aloger - algerp ) &
            / log ( real ( n, kind = 8 ) /real ( n - 2, kind = 8 ) )
    end if

    algerp = aloger
    write ( *, '(i4,e12.4,f11.2)' ) n, errmax, decay

  end do

  return
end
subroutine test02 ( )

!*****************************************************************************80
!
!! TEST02 demonstrates cubic Hermite interpolation on the Runge function.
!
!  Discussion:
!
!    This example comes from chapter IV of the reference.
!
!  Modified:
!
!    14 February 2007
!
!  Author:
!
!    Carl deBoor
!
!  Reference:
!
!    Carl DeBoor,
!    A Practical Guide to Splines,
!    Springer, 2001,
!    ISBN: 0387953663,
!    LC: QA1.A647.v27.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 20

  real ( kind = 8 ) algerp
  real ( kind = 8 ) aloger
  real ( kind = 8 ) c(4,n_max)
  real ( kind = 8 ) decay
  real ( kind = 8 ) divdf1
  real ( kind = 8 ) divdf3
  real ( kind = 8 ) dtau
  real ( kind = 8 ) dx
  real ( kind = 8 ) errmax
  real ( kind = 8 ) h
  integer ( kind = 4 ) i
  integer ( kind = 4 ) istep
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) n
  real ( kind = 8 ) pnatx
  real ( kind = 8 ), external :: runge
  real ( kind = 8 ), external :: rungep
  real ( kind = 8 ) step
  real ( kind = 8 ) tau(n_max)
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST02'
  write ( *, '(a)' ) '  Cubic Hermite interpolation appplied to Runge function,'
  write ( *, '(a)' ) '  equally spaced nodes.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Use polynomial interpolation of order N.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '   N   Max error   Decay exponent'
  write ( *, '(a)' ) ' '

  algerp = 0.0D+00
  decay = 0.0D+00
  istep = 20
  step = 20.0D+00

  do n = 2, n_max, 2
!
!  Choose interpolation points TAU equally spaced in (-1,1).
!
    do i = 1, n
      tau(i) =  ( real ( n - i,     kind = 8 ) * ( -1.0D+00 )   &
                + real (     i - 1, kind = 8 ) * (  1.0D+00 ) ) &
                / real ( n     - 1, kind = 8 )
    end do

    do i = 1, n
      c(1,i) = runge ( tau(i) )
      c(2,i) = rungep ( tau(i) )
    end do
!
!  Calculate the coefficients of the polynomial pieces.
!
    do i = 1, n-1
      dtau = tau(i+1) - tau(i)
      divdf1 = ( c(1,i+1) - c(1,i) ) / dtau
      divdf3 = c(2,i) + c(2,i+1) - 2.0D+00 * divdf1
      c(3,i) = ( divdf1 - c(2,i) - divdf3 ) / dtau
      c(4,i) = ( divdf3 / dtau ) / dtau
    end do
    c(3,n) = 0.0D+00
    c(4,n) = 0.0D+00
!
!  Estimate the maximum interpolation error on (-1,1).
!
    errmax = 0.0D+00
!
!  Consider the subinterval ( TAU(I-1), TAU(I) ).
!
    do i = 2, n

      dx = ( tau(i) - tau(i-1) ) / step
!
!  Sample ISTEP points in the subinterval.
!
      do j = 1, istep

        h = real ( j, kind = 8 ) * dx

        x = ( real ( istep - j, kind = 8 ) * tau(i-1) &
            + real (         j, kind = 8 ) * tau(i) ) &
            / real ( istep,     kind = 8 )

        pnatx = c(1,i-1) + h * ( &
                c(2,i-1) + h * ( &
                c(3,i-1) + h *   &
                c(4,i-1) ) )

!       errmax = max ( errmax, abs ( runge ( x ) - pnatx ) )
        errmax = max ( errmax, abs ( runge ( tau(i-1)+h ) - pnatx ) )

      end do
    end do

    aloger = log ( errmax )

    if ( 2 < n ) then
      decay = ( aloger - algerp ) &
            / log ( real ( n, kind = 8 ) /real ( n - 2, kind = 8 ) )
    end if

    algerp = aloger
    write ( *, '(i4,e12.4,f11.2)' ) n, errmax, decay

  end do

  return
end
subroutine test03 ( )

!*****************************************************************************80
!
!! TEST03 compares a B-spline representation to its values.
!
!  Discussion:
!
!    This example compares the B-spline representation of a cubic 
!    function with its values at knot averages.
!
!    This example comes from chapter IX of the reference.
!
!  Modified:
!
!    13 February 2007
!
!  Author:
!
!    Carl deBoor
!
!  Reference:
!
!    Carl DeBoor,
!    A Practical Guide to Splines,
!    Springer, 2001,
!    ISBN: 0387953663,
!    LC: QA1.A647.v27.
!
  implicit none

  real ( kind = 8 ) bcoef(23)
  real ( kind = 8 ), dimension ( 4 ) :: d0 = (/ &
    -162.0D+00, 99.0D+00, -18.0D+00, 1.0D+00 /)
  real ( kind = 8 ) dtip1
  real ( kind = 8 ) dtip2
  real ( kind = 8 ) f(23)
  real ( kind = 8 ) d(4)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) id
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jj
  integer ( kind = 4 ), parameter :: n = 13
  real ( kind = 8 ) t(27)
  real ( kind = 8 ) tave(23)
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST03'
  write ( *, '(a)' ) '  Compare B-spline coefficients to values'
  write ( *, '(a)' ) '  at knot averages.'
!
!  Set up the knot sequence in the array T.
!
  t(1:4) = 0.0D+00

  do i = 5, n
    t(i) = real ( i - 4, kind = 8 )
  end do

  t(n+1:n+4) = 10.0D+00

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  I   Tave(I)      F at Tave(I)      BCOEF(I)'
  write ( *, '(a)' ) ' '

  do i = 1, n
!
!  Use nested multiplication to get the Taylor coefficients D at
!  T(I+2) from those at 0.
!
    d(1:4) = d0(1:4)

    do j = 1, 3
      id = 4
      do jj = j, 3
        id = id - 1
        d(id) = d(id) + d(id+1) * t(i+2)
      end do
    end do
!
!  Compute the B spline coefficients by formula (9).
!
    dtip1 = t(i+2) - t(i+1)
    dtip2 = t(i+3) - t(i+2)

    bcoef(i) = d(1) &
           + ( d(2) * ( dtip2 - dtip1 ) &
             - d(3) * dtip1 * dtip2 ) / 3.0D+00
!
!  Evaluate F at the corresponding knot average.
!
    tave(i) = sum ( t(i+1:i+3) ) / 3.0D+00
    x = tave(i)
    f(i) = d0(1) + x * ( &
           d0(2) + x * ( &
           d0(3) + x *   &
           d0(4) ) )

    write ( *, '(i3,f10.5,f16.5,f16.5)' ) i, tave(i), f(i), bcoef(i)

  end do

  return
end
subroutine test04 ( )

!*****************************************************************************80
!
!! TEST04 "plots" some B splines.
!
!  Discussion:
!
!    This is example 1 from chapter X of the reference.
!
!  Modified:
!
!    13 February 2007
!
!  Author:
!
!    Carl deBoor
!
!  Reference:
!
!    Carl DeBoor,
!    A Practical Guide to Splines,
!    Springer, 2001,
!    ISBN: 0387953663,
!    LC: QA1.A647.v27.
!
  implicit none

  integer ( kind = 4 ), parameter :: k = 3
  integer ( kind = 4 ), parameter :: n = 7

  integer ( kind = 4 ) i
  integer ( kind = 4 ), parameter :: indx = 1
  integer ( kind = 4 ) left
  integer ( kind = 4 ) mflag
  integer ( kind = 4 ) :: npoint = 31
  real ( kind = 8 ), dimension (n+k) :: t = (/ &
    0.0D+00, 0.0D+00, 0.0D+00, 1.0D+00, 1.0D+00, &
    3.0D+00, 4.0D+00, 6.0D+00, 6.0D+00, 6.0D+00 /)
  real ( kind = 8 ) values(n)
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST04'
  write ( *, '(a)' ) '  INTERV locates a point X within a knot sequence;'
  write ( *, '(a)' ) '  BSPLVB evaluates a B-spline.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Use these routines to plot some B splines.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  B-spline order K = ', k
  write ( *, '(a,i8)' ) '  N = ', n

  call r8vec_print ( n+k, t, '  Knot sequence T:' )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) &
    '    X      LEFT     B1(X)       B2(X)       B3(X)       B4(X)       B5(X)'
  write ( *, '(a)' ) ' '
!
!  Initialize the B spline values.
!
  values(1:n) = 0.0D+00

  do i = 1, npoint
!
!  Set the next evaluation point X.
!
     x = ( real ( npoint - i,     kind = 8 ) * t(k) &
         + real (          i - 1, kind = 8 ) * t(n+1) ) &
         / real ( npoint     - 1, kind = 8 )
!
!  Locate X with respect to the knot array T.
!
     call interv ( t, n+1, x, left, mflag )
!
!  Get B(I,K)(X) in VALUES(1:N).  
!
!  K of these are supplied by BSPLVB:
!
!    B(LEFT-K+1,K)(X), 
!    ..., 
!    B(LEFT,K)(X).
!
!  All the others are zero.
!
     call bsplvb ( t, k, indx, x, left, values(left+1-k) )

     write ( *, '(f7.3, i8, 5f12.7)' ) x, left, values(3:7)
!
!  Zero out the values just computed in preparation for the next 
!  evaluation point.
!
    values(left+1-k:left) = 0.0D+00
 
  end do
 
  return
end
subroutine test05 ( )

!*****************************************************************************80
!
!! TEST05 plots the polynomials which make up a B spline.
!
!  Discussion:
!
!    This is example 2 from chapter X of the reference.
!
!  Modified:
!
!    13 February 2007
!
!  Author:
!
!    Carl deBoor
!
!  Reference:
!
!    Carl DeBoor,
!    A Practical Guide to Splines,
!    Springer, 2001,
!    ISBN: 0387953663,
!    LC: QA1.A647.v27.
!
  implicit none

  real ( kind = 8 ) biatx(4)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) left
  real ( kind = 8 ), dimension ( 11 ) :: t = (/ &
    0.0D+00, 0.0D+00, 0.0D+00, 0.0D+00, 1.0D+00, &
    3.0D+00, 4.0D+00, 6.0D+00, 6.0D+00, 6.0D+00, &
    6.0D+00 /)
  real ( kind = 8 ) values(4)
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST05'
  write ( *, '(a)' ) '  Polynomials that build splines.'
  write ( *, '(a)' ) ' '

  do i = 1, 40

    x = ( real ( 40 - i, kind = 8 ) * ( -1.0D+00 )   &
        + real (      i, kind = 8 ) * (  7.0D+00 ) ) &
        / real ( 40,     kind = 8 )

    do left = 4, 7
      call bsplvb ( t, 4, 1, x, left, biatx )
      values(left-3) = biatx(8-left)
    end do
!
!  According to the BSPLVB listing, BIATX now contains the value
!  at X of the polynomial which agrees on ( T(LEFT), T(LEFT+1) )
!  with the B-spline B(LEFT-4 + .,4,T).  
!
!  Hence, BIATX(8-LEFT) now contains the value of that polynomial for
!  B(LEFT-4+(8-LEFT),4,T), in other words, B(4,4,T).
!
!  Since this B-spline has support ( T(4), T(8) ), it makes sense to 
!  run LEFT = 4,...,7, storing the resulting values in VALUES(1:4),
!  for later printing.
!
    write ( *, '(2x,f10.1,2x,4f16.8)' ) x, values(1:4)

  end do

  return
end
subroutine test06 ( )

!*****************************************************************************80
!
!! TEST06 evaluates the piecewise polynomial form of a B spline.
!
!  Discussion:
!
!    This is example 3 from chapter X of the reference.
!
!  Modified:
!
!    16 February 2007
!
!  Author:
!
!    Carl deBoor
!
!  Reference:
!
!    Carl DeBoor,
!    A Practical Guide to Splines,
!    Springer, 2001,
!    ISBN: 0387953663,
!    LC: QA1.A647.v27.
!
  implicit none

  integer ( kind = 4 ), parameter :: k = 4
  integer ( kind = 4 ), parameter :: n = 7

  real ( kind = 8 ), dimension ( n) :: bcoef = (/ &
    0.0D+00, 0.0D+00, 0.0D+00, 1.0D+00, 0.0D+00, &
    0.0D+00, 0.0D+00 /)
  real ( kind = 8 ) break(5)
  real ( kind = 8 ) coef(4,4)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) l
  real ( kind = 8 ) ppvalu
  real ( kind = 8 ) scrtch(4,4)
  real ( kind = 8 ), dimension(n+k) :: t = (/ &
    0.0D+00, 0.0D+00, 0.0D+00, 0.0D+00, 1.0D+00, &
    3.0D+00, 4.0D+00, 6.0D+00, 6.0D+00, 6.0D+00, &
    6.0D+00 /)
  real ( kind = 8 ) value
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST06'
  write ( *, '(a)' ) '  Piecewise representation'
  write ( *, '(a)' ) ' '
!
!  Construct the piecewise polynomial representation.
!
  call bsplpp ( t, bcoef, n, k, scrtch, break, coef, l )
!
!  As a check, evaluate B(4,4,T) from its piecewise polynomial
!  representation on a fine mesh.
!
!  The values should agree with (some of) those generated in 
!  example 2 from chapter X.
!
  do i = 1, 40

    x = ( real ( 40 - i, kind = 8 ) * ( -1.0D+00 )   &
        + real (      i, kind = 8 ) * (  7.0D+00 ) ) &
        / real ( 40,     kind = 8 )

    value = ppvalu ( break, coef, l, 4, x, 0 )
    write ( *, '(g14.6,g14.6)' ) x, value

  end do

  return
end
subroutine test07 ( )

!*****************************************************************************80
!
!! TEST07 constructs a B spline using BVALUE.
!
!  Discussion:
!
!    This is example 4 from chapter X of the reference.
!
!  Modified:
!
!    13 February 2007
!
!  Author:
!
!    Carl deBoor
!
!  Reference:
!
!    Carl DeBoor,
!    A Practical Guide to Splines,
!    Springer, 2001,
!    ISBN: 0387953663,
!    LC: QA1.A647.v27.
!
  implicit none

  real ( kind = 8 ) bcoef(1)
  real ( kind = 8 ) bvalue
  integer ( kind = 4 ) i
  real ( kind = 8 ), dimension ( 5 ) :: t = (/ &
    0.0D+00, 1.0D+00, 3.0D+00, 4.0D+00, 6.0D+00 /)
  real ( kind = 8 ) value
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST07'
  write ( *, '(a)' ) '  Construct a spline via BVALUE'
  write ( *, '(a)' ) ' '

  bcoef(1) = 1.0D+00
!
!  Evaluate B(1,4,t) on a fine mesh.  On (0,6), the values should
!  coincide with those obtained in example 3.
!
  do i = 1, 40

    x = ( real ( 40 - i, kind = 8 ) * ( -1.0D+00 )   &
        + real (      i, kind = 8 ) * (  7.0D+00 ) ) &
        / real ( 40,     kind = 8 )

    value = bvalue ( t, bcoef, 1, 4, x, 0 )
    write ( *, '(g14.6,g14.6)' ) x, value

  end do

  return
end
subroutine test08 ( )

!*****************************************************************************80
!
!! TEST08 demonstrates cubic spline interpolation with good knots.
!
!  Discussion:
!
!    This is example 3 from chapter XII of the reference.
!
!  Modified:
!
!    13 February 2007
!
!  Author:
!
!    Carl deBoor
!
!  Reference:
!
!    Carl DeBoor,
!    A Practical Guide to Splines,
!    Springer, 2001,
!    ISBN: 0387953663,
!    LC: QA1.A647.v27.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 20

  real ( kind = 8 ) algerp
  real ( kind = 8 ) aloger
  real ( kind = 8 ) c(4,n_max)
  real ( kind = 8 ) decay
  real ( kind = 8 ) dx
  real ( kind = 8 ) errmax
  real ( kind = 8 ) g
  real ( kind = 8 ) h
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ibcbeg
  integer ( kind = 4 ) ibcend
  integer ( kind = 4 ) istep
  integer ( kind = 4 ) j
  integer ( kind = 4 ) n
  integer ( kind = 4 ) nhigh
  integer ( kind = 4 ) nlow
  real ( kind = 8 ) pnatx
  integer ( kind = 4 ) r
  integer ( kind = 4 ) run
  real ( kind = 8 ) step
  real ( kind = 8 ) tau(n_max)
  real ( kind = 8 ) x

  g(x) = sqrt ( x + 1.0D+00 )
 
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST08'
  write ( *, '(a)' ) '  Cubic spline with good knots'

  do run = 1, 3

    algerp = 0.0D+00
    decay = 0.0D+00
    istep = 20
    step = 20.0D+00

    if ( run == 1 ) then
      r = 8
      nlow = 4
      nhigh = 10
    else if ( run == 2 ) then
      r = 6
      nlow = 4
      nhigh = 20
    else if ( run == 3 ) then
      r = 4
      nlow = 4
      nhigh = 20
    end if

    write ( *, '(a)' ) ' '
    write ( *, '(a,i8)' ) '  Run number ', run
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  N   Max error   Decay exponent'
    write ( *, '(a)' ) ' '

    do n = nlow, nhigh, 2

      h = 1.0D+00 / real ( n - 1, kind = 8 )

      do i = 1, n
        tau(i) = &
          ( real ( ( n - 1 )**r - ( i - 1 )**r, kind = 8 ) * ( -1.0D+00 )   &
          + real (                ( i - 1 )**r, kind = 8 ) * ( +1.0D+00 ) ) &
          / real ( ( n - 1 )**r,                kind = 8 )
      end do

      do i = 1, n
        c(1,i) = g ( tau(i) )
      end do
!
!  Construct cubic spline interpolant.
!
      ibcbeg = 0
      ibcend = 0
      call cubspl ( tau, c, n, ibcbeg, ibcend )
!
!  Estimate maximum interpolation error on (-1,1).
!
      errmax = 0.0D+00

      do i = 1, n - 1

        dx = ( tau(i+1) - tau(i) ) / step

        do j = 1, istep

          h = real ( j, kind = 8 ) * dx

          pnatx = c(1,i) + h * ( &
                  c(2,i) + h * ( &
                  c(3,i) + h *   &
                  c(4,i) / 3.0D+00 ) / 2.0D+00 )

          x = tau(i) + h
          errmax = max ( errmax, abs ( g ( x ) - pnatx ) )

        end do

      end do

      aloger = log ( errmax )

      if ( nlow < n ) then
        decay = ( aloger - algerp ) &
              / log ( real ( n, kind = 8 ) / real ( n - 2, kind = 8 ) )
      end if

      algerp = aloger

      write ( *, '(i3,e12.4,f11.2)' ) n, errmax, decay

    end do

  end do

  return
end
subroutine test09 ( )

!*****************************************************************************80
!
!! TEST09 carries out cubic spline interpolation with good knots.
!
!  Discussion:
!
!    This is example 2 from chapter XII of the reference.
!
!  Modified:
!
!    13 February 2007
!
!  Author:
!
!    Carl deBoor
!
!  Reference:
!
!    Carl DeBoor,
!    A Practical Guide to Splines,
!    Springer, 2001,
!    ISBN: 0387953663,
!    LC: QA1.A647.v27.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 20

  real ( kind = 8 ) algerp
  real ( kind = 8 ) aloger
  real ( kind = 8 ) c(4,n_max)
  real ( kind = 8 ) decay
  real ( kind = 8 ) dx
  real ( kind = 8 ) errmax
  real ( kind = 8 ), external :: f02
  real ( kind = 8 ) h
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ibcbeg
  integer ( kind = 4 ) ibcend
  integer ( kind = 4 ) istep
  integer ( kind = 4 ) iter
  integer ( kind = 4 ) itermx
  integer ( kind = 4 ) j
  integer ( kind = 4 ) n
  integer ( kind = 4 ) nhigh
  integer ( kind = 4 ) nlow
  real ( kind = 8 ) pnatx
  integer ( kind = 4 ) run
  real ( kind = 8 ) scrtch(2,n_max)
  real ( kind = 8 ) step
  real ( kind = 8 ) tau(n_max)
  real ( kind = 8 ) taunew(n_max)
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST09'
  write ( *, '(a)' ) '  Cubic spline, good knots'

  do run = 1, 2

    algerp = 0.0D+00
    decay = 0.0D+00
    istep = 20
    step = 20.0D+00

    if ( run == 1 ) then
      itermx = 0
      nlow = 4
      nhigh = 20
    else
      itermx = 3
      nlow = 4
      nhigh = 20
    end if

    write ( *, '(a)' ) ' '
    write ( *, '(a,i8)' ) '  Run ', run
    write ( *, '(a)' ) ' '
    write ( *, '(a,i8)' ) '  ITERMX = ', itermx
    write ( *, '(a,i8)' ) '  NLOW =   ', nlow
    write ( *, '(a,i8)' ) '  NHIGH =  ', nhigh
    write ( *, '(a,i8,a)' ) '  Take ', itermx, ' cycles through NEWNOT'
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  N   Max error   Decay exponent'
    write ( *, '(a)' ) ' '
!
!  Loop over N = number of data points.
!
    do n = nlow, nhigh, 2
!
!  Knots are initially equispaced.
!
!  Note that the FORTRAN77 code uses what should be an equivalent
!  initialization of the knots TAU; however, at least for the
!  case with NLOW = 4 and ITERMX = 3, this results in very significant
!  difference in the output from NEWNOT, and hence a somewhat different
!  error calculation.  The differences between the F77 and F90 codes
!  seem to disappear after that.  This really seems to be something
!  about the way the algorithm is coded, not the particular language
!  used.
!
      h = 2.0D+00 / real ( n - 1, kind = 8 )

      do i = 1, n
        tau(i) = ( real ( n - i,     kind = 8 ) * ( - 1.0D+00 )   &
                 + real (     i - 1, kind = 8 ) * ( + 1.0D+00 ) ) &
                 / real ( n     - 1, kind = 8 )
      end do

      iter = 1
!
!  Construct the cubic spline interpolant. Then, ITERMX times,
!  determine new knots from it and find a new interpolant.
!
      do

        do i = 1, n
          c(1,i) = f02 ( tau(i) )
        end do

        ibcbeg = 0
        ibcend = 0
        call cubspl ( tau, c, n, ibcbeg, ibcend )

        if ( itermx < iter ) then
          exit
        end if

        iter = iter + 1
        call newnot ( tau, c, n-1, 4, taunew, n-1, scrtch )

        tau(1:n) = taunew(1:n)

      end do
!
!  Estimate the maximum interpolation error on (-1,1).
!
      errmax = 0.0D+00

      do i = 1, n - 1

        dx = ( tau(i+1) - tau(i) ) / step

        do j = 1, istep

           h = real ( j, kind = 8 ) * dx
           x = tau(i) + h

           pnatx = c(1,i) + h * ( &
                   c(2,i) + h * ( &
                   c(3,i) + h *   &
                   c(4,i) / 3.0D+00 ) / 2.0D+00 )

           errmax = max ( errmax, abs ( f02 ( x ) - pnatx ) )

        end do
      end do
!
!  Calculate the decay exponent.
!
      aloger = log ( errmax )

      if ( nlow < n ) then
        decay = ( aloger - algerp ) &
          / log ( real ( n, kind = 8 ) / real ( n - 2, kind = 8 ) )
      end if

      algerp = aloger
      write ( *, '(i3,e12.4,f11.2)' ) n, errmax, decay

    end do

  end do

  return
end
subroutine test10 ( )

!*****************************************************************************80
!
!! TEST10 demonstrates a quasi-interpolant with good knots.
!
!  Discussion:
!
!    This is example 4 from chapter XII of the reference.
!
!  Modified:
!
!    13 February 2007
!
!  Author:
!
!    Carl deBoor
!
!  Reference:
!
!    Carl DeBoor,
!    A Practical Guide to Splines,
!    Springer, 2001,
!    ISBN: 0387953663,
!    LC: QA1.A647.v27.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 20

  real ( kind = 8 ) algerp
  real ( kind = 8 ) aloger
  real ( kind = 8 ) bcoef(22)
  real ( kind = 8 ) break(n_max)
  real ( kind = 8 ) c(4,n_max)
  real ( kind = 8 ) decay
  real ( kind = 8 ) dg
  real ( kind = 8 ) ddg
  real ( kind = 8 ) dtip1
  real ( kind = 8 ) dtip2
  real ( kind = 8 ) dx
  real ( kind = 8 ) errmax
  real ( kind = 8 ) g
  real ( kind = 8 ) h
  integer ( kind = 4 ) i
  integer ( kind = 4 ) istep
  integer ( kind = 4 ) j
  integer ( kind = 4 ) l
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
  integer ( kind = 4 ) nlow
  integer ( kind = 4 ) nhigh
  real ( kind = 8 ) pnatx
  integer ( kind = 4 ) r
  integer ( kind = 4 ) run
  real ( kind = 8 ) scrtch(4,4)
  real ( kind = 8 ) step
  real ( kind = 8 ) t(26)
  real ( kind = 8 ) taui
  real ( kind = 8 ) x
!
!  G is the function to be approximated,
!  DG is its first, and
!  DDG its second derivative.
!
  g ( x ) = sqrt ( x + 1.0D+00 )
  dg ( x ) = 0.5D+00 / g ( x )
  ddg ( x ) = -0.5D+00 * dg ( x ) / ( x + 1.0D+00 )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST10'
  write ( *, '(a)' ) '  Quasi-interpolant'

  do run = 1, 2

    algerp = 0.0D+00
    decay = 0.0D+00
    istep = 20
    step = 20.0D+00

    if ( run == 1 ) then
      r = 8
      nlow = 4 
      nhigh = 10
    else if ( run == 2 ) then
      r = 6
      nlow = 4
      nhigh = 20
    end if

    write ( *, '(a)' ) ' '
    write ( *, '(a,i8)' ) '  Run ', run
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  N   Max error   Decay exponent'
    write ( *, '(a)' ) ' '
!
!  Loop over N = dim ( spline(4,T) ) - 2.
!
!  N is chosen as the parameter in order to afford 
!  comparison with examples in which cubic spline 
!  interpolation at N data points was used.
!
    do n = nlow, nhigh, 2

      h = 1.0D+00 / real ( n - 1, kind = 8 )
      m = n + 2
!
!  Interior knots are equidistributed with respect to the
!  function (X+1)**(1/R).
!
      t(1:4) = -1.0D+00

      do i = 5, m
        t(i) = &
          ( real ( ( m - 3 )**r - ( i - 4 )**r, kind = 8 ) * ( -1.0D+00 )   &
          + real (                ( i - 4 )**r, kind = 8 ) * ( +1.0D+00 ) ) &
          / real ( ( m - 3 )**r,                kind = 8 )
      end do

      t(m+1:m+4) = 1.0D+00
!
!  Construct quasi-interpolant.  
!  BCOEF(1) = G(-1.0) = 0.
!
      bcoef(1) = 0.0D+00
      dtip2 = t(5) - t(4)
      taui = t(5)
!
!  Special choice of TAU(2) to avoid infinite derivatives of G
!  at left endpoint.
!
      bcoef(2) = g ( taui ) - 2.0D+00 * dtip2 * dg ( taui ) / 3.0D+00 &
        + dtip2**2 * ddg ( taui ) / 6.0D+00

      do i = 3, m

        taui = t(i+2)
        dtip1 = dtip2
        dtip2 = t(i+3) - t(i+2)
!
!  Formula xii(30) of text is used.
!
        bcoef(i) = g ( taui ) + ( dtip2 - dtip1 ) * dg ( taui ) / 3.0D+00 &
          - dtip1 * dtip2 * ddg ( taui ) / 6.0D+00
      end do
!
!  Convert to piecewise polynomial representation.
!
      call bsplpp ( t, bcoef, m, 4, scrtch, break, c, l )
!
!  Estimate maximum interpolation error on (-1,1).
!
      errmax = 0.0D+00

      do i = 1, l

        dx = ( break(i+1) - break(i) ) / step

        do j = 1, istep

          h = real ( j, kind = 8 ) * dx

          pnatx = c(1,i) + h * ( &
                  c(2,i) + h * ( &
                  c(3,i) + h *   &
                  c(4,i) / 3.0D+00 ) / 2.0D+00 )

          errmax = max ( errmax, abs ( g(break(i)+h) - pnatx ) )

        end do

      end do
!
!  Calculate the decay exponent.
!
      aloger = log ( errmax )

      if ( nlow < n ) then
        decay = ( aloger - algerp ) &
              / log ( real ( n, kind = 8 ) / real ( n - 2, kind = 8 ) )
      end if

      algerp = aloger
      write ( *, '(i3,e12.4,f11.2)' ) n, errmax, decay

    end do

  end do

  return
end
subroutine test11 ( )

!*****************************************************************************80
!
!! TEST11 demonstrates that a large norm amplifies noise.
!
!  Discussion:
!
!    An initially uniform data point distribution of N points is
!    changed ITERMX times by moving the JCLOSE-th data point
!    toward its left neighbor, cutting the distance between the two 
!    by a factor of RATE each time.  
!
!    To demonstrate the corresponding increase in the norm of cubic 
!    spline interpolation at these data points, the data are taken 
!    from a cubic polynomial, which would be interpolated exactly,
!    but with noise of size SIZE added.  
!
!    The resulting noise or error in the interpolant, compared to the 
!    cubic, gives the norm or noise amplification factor and is 
!    printed out together with the diminishing distance H between the
!    two data points.
!
!    This is example 1 from chapter XIII of the reference.
!
!  Modified:
!
!    13 February 2007
!
!  Author:
!
!    Carl deBoor
!
!  Reference:
!
!    Carl DeBoor,
!    A Practical Guide to Splines,
!    Springer, 2001,
!    ISBN: 0387953663,
!    LC: QA1.A647.v27.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 200

  real ( kind = 8 ) amax
  real ( kind = 8 ) bcoef(n_max)
  real ( kind = 8 ) break(n_max)
  real ( kind = 8 ) coef(n_max*4)
  real ( kind = 8 ) dx
  real ( kind = 8 ) fx
  real ( kind = 8 ) g
  real ( kind = 8 ) gtau(n_max)
  real ( kind = 8 ) h
  integer ( kind = 4 ) i
  integer ( kind = 4 ) iflag
  integer ( kind = 4 ) istep
  integer ( kind = 4 ) iter
  integer ( kind = 4 ) itermx
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jclose
  integer ( kind = 4 ) l
  integer ( kind = 4 ) n
  real ( kind = 8 ) ppvalu
  real ( kind = 8 ) rate
  real ( kind = 8 ) round
  integer ( kind = 4 ) run
  real ( kind = 8 ) scrtch(n_max*7)
  real ( kind = 8 ) size
  real ( kind = 8 ) step
  real ( kind = 8 ) t(n_max+4)
  real ( kind = 8 ) tau(n_max)
  real ( kind = 8 ) x
!
!  The function to be interpolated.
!
  g ( x ) = 1.0D+00 + x + x**2 + x**3
 
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST11'
  write ( *, '(a)' ) '  A large norm amplifies noise.'
  write ( *, '(a)' ) ' '

  do run = 1, 2

    if ( run == 1 ) then
      n = 7
      itermx = 10
      jclose = 4
      size = 0.000001D+00
      rate = 2.0D+00
    else if ( run == 2 ) then
      n = 7
      itermx = 10
      jclose = 4
      size = 0.0D+00
      rate = 2.0D+00
    end if

    write ( *, '(a)' ) ' '
    write ( *, '(a,i8)' ) '  Run ', run
    write ( *, '(a,g14.6)' ) '  Size of noise  = ', size
    write ( *, '(a)' ) ' '

    istep = 20
    step = 20.0D+00

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '    H           Max error'
    write ( *, '(a)' ) ' '
!
!  Start with uniform data points.
!
    do i = 1, n
      tau(i) = real ( i - 1, kind = 8 ) &
             / real ( n - 1, kind = 8 )
    end do
!
!  Set up the knot sequence for not-a-knot end condition.
!  (Both TAU(2) and TAU(N-1) are "not a knot".)
!
    t(1:4) = tau(1)
    t(5:n) = tau(3:n-2)
    t(n+1:n+4) = tau(n)

    do iter = 1, itermx

      do i = 1, n
        gtau(i) = round ( g(tau(i)), size )
      end do

      call splint ( tau, gtau, t, n, 4, scrtch, bcoef, iflag )
 
      if ( iflag == 2 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) '  Error code IFLAG = 2 returned from SPLINT.'
        return
      end if
 
      call bsplpp ( t, bcoef, n, 4, scrtch, break, coef, l )
!
!  Calculate the maximum interpolation error.
!
      amax = 0.0D+00

      do i = 4, n

        dx = ( break(i-2) - break(i-3) ) / step

        do j = 2, istep

          x = break(i-2) - dx * real ( j - 1, kind = 8 )
          fx = ppvalu ( break, coef, l, 4, x, 0 )
          amax = max ( amax, abs ( fx - g ( x ) ) )

        end do
      end do

      h = tau(jclose) - tau(jclose-1)

      write ( *, '(e9.2,e15.3)' ) h, amax
!
!  Move TAU(JCLOSE) toward its left neighbor so as to cut
!  their distance by a factor of RATE.
!
      tau(jclose) = (          1.0D+00   * tau(jclose)     &
                    + ( rate - 1.0D+00 ) * tau(jclose-1) ) &
                    /   rate

    end do

  end do
 
  return
end
subroutine test12 ( )

!*****************************************************************************80
!
!! TEST12 shows a cubic spline interpolant at knot averages with good knots.
!
!  Discussion:
!
!    This is example 2 from chapter XIII of the reference.
!
!  Modified:
!
!    13 February 2007
!
!  Author:
!
!    Carl deBoor
!
!  Reference:
!
!    Carl DeBoor,
!    A Practical Guide to Splines,
!    Springer, 2001,
!    ISBN: 0387953663,
!    LC: QA1.A647.v27.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 20

  real ( kind = 8 ) algerp
  real ( kind = 8 ) aloger
  real ( kind = 8 ) bcoef(n_max+2)
  real ( kind = 8 ) break(n_max)
  real ( kind = 8 ) c(4,n_max)
  real ( kind = 8 ) decay
  real ( kind = 8 ) dx
  real ( kind = 8 ) errmax
  real ( kind = 8 ), external :: f02
  real ( kind = 8 ) gtau(n_max)
  real ( kind = 8 ) h
  integer ( kind = 4 ) i
  integer ( kind = 4 ) iflag
  integer ( kind = 4 ) istep
  integer ( kind = 4 ) iter
  integer ( kind = 4 ) itermx
  integer ( kind = 4 ) j
  integer ( kind = 4 ) l
  integer ( kind = 4 ) n
  integer ( kind = 4 ) nhigh
  integer ( kind = 4 ) nlow
  real ( kind = 8 ) pnatx
  integer ( kind = 4 ) run
  real ( kind = 8 ) scrtch(n_max*7)
  real ( kind = 8 ) step
  real ( kind = 8 ) t(n_max+6)
  real ( kind = 8 ) tau(n_max)
  real ( kind = 8 ) temp
  real ( kind = 8 ) tnew(n_max)
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST12'
  write ( *, '(a)' ) '  Interpolation at knot averages'

  do run = 1, 3

    algerp = 0.0D+00
    decay = 0.0D+00
    istep = 20
    step = 20.0D+00

    if ( run == 1 ) then
      itermx = 0
      nlow = 4
      nhigh = 20
    else if ( run == 2 ) then
      itermx = 3
      nlow = 4
      nhigh = 20
    else if ( run == 3 ) then
      itermx = 6
      nlow = 4
      nhigh = 20
    end if

    write ( *, '(a)' ) ' '
    write ( *, '(a,i8)' ) '  Run ', run
    write ( *, '(a,i8,a)' ) '  We will take ', itermx, ' cycles through NEWNOT.'
    write ( *, '(a,i8)' ) '  NLOW is  ', nlow
    write ( *, '(a,i8)' ) '  NHIGH is ', nhigh 
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  N   Max error   Decay exponent'
    write ( *, '(a)' ) ' '
!
!  Loop over N = number of data points.
!
    do n = nlow, nhigh, 2

      h = 2.0D+00 / real ( n - 3, kind = 8 )

      t(1:4) = -1.0D+00

      do i = 5, n
        t(i) = real ( i - 4, kind = 8 ) * h - 1.0D+00
      end do

      t(n+1:n+4) = 1.0D+00
!
!  Construct cubic spline interpolant.  Then, ITERMX times,
!  determine new knots from it and find a new interpolant.
!
      iter = 1

      do

        do i = 1, n
          tau(i) = sum ( t(i+1:i+3) ) / 3.0D+00
          gtau(i) = f02 ( tau(i) )
        end do

        call splint ( tau, gtau, t, n, 4, scrtch, bcoef, iflag )
 
        if ( iflag == 2 ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) '  Error code IFLAG = 2 from SPLINT.'
          return
        end if
 
        call bsplpp ( t, bcoef, n, 4, scrtch, break, c, l )

        if ( itermx < iter ) then
          exit
        end if

        iter = iter + 1
        call newnot ( break, c, l, 4, tnew, l, scrtch )

        t(5:3+l) = tnew(2:l)

      end do
!
!  Estimate maximum interpolation error on (-1,1).
!
      errmax = 0.0D+00

      do i = 1, l

        dx = ( break(i+1) - break(i) ) / step

        do j = 1, istep

          h = real ( j, kind = 8 ) * dx
          x = break(i) + h

          pnatx = c(1,i) + h * ( &
                  c(2,i) + h * ( &
                  c(3,i) + h *   &
                  c(4,i) / 3.0D+00 ) / 2.0D+00 )

          errmax = max ( errmax, abs ( f02 ( x ) - pnatx ) )

        end do
      end do
!
!  Calculate the decay exponent.
!
      aloger = log ( errmax )

      if ( nlow < n ) then
        temp = real ( n, kind = 8 ) / real ( n - 2, kind = 8 )
        decay = ( aloger - algerp ) / log ( temp )
      end if

      algerp = aloger

      write ( *, '(i3,e12.4,f11.2)' ) n, errmax, decay

    end do

  end do

  return
end
subroutine test13 ( )

!*****************************************************************************80
!
!! TEST13 is a cubic spline interpolant at knot averages with good knots.  
!  modified around label 4.
!
!  Discussion:
!
!    This is example 2m from chapter XIII of the reference.
!
!  Modified:
!
!    13 February 2007
!
!  Author:
!
!    Carl deBoor
!
!  Reference:
!
!    Carl DeBoor,
!    A Practical Guide to Splines,
!    Springer, 2001,
!    ISBN: 0387953663,
!    LC: QA1.A647.v27.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 20

  real ( kind = 8 ) algerp
  real ( kind = 8 ) aloger
  real ( kind = 8 ) bcoef(n_max+2)
  real ( kind = 8 ) break(n_max)
  real ( kind = 8 ) c(4,n_max)
  real ( kind = 8 ) decay
  real ( kind = 8 ) dx
  real ( kind = 8 ) errmax
  real ( kind = 8 ), external :: f02
  real ( kind = 8 ) gtau(n_max)
  real ( kind = 8 ) h
  integer ( kind = 4 ) i
  integer ( kind = 4 ) iflag
  integer ( kind = 4 ) istep
  integer ( kind = 4 ) iter
  integer ( kind = 4 ) itermx
  integer ( kind = 4 ) j
  integer ( kind = 4 ) l
  integer ( kind = 4 ) n
  integer ( kind = 4 ) nhigh
  integer ( kind = 4 ) nlow
  real ( kind = 8 ) pnatx
  integer ( kind = 4 ) run
  real ( kind = 8 ) scrtch(n_max*7)
  real ( kind = 8 ) step
  real ( kind = 8 ) t(n_max+6)
  real ( kind = 8 ) tau(n_max)
  real ( kind = 8 ) tnew(n_max)
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST13'
  write ( *, '(a)' ) '  Modified example 2'
  write ( *, '(a)' ) ' '

  do run = 1, 3

    algerp = 0.0D+00
    decay = 0.0D+00
    istep = 20
    step = 20.0D+00

    if ( run == 1 ) then
      itermx = 0
      nlow = 4
      nhigh = 20
    else if ( run == 2 ) then
      itermx = 1
      nlow = 4
      nhigh = 20
    else if ( run == 3 ) then
      itermx = 2
      nlow = 4
      nhigh = 20
    end if
 
    write ( *, '(a)' ) ' '
    write ( *, '(a,i8)' ) '  Run ', run
    write ( *, '(a,i8,a)' ) '  Take ', itermx, ' cycles through NEWNOT'
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  N   Max error   Decay exponent'
    write ( *, '(a)' ) ' '
!
!  Loop over number of data points.
!
    do n = nlow, nhigh, 2

      if ( nlow == n ) then

        h = 2.0D+00 / real ( n - 3, kind = 8 )

        t(1:4) = -1.0D+00

        do i = 5, n
          t(i) = real ( i - 4, kind = 8 ) * h - 1.0D+00
        end do

        t(n+1:n+4) = 1.0D+00

      else if ( nlow < n ) then

        call newnot ( break, c, l, 4, tnew, l+2, scrtch )

        l = l + 2

        t(5:l+3) = tnew(2:l)
        t(l+5) = 1.0D+00
        t(l+6) = 1.0D+00

      end if
!
!  Construct cubic spline interpolant. then,itermx  times,
!  determine new knots from it and find a new interpolant.
!
      iter = 1

      do

        do i = 1, n
          tau(i) = sum ( t(i+1:i+3) ) / 3.0D+00
          gtau(i) = f02 ( tau(i) )
        end do

        call splint ( tau, gtau, t, n, 4, scrtch, bcoef, iflag )

        call bsplpp ( t, bcoef, n, 4, scrtch, break, c, l )

        if ( itermx < iter ) then
          exit
        end if

        iter = iter + 1
        call newnot ( break, c, l, 4, tnew, l, scrtch )

        t(5:3+l) = tnew(2:l)

      end do
!
!  Estimate maximum interpolation error on (-1,1).
!
      errmax = 0.0D+00
      do i = 1, l
        dx = ( break(i+1) - break(i) ) / step

        do j = 1, istep

          h = real ( j, kind = 8 ) * dx
          x = break(i) + h

          pnatx = c(1,i) + h * ( &
                  c(2,i) + h * ( &
                  c(3,i) + h *   &
                  c(4,i) / 3.0D+00 ) / 2.0D+00 )

          errmax = max ( errmax, abs ( f02 ( x ) - pnatx ) )

        end do
      end do
!
!  Calculate the decay exponent.
!
      aloger = log ( errmax)

      if ( nlow < n ) then
        decay = ( aloger - algerp ) &
          / log ( real ( n, kind = 8 ) / real ( n - 2, kind = 8 ) )
      end if

      algerp = aloger

      write ( *, '(i3,e12.4,f11.2)' ) n, errmax, decay

    end do

  end do

  return
end
subroutine test14 ( )

!*****************************************************************************80
!
!! TEST14 tests optimal spline interpolation.
!
!  Discussion:
!
!    This is example 3 from chapter XIII of the reference.
!
!  Modified:
!
!    13 February 2007
!
!  Author:
!
!    Carl deBoor
!
!  Reference:
!
!    Carl DeBoor,
!    A Practical Guide to Splines,
!    Springer, 2001,
!    ISBN: 0387953663,
!    LC: QA1.A647.v27.
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 12
  integer ( kind = 4 ), parameter :: ntitan = 49
  integer ( kind = 4 ), parameter :: k = 5
  integer ( kind = 4 ), parameter :: lenscr = (n-k)*(2*k+3)+5*k+3

  real ( kind = 8 ) a(n)
  real ( kind = 8 ) bvalue
  real ( kind = 8 ) gtitan(ntitan)
  real ( kind = 8 ) gtau(ntitan)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) iflag
  integer ( kind = 4 ), dimension ( n ) :: ipick = (/ &
    1, 5, 11, 21, 27, 29, 31, 33, 35, 40, 45, 49 /)
  integer ( kind = 4 ) lx
  real ( kind = 8 ) scrtch(lenscr)
  real ( kind = 8 ) t(n+k)
  real ( kind = 8 ) tau(n)
  real ( kind = 8 ) x(ntitan)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST14'
  write ( *, '(a)' ) '  Optimal spline test'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  TITAND returns some test data for a material'
  write ( *, '(a)' ) '  property of titanium.'
!
!  Retrieve the titanium property values to be approximated.
!
  call titand ( x, gtitan, lx )
!
!  Select certain data points.
!
  tau(1:n) = x(ipick(1:n))
  gtau(1:n) = gtitan(ipick(1:n))
!
!  Evaluate the optimal spline interpolant.
!
  call splopt ( tau, n, k, scrtch, t, iflag )

  if ( 1 < iflag ) then
    return
  end if

  call splint ( tau, gtau, t, n, k, scrtch, a, iflag )

  if ( 1 < iflag ) then
    return
  end if

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  I  Data point  Data  Interpolant  Error'
  write ( *, '(a)' ) ' '

  do i = 1, lx
    gtau(i) = bvalue ( t, a, n, k, x(i), 0 )
    write ( *, '(i3,f8.0,f10.4,f9.4,e11.3)' ) &
      i, x(i), gtitan(i), gtau(i), gtitan(i) - gtau(i)
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Optimal knots:'
  write ( *, '(a)' ) ' '

  do i = 1, n - k
    write ( *, '(2x,i8,2x,g14.6)' ) i, t(k+i)
  end do

  return
end
subroutine test15 ( )

!*****************************************************************************80
!
!! TEST15 demonstrates the cubic smoothing spline.
!
!  Discussion:
!
!    Values from a cubic B spline are rounded to NDIGIT places
!    after the decimal point, then smoothed via SMOOTH for
!    various values of the control parameter S.
!
!    This is example 1 from chapter XIV of the reference.
!
!  Modified:
!
!    16 February 2007
!
!  Author:
!
!    Carl deBoor
!
!  Reference:
!
!    Carl DeBoor,
!    A Practical Guide to Splines,
!    Springer, 2001,
!    ISBN: 0387953663,
!    LC: QA1.A647.v27.
!
  implicit none

  integer ( kind = 4 ), parameter :: npoint = 61
  integer ( kind = 4 ), parameter :: ns = 7

  real ( kind = 8 ) a(npoint,4)
  real ( kind = 8 ), dimension ( 7 ) :: bcoef = (/ &
    0.0D+00, 0.0D+00, 0.0D+00, 1.0D+00, 0.0D+00, 0.0D+00, 0.0D+00 /)
  real ( kind = 8 ) break(5)
  real ( kind = 8 ) coef(4,4)
  real ( kind = 8 ) coefsm(4,60)
  real ( kind = 8 ) dely
  real ( kind = 8 ) dy(61)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) is
  integer ( kind = 4 ) j
  integer ( kind = 4 ) l
  integer ( kind = 4 ) ndigit
  real ( kind = 8 ) ppvalu
  real ( kind = 8 ) s(ns)
  real ( kind = 8 ) scrtch(427)
  real ( kind = 8 ) sfp
  real ( kind = 8 ) smooth
  real ( kind = 8 ), dimension ( 11 ) :: t = (/ &
    0.0D+00, 0.0D+00, 0.0D+00, 0.0D+00, 1.0D+00, &
    3.0D+00, 4.0D+00, 6.0D+00, 6.0D+00, 6.0D+00, &
    6.0D+00 /)
  real ( kind = 8 ) tenton
  real ( kind = 8 ) x(npoint)
  real ( kind = 8 ) y(npoint)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST15'
  write ( *, '(a)' ) '  SMOOTH defines a cubic smoothing spline'

  ndigit = 2

  do i = 1, 7
    s(i) = 6.0D+00 * 10.0D+00**(6-i)
  end do

  call bsplpp ( t, bcoef, 7, 4, scrtch, break, coef, l )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  The exact data values were rounded to ', ndigit
  write ( *, '(a)' ) '  digits after the decimal point.'

  tenton = 10.0D+00**ndigit
  dely = 0.5D+00 / tenton

  do i = 1, npoint
    x(i) = real ( i - 1, kind = 8 ) / 10.0D+00
    y(i) = ppvalu ( break, coef, l, 4, x(i), 0 )
  end do
!
!  Round the Y data.
!
  do i = 1, npoint
    y(i) = real ( int ( y(i) * tenton + 0.5D+00 ), kind = 8 ) / tenton
  end do
 
  do i = 1, npoint, 5
    do j = 1, 4
      a(i,j) = ppvalu ( break, coef, l, 4, x(i), j-1 )
    end do
  end do
 
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Value and derivatives of noisefree function ' &
    // 'at some points:'
  write ( *, '(a)' ) ' '

  do i = 1, npoint, 5
    write ( *, '(i4,4e15.7)' ) i, a(i,1:4)
  end do
!
!  Apply seven levels of smoothing S.
!
  dy(1:npoint) = dely

  do is = 1, ns

    sfp = smooth ( x, y, dy, npoint, s(is), scrtch, a )
 
    do i = 1, npoint - 1
      coefsm(1:4,i) = a(i,1:4)
    end do
 
    do i = 1, npoint, 5
      do j = 1, 4
        a(i,j) = ppvalu ( x, coefsm, npoint-1, 4, x(i), j-1 )
      end do
    end do
 
    write ( *, '(a)' ) ' '
    write ( *, '(a,g14.6)' ) '  Prescribed S =         ', s(is)
    write ( *, '(a,g14.6)' ) '  S(Smoothing spline)  = ', sfp
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Value and derivatives of smoothing spline at ' &
      // 'corresponding points:'
    write ( *, '(a)' ) ' '

    do i = 1, npoint, 5
      write ( *, '(i4,4e15.7)' ) i, a(i,1:4)
    end do

  end do
 
  return
end
subroutine test16 ( )

!*****************************************************************************80
!
!! TEST16 solves a second order boundary value problem.
!
!  Discussion:
!
!    Solution of a second order nonlinear two point boundary value
!    problem  on (0.,1.) by collocation with piecewise polynomial
!    functions having 4 pieces of order 6.  
!
!    Two passes through NEWNOT are to be made, without any knots
!    being added. 
!
!    Newton iteration is to be stopped when two iterates agree to 6  
!    decimal places.
!
!    This is an example from chapter XV of the reference.
!
!  Modified:
!
!    13 February 2007
!
!  Author:
!
!    Carl deBoor
!
!  Reference:
!
!    Carl DeBoor,
!    A Practical Guide to Splines,
!    Springer, 2001,
!    ISBN: 0387953663,
!    LC: QA1.A647.v27.
!
  implicit none

  real ( kind = 8 ) addbrk
  real ( kind = 8 ) aleft
  real ( kind = 8 ) aright
  integer ( kind = 4 ) iorder
  integer ( kind = 4 ) lbegin
  integer ( kind = 4 ) ntimes
  real ( kind = 8 ) relerr

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST16'
  write ( *, '(a)' ) '  Solution of second order two point'
  write ( *, '(a)' ) '  boundary value problem.'

  aleft = 0.0D+00
  aright = 1.0D+00
  lbegin = 4
  iorder = 6
  ntimes = 2
  addbrk = 0.0D+00
  relerr = 1.0D-06

  call colloc ( aleft, aright, lbegin, iorder, ntimes, addbrk, relerr )

  return
end
subroutine test17 ( )

!*****************************************************************************80
!
!! TEST17 demonstrates taut spline interpolation.
!
!  Discussion:
!
!    Taut spline interpolation to the titanium heat data.
!
!    This is example 2 from chapter XVI of the reference.
!
!  Modified:
!
!    13 February 2007
!
!  Author:
!
!    Carl deBoor
!
!  Reference:
!
!    Carl DeBoor,
!    A Practical Guide to Splines,
!    Springer, 2001,
!    ISBN: 0387953663,
!    LC: QA1.A647.v27.
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 12
  integer ( kind = 4 ), parameter :: npoint = 201

  real ( kind = 8 ) break(122)
  real ( kind = 8 ) coef(4,22)
  real ( kind = 8 ) :: gamma = 2.5D+00
  real ( kind = 8 ) gtau(49)
  real ( kind = 8 ) gtitan(49)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) iflag
  integer ( kind = 4 ), dimension ( n ) :: ipick = (/ &
    1, 5, 11, 21, 27, 29, 31, 33, 35, 40, 45, 49 /)
  integer ( kind = 4 ) ipicki
  integer ( kind = 4 ) :: k = 4
  integer ( kind = 4 ) l
  integer ( kind = 4 ) lx
  real ( kind = 8 ) plotf(npoint)
  real ( kind = 8 ) plotts(npoint)
  real ( kind = 8 ) plott(npoint)
  real ( kind = 8 ) ppvalu
  real ( kind = 8 ) scrtch(119)
  real ( kind = 8 ) step
  real ( kind = 8 ) tau(n)
  real ( kind = 8 ) x(49)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST17'
  write ( *, '(a)' ) '  Taut spline interpolation'

  call titand ( x, gtitan, lx )
 
  tau(1:n) = x(ipick(1:n))
  gtau(1:n) = gtitan(ipick(1:n))
 
  call tautsp ( tau, gtau, n, 0.0D+00, scrtch, break, coef, l, k, iflag )
 
  if ( 1 < iflag ) then
    return
  end if

  step = ( tau(n) - tau(1) ) / real ( npoint - 1, kind = 8 )
 
  do i = 1, npoint
    plott(i) = tau(1) + step * real ( i - 1, kind = 8 )
    plotf(i) = ppvalu ( break, coef, l, k, plott(i), 0 )
  end do
 
  if ( gamma < 0.0D+00 ) then
    return
  end if
 
  call tautsp ( tau, gtau, n, gamma, scrtch, break, coef, l, k, iflag )
 
  if ( 1 < iflag ) then
    return
  end if

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Cubic spline versus taut spline'
  write ( *, '(a,g14.6)' ) '  GAMMA = ', gamma
  write ( *, '(a)' ) ' '

  do i = 1, npoint

    plotts(i) = ppvalu ( break, coef, l, k, plott(i), 0 )

    write ( *, '(3g14.6)' ) plott(i), plotf(i), plotts(i)

  end do

  return
end
subroutine test18 ( )

!*****************************************************************************80
!
!! TEST18 demonstrates two parameterizations of some data.
!
!  Discussion:
!
!    This is example 3 from chapter XVI of the reference.
!
!  Modified:
!
!    13 February 2007
!
!  Author:
!
!    Carl deBoor
!
!  Reference:
!
!    Carl DeBoor,
!    A Practical Guide to Splines,
!    Springer, 2001,
!    ISBN: 0387953663,
!    LC: QA1.A647.v27.
!
  implicit none

  integer ( kind = 4 ), parameter :: k = 4
  integer ( kind = 4 ), parameter :: kpkm1 = k+k-1
  integer ( kind = 4 ), parameter :: n = 8
  integer ( kind = 4 ), parameter :: npiece = 6
  integer ( kind = 4 ), parameter :: npoint = 21

  real ( kind = 8 ) bcoef(n)
  real ( kind = 8 ) break(npiece)
  real ( kind = 8 ) ds
  integer ( kind = 4 ) i
  integer ( kind = 4 ) iflag
  integer ( kind = 4 ) l
  real ( kind = 8 ) ppvalu
  real ( kind = 8 ) q(n,kpkm1)
  integer ( kind = 4 ) run
  real ( kind = 8 ) s(n)
  real ( kind = 8 ) scrtch(k,k)
  real ( kind = 8 ) ss
  real ( kind = 8 ) t(n+k)
  real ( kind = 8 ), dimension ( n ) :: x = (/ &
    0.0D+00, 0.1D+00, 0.2D+00, 0.3D+00, 0.301D+00, &
    0.4D+00, 0.5D+00, 0.6D+00 /)
  real ( kind = 8 ) xcoef(k,npiece)
  real ( kind = 8 ) xx(npoint)
  real ( kind = 8 ) y(n)
  real ( kind = 8 ) ycoef(k,npiece)
  real ( kind = 8 ) yy(npoint)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST18'
  write ( *, '(a)' ) '  Two parameterizations'
!
!  Compute the Y component and set the natural parameterization.
!
  y(1:n) = ( x(1:n) - 0.3D+00 )**2
!
!  Convert X values to knots.  Note that second and second to
!  last X values are not knots.
!
  do run = 1, 2
 
    write ( *, '(a)' ) ' '

    if ( run == 1 ) then

      write ( *, '(a)' ) '  Using the "natural" parameterization.'

      s(1:n) = x(1:n)

    else if ( run == 2 ) then

      write ( *, '(a)' ) '  Using the "uniform" parameterization.'

      do i = 1, n
        s(i) = i
      end do

    end if

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '      X               Y'
    write ( *, '(a)' ) ' '

    t(1:k) = s(1) 
    t(k+1:n) = s(3:n-k+2)
    t(n+1:n+k) = s(n)
!
!  Interpolate to X component.
!
    call splint ( s, x, t, n, k, q, bcoef, iflag )

    call bsplpp ( t, bcoef, n, k, scrtch, break, xcoef, l )
!
!  Interpolate to Y component.  Since data abscissae and knots are
!  the same for both components, we only need to use 
!  backsubstitution
!
    bcoef(1:n) = y(1:n)
 
    call banslv ( q, kpkm1, n, k-1, k-1, bcoef )

    call bsplpp ( t, bcoef, n, k, scrtch, break, ycoef, l )
!
!  Evaluate the curve at some points near the potential trouble spot,
!  the fourth and fifth data points.
!
    ss = s(3)
    ds = ( s(6) - s(3) ) / real ( npoint - 1, kind = 8 )
 
    do i = 1, npoint

      xx(i) = ppvalu ( break, xcoef, l, k, ss, 0 )
      yy(i) = ppvalu ( break, ycoef, l, k, ss, 0 )
      ss = ss + ds

      write ( *, '(2x,g14.6,2x,g14.6)' ) xx(i), yy(i)

    end do
 
  end do

  return
end
subroutine test19 ( )

!*****************************************************************************80
!
!! TEST19 demonstrates bivariate spline interpolation.
!
!  Discussion:
!
!    A function G(X,Y) is interpolated on a rectangular mesh of data
!    points ( X(1:7), Y(1:6) ).
!
!    This is example 2 from chapter XVII of the reference.
!
!  Modified:
!
!    13 February 2007
!
!  Author:
!
!    Carl deBoor
!
!  Reference:
!
!    Carl DeBoor,
!    A Practical Guide to Splines,
!    Springer, 2001,
!    ISBN: 0387953663,
!    LC: QA1.A647.v27.
!
  implicit none

  integer ( kind = 4 ), parameter :: nx = 7
  integer ( kind = 4 ), parameter :: kx = 3
  integer ( kind = 4 ), parameter :: ny = 6
  integer ( kind = 4 ), parameter :: ky = 4
  integer ( kind = 4 ), parameter :: nq = (2*ky-1)*ny

  real ( kind = 8 ) bcoef(nx,ny)
  real ( kind = 8 ) bvalue
  real ( kind = 8 ), external :: f03
  real ( kind = 8 ) fdata(nx,ny)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) iflag
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jj
  integer ( kind = 4 ) lefty
  integer ( kind = 4 ) mflag
  real ( kind = 8 ) q(nq)
  real ( kind = 8 ) taux(nx)
  real ( kind = 8 ) tauy(ny)
  real ( kind = 8 ) tx(nx+kx)
  real ( kind = 8 ) ty(ny+ky)
  real ( kind = 8 ) work1(nx,ny)
  real ( kind = 8 ) work2(nx)
  real ( kind = 8 ) work3(nx,ny)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST19'
  write ( *, '(a)' ) '  Bivariate interpolation using B-splines.'
!
!  Set up the X data points TAUX and X spline knots TX.
!
  do i = 1, nx
    taux(i) = real ( i, kind = 8 )
  end do
 
  tx(1:kx) = taux(1)
 
  do i = kx+1, nx
    tx(i) = ( taux(i-kx+1) + taux(i-kx+2) ) / 2.0D+00
  end do

  tx(nx+1:nx+kx) = taux(nx)
!
!  Set up the Y data points TAUY and Y spline knots TY.
!
  do i = 1, ny
    tauy(i) = real ( i, kind = 8 )
  end do

  ty(1:ky) = tauy(1)

  do i = ky+1, ny
    ty(i) = tauy(i-ky+2)
  end do

  ty(ny+1:ny+ky) = tauy(ny)
!
!  Generate the function values at the mesh points.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The data to be interpolated:'
  write ( *, '(a)' ) ' '
  write ( *, '(6f12.1)' ) tauy(1:ny)
  write ( *, '(a)' ) ' '

  do i = 1, nx
    do j = 1, ny
      fdata(i,j) = f03 ( taux(i), tauy(j) )
    end do
    write ( *, '(f5.1,6e12.5)' ) taux(i), fdata(i,1:ny)
  end do
!
!  Construct the B spline coefficients of the interpolant.
!
  call spli2d ( taux, fdata, tx, nx, kx, ny, work2, q, work1, iflag )

  call spli2d ( tauy, work1, ty, ny, ky, nx, work2, q, bcoef, iflag )
!
!  Evaluate the interpolation error at the data points.
!
  do j = 1, ny

    call interv ( ty, ny+1, tauy(j), lefty, mflag )

    do i = 1, nx

      do jj = 1, ky
        work2(jj) = bvalue ( tx, bcoef(1,lefty-ky+jj), nx, kx, taux(i), 0 )
      end do

      work3(i,j) = bvalue ( ty(lefty-ky+1), work2, ky, ky, tauy(j), 0 )

    end do
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Interpolating function:'
  write ( *, '(a)' ) ' '
  write ( *, '(6f12.1)' ) tauy(1:ny)
  write ( *, '(a)' ) ' '

  do i = 1, nx
    write ( *, '(f5.1,6e12.5)' ) taux(i), work3(i,1:ny)
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Interpolation error:'
  write ( *, '(a)' ) ' '
  write ( *, '(6f12.1)' ) tauy(1:ny)
  write ( *, '(a)' ) ' '

  do i = 1, nx
    write ( *, '(f5.1,6e12.5)' ) taux(i), fdata(i,1:ny) - work3(i,1:ny)
  end do
 
  return
end
subroutine test20 ( )

!*****************************************************************************80
!
!! TEST20 demonstrates bivariate spline interpolation.
!
!  Discussion:
!
!    This is followed by conversion to piecewise polynomial representation 
!    and evaluation.
!
!    This is example 3 from chapter XVII of the reference.
!
!  Modified:
!
!    13 February 2007
!
!  Author:
!
!    Carl deBoor
!
!  Reference:
!
!    Carl DeBoor,
!    A Practical Guide to Splines,
!    Springer, 2001,
!    ISBN: 0387953663,
!    LC: QA1.A647.v27.
!
  implicit none

  integer ( kind = 4 ), parameter :: nx = 7
  integer ( kind = 4 ), parameter :: kx = 3
  integer ( kind = 4 ), parameter :: ny = 6
  integer ( kind = 4 ), parameter :: ky = 4
  integer ( kind = 4 ), parameter :: nq = (2*ky-1)*ny

  real ( kind = 8 ) bcoef(nx,ny)
  real ( kind = 8 ) breakx(6)
  real ( kind = 8 ) breaky(4)
  real ( kind = 8 ) coef(kx,5,ky,3)
  real ( kind = 8 ), external :: f03
  real ( kind = 8 ) fdata(nx,ny)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) iflag
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jj
  integer ( kind = 4 ) lefty
  integer ( kind = 4 ) lx
  integer ( kind = 4 ) ly
  integer ( kind = 4 ) mflag
  real ( kind = 8 ) ppvalu
  real ( kind = 8 ) q(nq)
  real ( kind = 8 ) taux(nx)
  real ( kind = 8 ) tauy(ny)
  real ( kind = 8 ) tx(nx+kx)
  real ( kind = 8 ) ty(ny+ky)
  real ( kind = 8 ) value
  real ( kind = 8 ) work1(nx,ny)
  real ( kind = 8 ) work2(nx)
  real ( kind = 8 ) work3(nx,ny)
  real ( kind = 8 ) work4(kx,kx,ny)
  real ( kind = 8 ) work5(ny,kx,5)
  real ( kind = 8 ) work6(ky,ky,21)
!
!  Note that, with the above parameters, lx = 5, ly = 3
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST20'
  write ( *, '(a)' ) '  Bivariate interpolation using piecewise polynomials.'
!
!  Set up data points and knots.
!
!  In X, interpolate between knots by parabolic splines, using
!  not-a-knot end condition
!
  do i = 1, nx
    taux(i) = real ( i, kind = 8 )
  end do

  tx(1:kx) = taux(1)

  do i = kx+1, nx
    tx(i) = ( taux(i-kx+1) + taux(i-kx+2) ) / 2.0D+00
  end do

  tx(nx+1:nx+kx) = taux(nx)
!
!  In Y, interpolate at knots by cubic splines, using not-a-knot
!  end condition.
!
  do i = 1, ny
    tauy(i) = i
  end do
 
  ty(1:ky) = tauy(1)
 
  do i = ky+1, ny
    ty(i) = tauy(i-ky+2)
  end do

  ty(ny+1:ny+ky) = tauy(ny)
!
!  Generate and print out function values.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The data to be interpolated:'
  write ( *, '(a)' ) ' '
  write ( *, '(6f12.1)' ) tauy(1:ny)
  write ( *, '(a)' ) ' '
 
  do i = 1, nx
    do j = 1, ny
      fdata(i,j) = f03 ( taux(i), tauy(j) )
    end do
    write ( *, '(f5.1,6e12.5)' ) taux(i), fdata(i,1:ny)
  end do
!
!  Construct the B spline coefficients of the interpolant.
!
  call spli2d ( taux, fdata, tx, nx, kx, ny, work2, q, work1, iflag )
 
  call spli2d ( tauy, work1, ty, ny, ky, nx, work2, q, bcoef, iflag )
!
!  Convert to piecewise polynomial representation.
!
  call bspp2d ( tx, bcoef, nx, kx, ny, work4, breakx, work5, lx )
 
  call bspp2d ( ty, work5, ny, ky, lx*kx, work6, breaky, coef, ly )
!
!  Evaluate the interpolation error at mesh points.
!
  do j = 1, ny

    call interv ( breaky, ly, tauy(j), lefty, mflag )

    do i = 1, nx
      do jj = 1, ky
        work2(jj) = ppvalu ( breakx, coef(1,1,jj,lefty), lx, kx, taux(i), 0 )
      end do

      work3(i,j) = ppvalu ( breaky(lefty), work2, 1, ky, tauy(j), 0 )

    end do

  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Interpolating function:'
  write ( *, '(a)' ) ' '
  write ( *, '(6f12.1)' ) tauy(1:ny)
  write ( *, '(a)' ) ' '

  do i = 1, nx
    write ( *, '(f5.1,6e12.5)' ) taux(i), work3(i,1:ny)
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Interpolation error:'
  write ( *, '(a)' ) ' '
  write ( *, '(6f12.1)' ) tauy(1:ny)
  write ( *, '(a)' ) ' '

  do i = 1, nx
    write ( *, '(f5.1,6e12.5)' ) taux(i), fdata(i,1:ny) - work3(i,1:ny)
  end do

  return
end
subroutine test21 ( )

!*****************************************************************************80
!
!! TEST21 demonstrates least-squares approximation by splines.
!
!  Discussion:
!
!    The program, though ostensibly written for L2-approximation, is 
!    typical for programs constructing a piecewise polynomial 
!    approximation to a function given in some sense.  The subprogram
!    L2APPR, for instance, could easily be replaced by one carrying 
!    out interpolation or some other form of approximation.
!
!    Input is expected in SETDAT, specifying both the data to
!    be approximated and the order and breakpoint sequence of the
!    piecewise polynomial approximating function to be used.  
!    Further, SETDAT is expected to terminate the run, for lack of 
!    further input or because ICOUNT has reached a critical value.
!
!    The number NTIMES is read in in the main program.  It 
!    specifies the number of passes through the knot improvement 
!    algorithm in NEWNOT to be made.  Also, ADDBRK is read in to 
!    specify that, on the average, ADDBRK knots are to be added per
!    pass through NEWNOT.  For example, ADDBRK = 0.34 would cause a knot 
!    to be added every third pass, as long as NTIMES < 50.
!
!    Printed output is governed by the three print control integers
!    PRBCO  = 1  gives printout of B spline coefficients of approximation
!    PRPCO  = 1  gives printout of piecewise polynomial representation
!                of the approximation.
!    PRFUN  = 1  gives printout of approximation and error at
!                every data point.
!    the order K, the number of pieces L, and the interior breakpoints
!    are always printed out as are (in L2ERR) the mean, mean square, and
!    maximum errors in the approximation.
!
!    ICOUNT provides communication with the data-input-and-
!    termination routine SETDAT.  It is initialized to 0 to
!    signal to SETDAT when it is being called for the first time.  After
!    that, it is up to SETDAT to use ICOUNT for keeping track of the
!    passes through SETDAT.
!
!    Information about the function to be approximated and order and
!    breakpoint sequence of the approximating piecewise polynomial 
!    functions is gathered by calling SETDAT.
!
!  Modified:
!
!    13 February 2007
!
!  Author:
!
!    Carl deBoor
!
!  Reference:
!
!    Carl DeBoor,
!    A Practical Guide to Splines,
!    Springer, 2001,
!    ISBN: 0387953663,
!    LC: QA1.A647.v27.
!
  implicit none

  integer ( kind = 4 ), parameter :: lpkmax = 100
  integer ( kind = 4 ), parameter :: ntmax = 200
  integer ( kind = 4 ), parameter :: ltkmax = 2000

  real ( kind = 8 ) addbrk
  real ( kind = 8 ) bcoef(lpkmax)
  real ( kind = 8 ) break
  real ( kind = 8 ) coef
  real ( kind = 8 ) gtau
  integer ( kind = 4 ) i
  integer ( kind = 4 ) icount
  integer ( kind = 4 ) ii
  integer ( kind = 4 ) inot
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) l
  integer ( kind = 4 ) lbegin
  integer ( kind = 4 ) ll
  integer ( kind = 4 ) lnew
  integer ( kind = 4 ) n
  integer ( kind = 4 ) nt
  integer ( kind = 4 ) ntau
  integer ( kind = 4 ) ntimes
  integer ( kind = 4 ) prbco
  integer ( kind = 4 ) prfun
  integer ( kind = 4 ) prpco
  real ( kind = 8 ) q(ltkmax)
  integer ( kind = 4 ) run
  real ( kind = 8 ) scrtch(ntmax)
  real ( kind = 8 ) t(ntmax)
  real ( kind = 8 ) tau
  real ( kind = 8 ) totalw
  real ( kind = 8 ) weight

  common / i4data / ntau
  common / r8data / tau(ntmax), gtau(ntmax), weight(ntmax), totalw
  common / approx / break(lpkmax), coef(ltkmax), l, k

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST21'
  write ( *, '(a)' ) '  Least squares approximation by splines.'

  do run = 1, 3

    write ( *, '(a)' ) ' '
    write ( *, '(a,i8)' ) '  Run number ', run

    icount = 0

    if ( run == 1 ) then
      ntimes = 1
      addbrk = 0.0D+00
      prbco = 0
      prpco = 0
      prfun = 0
      inot = 0
      call setdat1 ( icount )
    else if ( run == 2 ) then
      ntimes = 20
      addbrk = 1.0D+00
      prbco = 0
      prpco = 0
      prfun = 0
      inot = 1
      call setdat2 ( icount )
    else if ( run == 3 ) then
      ntimes = 4
      addbrk = 2.0D+00
      prbco = 0
      prpco = 0
      prfun = 0
      inot = 0
      call setdat3 ( icount )
    end if
!
!  Breakpoints are translated into knots, and the number N of
!  B splines to be used is obtained by calling L2KNTS.
!
    call l2knts ( break, l, k, t, n )
!
!  The integer NTIMES and the real ADDBRK are requested as well as 
!  the print controls  iprbco ,iprpco  and
!  iprfun.  ntimes  passes  are made through the routine new-
!  not, with an increase of  addbrk  knots for every pass.
!
    lbegin = l
    nt = 0
!
!  The B spline coefficients BCOEF of the L2-approximation are 
!  obtained.
!
    do

      call l2appr ( t, n, k, q, scrtch, bcoef )

      if ( prbco == 1 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) '  B spline coefficients:'
        write ( *, '(a)' ) ' '
        write ( *, '(2x,5g14.6)' ) bcoef(1:n)
      end if
!
!  Convert the B spline representation of the approximation to 
!  piecewise polynomial representation.
!
       call bsplpp ( t, bcoef, n, k, q, break, coef, l ) 

       write ( *, '(a)' ) ' '
       write ( *, '(a,i8)' ) '  Approximation by splines of order ',k
       write ( *, '(a,i8,a)' ) '  using ', l, ' intervals.'
       write ( *, '(a)' ) ' '
       write ( *, '(a)' ) '  Breakpoints:'
       write ( *, '(a)' ) ' '
       write ( *, '(2x,5g14.6)' ) break(2:l)

       if ( prpco == 1 ) then

         write ( *, '(a)' ) ' '
         write ( *, '(a)' ) '  Piecewise polynomial representation:'
         write ( *, '(a)' ) ' '

         do i = 1, l
           ii = ( i - 1 ) * k
           write ( *, '(f9.3,4e14.6/(9x,4e14.6))' ) &
             break(i), coef(ii+1:ii+k)
         end do

       end if
!
!  Compute and print out various error norms.
!
       call l2err ( prfun, scrtch, q )
!
!  If NEWNOT has been applied less than NTIMES times, try
!  it again to obtain, from the current approximation a possibly 
!  improved sequence of breakpoints with ADDBRK more breakpoints
!  (on the average) than the current approximation has.
!  If only an increase in breakpoints is wanted, without the
!  adjustment that NEWNOT provides, a fake NEWNOT routine could be
!  used here which merely returns the breakpoints for LNEW
!  equal intervals.
!
      if ( ntimes <= nt ) then
        exit
      end if

      lnew = lbegin + int ( real ( nt, kind = 8 ) * addbrk )

      if ( inot == 0 ) then
        call newnot ( break, coef, l, k, scrtch, lnew, t )
      else if ( inot == 1 ) then
        call evnnot ( break, coef, l, k, scrtch, lnew, t )
      end if

      call l2knts ( scrtch, lnew, k, t, n )
      nt = nt + 1

    end do

  end do

  return
end
subroutine setdat1 ( icount )

!*****************************************************************************80
!
!! SETDAT1 provides data for example 2 in chapter xiv. 
!
!  Discussion:
!
!    For a general purpose L2 approximation program, this routine
!    would have to be replaced by a subroutine reading in
!
!      NTAU, TAU(1:NTAU), GTAU(1:NTAU)
!
!    and reading in or setting
!
!      k, l, break(i),i = 1,...,l+1, and weight(i),i=1,...,ntau,
!
!    as well as totalw = sum(weight(i) ,i=1,...,ntau).
!
!    ICOUNT is equal to zero when setd is called in L2MAIN
!    for the first time.  After that, it is up to this routine to use ICOUNT
!    for keeping track of the number of calls.  This is important
!    since L2MAIN relies on this routine for termination.
!
!  Modified:
!
!    13 February 2007
!
!  Author:
!
!    Carl deBoor
!
!  Reference:
!
!    Carl DeBoor,
!    A Practical Guide to Splines,
!    Springer, 2001,
!    ISBN: 0387953663,
!    LC: QA1.A647.v27.
!
  implicit none

  integer ( kind = 4 ), parameter :: lpkmax = 100
  integer ( kind = 4 ), parameter :: ntmax = 200
  integer ( kind = 4 ), parameter :: ltkmax = 2000

  real ( kind = 8 ) break
  real ( kind = 8 ) coef
  real ( kind = 8 ) gtau
  integer ( kind = 4 ) i
  integer ( kind = 4 ) icount
  integer ( kind = 4 ) k
  integer ( kind = 4 ) l
  integer ( kind = 4 ) ntau
  real ( kind = 8 ) step
  real ( kind = 8 ) tau
  real ( kind = 8 ) totalw
  real ( kind = 8 ) weight

  common / i4data / ntau
  common / r8data / tau(ntmax), gtau(ntmax), weight(ntmax), totalw
  common / approx / break(lpkmax), coef(ltkmax), l, k

  icount = icount + 1
  ntau = 10

  do i = 1, ntau - 1
    tau(i) = 1.0D+00 - 0.5D+00**(i-1)
  end do
  tau(ntau) = 1.0D+00

  gtau(1:ntau) = tau(1:ntau)**2 + 1.0D+00

  weight(1:ntau) = 1.0D+00

  totalw = sum ( weight(1:ntau) )

  l = 6
  step = 1.0D+00 / real ( l, kind = 8 )
  k = 2

  do i = 1, l + 1
    break(i) = real ( i - 1, kind = 8 ) * step
  end do

  return
end
subroutine setdat2 ( icount )

!*****************************************************************************80
!
!! SETDAT2 provides data for example 3 in chapter xiv.
!
!  Modified:
!
!    14 February 2007
!
!  Author:
!
!    Carl deBoor
!
!  Reference:
!
!    Carl DeBoor,
!    A Practical Guide to Splines,
!    Springer, 2001,
!    ISBN: 0387953663,
!    LC: QA1.A647.v27.
!
  implicit none

  integer ( kind = 4 ), parameter :: lpkmax = 100
  integer ( kind = 4 ), parameter :: ntmax = 200
  integer ( kind = 4 ), parameter :: ltkmax = 2000

  real ( kind = 8 ) break
  real ( kind = 8 ) coef
  real ( kind = 8 ) gtau
  integer ( kind = 4 ) i
  integer ( kind = 4 ) icount
  integer ( kind = 4 ) k
  integer ( kind = 4 ) l
  integer ( kind = 4 ) ntau
  real ( kind = 8 ) roun
  real ( kind = 8 ) tau
  real ( kind = 8 ) totalw
  real ( kind = 8 ) weight
  real ( kind = 8 ) x

  common / i4data / ntau
  common / r8data / tau(ntmax), gtau(ntmax), weight(ntmax), totalw
  common / approx / break(lpkmax), coef(ltkmax), l, k

  roun ( x ) = real ( int ( x * 100.0D+00 ), kind = 8 ) / 100.0D+00

  icount = icount + 1
  ntau = 65

  do i = 1, ntau
    tau(i) = ( real ( ntau - i,     kind = 8 ) * 0.0D+00   &
             + real (        i - 1, kind = 8 ) * 3.0D+00 ) &
             / real ( ntau     - 1, kind = 8 )
  end do

  do i = 1, ntau
    gtau(i) = roun ( exp ( tau(i) ) )
  end do

  weight(1:ntau) = 1.0D+00

  totalw = sum ( weight(1:ntau) )

  l = 1
  break(1) = tau(1)
  break(2) = tau(ntau)
  k = 3

  return
end
subroutine setdat3 ( icount )

!*****************************************************************************80
!
!! SETDAT3 provides data for example 4 in chapter xiv.
!
!  Modified:
!
!    13 February 2007
!
!  Author:
!
!    Carl deBoor
!
!  Reference:
!
!    Carl DeBoor,
!    A Practical Guide to Splines,
!    Springer, 2001,
!    ISBN: 0387953663,
!    LC: QA1.A647.v27.
!
  implicit none

  integer ( kind = 4 ), parameter :: lpkmax = 100
  integer ( kind = 4 ), parameter :: n = 9
  integer ( kind = 4 ), parameter :: ntmax = 200
  integer ( kind = 4 ), parameter :: ltkmax = 2000

  real ( kind = 8 ) break
  real ( kind = 8 ), dimension ( n ) :: brkpic = (/ &
    595.000D+00,  730.985D+00,  794.414D+00,  844.476D+00,  880.060D+00, &
    907.814D+00,  938.001D+00,  976.752D+00, 1075.000D+00 /)
  real ( kind = 8 ) coef
  real ( kind = 8 ) gtau
  integer ( kind = 4 ) i
  integer ( kind = 4 ) icount
  integer ( kind = 4 ) k
  integer ( kind = 4 ) l
  integer ( kind = 4 ) ntau
  real ( kind = 8 ) tau
  real ( kind = 8 ) totalw
  real ( kind = 8 ) weight

  common / i4data / ntau
  common / r8data / tau(ntmax), gtau(ntmax), weight(ntmax), totalw
  common / approx / break(lpkmax), coef(ltkmax), l, k

  icount = icount + 1

  call titand ( tau, gtau, ntau )

  weight(1:ntau) = 1.0D+00

  totalw = sum ( weight(1:ntau) )

  l = n - 1
  k = 5

  break(1:n) = brkpic(1:n)

  return
end
subroutine test22 ( )

!*****************************************************************************80
!
!! TEST22 is the Runge example with cubic spline interpolation.
!
!  Modified:
!
!    13 February 2007
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 20

  real ( kind = 8 ) algerp
  real ( kind = 8 ) aloger
  real ( kind = 8 ) c(4,n_max)
  real ( kind = 8 ) decay
  real ( kind = 8 ) dx
  real ( kind = 8 ) errmax
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ibcbeg
  integer ( kind = 4 ) ibcend
  integer ( kind = 4 ) ii
  integer ( kind = 4 ) istep
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) n
  real ( kind = 8 ) pnatx
  real ( kind = 8 ) ppvalu
  real ( kind = 8 ), external :: runge
  real ( kind = 8 ), external :: rungep
  real ( kind = 8 ), external :: rungepp
  real ( kind = 8 ) tau(n_max)
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST22'
  write ( *, '(a)' ) '  The Runge example'
  write ( *, '(a)' ) '  Use cubic spline interpolation of order N.'
  write ( *, '(a)' ) ' '

  do ii = 0, 2

    ibcbeg = ii
    ibcend = ii

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Boundary conditions are '
    if ( ii == 0 ) then
      write ( *, '(a)' ) '  Not-a-knot.'
    else if ( ii == 1 ) then
      write ( *, '(a)' ) '  Derivatives at endpoints.'
    else if ( ii == 2 ) then
      write ( *, '(a)' ) '  Second derivatives at endpoints.'
    end if

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '   N   Max error   Decay exponent'
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) ' '

    algerp = 0.0D+00
    decay = 0.0D+00
    istep = 20

    do n = 2, n_max, 2
!
!  Choose interpolation points TAU equally spaced in (-1,1).
!
      do i = 1, n
        tau(i) = ( real ( n - i,     kind = 8 ) * ( -1.0D+00 ) &
                 + real (     i - 1, kind = 8 ) * ( +1.0D+00 ) ) &
                 / real ( n     - 1, kind = 8 )
      end do

      do i = 1, n
        c(1,i) = runge ( tau(i) )
      end do

      if ( ii == 0 ) then

      else if ( ii == 1 ) then
        c(2,1) = rungep ( tau(1) )
        c(2,n) = rungep ( tau(n) )
      else if ( ii == 2 ) then
        c(2,1) = rungepp ( tau(1) )
        c(2,2) = rungepp ( tau(n) )
      end if
!
!  Calculate the cubic spline.
!
      call cubspl ( tau, c, n, ibcbeg, ibcend )
!
!  Estimate the maximum interpolation error on (-1,1).
!
      errmax = 0.0D+00
!
!  Consider the subinterval ( TAU(I-1), TAU(I) ).
!
      do i = 2, n
!
!  Sample ISTEP points in the subinterval.
!
        do j = 1, istep

          x = ( real ( istep - j,     kind = 8 ) * tau(i-1) &
              + real (         j - 1, kind = 8 ) * tau(i) ) &
              / real ( istep     - 1, kind = 8 )

          pnatx = ppvalu ( tau, c, n-1, 4, x, 0 )

          errmax = max ( errmax, abs ( runge ( x ) - pnatx ) )

        end do
      end do

      aloger = log ( errmax )

      if ( 2 < n ) then
        decay = ( aloger - algerp ) &
          / log ( real ( n, kind = 8 ) / real ( n - 2, kind = 8 ) )
      end if

      algerp = aloger
      write ( *, '(i4,e12.4,f11.2)' ) n, errmax, decay

    end do

  end do

  return
end
function runge ( x )

!*****************************************************************************80
!
!! RUNGE evaluates the Runge function.
!
!  Modified:
!
!    13 February 2007
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument of the function.
!
!    Output, real ( kind = 8 ) RUNGE, the value of the function.
!
  implicit none

  real ( kind = 8 ) runge
  real ( kind = 8 ) x

  runge = 1.0D+00 / ( 1.0D+00 + 25.0D+00 * x * x )

  return
end
function rungep ( x )

!*****************************************************************************80
!
!! RUNGEP evaluates the derivative of the Runge function.
!
!  Modified:
!
!    13 February 2007
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument of the function.
!
!    Output, real ( kind = 8 ) RUNGEP, the value of the function.
!
  implicit none

  real ( kind = 8 ) rungep
  real ( kind = 8 ) x

  rungep = -50.0D+00 * x / ( 1.0D+00 + 25.0D+00 * x * x )**2

  return
end
function rungepp ( x )

!*****************************************************************************80
!
!! RUNGEPP evaluates the second derivative of the Runge function.
!
!  Modified:
!
!    13 February 2007
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument of the function.
!
!    Output, real ( kind = 8 ) RUNGEPP, the value of the function.
!
  implicit none

  real ( kind = 8 ) rungepp
  real ( kind = 8 ) u
  real ( kind = 8 ) up
  real ( kind = 8 ) v
  real ( kind = 8 ) vp
  real ( kind = 8 ) x

  u = -50.0D+00 * x
  up = -50.0D+00
  v = ( 1.0D+00 + 25.0D+00 * x * x )**2
  vp = 2.0D+00 * ( 1.0D+00 + 25.0D+00 * x * x ) * ( 50.0D+00 * x )

  rungepp = ( up * v - u * vp ) / v**2

  return
end
function f02 ( x )

!*****************************************************************************80
!
!! F02 evaluates the function SQRT(X+1).
!
!  Modified:
!
!    13 February 2007
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument of the function.
!
!    Output, real ( kind = 8 ) F02, the value of the function.
!
  implicit none

  real ( kind = 8 ) f02
  real ( kind = 8 ) x

  f02 = sqrt ( x + 1.0D+00 )

  return
end
function f03 ( x, y )

!*****************************************************************************80
!
!! F03 evaluates the function max(x-3.5,0.0)**2 + max(y-3.0,0.0)**3
!
!  Modified:
!
!    13 February 2007
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, Y, the arguments of the function.
!
!    Output, real ( kind = 8 ) F03, the value of the function.
!
  implicit none

  real ( kind = 8 ) f03
  real ( kind = 8 ) x
  real ( kind = 8 ) y

  f03 = max ( x - 3.5D+00, 0.0D+00 )**2 &
      + max ( y - 3.0D+00, 0.0D+00 )**3

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
