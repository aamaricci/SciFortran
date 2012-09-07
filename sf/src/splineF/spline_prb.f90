program main

!*****************************************************************************80
!
!! MAIN is the main program for SPLINE_PRB.
!
!  Discussion:
!
!    SPLINE_PRB calls the SPLINE tests.
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
  implicit none

  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'SPLINE_PRB'
  write ( *, '(a)' ) '  FORTRAN90 version:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Test the routines in the SPLINE library.'
 
  call test001
  call test002
  call test003
  call test004
  call test005
  call test006

  call test01
  call test02
  call test03
  call test04
  call test05
  call test06
  call test07
  call test08
  call test09
  call test10

  call test11
  call test115
  call test116
  call test12
  call test125
  call test126
  call test127
  call test13
  call test14
  call test143
  call test144
  call test145
  call test15
  call test16
  call test17
  call test18
  call test19
  call test20
  call test205

  call test21
  call test215
  call test22
  call test225
  call test23
  call test235
  call test24

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'SPLINE_PRB'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine test001

!***********************************************************************
!
!! TEST001 tests PARABOLA_VAL2.
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
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 1
  integer ( kind = 4 ), parameter :: ndata = 5

  integer ( kind = 4 ) i
  integer ( kind = 4 ) left
  real    ( kind = 8 ) xdata(ndata)
  real    ( kind = 8 ) xval
  real    ( kind = 8 ) ydata(dim_num,ndata)
  real    ( kind = 8 ) yval(dim_num)
  real    ( kind = 8 ) zdata(ndata)
  real    ( kind = 8 ) zval(dim_num)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST001'
  write ( *, '(a)' ) '  PARABOLA_VAL2 evaluates parabolas through'
  write ( *, '(a)' ) '    3 points in a table'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Our data tables will actually be parabolas:'
  write ( *, '(a)' ) '    Y: 2*x**2 + 3 * x + 1.'
  write ( *, '(a)' ) '    Z: 4*x**2 - 2 * x + 5.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '         I        X         Y           Z'
  write ( *, '(a)' ) ' '

  do i = 1, ndata

    xval = 2.0D+00 * real ( i, kind = 8 )
    xdata(i) = xval
    ydata(1,i) = 2.0D+00 * xval * xval + 3.0 * xval + 1.0D+00
    zdata(i) = 4.0D+00 * xval * xval - 2.0D+00 * xval + 5.0D+00
    write ( *, '(2x,i8,2x,f10.4,2x,f10.4,2x,f10.4)' ) &
      i, xdata(i), ydata(1,i), zdata(i)

  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Interpolated data:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '      LEFT        X         Y           Z'
  write ( *, '(a)' ) ' '

  do i = 1, 5

    xval = real ( 2 * i - 1, kind = 8 )
    left = i

    if ( ndata - 2 < left ) then
      left = ndata - 2
    end if

    if ( left < 1 ) then
      left = 1
    end if

    call parabola_val2 ( dim_num, ndata, xdata, ydata, left, xval, yval )

    call parabola_val2 ( dim_num, ndata, xdata, zdata, left, xval, zval )

    write ( *, '(2x,i8,2x,f10.4,2x,f10.4,2x,f10.4)' ) &
      left, xval, yval(1), zval(1)

  end do

  return
end
subroutine test002

!*****************************************************************************80
!
!! TEST002 tests R8VEC_BRACKET.
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
  implicit none

  integer ( kind = 4 ), parameter :: n = 10
  integer ( kind = 4 ), parameter :: test_num = 6

  integer ( kind = 4 ) left
  integer ( kind = 4 ) right
  integer ( kind = 4 ) test
  real    ( kind = 8 ) x(n)
  real    (  kind = 8 ), dimension ( test_num ) :: xtest = (/ &
   -10.0D+00, 1.0D+00, 4.5D+00, 5.0D+00, 10.0D+00, 12.0D+00 /)
  real    ( kind = 8 ) xval

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST002'
  write ( *, '(a)' ) '  R8VEC_BRACKET finds a pair of entries in a'
  write ( *, '(a)' ) '    sorted real array which bracket a value.'

  call r8vec_indicator ( n, x )
  x(6) = x(5)

  call r8vec_print ( n, x, '  Sorted array:' )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '    LEFT             RIGHT'
  write ( *, '(a)' ) '  X(LEFT)   XVAL   X(RIGHT)'
  write ( *, '(a)' ) ' '

  do test = 1, test_num

    xval = xtest(test)

    call r8vec_bracket ( n, x, xval, left, right )

    write ( *, '(2x,i14,14x,i14)' ) left, right

    if ( 1 <= left .and. 1 <= right ) then
      write ( *, '(2x,3g14.6)' ) x(left), xval, x(right)
    else if ( left < 1 .and. 1 <= right ) then
      write ( *, '(2x,14x,2g14.6)' )          xval, x(right)
    else if ( 1 <= left .and. right < 1 ) then
      write ( *, '(2x,2g14.6)' ) x(left), xval
    else if ( left < 1 .and. right < 1 ) then
      write ( *, '(2x,14x,g14.6)' )          xval
    end if

  end do

  return
end
subroutine test003

!*****************************************************************************80
!
!! TEST003 tests R8VEC_BRACKET3.
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
  implicit none

  integer ( kind = 4 ), parameter :: n = 10
  integer ( kind = 4 ), parameter :: test_num = 6

  integer ( kind = 4 ) left
  integer ( kind = 4 ) test
  real    ( kind = 8 ) x(n)
  real    (  kind = 8 ), dimension ( test_num ) :: xtest = (/ &
    -10.0D+00, 1.0D+00, 4.5D+00, 5.0D+00, 10.0D+00, 12.0D+00 /)
  real    ( kind = 8 ) xval

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST003'
  write ( *, '(a)' ) '  R8VEC_BRACKET3 finds a pair of entries in a'
  write ( *, '(a)' ) '    sorted real array which bracket a value.'

  call r8vec_indicator ( n, x )
  x(6) = x(5)

  call r8vec_print ( n, x, '  Sorted array:' )

  left = ( n + 1 ) / 2

  do test = 1, test_num

    xval = xtest(test)

    write ( *, '(a)' ) ' '
    write ( *, '(a,g14.6)' ) '  Search for XVAL = ', xval

    write ( *, '(a,i8)' ) '  Starting guess for interval is = ', left

    call r8vec_bracket3 ( n, x, xval, left )

    write ( *, '(a)' ) '  Nearest interval:'
    write ( *, '(2x,a,i8,a,g14.6)' ) '    X[', left,' ]= ', x(left)
    write ( *, '(2x,a,i8,a,g14.6)' ) '    X[', left+1, ' ]= ', x(left+1)

  end do

  return
end
subroutine test004

!*****************************************************************************80
!
!! TEST004 tests R8VEC_ORDER_TYPE.
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
  implicit none

  integer ( kind = 4 ), parameter :: n = 4
  integer ( kind = 4 ), parameter :: test_num = 6

  integer ( kind = 4 ) j
  integer ( kind = 4 ) order
  integer ( kind = 4 ) test
  real    ( kind = 8 ), dimension(n,test_num) :: x

  x(1,1) = 1.0D+00
  x(2,1) = 3.0D+00
  x(3,1) = 2.0D+00
  x(4,1) = 4.0D+00

  x(1,2) = 2.0D+00
  x(2,2) = 2.0D+00
  x(3,2) = 2.0D+00
  x(4,2) = 2.0D+00

  x(1,3) = 1.0D+00
  x(2,3) = 2.0D+00
  x(3,3) = 2.0D+00
  x(4,3) = 4.0D+00

  x(1,4) = 1.0D+00
  x(2,4) = 2.0D+00
  x(3,4) = 3.0D+00
  x(4,4) = 4.0D+00

  x(1,5) = 4.0D+00
  x(2,5) = 4.0D+00
  x(3,5) = 3.0D+00
  x(4,5) = 1.0D+00

  x(1,6) = 9.0D+00
  x(2,6) = 7.0D+00
  x(3,6) = 3.0D+00
  x(4,6) = 0.0D+00

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST004'
  write ( *, '(a)' ) '  R8VEC_ORDER_TYPE classifies a real vector as'
  write ( *, '(a)' ) '  -1: no order'
  write ( *, '(a)' ) '   0: all equal;'
  write ( *, '(a)' ) '   1: ascending;'
  write ( *, '(a)' ) '   2: strictly ascending;'
  write ( *, '(a)' ) '   3: descending;'
  write ( *, '(a)' ) '   4: strictly descending.'
  write ( *, '(a)' ) ' '

  do test = 1, test_num

    call r8vec_order_type ( n, x(1,test), order )

    write ( *, '(a)' ) ' '
    write ( *, '(a,i8)' ) '  The following vector has order type ', order
    write ( *, '(a)' ) ' '
    do j = 1, n
      write ( *, '(2x,i8,g14.6)' ) j, x(j,test)
    end do

  end do

  return
end
subroutine test005

!*****************************************************************************80
!
!! TEST005 tests R83_NP_FS.
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
  implicit none

  integer ( kind = 4 ), parameter :: n = 10

  real    ( kind = 8 ) a(3,n)
  real    ( kind = 8 ) b(n)
  integer ( kind = 4 ) :: seed = 123456789
  real    ( kind = 8 ) x(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST005'
  write ( *, '(a)' ) '  R83_NP_FS factors and solves a tridiagonal'
  write ( *, '(a)' ) '    linear system.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N = ', n
!
!  Set the matrix elements.
!
  call r83_uniform ( n, seed, a )
!
!  Set the desired solution.
!
  call r8vec_indicator ( n, x )
!
!  Compute b = A * x.
!
  call r83_mxv ( n, a, x, b )
!
!  Wipe out the solution.
!
  x(1:n) = 0.0D+00
!
!  Solve the system.
!
  call r83_np_fs ( n, a, b, x )

  call r8vec_print ( n, x, '  Solution:' )

  return
end
subroutine test006

!*****************************************************************************80
!
!! TEST006 tests DATA_TO_DIF and DIF_VAL.
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
  implicit none

  integer ( kind = 4 ), parameter :: maxtab = 8

  real    ( kind = 8 ) diftab(maxtab)
  real    ( kind = 8 ) err
  real    ( kind = 8 ) exact
  integer ( kind = 4 ) j
  integer ( kind = 4 ) ntab
  real    ( kind = 8 ) xtab(maxtab)
  real    ( kind = 8 ) xval
  real    ( kind = 8 ) ytab(maxtab)
  real    ( kind = 8 ) yval

  xval = 2.5D+00
  exact = exp ( xval )
 
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST006'
  write ( *, '(a,i8)' ) '  Approximate Y = EXP(X) using orders 1 to ', maxtab
  write ( *, '(a,g14.6)' ) '  Evaluate at X = ', xval
  write ( *, '(a,g14.6)' ) '  where EXP(X)=   ', exact

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '    Order  Approximate Y     Error'
  write ( *, '(a)' ) ' '
 
  do ntab = 1, maxtab
 
    do j = 1, ntab
      xtab(j) = real ( j - 1, kind = 8 )
      ytab(j) = exp ( xtab(j) )
    end do
 
    call data_to_dif ( ntab, xtab, ytab, diftab )

    call dif_val ( ntab, xtab, diftab, xval, yval )
 
    err = yval - exact
    write ( *, ' ( 2x, i8, 2g14.6 )' ) ntab, yval, err
 
  end do
 
  return
end
subroutine test01

!*****************************************************************************80
!
!! TEST01 tests BASIS_FUNCTION_B_VAL.
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
  implicit none

  integer ( kind = 4 ), parameter :: ndata = 5

  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jhi
  character mark
  integer ( kind = 4 ), parameter :: nsample = 4
  real    (  kind = 8 ), dimension ( ndata ) :: tdata = (/ &
    0.0D+00, 1.0D+00, 4.0D+00, 6.0D+00, 10.0D+00 /)
  real    ( kind = 8 ) thi
  real    ( kind = 8 ) tlo
  real    ( kind = 8 ) tval
  real    ( kind = 8 ) yval

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST01'
  write ( *, '(a)' ) '  BASIS_FUNCTION_B_VAL evaluates the '
  write ( *, '(a)' ) '    B spline basis function.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '    T            B(T)'
  write ( *, '(a)' ) ' '

  do i = 0, ndata

    if ( i == 0 ) then
      tlo = tdata(1) - 0.5D+00 * ( tdata(2) - tdata(1) )
      thi = tdata(1)
    else if ( i < ndata ) then
      tlo = tdata(i)
      thi = tdata(i+1)
    else if ( ndata <= i ) then
      tlo = tdata(ndata)
      thi = tdata(ndata) + 0.5D+00 * ( tdata(ndata) - tdata(ndata-1) )
    end if

    if ( i < ndata ) then
      jhi = nsample - 1
    else
      jhi = nsample
    end if

    do j = 0, jhi

      tval = ( real ( nsample - j, kind = 8 ) * tlo   &
             + real (           j, kind = 8 ) * thi ) &
             / real ( nsample,     kind = 8 )

      call basis_function_b_val ( tdata, tval, yval )

      if ( 0 < i .and. j == 0 ) then
        mark = '*'
      else
        mark = ' '
      end if

      write ( *, '(2x,a1,2g14.6)' ) mark, tval, yval

    end do

  end do

  return
end
subroutine test02

!*****************************************************************************80
!
!! TEST02 tests BASIS_FUNCTION_BETA_VAL.
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
  implicit none

  integer ( kind = 4 ), parameter :: ndata = 5

  real    ( kind = 8 ) beta1
  real    ( kind = 8 ) beta2
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jhi
  character mark
  integer ( kind = 4 ), parameter :: nsample = 4
  real    (  kind = 8 ), dimension ( ndata ) :: tdata = (/ &
    0.0D+00, 1.0D+00, 4.0D+00, 6.0D+00, 10.0D+00 /)
  real    ( kind = 8 ) thi
  real    ( kind = 8 ) tlo
  real    ( kind = 8 ) tval
  real    ( kind = 8 ) yval

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST02'
  write ( *, '(a)' ) '  BASIS_FUNCTION_BETA_VAL evaluates the '
  write ( *, '(a)' ) '    Beta spline basis function.'

  beta1 = 1.0D+00
  beta2 = 0.0D+00

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  BETA1 = ', beta1
  write ( *, '(a,g14.6)' ) '  BETA2 = ', beta2
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '       T            B(T)'
  write ( *, '(a)' ) ' '

  do i = 0, ndata

    if ( i == 0 ) then
      tlo = tdata(1) - 0.5D+00 * ( tdata(2) - tdata(1) )
      thi = tdata(1)
    else if ( i < ndata ) then
      tlo = tdata(i)
      thi = tdata(i+1)
    else if ( ndata <= i ) then
      tlo = tdata(ndata)
      thi = tdata(ndata) + 0.5D+00 * ( tdata(ndata) - tdata(ndata-1) )
    end if

    if ( i < ndata ) then
      jhi = nsample - 1
    else
      jhi = nsample
    end if

    do j = 0, jhi

      tval = ( real ( nsample - j, kind = 8 ) * tlo   &
             + real (           j, kind = 8 ) * thi ) &
             / real ( nsample,     kind = 8 )

      call basis_function_beta_val ( beta1, beta2, tdata, tval, yval )

      if ( 0 < i .and. j == 0 ) then
        mark = '*'
      else
        mark = ' '
      end if

      write ( *, '(2x,a1,2g14.6)' ) mark, tval, yval

    end do

  end do

  beta1 = 1.0D+00
  beta2 = 100.0D+00

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  BETA1 = ', beta1
  write ( *, '(a,g14.6)' ) '  BETA2 = ', beta2
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '       T            B(T)'
  write ( *, '(a)' ) ' '

  do i = 0, ndata

    if ( i == 0 ) then
      tlo = tdata(1) - 0.5D+00 * ( tdata(2) - tdata(1) )
      thi = tdata(1)
    else if ( i < ndata ) then
      tlo = tdata(i)
      thi = tdata(i+1)
    else if ( ndata <= i ) then
      tlo = tdata(ndata)
      thi = tdata(ndata) + 0.5D+00 * ( tdata(ndata) - tdata(ndata-1) )
    end if

    if ( i < ndata ) then
      jhi = nsample - 1
    else
      jhi = nsample
    end if

    do j = 0, jhi

      tval = ( real ( nsample - j, kind = 8 ) * tlo   &
             + real (           j, kind = 8 ) * thi ) &
             / real ( nsample,     kind = 8 )

      call basis_function_beta_val ( beta1, beta2, tdata, tval, yval )

      if ( 0 < i .and. j == 0 ) then
        mark = '*'
      else
        mark = ' '
      end if

      write ( *, '(2x,a1,2g14.6)' ) mark, tval, yval

    end do

  end do

  beta1 = 100.0D+00
  beta2 = 0.0D+00

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  BETA1 = ', beta1
  write ( *, '(a,g14.6)' ) '  BETA2 = ', beta2
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '       T            B(T)'
  write ( *, '(a)' ) ' '

  do i = 0, ndata

    if ( i == 0 ) then
      tlo = tdata(1) - 0.5D+00 * ( tdata(2) - tdata(1) )
      thi = tdata(1)
    else if ( i < ndata ) then
      tlo = tdata(i)
      thi = tdata(i+1)
    else if ( ndata <= i ) then
      tlo = tdata(ndata)
      thi = tdata(ndata) + 0.5D+00 * ( tdata(ndata) - tdata(ndata-1) )
    end if

    if ( i < ndata ) then
      jhi = nsample - 1
    else
      jhi = nsample
    end if

    do j = 0, jhi

      tval = ( real ( nsample - j, kind = 8 ) * tlo   &
             + real (           j, kind = 8 ) * thi ) &
             / real ( nsample,     kind = 8 )

      call basis_function_beta_val ( beta1, beta2, tdata, tval, yval )

      if ( 0 < i .and. j == 0 ) then
        mark = '*'
      else
        mark = ' '
      end if

      write ( *, '(2x,a1,2g14.6)' ) mark, tval, yval

    end do

  end do

  return
end
subroutine test03

!*****************************************************************************80
!
!! TEST03 tests BASIS_MATRIX_B_UNI and BASIS_MATRIX_TMP.
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
  implicit none

  integer ( kind = 4 ), parameter :: n = 4
  integer ( kind = 4 ), parameter :: ndata = 4

  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jhi
  integer ( kind = 4 ) left
  character mark
  real    ( kind = 8 ) mbasis(n,n)
  integer ( kind = 4 ), parameter :: nsample = 4
  real    (  kind = 8 ), dimension (ndata) ::  tdata = (/ &
    -1.0D+00, 0.0D+00, 1.0D+00, 2.0D+00 /)
  real    ( kind = 8 ) thi
  real    ( kind = 8 ) tlo
  real    ( kind = 8 ) tval
  real    (  kind = 8 ), dimension (ndata) ::  ydata = (/ &
    4.0D+00, 7.0D+00, 12.0D+00, 19.0D+00 /)
  real    ( kind = 8 ) yval

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST03'
  write ( *, '(a)' ) '  BASIS_MATRIX_B_UNI sets up the basis matrix'
  write ( *, '(a)' ) '    for the uniform B spline.'

  call basis_matrix_b_uni ( mbasis )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '       TDATA         YDATA'
  write ( *, '(a)' ) ' '
  do i = 1, ndata
    write ( *, '(2x,2g14.6)' ) tdata(i), ydata(i)
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '        T            Spline(T)'
  write ( *, '(a)' ) ' '

  left = 2

  do i = 0, ndata

    if ( i == 0 ) then
      tlo = tdata(1) - 0.5D+00 * ( tdata(2) - tdata(1) )
      thi = tdata(1)
    else if ( i < ndata ) then
      tlo = tdata(i)
      thi = tdata(i+1)
    else if ( ndata <= i ) then
      tlo = tdata(ndata)
      thi = tdata(ndata) + 0.5D+00 * ( tdata(ndata) - tdata(ndata-1) )
    end if

    if ( i < ndata ) then
      jhi = nsample - 1
    else
      jhi = nsample
    end if

    do j = 0, jhi

      tval = ( real ( nsample - j, kind = 8 ) * tlo   &
             + real (           j, kind = 8 ) * thi ) &
             / real ( nsample,     kind = 8 )

      call basis_matrix_tmp ( left, n, mbasis, ndata, tdata, ydata, tval, yval )

      if ( 0 < i .and. j == 0 ) then
        mark = '*'
      else
        mark = ' '
      end if

      write ( *, '(2x,a1,2g14.6)' ) mark, tval, yval

    end do

  end do

  return
end
subroutine test04

!*****************************************************************************80
!
!! TEST04 tests BASIS_MATRIX_BETA_UNI and BASIS_MATRIX_TMP.
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
  implicit none

  integer ( kind = 4 ), parameter :: n = 4
  integer ( kind = 4 ), parameter :: ndata = 4

  real    ( kind = 8 ) beta1
  real    ( kind = 8 ) beta2
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jhi
  integer ( kind = 4 ) left
  character mark
  real    ( kind = 8 ) mbasis(n,n)
  integer ( kind = 4 ), parameter :: nsample = 4
  real    ( kind = 8 ), dimension ( ndata ) :: tdata = (/ &
    -1.0D+00, 0.0D+00, 1.0D+00, 2.0D+00 /)
  real    ( kind = 8 ) thi
  real    ( kind = 8 ) tlo
  real    ( kind = 8 ) tval
  real    ( kind = 8 ), dimension ( ndata ) :: ydata = (/ &
    4.0D+00, 7.0D+00, 12.0D+00, 19.0D+00 /)
  real    ( kind = 8 ) yval

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST04'
  write ( *, '(a)' ) '  BASIS_MATRIX_BETA_UNI sets up the basis matrix'
  write ( *, '(a)' ) '    for the uniform beta spline.'
!
!  First test
!
  beta1 = 1.0D+00
  beta2 = 0.0D+00

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  BETA1 = ', beta1
  write ( *, '(a,g14.6)' ) '  BETA2 = ', beta2

  call basis_matrix_beta_uni ( beta1, beta2, mbasis )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '    TDATA, YDATA'
  write ( *, '(a)' ) ' '
  do i = 1, ndata
    write ( *, '(2g14.6)' ) tdata(i), ydata(i)
  end do

  left = 2

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '    T, Spline(T)'
  write ( *, '(a)' ) ' '

  do i = 0, ndata

    if ( i == 0 ) then
      tlo = tdata(1) - 0.5D+00 * ( tdata(2) - tdata(1) )
      thi = tdata(1)
    else if ( i < ndata ) then
      tlo = tdata(i)
      thi = tdata(i+1)
    else if ( ndata <= i ) then
      tlo = tdata(ndata)
      thi = tdata(ndata) + 0.5D+00 * ( tdata(ndata) - tdata(ndata-1) )
    end if

    if ( i < ndata ) then
      jhi = nsample - 1
    else
      jhi = nsample
    end if

    do j = 0, jhi

      tval = ( real ( nsample - j, kind = 8 ) * tlo   &
             + real (           j, kind = 8 ) * thi ) &
             / real ( nsample,     kind = 8 )

      call basis_matrix_tmp ( left, n, mbasis, ndata, tdata, ydata, tval, yval )

      if ( 0 < i .and. j == 0 ) then
        mark = '*'
      else
        mark = ' '
      end if

      write ( *, '(2x,a1,2g14.6)' ) mark, tval, yval

    end do

  end do
!
!  Second test
!
  beta1 = 1.0D+00
  beta2 = 100.0D+00

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  BETA1 = ', beta1
  write ( *, '(a,g14.6)' ) '  BETA2 = ', beta2

  call basis_matrix_beta_uni ( beta1, beta2, mbasis )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '    TDATA, YDATA'
  write ( *, '(a)' ) ' '
  do i = 1, ndata
    write ( *, '(2x,2g14.6)' ) tdata(i), ydata(i)
  end do

  left = 2

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '    T, Spline(T)'
  write ( *, '(a)' ) ' '

  do i = 0, ndata

    if ( i == 0 ) then
      tlo = tdata(1) - 0.5D+00 * ( tdata(2) - tdata(1) )
      thi = tdata(1)
    else if ( i < ndata ) then
      tlo = tdata(i)
      thi = tdata(i+1)
    else if ( ndata <= i ) then
      tlo = tdata(ndata)
      thi = tdata(ndata) + 0.5D+00 * ( tdata(ndata) - tdata(ndata-1) )
    end if

    if ( i < ndata ) then
      jhi = nsample - 1
    else
      jhi = nsample
    end if

    do j = 0, jhi

      tval = ( real ( nsample - j, kind = 8 ) * tlo   &
             + real (           j, kind = 8 ) * thi ) &
             / real ( nsample,     kind = 8 )

      call basis_matrix_tmp ( left, n, mbasis, ndata, tdata, ydata, tval, yval )

      if ( 0 < i .and. j == 0 ) then
        mark = '*'
      else
        mark = ' '
      end if

      write ( *, '(2x,a1,2g14.6)' ) mark, tval, yval

    end do

  end do
!
!  Third test
!
  beta1 = 100.0D+00
  beta2 = 0.0D+00

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  BETA1 = ', beta1
  write ( *, '(a,g14.6)' ) '  BETA2 = ', beta2

  call basis_matrix_beta_uni ( beta1, beta2, mbasis )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '     TDATA        YDATA'
  write ( *, '(a)' ) ' '
  do i = 1, ndata
    write ( *, '(2x,2g14.6)' ) tdata(i), ydata(i)
  end do

  left = 2

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '        T           Spline(T)'
  write ( *, '(a)' ) ' '

  do i = 0, ndata

    if ( i == 0 ) then
      tlo = tdata(1) - 0.5D+00 * ( tdata(2) - tdata(1) )
      thi = tdata(1)
    else if ( i < ndata ) then
      tlo = tdata(i)
      thi = tdata(i+1)
    else if ( ndata <= i ) then
      tlo = tdata(ndata)
      thi = tdata(ndata) + 0.5D+00 * ( tdata(ndata) - tdata(ndata-1) )
    end if

    if ( i < ndata ) then
      jhi = nsample - 1
    else
      jhi = nsample
    end if

    do j = 0, jhi

      tval = ( real ( nsample - j, kind = 8 ) * tlo   &
             + real (           j, kind = 8  ) * thi ) &
             / real ( nsample,     kind = 8  )

      call basis_matrix_tmp ( left, n, mbasis, ndata, tdata, ydata, tval, yval )

      if ( 0 < i .and. j == 0 ) then
        mark = '*'
      else
        mark = ' '
      end if

      write ( *, '(2x,a1,2g14.6)' ) mark, tval, yval

    end do

  end do

  return
end
subroutine test05

!*****************************************************************************80
!
!! TEST05 tests BASIS_MATRIX_BEZIER and BASIS_MATRIX_TMP.
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
  implicit none

  integer ( kind = 4 ), parameter :: n = 4
  integer ( kind = 4 ), parameter :: ndata = 4

  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jhi
  integer ( kind = 4 ) left
  character mark
  real    ( kind = 8 ) mbasis(n,n)
  integer ( kind = 4 ), parameter :: nsample = 4
  real    ( kind = 8 ), dimension ( ndata ) :: tdata = (/ &
    0.0D+00, 0.0D+00, 1.0D+00, 1.0D+00 /)
  real    ( kind = 8 ) thi
  real    ( kind = 8 ) tlo
  real    ( kind = 8 ) tval
  real    ( kind = 8 ), dimension ( ndata ) :: ydata = (/ &
    7.0D+00,  8.3333333D+00,   10.0D+00, 12.0D+00 /)
  real    ( kind = 8 ) yval

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST05'
  write ( *, '(a)' ) '  BASIS_MATRIX_BEZIER sets up the basis'
  write ( *, '(a)' ) '    matrix for the uniform Bezier spline.'

  call basis_matrix_bezier ( mbasis )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '       TDATA          YDATA'
  write ( *, '(a)' ) ' '
  do i = 1, ndata
    write ( *, '(2x,2g14.6)' ) tdata(i), ydata(i)
  end do

  left = 2

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '        T            Spline(T)'
  write ( *, '(a)' ) ' '

  do i = 0, ndata

    if ( i == 0 ) then
      tlo = tdata(1) - 0.5D+00 * ( tdata(2) - tdata(1) )
      thi = tdata(1)
    else if ( i < ndata ) then
      tlo = tdata(i)
      thi = tdata(i+1)
    else if ( ndata <= i ) then
      tlo = tdata(ndata)
      thi = tdata(ndata) + 0.5D+00 * ( tdata(ndata) - tdata(ndata-1) )
    end if

    if ( i < ndata ) then
      jhi = nsample - 1
    else
      jhi = nsample
    end if

    do j = 0, jhi

      tval = ( real ( nsample - j, kind = 8 ) * tlo   &
             + real (           j, kind = 8 ) * thi ) &
             / real ( nsample,     kind = 8 )

      call basis_matrix_tmp ( left, n, mbasis, ndata, tdata, ydata, tval, yval )

      if ( 0 < i .and. j == 0 ) then
        mark = '*'
      else
        mark = ' '
      end if

      write ( *, '(2x,a1,2g14.6)' ) mark, tval, yval

    end do

  end do

  return
end
subroutine test06

!*****************************************************************************80
!
!! TEST06 tests BASIS_MATRIX_HERMITE and BASIS_MATRIX_TMP.
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
  implicit none

  integer ( kind = 4 ), parameter :: n = 4
  integer ( kind = 4 ), parameter :: ndata = 4

  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jhi
  integer ( kind = 4 ) left
  character mark
  real    ( kind = 8 ) mbasis(n,n)
  integer ( kind = 4 ), parameter :: nsample = 4
  real    ( kind = 8 ), dimension ( ndata ) :: tdata = (/ &
    0.0D+00, 0.0D+00, 1.0D+00, 1.0D+00 /)
  real    ( kind = 8 ) thi
  real    ( kind = 8 ) tlo
  real    ( kind = 8 ) tval
  real    ( kind = 8 ), dimension ( ndata ) :: ydata = (/ &
    7.0D+00, 12.0D+00, 4.0D+00, 6.0D+00 /)
  real    ( kind = 8 ) yval

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST06'
  write ( *, '(a)' ) '  BASIS_MATRIX_HERMITE sets up the basis matrix'
  write ( *, '(a)' ) '    for the Hermite spline.'

  call basis_matrix_hermite ( mbasis )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '       TDATA        YDATA'
  write ( *, '(a)' ) ' '
  do i = 1, ndata
    write ( *, '(2x,2g14.6)' ) tdata(i), ydata(i)
  end do

  left = 2

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '        T           Spline(T)'
  write ( *, '(a)' ) ' '

  do i = 0, ndata

    if ( i == 0 ) then
      tlo = tdata(1) - 0.5D+00 * ( tdata(2) - tdata(1) )
      thi = tdata(1)
    else if ( i < ndata ) then
      tlo = tdata(i)
      thi = tdata(i+1)
    else if ( ndata <= i ) then
      tlo = tdata(ndata)
      thi = tdata(ndata) + 0.5D+00 * ( tdata(ndata) - tdata(ndata-1) )
    end if

    if ( i < ndata ) then
      jhi = nsample - 1
    else
      jhi = nsample
    end if

    do j = 0, jhi

      tval = ( real ( nsample - j, kind = 8 ) * tlo   &
             + real (           j, kind = 8 ) * thi ) &
             / real ( nsample,     kind = 8 )

      call basis_matrix_tmp ( left, n, mbasis, ndata, tdata, ydata, tval, yval )

      if ( 0 < i .and. j == 0 ) then
        mark = '*'
      else
        mark = ' '
      end if

      write ( *, '(2x,a1,2g14.6)' ) mark, tval, yval

    end do

  end do

  return
end
subroutine test07

!*****************************************************************************80
!
!! TEST07 tests BASIS_MATRIX_OVERHAUSER_UNI and BASIS_MATRIX_TMP.
!
!  Discussion:
!
!   YDATA(1:NDATA) = ( TDATA(1:NDATA) + 2 )**2 + 3
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
  implicit none

  integer ( kind = 4 ), parameter :: n = 4
  integer ( kind = 4 ), parameter :: ndata = 4

  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jhi
  integer ( kind = 4 ) left
  character mark
  real    ( kind = 8 ) mbasis(n,n)
  integer ( kind = 4 ), parameter :: nsample = 4
  real    ( kind = 8 ), dimension ( ndata ) :: tdata = (/ &
    -1.0D+00, 0.0D+00, 1.0D+00, 2.0D+00 /)
  real    ( kind = 8 ) thi
  real    ( kind = 8 ) tlo
  real    ( kind = 8 ) tval
  real    ( kind = 8 ), dimension ( ndata ) :: ydata = (/ &
    4.0D+00, 7.0D+00, 12.0D+00, 19.0D+00 /)
  real    ( kind = 8 ) yval

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST07'
  write ( *, '(a)' ) '  BASIS_MATRIX_OVERHAUSER_UNI sets up the basis'
  write ( *, '(a)' ) '    matrix for the uniform Overhauser spline.'

  call basis_matrix_overhauser_uni ( mbasis )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '       TDATA         YDATA'
  write ( *, '(a)' ) ' '
  do i = 1, ndata
    write ( *, '(2x,2g14.6)' ) tdata(i), ydata(i)
  end do

  left = 2

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '        T            Spline(T)'
  write ( *, '(a)' ) ' '

  do i = 0, ndata

    if ( i == 0 ) then
      tlo = tdata(1) - 0.5D+00 * ( tdata(2) - tdata(1) )
      thi = tdata(1)
    else if ( i < ndata ) then
      tlo = tdata(i)
      thi = tdata(i+1)
    else if ( ndata <= i ) then
      tlo = tdata(ndata)
      thi = tdata(ndata) + 0.5D+00 * ( tdata(ndata) - tdata(ndata-1) )
    end if

    if ( i < ndata ) then
      jhi = nsample - 1
    else
      jhi = nsample
    end if

    do j = 0, jhi

      tval = ( real ( nsample - j, kind = 8 ) * tlo   &
             + real (           j, kind = 8 ) * thi ) &
             / real ( nsample,     kind = 8 )

      call basis_matrix_tmp ( left, n, mbasis, ndata, tdata, ydata, tval, yval )

      if ( 0 < i .and. j == 0 ) then
        mark = '*'
      else
        mark = ' '
      end if

      write ( *, '(2x,a1,2g14.6)' ) mark, tval, yval

    end do

  end do

  return
end
subroutine test08

!*****************************************************************************80
!
!! TEST08 tests BASIS_MATRIX_OVERHAUSER_NONUNI and BASIS_MATRIX_TMP.
!
!  Discussion:
!
!    YDATA(1:NDATA) = ( TDATA(1:NDATA) - 2 )**2 + 3
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
  implicit none

  integer ( kind = 4 ), parameter :: n = 4
  integer ( kind = 4 ), parameter :: ndata = 4

  real    ( kind = 8 ) alpha
  real    ( kind = 8 ) beta
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jhi
  integer ( kind = 4 ) left
  character mark
  real    ( kind = 8 ) mbasis(n,n)
  integer ( kind = 4 ), parameter :: nsample = 4
  real    ( kind = 8 ), dimension ( ndata ) :: tdata
  real    ( kind = 8 ) thi
  real    ( kind = 8 ) tlo
  real    ( kind = 8 ) tval
  real    (  kind = 8 ), dimension ( ndata ) :: ydata
  real    ( kind = 8 ) yval

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST08'
  write ( *, '(a)' ) '  BASIS_MATRIX_OVERHAUSER_NONUNI sets up the'
  write ( *, '(a)' ) '    basis matrix for the nonuniform Overhauser'
  write ( *, '(a)' ) '    spline.'

  tdata = (/ 0.0D+00, 1.0D+00, 2.0D+00, 3.0D+00 /)

  alpha = ( tdata(3) - tdata(2) ) / ( tdata(3) - tdata(1) )
  beta = ( tdata(3) - tdata(2) ) / ( tdata(4) - tdata(2) )

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  ALPHA = ', alpha
  write ( *, '(a,g14.6)' ) '  BETA = ', beta

  call basis_matrix_overhauser_nonuni ( alpha, beta, mbasis )

  do i = 1, ndata
    ydata(i) = ( tdata(i) - 2.0D+00 )**2 + 3.0D+00
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '       TDATA         YDATA'
  write ( *, '(a)' ) ' '
  do i = 1, ndata
    write ( *, '(2x,2g14.6)' ) tdata(i), ydata(i)
  end do

  left = 2

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '       T           Spline(T)'
  write ( *, '(a)' ) ' '

  do i = 0, ndata

    if ( i == 0 ) then
      tlo = tdata(1) - 0.5D+00 * ( tdata(2) - tdata(1) )
      thi = tdata(1)
    else if ( i < ndata ) then
      tlo = tdata(i)
      thi = tdata(i+1)
    else if ( ndata <= i ) then
      tlo = tdata(ndata)
      thi = tdata(ndata) + 0.5D+00 * ( tdata(ndata) - tdata(ndata-1) )
    end if

    if ( i < ndata ) then
      jhi = nsample - 1
    else
      jhi = nsample
    end if

    do j = 0, jhi

      tval = ( real ( nsample - j, kind = 8 ) * tlo   &
             + real (           j, kind = 8 ) * thi ) &
             / real ( nsample,     kind = 8 )

      call basis_matrix_tmp ( left, n, mbasis, ndata, tdata, ydata, tval, yval )

      if ( 0 < i .and. j == 0 ) then
        mark = '*'
      else
        mark = ' '
      end if

      write ( *, '(2x,a1,2g14.6)' ) mark, tval, yval

    end do

  end do

  tdata(1:4) = (/ 0.0D+00, 1.0D+00, 2.0D+00, 5.0D+00 /)

  alpha = ( tdata(3) - tdata(2) ) / ( tdata(3) - tdata(1) )
  beta = ( tdata(3) - tdata(2) ) / ( tdata(4) - tdata(2) )

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  ALPHA = ', alpha
  write ( *, '(a,g14.6)' ) '  BETA = ', beta

  call basis_matrix_overhauser_nonuni ( alpha, beta, mbasis )

  do i = 1, ndata
    ydata(i) = ( tdata(i) - 2.0D+00 )**2 + 3.0D+00
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '    TDATA, YDATA'
  write ( *, '(a)' ) ' '
  do i = 1, ndata
    write ( *, '(2x,2g14.6)' ) tdata(i), ydata(i)
  end do

  left = 2

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '    T, Spline(T)'
  write ( *, '(a)' ) ' '

  do i = 0, ndata

    if ( i == 0 ) then
      tlo = tdata(1) - 0.5D+00 * ( tdata(2) - tdata(1) )
      thi = tdata(1)
    else if ( i < ndata ) then
      tlo = tdata(i)
      thi = tdata(i+1)
    else if ( ndata <= i ) then
      tlo = tdata(ndata)
      thi = tdata(ndata) + 0.5D+00 * ( tdata(ndata) - tdata(ndata-1) )
    end if

    if ( i < ndata ) then
      jhi = nsample - 1
    else
      jhi = nsample
    end if

    do j = 0, jhi

      tval = ( real ( nsample - j, kind = 8 ) * tlo   &
             + real (           j, kind = 8 ) * thi ) &
             / real ( nsample,     kind = 8 )

      call basis_matrix_tmp ( left, n, mbasis, ndata, tdata, ydata, tval, yval )

      if ( 0 < i .and. j == 0 ) then
        mark = '*'
      else
        mark = ' '
      end if

      write ( *, '(2x,a1,2g14.6)' ) mark, tval, yval

    end do

  end do

  tdata(1:4) = (/ 0.0D+00, 3.0D+00, 4.0D+00, 5.0D+00 /)

  alpha = ( tdata(3) - tdata(2) ) / ( tdata(3) - tdata(1) )
  beta = ( tdata(3) - tdata(2) ) / ( tdata(4) - tdata(2) )

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  ALPHA = ', alpha
  write ( *, '(a,g14.6)' ) '  BETA = ', beta

  call basis_matrix_overhauser_nonuni ( alpha, beta, mbasis )

  do i = 1, ndata
    ydata(i) = ( tdata(i) - 2.0D+00 )**2 + 3.0D+00
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '    TDATA, YDATA'
  write ( *, '(a)' ) ' '
  do i = 1, ndata
    write ( *, '(2x,2g14.6)' ) tdata(i), ydata(i)
  end do

  left = 2

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '    T, Spline(T)'
  write ( *, '(a)' ) ' '

  do i = 0, ndata

    if ( i == 0 ) then
      tlo = tdata(1) - 0.5D+00 * ( tdata(2) - tdata(1) )
      thi = tdata(1)
    else if ( i < ndata ) then
      tlo = tdata(i)
      thi = tdata(i+1)
    else if ( ndata <= i ) then
      tlo = tdata(ndata)
      thi = tdata(ndata) + 0.5D+00 * ( tdata(ndata) - tdata(ndata-1) )
    end if

    if ( i < ndata ) then
      jhi = nsample - 1
    else
      jhi = nsample
    end if

    do j = 0, jhi

      tval = ( real ( nsample - j, kind = 8 ) * tlo   &
             + real (           j, kind = 8 ) * thi ) &
             / real ( nsample,     kind = 8 )

      call basis_matrix_tmp ( left, n, mbasis, ndata, tdata, ydata, tval, yval )

      if ( 0 < i .and. j == 0 ) then
        mark = '*'
      else
        mark = ' '
      end if

      write ( *, '(2x,a1,2g14.6)' ) mark, tval, yval

    end do

  end do

  return
end
subroutine test09

!*****************************************************************************80
!
!! TEST09 tests BASIS_MATRIX_OVERHAUSER_NONUNI and BASIS_MATRIX_TMP.
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
  implicit none

  integer ( kind = 4 ), parameter :: n = 4
  integer ( kind = 4 ), parameter :: ndata = 4

  real    ( kind = 8 ) alpha
  real    ( kind = 8 ) beta
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jhi
  integer ( kind = 4 ) left
  character mark
  real    ( kind = 8 ) mbasis(n,n)
  integer ( kind = 4 ), parameter :: nsample = 4
  real    ( kind = 8 ) tdata(ndata)
  real    ( kind = 8 ) thi
  real    ( kind = 8 ) tlo
  real    ( kind = 8 ) tval
  real    ( kind = 8 ) ydata(ndata)
  real    ( kind = 8 ) yval

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST09'
  write ( *, '(a)' ) '  BASIS_MATRIX_OVERHAUSER_NONUNI sets up the'
  write ( *, '(a)' ) '    basis matrix for the nonuniform Overhauser '
  write ( *, '(a)' ) '    spline.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  First test that the nonuniform code can match'
  write ( *, '(a)' ) '  the uniform code.  Compare these results with'
  write ( *, '(a)' ) '  the uniform output.'
  write ( *, '(a)' ) ' '

  tdata(1:4) = (/ -1.0D+00, 0.0D+00, 1.0D+00, 2.0D+00 /)

  alpha = ( tdata(3) - tdata(2) ) / ( tdata(3) - tdata(1) )
  beta =  ( tdata(3) - tdata(2) ) / ( tdata(4) - tdata(2) )

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  ALPHA = ', alpha
  write ( *, '(a,g14.6)' ) '  BETA = ', beta

  call basis_matrix_overhauser_nonuni ( alpha, beta, mbasis )

  do i = 1, ndata
    ydata(i) = ( tdata(i) + 2.0D+00 )**2 + 3.0D+00
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '    TDATA, YDATA'
  write ( *, '(a)' ) ' '
  do i = 1, ndata
    write ( *, '(2x,2g14.6)' ) tdata(i), ydata(i)
  end do

  left = 2

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '    T, Spline(T)'
  write ( *, '(a)' ) ' '

  do i = 0, ndata

    if ( i == 0 ) then
      tlo = tdata(1) - 0.5D+00 * ( tdata(2) - tdata(1) )
      thi = tdata(1)
    else if ( i < ndata ) then
      tlo = tdata(i)
      thi = tdata(i+1)
    else if ( ndata <= i ) then
      tlo = tdata(ndata)
      thi = tdata(ndata) + 0.5D+00 * ( tdata(ndata) - tdata(ndata-1) )
    end if

    if ( i < ndata ) then
      jhi = nsample - 1
    else
      jhi = nsample
    end if

    do j = 0, jhi

      tval = ( real ( nsample - j, kind = 8 ) * tlo   &
             + real (           j, kind = 8 ) * thi ) &
             / real ( nsample,     kind = 8 )

      call basis_matrix_tmp ( left, n, mbasis, ndata, tdata, ydata, tval, yval )

      if ( 0 < i .and. j == 0 ) then
        mark = '*'
      else
        mark = ' '
      end if

      write ( *, '(2x,a1,2g14.6)' ) mark, tval, yval

    end do

  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Now test that the nonuniform code on a'
  write ( *, '(a)' ) '  nonuniform grid.'
  write ( *, '(a)' ) ' '

  tdata(1:4) = (/ -4.0D+00, -3.0D+00, -1.0D+00, 2.0D+00 /)

  alpha = ( tdata(3) - tdata(2) ) / ( tdata(3) - tdata(1) )
  beta =  ( tdata(3) - tdata(2) ) / ( tdata(4) - tdata(2) )

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  ALPHA = ', alpha
  write ( *, '(a,g14.6)' ) '  BETA =  ', beta

  call basis_matrix_overhauser_nonuni ( alpha, beta, mbasis )

  ydata(1:ndata) = ( tdata(1:ndata) + 2.0D+00 )**2 + 3.0D+00

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '       TDATA         YDATA'
  write ( *, '(a)' ) ' '
  do i = 1, ndata
    write ( *, '(2x,2g14.6)' ) tdata(i), ydata(i)
  end do

  left = 2

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '        T            Spline(T)'
  write ( *, '(a)' ) ' '

  do i = 0, ndata

    if ( i == 0 ) then
      tlo = tdata(1) - 0.5D+00 * ( tdata(2) - tdata(1) )
      thi = tdata(1)
    else if ( i < ndata ) then
      tlo = tdata(i)
      thi = tdata(i+1)
    else if ( ndata <= i ) then
      tlo = tdata(ndata)
      thi = tdata(ndata) + 0.5D+00 * ( tdata(ndata) - tdata(ndata-1) )
    end if

    if ( i < ndata ) then
      jhi = nsample - 1
    else
      jhi = nsample
    end if

    do j = 0, jhi

      tval = ( real ( nsample - j, kind = 8 ) * tlo   &
             + real (           j, kind = 8 ) * thi ) &
             / real ( nsample,     kind = 8 )

      call basis_matrix_tmp ( left, n, mbasis, ndata, tdata, ydata, tval, yval )

      if ( 0 < i .and. j == 0 ) then
        mark = '*'
      else
        mark = ' '
      end if

      write ( *, '(2x,a1,2g14.6)' ) mark, tval, yval

    end do

  end do

  return
end
subroutine test10

!*****************************************************************************80
!
!! TEST10 tests BC_VAL.
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
  implicit none

  integer ( kind = 4 ), parameter :: n = 2

  integer ( kind = 4 ) i
  integer ( kind = 4 ), parameter :: nsample = 101
  real    ( kind = 8 ) t
  real    (  kind = 8 ), dimension (0:n) :: xcon = (/ 0.0D+00, 0.75D+00, 1.0D+00 /)
  real    ( kind = 8 ) xval
  real    (  kind = 8 ), dimension (0:n) :: ycon = (/ 1.0D+00, 0.0D+00,  1.0D+00 /)
  real    ( kind = 8 ) yval

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST10'
  write ( *, '(a)' ) '  BC_VAL evaluates a general Bezier function.'
!
!  One point on the curve should be about (0.75, 0.536).
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '       T            X(T)          Y(T)'
  write ( *, '(a)' ) ' '

  do i = 1, nsample
    t = real (       i - 1, kind = 8 ) &
      / real ( nsample - 1, kind = 8 )
    call bc_val ( n, t, xcon, ycon, xval, yval )
    write ( *, '(2x,3g14.6)' ) t, xval, yval
  end do
 
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The point ( 0.75, 0.536 ) should be on the curve.'
 
  return
end
subroutine test11

!*****************************************************************************80
!
!! TEST11 tests BEZ_VAL.
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
  implicit none

  integer ( kind = 4 ), parameter :: n = 2

  real    ( kind = 8 ) :: a = 0.0D+00
  real    ( kind = 8 ) :: b = 1.0D+00
  real    ( kind = 8 ) bez_val
  integer ( kind = 4 ) i
  integer ( kind = 4 ) :: nsample = 21
  real    ( kind = 8 ) x
  real    (  kind = 8 ), dimension ( 0 : n ) :: y = (/ 1.0D+00, 0.0D+00, 1.0D+00 /)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST11'
  write ( *, '(a)' ) '  BEZ_VAL evaluates a Bezier function.'
!
!  One point on the curve should be (0.75, 20/32).
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '         T    X(T)          Y(T)'
  write ( *, '(a)' ) ' '

  do i = 1, nsample

    x = ( real ( nsample - i,     kind = 8 ) * a   &
        + real (           i - 1, kind = 8 ) * b ) &
        / real ( nsample     - 1, kind = 8 )

    write ( *, '(2x,i8,2g14.6)' ) i, x, bez_val ( n, x, a, b, y )

  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  When X = ', 0.75D+00
  write ( *, '(a,g14.6)' ) '  BEZ_VAL(X) should be ', 0.625D+00
 
  return
end
subroutine test115

!*****************************************************************************80
!
!! TEST115 tests BP01.
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
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 3

  real    ( kind = 8 ) :: a = 0.0D+00
  real    ( kind = 8 ) :: b = 1.0D+00
  real    ( kind = 8 ) bern(0:n_max)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) n
  integer ( kind = 4 ) :: nsample = 11
  real    ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST115'
  write ( *, '(a)' ) '  BP01 evaluates the Bernstein basis polynomials'
  write ( *, '(a)' ) '  for the interval [0,1].'

  do n = 0, n_max

    write ( *, '(a)' ) ' '
    write ( *, '(a,i8)' ) '  Degree N = ', n
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) &
      '   X         BERN(N,0,X)  BERN(N,1,X)  BERN(N,2,X)  BERN(N,3,X)'
    write ( *, '(a)' ) ' ' 

    do i = 1, nsample

      x = ( real ( nsample - i,     kind = 8 ) * a   &
          + real (           i - 1, kind = 8 ) * b ) &
          / real ( nsample     - 1, kind = 8 )

      call bp01 ( n, x, bern )

      write ( *, '(2x,f8.4,4x,5g14.6)' ) x, bern(0:n)

    end do

  end do

  return
end
subroutine test116

!*****************************************************************************80
!
!! TEST116 tests BPAB.
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
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 3

  real    (  kind = 8 ), parameter :: a = 1.0D+00
  real    (  kind = 8 ), parameter :: b = 3.0D+00
  real    ( kind = 8 ) bern(0:n_max)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) n
  real    ( kind = 8 ) x
  integer ( kind = 4 ), parameter :: nsample = 11

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST116'
  write ( *, '(a)' ) '  BPAB evaluates the Bernstein basis polynomials'
  write ( *, '(a)' ) '  for the interval [A,B].'
  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  A = ', a
  write ( *, '(a,g14.6)' ) '  B = ', b

  do n = 0, n_max

    write ( *, '(a)' ) ' '
    write ( *, '(a,i8)' ) '  Degree N = ', n
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) &
      '   X         BERN(N,0,X)  BERN(N,1,X)  BERN(N,2,X)  BERN(N,3,X)'
    write ( *, '(a)' ) ' ' 

    do i = 1, nsample

      x = ( real ( nsample - i,     kind = 8 ) * a   &
          + real (           i - 1, kind = 8 ) * b ) &
          / real ( nsample     - 1, kind = 8 )

      call bpab ( n, a, b, x, bern )

      write ( *, '(2x,f8.4,4x,5g14.6)' ) x, bern(0:n)

    end do

  end do

  return
end
subroutine test12

!*****************************************************************************80
!
!! TEST12 tests BP_APPROX.
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
  implicit none

  integer ( kind = 4 ), parameter :: maxdata = 10

  real    ( kind = 8 ) a
  real    ( kind = 8 ) b
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ndata
  integer ( kind = 4 ) nsample
  real    ( kind = 8 ) xdata(0:maxdata)
  real    ( kind = 8 ) xval
  real    ( kind = 8 ) ydata(0:maxdata)
  real    ( kind = 8 ) yval

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST12'
  write ( *, '(a)' ) '  BP_APPROX evaluates the Bernstein polynomial'
  write ( *, '(a)' ) '    approximant to a function F(X).'

  a = 1.0D+00
  b = 3.0D+00

  do ndata = 0, 9, 3

    do i = 0, ndata

      if ( ndata == 0 ) then
        xdata(i) = 0.5D+00 * ( a + b )
      else
        xdata(i) = ( real ( ndata - i, kind = 8 ) * a   &
                   + real (         i, kind = 8 ) * b ) &
                   / real ( ndata,     kind = 8 )
      end if

      ydata(i) = sin ( xdata(i) )

    end do

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '       XDATA        YDATA'
    write ( *, '(a)' ) ' '
    do i = 0, ndata
      write ( *, '(2x,2g14.6)' ) xdata(i), ydata(i)
    end do

    write ( *, '(a)' ) ' '
    write ( *, '(a,i8)' ) '  Bernstein approximant of degree N = ', ndata
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '       X            F(X)          BERN(X)        ERROR'
    write ( *, '(a)' ) ' '

    nsample = 2 * ndata + 1

    do i = 1, nsample

      if ( nsample == 1 ) then
        xval = 0.5D+00 * ( a + b )
      else
        xval = ( real ( nsample - i,     kind = 8 ) * a   &
               + real (           i - 1, kind = 8 ) * b ) &
               / real ( nsample     - 1, kind = 8 )
      end if

      call bp_approx ( ndata, a, b, ydata, xval, yval )

      write ( *, '(2x,4g14.6)' ) xval, sin(xval), yval, yval - sin(xval)

    end do

  end do

  return
end
subroutine test125

!*****************************************************************************80
!
!! TEST125 tests LEAST_SET_OLD and LEAST_VAL_OLD.
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
  implicit none

  integer ( kind = 4 ), parameter :: maxdeg = 6
  integer ( kind = 4 ), parameter :: ntab = 21

  real    ( kind = 8 ) b(1:maxdeg)
  real    ( kind = 8 ) c(0:maxdeg)
  real    ( kind = 8 ) d(2:maxdeg)
  real    ( kind = 8 ) eps
  real    ( kind = 8 ) error
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jhi
  integer ( kind = 4 ) ndeg
  real    ( kind = 8 ) ptab(ntab)
  real    ( kind = 8 ) xtab(ntab)
  real    ( kind = 8 ) xval
  real    ( kind = 8 ) ytab(ntab)
  real    ( kind = 8 ) ytrue
  real    ( kind = 8 ) yval

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST125'
  write ( *, '(a)' ) '  LEAST_SET_OLD sets a least squares polynomial,'
  write ( *, '(a)' ) '  LEAST_VAL_OLD evaluates it.'

  do i = 1, ntab
    xtab(i) = ( real ( ntab - i,     kind = 8 ) * ( -1.0D+00 )   &
              + real (        i - 1, kind = 8 ) * ( +1.0D+00 ) ) &
              / real ( ntab     - 1, kind = 8 )
    ytab(i) = real ( int ( exp ( xtab(i) ) * 100.0D+00 + 0.5D+00 ) ) / 100.0D+00
  end do

  write ( *, '(a)'    ) ' '
  write ( *, '(a)'    ) '  The data to be interpolated:'
  write ( *, '(a)'    ) ' '
  write ( *, '(a,i8)' ) '  Number of data values = ', ntab
  write ( *, '(a)'    ) ' '
  write ( *, '(a)'    ) '       X             Y'
  write ( *, '(a)'    ) ' '

  do i = 1, ntab
    write ( *, '(2x,2g14.6)' ) xtab(i), ytab(i)
  end do

  do ndeg = 1, maxdeg

    write ( *, '(a)' ) ' '
    write ( *, '(a,i8)' ) '  Using a polynomial of degree: ', ndeg
    write ( *, '(a)' ) ' '

    call least_set_old ( ntab, xtab, ytab, ndeg, ptab, b, c, d, eps, ierror )

    write ( *, '(a)' ) ' '
    write ( *, '(a,g14.6)' ) '  Total approximation error = ', eps
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '      X            F(X)          P(X)         Error'
    write ( *, '(a)' ) ' '

    do i = 1, ntab

      if ( i < ntab ) then
        jhi = 2
      else
        jhi = 0
      end if

      do j = 0, jhi

        if ( i < ntab ) then

          xval = ( real ( 3 - j, kind = 8 ) * xtab(i)     &
                 + real (     j, kind = 8 ) * xtab(i+1) ) &
                 / real ( 3,     kind = 8 )

        else

          xval = xtab(i)

        end if

        call least_val_old ( xval, ndeg, b, c, d, yval )

        ytrue = real ( int ( exp ( xval ) * 100.0D+00 + 0.5D+00 ), kind = 8 ) &
          / 100.0D+00

        error = yval - ytrue
        write ( *, '(2x,5g14.6)' ) xval, yval, ytrue, error
      end do

    end do

  end do

  return
end
subroutine test126

!*****************************************************************************80
!
!! TEST126 tests LEAST_SET and LEAST_VAL.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    09 December 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: point_num = 21
  integer ( kind = 4 ), parameter :: nterms = 4

  real    ( kind = 8 ) b(nterms)
  real    ( kind = 8 ) c(nterms)
  real    ( kind = 8 ) d(nterms)
  real    ( kind = 8 ) f(point_num)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) nterms2
  real    ( kind = 8 ) px
  real    ( kind = 8 ) w(point_num)
  real    ( kind = 8 ) x(point_num)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST126'
  write ( *, '(a)' ) '  LEAST_SET sets a least squares polynomial,'
  write ( *, '(a)' ) '  LEAST_VAL evaluates it.'

  w(1:point_num) = 1.0D+00

  do i = 1, point_num
    x(i) = - 1.0D+00 + real ( i - 1, kind = 8 ) / 10.0D+00
    f(i) = real ( int ( exp ( x(i) ) * 100.0D+00 + 0.5D+00 ), kind = 8 ) &
      / 100.0D+00
  end do

  call least_set ( point_num, x, f, w, nterms, b, c, d )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  X, F(X), P(X), Error'
  write ( *, '(a)' ) ' '

  do nterms2 = 1, nterms
    write ( *, '(a)' ) ' '
    write ( *, '(a,i8)' ) '  Using polynomial order = ', nterms2
    write ( *, '(a)' ) ' '
    do i = 1, point_num
      call least_val ( nterms2, b, c, d, x(i), px )
      write ( *, '(5g14.6)' ) x(i), f(i), px, px - f(i)
    end do
  end do

  return
end
subroutine test127

!*****************************************************************************80
!
!! TEST127 tests LEAST_SET and LEAST_VAL2.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    09 December 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: point_num = 21
  integer ( kind = 4 ), parameter :: nterms = 4

  real    ( kind = 8 ) b(nterms)
  real    ( kind = 8 ) c(nterms)
  real    ( kind = 8 ) d(nterms)
  real    ( kind = 8 ) f(point_num)
  real    ( kind = 8 ) fp(point_num)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) nterms2
  real    ( kind = 8 ) px
  real    ( kind = 8 ) pxp
  real    ( kind = 8 ) w(point_num)
  real    ( kind = 8 ) x(point_num)

  w(1:point_num) = 1.0D+00

  do i = 1, point_num
    x(i) = -1.0D+00 + real ( i - 1, kind = 8 ) / 10.0D+00
    f(i) = x(i)**2 - x(i) - 6.0D+00
    fp(i) = 2.0D+00 * x(i) - 1.0D+00
  end do

  call least_set ( point_num, x, f, w, nterms, b, c, d )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST127'
  write ( *, '(a)' ) '  LEAST_SET sets a least squares polynomial,'
  write ( *, '(a)' ) '  LEAST_VAL2 evaluates it.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  X, F(X), P(X), FP(X), PP(X)'
  write ( *, '(a)' ) ' '

  do nterms2 = 1, nterms
    write ( *, '(a)' ) ' '
    write ( *, '(a,i8)' ) '  Using polynomial order = ', nterms2
    write ( *, '(a)' ) ' '
    do i = 1, point_num
      call least_val2 ( nterms2, b, c, d, x(i), px, pxp )
      write ( *, '(5g14.6)' ) x(i), f(i), px, fp(i), pxp
    end do
  end do

  return
end
subroutine test13

!*****************************************************************************80
!
!! TEST13 tests SPLINE_B_VAL.
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
  implicit none

  integer ( kind = 4 ), parameter :: ndata = 11

  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jhi
  character mark
  integer ( kind = 4 ), parameter :: nsample = 4
  real    (  kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real    ( kind = 8 ) tdata(ndata)
  real    ( kind = 8 ) thi
  real    ( kind = 8 ) tlo
  real    ( kind = 8 ) tval
  real    ( kind = 8 ) ydata(ndata)
  real    ( kind = 8 ) yval

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST13'
  write ( *, '(a)' ) '  SPLINE_B_VAL evaluates the B spline.'

  do i = 1, ndata
    tdata(i) = real ( i - 1, kind = 8 )
    ydata(i) = sin ( 2.0D+00 * pi * tdata(i) / real ( ndata - 1, kind = 8 ) )
  end do

  write ( *, '(a)'    ) ' '
  write ( *, '(a)'    ) '  The data to be interpolated:'
  write ( *, '(a)'    ) ' '
  write ( *, '(a,i8)' ) '  Number of data values = ', ndata
  write ( *, '(a)'    ) ' '
  write ( *, '(a)'    ) '       T             Y'
  write ( *, '(a)'    ) ' '

  do i = 1, ndata
    write ( *, '(2x,2g14.6)' ) tdata(i), ydata(i)
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '       T           Spline(T)'
  write ( *, '(a)' ) ' '

  do i = 0, ndata

    if ( i == 0 ) then
      tlo = tdata(1) - 0.5D+00 * ( tdata(2) - tdata(1) )
      thi = tdata(1)
    else if ( i < ndata ) then
      tlo = tdata(i)
      thi = tdata(i+1)
    else if ( ndata <= i ) then
      tlo = tdata(ndata)
      thi = tdata(ndata) + 0.5D+00 * ( tdata(ndata) - tdata(ndata-1) )
    end if

    if ( i < ndata ) then
      jhi = nsample - 1
    else
      jhi = nsample
    end if

    do j = 0, jhi

      tval = ( real ( nsample - j, kind = 8 ) * tlo   &
             + real (           j, kind = 8 ) * thi ) &
             / real ( nsample,     kind = 8 )

      call spline_b_val ( ndata, tdata, ydata, tval, yval )

      if ( 0 < i .and. j == 0 ) then
        mark = '*'
      else
        mark = ' '
      end if

      write ( *, '(2x,a1,2g14.6)' ) mark, tval, yval

    end do

  end do

  return
end
subroutine test14

!*****************************************************************************80
!
!! TEST14 tests SPLINE_BETA_VAL.
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
  implicit none

  integer ( kind = 4 ), parameter :: ndata = 11

  real    ( kind = 8 ) beta1
  real    ( kind = 8 ) beta2
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jhi
  character mark
  integer ( kind = 4 ), parameter :: nsample = 4
  real    (  kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real    ( kind = 8 ) tdata(ndata)
  real    ( kind = 8 ) thi
  real    ( kind = 8 ) tlo
  real    ( kind = 8 ) tval
  real    ( kind = 8 ) ydata(ndata)
  real    ( kind = 8 ) yval

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST14'
  write ( *, '(a)' ) '  SPLINE_BETA_VAL evaluates the BETA spline.'

  do i = 1, ndata
    tdata(i) = real ( i - 1, kind = 8 )
    ydata(i) = sin ( 2.0D+00 * pi * tdata(i) / real ( ndata - 1, kind = 8 ) )
  end do

  write ( *, '(a)'    ) ' '
  write ( *, '(a)'    ) '  The data to be interpolated:'
  write ( *, '(a)'    ) ' '
  write ( *, '(a,i8)' ) '  Number of data values = ', ndata
  write ( *, '(a)'    ) ' '
  write ( *, '(a)'    ) '       T             Y'
  write ( *, '(a)'    ) ' '

  do i = 1, ndata
    write ( *, '(2x,2g14.6)' ) tdata(i), ydata(i)
  end do

  beta1 = 1.0D+00
  beta2 = 0.0D+00

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  BETA1 = ', beta1
  write ( *, '(a,g14.6)' ) '  BETA2 = ', beta2
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '    T, Spline(T)'
  write ( *, '(a)' ) ' '

  do i = 0, ndata

    if ( i == 0 ) then
      tlo = tdata(1) - 0.5D+00 * ( tdata(2) - tdata(1) )
      thi = tdata(1)
    else if ( i < ndata ) then
      tlo = tdata(i)
      thi = tdata(i+1)
    else if ( ndata <= i ) then
      tlo = tdata(ndata)
      thi = tdata(ndata) + 0.5D+00 * ( tdata(ndata) - tdata(ndata-1) )
    end if

    if ( i < ndata ) then
      jhi = nsample - 1
    else
      jhi = nsample
    end if

    do j = 0, jhi

      tval = ( real ( nsample - j, kind = 8 ) * tlo   &
             + real (           j, kind = 8 ) * thi ) &
             / real ( nsample,     kind = 8 )

      call spline_beta_val ( beta1, beta2, ndata, tdata, ydata, tval, yval )

      if ( 0 < i .and. j == 0 ) then
        mark = '*'
      else
        mark = ' '
      end if

      write ( *, '(2x,a1,2g14.6)' ) mark, tval, yval

    end do

  end do

  beta1 = 1.0D+00
  beta2 = 100.0D+00

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  BETA1 = ', beta1
  write ( *, '(a,g14.6)' ) '  BETA2 = ', beta2
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '    T, Spline(T)'
  write ( *, '(a)' ) ' '

  do i = 0, ndata

    if ( i == 0 ) then
      tlo = tdata(1) - 0.5D+00 * ( tdata(2) - tdata(1) )
      thi = tdata(1)
    else if ( i < ndata ) then
      tlo = tdata(i)
      thi = tdata(i+1)
    else if ( ndata <= i ) then
      tlo = tdata(ndata)
      thi = tdata(ndata) + 0.5D+00 * ( tdata(ndata) - tdata(ndata-1) )
    end if

    if ( i < ndata ) then
      jhi = nsample - 1
    else
      jhi = nsample
    end if

    do j = 0, jhi

      tval = ( real ( nsample - j, kind = 8 ) * tlo   &
             + real (           j, kind = 8 ) * thi ) &
             / real ( nsample,     kind = 8 )

      call spline_beta_val ( beta1, beta2, ndata, tdata, ydata, tval, yval )

      if ( 0 < i .and. j == 0 ) then
        mark = '*'
      else
        mark = ' '
      end if

      write ( *, '(2x,a1,2g14.6)' ) mark, tval, yval

    end do

  end do

  beta1 = 100.0D+00
  beta2 = 0.0D+00

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  BETA1 = ', beta1
  write ( *, '(a,g14.6)' ) '  BETA2 = ', beta2
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '    T, Spline(T)'
  write ( *, '(a)' ) ' '

  do i = 0, ndata

    if ( i == 0 ) then
      tlo = tdata(1) - 0.5D+00 * ( tdata(2) - tdata(1) )
      thi = tdata(1)
    else if ( i < ndata ) then
      tlo = tdata(i)
      thi = tdata(i+1)
    else if ( ndata <= i ) then
      tlo = tdata(ndata)
      thi = tdata(ndata) + 0.5D+00 * ( tdata(ndata) - tdata(ndata-1) )
    end if

    if ( i < ndata ) then
      jhi = nsample - 1
    else
      jhi = nsample
    end if

    do j = 0, jhi

      tval = ( real ( nsample - j, kind = 8 ) * tlo   &
             + real (           j, kind = 8 ) * thi ) &
             / real ( nsample,     kind = 8 )

      call spline_beta_val ( beta1, beta2, ndata, tdata, ydata, tval, yval )

      if ( 0 < i .and. j == 0 ) then
        mark = '*'
      else
        mark = ' '
      end if

      write ( *, '(2x,a1,2g14.6)' ) mark, tval, yval

    end do

  end do

  return
end
subroutine test143

!*****************************************************************************80
!
!! TEST143 tests SPLINE_BEZIER_VAL.
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
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 1
  integer ( kind = 4 ), parameter :: interval_num = 3

  real    ( kind = 8 ) data_val(dim_num,3*interval_num+1)
  real    ( kind = 8 ) dxdt
  integer ( kind = 4 ) interval
  integer ( kind = 4 ) j
  character mark
  integer ( kind = 4 ) point
  integer ( kind = 4 ) point_num
  real    ( kind = 8 ), allocatable, dimension ( : ) :: point_t
  real    ( kind = 8 ), allocatable, dimension ( :, : ) :: point_val
  real    ( kind = 8 ) t
  real    ( kind = 8 ) t_max
  real    ( kind = 8 ) t_min
  real    ( kind = 8 ) x
  real    ( kind = 8 ), parameter :: x_max = 2.0D+00 * 3.141592653589793D+00
  real    ( kind = 8 ), parameter :: x_min = 0.0D+00
  real    ( kind = 8 ) xdata(0:interval_num)
  real    ( kind = 8 ) ydata(0:interval_num)
  real    ( kind = 8 ) ypdata(0:interval_num)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST143'
  write ( *, '(a)' ) '  SPLINE_BEZIER_VAL evaluates a cubic Bezier spline.'
!
!  Construct the data.
!
  do interval = 0, interval_num

    x = ( real ( interval_num - interval, kind = 8 ) * x_min   &
        + real (                interval, kind = 8 ) * x_max ) &
        / real ( interval_num,            kind = 8 )
    
    xdata(interval) = x
    ydata(interval) =  sin ( x )
    ypdata(interval) = cos ( x )

  end do

  write ( *, '(a)'    ) ' '
  write ( *, '(a)'    ) '  The data to be interpolated:'
  write ( *, '(a)'    ) ' '
  write ( *, '(a,i8)' ) '  Number of intervals = ', interval_num
  write ( *, '(a)'    ) ' '
  write ( *, '(a)'    ) '       X             Y            dYdX'
  write ( *, '(a)'    ) ' '

  do interval = 0, interval_num
    write ( *, '(2x,3g14.6)' ) &
      xdata(interval), ydata(interval), ypdata(interval)
  end do
!
!  Construct the Bezier control data.
!
  dxdt = ( x_max - x_min ) / real ( interval_num, kind = 8 )
  j = 0

  do interval = 1, interval_num

    if ( interval == 1 ) then
      j = j + 1
      data_val(1,j) = ydata(interval-1)
    end if

    j = j + 1
    data_val(1,j) = ydata(interval-1) + ypdata(interval-1) * dxdt / 3.0D+00

    j = j + 1
    data_val(1,j) = ydata(interval)   - ypdata(interval) * dxdt  / 3.0D+00

    j = j + 1
    data_val(1,j) = ydata(interval)

  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The control points'
  write ( *, '(a)' ) '  Interpolation points are marked with a "*".'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '                     T            P(T)'
  write ( *, '(a)' ) ' '

  do j = 1, 3 * interval_num + 1

    t = real ( j - 1, kind = 8 ) / 3.0

    if ( abs ( nint ( t ) - t ) < 0.00001D+00 ) then
      mark = '*'
    else
      mark = ' '
    end if

    write ( *, '(2x,a,2x,i8,2g14.6)' ) mark, j, t, data_val(1,j)
  end do

  point_num = 6 * interval_num + 1
  allocate ( point_t(1:point_num) )
  allocate ( point_val(1:dim_num,1:point_num) )

  t_min = 0.0D+00
  t_max = real ( interval_num, kind = 8 )

  do point = 1, point_num

    t = ( real ( point_num - point,     kind = 8 ) * t_min   &
        + real (             point - 1, kind = 8 ) * t_max ) &
        / real ( point_num         - 1, kind = 8 )
    
    point_t(point) = t

  end do

  call spline_bezier_val ( dim_num, interval_num, data_val, point_num, &
    point_t, point_val )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The Bezier spline, sampled at various points.'
  write ( *, '(a)' ) '  Interpolation points are marked with a "*".'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '                   T        Spline(T)      F(X(T))'
  write ( *, '(a)' ) ' '

  do point = 1, point_num

    if ( abs ( nint ( point_t(point) ) - point_t(point) ) < 0.00001D+00 ) then
      mark = '*'
    else
      mark = ' '
    end if

    x = ( real ( point_num - point,     kind = 8 ) * x_min   &
        + real (             point - 1, kind = 8 ) * x_max ) &
        / real ( point_num         - 1, kind = 8 )

    write ( *, '(2x,a,2x,i8,2x,f10.6,2g14.6)' ) &
      mark, point, point_t(point), point_val(1,point), sin ( x )

  end do

  deallocate ( point_t )
  deallocate ( point_val )

  return
end
subroutine test144

!*****************************************************************************80
!
!! TEST144 tests SPLINE_BEZIER_VAL.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    18 June 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 1
  integer ( kind = 4 ), parameter :: interval_num = 3

  real    ( kind = 8 ) data_val(dim_num,3*interval_num+1)
  integer ( kind = 4 ) interval
  integer ( kind = 4 ) j
  character mark
  integer ( kind = 4 ) point
  integer ( kind = 4 ) point_num
  real    ( kind = 8 ), allocatable, dimension ( : ) :: point_t
  real    ( kind = 8 ), allocatable, dimension ( :, : ) :: point_val
  real    ( kind = 8 ) t
  real    ( kind = 8 ) t_max
  real    ( kind = 8 ) t_min
  real    ( kind = 8 ) x
  real    ( kind = 8 ), parameter :: x_max = 2.0D+00 * 3.141592653589793D+00
  real    ( kind = 8 ), parameter :: x_min = 0.0D+00
  real    ( kind = 8 ) xdata(0:3*interval_num)
  real    ( kind = 8 ) ydata(0:3*interval_num)
  real    ( kind = 8 ) y0
  real    ( kind = 8 ) y1
  real    ( kind = 8 ) y2
  real    ( kind = 8 ) y3
  real    ( kind = 8 ) ypdata(0:interval_num)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST144'
  write ( *, '(a)' ) '  SPLINE_BEZIER_VAL evaluates a cubic Bezier spline.'
  write ( *, '(a)' ) '  Normally, the "interior" points of a Bezier spline'
  write ( *, '(a)' ) '  are not interpolating points.  Instead, the'
  write ( *, '(a)' ) '  derivatives at the interval endpoints are used.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  This example shows, however, that it is possible'
  write ( *, '(a)' ) '  to start with only data values, and to "massage"'
  write ( *, '(a)' ) '  the data so that we can define a cubic Bezier spline'
  write ( *, '(a)' ) '  which interpolates ALL the data.'
!
!  Construct the data.
!
  do interval = 0, 3 * interval_num

    x = ( real ( 3 * interval_num - interval, kind = 8 ) * x_min   &
        + real (                    interval, kind = 8 ) * x_max ) &
        / real ( 3 * interval_num,            kind = 8 )
    
    xdata(interval) = x
    ydata(interval) = sin ( x )

  end do

  write ( *, '(a)'    ) ' '
  write ( *, '(a)'    ) '  The data to be interpolated:'
  write ( *, '(a)'    ) ' '
  write ( *, '(a,i8)' ) '  Number of intervals = ', interval_num
  write ( *, '(a)'    ) ' '
  write ( *, '(a)'    ) '       X             Y'
  write ( *, '(a)'    ) ' '

  do interval = 0, 3 * interval_num
    write ( *, '(2x,3g14.6)' ) xdata(interval), ydata(interval)
  end do
!
!  Construct the Bezier control data.
!  The middle points must be "massaged".
!
  j = 0

  do interval = 1, interval_num

    y0 = ydata(3*(interval-1))
    y1 = ydata(3*(interval-1)+1)
    y2 = ydata(3*(interval-1)+2)
    y3 = ydata(3*(interval-1)+3)

    if ( interval == 1 ) then
      j = j + 1
      data_val(1,j) = y0
    end if

    j = j + 1
    data_val(1,j) = &
      ( -5.0D+00 * y0 + 18.0D+00 * y1 - 9.0D+00 * y2 + 2.0D+00 * y3 ) / 6.0D+00 

    j = j + 1
    data_val(1,j) = &
      ( 2.0D+00 * y0 - 9.0D+00 * y1 + 18.0D+00 * y2 - 5.0D+00 * y3 ) / 6.0D+00

    j = j + 1
    data_val(1,j) = y3

  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The control points'
  write ( *, '(a)' ) '  ALL control points will be interpolation points!'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '                     T            P(T)'
  write ( *, '(a)' ) ' '

  do j = 1, 3 * interval_num + 1

    t = real ( j - 1, kind = 8 ) / 3.0

    mark = '*'

    write ( *, '(2x,a,2x,i8,2g14.6)' ) mark, j, t, data_val(1,j)
  end do

  point_num = 6 * interval_num + 1
  allocate ( point_t(1:point_num) )
  allocate ( point_val(1:dim_num,1:point_num) )

  t_min = 0.0D+00
  t_max = real ( interval_num, kind = 8 )

  do point = 1, point_num

    t = ( real ( point_num - point,     kind = 8 ) * t_min   &
        + real (             point - 1, kind = 8 ) * t_max ) &
        / real ( point_num         - 1, kind = 8 )
    
    point_t(point) = t

  end do

  call spline_bezier_val ( dim_num, interval_num, data_val, point_num, &
    point_t, point_val )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The Bezier spline, sampled at various points.'
  write ( *, '(a)' ) '  Interpolation points are marked with a "*".'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '                   T        Spline(T)      F(X(T))'
  write ( *, '(a)' ) ' '

  do point = 1, point_num

    if ( abs ( nint ( 3 * point_t(point) ) - 3 * point_t(point) ) &
      < 0.00001D+00 ) then
      mark = '*'
    else
      mark = ' '
    end if

    x = ( real ( point_num - point,     kind = 8 ) * x_min   &
        + real (             point - 1, kind = 8 ) * x_max ) &
        / real ( point_num         - 1, kind = 8 )

    write ( *, '(2x,a,2x,i8,2x,f10.6,2g14.6)' ) &
      mark, point, point_t(point), point_val(1,point), sin ( x )

  end do

  deallocate ( point_t )
  deallocate ( point_val )

  return
end
subroutine test145

!*****************************************************************************80
!
!! TEST145 tests SPLINE_CONSTANT_VAL.
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
  implicit none

  integer ( kind = 4 ), parameter :: ndata = 12
  integer ( kind = 4 ), parameter :: n_test = 20

  real    ( kind = 8 ) ahi
  real    ( kind = 8 ) alo
  real    ( kind = 8 ) frunge
  real    ( kind = 8 ) fval
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) :: seed = 123456789
  real    ( kind = 8 ) tdata(ndata-1)
  real    ( kind = 8 ) thi
  real    ( kind = 8 ) tlo
  real    ( kind = 8 ) t_test(n_test)
  real    ( kind = 8 ) tval
  real    ( kind = 8 ) ydata(ndata)
  real    ( kind = 8 ) yval

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST145'
  write ( *, '(a)' ) '  SPLINE_CONSTANT_VAL evaluates a piecewise '
  write ( *, '(a)' ) '  constant spline.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Runge''s function, evenly spaced knots.'
!
!  Set the data.
!
  tlo = -1.0D+00
  thi = +1.0D+00
  call r8vec_even ( ndata-1, tlo, thi, tdata )

  do i = 1, ndata

    if ( i == 1 ) then
      tval = tdata(1)
    else if ( i < ndata ) then
      tval = 0.5D+00 * ( tdata(i-1) + tdata(i) )
    else if ( i == ndata ) then
      tval = tdata(i-1)
    end if

    ydata(i) = frunge ( tval )

  end do

  write ( *, '(a)'    ) ' '
  write ( *, '(a)'    ) '  The data to be interpolated:'
  write ( *, '(a)'    ) ' '
  write ( *, '(a,i8)' ) '  Number of data values = ', ndata
  write ( *, '(a)'    ) ' '
  write ( *, '(a)'    ) '       T             Y'
  write ( *, '(a)'    ) ' '

  do i = 1, ndata
    write ( *, '(2x,a1,14x,g14.6)' ) '*', ydata(i)
    if ( i < ndata ) then
      write ( *, '(2x,a1, g14.6)' ) '*', tdata(i)
    end if
  end do
!
!  Sample the spline.
!
  write ( *, * ) 'DEBUG: TLO = ', tlo
  write ( *, * ) 'DEBUG: THI = ', thi

  alo = tlo - 1.0D+00
  ahi = thi + 1.0D+00

  call r8vec_uniform_01 ( n_test, seed, t_test )

  t_test(1:n_test) = ( 1.0D+00 - t_test(1:n_test) ) * alo &
                   +             t_test(1:n_test)   * ahi

  call r8vec_sort_bubble_a ( n_test, t_test )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '     T     Y(interp)    Y(exact)'
  write ( *, '(a)' ) ' '

  j = 0
  write ( *, '(2x,a1,14x,g14.6)' ) '*', ydata(j+1)
  j = j + 1

  do i = 1, n_test

    tval = t_test(i)

    call spline_constant_val ( ndata, tdata, ydata, tval, yval )

    if ( j <= ndata - 1 ) then
      do while ( tdata(j) <= tval )
        fval = frunge ( tdata(j) )
        write ( *, '(2x,a1,g14.6,14x,g14.6)' ) '*', tdata(j), fval
        write ( *, '(2x,a1,14x,g14.6)' ) '*', ydata(j+1)
        j = j + 1
        if ( ndata <= j ) then
          exit
        end if
      end do
    end if

    fval = frunge ( tval )

    write ( *, '(2x,a1,3g14.6)' ) ' ', tval, yval, fval

  end do

  return
end
subroutine test15

!*****************************************************************************80
!
!! TEST15 tests SPLINE_CUBIC_SET and SPLINE_CUBIC_VAL.
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
  implicit none

  integer ( kind = 4 ), parameter :: n = 11

  real    ( kind = 8 ) frunge
  real    ( kind = 8 ) fprunge
  real    ( kind = 8 ) fpprunge
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ibcbeg
  integer ( kind = 4 ) ibcend
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jhi
  integer ( kind = 4 ) k
  real    ( kind = 8 ) t(n)
  real    ( kind = 8 ) tval
  real    ( kind = 8 ) y(n)
  real    ( kind = 8 ) ybcbeg
  real    ( kind = 8 ) ybcend
  real    ( kind = 8 ) ypp(n)
  real    ( kind = 8 ) yppval
  real    ( kind = 8 ) ypval
  real    ( kind = 8 ) yval

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST15'
  write ( *, '(a)' ) '  SPLINE_CUBIC_SET sets up a cubic spline;'
  write ( *, '(a)' ) '  SPLINE_CUBIC_VAL evaluates it.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Runge''s function, evenly spaced knots.'

  do i = 1, n

    t(i) = ( real ( n - i,     kind = 8 ) * (-1.0D+00)   &
           + real (     i - 1, kind = 8 ) * (+1.0D+00) ) &
           / real ( n     - 1, kind = 8 )

    y(i) =  frunge ( t(i) )

  end do

  write ( *, '(a)'    ) ' '
  write ( *, '(a)'    ) '  The data to be interpolated:'
  write ( *, '(a)'    ) ' '
  write ( *, '(a,i8)' ) '  Number of data values = ', n
  write ( *, '(a)'    ) ' '
  write ( *, '(a)'    ) '       T             Y'
  write ( *, '(a)'    ) ' '

  do i = 1, n
    write ( *, '(2x,g14.6,2x,g14.6)' ) t(i), y(i)
  end do
!
!  Try boundary condition types 0, 1 and 2.
!
  do k = 0, 3

    if ( k == 0 ) then

      ibcbeg = 0
      ybcbeg = 0.0D+00

      ibcend = 0
      ybcend = 0.0D+00

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Boundary condition 0 at both ends:'
      write ( *, '(a)' ) '  Spline is quadratic in boundary intervals.'

    else if ( k == 1 ) then

      ibcbeg = 1
      ybcbeg = fprunge ( t(1) )

      ibcend = 1
      ybcend = fprunge ( t(n) )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Boundary condition 1 at both ends:'
      write ( *, '(a,g14.6)' ) '  Y''(left) =  ', ybcbeg
      write ( *, '(a,g14.6)' ) '  Y''(right) = ', ybcend

    else if ( k == 2 ) then

      ibcbeg = 2
      ybcbeg = fpprunge ( t(1) )

      ibcend = 2
      ybcend = fpprunge ( t(n) )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Boundary condition 2 at both ends:'
      write ( *, '(a,g14.6)' ) '  YP"(left) =  ', ybcbeg
      write ( *, '(a,g14.6)' ) '  YP"(right) = ', ybcend

    else if ( k == 3 ) then

      ibcbeg = 2
      ybcbeg = 0.0D+00

      ibcend = 2
      ybcend = 0.0D+00

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  "Natural" spline:'
      write ( *, '(a)' ) '  Boundary condition 2 at both ends:'
      write ( *, '(a,g14.6)' ) '  YP"(left) =  ', ybcbeg
      write ( *, '(a,g14.6)' ) '  YP"(right) = ', ybcend

    end if

    call spline_cubic_set ( n, t, y, ibcbeg, ybcbeg, ibcend, ybcend, ypp )

    if ( k == 3 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '       I      Y(I)          YPP(I)'
      write ( *, '(a)' ) ' '
      do i = 1, n
        write ( *, '(2x,i8,2g14.6)' ) i, y(i), ypp(i)
      end do
    end if

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '    SPLINE"(T)      F"(T):'
    write ( *, '(a)' ) ' '
    do i = 1, n
      write ( *, '(2x,2g14.6)' ) ypp(i), fpprunge(t(i))
    end do

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '      T            SPLINE(T)      F(T)'
    write ( *, '(a)' ) ' '

    do i = 0, n

      if ( i == 0 ) then
        jhi = 1
      else if ( i < n ) then
        jhi = 2
      else
        jhi = 2
      end if

      do j = 1, jhi

        if ( i == 0 ) then
          tval = t(1) - 1.0D+00
        else if ( i < n ) then
          tval = ( real ( jhi - j + 1, kind = 8 ) * t(i)     &
                 + real (       j - 1, kind = 8 ) * t(i+1) ) &
                 / real ( jhi,         kind = 8 )
        else
          if ( j == 1 ) then
            tval = t(n)
          else
            tval = t(n) + 1.0D+00
          end if
        end if

        call spline_cubic_val ( n, t, y, ypp, tval, yval, ypval, yppval )

        write ( *, '(2x,3g14.6)' ) tval, yval, frunge ( tval )

      end do
    end do

  end do

  return
end
subroutine test16

!*****************************************************************************80
!
!! TEST16 tests SPLINE_CUBIC_SET and SPLINE_CUBIC_VAL2.
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
  implicit none

  integer ( kind = 4 ), parameter :: n = 11

  real    ( kind = 8 ) frunge
  real    ( kind = 8 ) fprunge
  real    ( kind = 8 ) fpprunge
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ibcbeg
  integer ( kind = 4 ) ibcend
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jhi
  integer ( kind = 4 ) k
  integer ( kind = 4 ) left
  integer ( kind = 4 ) left_in
  real    ( kind = 8 ) t(n)
  real    ( kind = 8 ) tval
  real    ( kind = 8 ) y(n)
  real    ( kind = 8 ) ybcbeg
  real    ( kind = 8 ) ybcend
  real    ( kind = 8 ) ypp(n)
  real    ( kind = 8 ) yppval
  real    ( kind = 8 ) ypval
  real    ( kind = 8 ) yval

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST16'
  write ( *, '(a)' ) '  SPLINE_CUBIC_SET sets up a cubic spline;'
  write ( *, '(a)' ) '  SPLINE_CUBIC_VAL2 evaluates it.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Runge''s function, evenly spaced knots.'
  do i = 1, n

    t(i) = ( real ( n - i,     kind = 8 ) * (-1.0D+00)   &
           + real (     i - 1, kind = 8 ) * (+1.0D+00) ) &
           / real ( n     - 1, kind = 8 )

    y(i) =  frunge ( t(i) )

  end do

  write ( *, '(a)'    ) ' '
  write ( *, '(a)'    ) '  The data to be interpolated:'
  write ( *, '(a)'    ) ' '
  write ( *, '(a,i8)' ) '  Number of data values = ', n
  write ( *, '(a)'    ) ' '
  write ( *, '(a)'    ) '       T             Y'
  write ( *, '(a)'    ) ' '

  do i = 1, n
    write ( *, '(2x,2g14.6)' ) t(i), y(i)
  end do
!
!  Try boundary condition types 0, 1 and 2.
!
  do k = 0, 2

    if ( k == 0 ) then

      ibcbeg = 0
      ybcbeg = 0.0D+00

      ibcend = 0
      ybcend = 0.0D+00

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Boundary condition 0 at both ends:'
      write ( *, '(a)' ) '  Spline is quadratic in boundary intervals.'

    else if ( k == 1 ) then

      ibcbeg = 1
      ybcbeg = fprunge ( t(1) )

      ibcend = 1
      ybcend = fprunge ( t(n) )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Boundary condition 1 at both ends:'
      write ( *, '(a,g14.6)' ) '  Y''(left) =  ', ybcbeg
      write ( *, '(a,g14.6)' ) '  Y''(right) = ', ybcend

    else if ( k == 2 ) then

      ibcbeg = 2
      ybcbeg = fpprunge ( t(1) )

      ibcend = 2
      ybcend = fpprunge ( t(n) )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Boundary condition 2 at both ends:'
      write ( *, '(a,g14.6)' ) '  YP"(left) =  ', ybcbeg
      write ( *, '(a,g14.6)' ) '  YP"(right) = ', ybcend

    end if

    call spline_cubic_set ( n, t, y, ibcbeg, ybcbeg, ibcend, ybcend, ypp )

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  SPLINE"(T)      F"(T)'
    write ( *, '(a)' ) ' '
    do i = 1, n
      write ( *, '(2x,2g14.6)' ) ypp(i), fpprunge(t(i))
    end do

    left = 0

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) &
      '      T             SPLINE(T)       F(T)       LEFT_IN  LEFT_OUT'
    write ( *, '(a)' ) ' '

    do i = 0, n

      if ( i == 0 ) then
        jhi = 1
      else if ( i < n ) then
        jhi = 2
      else
        jhi = 2
      end if

      do j = 1, jhi

        if ( i == 0 ) then
          tval = t(1) - 1.0D+00
        else if ( i < n ) then
          tval = ( real ( jhi - j + 1, kind = 8 ) * t(i)     &
                 + real (       j - 1, kind = 8 ) * t(i+1) ) &
                 / real ( jhi,         kind = 8 )
        else
          if ( j == 1 ) then
            tval = t(n)
          else
            tval = t(n) + 1.0D+00
          end if
        end if

        left_in = left

        call spline_cubic_val2 ( n, t, y, ypp, left, tval, yval, ypval, yppval )

        write ( *, '(2x,3g14.6,2i8)' ) tval, yval, frunge ( tval ), &
          left_in, left

      end do
    end do

  end do

  return
end
subroutine test17

!*****************************************************************************80
!
!! TEST17 tests SPLINE_CUBIC_SET and SPLINE_CUBIC_VAL.
!
!  Discussion:
!
!    For boundary condition 0, the spline should come very close within
!    the interpolation interval.
!
!    For conditions 1 and 2, the spline should be essentially exactly equal
!    to the data, inside and outside the interpolation interval.
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
  implicit none

  integer ( kind = 4 ), parameter :: n = 11

  real    ( kind = 8 ) fcube
  real    ( kind = 8 ) fpcube
  real    ( kind = 8 ) fppcube
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ibcbeg
  integer ( kind = 4 ) ibcend
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jhi
  integer ( kind = 4 ) k
  real    ( kind = 8 ) t(n)
  real    ( kind = 8 ) tval
  real    ( kind = 8 ) y(n)
  real    ( kind = 8 ) ybcbeg
  real    ( kind = 8 ) ybcend
  real    ( kind = 8 ) ypp(n)
  real    ( kind = 8 ) yppval
  real    ( kind = 8 ) ypval
  real    ( kind = 8 ) yval

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST17'
  write ( *, '(a)' ) '  SPLINE_CUBIC_SET sets up a cubic spline;'
  write ( *, '(a)' ) '  SPLINE_CUBIC_VAL evaluates it.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Cubic data, unevenly spaced knots.'

  do i = 1, n
    t(i) = ( real ( i - 1, kind = 8 ) &
           / real ( n - 1, kind = 8 ) )**2
  end do

  do i = 1, n
    y(i) = fcube ( t(i) )
  end do

  write ( *, '(a)'    ) ' '
  write ( *, '(a)'    ) '  The data to be interpolated:'
  write ( *, '(a)'    ) ' '
  write ( *, '(a,i8)' ) '  Number of data values = ', n
  write ( *, '(a)'    ) ' '
  write ( *, '(a)'    ) '       T             Y'
  write ( *, '(a)'    ) ' '

  do i = 1, n
    write ( *, '(2x,g14.6,2x,g14.6)' ) t(i), y(i)
  end do
!
!  Try boundary condition types 0, 1 and 2.
!
  do k = 0, 2

    if ( k == 0 ) then

      ibcbeg = 0
      ybcbeg = 0.0D+00

      ibcend = 0
      ybcend = 0.0D+00

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Boundary condition 0 at both ends:'
      write ( *, '(a)' ) '  Spline is quadratic in boundary intervals.'

    else if ( k == 1 ) then

      ibcbeg = 1
      ybcbeg = fpcube ( t(1) )

      ibcend = 1
      ybcend = fpcube ( t(n) )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Boundary condition 1 at both ends:'
      write ( *, '(a,g14.6)' ) '  Y''(left) =  ', ybcbeg
      write ( *, '(a,g14.6)' ) '  Y''(right) = ', ybcend

    else if ( k == 2 ) then

      ibcbeg = 2
      ybcbeg = fppcube ( t(1) )

      ibcend = 2
      ybcend = fppcube ( t(n) )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Boundary condition 2 at both ends:'
      write ( *, '(a,g14.6)' ) '  YP"(left) =  ', ybcbeg
      write ( *, '(a,g14.6)' ) '  YP"(right) = ', ybcend

    end if

    call spline_cubic_set ( n, t, y, ibcbeg, ybcbeg, ibcend, ybcend, ypp )

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '       SPLINE"(T)    F"(T):'
    write ( *, '(a)' ) ' '
    do i = 1, n
      write ( *, '(2x,2g14.6)' ) ypp(i), fppcube(t(i))
    end do

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '       T           SPLINE(T)      F(T)'
    write ( *, '(a)' ) ' '

    do i = 0, n

      if ( i == 0 ) then
        jhi = 1
      else if ( i < n ) then
        jhi = 2
      else
        jhi = 2
      end if

      do j = 1, jhi

        if ( i == 0 ) then
          tval = t(1) - 1.0D+00
        else if ( i < n ) then
          tval = ( real ( jhi - j + 1, kind = 8 ) * t(i)     &
                 + real (       j - 1, kind = 8 ) * t(i+1) ) &
                 / real ( jhi,         kind = 8 )
        else
          if ( j == 1 ) then
            tval = t(n)
          else
            tval = t(n) + 1.0D+00
          end if
        end if

        call spline_cubic_val ( n, t, y, ypp, tval, yval, ypval, yppval )

        write ( *, '(2x,3g14.6)' ) tval, yval, fcube ( tval )

      end do
    end do

  end do

  return
end
subroutine test18

!*****************************************************************************80
!
!! TEST18 tests SPLINE_CUBIC_SET and SPLINE_CUBIC_VAL.
!
!  Discussion:
!
!    For boundary condition 0, the spline should come very close within
!    the interpolation interval.
!
!    For conditions 1 and 2, the spline should be essentially exactly equal
!    to the data, inside and outside the interpolation interval.
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
  implicit none

  integer ( kind = 4 ), parameter :: n = 11

  real    ( kind = 8 ) fcube
  real    ( kind = 8 ) fpcube
  real    ( kind = 8 ) fppcube
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ibcbeg
  integer ( kind = 4 ) ibcend
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jhi
  integer ( kind = 4 ) k
  real    ( kind = 8 ) t(n)
  real    ( kind = 8 ) tval
  real    ( kind = 8 ) y(n)
  real    ( kind = 8 ) ybcbeg
  real    ( kind = 8 ) ybcend
  real    ( kind = 8 ) ypp(n)
  real    ( kind = 8 ) yppval
  real    ( kind = 8 ) ypval
  real    ( kind = 8 ) yval

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST18'
  write ( *, '(a)' ) '  SPLINE_CUBIC_SET sets up a cubic spline;'
  write ( *, '(a)' ) '  SPLINE_CUBIC_VAL evaluates it.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Cubic data, evenly spaced knots.'
  do i = 1, n
    t(i) = real ( i - 1, kind = 8 ) / real ( n - 1, kind = 8 )
    y(i) =  fcube ( t(i) )
  end do

  write ( *, '(a)'    ) ' '
  write ( *, '(a)'    ) '  The data to be interpolated:'
  write ( *, '(a)'    ) ' '
  write ( *, '(a,i8)' ) '  Number of data values = ', n
  write ( *, '(a)'    ) ' '
  write ( *, '(a)'    ) '       T             Y'
  write ( *, '(a)'    ) ' '

  do i = 1, n
    write ( *, '(2x,g14.6,2x,g14.6)' ) t(i), y(i)
  end do
!
!  Try boundary condition types 0, 1 and 2.
!
  do k = 0, 2

    if ( k == 0 ) then

      ibcbeg = 0
      ybcbeg = 0.0D+00

      ibcend = 0
      ybcend = 0.0D+00

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Boundary condition 0 at both ends:'
      write ( *, '(a)' ) '  Spline is quadratic in boundary intervals.'

    else if ( k == 1 ) then

      ibcbeg = 1
      ybcbeg = fpcube ( t(1) )

      ibcend = 1
      ybcend = fpcube ( t(n) )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Boundary condition 1 at both ends:'
      write ( *, '(a,g14.6)' ) '  Y''(left) =  ', ybcbeg
      write ( *, '(a,g14.6)' ) '  Y''(right) = ', ybcend

    else if ( k == 2 ) then

      ibcbeg = 2
      ybcbeg = fppcube ( t(1) )

      ibcend = 2
      ybcend = fppcube ( t(n) )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Boundary condition 2 at both ends:'
      write ( *, '(a,g14.6)' ) '  YP"(left) =  ', ybcbeg
      write ( *, '(a,g14.6)' ) '  YP"(right) = ', ybcend

    end if

    call spline_cubic_set ( n, t, y, ibcbeg, ybcbeg, ibcend, ybcend, ypp )

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '   SPLINE"(T)      F"(T):'
    write ( *, '(a)' ) ' '
    do i = 1, n
      write ( *, '(2x,2g14.6)' ) ypp(i), fppcube(t(i))
    end do

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '        T      SPLINE(T)    F(T)'
    write ( *, '(a)' ) ' '

    do i = 0, n

      if ( i == 0 ) then
        jhi = 1
      else if ( i < n ) then
        jhi = 2
      else
        jhi = 2
      end if

      do j = 1, jhi

        if ( i == 0 ) then
          tval = t(1) - 1.0D+00
        else if ( i < n ) then
          tval = ( real ( jhi - j + 1, kind = 8 ) * t(i)     &
                 + real (       j - 1, kind = 8 ) * t(i+1) ) &
                 / real ( jhi,         kind = 8 )
        else
          if ( j == 1 ) then
            tval = t(n)
          else
            tval = t(n) + 1.0D+00
          end if
        end if

        call spline_cubic_val ( n, t, y, ypp, tval, yval, ypval, yppval )

        write ( *, '(2x,f10.4)' ) tval
        write ( *, '(2x,10x,2f10.4)' )   yval, fcube ( tval )
        write ( *, '(2x,10x,2f10.4)' )   ypval, fpcube ( tval )
        write ( *, '(2x,10x,2f10.4)' )   yppval, fppcube ( tval )

      end do
    end do

  end do

  return
end
subroutine test19

!*****************************************************************************80
!
!! TEST19 tests SPLINE_CUBIC_SET and SPLINE_CUBIC_VAL.
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
  implicit none

  integer ( kind = 4 ), parameter :: n = 2

  real    ( kind = 8 ) fcube
  real    ( kind = 8 ) fpcube
  real    ( kind = 8 ) fppcube
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ibcbeg
  integer ( kind = 4 ) ibcend
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jhi
  integer ( kind = 4 ) k1
  integer ( kind = 4 ) k2
  real    ( kind = 8 ) t(n)
  real    ( kind = 8 ) tval
  real    ( kind = 8 ) y(n)
  real    ( kind = 8 ) ybcbeg
  real    ( kind = 8 ) ybcend
  real    ( kind = 8 ) ypp(n)
  real    ( kind = 8 ) yppval
  real    ( kind = 8 ) ypval
  real    ( kind = 8 ) yval

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST19'
  write ( *, '(a)' ) '  SPLINE_CUBIC_SET sets up a cubic spline;'
  write ( *, '(a)' ) '  SPLINE_CUBIC_VAL evaluates it.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Cubic data, evenly spaced knots.'
  write ( *, '(a)' ) '  ONLY TWO KNOTS!'

  do i = 1, n
    t(i) = real ( i - 1, kind = 8 ) / real ( n - 1, kind = 8 )
    y(i) = fcube ( t(i) )
  end do

  write ( *, '(a)'    ) ' '
  write ( *, '(a)'    ) '  The data to be interpolated:'
  write ( *, '(a)'    ) ' '
  write ( *, '(a,i8)' ) '  Number of data values = ', n
  write ( *, '(a)'    ) ' '
  write ( *, '(a)'    ) '       T             Y'
  write ( *, '(a)'    ) ' '

  do i = 1, n
    write ( *, '(2x,g14.6,2x,g14.6)' ) t(i), y(i)
  end do
!
!  Try all 9 pairs of boundary condition types 0, 1 and 2.
!
  do k1 = 0, 2

    if ( k1 == 0 ) then

      ibcbeg = 0
      ybcbeg = 0.0D+00

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Boundary condition 0 at left end.'

    else if ( k1 == 1 ) then

      ibcbeg = 1
      ybcbeg = fpcube ( t(1) )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Boundary condition 1 at left end.'
      write ( *, '(a,g14.6)' ) '  Y''(left) =  ', ybcbeg

    else if ( k1 == 2 ) then

      ibcbeg = 2
      ybcbeg = fppcube ( t(1) )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Boundary condition 2 at left ends:'
      write ( *, '(a,g14.6)' ) '  YP"(left) =  ', ybcbeg

    end if

    do k2 = 0, 2

      if ( k2 == 0 ) then

        ibcend = 0
        ybcend = 0.0D+00

        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) '  Boundary condition 0 at right end.'

      else if ( k2 == 1 ) then

        ibcend = 1
        ybcend = fpcube ( t(n) )

        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) '  Boundary condition 1 at right end.'
        write ( *, '(a,g14.6)' ) '  Y''(right) = ', ybcend

      else if ( k2 == 2 ) then

        ibcend = 2
        ybcend = fppcube ( t(n) )

        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) '  Boundary condition 2 at right end.'
        write ( *, '(a,g14.6)' ) '  YP"(right) = ', ybcend

      end if

      call spline_cubic_set ( n, t, y, ibcbeg, ybcbeg, ibcend, ybcend, ypp )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  SPLINE"(T), F"(T):'
      write ( *, '(a)' ) ' '
      do i = 1, n
        write ( *, '(2x,2g14.6)' ) ypp(i), fppcube(t(i))
      end do

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  T, SPLINE(T), F(T)'
      write ( *, '(a)' ) ' '

      do i = 0, n

        if ( i == 0 ) then
          jhi = 1
        else if ( i < n ) then
          jhi = 2
        else
          jhi = 2
        end if

        do j = 1, jhi

          if ( i == 0 ) then
            tval = t(1) - 1.0D+00
          else if ( i < n ) then
            tval = ( real ( jhi - j + 1, kind = 8 ) * t(i)     &
                   + real (       j - 1, kind = 8 ) * t(i+1) ) &
                   / real ( jhi,         kind = 8 )
          else
            if ( j == 1 ) then
              tval = t(n)
            else
              tval = t(n) + 1.0D+00
            end if
          end if

          call spline_cubic_val ( n, t, y, ypp, tval, yval, ypval, yppval )

          write ( *, '(2x,f10.4)' ) tval
          write ( *, '(2x,10x,2f10.4)' )   yval, fcube ( tval )
          write ( *, '(2x,10x,2f10.4)' )   ypval, fpcube ( tval )
          write ( *, '(2x,10x,2f10.4)' )   yppval, fppcube ( tval )

        end do
      end do

    end do

  end do

  return
end
subroutine test20

!*****************************************************************************80
!
!! TEST20 tests SPLINE_HERMITE_SET and SPLINE_HERMITE_VAL.
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
  implicit none

  integer ( kind = 4 ), parameter :: ndata = 4

  real    ( kind = 8 ) c(4,ndata)
  real    ( kind = 8 ) fpval
  real    ( kind = 8 ) fval
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jhi
  character mark
  real    (  kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real    ( kind = 8 ) tdata(ndata)
  real    ( kind = 8 ) tval
  real    ( kind = 8 ) ydata(ndata)
  real    ( kind = 8 ) ypdata(ndata)
  real    ( kind = 8 ) ypval
  real    ( kind = 8 ) yval

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST20'
  write ( *, '(a)' ) '  SPLINE_HERMITE_SET sets up a Hermite spline;'
  write ( *, '(a)' ) '  SPLINE_HERMITE_VAL evaluates it.'
!
!  Set the data.
!
  do i = 1, ndata
    tdata(i) = ( real ( ndata - i,     kind = 8 ) *   0.0D+00          &
               + real (         i - 1, kind = 8 ) * ( 0.5D+00 * pi ) ) &
               / real ( ndata     - 1, kind = 8 )
    ydata(i) = sin ( tdata(i) )
    ypdata(i) = cos ( tdata(i) )
  end do

  write ( *, '(a)'    ) ' '
  write ( *, '(a)'    ) '  The data to be interpolated:'
  write ( *, '(a)'    ) ' '
  write ( *, '(a,i8)' ) '  Number of data values = ', ndata
  write ( *, '(a)'    ) ' '
  write ( *, '(a)'    ) '       T             Y             Y'''
  write ( *, '(a)'    ) ' '
  do i = 1, ndata
    write ( *, '(2x,3g14.6)' ) tdata(i), ydata(i), ypdata(i)
  end do
!
!  Set up the spline.
!
  call spline_hermite_set ( ndata, tdata, ydata, ypdata, c )
!
!  Now evaluate the spline all over the place.
!
  write ( *, '(a)'    ) ' '
  write ( *, '(a)'    ) '       T        Y(hermite)     ' // &
    'Y(exact)      Y''(hermite)   Y''(exact)'
  write ( *, '(a)'    ) ' '

  do i = 1, ndata

    if ( i == ndata ) then
      jhi = 0
    else
      jhi = 2
    end if

    do j = 0, jhi

      tval = real ( 3 * ( i - 1 ) + j, kind = 8 ) * ( 0.5D+00 * pi ) &
           / real ( 3 * ( ndata - 1 ), kind = 8 )

      fval = sin ( tval )
      fpval = cos ( tval )
 
      call spline_hermite_val ( ndata, tdata, c, tval, yval, ypval )

      if (j == 0 ) then
        mark = '*'
      else
        mark = ' '
      end if

      write ( *, '(2x,a1,5g14.6)' ) mark, tval, yval, fval, ypval, fpval

    end do

  end do

  return
end
subroutine test205

!*****************************************************************************80
!
!! TEST205 tests SPLINE_LINEAR_INT and SPLINE_LINEAR_INTSET.
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
  implicit none

  integer ( kind = 4 ), parameter :: n = 4

  real    ( kind = 8 ) a
  real    ( kind = 8 ) b
  real    ( kind = 8 ) data_x(n)
  real    ( kind = 8 ) data_y(n)
  integer ( kind = 4 ) i
  real    ( kind = 8 ) int_x(n+1)
  real    ( kind = 8 ) int_v(n)
  real    ( kind = 8 ) value

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST205'
  write ( *, '(a)' ) '  SPLINE_LINEAR_INTSET is given some interval endpoints,'
  write ( *, '(a)' ) '  and a value associated with each interval.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  It determines a linear spline, with breakpoints'
  write ( *, '(a)' ) '  at the centers of each interval, whose integral'
  write ( *, '(a)' ) '  over each interval is equal to the given value.'

  int_x(1:n+1) = (/ 0.0D+00, 1.0D+00, 4.0D+00, 5.0D+00, 10.0D+00 /)
  int_v(1:n) = (/ 10.0D+00, 2.0D+00, 8.0D+00, 27.5D+00 /)

  call r8vec_print ( n+1, int_x, '  The interval end points:' )

  call r8vec_print ( n, int_v, '  The desired interval integral values:' )

  call spline_linear_intset ( n, int_x, int_v, data_x, data_y )

  call r8vec_print ( n, data_x, '  The spline break points:' )
  call r8vec_print ( n, data_y, '  The spline data values: ' )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  As a check, call SPLINE_LINEAR_INT to compute'
  write ( *, '(a)' ) '  the integral of the spline over each interval,'
  write ( *, '(a)' ) '  and compare to the desired value.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '     A       B    Desired      Computed'
  write ( *, '(a)' ) ' '

  do i = 1, n
    a = int_x(i)
    b = int_x(i+1)
    call spline_linear_int ( n, data_x, data_y, a, b, value )
    write ( *, '(2x,2f8.2,2g14.6)' ) a, b, int_v(i), value
  end do

  return
end
subroutine test21

!*****************************************************************************80
!
!! TEST21 tests SPLINE_LINEAR_VAL.
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
  implicit none

  integer ( kind = 4 ), parameter :: n = 11

  real    ( kind = 8 ) frunge
  real    ( kind = 8 ) fval
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jhi
  real    ( kind = 8 ) t(n)
  real    ( kind = 8 ) tval
  real    ( kind = 8 ) y(n)
  real    ( kind = 8 ) ypval
  real    ( kind = 8 ) yval

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST21'
  write ( *, '(a)' ) '  SPLINE_LINEAR_VAL evaluates a linear spline.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Runge''s function, evenly spaced knots.'

  do i = 1, n

    t(i) = ( real ( n - i,     kind = 8 ) * (-1.0D+00)   &
           + real (     i - 1, kind = 8 ) * (+1.0D+00) ) &
           / real ( n     - 1, kind = 8 )

    y(i) =  frunge ( t(i) )

  end do

  write ( *, '(a)'    ) ' '
  write ( *, '(a)'    ) '  The data to be interpolated:'
  write ( *, '(a)'    ) ' '
  write ( *, '(a,i8)' ) '  Number of data values = ', n
  write ( *, '(a)'    ) ' '
  write ( *, '(a)'    ) '       T             Y'
  write ( *, '(a)'    ) ' '
  do i = 1, n
    write ( *, '(2x,2g14.6)' ) t(i), y(i)
  end do

  write ( *, '(a)'    ) ' '
  write ( *, '(a)'    ) '  Interpolation:'
  write ( *, '(a)'    ) ' '
  write ( *, '(a)'    ) '       T             Y            Yexact'
  write ( *, '(a)'    ) ' '

  do i = 0, n

    if ( i == 0 ) then
      jhi = 1
    else if ( i < n ) then
      jhi = 2
    else
      jhi = 2
    end if

    do j = 1, jhi

      if ( i == 0 ) then
        tval = t(1) - 1.0D+00
      else if ( i < n ) then
        tval = ( real ( jhi - j + 1, kind = 8 ) * t(i)     &
               + real (       j - 1, kind = 8 ) * t(i+1) ) &
               / real ( jhi,         kind = 8 )
      else
        if ( j == 1 ) then
          tval = t(n)
        else
          tval = t(n) + 1.0D+00
        end if
      end if

      call spline_linear_val ( n, t, y, tval, yval, ypval )

      fval = frunge ( tval )

      write ( *, '(2x,3g14.6)' ) tval, yval, fval

    end do

  end do

  return
end
subroutine test215

!*****************************************************************************80
!
!! TEST215 tests SPLINE_LINEAR_INT.
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
  implicit none

  integer ( kind = 4 ), parameter :: n = 3

  real    ( kind = 8 ) a
  real    ( kind = 8 ) b
  integer ( kind = 4 ) i
  real    ( kind = 8 ) int_val
  real    (  kind = 8 ), dimension ( n ) :: t = (/ 2.0D+00, 4.5D+00, 7.5D+00 /)
  real    (  kind = 8 ), dimension ( n ) :: y = (/ 3.0D+00, 3.75D+00, 5.5D+00 /)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST215'
  write ( *, '(a)' ) '  SPLINE_LINEAR_INT computes the integral '
  write ( *, '(a)' ) '  of a linear spline.'
  write ( *, '(a)' ) ' '

  write ( *, '(a)'    ) ' '
  write ( *, '(a)'    ) '  The data to be interpolated:'
  write ( *, '(a)'    ) ' '
  write ( *, '(a,i8)' ) '  Number of data values = ', n
  write ( *, '(a)'    ) ' '
  write ( *, '(a)'    ) '       T             Y'
  write ( *, '(a)'    ) ' '
  do i = 1, n
    write ( *, '(2x,2g14.6)' ) t(i), y(i)
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '    A             B           Integral'
  write ( *, '(a)' ) ' '

  do i = 1, 5

    if ( i == 1 ) then
      a = 0.0D+00
      b = 4.0D+00
    else if ( i == 2 ) then
      a = 4.0D+00
      b = 5.0D+00
    else if ( i == 3 ) then
      a = 5.0D+00
      b = 10.0D+00
    else if ( i == 4 ) then
      a = 0.0D+00
      b = 10.0D+00
    else
      a = 10.0D+00
      b = 0.0D+00
    end if

    call spline_linear_int ( n, t, y, a, b, int_val )

    write ( *, '(2x,3g14.6)' ) a, b, int_val

  end do

  return
end
subroutine test22

!*****************************************************************************80
!
!! TEST22 tests SPLINE_OVERHAUSER_UNI_VAL.
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
  implicit none

  integer ( kind = 4 ), parameter :: ndata = 11

  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jhi
  character mark
  integer ( kind = 4 ), parameter :: nsample = 4
  real    (  kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real    ( kind = 8 ) tdata(ndata)
  real    ( kind = 8 ) thi
  real    ( kind = 8 ) tlo
  real    ( kind = 8 ) tval
  real    ( kind = 8 ) ydata(ndata)
  real    ( kind = 8 ) yval

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST22'
  write ( *, '(a)' ) '  SPLINE_OVERHAUSER_UNI_VAL evaluates the'
  write ( *, '(a)' ) '    uniform Overhauser spline.'

  do i = 1, ndata
    tdata(i) = real ( i - 1, kind = 8 )
    ydata(i) = sin ( 2.0D+00 * pi * tdata(i) / real ( ndata - 1, kind = 8 ) )
  end do

  write ( *, '(a)'    ) ' '
  write ( *, '(a)'    ) '  The data to be interpolated:'
  write ( *, '(a)'    ) ' '
  write ( *, '(a,i8)' ) '  Number of data values = ', ndata
  write ( *, '(a)'    ) ' '
  write ( *, '(a)'    ) '         T             Y'
  write ( *, '(a)'    ) ' '
  do i = 1, ndata
    write ( *, '(2x,2g14.6)' ) tdata(i), ydata(i)
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '      T, Spline(T)'
  write ( *, '(a)' ) ' '

  do i = 0, ndata

    if ( i == 0 ) then
      tlo = tdata(1) - 0.5D+00 * ( tdata(2) - tdata(1) )
      thi = tdata(1)
    else if ( i < ndata ) then
      tlo = tdata(i)
      thi = tdata(i+1)
    else if ( ndata <= i ) then
      tlo = tdata(ndata)
      thi = tdata(ndata) + 0.5D+00 * ( tdata(ndata) - tdata(ndata-1) )
    end if

    if ( i < ndata ) then
      jhi = nsample - 1
    else
      jhi = nsample
    end if

    do j = 0, jhi

      tval = ( real ( nsample - j, kind = 8 ) * tlo   &
             + real (           j, kind = 8 ) * thi ) &
             / real ( nsample,     kind = 8 )

      call spline_overhauser_uni_val ( ndata, tdata, ydata, tval, yval )

      if ( 0 < i .and. j == 0 ) then
        mark = '*'
      else
        mark = ' '
      end if

      write ( *, '(2x,a1,2g14.6)' ) mark, tval, yval

    end do

  end do

  return
end
subroutine test225

!*****************************************************************************80
!
!! TEST225 tests SPLINE_OVERHAUSER_NONUNI_VAL.
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
  implicit none

  integer ( kind = 4 ), parameter :: ndata = 11

  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jhi
  character mark
  integer ( kind = 4 ), parameter :: nsample = 4
  real    (  kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real    ( kind = 8 ) tdata(ndata)
  real    ( kind = 8 ) thi
  real    ( kind = 8 ) tlo
  real    ( kind = 8 ) tval
  real    ( kind = 8 ) ydata(ndata)
  real    ( kind = 8 ) yval

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST225'
  write ( *, '(a)' ) '  SPLINE_OVERHAUSER_NONUNI_VAL evaluates the'
  write ( *, '(a)' ) '    nonuniform Overhauser spline.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  In this initial draft of a test, we simply'
  write ( *, '(a)' ) '  use uniform nodes.'

  do i = 1, ndata
    tdata(i) = real ( i - 1, kind = 8 )
    ydata(i) = sin ( 2.0D+00 * pi * tdata(i) / real ( ndata - 1, kind = 8 ) )
  end do

  write ( *, '(a)'    ) ' '
  write ( *, '(a)'    ) '  The data to be interpolated:'
  write ( *, '(a)'    ) ' '
  write ( *, '(a,i8)' ) '  Number of data values = ', ndata
  write ( *, '(a)'    ) ' '
  write ( *, '(a)'    ) '         T             Y'
  write ( *, '(a)'    ) ' '
  do i = 1, ndata
    write ( *, '(2x,2g14.6)' ) tdata(i), ydata(i)
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '      T, Spline(T)'
  write ( *, '(a)' ) ' '

  do i = 0, ndata

    if ( i == 0 ) then
      tlo = tdata(1) - 0.5D+00 * ( tdata(2) - tdata(1) )
      thi = tdata(1)
    else if ( i < ndata ) then
      tlo = tdata(i)
      thi = tdata(i+1)
    else if ( ndata <= i ) then
      tlo = tdata(ndata)
      thi = tdata(ndata) + 0.5D+00 * ( tdata(ndata) - tdata(ndata-1) )
    end if

    if ( i < ndata ) then
      jhi = nsample - 1
    else
      jhi = nsample
    end if

    do j = 0, jhi

      tval = ( real ( nsample - j, kind = 8 ) * tlo   &
             + real (           j, kind = 8 ) * thi ) &
             / real ( nsample,     kind = 8 )

      call spline_overhauser_nonuni_val ( ndata, tdata, ydata, tval, yval )

      if ( 0 < i .and. j == 0 ) then
        mark = '*'
      else
        mark = ' '
      end if

      write ( *, '(2x,a1,2g14.6)' ) mark, tval, yval

    end do

  end do

  return
end
subroutine test23

!*****************************************************************************80
!
!! TEST23 tests SPLINE_OVERHAUSER_VAL.
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
  implicit none

  integer ( kind = 4 ), parameter :: ndata = 4
  integer ( kind = 4 ), parameter :: dim_num = 2

  integer ( kind = 4 ) i
  real    ( kind = 8 ) tdata(ndata)
  real    ( kind = 8 ) tval
  real    ( kind = 8 ) ydata(dim_num,ndata)
  real    ( kind = 8 ) yval(dim_num)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST23'
  write ( *, '(a)' ) '  SPLINE_OVERHAUSER_VAL evaluates the'
  write ( *, '(a)' ) '    Overhauser spline.'
!
!  Set the data.
!
  tdata(1) = 1.0D+00
  ydata(1,1) =   0.0D+00
  ydata(2,1) =   0.0D+00

  tdata(2) = 2.0D+00
  ydata(1,2) =   1.0D+00
  ydata(2,2) =   1.0D+00

  tdata(3) = 3.0D+00
  ydata(1,3) =   2.0D+00
  ydata(2,3) = - 1.0D+00

  tdata(4) = 4.0D+00
  ydata(1,4) =   3.0D+00
  ydata(2,4) =   0.0D+00

  write ( *, '(a)'    ) ' '
  write ( *, '(a)'    ) '  The data to be interpolated:'
  write ( *, '(a)'    ) ' '
  write ( *, '(a,i8)' ) '  Number of data values = ', ndata
  write ( *, '(a)'    ) ' '
  write ( *, '(a)'    ) '       T             Y'
  write ( *, '(a)'    ) ' '

  do i = 1, ndata
    write ( *, '(2x,3g14.6)' ) tdata(i), ydata(1:dim_num,i)
  end do
!
!  Now evaluate the spline all over the place.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  T, Spline value'
  write ( *, '(a)' ) ' '

  do i = 0, 6 * ndata + 3

    tval = real ( i, kind = 8 ) / 6.0D+00
    call spline_overhauser_val ( dim_num, ndata, tdata, ydata, tval, yval )
    write ( *, '(2x,3g14.6)' ) tval, yval(1:dim_num)

  end do

  return
end
subroutine test235

!*****************************************************************************80
!
!! TEST235 tests SPLINE_PCHIP_SET and SPLINE_PCHIP_VAL.
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
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 21
  integer ( kind = 4 ), parameter :: ne = 101

  real    ( kind = 8 ) d(n)
  real    ( kind = 8 ) diff
  real    ( kind = 8 ) f(n)
  real    ( kind = 8 ) fd(ne)
  real    ( kind = 8 ) fe(ne)
  integer ( kind = 4 ) i
  real    ( kind = 8 ), external :: frunge
  real    ( kind = 8 ) x(n)
  real    ( kind = 8 ) xe(ne)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST235'
  write ( *, '(a)' ) '  SPLINE_PCHIP_SET carries out piecewise cubic '
  write ( *, '(a)' ) '    Hermite interpolation.'
  write ( *, '(a)' ) '  SPLINE_PCHIP_VAL evaluates the interpolant.'
  write ( *, '(a)' ) ' '
!
!  Compute Runge's function at N points in [-1,1].
!
  do i = 1, n
    x(i) = -1.0D+00 + real ( i - 1, kind = 8 ) / 10.0D+00
    f(i) = frunge ( x(i) )
  end do
!
!  SPLINE_PCHIP_SET takes the data in X and F, and constructs a table in D
!  that defines the interpolant.
!
  call spline_pchip_set ( n, x, f, d )
!
!  Evaluate the interpolant and derivative at NE points from -1 to 0.
!
  do i = 1, ne
    xe(i) = -1.0D+00 + real ( i - 1, kind = 8 ) / real ( ne - 1, kind = 8 )
  end do

  call spline_pchip_val ( n, x, f, d, ne, xe, fe )
!
!  Print the table of X, F(exact) and F(interpolated)
!
  do i = 1, ne
    diff = fe(i) - frunge ( xe(i) )
    write ( *, '(2x,f8.4,2x,f10.6,2x,f10.6,2x,g14.6)' ) &
      xe(i), frunge ( xe(i) ), fe(i), diff
  end do

  return
end
subroutine test24

!*****************************************************************************80
!
!! TEST24 tests SPLINE_QUADRATIC_VAL.
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
  implicit none

  integer ( kind = 4 ), parameter :: n = 11

  real    ( kind = 8 ) frunge
  real    ( kind = 8 ) fval
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jhi
  real    ( kind = 8 ) t(n)
  real    ( kind = 8 ) tval
  real    ( kind = 8 ) y(n)
  real    ( kind = 8 ) ypval
  real    ( kind = 8 ) yval

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST24'
  write ( *, '(a)' ) '  SPLINE_QUADRATIC_VAL evaluates a '
  write ( *, '(a)' ) '    quadratic spline.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Runge''s function, evenly spaced knots.'

  do i = 1, n

    t(i) = ( real ( n - i,     kind = 8 ) * (-1.0D+00)   &
           + real (     i - 1, kind = 8 ) * (+1.0D+00) ) &
           / real ( n     - 1, kind = 8 )

    y(i) =  frunge ( t(i) )

  end do

  write ( *, '(a)'    ) ' '
  write ( *, '(a)'    ) '  The data to be interpolated:'
  write ( *, '(a)'    ) ' '
  write ( *, '(a,i8)' ) '  Number of data values = ', n
  write ( *, '(a)'    ) ' '
  write ( *, '(a)'    ) '         T             Y'
  write ( *, '(a)'    ) ' '
  do i = 1, n
    write ( *, '(2x,2g14.6)' ) t(i), y(i)
  end do

  write ( *, '(a)'    ) ' '
  write ( *, '(a)'    ) '  Interpolated values'
  write ( *, '(a)'    ) ' '
  write ( *, '(a)'    ) '       T             Y           Y(exact)'
  write ( *, '(a)'    ) ' '

  do i = 0, n

    if ( i == 0 ) then
      jhi = 1
    else if ( i < n ) then
      jhi = 2
    else
      jhi = 2
    end if

    do j = 1, jhi

      if ( i == 0 ) then
        tval = t(1) - 1.0D+00
      else if ( i < n ) then
        tval = ( real ( jhi - j + 1, kind = 8 ) * t(i)     &
               + real (       j - 1, kind = 8 ) * t(i+1) ) &
               / real ( jhi,         kind = 8 )
      else
        if ( j == 1 ) then
          tval = t(n)
        else
          tval = t(n) + 1.0D+00
        end if
      end if

      call spline_quadratic_val ( n, t, y, tval, yval, ypval )

      fval = frunge ( tval )

      write ( *, '(2x,3g14.6)' ) tval, yval, fval

    end do

  end do

  return
end
function fcube ( x )

!*****************************************************************************80
!
!! FCUBE evaluates a cubic function.
!
!  Discussion:
!
!    Y(X) = ( ( X + 2 ) * X + 3 ) * X + 4
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
!    Input, real X, the argument.
!
!    Output, real FCUBE, the value of the function.
!
  implicit none

  real    ( kind = 8 ) fcube
  real    ( kind = 8 ) x

  fcube = ( ( (       1.0D+00 ) &
                * x + 2.0D+00 ) &
                * x + 3.0D+00 ) &
                * x + 4.0D+00

  return
end
function fpcube ( x )

!*****************************************************************************80
!
!! FPCUBE evaluates the derivative of a cubic function.
!
!  Discussion:
!
!    Y(X) = ( ( X + 2 ) * X + 3 ) * X + 4
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
!    Input, real X, the argument.
!
!    Output, real FPCUBE, the value of the derivative of the cubic function.
!
  implicit none

  real    ( kind = 8 ) fpcube
  real    ( kind = 8 ) x

  fpcube = ( 3.0D+00 * x + 4.0D+00 ) * x + 3.0D+00

  return
end
function fppcube ( x )

!*****************************************************************************80
!
!! FPPCUBE evaluates the second derivative of a cubic function.
!
!  Discussion:
!
!    Y(X) = ( ( X + 2 ) * X + 3 ) * X + 4
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
!    Input, real X, the argument.
!
!    Output, real FPPCUBE, the second derivative of the cubic function.
!
  implicit none

  real    ( kind = 8 ) fppcube
  real    ( kind = 8 ) x

  fppcube = 6.0D+00 * x + 4.0D+00

  return
end
function frunge ( x )

!*****************************************************************************80
!
!! FRUNGE evaluates the Runge function.
!
!  Discussion:
!
!    Interpolation of the Runge function at evenly spaced points in [-1,1]
!    is a common test.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    24 January 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real X, the argument.
!
!    Output, real FRUNGE, the value of the function.
!
  implicit none

  real    ( kind = 8 ) frunge
  real    ( kind = 8 ) x

  frunge = 1.0D+00 / ( 1.0D+00 + 25.0D+00 * x * x )

  return
end
function fprunge ( x )

!*****************************************************************************80
!
!! FPRUNGE evaluates the derivative of the Runge function.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    24 January 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real X, the argument.
!
!    Output, real FPRUNGE, the value of the derivative of the Runge function.
!
  implicit none

  real    ( kind = 8 ) fprunge
  real    ( kind = 8 ) x

  fprunge = - 50.0D+00 * x / ( 1.0D+00 + 25.0D+00 * x * x )**2

  return
end
function fpprunge ( x )

!*****************************************************************************80
!
!! FPPRUNGE evaluates the second derivative of the Runge function.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    24 January 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real X, the argument.
!
!    Output, real FPPRUNGE, the value of the second derivative of
!    the Runge function.
!
  implicit none

  real    ( kind = 8 ) fpprunge
  real    ( kind = 8 ) x

  fpprunge = ( - 50.0D+00 + 3750.0D+00 * x * x ) &
    / ( 1.0D+00 + 25.0D+00 * x * x )**3

  return
end
