program main

!*****************************************************************************80
!
!! MAIN is the main program for INTERP_PRB.
!
!  Discussion:
!
!    INTERP_PRB calls the INTERP tests.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    03 December 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) data_num

  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'INTERP_PRB'
  write ( *, '(a)' ) '  FORTRAN90 version:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Test the routines in the INTERP library.'
 
  call test01

  data_num = 6
  call test02 ( data_num )

  data_num = 11
  call test02 ( data_num )

  data_num = 6
  call test03 ( data_num )

  data_num = 11
  call test03 ( data_num )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'INTERP_PRB'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine test01

!*****************************************************************************80
!
!! TEST01 tests INTERP_LINEAR on 1-dimensional data.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    03 December 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer( kind = 4 ), parameter :: data_num = 11
  integer( kind = 4 ), parameter :: dim_num = 1

  integer ( kind = 4 ) after
  integer ( kind = 4 ) before
  integer ( kind = 4 ) fat
  integer ( kind = 4 ) i
  integer ( kind = 4 ) interp
  integer ( kind = 4 ) interp_num
  real    ( kind = 8 ) p
  real    ( kind = 8 ) p_data(dim_num,data_num)
  real    ( kind = 8 ), allocatable, dimension ( :, : ) :: p_interp
  real    ( kind = 8 ), allocatable, dimension ( : ) :: p_value
  real    ( kind = 8 ) t
  real    ( kind = 8 ) t_data(data_num)
  real    ( kind = 8 ), allocatable, dimension ( : ) :: t_interp
  real    ( kind = 8 ) t_max
  real    ( kind = 8 ) t_min

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST01'
  write ( *, '(a)' ) '  INTERP_LINEAR evaluates a piecewise linear spline.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  In this example, the function we are interpolating is'
  write ( *, '(a)' ) '  Runge''s function, with evenly spaced knots.'

  t_min = -1.0D+00
  t_max = +1.0D+00

  call ncc_abscissas_ab ( t_min, t_max, data_num, t_data )

  call f_runge ( data_num, t_data, p_data )

  write ( *, '(a)'    ) ' '
  write ( *, '(a)'    ) '  The data to be interpolated:'
  write ( *, '(a)'    ) ' '
  write ( *, '(a,i8)' ) '  Data dimension =        ', dim_num
  write ( *, '(a,i8)' ) '  Number of data values = ', data_num
  write ( *, '(a)'    ) ' '
  write ( *, '(a)'    ) '       T_data        P_data'
  write ( *, '(a)'    ) ' '
  do i = 1, data_num
    write ( *, '(2x,2g14.6)' ) t_data(i), p_data(1,i)
  end do
!
!  Our interpolation values will include the original T values, plus
!  3 new values in between each pair of original values.
!
  before = 4
  fat = 3
  after = 2

  interp_num = before + 1 + ( data_num - 1 ) * ( fat + 1 ) + after

  allocate ( t_interp(interp_num) )
  allocate ( p_interp(dim_num,interp_num) )

  call r8vec_expand_linear2 ( data_num, t_data, before, fat, after, t_interp )

  call interp_linear ( dim_num, data_num, t_data, p_data, interp_num, &
    t_interp, p_interp )

  allocate ( p_value(1:interp_num) )

  call f_runge ( interp_num, t_interp, p_value )

  write ( *, '(a)'    ) ' '
  write ( *, '(a)'    ) '  Interpolation:'
  write ( *, '(a)'    ) ' '
  write ( *, '(a)'    ) '    T_interp      P_interp        P_exact        Error'
  write ( *, '(a)'    ) ' '

  do interp = 1, interp_num

    write ( *, '(2x,f10.4,2x,g14.6,2x,g14.6,2x,g10.2)' ) &
      t_interp(interp), p_interp(1,interp), p_value(interp), &
      p_interp(1,interp) - p_value(interp)

  end do

  deallocate ( p_interp )
  deallocate ( p_value )
  deallocate ( t_interp )

  return
end
subroutine test02 ( data_num )

!*****************************************************************************80
!
!! TEST02 tests INTERP_LAGRANGE on 1-dimensional data, equally spaced data.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    03 December 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer( kind = 4 ) data_num
  integer( kind = 4 ), parameter :: dim_num = 1

  integer ( kind = 4 ) after
  integer ( kind = 4 ) before
  integer ( kind = 4 ) fat
  integer ( kind = 4 ) i
  integer ( kind = 4 ) interp
  integer ( kind = 4 ) interp_num
  real    ( kind = 8 ) p
  real    ( kind = 8 ) p_data(dim_num,data_num)
  real    ( kind = 8 ), allocatable, dimension ( :, : ) :: p_interp
  real    ( kind = 8 ), allocatable, dimension ( : ) :: p_value
  real    ( kind = 8 ) t
  real    ( kind = 8 ) t_data(data_num)
  real    ( kind = 8 ), allocatable, dimension ( : ) :: t_interp
  real    ( kind = 8 ) t_max
  real    ( kind = 8 ) t_min

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST02'
  write ( *, '(a)' ) '  INTERP_LAGRANGE evaluates a polynomial interpolant.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  In this example, the function we are interpolating is'
  write ( *, '(a)' ) '  Runge''s function, with evenly spaced knots.'

  t_min = -1.0D+00
  t_max = +1.0D+00

  call ncc_abscissas_ab ( t_min, t_max, data_num, t_data )

  call f_runge ( data_num, t_data, p_data )

  write ( *, '(a)'    ) ' '
  write ( *, '(a)'    ) '  The data to be interpolated:'
  write ( *, '(a)'    ) ' '
  write ( *, '(a,i8)' ) '  Data dimension =        ', dim_num
  write ( *, '(a,i8)' ) '  Number of data values = ', data_num
  write ( *, '(a)'    ) ' '
  write ( *, '(a)'    ) '       T_data        P_data'
  write ( *, '(a)'    ) ' '
  do i = 1, data_num
    write ( *, '(2x,2g14.6)' ) t_data(i), p_data(1,i)
  end do
!
!  Our interpolation values will include the original T values, plus
!  3 new values in between each pair of original values.
!
  before = 4
  fat = 3
  after = 2

  interp_num = before + 1 + ( data_num - 1 ) * ( fat + 1 ) + after

  allocate ( t_interp(interp_num) )
  allocate ( p_interp(dim_num,interp_num) )

  call r8vec_expand_linear2 ( data_num, t_data, before, fat, after, t_interp )

  call interp_lagrange ( dim_num, data_num, t_data, p_data, interp_num, &
    t_interp, p_interp )

  allocate ( p_value(1:interp_num) )

  call f_runge ( interp_num, t_interp, p_value )

  write ( *, '(a)'    ) ' '
  write ( *, '(a)'    ) '  Interpolation:'
  write ( *, '(a)'    ) ' '
  write ( *, '(a)'    ) '    T_interp      P_interp        P_exact        Error'
  write ( *, '(a)'    ) ' '

  do interp = 1, interp_num

    write ( *, '(2x,f10.4,2x,g14.6,2x,g14.6,2x,g10.2)' ) &
      t_interp(interp), p_interp(1,interp), p_value(interp), &
      p_interp(1,interp) - p_value(interp)

  end do

  deallocate ( p_interp )
  deallocate ( p_value )
  deallocate ( t_interp )

  return
end
subroutine test03 ( data_num )

!*****************************************************************************80
!
!! TEST03 tests INTERP_LAGRANGE on 1-dimensional data, Clenshaw-Curtis data.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    04 December 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer( kind = 4 ) data_num
  integer( kind = 4 ), parameter :: dim_num = 1

  integer ( kind = 4 ) after
  integer ( kind = 4 ) before
  integer ( kind = 4 ) fat
  integer ( kind = 4 ) i
  integer ( kind = 4 ) interp
  integer ( kind = 4 ) interp_num
  real    ( kind = 8 ) p
  real    ( kind = 8 ) p_data(dim_num,data_num)
  real    ( kind = 8 ), allocatable, dimension ( :, : ) :: p_interp
  real    ( kind = 8 ), allocatable, dimension ( : ) :: p_value
  real    ( kind = 8 ) t
  real    ( kind = 8 ) t_data(data_num)
  real    ( kind = 8 ), allocatable, dimension ( : ) :: t_interp
  real    ( kind = 8 ) t_max
  real    ( kind = 8 ) t_min

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST03'
  write ( *, '(a)' ) '  INTERP_LAGRANGE evaluates a polynomial interpolant.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  In this example, the function we are interpolating is'
  write ( *, '(a)' ) '  Runge''s function, with Clenshaw Curtis knots.'

  t_min = -1.0D+00
  t_max = +1.0D+00

  call cc_abscissas_ab ( t_min, t_max, data_num, t_data )

  call f_runge ( data_num, t_data, p_data )

  write ( *, '(a)'    ) ' '
  write ( *, '(a)'    ) '  The data to be interpolated:'
  write ( *, '(a)'    ) ' '
  write ( *, '(a,i8)' ) '  Data dimension =        ', dim_num
  write ( *, '(a,i8)' ) '  Number of data values = ', data_num
  write ( *, '(a)'    ) ' '
  write ( *, '(a)'    ) '       T_data        P_data'
  write ( *, '(a)'    ) ' '
  do i = 1, data_num
    write ( *, '(2x,2g14.6)' ) t_data(i), p_data(1,i)
  end do
!
!  Our interpolation values will include the original T values, plus
!  3 new values in between each pair of original values.
!
  before = 4
  fat = 3
  after = 2

  interp_num = before + 1 + ( data_num - 1 ) * ( fat + 1 ) + after

  allocate ( t_interp(interp_num) )
  allocate ( p_interp(dim_num,interp_num) )

  call r8vec_expand_linear2 ( data_num, t_data, before, fat, after, t_interp )

  call interp_lagrange ( dim_num, data_num, t_data, p_data, interp_num, &
    t_interp, p_interp )

  allocate ( p_value(1:interp_num) )

  call f_runge ( interp_num, t_interp, p_value )

  write ( *, '(a)'    ) ' '
  write ( *, '(a)'    ) '  Interpolation:'
  write ( *, '(a)'    ) ' '
  write ( *, '(a)'    ) '    T_interp      P_interp        P_exact        Error'
  write ( *, '(a)'    ) ' '

  do interp = 1, interp_num

    write ( *, '(2x,f10.4,2x,g14.6,2x,g14.6,2x,g10.2)' ) &
      t_interp(interp), p_interp(1,interp), p_value(interp), &
      p_interp(1,interp) - p_value(interp)

  end do

  deallocate ( p_interp )
  deallocate ( p_value )
  deallocate ( t_interp )

  return
end
subroutine test04 ( data_num )

!*****************************************************************************80
!
!! TEST04 tests INTERP_LAGRANGE on 2-dimensional, equally spaced T-data.
!
!  Discussion:
!
!    The independent variable T has dimension 2.
!    The   dependent variable P has dimension 1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    20 January 2008
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer( kind = 4 ) data_num
  integer( kind = 4 ), parameter :: dim_num = 2

  integer ( kind = 4 ) after
  integer ( kind = 4 ) before
  integer ( kind = 4 ) fat
  integer ( kind = 4 ) i
  integer ( kind = 4 ) interp
  integer ( kind = 4 ) interp_num
  real    ( kind = 8 ) p
  real    ( kind = 8 ), allocatable, dimension ( :    ) :: p_data
  real    ( kind = 8 ), allocatable, dimension ( :, : ) :: p_interp
  real    ( kind = 8 ), allocatable, dimension ( :    ) :: p_value
  real    ( kind = 8 ) t
  real    ( kind = 8 ), allocatable, dimension ( :, : ) :: t_data
  real    ( kind = 8 ), allocatable, dimension ( :, : ) :: t_interp
  real    ( kind = 8 ) t_max(dim_num)
  real    ( kind = 8 ) t_min(dim_num)

  allocate ( p_data(1:data_num) )
  allocate ( t_data(1:dim_num,1:data_num) )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST04'
  write ( *, '(a)' ) '  INTERP_LAGRANGE evaluates a polynomial interpolant.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  In this example, the function we are interpolating is'
  write ( *, '(a)' ) '  Runge''s function, with evenly spaced knots.'

  t_min(1:dim_num) = -1.0D+00
  t_max(1:dim_num) = +1.0D+00

  do dim = 1, dim_num
    call ncc_abscissas_ab ( t_min(dim), t_max(dim), data_num, &
    t_data(dim,1:data_num) )
  end do

  call f_runge_2d ( data_num, t_data, p_data )

  write ( *, '(a)'    ) ' '
  write ( *, '(a)'    ) '  The data to be interpolated:'
  write ( *, '(a)'    ) ' '
  write ( *, '(a,i8)' ) '  Data dimension =        ', dim_num
  write ( *, '(a,i8)' ) '  Number of data values = ', data_num
  write ( *, '(a)'    ) ' '
  write ( *, '(a)'    ) '         T_data                 P_data'
  write ( *, '(a)'    ) ' '
  write ( *, '(a)'    ) '      X           Y             F(X,Y)'
  write ( *, '(a)'    ) '  ----------  ----------    --------------'
  write ( *, '(a)'    ) ' '
  do i = 1, data_num
    write ( *, '(2x,f10.4,2x,f10.4,4x,g14.6)' ) t_data(1:2,i), p_data(1,i)
  end do
!
!  I HOPE THIS IS THE RIGHT WAY TO EXTEND THIS TO 2D...
!
!  Our interpolation values will include the original T values, plus
!  3 new values in between each pair of original values.
!
  before = 4
  fat = 3
  after = 2

  interp_num = before + 1 + ( data_num - 1 ) * ( fat + 1 ) + after

  allocate ( t_interp(1:dim_num,1:interp_num) )
  allocate ( p_interp(1:interp_num) )

  do dim = 1, dim_num
    call r8vec_expand_linear2 ( data_num, t_data(dim,1:data_num), before, 
    fat, after, t_interp(dim,1:interp_num) )
  end do
!
!  DO YOU REALLY THINK THIS WILL WORK?  NOPE!
!
  call interp_lagrange ( dim_num, data_num, t_data, p_data, interp_num, &
    t_interp, p_interp )

  allocate ( p_value(1:interp_num) )

  call f_runge_2d ( interp_num, t_interp, p_value )

  write ( *, '(a)'    ) ' '
  write ( *, '(a)'    ) '  Interpolation:'
  write ( *, '(a)'    ) ' '
  write ( *, '(a)'    ) '        T_interp              P_interp        P_exact        Error'
  write ( *, '(a)'    ) '      X           Y '
  write ( *, '(a)'    ) '  ----------  ----------'
  write ( *, '(a)'    ) ' '

  do interp = 1, interp_num

    write ( *, '(2x,f10.4,2x,f10.4,2x,g14.6,2x,g14.6,2x,g10.2)' ) &
      t_interp(1:2interp), p_interp(1,interp), p_value(interp), &
      p_interp(1,interp) - p_value(interp)

  end do

  deallocate ( p_data )
  deallocate ( p_interp )
  deallocate ( p_value )
  deallocate ( t_data )
  deallocate ( t_interp )

  return
end
subroutine f_runge ( n, x, f )

!*****************************************************************************80
!
!! F_RUNGE evaluates the Runge function.
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
!    03 December 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of evaluation points.
!
!    Input, real X(N), the evaluation points.
!
!    Output, real F(N), the function values.
!
  implicit none

  integer ( kind = 4 ) n

  real    ( kind = 8 ) f(n)
  real    ( kind = 8 ) x(n)
  
  f(1:n) = 1.0D+00 / ( 1.0D+00 + 25.0D+00 * x(1:n)**2 )

  return
end
subroutine f_runge_2d ( n, x, f )

!*****************************************************************************80
!
!! F_RUNGE_2D evaluates the Runge function in 2D.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    20 January 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of evaluation points.
!
!    Input, real X(2,N), the evaluation points.
!
!    Output, real F(N), the function values.
!
  implicit none

  integer ( kind = 4 ) n

  real    ( kind = 8 ) f(n)
  real    ( kind = 8 ) x(2,n)
  
  f(1:n) = 1.0D+00 / ( 1.0D+00 + 25.0D+00 * ( x(1,1:n)**2 + x(1,1:n)**2 ) )

  return
end
