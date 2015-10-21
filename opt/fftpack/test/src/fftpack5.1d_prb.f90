program main

!*****************************************************************************80
!
!! MAIN is the main program for FFTPACK5.1D_PRB.
!
!  Discussion:
!
!    FFTPACK5.1D_PRB tests the FFTPACK5.1D library.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 November 2011
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'FFTPACK5.1D_PRB'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test the FFTPACK5.1D library.'

  call test01 ( )
  call test02 ( )
  call test03 ( )
  call test04 ( )
  call test05 ( )
  call test06 ( )
  call test07 ( )
  call test09 ( )

  call test10 ( )
  call test11 ( )
  call test12 ( )
  call test13 ( )
  call test14 ( )
  call test15 ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'FFTPACK5.1D_PRB'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine test01 ( )

!*****************************************************************************80
!
!! TEST01 tests CFFT1B, CFFT1F and CFFT1I.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 November 2011
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 4096

  complex ( kind = 8 ) c(n)
  integer ( kind = 4 ) ier
  integer ( kind = 4 ) inc
  integer ( kind = 4 ) lenc
  integer ( kind = 4 ) lensav
  integer ( kind = 4 ) lenwrk
  integer ( kind = 4 ) seed
  real ( kind = 8 ), allocatable, dimension ( : ) :: work
  real ( kind = 8 ), allocatable, dimension ( : ) :: wsave

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST01'
  write ( *, '(a)' ) '  For complex double precision fast Fourier transforms, 1D,'
  write ( *, '(a)' ) '  CFFT1I initializes the transform,'
  write ( *, '(a)' ) '  CFFT1F does a forward transform;'
  write ( *, '(a)' ) '  CFFT1B does a backward transform.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  The number of data items is N = ', n
!
!  Allocate the work arrays.
!
  lenwrk = 2 * n
  lensav = 2 * n + int ( log ( real ( n, kind = 8 ) ) / log ( 2.0D+00 ) ) + 4

  write ( *, '(a,i8)' ) '  LENSAV = ', lensav
  write ( *, '(a,i8)' ) '  LENWRK = ', lenwrk

  allocate ( work(1:lenwrk) )
  allocate ( wsave(1:lensav) )

  call cfft1i ( n, wsave, lensav, ier )
!
!  Set the data values.
!
  seed = 1973

  call c8vec_uniform_01 ( n, seed, c )

  call c8vec_print_part ( n, c, 10, '  The original data:' ) 
!
!  Compute the FFT coefficients.
!
  inc = 1
  lenc = n

  call cfft1f ( n, inc, c, lenc, wsave, lensav, work, lenwrk, ier )

  call c8vec_print_part ( n, c, 10, '  The FFT coefficients:' )
!
!  Compute inverse FFT of coefficients.  Should get back the
!  original data.
!
  call cfft1b ( n, inc, c, lenc, wsave, lensav, work, lenwrk, ier )

  call c8vec_print_part ( n, c, 10, '  The retrieved data:' )

  deallocate ( work )
  deallocate ( wsave )

  return
end
subroutine test02 ( )

!*****************************************************************************80
!
!! TEST02 tests CFFT2B, CFFT2F and CFFT2I.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 November 2011
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: l = 32
  integer ( kind = 4 ), parameter :: m = 64

  complex ( kind = 8 ) c(l,m)
  integer ( kind = 4 ) ier
  integer ( kind = 4 ) ldim
  integer ( kind = 4 ) lensav
  integer ( kind = 4 ) lenwrk
  integer ( kind = 4 ) seed
  real ( kind = 8 ), allocatable, dimension ( : ) :: work
  real ( kind = 8 ), allocatable, dimension ( : ) :: wsave

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST02'
  write ( *, '(a)' ) '  For complex double precision fast Fourier transforms, 2D,'
  write ( *, '(a)' ) '  CFFT2I initializes the transform,'
  write ( *, '(a)' ) '  CFFT2F does a forward transform;'
  write ( *, '(a)' ) '  CFFT2B does a backward transform.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The data is stored in an L by M array, with'
  write ( *, '(a,i8)' ) '  L = ', l
  write ( *, '(a,i8)' ) '  M = ', m
!
!  Allocate work arrays.
!
  lenwrk = 2 * l * m

  lensav = 2 * l + int ( log ( real ( l, kind = 8 ) ) / log ( 2.0D+00 ) ) &
    + 2 * m + int ( log ( real ( m, kind = 8 ) ) / log ( 2.0D+00 ) ) &
    + 8 

  write ( *, '(a,i8)' ) '  LENSAV = ', lensav
  write ( *, '(a,i8)' ) '  LENWRK = ', lenwrk

  allocate ( work(1:lenwrk) )
  allocate ( wsave(1:lensav) )

  call cfft2i ( l, m, wsave, lensav, ier )
!
!  Set the data values.
!
  seed = 1973

  call c8mat_uniform_01 ( l, m, seed, c )

  call c8mat_print_some ( l, m, c, 1, 1, 5, 5, &
    '  Part of the original data:' )
!
!  Compute the FFT coefficients.
!
  ldim = l

  call cfft2f ( ldim, l, m, c, wsave, lensav, work, lenwrk, ier )

  call c8mat_print_some ( l, m, c, 1, 1, 5, 5, &
    '  Part of the FFT coefficients:' )
!
!  Compute inverse FFT of coefficients.  Should get back the
!  original data.
!
  call cfft2b ( ldim, l, m, c, wsave, lensav, work, lenwrk, ier )

  call c8mat_print_some ( l, m, c, 1, 1, 5, 5, '  Part of the retrieved data:' )

  deallocate ( work )
  deallocate ( wsave )

  return
end
subroutine test03 ( )

!*****************************************************************************80
!
!! TEST03 tests CFFTMB, CFFTMF and CFFTMI.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 November 2011
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 32
  integer ( kind = 4 ), parameter :: lot = 6

  complex ( kind = 8 ), allocatable, dimension ( : ) :: c
  integer ( kind = 4 ) ier
  integer ( kind = 4 ) inc
  integer ( kind = 4 ) jump
  integer ( kind = 4 ) lenc
  integer ( kind = 4 ) lensav
  integer ( kind = 4 ) lenwrk
  integer ( kind = 4 ) seed
  real ( kind = 8 ), allocatable, dimension ( : ) :: work
  real ( kind = 8 ), allocatable, dimension ( : ) :: wsave

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST03'
  write ( *, '(a)' ) '  For complex double precision fast Fourier transforms, 1D, multiple'
  write ( *, '(a)' ) '  CFFTMI initializes the transform,'
  write ( *, '(a)' ) '  CFFTMF does a forward transform;'
  write ( *, '(a)' ) '  CFFTMB does a backward transform.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  The number of sequences is LOT = ', lot
  write ( *, '(a,i8)' ) '  The length of each sequence is N = ', n
!
!  Set work vectors.
!
  lenc = n * lot
  lenwrk = 2 * lot * n
  lensav = 2 * n + int ( log ( real ( n, kind = 8 ) ) / log ( 2.0D+00 ) ) + 4

  write ( *, '(a,i8)' ) '  LENC   = ', lenc
  write ( *, '(a,i8)' ) '  LENSAV = ', lensav
  write ( *, '(a,i8)' ) '  LENWRK = ', lenwrk

  allocate ( c(1:lenc) )
  allocate ( wsave(1:lensav) )
  allocate ( work(1:lenwrk) )

  call cfftmi ( n, wsave, lensav, ier )
!
!  Set the data values.
!
  seed = 1973

  call c8mat_uniform_01 ( n, lot, seed, c )

  call c8mat_print_some ( n, lot, c, 1, 1, 5, 5, &
    '  Part of the original data:' )
!
!  Compute the FFT coefficients.
!
  jump = n
  inc = 1

  call cfftmf ( lot, jump, n, inc, c, lenc, wsave, lensav, work, lenwrk, ier )

  call c8mat_print_some ( n, lot, c, 1, 1, 5, 5, &
    '  Part of the FFT coefficients:' )
!
!  Compute inverse FFT of coefficients.  Should get back the
!  original data.
!
  call cfftmb ( lot, jump, n, inc, c, lenc, wsave, lensav, work, lenwrk, ier )

  call c8mat_print_some ( n, lot, c, 1, 1, 5, 5, &
    '  Part of the retrieved data:' )

  deallocate ( c )
  deallocate ( wsave )
  deallocate ( work )

  return
end
subroutine test04 ( )

!*****************************************************************************80
!
!! TEST04 tests COSQ1B, COSQ1F and COSQ1I.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 November 2011
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 4096

  integer ( kind = 4 ) ier
  integer ( kind = 4 ) inc
  integer ( kind = 4 ) lenr
  integer ( kind = 4 ) lensav
  integer ( kind = 4 ) lenwrk
  real ( kind = 8 ) r(n)
  integer ( kind = 4 ) seed
  real ( kind = 8 ), allocatable, dimension ( : ) :: work
  real ( kind = 8 ), allocatable, dimension ( : ) :: wsave

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST04'
  write ( *, '(a)' ) '  For real double precision fast cosine transforms, 1D,'
  write ( *, '(a)' ) '  COSQ1I initializes the transform,'
  write ( *, '(a)' ) '  COSQ1F does a forward transform;'
  write ( *, '(a)' ) '  COSQ1B does a backward transform.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  The number of data items is N = ', n
!
!  Set work vectors.
!
  lensav = 2 * n + int ( log ( real ( n, kind = 8 ) ) / log ( 2.0D+00 ) ) + 4
  lenwrk = n

  write ( *, '(a,i8)' ) '  LENSAV = ', lensav
  write ( *, '(a,i8)' ) '  LENWRK = ', lenwrk

  allocate ( work(1:lenwrk) )
  allocate ( wsave(1:lensav) )

  call cosq1i ( n, wsave, lensav, ier )
!
!  Set the data values.
!
  seed = 1973

  call r8vec_uniform_01 ( n, seed, r )

  call r8vec_print_part ( n, r, 10, '  The original data:' )
!
!  Compute the FFT coefficients.
!
  inc = 1
  lenr = n

  call cosq1f ( n, inc, r, lenr, wsave, lensav, work, lenwrk, ier )

  call r8vec_print_part ( n, r, 10, '  The FFT coefficients:' )
!
!  Compute inverse FFT of coefficients.  Should get back the
!  original data.
!
  call cosq1b ( n, inc, r, lenr, wsave, lensav, work, lenwrk, ier )

  call r8vec_print_part ( n, r, 10, '  The retrieved data:' )

  deallocate ( work )
  deallocate ( wsave )

  return
end
subroutine test05 ( )

!*****************************************************************************80
!
!! TEST05 tests COSQMB, COSQMF and COSQMI.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 November 2011
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 32
  integer ( kind = 4 ), parameter :: lot = 6    

  integer ( kind = 4 ) ier
  integer ( kind = 4 ) inc
  integer ( kind = 4 ) jump
  integer ( kind = 4 ) lenr
  integer ( kind = 4 ) lensav
  integer ( kind = 4 ) lenwrk
  real ( kind = 8 ), allocatable, dimension ( : ) :: r
  integer ( kind = 4 ) seed
  real ( kind = 8 ), allocatable, dimension ( : ) :: work
  real ( kind = 8 ), allocatable, dimension ( : ) :: wsave

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST05'
  write ( *, '(a)' ) '  For real double precision fast cosine transform, 1D, multiple'
  write ( *, '(a)' ) '  COSQMI initializes the transform,'
  write ( *, '(a)' ) '  COSQMF does a forward transform;'
  write ( *, '(a)' ) '  COSQMB does a backward transform.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  The number of sequences is LOT =   ', lot
  write ( *, '(a,i8)' ) '  The length of each sequence is N = ', n
!
!  Set work arrays.
!
  lenr = n * lot
  lenwrk = lot * n
  lensav = 2 * n + int ( log ( real ( n, kind = 8 ) ) / log ( 2.0D+00 ) ) + 4

  write ( *, '(a,i8)' ) '  LENR   = ', lenr
  write ( *, '(a,i8)' ) '  LENSAV = ', lensav
  write ( *, '(a,i8)' ) '  LENWRK = ', lenwrk

  allocate ( r(lenr) )
  allocate ( work(1:lenwrk) )
  allocate ( wsave(1:lensav) )

  call cosqmi ( n, wsave, lensav, ier )
!
!  Set the data values.
!
  seed = 1973

  call r8mat_uniform_01 ( n, lot, seed, r )

  call r8mat_print_some ( n, lot, r, 1, 1, 5, 5, &
    '  Part of the original data:' )
!
!  Compute the FFT coefficients.
!
  jump = n
  inc = 1

  call cosqmf ( lot, jump, n, inc, r, lenr, wsave, lensav, work, lenwrk, ier )

  call r8mat_print_some ( n, lot, r, 1, 1, 5, 5, &
    '  Part of the FFT coefficients:' )
!
!  Compute inverse FFT of coefficients.  Should get back the
!  original data.
!
  call cosqmb ( lot, jump, n, inc, r, lenr, wsave, lensav, work, lenwrk, ier )

  call r8mat_print_some ( n, lot, r, 1, 1, 5, 5, &
    '  Part of the retrieved data:' )

  deallocate ( r )
  deallocate ( work )
  deallocate ( wsave )

  return
end
subroutine test06 ( )

!*****************************************************************************80
!
!! TEST06 tests COST1B, COST1F and COST1I.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 November 2011
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 4096 

  integer ( kind = 4 ) ier
  integer ( kind = 4 ) inc
  integer ( kind = 4 ) lenr
  integer ( kind = 4 ) lensav
  integer ( kind = 4 ) lenwrk
  real ( kind = 8 ) r(n)
  integer ( kind = 4 ) seed
  real ( kind = 8 ), allocatable, dimension ( : ) :: work
  real ( kind = 8 ), allocatable, dimension ( : ) :: wsave

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST06'
  write ( *, '(a)' ) '  For real double precision fast cosine transforms, 1D,'
  write ( *, '(a)' ) '  COST1I initializes the transform,'
  write ( *, '(a)' ) '  COST1F does a forward transform;'
  write ( *, '(a)' ) '  COST1B does a backward transform.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  The number of data items is N = ', n
!
!  Set work arrays.
!
  lenwrk = n - 1 
  lensav = 2 * n + int ( log ( real ( n, kind = 8 ) ) / log ( 2.0D+00 ) ) + 4 

  write ( *, '(a,i8)' ) '  LENSAV = ', lensav
  write ( *, '(a,i8)' ) '  LENWRK = ', lenwrk

  allocate ( work(1:lenwrk) )
  allocate ( wsave(1:lensav) )

  call cost1i ( n, wsave, lensav, ier )
!
!  Set the data values.
!
  seed = 1973

  call r8vec_uniform_01 ( n, seed, r )

  call r8vec_print_part ( n, r, 10, '  The original data:' )
!
!  Compute the FFT coefficients.
!
  inc = 1
  lenr = n

  call cost1f ( n, inc, r, lenr, wsave, lensav, work, lenwrk, ier )

  call r8vec_print_part ( n, r, 10, '  The FFT coefficients:' )
!
!  Compute inverse FFT of coefficients.
!
  call cost1b ( n, inc, r, lenr, wsave, lensav, work, lenwrk, ier )

  call r8vec_print_part ( n, r, 10, '  The retrieved data:' )

  deallocate ( work )
  deallocate ( wsave )

  return
end
subroutine test07 ( )

!*****************************************************************************80
!
!! TEST07 tests COSTMB, COSTMF and COSTMI.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 November 2011
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 32
  integer ( kind = 4 ), parameter :: lot = 6

  integer ( kind = 4 ) ier
  integer ( kind = 4 ) inc
  integer ( kind = 4 ) jump
  integer ( kind = 4 ) lenr
  integer ( kind = 4 ) lensav
  integer ( kind = 4 ) lenwrk
  real ( kind = 8 ), allocatable, dimension ( : ) :: r
  integer ( kind = 4 ) seed
  real ( kind = 8 ), allocatable, dimension ( : ) :: work
  real ( kind = 8 ), allocatable, dimension ( : ) :: wsave

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST07'
  write ( *, '(a)' ) '  For real double precision fast cosine transforms, 1D, multiple'
  write ( *, '(a)' ) '  COSTMI initializes the transform,'
  write ( *, '(a)' ) '  COSTMF does a forward transform;'
  write ( *, '(a)' ) '  COSTMB does a backward transform.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  The number of sequences is LOT =   ', lot
  write ( *, '(a,i8)' ) '  The length of each sequence is N = ', n
!
!  Set work arrays.
!
  lenr = n * lot
  lensav = 2 * n + int ( log ( real ( n, kind = 8 ) ) / log ( 2.0D+00 ) ) + 4
  lenwrk = lot * ( n + 1 )

  write ( *, '(a,i8)' ) '  LENR   = ', lenr
  write ( *, '(a,i8)' ) '  LENSAV = ', lensav
  write ( *, '(a,i8)' ) '  LENWRK = ', lenwrk

  allocate ( r(1:lenr) )
  allocate ( wsave(1:lensav) )
  allocate ( work(1:lenwrk) )

  call costmi ( n, wsave, lensav, ier )
!
!  Set the data values.
!
  seed = 1973

  call r8mat_uniform_01 ( n, lot, seed, r )

  call r8mat_print_some ( n, lot, r, 1, 1, 5, 5, &
    '  Part of the original data:' )
!
!  Compute the FFT coefficients.
!
  jump = n
  inc = 1

  call costmf ( lot, jump, n, inc, r, lenr, wsave, lensav, work, lenwrk, ier )

  call r8mat_print_some ( n, lot, r, 1, 1, 5, 5, &
    '  Part of the FFT coefficients:' )
!
!  Compute inverse FFT of coefficients.  Should get back the
!  original data.
!
  call costmb ( lot, jump, n, inc, r, lenr, wsave, lensav, work, lenwrk, ier )

  call r8mat_print_some ( n, lot, r, 1, 1, 5, 5, &
    '  Part of the retrieved data:' )

  deallocate ( r )
  deallocate ( work )
  deallocate ( wsave )

  return
end
subroutine test09 ( )

!*****************************************************************************80
!
!! TEST09 tests RFFT1B, RFFT1F and RFFT1I.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 November 2011
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 4096

  integer ( kind = 4 ) ier
  integer ( kind = 4 ) inc
  integer ( kind = 4 ) lenr
  integer ( kind = 4 ) lensav
  integer ( kind = 4 ) lenwrk
  real ( kind = 8 ) r(n)
  integer ( kind = 4 ) seed
  real ( kind = 8 ), allocatable, dimension ( : ) :: work
  real ( kind = 8 ), allocatable, dimension ( : ) :: wsave

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST09'
  write ( *, '(a)' ) '  For real double precision fast cosine transforms, 1D,'
  write ( *, '(a)' ) '  RFFT1I initializes the transform,'
  write ( *, '(a)' ) '  RFFT1F does a forward transform;'
  write ( *, '(a)' ) '  RFFT1B does a backward transform.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  The number of data items is N = ', n
!
!  Set work vectors.
!
  lensav = n + int ( log ( real ( n, kind = 8 ) ) / log ( 2.0D+00 ) ) + 4
  lenwrk = n

  write ( *, '(a,i8)' ) '  LENSAV = ', lensav
  write ( *, '(a,i8)' ) '  LENWRK = ', lenwrk

  allocate ( work(1:lenwrk) )
  allocate ( wsave(1:lensav) )

  call rfft1i ( n, wsave, lensav, ier )
!
!  Set the data values.
!
  seed = 1973

  call r8vec_uniform_01 ( n, seed, r )

  call r8vec_print_part ( n, r, 10, '  The original data:' )
!
!  Compute the FFT coefficients.
!
  inc = 1
  lenr = n

  call rfft1f ( n, inc, r, lenr, wsave, lensav, work, lenwrk, ier )

  call r8vec_print_part ( n, r, 10, '  The FFT coefficients:' )
!
!  Compute inverse FFT of coefficients.  Should get back the
!  original data.
!
  call rfft1b ( n, inc, r, lenr, wsave, lensav, work, lenwrk, ier )

  call r8vec_print_part ( n, r, 10, '  The retrieved data:' )

  deallocate ( work )
  deallocate ( wsave )

  return
end
subroutine test10 ( )

!*****************************************************************************80
!
!! TEST10 tests RFFT2B, RFFT2F and RFFT2I.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 November 2011
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: l = 32
  integer ( kind = 4 ), parameter :: m = 64

  integer ( kind = 4 ), parameter :: ldim = 2 * ( l / 2 + 1 )

  integer ( kind = 4 ) ier
  integer ( kind = 4 ) lensav
  integer ( kind = 4 ) lenwrk
  real ( kind = 8 ) r(ldim,m)
  integer ( kind = 4 ) seed
  real ( kind = 8 ), allocatable, dimension ( : ) :: work
  real ( kind = 8 ), allocatable, dimension ( : ) :: wsave

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST10'
  write ( *, '(a)' ) '  For real double precision fast Fourier transform, 2D,'
  write ( *, '(a)' ) '  RFFT2I initializes the transform,'
  write ( *, '(a)' ) '  RFFT2F does a forward transform;'
  write ( *, '(a)' ) '  RFFT2B does a backward transform.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The L by M data is stored in an LDIM by M array, with'
  write ( *, '(a,i8)' ) '  L =    ', l
  write ( *, '(a,i8)' ) '  LDIM = ', ldim
  write ( *, '(a,i8)' ) '  M =    ', m
!
!  Set work arrays.
!
  lenwrk = 2 * ldim * m

  lensav = &
          l + int ( log ( real ( l, kind = 8 ) ) / log ( 2.0D+00 ) ) + 4 &
    + 2 * m + int ( log ( real ( m, kind = 8 ) ) / log ( 2.0D+00 ) ) + 4 &
    +     m + int ( log ( real ( m, kind = 8 ) ) / log ( 2.0D+00 ) ) + 4

  write ( *, '(a,i8)' ) '  LENSAV = ', lensav
  write ( *, '(a,i8)' ) '  LENWRK = ', lenwrk

  allocate ( work(1:lenwrk) )
  allocate ( wsave(1:lensav) )

  call rfft2i ( l, m, wsave, lensav, ier )
!
!  Set the data values.
!
  seed = 1973

  call r8mat_uniform_01 ( ldim, m, seed, r )

  call r8mat_print_some ( ldim, m, r, 1, 1, 5, 5, &
    '  Part of the original data:' )
!
!  Compute the FFT coefficients.
!
  call rfft2f ( ldim, l, m, r, wsave, lensav, work, lenwrk, ier )

  call r8mat_print_some ( ldim, m, r, 1, 1, 5, 5, &
    '  Part of the FFT coefficients:' )
!
!  Compute inverse FFT of coefficients.  Should get back the
!  original data.
!
  call rfft2b ( ldim, l, m, r, wsave, lensav, work, lenwrk, ier )

  call r8mat_print_some ( ldim, m, r, 1, 1, 5, 5, &
    '  Part of the retrieved data:' )

  deallocate ( work )
  deallocate ( wsave )

  return
end
subroutine test11 ( )

!*****************************************************************************80
!
!! TEST11 tests RFFTMB, RFFTMF and RFFTMI.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 November 2011
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 32
  integer ( kind = 4 ), parameter :: lot = 6

  integer ( kind = 4 ) ier
  integer ( kind = 4 ) inc
  integer ( kind = 4 ) jump
  integer ( kind = 4 ) lenr
  integer ( kind = 4 ) lensav
  integer ( kind = 4 ) lenwrk
  real ( kind = 8 ), allocatable, dimension ( : ) :: r
  integer ( kind = 4 ) seed
  real ( kind = 8 ), allocatable, dimension ( : ) :: work
  real ( kind = 8 ), allocatable, dimension ( : ) :: wsave

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST11'
  write ( *, '(a)' ) '  For real double precision fast Fourier transform, 1D, multiple'
  write ( *, '(a)' ) '  RFFTMI initializes the transform,'
  write ( *, '(a)' ) '  RFFTMF does a forward transform;'
  write ( *, '(a)' ) '  RFFTMB does a backward transform.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  The number of sequences is LOT =   ', lot
  write ( *, '(a,i8)' ) '  The length of each sequence is N = ', n
!
!  Set work arrays.
!
  lenr = n * lot
  lenwrk = lot * n
  lensav = n + int ( log ( real ( n, kind = 8 ) ) / log ( 2.0D+00 ) ) + 4

  write ( *, '(a,i8)' ) '  LENR =   ', lenr
  write ( *, '(a,i8)' ) '  LENSAV = ', lensav
  write ( *, '(a,i8)' ) '  LENWRK = ', lenwrk

  allocate ( r(1:lenr) )
  allocate ( wsave(1:lensav) )
  allocate ( work(1:lenwrk) )

  call rfftmi ( n, wsave, lensav, ier )
!
!  Set the data values.
!
  seed = 1973

  call r8mat_uniform_01 ( n, lot, seed, r )

  call r8mat_print_some ( n, lot, r, 1, 1, 5, 5, '  Part of the original data:' )
!
!  Compute the FFT coefficients.
!
  jump = n
  inc = 1

  call rfftmf ( lot, jump, n, inc, r, lenr, wsave, lensav, work, lenwrk, ier )

  call r8mat_print_some ( n, lot, r, 1, 1, 5, 5, &
    '  Part of the FFT coefficients:' )
!
!  Compute inverse FFT of coefficients.  Should get back the
!  original data.
!
  call rfftmb ( lot, jump, n, inc, r, lenr, wsave, lensav, work, lenwrk, ier )

  call r8mat_print_some ( n, lot, r, 1, 1, 5, 5, &
    '  Part of the retrieved data:' )

  deallocate ( r )
  deallocate ( work )
  deallocate ( wsave )

  return
end
subroutine test12 ( )

!*****************************************************************************80
!
!! TEST12 tests SINQ1B, SINQ1F and SINQ1I.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 November 2011
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 4096

  integer ( kind = 4 ) ier
  integer ( kind = 4 ) inc
  integer ( kind = 4 ) lenr
  integer ( kind = 4 ) lensav
  integer ( kind = 4 ) lenwrk
  real ( kind = 8 ) r(n)
  integer ( kind = 4 ) seed
  real ( kind = 8 ), allocatable, dimension ( : ) :: work
  real ( kind = 8 ), allocatable, dimension ( : ) :: wsave

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST12'
  write ( *, '(a)' ) '  For real double precision fast sine transforms, 1D,'
  write ( *, '(a)' ) '  SINQ1I initializes the transform,'
  write ( *, '(a)' ) '  SINQ1F does a forward transform;'
  write ( *, '(a)' ) '  SINQ1B does a backward transform.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  The number of data items is N = ', n
!
!  Set work arrays.
!
  lenwrk = n
  lensav = 2 * n + int ( log ( real ( n, kind = 8 ) ) / log ( 2.0D+00 ) ) + 4

  write ( *, '(a,i8)' ) '  LENSAV = ', lensav
  write ( *, '(a,i8)' ) '  LENWRK = ', lenwrk

  allocate ( work(1:lenwrk) )
  allocate ( wsave(1:lensav) )

  call sinq1i ( n, wsave, lensav, ier )
!
!  Set the data values.
!
  seed = 1973

  call r8vec_uniform_01 ( n, seed, r )

  call r8vec_print_part ( n, r, 10, '  The original data:' )
!
!  Compute the FFT coefficients.
!
  inc = 1
  lenr = n

  call sinq1f ( n, inc, r, lenr, wsave, lensav, work, lenwrk, ier )

  call r8vec_print_part ( n, r, 10, '  The FFT coefficients:' )
!
!  Compute inverse FFT of coefficients.  Should get back the
!  original data.
!
  call sinq1b ( n, inc, r, lenr, wsave, lensav, work, lenwrk, ier )

  call r8vec_print_part ( n, r, 10, '  The retrieved data:' )

  deallocate ( work )
  deallocate ( wsave )

  return
end
subroutine test13 ( )

!*****************************************************************************80
!
!! TEST13 tests SINQMB, SINQMF and SINQMI.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 November 2011
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 32
  integer ( kind = 4 ), parameter :: lot = 6

  integer ( kind = 4 ) ier
  integer ( kind = 4 ) inc
  integer ( kind = 4 ) jump
  integer ( kind = 4 ) lenr
  integer ( kind = 4 ) lensav
  integer ( kind = 4 ) lenwrk
  real ( kind = 8 ), allocatable, dimension ( : ) :: r
  integer ( kind = 4 ) seed
  real ( kind = 8 ), allocatable, dimension ( : ) :: work
  real ( kind = 8 ), allocatable, dimension ( : ) :: wsave

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST13'
  write ( *, '(a)' ) '  For real double precision fast sine transforms, 1D, multiple'
  write ( *, '(a)' ) '  SINQMI initializes the transform,'
  write ( *, '(a)' ) '  SINQMF does a forward transform;'
  write ( *, '(a)' ) '  SINQMB does a backward transform.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  The number of sequences is LOT =   ', lot
  write ( *, '(a,i8)' ) '  The length of each sequence is N = ', n
!
!  Set work arrays.
!
  lenr = n * lot
  lenwrk = lot * n
  lensav = 2 * n + int ( log ( real ( n, kind = 8 ) ) / log ( 2.0D+00 ) ) + 4

  write ( *, '(a,i8)' ) '  LENR   = ', lenr
  write ( *, '(a,i8)' ) '  LENSAV = ', lensav
  write ( *, '(a,i8)' ) '  LENWRK = ', lenwrk

  allocate ( r(1:lenr) )
  allocate ( wsave(1:lensav) )
  allocate ( work(1:lenwrk) )

  call sinqmi ( n, wsave, lensav, ier )
!
!  Set the data values.
!
  seed = 1973

  call r8mat_uniform_01 ( n, lot, seed, r )

  call r8mat_print_some ( n, lot, r, 1, 1, 5, 5, &
    '  Part of the original data:' )
!
!  Compute the FFT coefficients.
!
  jump = n
  inc = 1

  call sinqmf ( lot, jump, n, inc, r, lenr, wsave, lensav, work, lenwrk, ier )

  call r8mat_print_some ( n, lot, r, 1, 1, 5, 5, &
    '  Part of the FFT coefficients:' )
!
!  Compute inverse FFT of coefficients.  Should get back the
!  original data.
!
  call sinqmb ( lot, jump, n, inc, r, lenr, wsave, lensav, work, lenwrk, ier )

  call r8mat_print_some ( n, lot, r, 1, 1, 5, 5, &
    '  Part of the retrieved data:' )

  deallocate ( r )
  deallocate ( work )
  deallocate ( wsave )

  return
end
subroutine test14 ( )

!*****************************************************************************80
!
!! TEST14 tests SINT1B, SINT1F and SINT1I.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 November 2011
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 4096

  integer ( kind = 4 ) ier
  integer ( kind = 4 ) inc
  integer ( kind = 4 ) lenr
  integer ( kind = 4 ) lensav
  integer ( kind = 4 ) lenwrk
  real ( kind = 8 ) r(n)
  integer ( kind = 4 ) seed
  real ( kind = 8 ), allocatable, dimension ( : ) :: work
  real ( kind = 8 ), allocatable, dimension ( : ) :: wsave

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST14'
  write ( *, '(a)' ) '  For real double precision fast sine transforms, 1D,'
  write ( *, '(a)' ) '  SINT1I initializes the transform,'
  write ( *, '(a)' ) '  SINT1F does a forward transform;'
  write ( *, '(a)' ) '  SINT1B does a backward transform.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  The number of data items is N = ', n
!
!  Set work arrays.
!
  lenwrk = 2 * ( n + 1 )
  lensav = n / 2 + n + int ( log ( real ( n, kind = 8 ) ) / log ( 2.0D+00 ) ) + 4

  write ( *, '(a,i8)' ) '  LENSAV = ', lensav
  write ( *, '(a,i8)' ) '  LENWRK = ', lenwrk

  allocate ( wsave(1:lensav) )
  allocate ( work(1:lenwrk) )

  call sint1i ( n, wsave, lensav, ier )
!
!  Set the data values.
!
  seed = 1973

  call r8vec_uniform_01 ( n, seed, r )

  call r8vec_print_part ( n, r, 10, '  The original data:' )
!
!  Compute the FFT coefficients.
!
  inc = 1
  lenr = n

  call sint1f ( n, inc, r, lenr, wsave, lensav, work, lenwrk, ier )

  call r8vec_print_part ( n, r, 10, '  The FFT coefficients:' )
!
!  Compute inverse FFT of coefficients.  Should get back the
!  original data.
!
  call sint1b ( n, inc, r, lenr, wsave, lensav, work, lenwrk, ier )

  call r8vec_print_part ( n, r, 10, '  The retrieved data:' )

  deallocate ( work )
  deallocate ( wsave )

  return
end
subroutine test15 ( )

!*****************************************************************************80
!
!! TEST15 tests SINTMB, SINTMF and SINTMI.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 November 2011
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 32
  integer ( kind = 4 ), parameter :: lot = 6

  integer ( kind = 4 ) ier
  integer ( kind = 4 ) inc
  integer ( kind = 4 ) jump
  integer ( kind = 4 ) lenr
  integer ( kind = 4 ) lensav
  integer ( kind = 4 ) lenwrk
  real ( kind = 8 ), allocatable, dimension ( : ) :: r
  integer ( kind = 4 ) seed
  real ( kind = 8 ), allocatable, dimension ( : ) :: work
  real ( kind = 8 ), allocatable, dimension ( : ) :: wsave

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST15'
  write ( *, '(a)' ) '  For real double precision fast sine transforms, 1D, multiple'
  write ( *, '(a)' ) '  SINTMI initializes the transform,'
  write ( *, '(a)' ) '  SINTMF does a forward transform;'
  write ( *, '(a)' ) '  SINTMB does a backward transform.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  The number of sequences is LOT =   ', lot
  write ( *, '(a,i8)' ) '  The length of each sequence is N = ', n
!
!  Set work arrays.
!
  lenr = n * lot
  lenwrk = lot * 2 * ( n + 2 )
  lensav = n / 2 + n + int ( log ( real ( n, kind = 8 ) ) / log ( 2.0D+00 ) ) + 4

  write ( *, '(a,i8)' ) '  LENR =   ', lenr
  write ( *, '(a,i8)' ) '  LENSAV = ', lensav
  write ( *, '(a,i8)' ) '  LENWRK = ', lenwrk

  allocate ( r(1:lenr) )
  allocate ( wsave(1:lensav) )
  allocate ( work(1:lenwrk) )

  call sintmi ( n, wsave, lensav, ier )
!
!  Set the data values.
!
  seed = 1973

  call r8mat_uniform_01 ( n, lot, seed, r )

  call r8mat_print_some ( n, lot, r, 1, 1, 5, 5, &
    '  Part of the original data:' )
!
!  Compute the FFT coefficients.
!
  jump = n
  inc = 1

  call sintmf ( lot, jump, n, inc, r, lenr, wsave, lensav, &
    work, lenwrk, ier )

  call r8mat_print_some ( n, lot, r, 1, 1, 5, 5, &
    '  Part of the FFT coefficients:' )
!
!  Compute inverse FFT of coefficients.  Should get back the
!  original data.
!
  call sintmb ( lot, jump, n, inc, r, lenr, wsave, lensav, &
    work, lenwrk, ier )

  call r8mat_print_some ( n, lot, r, 1, 1, 5, 5, &
    '  Part of the retrieved data:' )

  deallocate ( r )
  deallocate ( work )
  deallocate ( wsave )

  return
end
subroutine c8mat_print_some ( m, n, a, ilo, jlo, ihi, jhi, title )

!*****************************************************************************80
!
!! C8MAT_PRINT_SOME prints some of a C8MAT.
!
!  Discussion:
!
!    A C8MAT is a matrix of C8's.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    23 March 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns 
!    in the matrix.
!
!    Input, complex ( kind = 8 ) A(M,N), the matrix.
!
!    Input, integer ( kind = 4 ) ILO, JLO, IHI, JHI, the first row and
!    column, and the last row and column to be printed.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer ( kind = 4 ), parameter :: incx = 4
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  complex ( kind = 8 ) a(m,n)
  character ( len = 20 ) ctemp(incx)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i2hi
  integer ( kind = 4 ) i2lo
  integer ( kind = 4 ) ihi
  integer ( kind = 4 ) ilo
  integer ( kind = 4 ) inc
  integer ( kind = 4 ) j
  integer ( kind = 4 ) j2
  integer ( kind = 4 ) j2hi
  integer ( kind = 4 ) j2lo
  integer ( kind = 4 ) jhi
  integer ( kind = 4 ) jlo
  character ( len = * )  title
  complex ( kind = 8 ) zero

  zero = cmplx ( 0.0D+00, 0.0D+00, kind = 8 )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )

  if ( m <= 0 .or. n <= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  (None)' 
    return
  end if
!
!  Print the columns of the matrix, in strips of INCX.
!
  do j2lo = jlo, min ( jhi, n ), incx

    j2hi = j2lo + incx - 1
    j2hi = min ( j2hi, n )
    j2hi = min ( j2hi, jhi )

    inc = j2hi + 1 - j2lo

    write ( *, '(a)' ) ' '

    do j = j2lo, j2hi
      j2 = j + 1 - j2lo
      write ( ctemp(j2), '(i10,10x)' ) j
    end do

    write ( *, '(a,4a20)' ) '  Col: ', ( ctemp(j2), j2 = 1, inc )
    write ( *, '(a)' ) '  Row'
    write ( *, '(a)' ) '  ---'
!
!  Determine the range of the rows in this strip.
!
    i2lo = max ( ilo, 1 )
    i2hi = min ( ihi, m )

    do i = i2lo, i2hi
!
!  Print out (up to) INCX entries in row I, that lie in the current strip.
!
      do j2 = 1, inc

        j = j2lo - 1 + j2

        if ( a(i,j) == zero ) then
          ctemp(j2) = '       0.0          '
        else if ( imag ( a(i,j) ) == 0.0D+00 ) then
          write ( ctemp(j2), '(g10.3,10x)' ) real ( a(i,j), kind = 8 )
        else
          write ( ctemp(j2), '(2g10.3)' ) a(i,j)
        end if

      end do

      write ( *, '(i5,a1,4a20)' ) i, ':', ( ctemp(j2), j2 = 1, inc )

    end do

  end do

  return
end
subroutine c8mat_uniform_01 ( m, n, seed, c )

!*****************************************************************************80
!
!! C8MAT_UNIFORM_01 returns a unit pseudorandom C8MAT.
!
!  Discussion:
!
!    A C8MAT is a matrix of C8's.
!
!    The angles should be uniformly distributed between 0 and 2 * PI,
!    the square roots of the radius uniformly distributed between 0 and 1.
!
!    This results in a uniform distribution of values in the unit circle.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    15 March 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns
!    in the matrix.
!
!    Input/output, integer ( kind = 4 ) SEED, the "seed" value, which 
!    should NOT be 0.  On output, SEED has been updated.
!
!    Output, complex ( kind = 8 ) C(M,N), the pseudorandom complex matrix.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  complex ( kind = 8 ) c(m,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ), parameter :: i4_huge = 2147483647
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  real ( kind = 8 ) r
  real ( kind = 8 ), parameter :: r8_pi = 3.141592653589793D+00
  integer ( kind = 4 ) seed
  real ( kind = 8 ) theta

  do j = 1, n
    do i = 1, m

      k = seed / 127773

      seed = 16807 * ( seed - k * 127773 ) - k * 2836

      if ( seed < 0 ) then
        seed = seed + i4_huge
      end if

      r = sqrt ( real ( seed, kind = 8 ) * 4.656612875D-10 )

      k = seed / 127773

      seed = 16807 * ( seed - k * 127773 ) - k * 2836

      if ( seed < 0 ) then
        seed = seed + i4_huge
      end if

      theta = 2.0D+00 * r8_pi * ( real ( seed, kind = 8 ) * 4.656612875D-10 )

      c(i,j) = r * cmplx ( cos ( theta ), sin ( theta ), kind = 8 )

    end do

  end do

  return
end
subroutine c8vec_print_part ( n, a, max_print, title )

!*****************************************************************************80
!
!! C8VEC_PRINT_PART prints "part" of a C8VEC.
!
!  Discussion:
!
!    The user specifies MAX_PRINT, the maximum number of lines to print.
!
!    If N, the size of the vector, is no more than MAX_PRINT, then
!    the entire vector is printed, one entry per line.
!
!    Otherwise, if possible, the first MAX_PRINT-2 entries are printed,
!    followed by a line of periods suggesting an omission,
!    and the last entry.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    22 June 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries of the vector.
!
!    Input, complex ( kind = 8 ) A(N), the vector to be printed.
!
!    Input, integer ( kind = 4 ) MAX_PRINT, the maximum number of lines
!    to print.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer ( kind = 4 ) n

  complex ( kind = 8 ) a(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) max_print
  character ( len = * )  title

  if ( max_print <= 0 ) then
    return
  end if

  if ( n <= 0 ) then
    return
  end if

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )
  write ( *, '(a)' ) ' '

  if ( n <= max_print ) then

    do i = 1, n
      write ( *, '(2x,i8,a,1x,g14.6,2x,g14.6)' ) i, ':', a(i)
    end do

  else if ( 3 <= max_print ) then

    do i = 1, max_print - 2
      write ( *, '(2x,i8,a,1x,g14.6,2x,g14.6)' ) i, ':', a(i)
    end do
    write ( *, '(a)' ) '  ........  ..............  ..............'
    i = n
    write ( *, '(2x,i8,a,1x,g14.6,2x,g14.6)' ) i, ':', a(i)

  else

    do i = 1, max_print - 1
      write ( *, '(2x,i8,a,1x,g14.6,2x,g14.6)' ) i, ':', a(i)
    end do
    i = max_print
    write ( *, '(2x,i8,a,1x,g14.6,2x,g14.6,2x,a)' ) i, ':', a(i), &
      '...more entries...'

  end if

  return
end
subroutine c8vec_uniform_01 ( n, seed, c )

!*****************************************************************************80
!
!! C8VEC_UNIFORM_01 returns a unit pseudorandom C8VEC.
!
!  Discussion:
!
!    A C8VEC is a vector of C8's.
!
!    The angles should be uniformly distributed between 0 and 2 * PI,
!    the square roots of the radius uniformly distributed between 0 and 1.
!
!    This results in a uniform distribution of values in the unit circle.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    15 March 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of values to compute.
!
!    Input/output, integer ( kind = 4 ) SEED, the "seed" value, which should 
!    NOT be 0.  On output, SEED has been updated.
!
!    Output, complex ( kind = 8 ) C(N), the pseudorandom complex vector.
!
  implicit none

  integer ( kind = 4 ) n

  complex ( kind = 8 ) c(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ), parameter :: i4_huge = 2147483647
  integer ( kind = 4 ) k
  real    ( kind = 8 ) r
  real    ( kind = 8 ), parameter :: r8_pi = 3.141592653589793D+00
  integer ( kind = 4 ) seed
  real    ( kind = 8 ) theta

  do i = 1, n

    k = seed / 127773

    seed = 16807 * ( seed - k * 127773 ) - k * 2836

    if ( seed < 0 ) then
      seed = seed + i4_huge
    end if

    r = sqrt ( real ( seed, kind = 8 ) * 4.656612875D-10 )

    k = seed / 127773

    seed = 16807 * ( seed - k * 127773 ) - k * 2836

    if ( seed < 0 ) then
      seed = seed + i4_huge
    end if

    theta = 2.0D+00 * r8_pi * ( real ( seed, kind = 8 ) * 4.656612875D-10 )

    c(i) = r * cmplx ( cos ( theta ), sin ( theta ), kind = 8 )

  end do

  return
end
subroutine r8mat_print ( m, n, a, title )

!*****************************************************************************80
!
!! R8MAT_PRINT prints an R8MAT.
!
!  Discussion:
!
!    An R8MAT is an array of R8 values.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 September 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of rows in A.
!
!    Input, integer ( kind = 4 ) N, the number of columns in A.
!
!    Input, real ( kind = 8 ) A(M,N), the matrix.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(m,n)
  character ( len = * )  title

  call r8mat_print_some ( m, n, a, 1, 1, m, n, title )

  return
end
subroutine r8mat_print_some ( m, n, a, ilo, jlo, ihi, jhi, title )

!*****************************************************************************80
!
!! R8MAT_PRINT_SOME prints some of an R8MAT.
!
!  Discussion:
!
!    An R8MAT is an array of R8 values.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 September 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns.
!
!    Input, real ( kind = 8 ) A(M,N), an M by N matrix to be printed.
!
!    Input, integer ( kind = 4 ) ILO, JLO, the first row and column to print.
!
!    Input, integer ( kind = 4 ) IHI, JHI, the last row and column to print.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer ( kind = 4 ), parameter :: incx = 5
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(m,n)
  character ( len = 14 ) ctemp(incx)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i2hi
  integer ( kind = 4 ) i2lo
  integer ( kind = 4 ) ihi
  integer ( kind = 4 ) ilo
  integer ( kind = 4 ) inc
  integer ( kind = 4 ) j
  integer ( kind = 4 ) j2
  integer ( kind = 4 ) j2hi
  integer ( kind = 4 ) j2lo
  integer ( kind = 4 ) jhi
  integer ( kind = 4 ) jlo
  character ( len = * ) title

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )

  if ( m <= 0 .or. n <= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  (None)'
    return
  end if

  do j2lo = max ( jlo, 1 ), min ( jhi, n ), incx

    j2hi = j2lo + incx - 1
    j2hi = min ( j2hi, n )
    j2hi = min ( j2hi, jhi )

    inc = j2hi + 1 - j2lo

    write ( *, '(a)' ) ' '

    do j = j2lo, j2hi
      j2 = j + 1 - j2lo
      write ( ctemp(j2), '(i8,6x)' ) j
    end do

    write ( *, '(''  Col   '',5a14)' ) ctemp(1:inc)
    write ( *, '(a)' ) '  Row'
    write ( *, '(a)' ) ' '

    i2lo = max ( ilo, 1 )
    i2hi = min ( ihi, m )

    do i = i2lo, i2hi

      do j2 = 1, inc

        j = j2lo - 1 + j2

        if ( a(i,j) == real ( int ( a(i,j) ), kind = 8 ) ) then
          write ( ctemp(j2), '(f8.0,6x)' ) a(i,j)
        else
          write ( ctemp(j2), '(g14.6)' ) a(i,j)
        end if

      end do

      write ( *, '(i5,a,5a14)' ) i, ':', ( ctemp(j), j = 1, inc )

    end do

  end do

  return
end
subroutine r8mat_uniform_01 ( m, n, seed, r )

!*****************************************************************************80
!
!! R8MAT_UNIFORM_01 fills an R8MAT with unit pseudorandom numbers.
!
!  Discussion:
!
!    An R8MAT is an array of R8 values.
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
!    Peter Lewis, Allen Goodman, James Miller,
!    A Pseudo-Random Number Generator for the System/360,
!    IBM Systems Journal,
!    Volume 8, pages 136-143, 1969.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns in
!    the array.
!
!    Input/output, integer ( kind = 4 ) SEED, the "seed" value, which
!    should NOT be 0.  On output, SEED has been updated.
!
!    Output, real ( kind = 8 ) R(M,N), the array of pseudorandom values.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  integer ( kind = 4 ), parameter :: i4_huge = 2147483647
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) seed
  real ( kind = 8 ) r(m,n)

  do j = 1, n

    do i = 1, m

      k = seed / 127773

      seed = 16807 * ( seed - k * 127773 ) - k * 2836

      if ( seed < 0 ) then
        seed = seed + i4_huge
      end if

      r(i,j) = real ( seed, kind = 8 ) * 4.656612875D-10

    end do
  end do

  return
end
subroutine r8vec_print_part ( n, a, max_print, title )

!*****************************************************************************80
!
!! R8VEC_PRINT_PART prints "part" of an R8VEC.
!
!  Discussion:
!
!    The user specifies MAX_PRINT, the maximum number of lines to print.
!
!    If N, the size of the vector, is no more than MAX_PRINT, then
!    the entire vector is printed, one entry per line.
!
!    Otherwise, if possible, the first MAX_PRINT-2 entries are printed,
!    followed by a line of periods suggesting an omission,
!    and the last entry.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 December 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries of the vector.
!
!    Input, real ( kind = 8 ) A(N), the vector to be printed.
!
!    Input, integer ( kind = 4 ) MAX_PRINT, the maximum number of lines
!    to print.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) max_print
  character ( len = * ) title

  if ( max_print <= 0 ) then
    return
  end if

  if ( n <= 0 ) then
    return
  end if

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )
  write ( *, '(a)' ) ' '

  if ( n <= max_print ) then

    do i = 1, n
      write ( *, '(2x,i8,a,1x,g14.6)' ) i, ':', a(i)
    end do

  else if ( 3 <= max_print ) then

    do i = 1, max_print - 2
      write ( *, '(2x,i8,a,1x,g14.6)' ) i, ':', a(i)
    end do
    write ( *, '(a)' ) '  ........  ..............'
    i = n
    write ( *, '(2x,i8,a,1x,g14.6)' ) i, ':', a(i)

  else

    do i = 1, max_print - 1
      write ( *, '(2x,i8,a,1x,g14.6)' ) i, ':', a(i)
    end do
    i = max_print
    write ( *, '(2x,i8,a,1x,g14.6,2x,a)' ) i, ':', a(i), '...more entries...'

  end if

  return
end
subroutine r8vec_uniform_01 ( n, seed, r )

!*****************************************************************************80
!
!! R8VEC_UNIFORM_01 returns a unit pseudorandom R8VEC.
!
!  Discussion:
!
!    An R8VEC is a vector of R8's.
!
!    For now, the input quantity SEED is an integer ( kind = 4 ) variable.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 July 2006
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
!    Peter Lewis, Allen Goodman, James Miller
!    A Pseudo-Random Number Generator for the System/360,
!    IBM Systems Journal,
!    Volume 8, pages 136-143, 1969.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in the vector.
!
!    Input/output, integer ( kind = 4 ) SEED, the "seed" value, which
!    should NOT be 0.  On output, SEED has been updated.
!
!    Output, real ( kind = 8 ) R(N), the vector of pseudorandom values.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  integer ( kind = 4 ) k
  integer ( kind = 4 ) seed
  real ( kind = 8 ) r(n)

  if ( seed == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8VEC_UNIFORM_01 - Fatal error!'
    write ( *, '(a)' ) '  Input value of SEED = 0.'
    stop
  end if

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
!    18 May 2013
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
