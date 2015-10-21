subroutine cfft2i ( l, m, wsave, lensav, ier )

!*****************************************************************************80
!
!! CFFT2I: initialization for CFFT2B and CFFT2F.
!
!  Discussion:
!
!    CFFT2I initializes real array WSAVE for use in its companion
!    routines CFFT2F and CFFT2B for computing two-dimensional fast
!    Fourier transforms of complex data.  Prime factorizations of L and M,
!    together with tabulations of the trigonometric functions, are
!    computed and stored in array WSAVE.
!
!    On 10 May 2010, this code was modified by changing the value
!    of an index into the WSAVE array.
!
!  License:
!
!    Licensed under the GNU General Public License (GPL).
!    Copyright (C) 1995-2004, Scientific Computing Division,
!    University Corporation for Atmospheric Research
!
!  Modified:
!
!    15 November 2011
!
!  Author:
!
!    Original FORTRAN77 version by Paul Swarztrauber, Richard Valent.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Paul Swarztrauber,
!    Vectorizing the Fast Fourier Transforms,
!    in Parallel Computations,
!    edited by G. Rodrigue,
!    Academic Press, 1982.
!
!    Paul Swarztrauber,
!    Fast Fourier Transform Algorithms for Vector Computers,
!    Parallel Computing, pages 45-63, 1984.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) L, the number of elements to be transformed
!    in the first dimension.  The transform is most efficient when L is a
!    product of small primes.
!
!    Input, integer ( kind = 4 ) M, the number of elements to be transformed
!    in the second dimension.  The transform is most efficient when M is a
!    product of small primes.
!
!    Input, integer ( kind = 4 ) LENSAV, the dimension of the WSAVE array.
!    LENSAV must be at least 2*(L+M) + INT(LOG(REAL(L)))
!    + INT(LOG(REAL(M))) + 8.
!
!    Output, real ( kind = 8 ) WSAVE(LENSAV), contains the prime factors of L
!    and M, and also certain trigonometric values which will be used in
!    routines CFFT2B or CFFT2F.
!
!    Output, integer ( kind = 4 ) IER, error flag.
!    0, successful exit;
!    2, input parameter LENSAV not big enough;
!    20, input error returned by lower level routine.
!
  implicit none

  integer ( kind = 4 ) lensav

  integer ( kind = 4 ) ier
  integer ( kind = 4 ) ier1
  integer ( kind = 4 ) iw
  integer ( kind = 4 ) l
  integer ( kind = 4 ) m
  real ( kind = 8 ) wsave(lensav)

  ier = 0

  if ( lensav < 2 * l + int ( log ( real ( l, kind = 8 ) ) &
    / log ( 2.0D+00 ) ) + 2 * m + int ( log ( real ( m, kind = 8 ) ) &
    / log ( 2.0D+00 ) ) + 8 ) then
    ier = 2
    call xerfft ('cfft2i', 4)
    return
  end if

  call cfftmi ( l, wsave(1), 2 * l + int ( log ( real ( l, kind = 8 ) ) &
    / log ( 2.0D+00 ) ) + 4, ier1 )

  if ( ier1 /= 0) then
    ier = 20
    call xerfft ('cfft2i',-5)
    return
  end if

  call cfftmi ( m, &
    wsave(2*l+int(log( real ( l, kind = 8 ) )/log( 2.0D+00 )) + 3), &
    2*m + int(log( real ( m, kind = 8 ) )/log( 2.0D+00 )) + 4, ier1)

  if (ier1 /= 0) then
    ier = 20
    call xerfft ('cfft2i',-5)
  end if

  return
end
