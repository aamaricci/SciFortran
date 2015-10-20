subroutine cosq1b ( n, inc, x, lenx, wsave, lensav, work, lenwrk, ier )

!*****************************************************************************80
!
!! COSQ1B: real double precision backward cosine quarter wave transform, 1D.
!
!  Discussion:
!
!    COSQ1B computes the one-dimensional Fourier transform of a sequence
!    which is a cosine series with odd wave numbers.  This transform is
!    referred to as the backward transform or Fourier synthesis, transforming
!    the sequence from spectral to physical space.
!
!    This transform is normalized since a call to COSQ1B followed
!    by a call to COSQ1F (or vice-versa) reproduces the original
!    array  within roundoff error.
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
!    Input, integer ( kind = 4 ) N, the number of elements to be transformed
!    in the sequence.  The transform is most efficient when N is a
!    product of small primes.
!
!    Input, integer ( kind = 4 ) INC, the increment between the locations,
!    in array R, of two consecutive elements within the sequence.
!
!    Input/output, real ( kind = 8 ) R(LENR); on input, containing the sequence
!    to be transformed, and on output, containing the transformed sequence.
!
!    Input, integer ( kind = 4 ) LENR, the dimension of the R array.
!    LENR must be at least INC*(N-1)+ 1.
!
!    Input, real ( kind = 8 ) WSAVE(LENSAV).  WSAVE's contents must be
!    initialized with a call to COSQ1I before the first call to routine COSQ1F
!    or COSQ1B for a given transform length N.  WSAVE's contents may be
!    re-used for subsequent calls to COSQ1F and COSQ1B with the same N.
!
!    Input, integer ( kind = 4 ) LENSAV, the dimension of the WSAVE array.
!    LENSAV must be at least 2*N + INT(LOG(REAL(N))) + 4.
!
!    Workspace, real ( kind = 8 ) WORK(LENWRK).
!
!    Input, integer ( kind = 4 ) LENWRK, the dimension of the WORK array.
!    LENWRK must be at least N.
!
!    Output, integer ( kind = 4 ) IER, error flag.
!    0, successful exit;
!    1, input parameter LENR   not big enough;
!    2, input parameter LENSAV not big enough;
!    3, input parameter LENWRK not big enough;
!    20, input error returned by lower level routine.
!
  implicit none

  integer ( kind = 4 ) inc
  integer ( kind = 4 ) lensav
  integer ( kind = 4 ) lenwrk

  integer ( kind = 4 ) ier
  integer ( kind = 4 ) ier1
  integer ( kind = 4 ) lenx
  integer ( kind = 4 ) n
  real ( kind = 8 ) ssqrt2
  real ( kind = 8 ) work(lenwrk)
  real ( kind = 8 ) wsave(lensav)
  real ( kind = 8 ) x(inc,*)
  real ( kind = 8 ) x1

  ier = 0

  if (lenx < inc*(n-1) + 1) then
    ier = 1
    call xerfft ('cosq1b', 6)
    return
  else if (lensav < &
    2*n + int(log( real ( n, kind = 8 ) )/log( 2.0D+00 )) +4) then
    ier = 2
    call xerfft ('cosq1b', 8)
    return
  else if (lenwrk < n) then
    ier = 3
    call xerfft ('cosq1b', 10)
    return
  end if

      if (n-2) 300,102,103

 102  ssqrt2 = 1.0D+00 / sqrt ( 2.0D+00 )
      x1 = x(1,1)+x(1,2)
      x(1,2) = ssqrt2*(x(1,1)-x(1,2))
      x(1,1) = x1
      return

  103 call cosqb1 (n,inc,x,wsave,work,ier1)

      if (ier1 /= 0) then
        ier = 20
        call xerfft ('cosq1b',-5)
      end if

  300 continue

  return
end
