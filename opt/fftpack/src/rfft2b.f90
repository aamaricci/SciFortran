subroutine rfft2b ( ldim, l, m, r, wsave, lensav, work, lenwrk, ier )

!*****************************************************************************80
!
!! RFFT2B: real double precision backward fast Fourier transform, 2D.
!
!  Discussion:
!
!    RFFT2B computes the two-dimensional discrete Fourier transform of the
!    complex Fourier coefficients a real periodic array.  This transform is
!    known as the backward transform or Fourier synthesis, transforming from
!    spectral to physical space.  Routine RFFT2B is normalized: a call to
!    RFFT2B followed by a call to RFFT2F (or vice-versa) reproduces the
!    original array within roundoff error.
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
!    Input, integer ( kind = 4 ) LDIM, the first dimension of the 2D real
!    array R, which must be at least 2*(L/2+1).
!
!    Input, integer ( kind = 4 ) L, the number of elements to be transformed
!    in the first dimension of the two-dimensional real array R.  The value of
!    L must be less than or equal to that of LDIM.  The transform is most
!    efficient when L is a product of small primes.
!
!    Input, integer ( kind = 4 ) M, the number of elements to be transformed
!    in the second dimension of the two-dimensional real array R.  The transform
!    is most efficient when M is a product of small primes.
!
!    Input/output, real ( kind = 8 ) R(LDIM,M), the real array of two
!    dimensions.  On input, R contains the L/2+1-by-M complex subarray of
!    spectral coefficients, on output, the physical coefficients.
!
!    Input, real ( kind = 8 ) WSAVE(LENSAV).  WSAVE's contents must be
!    initialized with a call to RFFT2I before the first call to routine RFFT2F
!    or RFFT2B with lengths L and M.  WSAVE's contents may be re-used for
!    subsequent calls to RFFT2F and RFFT2B with the same transform lengths
!    L and M.
!
!    Input, integer ( kind = 4 ) LENSAV, the number of elements in the WSAVE
!    array.  LENSAV must be at least L + M + INT(LOG(REAL(L)))
!    + INT(LOG(REAL(M))) + 8.
!
!    Workspace, real ( kind = 8 ) WORK(LENWRK).  WORK provides workspace, and
!    its contents need not be saved between calls to routines RFFT2B and RFFT2F.
!
!    Input, integer ( kind = 4 )  LENWRK, the number of elements in the WORK
!    array.  LENWRK must be at least LDIM*M.
!
!    Output, integer ( kind = 4 ) IER, the error flag.
!    0, successful exit;
!    2, input parameter LENSAV not big enough;
!    3, input parameter LENWRK not big enough;
!    6, input parameter LDIM < 2*(L/2+1);
!    20, input error returned by lower level routine.
!
  implicit none

  integer ( kind = 4 ) ldim
  integer ( kind = 4 ) lensav
  integer ( kind = 4 ) lenwrk
  integer ( kind = 4 ) m

  integer ( kind = 4 ) i
  integer ( kind = 4 ) ier
  integer ( kind = 4 ) ier1
  integer ( kind = 4 ) j
  integer ( kind = 4 ) l
  integer ( kind = 4 ) ldh
  integer ( kind = 4 ) ldw
  integer ( kind = 4 ) ldx
  integer ( kind = 4 ) lwsav
  integer ( kind = 4 ) mmsav
  integer ( kind = 4 ) modl
  integer ( kind = 4 ) modm
  integer ( kind = 4 ) mwsav
  real ( kind = 8 ) work(lenwrk)
  real ( kind = 8 ) wsave(lensav)
  real ( kind = 8 ) r(ldim,m)

  ier = 0
!
!  verify lensav
!
  lwsav =   l+int(log( real ( l, kind = 8 ) )/log( 2.0D+00 ))+4
  mwsav =   2*m+int(log( real ( m, kind = 8 ) )/log( 2.0D+00 ))+4
  mmsav =   m+int(log( real ( m, kind = 8 ) )/log( 2.0D+00 ))+4
  modl = mod(l,2)
  modm = mod(m,2)

  if (lensav < lwsav+mwsav+mmsav) then
    ier = 2
    call xerfft ('rfft2f', 6)
    return
  end if
!
! verify lenwrk
!
  if (lenwrk < (l+1)*m) then
    ier = 3
    call xerfft ('rfft2f', 8)
    return
  end if
!
! verify ldim is as big as l
!
  if (ldim < l) then
    ier = 5
    call xerfft ('rfft2f', -6)
    return
  end if
!
! transform second dimension of array
!
      do j=2,2*((m+1)/2)-1
      r(1,j) = r(1,j)+r(1,j)
      end do
      do j=3,m,2
      r(1,j) = -r(1,j)
      end do
      call rfftmb(1,1,m,ldim,r,m*ldim, &
           wsave(lwsav+mwsav+1),mmsav,work,lenwrk,ier1)
      ldh = int((l+1)/2)
      if( 1 < ldh ) then
      ldw = ldh+ldh
!
!  r and work are switched because the the first dimension
!  of the input to complex cfftmf must be even.
!
      call r2w(ldim,ldw,l,m,r,work)
      call cfftmb(ldh-1,1,m,ldh,work(2),ldh*m, &
           wsave(lwsav+1),mwsav,r,l*m, ier1)

      if(ier1/=0) then
         ier=20
         call xerfft('rfft2b',-5)
         return
      end if

      call w2r(ldim,ldw,l,m,r,work)
      end if

      if(modl == 0) then
      do j=2,2*((m+1)/2)-1
      r(l,j) = r(l,j)+r(l,j)
      end do
      do j=3,m,2
      r(l,j) = -r(l,j)
      end do
      call rfftmb(1,1,m,ldim,r(l,1),m*ldim, &
           wsave(lwsav+mwsav+1),mmsav,work,lenwrk,ier1)
      end if
!
!  transform first dimension of array
!
      ldx = 2*int((l+1)/2)-1
      do i=2,ldx
      do j=1,m
      r(i,j) = r(i,j)+r(i,j)
      end do
      end do
      do j=1,m
      do i=3,ldx,2
      r(i,j) = -r(i,j)
      end do
      end do
      call rfftmb(m,ldim,l,1,r,m*ldim,wsave(1), &
           l+int(log( real ( l, kind = 8 ) )/log( 2.0D+00 ))+4,work,lenwrk,ier1)

      if(ier1/=0) then
         ier=20
         call xerfft('rfft2f',-5)
         return
      end if

      if(ier1/=0) then
         ier=20
         call xerfft('rfft2f',-5)
         return
      end if

  100 continue

  return
end
