subroutine sinqmb ( lot, jump, n, inc, x, lenx, wsave, lensav, &
  work, lenwrk, ier )

!*****************************************************************************80
!
!! SINQMB: real double precision backward sine quarter wave, multiple vectors.
!
!  Discussion:
!
!    SINQMB computes the one-dimensional Fourier transform of multiple
!    sequences within a real array, where each of the sequences is a
!    sine series with odd wave numbers.  This transform is referred to as
!    the backward transform or Fourier synthesis, transforming the
!    sequences from spectral to physical space.
!
!    This transform is normalized since a call to SINQMB followed
!    by a call to SINQMF (or vice-versa) reproduces the original
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
!    Input, integer ( kind = 4 ) LOT, the number of sequences to be transformed
!    within array R.
!
!    Input, integer ( kind = 4 ) JUMP, the increment between the locations, in
!    array R, of the first elements of two consecutive sequences to be
!    transformed.
!
!    Input, integer ( kind = 4 ) N, the length of each sequence to be
!    transformed.  The transform is most efficient when N is a product of
!    small primes.
!
!    Input, integer ( kind = 4 ) INC, the increment between the locations, in
!    array R, of two consecutive elements within the same sequence.
!
!    Input/output, real ( kind = 8 ) R(LENR), containing LOT sequences, each
!    having length N.  R can have any number of dimensions, but the total
!    number of locations must be at least LENR.  On input, R contains the data
!    to be transformed, and on output the transformed data.
!
!    Input, integer ( kind = 4 ) LENR, the dimension of the R array.
!    LENR must be at least (LOT-1)*JUMP + INC*(N-1)+ 1.
!
!    Input, real ( kind = 8 ) WSAVE(LENSAV).  WSAVE's contents must be
!    initialized with a call to SINQMI before the first call to routine SINQMF
!    or SINQMB for a given transform length N.  WSAVE's contents may be re-used
!    for subsequent calls to SINQMF and SINQMB with the same N.
!
!    Input, integer ( kind = 4 ) LENSAV, the dimension of the WSAVE array.
!    LENSAV must be at least 2*N + INT(LOG(REAL(N))) + 4.
!
!    Workspace, real ( kind = 8 ) WORK(LENWRK).
!
!    Input, integer ( kind = 4 ) LENWRK, the dimension of the WORK array.
!    LENWRK must be at least LOT*N.
!
!    Output, integer ( kind = 4 ) IER, error flag.
!    0, successful exit;
!    1, input parameter LENR not big enough;
!    2, input parameter LENSAV not big enough;
!    3, input parameter LENWRK not big enough;
!    4, input parameters INC,JUMP,N,LOT are not consistent;
!    20, input error returned by lower level routine.
!
  implicit none

  integer ( kind = 4 ) inc
  integer ( kind = 4 ) lensav
  integer ( kind = 4 ) lenwrk

  integer ( kind = 4 ) ier
  integer ( kind = 4 ) ier1
  integer ( kind = 4 ) jump
  integer ( kind = 4 ) k
  integer ( kind = 4 ) kc
  integer ( kind = 4 ) lenx
  integer ( kind = 4 ) lj
  integer ( kind = 4 ) lot
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
  integer ( kind = 4 ) ns2
  real ( kind = 8 ) work(lenwrk)
  real ( kind = 8 ) wsave(lensav)
  real ( kind = 8 ) x(inc,*)
  logical xercon
  real ( kind = 8 ) xhold

      ier = 0

      if (lenx < (lot-1)*jump + inc*(n-1) + 1) then
        ier = 1
        call xerfft ('sinqmb', 6)
      else if (lensav < 2*n + int(log( real ( n, kind = 8 ) ) &
        /log( 2.0D+00 )) +4) then
        ier = 2
        call xerfft ('sinqmb', 8)
      else if (lenwrk < lot*n) then
        ier = 3
        call xerfft ('sinqmb', 10)
      else if (.not. xercon(inc,jump,n,lot)) then
        ier = 4
        call xerfft ('sinqmb', -1)
      end if

      lj = (lot-1)*jump+1
      if (1 < n ) go to 101

      do m=1,lj,jump
         x(m,1) =  4.0D+00 * x(m,1)
      end do

      return
  101 ns2 = n/2

    do k=2,n,2
       do m=1,lj,jump
         x(m,k) = -x(m,k)
      end do
    end do

      call cosqmb (lot,jump,n,inc,x,lenx,wsave,lensav,work,lenwrk,ier1)
      if (ier1 /= 0) then
        ier = 20
        call xerfft ('sinqmb',-5)
        go to 300
      end if
      do 103 k=1,ns2
         kc = n-k
         do 203 m=1,lj,jump
         xhold = x(m,k)
         x(m,k) = x(m,kc+1)
         x(m,kc+1) = xhold
 203     continue
  103 continue
  300 continue

  return
end
