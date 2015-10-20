subroutine mcsqf1 (lot,jump,n,inc,x,wsave,work,ier)

!*****************************************************************************80
!
!! MCSQF1 is an FFTPACK5.1 auxilliary function.
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
  implicit none

  integer ( kind = 4 ) inc
  integer ( kind = 4 ) lot

  integer ( kind = 4 ) i
  integer ( kind = 4 ) ier
  integer ( kind = 4 ) ier1
  integer ( kind = 4 ) jump
  integer ( kind = 4 ) k
  integer ( kind = 4 ) kc
  integer ( kind = 4 ) lenx
  integer ( kind = 4 ) lnsv
  integer ( kind = 4 ) lnwk
  integer ( kind = 4 ) lj
  integer ( kind = 4 ) m
  integer ( kind = 4 ) m1
  integer ( kind = 4 ) modn
  integer ( kind = 4 ) n
  integer ( kind = 4 ) np2
  integer ( kind = 4 ) ns2
  real ( kind = 8 ) work(lot,*)
  real ( kind = 8 ) wsave(*)
  real ( kind = 8 ) x(inc,*)
  real ( kind = 8 ) xim1

  ier = 0

  lj = (lot-1)*jump+1
  ns2 = (n+1)/2
  np2 = n+2

  do k=2,ns2
    kc = np2-k
    m1 = 0
    do m=1,lj,jump
      m1 = m1 + 1
      work(m1,k)  = x(m,k)+x(m,kc)
      work(m1,kc) = x(m,k)-x(m,kc)
    end do
  end do

  modn = mod(n,2)

  if (modn == 0) then
    m1 = 0
    do m=1,lj,jump
      m1 = m1 + 1
      work(m1,ns2+1) = x(m,ns2+1)+x(m,ns2+1)
    end do
  end if

       do 102 k=2,ns2
         kc = np2-k
         m1 = 0
         do 302 m=1,lj,jump
         m1 = m1 + 1
         x(m,k)  = wsave(k-1)*work(m1,kc)+wsave(kc-1)*work(m1,k)
         x(m,kc) = wsave(k-1)*work(m1,k) -wsave(kc-1)*work(m1,kc)
 302     continue
  102 continue

      if (modn /= 0) go to 303
      m1 = 0
      do 304 m=1,lj,jump
         m1 = m1 + 1
         x(m,ns2+1) = wsave(ns2)*work(m1,ns2+1)
 304  continue
 303  continue

      lenx = (lot-1)*jump + inc*(n-1)  + 1
      lnsv = n + int(log( real ( n, kind = 8 ) )/log( 2.0D+00 )) + 4
      lnwk = lot*n

      call rfftmf(lot,jump,n,inc,x,lenx,wsave(n+1),lnsv,work,lnwk,ier1)
      if (ier1 /= 0) then
        ier = 20
        call xerfft ('mcsqf1',-5)
        go to 400
      end if

      do 103 i=3,n,2
         do 203 m=1,lj,jump
            xim1 = 0.5D+00 * (x(m,i-1)+x(m,i))
            x(m,i) = 0.5D+00 * (x(m,i-1)-x(m,i))
            x(m,i-1) = xim1
 203     continue
  103 continue

  400 continue

  return
end
