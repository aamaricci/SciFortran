subroutine msntb1(lot,jump,n,inc,x,wsave,dsum,xh,work,ier)

!*****************************************************************************80
!
!! MSNTB1 is an FFTPACK5.1 auxilliary function.
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

  real ( kind = 8 ) dsum(*)
  real ( kind = 8 ) fnp1s4
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ier
  integer ( kind = 4 ) ier1
  integer ( kind = 4 ) jump
  integer ( kind = 4 ) k
  integer ( kind = 4 ) kc
  integer ( kind = 4 ) lj
  integer ( kind = 4 ) lnsv
  integer ( kind = 4 ) lnwk
  integer ( kind = 4 ) lnxh
  integer ( kind = 4 ) m
  integer ( kind = 4 ) m1
  integer ( kind = 4 ) modn
  integer ( kind = 4 ) n
  integer ( kind = 4 ) np1
  integer ( kind = 4 ) ns2
  real ( kind = 8 ) srt3s2
  real ( kind = 8 ) t1
  real ( kind = 8 ) t2
  real ( kind = 8 ) work(*)
  real ( kind = 8 ) wsave(*)
  real ( kind = 8 ) x(inc,*)
  real ( kind = 8 ) xh(lot,*)
  real ( kind = 8 ) xhold

  ier = 0
  lj = (lot-1)*jump+1

      if (n-2) 200,102,103
 102  srt3s2 = sqrt( 3.0D+00 )/ 2.0D+00

      do m=1,lj,jump
         xhold = srt3s2*(x(m,1)+x(m,2))
         x(m,2) = srt3s2*(x(m,1)-x(m,2))
         x(m,1) = xhold
      end do

      go to 200
  103 np1 = n+1
      ns2 = n/2
      do 104 k=1,ns2
         kc = np1-k
         m1 = 0
         do 114 m=1,lj,jump
         m1 = m1+1
         t1 = x(m,k)-x(m,kc)
         t2 = wsave(k)*(x(m,k)+x(m,kc))
         xh(m1,k+1) = t1+t2
         xh(m1,kc+1) = t2-t1
  114    continue
  104 continue
      modn = mod(n,2)
      if (modn == 0) go to 124
      m1 = 0
      do 123 m=1,lj,jump
         m1 = m1+1
         xh(m1,ns2+2) =  4.0D+00 * x(m,ns2+1)
  123 continue
  124 do m=1,lot
         xh(m,1) = 0.0D+00
      end do

      lnxh = lot-1 + lot*(np1-1) + 1
      lnsv = np1 + int(log( real ( np1, kind = 8 ))/log( 2.0D+00 )) + 4
      lnwk = lot*np1

      call rfftmf(lot,1,np1,lot,xh,lnxh,wsave(ns2+1),lnsv,work,lnwk,ier1)

      if (ier1 /= 0) then
        ier = 20
        call xerfft ('msntb1',-5)
        go to 200
      end if

      if(mod(np1,2) /= 0) go to 30
      do m=1,lot
        xh(m,np1) = xh(m,np1)+xh(m,np1)
      end do
 30   fnp1s4 = real ( np1 ) / 4.0D+00
      m1 = 0
      do 125 m=1,lj,jump
         m1 = m1+1
         x(m,1) = fnp1s4*xh(m1,1)
         dsum(m1) = x(m,1)
  125 continue
      do 105 i=3,n,2
         m1 = 0
         do 115 m=1,lj,jump
            m1 = m1+1
            x(m,i-1) = fnp1s4*xh(m1,i)
            dsum(m1) = dsum(m1)+fnp1s4*xh(m1,i-1)
            x(m,i) = dsum(m1)
  115    continue
  105 continue
      if (modn /= 0) go to 200
      m1 = 0
      do 116 m=1,lj,jump
         m1 = m1+1
         x(m,n) = fnp1s4*xh(m1,n+1)
  116 continue

  200 continue

  return
end
