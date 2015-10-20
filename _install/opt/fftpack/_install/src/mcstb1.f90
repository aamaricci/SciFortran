subroutine mcstb1(lot,jump,n,inc,x,wsave,dsum,work,ier)

!*****************************************************************************80
!
!! MCSTB1 is an FFTPACK5.1 auxilliary function.
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

  real ( kind = 8 ) dsum(*)
  real ( kind = 8 ) fnm1s2
  real ( kind = 8 ) fnm1s4
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ier
  integer ( kind = 4 ) ier1
  integer ( kind = 4 ) jump
  integer ( kind = 4 ) k
  integer ( kind = 4 ) kc
  integer ( kind = 4 ) lenx
  integer ( kind = 4 ) lj
  integer ( kind = 4 ) lnsv
  integer ( kind = 4 ) lnwk
  integer ( kind = 4 ) lot
  integer ( kind = 4 ) m
  integer ( kind = 4 ) m1
  integer ( kind = 4 ) modn
  integer ( kind = 4 ) n
  integer ( kind = 4 ) nm1
  integer ( kind = 4 ) np1
  integer ( kind = 4 ) ns2
  real ( kind = 8 ) t1
  real ( kind = 8 ) t2
  real ( kind = 8 ) work(*)
  real ( kind = 8 ) wsave(*)
  real ( kind = 8 ) x(inc,*)
  real ( kind = 8 ) x1h
  real ( kind = 8 ) x1p3
  real ( kind = 8 ) x2
  real ( kind = 8 ) xi

  ier = 0

      nm1 = n-1
      np1 = n+1
      ns2 = n/2
      lj = (lot-1)*jump+1
      if (n-2) 106,101,102

  101 do 111 m=1,lj,jump
         x1h = x(m,1)+x(m,2)
         x(m,2) = x(m,1)-x(m,2)
         x(m,1) = x1h
  111 continue

      return

  102 if ( 3 < n ) go to 103
      do 112 m=1,lj,jump
         x1p3 = x(m,1)+x(m,3)
         x2 = x(m,2)
         x(m,2) = x(m,1)-x(m,3)
         x(m,1) = x1p3+x2
         x(m,3) = x1p3-x2
  112 continue

      return

 103  do m=1,lj,jump
        x(m,1) = x(m,1)+x(m,1)
        x(m,n) = x(m,n)+x(m,n)
      end do

      m1 = 0

      do m=1,lj,jump
         m1 = m1+1
         dsum(m1) = x(m,1)-x(m,n)
         x(m,1) = x(m,1)+x(m,n)
      end do

      do 104 k=2,ns2
         m1 = 0
         do 114 m=1,lj,jump
           m1 = m1+1
           kc = np1-k
           t1 = x(m,k)+x(m,kc)
           t2 = x(m,k)-x(m,kc)
           dsum(m1) = dsum(m1)+wsave(kc)*t2
           t2 = wsave(k)*t2
           x(m,k) = t1-t2
           x(m,kc) = t1+t2
  114    continue
  104 continue

      modn = mod(n,2)
      if (modn == 0) go to 124
         do 123 m=1,lj,jump
         x(m,ns2+1) = x(m,ns2+1)+x(m,ns2+1)
  123    continue
 124  continue
      lenx = (lot-1)*jump + inc*(nm1-1)  + 1
      lnsv = nm1 + int(log( real ( nm1, kind = 8 ))/log( 2.0D+00 )) + 4
      lnwk = lot*nm1

      call rfftmf(lot,jump,nm1,inc,x,lenx,wsave(n+1),lnsv,work,lnwk,ier1)

      if (ier1 /= 0) then
        ier = 20
        call xerfft ('mcstb1',-5)
        go to 106
      end if

      fnm1s2 = real ( nm1, kind = 8 ) /  2.0D+00
      m1 = 0
      do 10 m=1,lj,jump
      m1 = m1+1
      dsum(m1) = 0.5D+00 * dsum(m1)
      x(m,1) = fnm1s2 * x(m,1)
   10 continue
      if(mod(nm1,2) /= 0) go to 30
      do 20 m=1,lj,jump
      x(m,nm1) = x(m,nm1)+x(m,nm1)
   20 continue
 30   fnm1s4 = real ( nm1, kind = 8 ) / 4.0D+00
      do 105 i=3,n,2
         m1 = 0
         do 115 m=1,lj,jump
            m1 = m1+1
            xi = fnm1s4*x(m,i)
            x(m,i) = fnm1s4*x(m,i-1)
            x(m,i-1) = dsum(m1)
            dsum(m1) = dsum(m1)+xi
  115 continue
  105 continue
      if (modn /= 0) return
      m1 = 0
      do m=1,lj,jump
         m1 = m1+1
         x(m,n) = dsum(m1)
      end do
  106 continue

  return
end
