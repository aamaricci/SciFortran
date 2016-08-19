subroutine mrftf1 (m,im,n,in,c,ch,wa,fac)

!*****************************************************************************80
!
!! MRFTF1 is an FFTPACK5.1 auxilliary function.
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

  integer ( kind = 4 ) in
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) c(in,*)
  real ( kind = 8 ) ch(m,*)
  real ( kind = 8 ) fac(15)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) idl1
  integer ( kind = 4 ) ido
  integer ( kind = 4 ) im
  integer ( kind = 4 ) ip
  integer ( kind = 4 ) iw
  integer ( kind = 4 ) ix2
  integer ( kind = 4 ) ix3
  integer ( kind = 4 ) ix4
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k1
  integer ( kind = 4 ) kh
  integer ( kind = 4 ) l1
  integer ( kind = 4 ) l2
  integer ( kind = 4 ) m2
  integer ( kind = 4 ) modn
  integer ( kind = 4 ) na
  integer ( kind = 4 ) nf
  integer ( kind = 4 ) nl
  real ( kind = 8 ) sn
  real ( kind = 8 ) tsn
  real ( kind = 8 ) tsnm
  real ( kind = 8 ) wa(n)

  nf = fac(2)
  na = 1
  l2 = n
  iw = n

      do 111 k1=1,nf
         kh = nf-k1
         ip = fac(kh+3)
         l1 = l2/ip
         ido = n/l2
         idl1 = ido*l1
         iw = iw-(ip-1)*ido
         na = 1-na
         if (ip /= 4) go to 102
         ix2 = iw+ido
         ix3 = ix2+ido
         if (na /= 0) go to 101
	     call mradf4 (m,ido,l1,c,im,in,ch,1,m,wa(iw),wa(ix2),wa(ix3))
         go to 110
  101    call mradf4 (m,ido,l1,ch,1,m,c,im,in,wa(iw),wa(ix2),wa(ix3))
         go to 110
  102    if (ip /= 2) go to 104
         if (na /= 0) go to 103
	     call mradf2 (m,ido,l1,c,im,in,ch,1,m,wa(iw))
         go to 110
  103    call mradf2 (m,ido,l1,ch,1,m,c,im,in,wa(iw))
         go to 110
  104    if (ip /= 3) go to 106
         ix2 = iw+ido
         if (na /= 0) go to 105
	     call mradf3 (m,ido,l1,c,im,in,ch,1,m,wa(iw),wa(ix2))
         go to 110
  105    call mradf3 (m,ido,l1,ch,1,m,c,im,in,wa(iw),wa(ix2))
         go to 110
  106    if (ip /= 5) go to 108
         ix2 = iw+ido
         ix3 = ix2+ido
         ix4 = ix3+ido
         if (na /= 0) go to 107
         call mradf5(m,ido,l1,c,im,in,ch,1,m,wa(iw),wa(ix2),wa(ix3),wa(ix4))
         go to 110
  107    call mradf5(m,ido,l1,ch,1,m,c,im,in,wa(iw),wa(ix2),wa(ix3),wa(ix4))
         go to 110
  108    if (ido == 1) na = 1-na
         if (na /= 0) go to 109
         call mradfg (m,ido,ip,l1,idl1,c,c,c,im,in,ch,ch,1,m,wa(iw))
         na = 1
         go to 110
  109    call mradfg (m,ido,ip,l1,idl1,ch,ch,ch,1,m,c,c,im,in,wa(iw))
         na = 0
  110    l2 = l1
  111 continue

      sn = 1.0D+00 / real ( n, kind = 8 )
      tsn =  2.0D+00 / real ( n, kind = 8 )
      tsnm = -tsn
      modn = mod(n,2)
      nl = n-2
      if(modn /= 0) nl = n-1
      if (na /= 0) go to 120
      m2 = 1-im

      do i=1,m
        m2 = m2+im
        c(m2,1) = sn*ch(i,1)
      end do

      do j=2,nl,2
        m2 = 1-im
        do i=1,m
          m2 = m2+im
	      c(m2,j) = tsn*ch(i,j)
	      c(m2,j+1) = tsnm*ch(i,j+1)
        end do
      end do

      if(modn /= 0) return
      m2 = 1-im
      do 119 i=1,m
         m2 = m2+im
         c(m2,n) = sn*ch(i,n)
  119 continue
      return
  120 m2 = 1-im
      do 121 i=1,m
         m2 = m2+im
         c(m2,1) = sn*c(m2,1)
  121 continue
      do 122 j=2,nl,2
      m2 = 1-im
      do 122 i=1,m
         m2 = m2+im
	 c(m2,j) = tsn*c(m2,j)
	 c(m2,j+1) = tsnm*c(m2,j+1)
  122 continue
      if(modn /= 0) return
      m2 = 1-im
      do i=1,m
         m2 = m2+im
         c(m2,n) = sn*c(m2,n)
      end do

  return
end
