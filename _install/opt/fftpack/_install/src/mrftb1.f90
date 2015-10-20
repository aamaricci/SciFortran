subroutine mrftb1 (m,im,n,in,c,ch,wa,fac)

!*****************************************************************************80
!
!! MRFTB1 is an FFTPACK5.1 auxilliary function.
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
  real ( kind = 8 ) half
  real ( kind = 8 ) halfm
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
  integer ( kind = 4 ) l1
  integer ( kind = 4 ) l2
  integer ( kind = 4 ) m2
  integer ( kind = 4 ) modn
  integer ( kind = 4 ) na
  integer ( kind = 4 ) nf
  integer ( kind = 4 ) nl
  real ( kind = 8 ) wa(n)

  nf = fac(2)
  na = 0

  do k1=1,nf
      ip = fac(k1+2)
      na = 1-na
      if(ip <= 5) go to 10
      if(k1 == nf) go to 10
      na = 1-na
   10 continue
  end do

      half = 0.5D+00
      halfm = -0.5D+00
      modn = mod(n,2)
      nl = n-2
      if(modn /= 0) nl = n-1
      if (na == 0) go to 120
      m2 = 1-im
      do 117 i=1,m
      m2 = m2+im
      ch(i,1) = c(m2,1)
      ch(i,n) = c(m2,n)
  117 continue
      do 118 j=2,nl,2
      m2 = 1-im
      do 118 i=1,m
         m2 = m2+im
	 ch(i,j) = half*c(m2,j)
	 ch(i,j+1) = halfm*c(m2,j+1)
  118 continue
      go to 124
  120 continue
      do 122 j=2,nl,2
      m2 = 1-im
      do 122 i=1,m
         m2 = m2+im
	 c(m2,j) = half*c(m2,j)
	 c(m2,j+1) = halfm*c(m2,j+1)
  122 continue
  124 l1 = 1
      iw = 1
      do 116 k1=1,nf
         ip = fac(k1+2)
         l2 = ip*l1
         ido = n/l2
         idl1 = ido*l1
         if (ip /= 4) go to 103
         ix2 = iw+ido
         ix3 = ix2+ido
         if (na /= 0) go to 101
         call mradb4 (m,ido,l1,c,im,in,ch,1,m,wa(iw),wa(ix2),wa(ix3))
         go to 102
  101    call mradb4 (m,ido,l1,ch,1,m,c,im,in,wa(iw),wa(ix2),wa(ix3))
  102    na = 1-na
         go to 115
  103    if (ip /= 2) go to 106
         if (na /= 0) go to 104
	 call mradb2 (m,ido,l1,c,im,in,ch,1,m,wa(iw))
         go to 105
  104    call mradb2 (m,ido,l1,ch,1,m,c,im,in,wa(iw))
  105    na = 1-na
         go to 115
  106    if (ip /= 3) go to 109
         ix2 = iw+ido
         if (na /= 0) go to 107
	 call mradb3 (m,ido,l1,c,im,in,ch,1,m,wa(iw),wa(ix2))
         go to 108
  107    call mradb3 (m,ido,l1,ch,1,m,c,im,in,wa(iw),wa(ix2))
  108    na = 1-na
         go to 115
  109    if (ip /= 5) go to 112
         ix2 = iw+ido
         ix3 = ix2+ido
         ix4 = ix3+ido
         if (na /= 0) go to 110
         call mradb5 (m,ido,l1,c,im,in,ch,1,m,wa(iw),wa(ix2),wa(ix3),wa(ix4))
         go to 111
  110    call mradb5 (m,ido,l1,ch,1,m,c,im,in,wa(iw),wa(ix2),wa(ix3),wa(ix4))
  111    na = 1-na
         go to 115
  112    if (na /= 0) go to 113
	 call mradbg (m,ido,ip,l1,idl1,c,c,c,im,in,ch,ch,1,m,wa(iw))
         go to 114
  113    call mradbg (m,ido,ip,l1,idl1,ch,ch,ch,1,m,c,c,im,in,wa(iw))
  114    if (ido == 1) na = 1-na
  115    l1 = l2
         iw = iw+(ip-1)*ido
  116 continue

  return
end
