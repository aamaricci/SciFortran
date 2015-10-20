subroutine mradbg (m,ido,ip,l1,idl1,cc,c1,c2,im1,in1,ch,ch2,im2,in2,wa)

!*****************************************************************************80
!
!! MRADBG is an FFTPACK5.1 auxilliary function.
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

  integer ( kind = 4 ) idl1
  integer ( kind = 4 ) ido
  integer ( kind = 4 ) in1
  integer ( kind = 4 ) in2
  integer ( kind = 4 ) ip
  integer ( kind = 4 ) l1

  real ( kind = 8 ) ai1
  real ( kind = 8 ) ai2
  real ( kind = 8 ) ar1
  real ( kind = 8 ) ar1h
  real ( kind = 8 ) ar2
  real ( kind = 8 ) ar2h
  real ( kind = 8 ) arg
  real ( kind = 8 ) c1(in1,ido,l1,ip)
  real ( kind = 8 ) c2(in1,idl1,ip)
  real ( kind = 8 ) cc(in1,ido,ip,l1)
  real ( kind = 8 ) ch(in2,ido,l1,ip)
  real ( kind = 8 ) ch2(in2,idl1,ip)
  real ( kind = 8 ) dc2
  real ( kind = 8 ) dcp
  real ( kind = 8 ) ds2
  real ( kind = 8 ) dsp
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ic
  integer ( kind = 4 ) idij
  integer ( kind = 4 ) idp2
  integer ( kind = 4 ) ik
  integer ( kind = 4 ) im1
  integer ( kind = 4 ) im2
  integer ( kind = 4 ) ipp2
  integer ( kind = 4 ) ipph
  integer ( kind = 4 ) is
  integer ( kind = 4 ) j
  integer ( kind = 4 ) j2
  integer ( kind = 4 ) jc
  integer ( kind = 4 ) k
  integer ( kind = 4 ) l
  integer ( kind = 4 ) lc
  integer ( kind = 4 ) m
  integer ( kind = 4 ) m1
  integer ( kind = 4 ) m1d
  integer ( kind = 4 ) m2
  integer ( kind = 4 ) m2s
  integer ( kind = 4 ) nbd
  real ( kind = 8 ) tpi
  real ( kind = 8 ) wa(ido)

  m1d = ( m - 1 ) * im1 + 1
  m2s = 1 - im2
  tpi = 2.0D+00 * 4.0D+00 * atan ( 1.0D+00 )
  arg = tpi / real ( ip, kind = 8 )
  dcp = cos ( arg )
  dsp = sin ( arg )
  idp2 = ido + 2
  nbd = (ido-1)/2
  ipp2 = ip+2
  ipph = (ip+1)/2

      if (ido < l1) go to 103

      do k = 1, l1
         do i=1,ido
            m2 = m2s
            do m1=1,m1d,im1
            m2 = m2+im2
            ch(m2,i,k,1) = cc(m1,i,1,k)
            end do
      end do
  end do

      go to 106

  103 do 105 i=1,ido
         do 104 k = 1, l1
            m2 = m2s
            do 1004 m1=1,m1d,im1
            m2 = m2+im2
            ch(m2,i,k,1) = cc(m1,i,1,k)
 1004       continue
  104    continue
  105 continue
  106 do 108 j=2,ipph
         jc = ipp2-j
         j2 = j+j
         do 107 k = 1, l1
            m2 = m2s
            do 1007 m1=1,m1d,im1
            m2 = m2+im2
            ch(m2,1,k,j) = cc(m1,ido,j2-2,k)+cc(m1,ido,j2-2,k)
            ch(m2,1,k,jc) = cc(m1,1,j2-1,k)+cc(m1,1,j2-1,k)
 1007       continue
  107    continue
  108 continue
      if (ido == 1) go to 116
      if (nbd < l1) go to 112
      do 111 j=2,ipph
         jc = ipp2-j
         do 110 k = 1, l1
            do 109 i=3,ido,2
               ic = idp2-i
               m2 = m2s
               do 1009 m1=1,m1d,im1
               m2 = m2+im2
               ch(m2,i-1,k,j) = cc(m1,i-1,2*j-1,k)+cc(m1,ic-1,2*j-2,k)
               ch(m2,i-1,k,jc) = cc(m1,i-1,2*j-1,k)-cc(m1,ic-1,2*j-2,k)
               ch(m2,i,k,j) = cc(m1,i,2*j-1,k)-cc(m1,ic,2*j-2,k)
               ch(m2,i,k,jc) = cc(m1,i,2*j-1,k)+cc(m1,ic,2*j-2,k)
 1009          continue
  109       continue
  110    continue
  111 continue
      go to 116
  112 do 115 j=2,ipph
         jc = ipp2-j
         do 114 i=3,ido,2
            ic = idp2-i
            do 113 k = 1, l1
               m2 = m2s
               do 1013 m1=1,m1d,im1
               m2 = m2+im2
               ch(m2,i-1,k,j) = cc(m1,i-1,2*j-1,k)+cc(m1,ic-1,2*j-2,k)
               ch(m2,i-1,k,jc) = cc(m1,i-1,2*j-1,k)-cc(m1,ic-1,2*j-2,k)
               ch(m2,i,k,j) = cc(m1,i,2*j-1,k)-cc(m1,ic,2*j-2,k)
               ch(m2,i,k,jc) = cc(m1,i,2*j-1,k)+cc(m1,ic,2*j-2,k)
 1013          continue
  113       continue
  114    continue
  115 continue
  116 ar1 = 1.0D+00
      ai1 = 0.0D+00
      do 120 l=2,ipph
         lc = ipp2-l
         ar1h = dcp*ar1-dsp*ai1
         ai1 = dcp*ai1+dsp*ar1
         ar1 = ar1h
         do 117 ik=1,idl1
            m2 = m2s
            do 1017 m1=1,m1d,im1
            m2 = m2+im2
            c2(m1,ik,l) = ch2(m2,ik,1)+ar1*ch2(m2,ik,2)
            c2(m1,ik,lc) = ai1*ch2(m2,ik,ip)
 1017       continue
  117    continue
         dc2 = ar1
         ds2 = ai1
         ar2 = ar1
         ai2 = ai1
         do 119 j=3,ipph
            jc = ipp2-j
            ar2h = dc2*ar2-ds2*ai2
            ai2 = dc2*ai2+ds2*ar2
            ar2 = ar2h
            do 118 ik=1,idl1
               m2 = m2s
               do 1018 m1=1,m1d,im1
               m2 = m2+im2
               c2(m1,ik,l) = c2(m1,ik,l)+ar2*ch2(m2,ik,j)
               c2(m1,ik,lc) = c2(m1,ik,lc)+ai2*ch2(m2,ik,jc)
 1018          continue
  118       continue
  119    continue
  120 continue
      do 122 j=2,ipph
         do 121 ik=1,idl1
            m2 = m2s
            do 1021 m1=1,m1d,im1
            m2 = m2+im2
            ch2(m2,ik,1) = ch2(m2,ik,1)+ch2(m2,ik,j)
 1021       continue
  121    continue
  122 continue
      do 124 j=2,ipph
         jc = ipp2-j
         do 123 k = 1, l1
            m2 = m2s
            do 1023 m1=1,m1d,im1
            m2 = m2+im2
            ch(m2,1,k,j) = c1(m1,1,k,j)-c1(m1,1,k,jc)
            ch(m2,1,k,jc) = c1(m1,1,k,j)+c1(m1,1,k,jc)
 1023       continue
  123    continue
  124 continue
      if (ido == 1) go to 132
      if (nbd < l1) go to 128
      do 127 j=2,ipph
         jc = ipp2-j
         do 126 k = 1, l1
            do 125 i=3,ido,2
               m2 = m2s
               do 1025 m1=1,m1d,im1
               m2 = m2+im2
               ch(m2,i-1,k,j) = c1(m1,i-1,k,j)-c1(m1,i,k,jc)
               ch(m2,i-1,k,jc) = c1(m1,i-1,k,j)+c1(m1,i,k,jc)
               ch(m2,i,k,j) = c1(m1,i,k,j)+c1(m1,i-1,k,jc)
               ch(m2,i,k,jc) = c1(m1,i,k,j)-c1(m1,i-1,k,jc)
 1025          continue
  125       continue
  126    continue
  127 continue
      go to 132
  128 do 131 j=2,ipph
         jc = ipp2-j
         do 130 i=3,ido,2
            do 129 k = 1, l1
               m2 = m2s
               do 1029 m1=1,m1d,im1
               m2 = m2+im2
               ch(m2,i-1,k,j) = c1(m1,i-1,k,j)-c1(m1,i,k,jc)
               ch(m2,i-1,k,jc) = c1(m1,i-1,k,j)+c1(m1,i,k,jc)
               ch(m2,i,k,j) = c1(m1,i,k,j)+c1(m1,i-1,k,jc)
               ch(m2,i,k,jc) = c1(m1,i,k,j)-c1(m1,i-1,k,jc)
 1029          continue
  129       continue
  130    continue
  131 continue
  132 continue
      if (ido == 1) return
      do 133 ik=1,idl1
         m2 = m2s
         do 1033 m1=1,m1d,im1
         m2 = m2+im2
         c2(m1,ik,1) = ch2(m2,ik,1)
 1033    continue
  133 continue
      do 135 j=2,ip
         do 134 k = 1, l1
            m2 = m2s
            do 1034 m1=1,m1d,im1
            m2 = m2+im2
            c1(m1,1,k,j) = ch(m2,1,k,j)
 1034       continue
  134    continue
  135 continue
      if (l1 < nbd ) go to 139
      is = -ido
      do 138 j=2,ip
         is = is+ido
         idij = is
         do 137 i=3,ido,2
            idij = idij+2
            do 136 k = 1, l1
               m2 = m2s
               do 1036 m1=1,m1d,im1
               m2 = m2+im2
               c1(m1,i-1,k,j) = wa(idij-1)*ch(m2,i-1,k,j)-wa(idij)* ch(m2,i,k,j)
               c1(m1,i,k,j) = wa(idij-1)*ch(m2,i,k,j)+wa(idij)* ch(m2,i-1,k,j)
 1036          continue
  136       continue
  137    continue
  138 continue
      go to 143
  139 is = -ido
      do 142 j=2,ip
         is = is+ido
         do 141 k = 1, l1
            idij = is
            do 140 i=3,ido,2
               idij = idij+2
               m2 = m2s
               do 1040 m1=1,m1d,im1
               m2 = m2+im2
               c1(m1,i-1,k,j) = wa(idij-1)*ch(m2,i-1,k,j)-wa(idij)*ch(m2,i,k,j)
               c1(m1,i,k,j) = wa(idij-1)*ch(m2,i,k,j)+wa(idij)*ch(m2,i-1,k,j)
 1040          continue
  140       continue
  141    continue
  142 continue
  143 continue

  return
end
