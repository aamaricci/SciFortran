subroutine mradfg (m,ido,ip,l1,idl1,cc,c1,c2,im1,in1,ch,ch2,im2,in2,wa)

!*****************************************************************************80
!
!! MRADFG is an FFTPACK5.1 auxilliary function.
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

  m1d = (m-1)*im1+1
  m2s = 1-im2
  tpi= 2.0D+00 * 4.0D+00 * atan ( 1.0D+00 )
  arg = tpi / real ( ip, kind = 8 )
  dcp = cos(arg)
  dsp = sin(arg)
  ipph = (ip+1)/2
  ipp2 = ip+2
  idp2 = ido+2
  nbd = (ido-1)/2

      if (ido == 1) go to 119

      do 101 ik=1,idl1
         m2 = m2s
         do 1001 m1=1,m1d,im1
         m2 = m2+im2
         ch2(m2,ik,1) = c2(m1,ik,1)
 1001    continue
  101 continue

      do 103 j=2,ip
         do 102 k = 1, l1
            m2 = m2s
            do 1002 m1=1,m1d,im1
            m2 = m2+im2
            ch(m2,1,k,j) = c1(m1,1,k,j)
 1002       continue
  102    continue
  103 continue
      if ( l1 < nbd ) go to 107
      is = -ido
      do 106 j=2,ip
         is = is+ido
         idij = is
         do 105 i=3,ido,2
            idij = idij+2
            do 104 k = 1, l1
               m2 = m2s
               do 1004 m1=1,m1d,im1
               m2 = m2+im2
               ch(m2,i-1,k,j) = wa(idij-1)*c1(m1,i-1,k,j)+wa(idij)*c1(m1,i,k,j)
               ch(m2,i,k,j) = wa(idij-1)*c1(m1,i,k,j)-wa(idij)*c1(m1,i-1,k,j)
 1004          continue
  104       continue
  105    continue
  106 continue
      go to 111
  107 is = -ido
      do 110 j=2,ip
         is = is+ido
         do 109 k = 1, l1
            idij = is
            do 108 i=3,ido,2
               idij = idij+2
               m2 = m2s
               do 1008 m1=1,m1d,im1
               m2 = m2+im2
               ch(m2,i-1,k,j) = wa(idij-1)*c1(m1,i-1,k,j)+wa(idij)*c1(m1,i,k,j)
               ch(m2,i,k,j) = wa(idij-1)*c1(m1,i,k,j)-wa(idij)*c1(m1,i-1,k,j)
 1008          continue
  108       continue
  109    continue
  110 continue
  111 if (nbd < l1) go to 115
      do 114 j=2,ipph
         jc = ipp2-j
         do 113 k = 1, l1
            do 112 i=3,ido,2
               m2 = m2s
               do 1012 m1=1,m1d,im1
               m2 = m2+im2
               c1(m1,i-1,k,j) = ch(m2,i-1,k,j)+ch(m2,i-1,k,jc)
               c1(m1,i-1,k,jc) = ch(m2,i,k,j)-ch(m2,i,k,jc)
               c1(m1,i,k,j) = ch(m2,i,k,j)+ch(m2,i,k,jc)
               c1(m1,i,k,jc) = ch(m2,i-1,k,jc)-ch(m2,i-1,k,j)
 1012          continue
  112       continue
  113    continue
  114 continue
      go to 121
  115 do 118 j=2,ipph
         jc = ipp2-j
         do 117 i=3,ido,2
            do 116 k = 1, l1
               m2 = m2s
               do 1016 m1=1,m1d,im1
               m2 = m2+im2
               c1(m1,i-1,k,j) = ch(m2,i-1,k,j)+ch(m2,i-1,k,jc)
               c1(m1,i-1,k,jc) = ch(m2,i,k,j)-ch(m2,i,k,jc)
               c1(m1,i,k,j) = ch(m2,i,k,j)+ch(m2,i,k,jc)
               c1(m1,i,k,jc) = ch(m2,i-1,k,jc)-ch(m2,i-1,k,j)
 1016          continue
  116       continue
  117    continue
  118 continue
      go to 121
  119 do 120 ik=1,idl1
         m2 = m2s
         do 1020 m1=1,m1d,im1
         m2 = m2+im2
         c2(m1,ik,1) = ch2(m2,ik,1)
 1020    continue
  120 continue
  121 do 123 j=2,ipph
         jc = ipp2-j
         do 122 k = 1, l1
            m2 = m2s
            do 1022 m1=1,m1d,im1
            m2 = m2+im2
            c1(m1,1,k,j) = ch(m2,1,k,j)+ch(m2,1,k,jc)
            c1(m1,1,k,jc) = ch(m2,1,k,jc)-ch(m2,1,k,j)
 1022       continue
  122    continue
  123 continue

      ar1 = 1.0D+00
      ai1 = 0.0D+00
      do 127 l=2,ipph
         lc = ipp2-l
         ar1h = dcp*ar1-dsp*ai1
         ai1 = dcp*ai1+dsp*ar1
         ar1 = ar1h
         do 124 ik=1,idl1
            m2 = m2s
            do 1024 m1=1,m1d,im1
            m2 = m2+im2
            ch2(m2,ik,l) = c2(m1,ik,1)+ar1*c2(m1,ik,2)
            ch2(m2,ik,lc) = ai1*c2(m1,ik,ip)
 1024       continue
  124    continue
         dc2 = ar1
         ds2 = ai1
         ar2 = ar1
         ai2 = ai1
         do 126 j=3,ipph
            jc = ipp2-j
            ar2h = dc2*ar2-ds2*ai2
            ai2 = dc2*ai2+ds2*ar2
            ar2 = ar2h
            do 125 ik=1,idl1
               m2 = m2s
               do 1025 m1=1,m1d,im1
               m2 = m2+im2
               ch2(m2,ik,l) = ch2(m2,ik,l)+ar2*c2(m1,ik,j)
               ch2(m2,ik,lc) = ch2(m2,ik,lc)+ai2*c2(m1,ik,jc)
 1025          continue
  125       continue
  126    continue
  127 continue
      do 129 j=2,ipph
         do 128 ik=1,idl1
            m2 = m2s
            do 1028 m1=1,m1d,im1
            m2 = m2+im2
            ch2(m2,ik,1) = ch2(m2,ik,1)+c2(m1,ik,j)
 1028       continue
  128    continue
  129 continue

      if (ido < l1) go to 132
      do 131 k = 1, l1
         do 130 i=1,ido
            m2 = m2s
            do 1030 m1=1,m1d,im1
            m2 = m2+im2
            cc(m1,i,1,k) = ch(m2,i,k,1)
 1030       continue
  130    continue
  131 continue
      go to 135
  132 do 134 i=1,ido
         do 133 k = 1, l1
            m2 = m2s
            do 1033 m1=1,m1d,im1
            m2 = m2+im2
            cc(m1,i,1,k) = ch(m2,i,k,1)
 1033       continue
  133    continue
  134 continue
  135 do 137 j=2,ipph
         jc = ipp2-j
         j2 = j+j
         do 136 k = 1, l1
            m2 = m2s
            do 1036 m1=1,m1d,im1
            m2 = m2+im2
            cc(m1,ido,j2-2,k) = ch(m2,1,k,j)
            cc(m1,1,j2-1,k) = ch(m2,1,k,jc)
 1036       continue
  136    continue
  137 continue
      if (ido == 1) return
      if (nbd < l1) go to 141
      do 140 j=2,ipph
         jc = ipp2-j
         j2 = j+j
         do 139 k = 1, l1
            do 138 i=3,ido,2
               ic = idp2-i
               m2 = m2s
               do 1038 m1=1,m1d,im1
               m2 = m2+im2
               cc(m1,i-1,j2-1,k) = ch(m2,i-1,k,j)+ch(m2,i-1,k,jc)
               cc(m1,ic-1,j2-2,k) = ch(m2,i-1,k,j)-ch(m2,i-1,k,jc)
               cc(m1,i,j2-1,k) = ch(m2,i,k,j)+ch(m2,i,k,jc)
               cc(m1,ic,j2-2,k) = ch(m2,i,k,jc)-ch(m2,i,k,j)
 1038          continue
  138       continue
  139    continue
  140 continue
      return
  141 do 144 j=2,ipph
         jc = ipp2-j
         j2 = j+j
         do 143 i=3,ido,2
            ic = idp2-i
            do 142 k = 1, l1
               m2 = m2s
               do 1042 m1=1,m1d,im1
               m2 = m2+im2
               cc(m1,i-1,j2-1,k) = ch(m2,i-1,k,j)+ch(m2,i-1,k,jc)
               cc(m1,ic-1,j2-2,k) = ch(m2,i-1,k,j)-ch(m2,i-1,k,jc)
               cc(m1,i,j2-1,k) = ch(m2,i,k,j)+ch(m2,i,k,jc)
               cc(m1,ic,j2-2,k) = ch(m2,i,k,jc)-ch(m2,i,k,j)
 1042          continue
  142       continue
  143    continue
  144 continue

  return
end
