subroutine cmfgkf ( lot, ido, ip, l1, lid, na, cc, cc1, im1, in1, &
  ch, ch1, im2, in2, wa )

!*****************************************************************************80
!
!! CMFGKF is an FFTPACK5.1 auxiliary routine.
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

  integer ( kind = 4 ) ido
  integer ( kind = 4 ) in1
  integer ( kind = 4 ) in2
  integer ( kind = 4 ) ip
  integer ( kind = 4 ) l1
  integer ( kind = 4 ) lid

  real ( kind = 8 ) cc(2,in1,l1,ip,ido)
  real ( kind = 8 ) cc1(2,in1,lid,ip)
  real ( kind = 8 ) ch(2,in2,l1,ido,ip)
  real ( kind = 8 ) ch1(2,in2,lid,ip)
  real ( kind = 8 ) chold1
  real ( kind = 8 ) chold2
  integer ( kind = 4 ) i
  integer ( kind = 4 ) idlj
  integer ( kind = 4 ) im1
  integer ( kind = 4 ) im2
  integer ( kind = 4 ) ipp2
  integer ( kind = 4 ) ipph
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jc
  integer ( kind = 4 ) k
  integer ( kind = 4 ) ki
  integer ( kind = 4 ) l
  integer ( kind = 4 ) lc
  integer ( kind = 4 ) lot
  integer ( kind = 4 ) m1
  integer ( kind = 4 ) m1d
  integer ( kind = 4 ) m2
  integer ( kind = 4 ) m2s
  integer ( kind = 4 ) na
  real ( kind = 8 ) sn
  real ( kind = 8 ) wa(ido,ip-1,2)
  real ( kind = 8 ) wai
  real ( kind = 8 ) war

      m1d = (lot-1)*im1+1
      m2s = 1-im2
      ipp2 = ip+2
      ipph = (ip+1)/2
      do 110 ki=1,lid
         m2 = m2s
         do 110 m1=1,m1d,im1
         m2 = m2+im2
         ch1(1,m2,ki,1) = cc1(1,m1,ki,1)
         ch1(2,m2,ki,1) = cc1(2,m1,ki,1)
  110 continue

      do 111 j=2,ipph
         jc = ipp2-j
         do 112 ki=1,lid
         m2 = m2s
         do 112 m1=1,m1d,im1
         m2 = m2+im2
            ch1(1,m2,ki,j) =  cc1(1,m1,ki,j)+cc1(1,m1,ki,jc)
            ch1(1,m2,ki,jc) = cc1(1,m1,ki,j)-cc1(1,m1,ki,jc)
            ch1(2,m2,ki,j) =  cc1(2,m1,ki,j)+cc1(2,m1,ki,jc)
            ch1(2,m2,ki,jc) = cc1(2,m1,ki,j)-cc1(2,m1,ki,jc)
  112    continue
  111 continue
      do 118 j=2,ipph
         do 117 ki=1,lid
         m2 = m2s
         do 117 m1=1,m1d,im1
         m2 = m2+im2
            cc1(1,m1,ki,1) = cc1(1,m1,ki,1)+ch1(1,m2,ki,j)
            cc1(2,m1,ki,1) = cc1(2,m1,ki,1)+ch1(2,m2,ki,j)
  117    continue
  118 continue
      do 116 l=2,ipph
         lc = ipp2-l
         do 113 ki=1,lid
         m2 = m2s
         do 113 m1=1,m1d,im1
         m2 = m2+im2
            cc1(1,m1,ki,l) = ch1(1,m2,ki,1)+wa(1,l-1,1)*ch1(1,m2,ki,2)
            cc1(1,m1,ki,lc) = -wa(1,l-1,2)*ch1(1,m2,ki,ip)
            cc1(2,m1,ki,l) = ch1(2,m2,ki,1)+wa(1,l-1,1)*ch1(2,m2,ki,2)
            cc1(2,m1,ki,lc) = -wa(1,l-1,2)*ch1(2,m2,ki,ip)
  113    continue
         do 115 j=3,ipph
            jc = ipp2-j
            idlj = mod((l-1)*(j-1),ip)
            war = wa(1,idlj,1)
            wai = -wa(1,idlj,2)
            do 114 ki=1,lid
               m2 = m2s
               do 114 m1=1,m1d,im1
               m2 = m2+im2
               cc1(1,m1,ki,l) = cc1(1,m1,ki,l)+war*ch1(1,m2,ki,j)
               cc1(1,m1,ki,lc) = cc1(1,m1,ki,lc)+wai*ch1(1,m2,ki,jc)
               cc1(2,m1,ki,l) = cc1(2,m1,ki,l)+war*ch1(2,m2,ki,j)
               cc1(2,m1,ki,lc) = cc1(2,m1,ki,lc)+wai*ch1(2,m2,ki,jc)
  114       continue
  115    continue
  116 continue
      if ( 1 < ido ) go to 136
      sn = 1.0D+00 / real ( ip * l1, kind = 8 )
      if (na == 1) go to 146
      do 149 ki=1,lid
         m2 = m2s
         do 149 m1=1,m1d,im1
         m2 = m2+im2
         cc1(1,m1,ki,1) = sn*cc1(1,m1,ki,1)
         cc1(2,m1,ki,1) = sn*cc1(2,m1,ki,1)
  149 continue
      do 120 j=2,ipph
         jc = ipp2-j
         do 119 ki=1,lid
         do 119 m1=1,m1d,im1
            chold1 = sn*(cc1(1,m1,ki,j)-cc1(2,m1,ki,jc))
            chold2 = sn*(cc1(1,m1,ki,j)+cc1(2,m1,ki,jc))
            cc1(1,m1,ki,j) = chold1
            cc1(2,m1,ki,jc) = sn*(cc1(2,m1,ki,j)-cc1(1,m1,ki,jc))
            cc1(2,m1,ki,j) = sn*(cc1(2,m1,ki,j)+cc1(1,m1,ki,jc))
            cc1(1,m1,ki,jc) = chold2
  119    continue
  120 continue

      return

  146 do 147 ki=1,lid
         m2 = m2s
         do 147 m1=1,m1d,im1
         m2 = m2+im2
         ch1(1,m2,ki,1) = sn*cc1(1,m1,ki,1)
         ch1(2,m2,ki,1) = sn*cc1(2,m1,ki,1)
  147 continue
      do 145 j=2,ipph
         jc = ipp2-j
         do 144 ki=1,lid
         m2 = m2s
         do 144 m1=1,m1d,im1
         m2 = m2+im2
            ch1(1,m2,ki,j) = sn*(cc1(1,m1,ki,j)-cc1(2,m1,ki,jc))
            ch1(2,m2,ki,j) = sn*(cc1(2,m1,ki,j)+cc1(1,m1,ki,jc))
            ch1(1,m2,ki,jc) = sn*(cc1(1,m1,ki,j)+cc1(2,m1,ki,jc))
            ch1(2,m2,ki,jc) = sn*(cc1(2,m1,ki,j)-cc1(1,m1,ki,jc))
  144    continue
  145 continue

      return

  136 do 137 ki=1,lid
         m2 = m2s
         do 137 m1=1,m1d,im1
         m2 = m2+im2
         ch1(1,m2,ki,1) = cc1(1,m1,ki,1)
         ch1(2,m2,ki,1) = cc1(2,m1,ki,1)
  137 continue
      do 135 j=2,ipph
         jc = ipp2-j
         do 134 ki=1,lid
         m2 = m2s
         do 134 m1=1,m1d,im1
         m2 = m2+im2
            ch1(1,m2,ki,j) = cc1(1,m1,ki,j)-cc1(2,m1,ki,jc)
            ch1(2,m2,ki,j) = cc1(2,m1,ki,j)+cc1(1,m1,ki,jc)
            ch1(1,m2,ki,jc) = cc1(1,m1,ki,j)+cc1(2,m1,ki,jc)
            ch1(2,m2,ki,jc) = cc1(2,m1,ki,j)-cc1(1,m1,ki,jc)
  134    continue
  135 continue
      do 131 i=1,ido
         do 130 k = 1, l1
         m2 = m2s
         do 130 m1=1,m1d,im1
         m2 = m2+im2
            cc(1,m1,k,1,i) = ch(1,m2,k,i,1)
            cc(2,m1,k,1,i) = ch(2,m2,k,i,1)
  130    continue
  131 continue
      do 123 j=2,ip
         do 122 k = 1, l1
         m2 = m2s
         do 122 m1=1,m1d,im1
         m2 = m2+im2
            cc(1,m1,k,j,1) = ch(1,m2,k,1,j)
            cc(2,m1,k,j,1) = ch(2,m2,k,1,j)
  122    continue
  123 continue
      do 126 j=2,ip
         do 125 i = 2, ido
            do 124 k = 1, l1
               m2 = m2s
               do 124 m1=1,m1d,im1
               m2 = m2+im2
               cc(1,m1,k,j,i) = wa(i,j-1,1)*ch(1,m2,k,i,j) &
                            +wa(i,j-1,2)*ch(2,m2,k,i,j)
               cc(2,m1,k,j,i) = wa(i,j-1,1)*ch(2,m2,k,i,j) &
                            -wa(i,j-1,2)*ch(1,m2,k,i,j)
  124       continue
  125    continue
  126 continue

  return
end
