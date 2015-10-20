subroutine mradf3 (m,ido,l1,cc,im1,in1,ch,im2,in2,wa1,wa2)

!*****************************************************************************80
!
!! MRADF3 is an FFTPACK5.1 auxilliary function.
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
  integer ( kind = 4 ) l1

  real ( kind = 8 ) arg
  real ( kind = 8 ) cc(in1,ido,l1,3)
  real ( kind = 8 ) ch(in2,ido,3,l1)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ic
  integer ( kind = 4 ) idp2
  integer ( kind = 4 ) im1
  integer ( kind = 4 ) im2
  integer ( kind = 4 ) k
  integer ( kind = 4 ) m
  integer ( kind = 4 ) m1
  integer ( kind = 4 ) m1d
  integer ( kind = 4 ) m2
  integer ( kind = 4 ) m2s
  real ( kind = 8 ) taui
  real ( kind = 8 ) taur
  real ( kind = 8 ) wa1(ido)
  real ( kind = 8 ) wa2(ido)

  m1d = (m-1)*im1+1
  m2s = 1-im2
  arg= 2.0D+00 * 4.0D+00 * atan ( 1.0D+00 ) / 3.0D+00
  taur=cos(arg)
  taui=sin(arg)

      do 101 k = 1, l1
         m2 = m2s
         do 1001 m1=1,m1d,im1
         m2 = m2+im2
         ch(m2,1,1,k) = cc(m1,1,k,1)+(cc(m1,1,k,2)+cc(m1,1,k,3))
         ch(m2,1,3,k) = taui*(cc(m1,1,k,3)-cc(m1,1,k,2))
         ch(m2,ido,2,k) = cc(m1,1,k,1)+taur*(cc(m1,1,k,2)+cc(m1,1,k,3))
 1001    continue
  101 continue

  if (ido == 1) then
    return
  end if

      idp2 = ido+2
      do 103 k = 1, l1
         do 102 i=3,ido,2
            ic = idp2-i
            m2 = m2s
            do 1002 m1=1,m1d,im1
            m2 = m2+im2
            ch(m2,i-1,1,k) = cc(m1,i-1,k,1)+((wa1(i-2)*cc(m1,i-1,k,2)+ &
             wa1(i-1)*cc(m1,i,k,2))+(wa2(i-2)*cc(m1,i-1,k,3)+wa2(i-1)* &
             cc(m1,i,k,3)))

            ch(m2,i,1,k) = cc(m1,i,k,1)+((wa1(i-2)*cc(m1,i,k,2)- &
             wa1(i-1)*cc(m1,i-1,k,2))+(wa2(i-2)*cc(m1,i,k,3)-wa2(i-1)* &
             cc(m1,i-1,k,3)))

            ch(m2,i-1,3,k) = (cc(m1,i-1,k,1)+taur*((wa1(i-2)* &
             cc(m1,i-1,k,2)+wa1(i-1)*cc(m1,i,k,2))+(wa2(i-2)* &
             cc(m1,i-1,k,3)+wa2(i-1)*cc(m1,i,k,3))))+(taui*((wa1(i-2)* &
             cc(m1,i,k,2)-wa1(i-1)*cc(m1,i-1,k,2))-(wa2(i-2)* &
             cc(m1,i,k,3)-wa2(i-1)*cc(m1,i-1,k,3))))

            ch(m2,ic-1,2,k) = (cc(m1,i-1,k,1)+taur*((wa1(i-2)* &
             cc(m1,i-1,k,2)+wa1(i-1)*cc(m1,i,k,2))+(wa2(i-2)* &
             cc(m1,i-1,k,3)+wa2(i-1)*cc(m1,i,k,3))))-(taui*((wa1(i-2)* &
             cc(m1,i,k,2)-wa1(i-1)*cc(m1,i-1,k,2))-(wa2(i-2)* &
             cc(m1,i,k,3)-wa2(i-1)*cc(m1,i-1,k,3))))

            ch(m2,i,3,k) = (cc(m1,i,k,1)+taur*((wa1(i-2)*cc(m1,i,k,2)- &
             wa1(i-1)*cc(m1,i-1,k,2))+(wa2(i-2)*cc(m1,i,k,3)-wa2(i-1)* &
             cc(m1,i-1,k,3))))+(taui*((wa2(i-2)*cc(m1,i-1,k,3)+wa2(i-1)* &
             cc(m1,i,k,3))-(wa1(i-2)*cc(m1,i-1,k,2)+wa1(i-1)* &
             cc(m1,i,k,2))))

            ch(m2,ic,2,k) = (taui*((wa2(i-2)*cc(m1,i-1,k,3)+wa2(i-1)* &
             cc(m1,i,k,3))-(wa1(i-2)*cc(m1,i-1,k,2)+wa1(i-1)* &
             cc(m1,i,k,2))))-(cc(m1,i,k,1)+taur*((wa1(i-2)*cc(m1,i,k,2)- &
             wa1(i-1)*cc(m1,i-1,k,2))+(wa2(i-2)*cc(m1,i,k,3)-wa2(i-1)* &
             cc(m1,i-1,k,3))))
 1002       continue
  102    continue
  103 continue

  return
end
