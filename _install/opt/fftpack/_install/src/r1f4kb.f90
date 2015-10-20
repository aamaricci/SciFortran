subroutine r1f4kb (ido,l1,cc,in1,ch,in2,wa1,wa2,wa3)

!*****************************************************************************80
!
!! R1F4KB is an FFTPACK5.1 auxilliary function.
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

  real ( kind = 8 ) cc(in1,ido,4,l1)
  real ( kind = 8 ) ch(in2,ido,l1,4)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ic
  integer ( kind = 4 ) idp2
  integer ( kind = 4 ) k
  real ( kind = 8 ) sqrt2
  real ( kind = 8 ) wa1(ido)
  real ( kind = 8 ) wa2(ido)
  real ( kind = 8 ) wa3(ido)

  sqrt2=sqrt( 2.0D+00 )

  do k = 1, l1
    ch(1,1,k,3) = (cc(1,1,1,k)+cc(1,ido,4,k)) &
      -(cc(1,ido,2,k)+cc(1,ido,2,k))
    ch(1,1,k,1) = (cc(1,1,1,k)+cc(1,ido,4,k)) &
      +(cc(1,ido,2,k)+cc(1,ido,2,k))
    ch(1,1,k,4) = (cc(1,1,1,k)-cc(1,ido,4,k)) &
      +(cc(1,1,3,k)+cc(1,1,3,k))
    ch(1,1,k,2) = (cc(1,1,1,k)-cc(1,ido,4,k)) &
      -(cc(1,1,3,k)+cc(1,1,3,k))
  end do

      if (ido-2) 107,105,102

  102 idp2 = ido+2
      do 104 k = 1, l1
         do 103 i=3,ido,2
            ic = idp2-i
        ch(1,i-1,k,1) = (cc(1,i-1,1,k)+cc(1,ic-1,4,k)) &
        +(cc(1,i-1,3,k)+cc(1,ic-1,2,k))
        ch(1,i,k,1) = (cc(1,i,1,k)-cc(1,ic,4,k)) &
        +(cc(1,i,3,k)-cc(1,ic,2,k))
        ch(1,i-1,k,2)=wa1(i-2)*((cc(1,i-1,1,k)-cc(1,ic-1,4,k)) &
        -(cc(1,i,3,k)+cc(1,ic,2,k)))-wa1(i-1) &
        *((cc(1,i,1,k)+cc(1,ic,4,k))+(cc(1,i-1,3,k)-cc(1,ic-1,2,k)))
        ch(1,i,k,2)=wa1(i-2)*((cc(1,i,1,k)+cc(1,ic,4,k)) &
        +(cc(1,i-1,3,k)-cc(1,ic-1,2,k)))+wa1(i-1) &
        *((cc(1,i-1,1,k)-cc(1,ic-1,4,k))-(cc(1,i,3,k)+cc(1,ic,2,k)))
        ch(1,i-1,k,3)=wa2(i-2)*((cc(1,i-1,1,k)+cc(1,ic-1,4,k)) &
        -(cc(1,i-1,3,k)+cc(1,ic-1,2,k)))-wa2(i-1) &
        *((cc(1,i,1,k)-cc(1,ic,4,k))-(cc(1,i,3,k)-cc(1,ic,2,k)))
        ch(1,i,k,3)=wa2(i-2)*((cc(1,i,1,k)-cc(1,ic,4,k)) &
        -(cc(1,i,3,k)-cc(1,ic,2,k)))+wa2(i-1) &
        *((cc(1,i-1,1,k)+cc(1,ic-1,4,k))-(cc(1,i-1,3,k) &
        +cc(1,ic-1,2,k)))
        ch(1,i-1,k,4)=wa3(i-2)*((cc(1,i-1,1,k)-cc(1,ic-1,4,k)) &
        +(cc(1,i,3,k)+cc(1,ic,2,k)))-wa3(i-1) &
       *((cc(1,i,1,k)+cc(1,ic,4,k))-(cc(1,i-1,3,k)-cc(1,ic-1,2,k)))
        ch(1,i,k,4)=wa3(i-2)*((cc(1,i,1,k)+cc(1,ic,4,k)) &
        -(cc(1,i-1,3,k)-cc(1,ic-1,2,k)))+wa3(i-1) &
        *((cc(1,i-1,1,k)-cc(1,ic-1,4,k))+(cc(1,i,3,k)+cc(1,ic,2,k)))
  103    continue
  104 continue
      if (mod(ido,2) == 1) return
  105 continue
      do 106 k = 1, l1
         ch(1,ido,k,1) = (cc(1,ido,1,k)+cc(1,ido,3,k)) &
         +(cc(1,ido,1,k)+cc(1,ido,3,k))
         ch(1,ido,k,2) = sqrt2*((cc(1,ido,1,k)-cc(1,ido,3,k)) &
         -(cc(1,1,2,k)+cc(1,1,4,k)))
         ch(1,ido,k,3) = (cc(1,1,4,k)-cc(1,1,2,k)) &
         +(cc(1,1,4,k)-cc(1,1,2,k))
         ch(1,ido,k,4) = -sqrt2*((cc(1,ido,1,k)-cc(1,ido,3,k)) &
         +(cc(1,1,2,k)+cc(1,1,4,k)))
  106 continue
  107 continue

  return
end
