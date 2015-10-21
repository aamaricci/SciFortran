subroutine cmfm1b ( lot, jump, n, inc, c, ch, wa, fnf, fac )

!*****************************************************************************80
!
!! CMFM1B is an FFTPACK5.1 auxiliary routine.
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

  complex ( kind = 8 ) c(*)
  real ( kind = 8 ) ch(*)
  real ( kind = 8 ) fac(*)
  real ( kind = 8 ) fnf
  integer ( kind = 4 ) ido
  integer ( kind = 4 ) inc
  integer ( kind = 4 ) ip
  integer ( kind = 4 ) iw
  integer ( kind = 4 ) jump
  integer ( kind = 4 ) k1
  integer ( kind = 4 ) l1
  integer ( kind = 4 ) l2
  integer ( kind = 4 ) lid
  integer ( kind = 4 ) lot
  integer ( kind = 4 ) n
  integer ( kind = 4 ) na
  integer ( kind = 4 ) nbr
  integer ( kind = 4 ) nf
  real ( kind = 8 ) wa(*)

      nf = int ( fnf )
      na = 0
      l1 = 1
      iw = 1
      do 125 k1=1,nf
         ip = fac(k1)
         l2 = ip*l1
         ido = n/l2
         lid = l1*ido
         nbr = 1+na+2*min(ip-2,4)
         go to (52,62,53,63,54,64,55,65,56,66),nbr
   52    call cmf2kb (lot,ido,l1,na,c,jump,inc,ch,1,lot,wa(iw))
         go to 120
   62    call cmf2kb (lot,ido,l1,na,ch,1,lot,c,jump,inc,wa(iw))
         go to 120
   53    call cmf3kb (lot,ido,l1,na,c,jump,inc,ch,1,lot,wa(iw))
         go to 120
   63    call cmf3kb (lot,ido,l1,na,ch,1,lot,c,jump,inc,wa(iw))
         go to 120
   54    call cmf4kb (lot,ido,l1,na,c,jump,inc,ch,1,lot,wa(iw))
         go to 120
   64    call cmf4kb (lot,ido,l1,na,ch,1,lot,c,jump,inc,wa(iw))
         go to 120
   55    call cmf5kb (lot,ido,l1,na,c,jump,inc,ch,1,lot,wa(iw))
         go to 120
   65    call cmf5kb (lot,ido,l1,na,ch,1,lot,c,jump,inc,wa(iw))
         go to 120
   56    call cmfgkb (lot,ido,ip,l1,lid,na,c,c,jump,inc,ch,ch,1,lot,wa(iw))
         go to 120
   66    call cmfgkb (lot,ido,ip,l1,lid,na,ch,ch,1,lot,c,c, &
           jump,inc,wa(iw))
  120    l1 = l2
         iw = iw+(ip-1)*(ido+ido)
         if(ip <= 5) na = 1-na
  125 continue

  return
end
