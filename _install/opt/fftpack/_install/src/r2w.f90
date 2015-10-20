subroutine r2w ( ldr, ldw, l, m, r, w )

!*****************************************************************************80
!
!! R2W copies a 2D array, allowing for different leading dimensions.
!
!  License:
!
!    Licensed under the GNU General Public License (GPL).
!    Copyright (C) 1995-2004, Scientific Computing Division,
!    University Corporation for Atmospheric Research
!
!  Modified:
!
!    13 May 2013
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
!    Input, integer ( kind = 4 ) LDR, the leading dimension of R.
!
!    Input, integer ( kind = 4 ) LDW, the leading dimension of W.
!
!    Input, integer ( kind = 4 ) L, M, the number of rows and columns
!    of information to be copied.
!
!    Input, real ( kind = 8 ) R(LDR,M), the copied information.
!
!    Output, real (  kind = 8 ) W(LDW,M), the original information.
!
  implicit none

  integer ( kind = 4 ) ldr
  integer ( kind = 4 ) ldw
  integer ( kind = 4 ) m

  integer ( kind = 4 ) l
  real ( kind = 8 ) r(ldr,m)
  real ( kind = 8 ) w(ldw,m)

  w(1:l,1:m) = r(1:l,1:m)

  return
end
