subroutine xerfft ( srname, info )

!*****************************************************************************80
!
!! XERFFT is an error handler for the FFTPACK routines.
!
!  Discussion:
!
!    XERFFT is an error handler for FFTPACK version 5.1 routines.
!    It is called by an FFTPACK 5.1 routine if an input parameter has an
!    invalid value.  A message is printed and execution stops.
!
!    Installers may consider modifying the stop statement in order to
!    call system-specific exception-handling facilities.
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
!    Input, character ( len = * ) SRNAME, the name of the calling routine.
!
!    Input, integer ( kind = 4 ) INFO, an error code.  When a single invalid
!    parameter in the parameter list of the calling routine has been detected,
!    INFO is the position of that parameter.  In the case when an illegal
!    combination of LOT, JUMP, N, and INC has been detected, the calling
!    subprogram calls XERFFT with INFO = -1.
!
  implicit none

  integer ( kind = 4 ) info
  character ( len = * ) srname

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'XERFFT - Fatal error!'

  if ( 1 <= info ) then
    write ( *, '(a,a,a,i3,a)' ) '  On entry to ', trim ( srname ), &
      ' parameter number ', info, ' had an illegal value.'
  else if ( info == -1 ) then
    write ( *, '(a,a,a,a)' ) '  On entry to ', trim ( srname ), &
      ' parameters LOT, JUMP, N and INC are inconsistent.'
  else if ( info == -2 ) then
    write ( *, '(a,a,a,a)' ) '  On entry to ', trim ( srname ), &
      ' parameter L is greater than LDIM.'
  else if ( info == -3 ) then
    write ( *, '(a,a,a,a)' ) '  On entry to ', trim ( srname ), &
      ' parameter M is greater than MDIM.'
  else if ( info == -5 ) then
    write ( *, '(a,a,a,a)' ) '  Within ', trim ( srname ), &
      ' input error returned by lower level routine.'
  else if ( info == -6 ) then
    write ( *, '(a,a,a,a)' ) '  On entry to ', trim ( srname ), &
      ' parameter LDIM is less than 2*(L/2+1).'
  end if

  stop
end
