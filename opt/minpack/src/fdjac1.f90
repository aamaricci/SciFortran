subroutine fdjac1 ( fcn, n, x, fvec, fjac, ldfjac, iflag, ml, mu, epsfcn )

!*****************************************************************************80
!
!! FDJAC1 estimates an N by N jacobian matrix using forward differences.
!
!  Discussion:
!
!    This subroutine computes a forward-difference approximation
!    to the N by N jacobian matrix associated with a specified
!    problem of N functions in N variables. If the jacobian has
!    a banded form, then function evaluations are saved by only
!    approximating the nonzero terms.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 April 2010
!
!  Author:
!
!    Original FORTRAN77 version by Jorge More, Burton Garbow, Kenneth Hillstrom.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Jorge More, Burton Garbow, Kenneth Hillstrom,
!    User Guide for MINPACK-1,
!    Technical Report ANL-80-74,
!    Argonne National Laboratory, 1980.
!
!  Parameters:
!
!    Input, external FCN, the name of the user-supplied subroutine which
!    calculates the functions.  The routine should have the form:
!
!      subroutine fcn ( n, x, fvec, iflag )
!
!      integer ( kind = 4 ) n
!
!      real fvec(n)
!      integer ( kind = 4 ) iflag
!      real x(n)
!
!    The value of IFLAG should not be changed by FCN unless
!    the user wants to terminate execution of the routine.
!    In this case set IFLAG to a negative integer.
!
!    Input, integer ( kind = 4 ) N, the number of functions and variables.
!
!    Input, real ( kind = 8 ) X(N), the point where the jacobian is evaluated.
!
!    Input, real ( kind = 8 ) FVEC(N), the functions evaluated at X.
!
!    Output, real ( kind = 8 ) FJAC(LDFJAC,N), the N by N approximate
!    jacobian matrix.
!
!    Input, integer ( kind = 4 ) LDFJAC, the leading dimension of FJAC, which must
!    not be less than N.
!
!    Output, integer ( kind = 4 ) IFLAG, is an error flag returned by FCN.  If FCN
!    returns a nonzero value of IFLAG, then this routine returns immediately
!    to the calling program, with the value of IFLAG.
!
!    Input, integer ( kind = 4 ) ML, MU, specify the number of subdiagonals and
!    superdiagonals within the band of the jacobian matrix.  If the
!    jacobian is not banded, set ML and MU to N-1.
!
!    Input, real ( kind = 8 ) EPSFCN, is used in determining a suitable step
!    length for the forward-difference approximation.  This approximation
!    assumes that the relative errors in the functions are of the order of
!    EPSFCN.  If EPSFCN is less than the machine precision, it is assumed that
!    the relative errors in the functions are of the order of the machine
!    precision.
!
  implicit none

  integer ( kind = 4 ) ldfjac
  integer ( kind = 4 ) n

  real    ( kind = 8 ) eps
  real    ( kind = 8 ) epsfcn
  real    ( kind = 8 ) epsmch
  external             fcn
  real    ( kind = 8 ) fjac(ldfjac,n)
  real    ( kind = 8 ) fvec(n)
  real    ( kind = 8 ) h
  integer ( kind = 4 ) i
  integer ( kind = 4 ) iflag
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) ml
  integer ( kind = 4 ) msum
  integer ( kind = 4 ) mu
  real    ( kind = 8 ) temp
  real    ( kind = 8 ) wa1(n)
  real    ( kind = 8 ) wa2(n)
  real    ( kind = 8 ) x(n)

  epsmch = epsilon ( epsmch )

  eps = sqrt ( max ( epsfcn, epsmch ) )
  msum = ml + mu + 1
!
!  Computation of dense approximate jacobian.
!
  if ( n <= msum ) then

     do j = 1, n

        temp = x(j)
        h = eps * abs ( temp )
        if ( h == 0.0D+00 ) then
          h = eps
        end if

        x(j) = temp + h
        call fcn ( n, x, wa1, iflag )

        if ( iflag < 0 ) then
          exit
        end if

        x(j) = temp
        fjac(1:n,j) = ( wa1(1:n) - fvec(1:n) ) / h

     end do

  else
!
!  Computation of banded approximate jacobian.
!
     do k = 1, msum

        do j = k, n, msum
          wa2(j) = x(j)
          h = eps * abs ( wa2(j) )
          if ( h == 0.0D+00 ) then
            h = eps
          end if
          x(j) = wa2(j) + h
        end do

        call fcn ( n, x, wa1, iflag )

        if ( iflag < 0 ) then
          exit
        end if

        do j = k, n, msum

           x(j) = wa2(j)

           h = eps * abs ( wa2(j) )
           if ( h == 0.0D+00 ) then
             h = eps
           end if

           fjac(1:n,j) = 0.0D+00

           do i = 1, n
             if ( j - mu <= i .and. i <= j + ml ) then
               fjac(i,j) = ( wa1(i) - fvec(i) ) / h
             end if
           end do

        end do

     end do

  end if

  return
end
