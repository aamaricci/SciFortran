subroutine rwupdt ( n, r, ldr, w, b, alpha, c, s )

!*****************************************************************************80
!
!! RWUPDT computes the decomposition of a triangular matrix augmented by one row.
!
!  Discussion:
!
!    Given an N by N upper triangular matrix R, this subroutine
!    computes the QR decomposition of the matrix formed when a row
!    is added to R.  If the row is specified by the vector W, then
!    RWUPDT determines an orthogonal matrix Q such that when the
!    N+1 by N matrix composed of R augmented by W is premultiplied
!    by Q', the resulting matrix is upper trapezoidal.
!    The matrix Q' is the product of N transformations
!
!      G(N)*G(N-1)* ... *G(1)
!
!    where G(I) is a Givens rotation in the (I,N+1) plane which eliminates
!    elements in the (N+1)-st plane.  RWUPDT also computes the product
!    Q'*C where C is the (N+1)-vector (B,ALPHA).  Q itself is not
!    accumulated, rather the information to recover the G rotations is
!    supplied.
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
!    Input, integer ( kind = 4 ) N, the order of R.
!
!    Input/output, real ( kind = 8 ) R(LDR,N).  On input the upper triangular
!    part of R must contain the matrix to be updated.  On output R contains the
!    updated triangular matrix.
!
!    Input, integer ( kind = 4 ) LDR, the leading dimension of the array R.
!    LDR must not be less than N.
!
!    Input, real ( kind = 8 ) W(N), the row vector to be added to R.
!
!    Input/output, real ( kind = 8 ) B(N).  On input, the first N elements
!    of the vector C.  On output the first N elements of the vector Q'*C.
!
!    Input/output, real ( kind = 8 ) ALPHA.  On input, the (N+1)-st element
!    of the vector C.  On output the (N+1)-st element of the vector Q'*C.
!
!    Output, real ( kind = 8 ) C(N), S(N), the cosines and sines of the
!    transforming Givens rotations.
!
  implicit none

  integer ( kind = 4 ) ldr
  integer ( kind = 4 ) n

  real    ( kind = 8 ) alpha
  real    ( kind = 8 ) b(n)
  real    ( kind = 8 ) c(n)
  real    ( kind = 8 ) cotan
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real    ( kind = 8 ) r(ldr,n)
  real    ( kind = 8 ) rowj
  real    ( kind = 8 ) s(n)
  real    ( kind = 8 ) tan
  real    ( kind = 8 ) temp
  real    ( kind = 8 ) w(n)

  do j = 1, n

    rowj = w(j)
!
!  Apply the previous transformations to R(I,J), I=1,2,...,J-1, and to W(J).
!
    do i = 1, j-1
      temp =   c(i) * r(i,j) + s(i) * rowj
      rowj = - s(i) * r(i,j) + c(i) * rowj
      r(i,j) = temp
    end do
!
!  Determine a Givens rotation which eliminates W(J).
!
    c(j) = 1.0D+00
    s(j) = 0.0D+00

    if ( rowj /= 0.0D+00 ) then

      if ( abs ( r(j,j) ) < abs ( rowj ) ) then
        cotan = r(j,j) / rowj
        s(j) = 0.5D+00 / sqrt ( 0.25D+00 + 0.25D+00 * cotan**2 )
        c(j) = s(j) * cotan
      else
        tan = rowj / r(j,j)
        c(j) = 0.5D+00 / sqrt ( 0.25D+00 + 0.25D+00 * tan**2 )
        s(j) = c(j) * tan
      end if
!
!  Apply the current transformation to R(J,J), B(J), and ALPHA.
!
      r(j,j) =  c(j) * r(j,j) + s(j) * rowj
      temp =    c(j) * b(j)   + s(j) * alpha
      alpha = - s(j) * b(j)   + c(j) * alpha
      b(j) = temp

    end if

  end do

  return
end
