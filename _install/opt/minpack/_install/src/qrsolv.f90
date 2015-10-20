subroutine qrsolv ( n, r, ldr, ipvt, diag, qtb, x, sdiag )

!*****************************************************************************80
!
!! QRSOLV solves a rectangular linear system A*x=b in the least squares sense.
!
!  Discussion:
!
!    Given an M by N matrix A, an N by N diagonal matrix D,
!    and an M-vector B, the problem is to determine an X which
!    solves the system
!
!      A*X = B
!      D*X = 0
!
!    in the least squares sense.
!
!    This subroutine completes the solution of the problem
!    if it is provided with the necessary information from the
!    QR factorization, with column pivoting, of A.  That is, if
!    Q*P = Q*R, where P is a permutation matrix, Q has orthogonal
!    columns, and R is an upper triangular matrix with diagonal
!    elements of nonincreasing magnitude, then QRSOLV expects
!    the full upper triangle of R, the permutation matrix p,
!    and the first N components of Q'*B.
!
!    The system is then equivalent to
!
!      R*Z = Q'*B
!      P'*D*P*Z = 0
!
!    where X = P*Z.  If this system does not have full rank,
!    then a least squares solution is obtained.  On output QRSOLV
!    also provides an upper triangular matrix S such that
!
!      P'*(A'*A + D*D)*P = S'*S.
!
!    S is computed within QRSOLV and may be of separate interest.
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
!    Input/output, real ( kind = 8 ) R(LDR,N), the N by N matrix.
!    On input the full upper triangle must contain the full upper triangle
!    of the matrix R.  On output the full upper triangle is unaltered, and
!    the strict lower triangle contains the strict upper triangle
!    (transposed) of the upper triangular matrix S.
!
!    Input, integer ( kind = 4 ) LDR, the leading dimension of R, which must be
!    at least N.
!
!    Input, integer ( kind = 4 ) IPVT(N), defines the permutation matrix P such that
!    A*P = Q*R.  Column J of P is column IPVT(J) of the identity matrix.
!
!    Input, real ( kind = 8 ) DIAG(N), the diagonal elements of the matrix D.
!
!    Input, real ( kind = 8 ) QTB(N), the first N elements of the vector Q'*B.
!
!    Output, real ( kind = 8 ) X(N), the least squares solution.
!
!    Output, real ( kind = 8 ) SDIAG(N), the diagonal elements of the upper
!    triangular matrix S.
!
  implicit none

  integer ( kind = 4 ) ldr
  integer ( kind = 4 ) n

  real    ( kind = 8 ) c
  real    ( kind = 8 ) cotan
  real    ( kind = 8 ) diag(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ipvt(n)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) l
  integer ( kind = 4 ) nsing
  real    ( kind = 8 ) qtb(n)
  real    ( kind = 8 ) qtbpj
  real    ( kind = 8 ) r(ldr,n)
  real    ( kind = 8 ) s
  real    ( kind = 8 ) sdiag(n)
  real    ( kind = 8 ) sum2
  real    ( kind = 8 ) t
  real    ( kind = 8 ) temp
  real    ( kind = 8 ) wa(n)
  real    ( kind = 8 ) x(n)
!
!  Copy R and Q'*B to preserve input and initialize S.
!
!  In particular, save the diagonal elements of R in X.
!
  do j = 1, n
    r(j:n,j) = r(j,j:n)
    x(j) = r(j,j)
  end do

  wa(1:n) = qtb(1:n)
!
!  Eliminate the diagonal matrix D using a Givens rotation.
!
  do j = 1, n
!
!  Prepare the row of D to be eliminated, locating the
!  diagonal element using P from the QR factorization.
!
    l = ipvt(j)

    if ( diag(l) /= 0.0D+00 ) then

      sdiag(j:n) = 0.0D+00
      sdiag(j) = diag(l)
!
!  The transformations to eliminate the row of D
!  modify only a single element of Q'*B
!  beyond the first N, which is initially zero.
!
      qtbpj = 0.0D+00

      do k = j, n
!
!  Determine a Givens rotation which eliminates the
!  appropriate element in the current row of D.
!
        if ( sdiag(k) /= 0.0D+00 ) then

          if ( abs ( r(k,k) ) < abs ( sdiag(k) ) ) then
            cotan = r(k,k) / sdiag(k)
            s = 0.5D+00 / sqrt ( 0.25D+00 + 0.25D+00 * cotan**2 )
            c = s * cotan
          else
            t = sdiag(k) / r(k,k)
            c = 0.5D+00 / sqrt ( 0.25D+00 + 0.25D+00 * t**2 )
            s = c * t
          end if
!
!  Compute the modified diagonal element of R and
!  the modified element of (Q'*B,0).
!
          r(k,k) = c * r(k,k) + s * sdiag(k)
          temp = c * wa(k) + s * qtbpj
          qtbpj = - s * wa(k) + c * qtbpj
          wa(k) = temp
!
!  Accumulate the tranformation in the row of S.
!
          do i = k+1, n
            temp = c * r(i,k) + s * sdiag(i)
            sdiag(i) = - s * r(i,k) + c * sdiag(i)
            r(i,k) = temp
          end do

        end if

      end do

    end if
!
!  Store the diagonal element of S and restore
!  the corresponding diagonal element of R.
!
    sdiag(j) = r(j,j)
    r(j,j) = x(j)

  end do
!
!  Solve the triangular system for Z.  If the system is
!  singular, then obtain a least squares solution.
!
  nsing = n

  do j = 1, n

    if ( sdiag(j) == 0.0D+00 .and. nsing == n ) then
      nsing = j - 1
    end if

    if ( nsing < n ) then
      wa(j) = 0.0D+00
    end if

  end do

  do j = nsing, 1, -1
    sum2 = dot_product ( wa(j+1:nsing), r(j+1:nsing,j) )
    wa(j) = ( wa(j) - sum2 ) / sdiag(j)
  end do
!
!  Permute the components of Z back to components of X.
!
  do j = 1, n
    l = ipvt(j)
    x(l) = wa(j)
  end do

  return
end
