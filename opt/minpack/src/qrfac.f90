subroutine qrfac ( m, n, a, lda, pivot, ipvt, lipvt, rdiag, acnorm )

!*****************************************************************************80
!
!! QRFAC computes a QR factorization using Householder transformations.
!
!  Discussion:
!
!    This subroutine uses Householder transformations with column
!    pivoting (optional) to compute a QR factorization of the
!    M by N matrix A.  That is, QRFAC determines an orthogonal
!    matrix Q, a permutation matrix P, and an upper trapezoidal
!    matrix R with diagonal elements of nonincreasing magnitude,
!    such that A*P = Q*R.  The Householder transformation for
!    column K, K = 1,2,...,min(M,N), is of the form
!
!      I - ( 1 / U(K) ) * U * U'
!
!    where U has zeros in the first K-1 positions.  The form of
!    this transformation and the method of pivoting first
!    appeared in the corresponding LINPACK subroutine.
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
!    Input, integer ( kind = 4 ) M, the number of rows of A.
!
!    Input, integer ( kind = 4 ) N, the number of columns of A.
!
!    Input/output, real ( kind = 8 ) A(LDA,N), the M by N array.
!    On input, A contains the matrix for which the QR factorization is to
!    be computed.  On output, the strict upper trapezoidal part of A contains
!    the strict upper trapezoidal part of R, and the lower trapezoidal
!    part of A contains a factored form of Q (the non-trivial elements of
!    the U vectors described above).
!
!    Input, integer ( kind = 4 ) LDA, the leading dimension of A, which must
!    be no less than M.
!
!    Input, logical PIVOT, is TRUE if column pivoting is to be carried out.
!
!    Output, integer ( kind = 4 ) IPVT(LIPVT), defines the permutation matrix P such
!    that A*P = Q*R.  Column J of P is column IPVT(J) of the identity matrix.
!    If PIVOT is false, IPVT is not referenced.
!
!    Input, integer ( kind = 4 ) LIPVT, the dimension of IPVT, which should be N if
!    pivoting is used.
!
!    Output, real ( kind = 8 ) RDIAG(N), contains the diagonal elements of R.
!
!    Output, real ( kind = 8 ) ACNORM(N), the norms of the corresponding
!    columns of the input matrix A.  If this information is not needed,
!    then ACNORM can coincide with RDIAG.
!
  implicit none

  integer ( kind = 4 ) lda
  integer ( kind = 4 ) lipvt
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real    ( kind = 8 ) a(lda,n)
  real    ( kind = 8 ) acnorm(n)
  real    ( kind = 8 ) ajnorm
  real    ( kind = 8 ) enorm
  real    ( kind = 8 ) epsmch
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i4_temp
  integer ( kind = 4 ) ipvt(lipvt)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) kmax
  integer ( kind = 4 ) minmn
  logical              pivot
  real    ( kind = 8 ) r8_temp(m)
  real    ( kind = 8 ) rdiag(n)
  real    ( kind = 8 ) temp
  real    ( kind = 8 ) wa(n)

  epsmch = epsilon ( epsmch )
!
!  Compute the initial column norms and initialize several arrays.
!
  do j = 1, n
    acnorm(j) = enorm ( m, a(1:m,j) )
  end do

  rdiag(1:n) = acnorm(1:n)
  wa(1:n) = acnorm(1:n)

  if ( pivot ) then
    do j = 1, n
      ipvt(j) = j
    end do
  end if
!
!  Reduce A to R with Householder transformations.
!
  minmn = min ( m, n )

  do j = 1, minmn
!
!  Bring the column of largest norm into the pivot position.
!
    if ( pivot ) then

      kmax = j

      do k = j, n
        if ( rdiag(k) > rdiag(kmax) ) then
          kmax = k
        end if
      end do

      if ( kmax /= j ) then

        r8_temp(1:m) = a(1:m,j)
        a(1:m,j)     = a(1:m,kmax)
        a(1:m,kmax)  = r8_temp(1:m)

        rdiag(kmax) = rdiag(j)
        wa(kmax) = wa(j)

        i4_temp    = ipvt(j)
        ipvt(j)    = ipvt(kmax)
        ipvt(kmax) = i4_temp

      end if

    end if
!
!  Compute the Householder transformation to reduce the
!  J-th column of A to a multiple of the J-th unit vector.
!
    ajnorm = enorm ( m-j+1, a(j,j) )

    if ( ajnorm /= 0.0D+00 ) then

      if ( a(j,j) < 0.0D+00 ) then
        ajnorm = -ajnorm
      end if

      a(j:m,j) = a(j:m,j) / ajnorm
      a(j,j) = a(j,j) + 1.0D+00
!
!  Apply the transformation to the remaining columns and update the norms.
!
      do k = j+1, n

        temp = dot_product ( a(j:m,j), a(j:m,k) ) / a(j,j)

        a(j:m,k) = a(j:m,k) - temp * a(j:m,j)

        if ( pivot .and. rdiag(k) /= 0.0D+00 ) then

          temp = a(j,k) / rdiag(k)
          rdiag(k) = rdiag(k) * sqrt ( max ( 0.0D+00, 1.0D+00-temp**2 ) )

          if ( 0.05D+00 * ( rdiag(k) / wa(k) )**2 <= epsmch ) then
            rdiag(k) = enorm ( m-j, a(j+1,k) )
            wa(k) = rdiag(k)
          end if

        end if

      end do

    end if

    rdiag(j) = -ajnorm

  end do

  return
end
