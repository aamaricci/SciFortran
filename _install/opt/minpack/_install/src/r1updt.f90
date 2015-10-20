subroutine r1updt ( m, n, s, ls, u, v, w, sing )

!*****************************************************************************80
!
!! R1UPDT re-triangularizes a matrix after a rank one update.
!
!  Discussion:
!
!    Given an M by N lower trapezoidal matrix S, an M-vector U, and an
!    N-vector V, the problem is to determine an orthogonal matrix Q such that
!
!      (S + U * V' ) * Q
!
!    is again lower trapezoidal.
!
!    This subroutine determines Q as the product of 2 * (N - 1)
!    transformations
!
!      GV(N-1)*...*GV(1)*GW(1)*...*GW(N-1)
!
!    where GV(I), GW(I) are Givens rotations in the (I,N) plane
!    which eliminate elements in the I-th and N-th planes,
!    respectively.  Q itself is not accumulated, rather the
!    information to recover the GV and GW rotations is returned.
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
!    Input, integer ( kind = 4 ) M, the number of rows of S.
!
!    Input, integer ( kind = 4 ) N, the number of columns of S.  N must not exceed M.
!
!    Input/output, real ( kind = 8 ) S(LS).  On input, the lower trapezoidal
!    matrix S stored by columns.  On output S contains the lower trapezoidal
!    matrix produced as described above.
!
!    Input, integer ( kind = 4 ) LS, the length of the S array.  LS must be at least
!    (N*(2*M-N+1))/2.
!
!    Input, real ( kind = 8 ) U(M), the U vector.
!
!    Input/output, real ( kind = 8 ) V(N).  On input, V must contain the vector V.
!    On output V contains the information necessary to recover the Givens
!    rotations GV described above.
!
!    Output, real ( kind = 8 ) W(M), contains information necessary to
!    recover the Givens rotations GW described above.
!
!    Output, logical SING, is set to TRUE if any of the diagonal elements
!    of the output S are zero.  Otherwise SING is set FALSE.
!
  implicit none

  integer ( kind = 4 ) ls
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real    ( kind = 8 ) cos
  real    ( kind = 8 ) cotan
  real    ( kind = 8 ) giant
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jj
  integer ( kind = 4 ) l
  real    ( kind = 8 ) s(ls)
  real    ( kind = 8 ) sin
  logical              sing
  real    ( kind = 8 ) tan
  real    ( kind = 8 ) tau
  real    ( kind = 8 ) temp
  real    ( kind = 8 ) u(m)
  real    ( kind = 8 ) v(n)
  real    ( kind = 8 ) w(m)
!
!  GIANT is the largest magnitude.
!
  giant = huge ( giant )
!
!  Initialize the diagonal element pointer.
!
  jj = ( n * ( 2 * m - n + 1 ) ) / 2 - ( m - n )
!
!  Move the nontrivial part of the last column of S into W.
!
  l = jj
  do i = n, m
    w(i) = s(l)
    l = l + 1
  end do
!
!  Rotate the vector V into a multiple of the N-th unit vector
!  in such a way that a spike is introduced into W.
!
  do j = n-1, 1, -1

    jj = jj - ( m - j + 1 )
    w(j) = 0.0D+00

    if ( v(j) /= 0.0D+00 ) then
!
!  Determine a Givens rotation which eliminates the
!  J-th element of V.
!
      if ( abs ( v(n) ) < abs ( v(j) ) ) then
        cotan = v(n) / v(j)
        sin = 0.5D+00 / sqrt ( 0.25D+00 + 0.25D+00 * cotan**2 )
        cos = sin * cotan
        tau = 1.0D+00
        if ( abs ( cos ) * giant > 1.0D+00 ) then
          tau = 1.0D+00 / cos
        end if
      else
        tan = v(j) / v(n)
        cos = 0.5D+00 / sqrt ( 0.25D+00 + 0.25D+00 * tan**2 )
        sin = cos * tan
        tau = sin
      end if
!
!  Apply the transformation to V and store the information
!  necessary to recover the Givens rotation.
!
      v(n) = sin * v(j) + cos * v(n)
      v(j) = tau
!
!  Apply the transformation to S and extend the spike in W.
!
      l = jj
      do i = j, m
        temp = cos * s(l) - sin * w(i)
        w(i) = sin * s(l) + cos * w(i)
        s(l) = temp
        l = l + 1
      end do

    end if

  end do
!
!  Add the spike from the rank 1 update to W.
!
   w(1:m) = w(1:m) + v(n) * u(1:m)
!
!  Eliminate the spike.
!
  sing = .false.

  do j = 1, n-1

    if ( w(j) /= 0.0D+00 ) then
!
!  Determine a Givens rotation which eliminates the
!  J-th element of the spike.
!
      if ( abs ( s(jj) ) < abs ( w(j) ) ) then

        cotan = s(jj) / w(j)
        sin = 0.5D+00 / sqrt ( 0.25D+00 + 0.25D+00 * cotan**2 )
        cos = sin * cotan

        if ( 1.0D+00 < abs ( cos ) * giant ) then
          tau = 1.0D+00 / cos
        else
          tau = 1.0D+00
        end if

      else

        tan = w(j) / s(jj)
        cos = 0.5D+00 / sqrt ( 0.25D+00 + 0.25D+00 * tan**2 )
        sin = cos * tan
        tau = sin

      end if
!
!  Apply the transformation to S and reduce the spike in W.
!
      l = jj
      do i = j, m
        temp = cos * s(l) + sin * w(i)
        w(i) = - sin * s(l) + cos * w(i)
        s(l) = temp
        l = l + 1
      end do
!
!  Store the information necessary to recover the Givens rotation.
!
      w(j) = tau

    end if
!
!  Test for zero diagonal elements in the output S.
!
    if ( s(jj) == 0.0D+00 ) then
      sing = .true.
    end if

    jj = jj + ( m - j + 1 )

  end do
!
!  Move W back into the last column of the output S.
!
  l = jj
  do i = n, m
    s(l) = w(i)
    l = l + 1
  end do

  if ( s(jj) == 0.0D+00 ) then
    sing = .true.
  end if

  return
end
