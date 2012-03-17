subroutine banfac ( w, nroww, nrow, nbandl, nbandu, iflag )

!*****************************************************************************80
!
!! BANFAC factors a banded matrix without pivoting.
!
!  Discussion:
!
!    BANFAC returns in W the LU-factorization, without pivoting, of 
!    the banded matrix A of order NROW with (NBANDL+1+NBANDU) bands 
!    or diagonals in the work array W.
! 
!    Gauss elimination without pivoting is used.  The routine is 
!    intended for use with matrices A which do not require row 
!    interchanges during factorization, especially for the totally 
!    positive matrices which occur in spline calculations.
!
!    The matrix storage mode used is the same one used by LINPACK 
!    and LAPACK, and results in efficient innermost loops.
! 
!    Explicitly, A has 
! 
!      NBANDL bands below the diagonal
!      1     main diagonal
!      NBANDU bands above the diagonal
!
!    and thus, with MIDDLE=NBANDU+1,
!    A(I+J,J) is in W(I+MIDDLE,J) for I=-NBANDU,...,NBANDL, J=1,...,NROW.
!
!    For example, the interesting entries of a banded matrix
!    matrix of order 9, with NBANDL=1, NBANDU=2:
!
!      11 12 13  0  0  0  0  0  0
!      21 22 23 24  0  0  0  0  0
!       0 32 33 34 35  0  0  0  0
!       0  0 43 44 45 46  0  0  0
!       0  0  0 54 55 56 57  0  0
!       0  0  0  0 65 66 67 68  0
!       0  0  0  0  0 76 77 78 79
!       0  0  0  0  0  0 87 88 89
!       0  0  0  0  0  0  0 98 99
!
!    would appear in the first 1+1+2=4 rows of W as follows:
!
!       0  0 13 24 35 46 57 68 79
!       0 12 23 34 45 56 67 78 89
!      11 22 33 44 55 66 77 88 99
!      21 32 43 54 65 76 87 98  0
! 
!    All other entries of W not identified in this way with an
!    entry of A are never referenced.
! 
!    This routine makes it possible to solve any particular linear system 
!    A*X=B for X by the call
!
!      call banslv ( w, nroww, nrow, nbandl, nbandu, b )
!
!    with the solution X contained in B on return.
! 
!    If IFLAG=2, then one of NROW-1, NBANDL, NBANDU failed to be nonnegative, 
!    or else one of the potential pivots was found to be zero 
!    indicating that A does not have an LU-factorization.  This 
!    implies that A is singular in case it is totally positive.
!
!  Modified:
!
!    14 February 2007
!
!  Author:
!
!    Carl DeBoor
!
!  Reference:
!
!    Carl DeBoor,
!    A Practical Guide to Splines,
!    Springer, 2001,
!    ISBN: 0387953663,
!    LC: QA1.A647.v27.
!
!  Parameters:
! 
!    Input/output, real ( kind = 8 ) W(NROWW,NROW).
!    On input, W contains the "interesting" part of a banded 
!    matrix A, with the diagonals or bands of A stored in the
!    rows of W, while columns of A correspond to columns of W. 
!    On output, W contains the LU-factorization of A into a unit 
!    lower triangular matrix L and an upper triangular matrix U 
!    (both banded) and stored in customary fashion over the 
!    corresponding entries of A.  
!
!    Input, integer ( kind = 4 ) NROWW, the row dimension of the work array W.
!    NROWW must be at least NBANDL+1 + NBANDU.
! 
!    Input, integer ( kind = 4 ) NROW, the number of rows in A.
!
!    Input, integer ( kind = 4 ) NBANDL, the number of bands of A below 
!    the main diagonal.
! 
!    Input, integer ( kind = 4 ) NBANDU, the number of bands of A above 
!    the main diagonal.
! 
!    Output, integer ( kind = 4 ) IFLAG, error flag.
!    1, success.
!    2, failure, the matrix was not factored.
!
  implicit none

  integer ( kind = 4 ) nrow
  integer ( kind = 4 ) nroww

  real    ( kind = 8 ) factor
  integer ( kind = 4 ) i
  integer ( kind = 4 ) iflag
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) middle
  integer ( kind = 4 ) nbandl
  integer ( kind = 4 ) nbandu
  real    ( kind = 8 ) pivot
  real    ( kind = 8 ) w(nroww,nrow)

  iflag = 1

  if ( nrow < 1 ) then
    iflag = 2
    return
  end if
!
!  W(MIDDLE,*) contains the main diagonal of A.
!
  middle = nbandu + 1
  
  if ( nrow == 1 ) then
    if ( w(middle,nrow) == 0.0D+00 ) then
      iflag = 2
    end if
    return
  end if
!
!  A is upper triangular.  Check that the diagonal is nonzero.
!
  if ( nbandl <= 0 ) then

    do i = 1, nrow-1
      if ( w(middle,i) == 0.0D+00 ) then
        iflag = 2
        return
      end if
    end do

    if ( w(middle,nrow) == 0.0D+00 ) then
      iflag = 2
    end if

    return
!
!  A is lower triangular.  Check that the diagonal is nonzero and
!  divide each column by its diagonal.
!
  else if ( nbandu <= 0 ) then

    do i = 1, nrow-1

      pivot = w(middle,i)

      if ( pivot == 0.0D+00 ) then
        iflag = 2
        return
      end if

      do j = 1, min ( nbandl, nrow-i )
        w(middle+j,i) = w(middle+j,i) / pivot
      end do

    end do

    return

  end if
!
!  A is not just a triangular matrix.  
!  Construct the LU factorization.
!
  do i = 1, nrow-1
!
!  W(MIDDLE,I) is the pivot for the I-th step.
!
    if ( w(middle,i) == 0.0D+00 ) then
      iflag = 2
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'BANFAC - Fatal error!'
      write ( *, '(a,i8)' ) '  Zero pivot encountered in column ', i
      stop
    end if
!
!  Divide each entry in column I below the diagonal by PIVOT.
!
    do j = 1, min ( nbandl, nrow-i )
      w(middle+j,i) = w(middle+j,i) / w(middle,i)
    end do
!
!  Subtract A(I,I+K)*(I-th column) from (I+K)-th column (below row I).
!
    do k = 1, min ( nbandu, nrow-i )
      factor = w(middle-k,i+k)
      do j = 1, min ( nbandl, nrow-i )
        w(middle-k+j,i+k) = w(middle-k+j,i+k) - w(middle+j,i) * factor
      end do
    end do
 
  end do
!
!  Check the last diagonal entry.
!
  if ( w(middle,nrow) == 0.0D+00 ) then
    iflag = 2
  end if

  return
end
subroutine banslv ( w, nroww, nrow, nbandl, nbandu, b )

!*****************************************************************************80
!
!! BANSLV solves a banded linear system A * X = B factored by BANFAC.
!
!  Modified:
!
!    14 February 2007
!
!  Author:
!
!    Carl DeBoor
!
!  Reference:
!
!    Carl DeBoor,
!    A Practical Guide to Splines,
!    Springer, 2001,
!    ISBN: 0387953663,
!    LC: QA1.A647.v27.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) W(NROWW,NROW).  W contains the banded matrix,
!    after it has been factored by BANFAC.
!
!    Input, integer ( kind = 4 ) NROWW, the row dimension of the work array W.
!    NROWW must be at least NBANDL+1 + NBANDU.
! 
!    Input, integer ( kind = 4 ) NROW, the number of rows in A.
!
!    Input, integer ( kind = 4 ) NBANDL, the number of bands of A below the 
!    main diagonal.
! 
!    Input, integer ( kind = 4 ) NBANDU, the number of bands of A above the 
!    main diagonal.
! 
!    Input/output, real ( kind = 8 ) B(NROW).
!    On input, B contains the right hand side of the system to be solved.
!    On output, B contains the solution, X.
!
  implicit none

  integer ( kind = 4 ) nrow
  integer ( kind = 4 ) nroww

  real ( kind = 8 ) b(nrow)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jmax
  integer ( kind = 4 ) middle
  integer ( kind = 4 ) nbandl
  integer ( kind = 4 ) nbandu
  real ( kind = 8 ) w(nroww,nrow)

  middle = nbandu + 1

  if ( nrow == 1 ) then
    b(1) = b(1) / w(middle,1)
    return
  end if
!
!  Forward pass:
!
!  For I = 1, 2, ..., NROW-1, subtract RHS(I)*(I-th column of L) 
!  from the right hand side, below the I-th row.
!
  if ( 0 < nbandl ) then
    do i = 1, nrow-1
      jmax = min ( nbandl, nrow-i )
      do j = 1, jmax
        b(i+j) = b(i+j) - b(i) * w(middle+j,i)
      end do
    end do
  end if
!
!  Backward pass:
!
!  For I=NROW, NROW-1,...,1, divide RHS(I) by 
!  the I-th diagonal entry of U, then subtract 
!  RHS(I)*(I-th column of U) from right hand side, above the I-th row.
!
  do i = nrow, 2, -1
   
    b(i) = b(i) / w(middle,i)

    do j = 1, min ( nbandu, i-1 )
      b(i-j) = b(i-j) - b(i) * w(middle-j,i)
    end do

  end do

  b(1) = b(1) / w(middle,1)

  return
end
subroutine bchfac ( w, nbands, nrow, diag )

!*****************************************************************************80
!
!! BCHFAC constructs a Cholesky factorization of a matrix.
!
!  Discussion:
!
!    The factorization has the form
!
!      C = L * D * L'
!  
!    with L unit lower triangular and D diagonal, for a given matrix C of 
!    order NROW, where C is symmetric positive semidefinite and banded, 
!    having NBANDS diagonals at and below the main diagonal.
! 
!    Gauss elimination is used, adapted to the symmetry and bandedness of C.
! 
!    Near-zero pivots are handled in a special way.  The diagonal 
!    element C(N,N) = W(1,N) is saved initially in DIAG(N), all N. 
! 
!    At the N-th elimination step, the current pivot element, W(1,N), 
!    is compared with its original value, DIAG(N).  If, as the result 
!    of prior elimination steps, this element has been reduced by about 
!    a word length, that is, if W(1,N) + DIAG(N) <= DIAG(N), then the pivot 
!    is declared to be zero, and the entire N-th row is declared to
!    be linearly dependent on the preceding rows.  This has the effect 
!    of producing X(N) = 0 when solving C * X = B for X, regardless of B.
! 
!    Justification for this is as follows.  In contemplated applications 
!    of this program, the given equations are the normal equations for 
!    some least-squares approximation problem, DIAG(N) = C(N,N) gives 
!    the norm-square of the N-th basis function, and, at this point, 
!    W(1,N) contains the norm-square of the error in the least-squares 
!    approximation to the N-th basis function by linear combinations 
!    of the first N-1.  
!
!    Having W(1,N)+DIAG(N) <= DIAG(N) signifies that the N-th function 
!    is linearly dependent to machine accuracy on the first N-1 
!    functions, therefore can safely be left out from the basis of 
!    approximating functions.
!
!    The solution of a linear system C * X = B is effected by the 
!    succession of the following two calls:
! 
!      call bchfac ( w, nbands, nrow, diag )
!
!      call bchslv ( w, nbands, nrow, b, x )
!
!  Modified:
!
!    14 February 2007
!
!  Author:
!
!    Carl DeBoor
!
!  Reference:
!
!    Carl DeBoor,
!    A Practical Guide to Splines,
!    Springer, 2001,
!    ISBN: 0387953663,
!    LC: QA1.A647.v27.
!
!  Parameters:
!
!    Input/output, real ( kind = 8 ) W(NBANDS,NROW).
!    On input, W contains the NBANDS diagonals in its rows, 
!    with the main diagonal in row 1.  Precisely, W(I,J) 
!    contains C(I+J-1,J), I=1,...,NBANDS, J=1,...,NROW.
!    For example, the interesting entries of a seven diagonal
!    symmetric matrix C of order 9 would be stored in W as
! 
!      11 22 33 44 55 66 77 88 99
!      21 32 43 54 65 76 87 98  *
!      31 42 53 64 75 86 97  *  *
!      41 52 63 74 85 96  *  *  *
!
!    Entries of the array not associated with an
!    entry of C are never referenced.
!    On output, W contains the Cholesky factorization 
!    C = L*D*L', with W(1,I) containing 1/D(I,I) and W(I,J) 
!    containing L(I-1+J,J), I=2,...,NBANDS.
!
!    Input, integer ( kind = 4 ) NBANDS, indicates the bandwidth of the
!    matrix C, that is, C(I,J) = 0 for NBANDS < abs(I-J).
! 
!    Input, integer ( kind = 4 ) NROW, is the order of the matrix C.
! 
!    Work array, real ( kind = 8 ) DIAG(NROW).
!
  implicit none

  integer ( kind = 4 ) nbands
  integer ( kind = 4 ) nrow

  real ( kind = 8 ) diag(nrow)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) imax
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jmax
  integer ( kind = 4 ) n
  real ( kind = 8 ) ratio
  real ( kind = 8 ) w(nbands,nrow)

  if ( nrow <= 1 ) then
    if ( 0.0D+00 < w(1,1) ) then
      w(1,1) = 1.0D+00 / w(1,1)
    end if
    return
  end if
!
!  Store the diagonal.
!
  diag(1:nrow) = w(1,1:nrow)
!
!  Factorization.
!
  do n = 1, nrow
 
    if ( w(1,n) + diag(n) <= diag(n) ) then
      w(1:nbands,n) = 0.0D+00
    else
 
      w(1,n) = 1.0D+00 / w(1,n)
 
      imax = min ( nbands - 1, nrow - n )
 
      jmax = imax
 
      do i = 1, imax
 
        ratio = w(i+1,n) * w(1,n)
 
        do j = 1, jmax
          w(j,n+i) = w(j,n+i) - w(j+i,n) * ratio
        end do
 
        jmax = jmax-1
        w(i+1,n) = ratio
 
      end do
 
    end if
 
  end do
 
  return
end
subroutine bchslv ( w, nbands, nrow, b )

!*****************************************************************************80
!
!! BCHSLV solves a banded symmetric positive definite system.
!
!  Discussion:
!
!    The system is of the form:
!
!      C * X = B 
!  
!    and the Cholesky factorization of C has been constructed 
!    by BCHFAC.
! 
!    With the factorization 
!
!      C = L * D * L'
!
!    available, where L is unit lower triangular and D is diagonal, 
!    the triangular system 
!
!      L * Y = B 
!
!    is solved for Y (forward substitution), Y is stored in B, the 
!    vector D**(-1)*Y is computed and stored in B, then the 
!    triangular system L'*X = D**(-1)*Y is solved for X 
!    (back substitution).
!
!  Modified:
!
!    14 February 2007
!
!  Author:
!
!    Carl DeBoor
!
!  Reference:
!
!    Carl DeBoor,
!    A Practical Guide to Splines,
!    Springer, 2001,
!    ISBN: 0387953663,
!    LC: QA1.A647.v27.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) W(NBANDS,NROW), the Cholesky factorization for C, 
!    as computed by BCHFAC.
! 
!    Input, integer ( kind = 4 ) NBANDS, the bandwidth of C.
!
!    Input, integer ( kind = 4 ) NROW, the order of the matrix C.
! 
!    Input/output, real ( kind = 8 ) B(NROW).
!    On input, the right hand side.
!    On output, the solution.
!
  implicit none

  integer ( kind = 4 ) nbands
  integer ( kind = 4 ) nrow

  real ( kind = 8 ) b(nrow)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) n
  real ( kind = 8 ) w(nbands,nrow)

  if ( nrow <= 1 ) then
    b(1) = b(1) * w(1,1)
    return
  end if
!
!  Forward substitution. 
!  Solve L*Y = B.
!
  do n = 1, nrow

    do j = 1, min ( nbands - 1, nrow - n )
      b(j+n) = b(j+n) - w(j+1,n) * b(n)
    end do

  end do
!
!  Back substitution. 
!  Solve L'*X = D**(-1)*Y.
!
  do n = nrow, 1, -1

    b(n) = b(n) * w(1,n)

    do j = 1, min ( nbands - 1, nrow - n )
      b(n) = b(n) - w(j+1,n) * b(j+n)
    end do

  end do

  return
end
subroutine bsplpp ( t, bcoef, n, k, scrtch, break, coef, l )

!*****************************************************************************80
!
!! BSPLPP converts from B-spline to piecewise polynomial form.
!
!  Discussion:
!
!    The B-spline representation of a spline is ( T, BCOEF, N, K ),
!    while the piecewise polynomial representation is 
!    ( BREAK, COEF, L, K ).
!
!    For each breakpoint interval, the K relevant B-spline coefficients 
!    of the spline are found and then differenced repeatedly to get the 
!    B-spline coefficients of all the derivatives of the spline on that 
!    interval. 
!
!    The spline and its first K-1 derivatives are then evaluated at the 
!    left end point of that interval, using BSPLVB repeatedly to obtain 
!    the values of all B-splines of the appropriate order at that point.
!
!  Modified:
!
!    14 February 2007
!
!  Author:
!
!    Carl DeBoor
!
!  Reference:
!
!    Carl DeBoor,
!    A Practical Guide to Splines,
!    Springer, 2001,
!    ISBN: 0387953663,
!    LC: QA1.A647.v27.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) T(N+K), the knot sequence.
! 
!    Input, real ( kind = 8 ) BCOEF(N), the B spline coefficient sequence.
! 
!    Input, integer ( kind = 4 ) N, the number of B spline coefficients.
! 
!    Input, integer ( kind = 4 ) K, the order of the spline.
! 
!    Work array, real ( kind = 8 ) SCRTCH(K,K).
! 
!    Output, real ( kind = 8 ) BREAK(L+1), the piecewise polynomial breakpoint 
!    sequence.  BREAK contains the distinct points in the sequence T(K:N+1)
! 
!    Output, real ( kind = 8 ) COEF(K,N), with COEF(I,J) = (I-1)st derivative 
!    of the spline at BREAK(J) from the right.
! 
!    Output, integer ( kind = 4 ) L, the number of polynomial pieces which 
!    make up the spline in the interval ( T(K), T(N+1) ).
!
  implicit none

  integer ( kind = 4 ) k
  integer ( kind = 4 ) l
  integer ( kind = 4 ) n

  real ( kind = 8 ) bcoef(n)
  real ( kind = 8 ) biatx(k)
  real ( kind = 8 ) break(*)
  real ( kind = 8 ) coef(k,n)
  real ( kind = 8 ) diff
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jp1
  integer ( kind = 4 ) left
  integer ( kind = 4 ) lsofar
  real ( kind = 8 ) scrtch(k,k)      
  real ( kind = 8 ) sum1
  real ( kind = 8 ) t(n+k)

  lsofar = 0
  break(1) = t(k)
  
  do left = k, n
!
!  Find the next nontrivial knot interval.
!
    if ( t(left+1) == t(left) ) then
      cycle
    end if

    lsofar = lsofar + 1
    break(lsofar+1) = t(left+1)

    if ( k <= 1 ) then
      coef(1,lsofar) = bcoef(left)
      cycle
    end if
!
!  Store the K B-spline coefficients relevant to current knot 
!  interval in SCRTCH(*,1).
!
    do i = 1, k
      scrtch(i,1) = bcoef(left-k+i)
    end do
!
!  For J=1,...,K-1, compute the  K-J  B-spline coefficients relevant to
!  the current knot interval for the J-th derivative by differencing
!  those for the (J-1)st derivative, and store in SCRTCH(.,J+1).
!
    do jp1 = 2, k
      j = jp1 - 1
      do i = 1, k - j
        diff = t(left+i) - t(left+i-(k-j))
        if ( 0.0D+00 < diff ) then
          scrtch(i,jp1) = ( ( scrtch(i+1,j) - scrtch(i,j) ) / diff ) &
            * real ( k - j, kind = 8 )
        end if
      end do
    end do
!
!  For J=0, ..., K-1, find the values at T(left) of the J+1
!  B-splines of order J+1 whose support contains the current
!  knot interval from those of order J (in  BIATX ), then combine
!  with the B-spline coefficients (in SCRTCH(.,K-J) ) found earlier
!  to compute the (K-J-1)st derivative at  T(LEFT) of the given
!  spline.
!
    call bsplvb ( t, 1, 1, t(left), left, biatx )

    coef(k,lsofar) = scrtch(1,k)
    
    do jp1 = 2, k
    
      call bsplvb ( t, jp1, 2, t(left), left, biatx )

      coef(k+1-jp1,lsofar) = dot_product ( biatx(1:jp1), scrtch(1:jp1,k+1-jp1) )
      
    end do

  end do
   
  l = lsofar

  return
end
subroutine bsplvb ( t, jhigh, index, x, left, biatx )

!*****************************************************************************80
!
!! BSPLVB evaluates B-splines at a point X with a given knot sequence.
!
!  Discusion:
!
!    BSPLVB evaluates all possibly nonzero B-splines at X of order
!
!      JOUT = MAX ( JHIGH, (J+1)*(INDEX-1) ) 
!  
!    with knot sequence T.
! 
!    The recurrence relation
! 
!                     X - T(I)               T(I+J+1) - X
!    B(I,J+1)(X) = ----------- * B(I,J)(X) + --------------- * B(I+1,J)(X)
!                  T(I+J)-T(I)               T(I+J+1)-T(I+1)
! 
!    is used to generate B(LEFT-J:LEFT,J+1)(X) from B(LEFT-J+1:LEFT,J)(X)
!    storing the new values in BIATX over the old. 
!
!    The facts that 
!
!      B(I,1)(X) = 1  if  T(I) <= X < T(I+1)
!
!    and that 
!
!      B(I,J)(X) = 0  unless  T(I) <= X < T(I+J)
!
!    are used. 
!
!    The particular organization of the calculations follows 
!    algorithm 8 in chapter X of the text.
!
!  Modified:
!
!    14 February 2007
!
!  Author:
!
!    Carl DeBoor
!
!  Reference:
!
!    Carl DeBoor,
!    A Practical Guide to Splines,
!    Springer, 2001,
!    ISBN: 0387953663,
!    LC: QA1.A647.v27.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) T(LEFT+JOUT), the knot sequence.  T is assumed to 
!    be nondecreasing, and also, T(LEFT) must be strictly less than 
!    T(LEFT+1).
! 
!    Input, integer ( kind = 4 ) JHIGH, INDEX, determine the order 
!    JOUT = max ( JHIGH, (J+1)*(INDEX-1) )  
!    of the B-splines whose values at X are to be returned.  
!    INDEX is used to avoid recalculations when several 
!    columns of the triangular array of B-spline values are
!    needed, for example, in BVALUE or in BSPLVD.
!    If INDEX = 1, the calculation starts from scratch and the entire 
!    triangular array of B-spline values of orders
!    1, 2, ...,JHIGH is generated order by order, that is, 
!    column by column.
!    If INDEX = 2, only the B-spline values of order J+1, J+2, ..., JOUT  
!    are generated, the assumption being that BIATX, J, 
!    DELTAL, DELTAR are, on entry, as they were on exit 
!    at the previous call.  In particular, if JHIGH = 0, 
!    then JOUT = J+1, that is, just the next column of B-spline 
!    values is generated.
!    Warning: the restriction  JOUT <= JMAX (= 20) is
!    imposed arbitrarily by the dimension statement for DELTAL
!    and DELTAR, but is nowhere checked for.
! 
!    Input, real ( kind = 8 ) X, the point at which the B-splines 
!    are to be evaluated.
! 
!    Input, integer ( kind = 4 ) LEFT, an integer chosen so that 
!    T(LEFT) <= X <= T(LEFT+1).
! 
!    Output, real ( kind = 8 ) BIATX(JOUT), with BIATX(I) containing the
!    value at X of the polynomial of order JOUT which agrees 
!    with the B-spline B(LEFT-JOUT+I,JOUT,T) on the interval 
!    (T(LEFT),T(LEFT+1)).
!
  implicit none

  integer ( kind = 4 ), parameter :: jmax = 20

  integer ( kind = 4 ) jhigh

  real ( kind = 8 ) biatx(jhigh)
  real ( kind = 8 ), save, dimension ( jmax ) :: deltal
  real ( kind = 8 ), save, dimension ( jmax ) :: deltar
  integer ( kind = 4 ) i
  integer ( kind = 4 ) index
  integer ( kind = 4 ), save :: j = 1
  integer ( kind = 4 ) left
  real ( kind = 8 ) saved
  real ( kind = 8 ) t(left+jhigh)
  real ( kind = 8 ) term
  real ( kind = 8 ) x

  if ( index == 1 ) then 
    j = 1
    biatx(1) = 1.0D+00
    if ( jhigh <= j ) then
      return
    end if
  end if

  if ( t(left+1) <= t(left) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'BSPLVB - Fatal error!'
    write ( *, '(a)' ) '  It is required that T(LEFT) < T(LEFT+1).'
    write ( *, '(a,i8)' ) '  But LEFT = ', left
    write ( *, '(a,g14.6)' ) '  T(LEFT) =   ', t(left)
    write ( *, '(a,g14.6)' ) '  T(LEFT+1) = ', t(left+1)
    stop
  end if

  do
   
    deltar(j) = t(left+j) - x
    deltal(j) = x - t(left+1-j)

    saved = 0.0D+00
    do i = 1, j
      term = biatx(i) / ( deltar(i) + deltal(j+1-i) )
      biatx(i) = saved + deltar(i) * term
      saved = deltal(j+1-i) * term
    end do

    biatx(j+1) = saved
    j = j + 1

    if ( jhigh <= j ) then
      exit
    end if

  end do

  return
end
subroutine bsplvd ( t, k, x, left, a, dbiatx, nderiv )

!*****************************************************************************80
!
!! BSPLVD calculates the nonvanishing B-splines and derivatives at X.
!
!  Discussion:
!
!    Values at X of all the relevant B-splines of order K:K+1-NDERIV 
!    are generated via BSPLVB and stored temporarily in DBIATX.  
!
!    Then the B-spline coefficients of the required derivatives 
!    of the B-splines of interest are generated by differencing, 
!    each from the preceding one of lower order, and combined with 
!    the values of B-splines of corresponding order in DBIATX 
!    to produce the desired values.
!
!  Modified:
!
!    14 February 2007
!
!  Author:
!
!    Carl DeBoor
!
!  Reference:
!
!    Carl DeBoor,
!    A Practical Guide to Splines,
!    Springer, 2001,
!    ISBN: 0387953663,
!    LC: QA1.A647.v27.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) T(LEFT+K), the knot sequence.  It is assumed that 
!    T(LEFT) < T(LEFT+1).  Also, the output is correct only if 
!    T(LEFT) <= X <= T(LEFT+1).
!
!    Input, integer ( kind = 4 ) K, the order of the B-splines to be evaluated.
! 
!    Input, real ( kind = 8 ) X, the point at which these values are sought.
! 
!    Input, integer ( kind = 4 ) LEFT, indicates the left endpoint of the 
!    interval of interest.  The K B-splines whose support contains the interval 
!    ( T(LEFT), T(LEFT+1) ) are to be considered.
! 
!    Workspace, real ( kind = 8 ) A(K,K).
! 
!    Output, real ( kind = 8 ) DBIATX(K,NDERIV).  DBIATX(I,M) contains 
!    the value of the (M-1)st derivative of the (LEFT-K+I)-th B-spline 
!    of order K for knot sequence T, I=M,...,K, M=1,...,NDERIV.
!
!    Input, integer ( kind = 4 ) NDERIV, indicates that values of B-splines and 
!    their derivatives up to but not including the NDERIV-th are asked for. 
!
  implicit none

  integer ( kind = 4 ) k
  integer ( kind = 4 ) left
  integer ( kind = 4 ) nderiv

  real    ( kind = 8 ) a(k,k)
  real    ( kind = 8 ) dbiatx(k,nderiv)
  real    ( kind = 8 ) factor
  real    ( kind = 8 ) fkp1mm
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ideriv
  integer ( kind = 4 ) il
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jlow
  integer ( kind = 4 ) jp1mid
  integer ( kind = 4 ) ldummy
  integer ( kind = 4 ) m
  integer ( kind = 4 ) mhigh
  real    ( kind = 8 ) sum1
  real    ( kind = 8 ) t(left+k)
  real    ( kind = 8 ) x

  mhigh = max ( min ( nderiv, k ), 1 )
!
!  MHIGH is usually equal to NDERIV.
!
  call bsplvb ( t, k+1-mhigh, 1, x, left, dbiatx )
  
  if ( mhigh == 1 ) then
    return
  end if
!
!  The first column of DBIATX always contains the B-spline values
!  for the current order.  These are stored in column K+1-current
!  order before BSPLVB is called to put values for the next
!  higher order on top of it.
!
  ideriv = mhigh
  do m = 2, mhigh
    jp1mid = 1
    do j = ideriv, k
      dbiatx(j,ideriv) = dbiatx(jp1mid,1)
      jp1mid = jp1mid + 1
    end do
    ideriv = ideriv - 1
    call bsplvb ( t, k+1-ideriv, 2, x, left, dbiatx )
  end do
!
!  At this point, B(LEFT-K+I, K+1-J)(X) is in DBIATX(I,J) for
!  I=J,...,K and J=1,...,MHIGH ('=' NDERIV). 
!
!  In particular, the first column of DBIATX is already in final form. 
!
!  To obtain corresponding derivatives of B-splines in subsequent columns, 
!  generate their B-representation by differencing, then evaluate at X.
!
  jlow = 1
  do i = 1, k
    a(jlow:k,i) = 0.0D+00
    jlow = i
    a(i,i) = 1.0D+00
  end do
!
!  At this point, A(.,J) contains the B-coefficients for the J-th of the
!  K B-splines of interest here.
!
  do m = 2, mhigh
  
    fkp1mm = real ( k + 1 - m, kind = 8 )
    il = left
    i = k
!
!  For J = 1,...,K, construct B-coefficients of (M-1)st derivative of
!  B-splines from those for preceding derivative by differencing
!  and store again in  A(.,J).  The fact that  A(I,J) = 0 for
!  I < J is used.
!
    do ldummy = 1, k+1-m
    
      factor = fkp1mm / ( t(il+k+1-m) - t(il) )
!
!  The assumption that T(LEFT) < T(LEFT+1) makes denominator
!  in FACTOR nonzero.
!
      a(i,1:i) = ( a(i,1:i) - a(i-1,1:i) ) * factor

      il = il - 1
      i = i - 1
      
    end do
!
!  For I = 1,...,K, combine B-coefficients A(.,I) with B-spline values
!  stored in DBIATX(.,M) to get value of (M-1)st derivative of
!  I-th B-spline (of interest here) at X, and store in DBIATX(I,M).
!
!  Storage of this value over the value of a B-spline
!  of order M there is safe since the remaining B-spline derivatives
!  of the same order do not use this value due to the fact
!  that  A(J,I) = 0  for J < I.
!
    do i = 1, k
 
      jlow = max ( i, m )
 
      dbiatx(i,m) = dot_product ( a(jlow:k,i), dbiatx(jlow:k,m) )

    end do
    
  end do
 
  return
end
subroutine bspp2d ( t, bcoef, n, k, m, scrtch, break, coef, l )

!*****************************************************************************80
!
!! BSPP2D converts from B-spline to piecewise polynomial representation.
!
!  Discussion:
!
!    The B-spline representation
! 
!      T, BCOEF(.,J), N, K 
!
!    is converted to its piecewise polynomial representation 
!
!      BREAK, COEF(J,.,.), L, K, J=1, ..., M.
!
!    This is an extended version of BSPLPP for use with tensor products.
!
!    For each breakpoint interval, the K relevant B-spline 
!    coefficients of the spline are found and then differenced 
!    repeatedly to get the B-spline coefficients of all the 
!    derivatives of the spline on that interval. 
!
!    The spline and its first K-1 derivatives are then evaluated 
!    at the left endpoint of that interval, using BSPLVB 
!    repeatedly to obtain the values of all B-splines of the 
!    appropriate order at that point.
!
!  Modified:
!
!    14 February 2007
!
!  Author:
!
!    Carl DeBoor
!
!  Reference:
!
!    Carl DeBoor,
!    A Practical Guide to Splines,
!    Springer, 2001,
!    ISBN: 0387953663,
!    LC: QA1.A647.v27.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) T(N+K), the knot sequence.
! 
!    Input, real ( kind = 8 ) BCOEF(N,M).  For each J, B(*,J) is the 
!    B-spline coefficient sequence, of length N.
! 
!    Input, integer ( kind = 4 ) N, the length of BCOEF.
! 
!    Input, integer ( kind = 4 ) K, the order of the spline.
! 
!    Input, integer ( kind = 4 ) M, the number of data sets.
! 
!    Work array, real ( kind = 8 ) SCRTCH(K,K,M).
! 
!    Output, real ( kind = 8 ) BREAK(L+1), the breakpoint sequence
!    containing the distinct points in the sequence T(K),...,T(N+1)
! 
!    Output, real ( kind = 8 ) COEF(M,K,N), with COEF(MM,I,J) = the (I-1)st
!    derivative of the MM-th spline at BREAK(J) from the right, MM=1, ..., M.
! 
!    Output, integer ( kind = 4 ) L, the number of polynomial pieces which 
!    make up the spline in the interval (T(K), T(N+1)).
!
  implicit none

  integer ( kind = 4 ) k
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real    ( kind = 8 ) bcoef(n,m)
  real    ( kind = 8 ) biatx(k)
  real    ( kind = 8 ) break(*)
  real    ( kind = 8 ) coef(m,k,*)
  real    ( kind = 8 ) diff
  real    ( kind = 8 ) fkmj
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jp1
  integer ( kind = 4 ) kmj
  integer ( kind = 4 ) l
  integer ( kind = 4 ) left
  integer ( kind = 4 ) lsofar
  integer ( kind = 4 ) mm
  real    ( kind = 8 ) scrtch(k,k,m)
  real    ( kind = 8 ) sum1
  real    ( kind = 8 ) t(n+k)

  lsofar = 0
  break(1) = t(k)
  
  do left = k, n
!
!  Find the next nontrivial knot interval.
!
    if ( t(left+1) == t(left) ) then
      cycle
    end if

    lsofar = lsofar+1
    break(lsofar+1) = t(left+1)

    if ( k <= 1 ) then
    
      coef(1:m,1,lsofar) = bcoef(left,1:m)

      cycle

    end if
!
!  Store the K B-spline coefficients relevant to current knot interval
!  in SCRTCH(.,1).
!
    do i = 1, k
      scrtch(i,1,1:m) = bcoef(left-k+i,1:m)
    end do
!
!  For J = 1,...,K-1, compute the ( K - J ) B-spline coefficients relevant to
!  current knot interval for the J-th derivative by differencing
!  those for the (J-1)st derivative, and store in SCRTCH(.,J+1).
!
    do jp1 = 2, k
    
      j = jp1 - 1
      kmj = k - j
      fkmj = real ( k - j, kind = 8 )
      
      do i = 1, k - j
      
        diff = ( t(left+i) - t(left+i-kmj) ) / fkmj
        
        if ( 0.0D+00 < diff ) then
        
          scrtch(i,jp1,1:m) = ( scrtch(i+1,j,1:m) - scrtch(i,j,1:m) ) / diff
          
        end if
        
      end do
        
    end do
!
!  For  J = 0, ..., K-1, find the values at T(LEFT) of the J+1
!  B-splines of order J+1 whose support contains the current
!  knot interval from those of order J (in  BIATX ), then combine
!  with the B-spline coefficients (in SCRTCH(.,K-J) ) found earlier
!  to compute the (K-J-1)st derivative at T(LEFT) of the given spline.
!
    call bsplvb ( t, 1, 1, t(left), left, biatx )
    
    coef(1:m,k,lsofar) = scrtch(1,k,1:m)
    
    do jp1 = 2, k
    
      call bsplvb ( t, jp1, 2, t(left), left, biatx )
      kmj = k + 1 - jp1
      
      do mm = 1, m
      
        sum1 = 0.0D+00
        do i = 1, jp1
          sum1 = sum1 + biatx(i) * scrtch(i,kmj,mm)
        end do
        
        coef(mm,kmj,lsofar) = sum1
        
      end do
      
    end do
   
  end do
   
  l = lsofar
  
  return
end
function bvalue ( t, bcoef, n, k, x, jderiv )

!*****************************************************************************80
!
!! BVALUE evaluates a derivative of a spline from its B-spline representation.  
!
!  Discussion:
!
!    The spline is taken to be continuous from the right.
! 
!    The nontrivial knot interval (T(I),T(I+1)) containing X is 
!    located with the aid of INTERV.  The K B-spline coefficients 
!    of F relevant for this interval are then obtained from BCOEF, 
!    or are taken to be zero if not explicitly available, and are 
!    then differenced JDERIV times to obtain the B-spline 
!    coefficients of (D**JDERIV)F relevant for that interval.  
!
!    Precisely, with J = JDERIV, we have from X.(12) of the text that:
! 
!      (D**J)F = sum ( BCOEF(.,J)*B(.,K-J,T) )
! 
!    where
!                      / BCOEF(.),                    if J == 0
!                     /
!       BCOEF(.,J) = / BCOEF(.,J-1) - BCOEF(.-1,J-1)
!                   / -----------------------------,  if 0 < J
!                  /    (T(.+K-J) - T(.))/(K-J)
! 
!    Then, we use repeatedly the fact that
! 
!      sum ( A(.) * B(.,M,T)(X) ) = sum ( A(.,X) * B(.,M-1,T)(X) )
! 
!    with
!                   (X - T(.))*A(.) + (T(.+M-1) - X)*A(.-1)
!      A(.,X) =   ---------------------------------------
!                   (X - T(.))      + (T(.+M-1) - X)
! 
!    to write (D**J)F(X) eventually as a linear combination of 
!    B-splines of order 1, and the coefficient for B(I,1,T)(X) 
!    must then be the desired number (D**J)F(X).
!    See Chapter X, (17)-(19) of text.
!
!  Modified:
!
!    14 February 2007
!
!  Author:
!
!    Carl DeBoor
!
!  Reference:
!
!    Carl DeBoor,
!    A Practical Guide to Splines,
!    Springer, 2001,
!    ISBN: 0387953663,
!    LC: QA1.A647.v27.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) T(N+K), the knot sequence.  T is assumed 
!    to be nondecreasing.
!
!    Input, real ( kind = 8 ) BCOEF(N), B-spline coefficient sequence.
! 
!    Input, integer ( kind = 4 ) N, the length of BCOEF.
! 
!    Input, integer ( kind = 4 ) K, the order of the spline.
!
!    Input, real ( kind = 8 ) X, the point at which to evaluate.
! 
!    Input, integer ( kind = 4 ) JDERIV, the order of the derivative to 
!    be evaluated.  JDERIV is assumed to be zero or positive.
! 
!    Output, real ( kind = 8 ) BVALUE, the value of the (JDERIV)-th 
!    derivative of the spline at X.
!
  implicit none

  integer ( kind = 4 ) k
  integer ( kind = 4 ) n

  real ( kind = 8 ) aj(k)
  real ( kind = 8 ) bcoef(n)
  real ( kind = 8 ) bvalue
  real ( kind = 8 ) dl(k)
  real ( kind = 8 ) dr(k)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ilo
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jc
  integer ( kind = 4 ) jcmax
  integer ( kind = 4 ) jcmin
  integer ( kind = 4 ) jderiv
  integer ( kind = 4 ) jj
  integer ( kind = 4 ) mflag
  real ( kind = 8 ) t(n+k)
  real ( kind = 8 ) x

  bvalue = 0.0D+00
  
  if ( k <= jderiv ) then
    return
  end if
!
!  Find I so that 1 <= I < N+K and T(I) < T(I+1) and T(I) <= X < T(I+1). 
!
!  If no such I can be found, X lies outside the support of the 
!  spline F and  BVALUE = 0.  The asymmetry in this choice of I makes F 
!  right continuous.
!
  call interv ( t, n+k, x, i, mflag )
  
  if ( mflag /= 0 ) then
    return
  end if
!
!  If K = 1 (and JDERIV = 0), BVALUE = BCOEF(I).
!
  if ( k <= 1 ) then
    bvalue = bcoef(i)
    return
  end if
!
!  Store the K B-spline coefficients relevant for the knot interval
!  ( T(I),T(I+1) ) in AJ(1),...,AJ(K) and compute DL(J) = X - T(I+1-J),
!  DR(J) = T(I+J)-X, J=1,...,K-1.  Set any of the AJ not obtainable
!  from input to zero.  
!
!  Set any T's not obtainable equal to T(1) or to T(N+K) appropriately.
!
  jcmin = 1
  
  if ( k <= i ) then
  
    do j = 1, k-1
      dl(j) = x - t(i+1-j)
    end do
    
  else
  
    jcmin = 1 - ( i - k )
 
    do j = 1, i
      dl(j) = x - t(i+1-j)
    end do
 
    do j = i, k-1
      aj(k-j) = 0.0D+00
      dl(j) = dl(i)
    end do
  
  end if
 
  jcmax = k

  if ( n < i ) then
 
    jcmax = k + n - i
    do j = 1, k + n - i
      dr(j) = t(i+j) - x
    end do
 
    do j = k+n-i, k-1
      aj(j+1) = 0.0D+00
      dr(j) = dr(k+n-i)
    end do
 
  else
 
    do j = 1, k-1
      dr(j) = t(i+j) - x
    end do
 
  end if
 
  do jc = jcmin, jcmax
    aj(jc) = bcoef(i-k+jc)
  end do
!
!  Difference the coefficients JDERIV times.
!
  do j = 1, jderiv
 
    ilo = k - j
    do jj = 1, k - j
      aj(jj) = ( ( aj(jj+1) - aj(jj) ) / ( dl(ilo) + dr(jj) ) ) &
        * real ( k - j, kind = 8 )
      ilo = ilo - 1
    end do
 
  end do
!
!  Compute value at X in (T(I),T(I+1)) of JDERIV-th derivative,
!  given its relevant B-spline coefficients in AJ(1),...,AJ(K-JDERIV).
!
  do j = jderiv+1, k-1
    ilo = k-j
    do jj = 1, k-j
      aj(jj) = ( aj(jj+1) * dl(ilo) + aj(jj) * dr(jj) ) &
        / ( dl(ilo) + dr(jj) )
      ilo = ilo - 1
    end do
  end do
  
  bvalue = aj(1)
 
  return
end
subroutine chol1d ( p, v, qty, npoint, ncol, u, qu )

!*****************************************************************************80
!
!! CHOL1D sets up and solves linear systems needed by SMOOTH.
!
!  Discussion:
!
!    This routine constructs the upper three diagonals of
!
!      V(I,J), I = 2 to NPOINT-1, J=1,3, 
!
!    of the matrix 
!
!      6 * (1-P) * Q' * (D**2) * Q + P * R.
! 
!    It then computes its L*L' decomposition and stores it also
!    in V, then applies forward and back substitution to the right hand side 
!
!      Q'*Y 
!
!    in QTY to obtain the solution in U.
!
!  Modified:
!
!    16 February 2007
!
!  Author:
!
!    Carl DeBoor
!
!  Reference:
!
!    Carl DeBoor,
!    A Practical Guide to Splines,
!    Springer, 2001,
!    ISBN: 0387953663,
!    LC: QA1.A647.v27.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) P, the smoothing parameter that defines 
!    the linear system.
!
!    Input/output, real ( kind = 8 ) V(NPOINT,7), contains data used
!    to define the linear system, some of which is determined by
!    routine SETUPQ.
!
!    Input, real ( kind = 8 ) QTY(NPOINT), the value of Q' * Y.
!
!    Input, integer ( kind = 4 ) NPOINT, the number of equations.
!
!    Input, integer ( kind = 4 ) NCOL, an unused parameter, which may be 
!    set to 1.
!
!    Output, real ( kind = 8 ) U(NPOINT), the solution.
!
!    Output, real ( kind = 8 ) QU(NPOINT), the value of Q * U.
!
  implicit none

  integer ( kind = 4 ) npoint

  integer ( kind = 4 ) i
  integer ( kind = 4 ) ncol
  real    ( kind = 8 ) p
  real    ( kind = 8 ) qty(npoint)
  real    ( kind = 8 ) qu(npoint)
  real    ( kind = 8 ) u(npoint)
  real    ( kind = 8 ) v(npoint,7)
  real    ( kind = 8 ) prev
  real    ( kind = 8 ) ratio
  real    ( kind = 8 ) six1mp
  real    ( kind = 8 ) twop
!
!  Construct 6*(1-P)*Q'*(D**2)*Q + P*R.
!
  six1mp = 6.0D+00 * ( 1.0D+00 - p )
  twop = 2.0D+00 * p

  v(2:npoint-1,1) = six1mp * v(2:npoint-1,5) &
    + twop * ( v(1:npoint-2,4) + v(2:npoint-1,4) )
  v(2:npoint-1,2) = six1mp * v(2:npoint-1,6) + p * v(2:npoint-1,4)
  v(2:npoint-1,3) = six1mp * v(2:npoint-1,7)
 
  if ( npoint < 4 ) then

    u(1) = 0.0D+00
    u(2) = qty(2) / v(2,1)
    u(3) = 0.0D+00
!
!  Factorization.
!
  else

    do i = 2, npoint-2
      ratio = v(i,2) / v(i,1)
      v(i+1,1) = v(i+1,1) - ratio * v(i,2)
      v(i+1,2) = v(i+1,2) - ratio * v(i,3)
      v(i,2) = ratio
      ratio = v(i,3) / v(i,1)
      v(i+2,1) = v(i+2,1) - ratio * v(i,3)
      v(i,3) = ratio
    end do
!
!  Forward substitution
!
    u(1) = 0.0D+00
    v(1,3) = 0.0D+00
    u(2) = qty(2)
    do i = 2, npoint-2
      u(i+1) = qty(i+1) - v(i,2) * u(i) - v(i-1,3) * u(i-1)
    end do
!
!  Back substitution.
!
    u(npoint) = 0.0D+00
    u(npoint-1) = u(npoint-1) / v(npoint-1,1)

    do i = npoint-2, 2, -1
      u(i) = u(i) / v(i,1) - u(i+1) * v(i,2) - u(i+2) * v(i,3)
    end do

  end if
!
!  Construct Q * U.
!
  prev = 0.0D+00
  do i = 2, npoint
    qu(i) = ( u(i) - u(i-1) ) / v(i-1,4)
    qu(i-1) = qu(i) - prev
    prev = qu(i)
  end do

  qu(npoint) = -qu(npoint)
    
  return
end
subroutine colloc ( aleft, aright, lbegin, iorder, ntimes, addbrk, relerr )

!*****************************************************************************80
!
!! COLLOC solves an ordinary differential equation by collocation.
!
!  Method:
!
!    The M-th order ordinary differential equation with M side 
!    conditions, to be specified in subroutine DIFEQU, is solved 
!    approximately by collocation.
!
!    The approximation F to the solution G is piecewise polynomial of order 
!    K+M with L pieces and M-1 continuous derivatives.   F is determined by 
!    the requirement that it satisfy the differential equation at K points 
!    per interval (to be specified in COLPNT ) and the M side conditions.
!
!    This usually nonlinear system of equations for F is solved by
!    Newton's method. the resulting linear system for the B-coefficients of an
!    iterate is constructed appropriately in EQBLOK and then solved
!    in SLVBLK, a program designed to solve almost block
!    diagonal linear systems efficiently.
!
!    There is an opportunity to attempt improvement of the breakpoint
!    sequence, both in number and location, through the use of NEWNOT.
!
!    Printed output consists of the piecewise polynomial representation 
!    of the approximate solution, and of the error at selected points.
!
!  Modified:
!
!    14 February 2007
!
!  Author:
!
!    Carl DeBoor
!
!  Reference:
!
!    Carl DeBoor,
!    A Practical Guide to Splines,
!    Springer, 2001,
!    ISBN: 0387953663,
!    LC: QA1.A647.v27.
!
!  Parameters: 
!
!    Input, real ( kind = 8 ) ALEFT, ARIGHT, the endpoints of the interval.
!
!    Input, integer ( kind = 4 ) LBEGIN, the initial number of polynomial pieces 
!    in the approximation.  A uniform breakpoint sequence will be chosen.
!
!    Input, integer ( kind = 4 ) IORDER, the order of the polynomial pieces to be
!    used in the approximation
!
!    Input, integer ( kind = 4 ) NTIMES, the number of passes to be made 
!    through NEWNOT.
!
!    Input, real ( kind = 8 ) ADDBRK, the number, possibly fractional, of 
!    breaks to be added per pass through NEWNOT.  For instance, if 
!    ADDBRK = 0.33334, then a breakpoint will be added at every third pass 
!    through NEWNOT.
!
!    Input, real ( kind = 8 ) RELERR, a tolerance.  Newton iteration is 
!    stopped if the difference between the B-coefficients of two successive 
!    iterates is no more than RELERR*(absolute largest B-coefficient).
!
  implicit none

  integer ( kind = 4 ), parameter :: npiece = 100
  integer ( kind = 4 ), parameter :: ndim = 200
  integer ( kind = 4 ), parameter :: ncoef = 2000
  integer ( kind = 4 ), parameter :: lenblk = 2000

  real ( kind = 8 ) a(ndim)
  real ( kind = 8 ) addbrk
  real ( kind = 8 ) aleft
  real ( kind = 8 ) amax
  real ( kind = 8 ) aright
  real ( kind = 8 ) asave(ndim)
  real ( kind = 8 ) b(ndim)
  real ( kind = 8 ) bloks(lenblk)
  real ( kind = 8 ) break
  real ( kind = 8 ) coef
  real ( kind = 8 ) dx
  real ( kind = 8 ) err
  integer ( kind = 4 ) i
  integer ( kind = 4 ) iflag
  integer ( kind = 4 ) ii
  integer ( kind = 4 ) integs(3,npiece)
  integer ( kind = 4 ) iorder
  integer ( kind = 4 ) iside
  integer ( kind = 4 ) itemps(ndim)
  integer ( kind = 4 ) iter
  integer ( kind = 4 ) itermx
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) kpm
  integer ( kind = 4 ) l
  integer ( kind = 4 ) lbegin
  integer ( kind = 4 ) lnew
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
  integer ( kind = 4 ) nbloks
  integer ( kind = 4 ) nt
  integer ( kind = 4 ) ntimes
  real ( kind = 8 ) relerr
  real ( kind = 8 ) rho
  real ( kind = 8 ) t(ndim)
  real ( kind = 8 ) templ(lenblk)
  real ( kind = 8 ) temps(ndim)
  real ( kind = 8 ) xside

  equivalence ( bloks, templ )

  save / approx /
  save / other /
  save / side /

  common / approx / break(npiece), coef(ncoef), l, kpm
  common / other / itermx, k, rho(19)
  common / side / m, iside, xside(10)

  kpm = iorder

  if ( ncoef < lbegin * kpm ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'COLLOC - Fatal error!'
    write ( *, '(a)' ) '  The assigned dimension for COEF is too small.'
    stop
  end if
!
!  Set the various parameters concerning the particular differential
!  equation, including a first approximation in case the differential
!  equation is to be solved by iteration ( 0 < ITERMX ).
!
  call difequ ( 1, temps(1), temps )
!
!  Obtain the K collocation points for the standard interval.
!
  k = kpm - m
  call colpnt ( k, rho )
!
!  The following five statements could be replaced by a read in
!  order to obtain a nonuniform spacing of the breakpoints.
!
  dx = ( aright - aleft ) / real ( lbegin, kind = 8 )
 
  temps(1) = aleft
  do i = 2, lbegin
    temps(i) = temps(i-1) + dx
  end do
  temps(lbegin+1) = aright
!
!  Generate the required knots T(1:N+KPM).
!
  call knots ( temps, lbegin, kpm, m, t, n )
  nt = 1
!
!  Generate the almost block diagonal coefficient matrix BLOKS and
!  right hand side B from collocation equations and side conditions.
!
!  Then solve via SLVBLK, obtaining the B-representation of the 
!  approximation in T, A, N, KPM.
!
  do

    call eqblok ( t, n, kpm, temps, a, bloks, lenblk, integs, nbloks, b )

    call slvblk ( bloks, integs, nbloks, b, itemps, a, iflag )
!
!  Save B-spline coefficients of current approximation in ASAVE, then 
!  get new approximation and compare with old. 
!
!  If coefficients are more than RELERR apart (relatively) or if number 
!  of iterations is less than ITERMX, continue iterating.
!
    do iter = 1, itermx

      call bsplpp ( t, a, n, kpm, templ, break, coef, l )
 
      asave(1:n) = a(1:n)
 
      call eqblok ( t, n, kpm, temps, a, bloks, lenblk, integs, nbloks, b )

      call slvblk ( bloks, integs, nbloks, b, itemps, a, iflag )

      amax = maxval ( abs ( a(1:n) ) )
      err = maxval ( abs ( a(1:n) - asave(1:n) ) )
 
      if ( err <= relerr * amax ) then
        exit
      end if

    end do
!
!  Iteration (if any) completed.  Print out approximation based on current
!  breakpoint sequence, then try to improve the sequence.
!
    write ( *, '(a)' ) ' '
    write ( *,'(a,i3,a,i3,a)' ) &
    '  Approximation from a space of splines of order ', kpm, &
    ' on ', l, ' intervals'
    write ( *, '(a,i4)' ) '  of dimension ', n
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Breakpoints:'
    write ( *, '(a)' ) ' '
    write ( *, '(5g14.6)' ) break(2:l)

    if ( 0 < itermx ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' ) '  Results on interation ', iter
    end if

    call bsplpp ( t, a, n, kpm, templ, break, coef, l )
  
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) &
      '  The piecewise polynomial representation of the approximation:'
    write ( *, '(a)' ) ' '
   
    do i = 1, l
      ii = ( i - 1 ) * kpm
      write ( *, '(f9.3,2x,e12.4,10e11.3)' ) break(i), coef(ii+1:ii+kpm)
    end do
!
!  The following call is provided here for possible further analysis
!  of the approximation specific to the problem being solved.
!  It is, of course, easily omitted.
!
    call difequ ( 4, temps(1), temps )
 
    if ( ntimes < nt ) then
      exit
    end if
!
!  From the piecewise polynomial representation of the current approximation, 
!  obtain in NEWNOT a new, and possibly better, sequence of breakpoints, 
!  adding, on average, ADDBRK breakpoints per pass through NEWNOT.
!
    lnew = lbegin + int ( real ( nt, kind = 8 ) * addbrk )

    if ( ncoef < lnew * kpm ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'COLLOC - Fatal error!'
      write ( *, '(a)' ) '  The assigned dimension for COEF is too small.'
      stop
    end if

    call newnot ( break, coef, l, kpm, temps, lnew, templ )

    call knots ( temps, lnew, kpm, m, t, n )
    nt = nt + 1

  end do

  return  
end
subroutine colpnt ( k, rho )

!*****************************************************************************80
!
!! COLPNT supplies collocation points.
!
!  Discussion:
!
!    The collocation points are for the standard interval (-1,1) as the 
!    zeros of the Legendre polynomial of degree K, provided K <= 8.  
!
!    Otherwise, uniformly spaced points are given.
!
!  Modified:
!
!    14 February 2007
!
!  Author:
!
!    Carl DeBoor
!
!  Reference:
!
!    Carl DeBoor,
!    A Practical Guide to Splines,
!    Springer, 2001,
!    ISBN: 0387953663,
!    LC: QA1.A647.v27.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) K, the number of collocation points desired.
!
!    Output, real ( kind = 8 ) RHO(K), the collocation points.
!
  implicit none

  integer ( kind = 4 ) k

  integer ( kind = 4 ) j
  real ( kind = 8 ) rho(k)

  if ( k == 1 ) then

    rho(1) = 0.0D+00

  else if ( k == 2 ) then

    rho(1) = -0.577350269189626D+00
    rho(2) =  0.577350269189626D+00

  else if ( k == 3 ) then

    rho(1) = -0.774596669241483D+00
    rho(2) =  0.0D+00
    rho(3) =  0.774596669241483D+00

  else if ( k == 4 ) then

    rho(1) = -0.861136311594053D+00
    rho(2) = -0.339981043584856D+00
    rho(3) =  0.339981043584856D+00
    rho(4) =  0.861136311594053D+00

  else if ( k == 5 ) then

    rho(1) = -0.906179845938664D+00
    rho(2) = -0.538469310105683D+00
    rho(3) =  0.0D+00
    rho(4) =  0.538469310105683D+00
    rho(5) =  0.906179845938664D+00

  else if ( k == 6 ) then

    rho(1) = -0.932469514203152D+00
    rho(2) = -0.661209386466265D+00
    rho(3) = -0.238619186083197D+00
    rho(4) =  0.238619186083197D+00
    rho(5) =  0.661209386466265D+00
    rho(6) =  0.932469514203152D+00

  else if ( k == 7 ) then

    rho(1) = -0.949107912342759D+00
    rho(2) = -0.741531185599394D+00
    rho(3) = -0.405845151377397D+00
    rho(4) =  0.0D+00
    rho(5) =  0.405845151377397D+00
    rho(6) =  0.741531185599394D+00
    rho(7) =  0.949107912342759D+00

  else if ( k == 8 ) then

    rho(1) = -0.960289856497536D+00
    rho(2) = -0.796666477413627D+00
    rho(3) = -0.525532409916329D+00
    rho(4) = -0.183434642495650D+00
    rho(5) =  0.183434642495650D+00
    rho(6) =  0.525532409916329D+00
    rho(7) =  0.796666477413627D+00
    rho(8) =  0.960289856497536D+00

  else

    write ( *, '(a)' ) ' '
    write ( *, '(a)' )'COLPNT - Warning!'
    write ( *, '(a)' )'  Equispaced collocation points will be used,'
    write ( *, '(a,i8)' ) '  because K = ', k
 
    do j = 1, k
      rho(j) = ( real ( k - j,     kind = 8 ) * ( -1.0D+00 )   &
               + real (     j - 1, kind = 8 ) * ( +1.0D+00 ) ) &
               / real ( k     - 1, kind = 8 )
    end do
    
  end if
 
  return
end
subroutine cubspl ( tau, c, n, ibcbeg, ibcend )

!*****************************************************************************80
!
!! CUBSPL defines an interpolatory cubic spline.
!
!  Discussion:
!
!    A tridiagonal linear system for the unknown slopes S(I) of
!    F at TAU(I), I=1,..., N, is generated and then solved by Gauss
!    elimination, with S(I) ending up in C(2,I), for all I.
!
!  Modified:
!
!    14 February 2007
!
!  Author:
!
!    Carl DeBoor
!
!  Reference:
!
!    Carl DeBoor,
!    A Practical Guide to Splines,
!    Springer, 2001,
!    ISBN: 0387953663,
!    LC: QA1.A647.v27.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) TAU(N), the abscissas or X values of
!    the data points.  The entries of TAU are assumed to be
!    strictly increasing.
!
!    Input, integer ( kind = 4 ) N, the number of data points.  N is
!    assumed to be at least 2.
!
!    Input/output, real ( kind = 8 ) C(4,N).
!    On input, if IBCBEG or IBCBEG is 1 or 2, then C(2,1)
!    or C(2,N) should have been set to the desired derivative
!    values, as described further under IBCBEG and IBCEND.
!    On output, C contains the polynomial coefficients of
!    the cubic interpolating spline with interior knots
!    TAU(2) through TAU(N-1).
!    In the interval interval (TAU(I), TAU(I+1)), the spline
!    F is given by
!      F(X) = 
!        C(1,I) + 
!        C(2,I) * H +
!        C(3,I) * H**2 / 2 + 
!        C(4,I) * H**3 / 6.
!    where H=X-TAU(I).  The routine PPVALU may be used to
!    evaluate F or its derivatives from TAU, C, L=N-1,
!    and K=4.
!
!    Input, integer ( kind = 4 ) IBCBEG, IBCEND, boundary condition indicators.
!    IBCBEG = 0 means no boundary condition at TAU(1) is given.
!    In this case, the "not-a-knot condition" is used.  That
!    is, the jump in the third derivative across TAU(2) is
!    forced to zero.  Thus the first and the second cubic
!    polynomial pieces are made to coincide.
!    IBCBEG = 1 means the slope at TAU(1) is to equal the
!    input value C(2,1).
!    IBCBEG = 2 means the second derivative at TAU(1) is
!    to equal C(2,1).
!    IBCEND = 0, 1, or 2 has analogous meaning concerning the
!    boundary condition at TAU(N), with the additional
!    information taken from C(2,N).
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) c(4,n)
  real ( kind = 8 ) divdf1
  real ( kind = 8 ) divdf3
  real ( kind = 8 ) dtau
  real ( kind = 8 ) g
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ibcbeg
  integer ( kind = 4 ) ibcend
  real ( kind = 8 ) tau(n)
!
!  C(3,*) and C(4,*) are used initially for temporary storage.
!
!  Store first differences of the TAU sequence in C(3,*).
!
!  Store first divided difference of data in C(4,*).
!
  do i = 2, n
    c(3,i) = tau(i) - tau(i-1)
  end do

  do i = 2, n 
    c(4,i) = ( c(1,i) - c(1,i-1) ) / ( tau(i) - tau(i-1) )
  end do
!
!  Construct the first equation from the boundary condition
!  at the left endpoint, of the form:
!
!    C(4,1) * S(1) + C(3,1) * S(2) = C(2,1)
!
!  IBCBEG = 0: Not-a-knot
!
  if ( ibcbeg == 0 ) then

    if ( n <= 2 ) then
      c(4,1) = 1.0D+00
      c(3,1) = 1.0D+00
      c(2,1) = 2.0D+00 * c(4,2)
      go to 120
    end if

    c(4,1) = c(3,3)
    c(3,1) = c(3,2) + c(3,3)
    c(2,1) = ( ( c(3,2) + 2.0D+00 * c(3,1) ) * c(4,2) * c(3,3) &
      + c(3,2)**2 * c(4,3) ) / c(3,1)
!
!  IBCBEG = 1: derivative specified.
!
  else if ( ibcbeg == 1 ) then

    c(4,1) = 1.0D+00
    c(3,1) = 0.0D+00

    if ( n == 2 ) then
      go to 120
    end if
!
!  Second derivative prescribed at left end.
!
  else

    c(4,1) = 2.0D+00
    c(3,1) = 1.0D+00
    c(2,1) = 3.0D+00 * c(4,2) - c(3,2) / 2.0D+00 * c(2,1)

    if ( n == 2 ) then
      go to 120
    end if

  end if
!
!  If there are interior knots, generate the corresponding
!  equations and carry out the forward pass of Gauss elimination,
!  after which the I-th equation reads:
!
!    C(4,I) * S(I) + C(3,I) * S(I+1) = C(2,I).
!
  do i = 2, n-1
    g = -c(3,i+1) / c(4,i-1)
    c(2,i) = g * c(2,i-1) + 3.0D+00 * ( c(3,i) * c(4,i+1) + c(3,i+1) * c(4,i) )
    c(4,i) = g * c(3,i-1) + 2.0D+00 * ( c(3,i) + c(3,i+1))
  end do
!
!  Construct the last equation from the second boundary condition, of
!  the form
!
!    -G * C(4,N-1) * S(N-1) + C(4,N) * S(N) = C(2,N)
!
!  If slope is prescribed at right end, one can go directly to
!  back-substitution, since the C array happens to be set up just
!  right for it at this point.
!
  if ( ibcend == 1 ) then
    go to 160
  end if

  if ( 1 < ibcend ) then
    go to 110
  end if
 
90    continue
!
!  Not-a-knot and 3 <= N, and either 3 < N or also not-a-knot
!  at left end point.
!
  if ( n /= 3 .or. ibcbeg /= 0 ) then
    g = c(3,n-1) + c(3,n)
    c(2,n) = ( ( c(3,n) + 2.0D+00 * g ) * c(4,n) * c(3,n-1) + c(3,n)**2 &
      * ( c(1,n-1) - c(1,n-2) ) / c(3,n-1) ) / g
    g = - g / c(4,n-1)
    c(4,n) = c(3,n-1)
    c(4,n) = c(4,n) + g * c(3,n-1)
    c(2,n) = ( g * c(2,n-1) + c(2,n) ) / c(4,n)
    go to 160
  end if
!
!  N = 3 and not-a-knot also at left.
!
100   continue
 
  c(2,n) = 2.0D+00 * c(4,n)
  c(4,n) = 1.0D+00
  g = -1.0D+00 / c(4,n-1)
  c(4,n) = c(4,n) - c(3,n-1) / c(4,n-1)
  c(2,n) = ( g * c(2,n-1) + c(2,n) ) / c(4,n)
  go to 160
!
!  IBCEND = 2: Second derivative prescribed at right endpoint.
!
110   continue
 
  c(2,n) = 3.0D+00 * c(4,n) + c(3,n) / 2.0D+00 * c(2,n)
  c(4,n) = 2.0D+00
  g = -1.0D+00 / c(4,n-1)
  c(4,n) = c(4,n) - c(3,n-1) / c(4,n-1)
  c(2,n) = ( g * c(2,n-1) + c(2,n) ) / c(4,n)
  go to 160
!
!  N = 2.
!
120   continue
  
  if ( ibcend == 2  ) then

    c(2,n) = 3.0D+00 * c(4,n) + c(3,n) / 2.0D+00 * c(2,n)
    c(4,n) = 2.0D+00
    g = -1.0D+00 / c(4,n-1)
    c(4,n) = c(4,n) - c(3,n-1) / c(4,n-1)
    c(2,n) = ( g * c(2,n-1) + c(2,n) ) / c(4,n)
 
  else if ( ibcend == 0 .and. ibcbeg /= 0 ) then

    c(2,n) = 2.0D+00 * c(4,n)
    c(4,n) = 1.0D+00
    g = -1.0D+00 / c(4,n-1)
    c(4,n) = c(4,n) - c(3,n-1) / c(4,n-1)
    c(2,n) = ( g * c(2,n-1) + c(2,n) ) / c(4,n)

  else if ( ibcend == 0 .and. ibcbeg == 0 ) then

    c(2,n) = c(4,n)

  end if
!
!  Back solve the upper triangular system 
!
!    C(4,I) * S(I) + C(3,I) * S(I+1) = B(I)
!
!  for the slopes C(2,I), given that S(N) is already known.
!
160   continue
 
  do i = n-1, 1, -1
    c(2,i) = ( c(2,i) - c(3,i) * c(2,i+1) ) / c(4,i)
  end do
!
!  Generate cubic coefficients in each interval, that is, the
!  derivatives at its left endpoint, from value and slope at its
!  endpoints.
!
  do i = 2, n
    dtau = c(3,i)
    divdf1 = ( c(1,i) - c(1,i-1) ) / dtau
    divdf3 = c(2,i-1) + c(2,i) - 2.0D+00 * divdf1
    c(3,i-1) = 2.0D+00 * ( divdf1 - c(2,i-1) - divdf3 ) / dtau
    c(4,i-1) = 6.0D+00 * divdf3 / dtau**2
  end do
 
  return
end
subroutine cwidth ( w, b, nequ, ncols, integs, nbloks, d, x, iflag )

!*****************************************************************************80
!
!! CWIDTH solves an almost block diagonal linear system.
!
!  Discussion:
!
!    This routine is a variation of the theme in the algorithm
!    by Martin and Wilkinson.  It solves the linear system
!
!      A * X = B
!
!    of NEQU equations in case A is almost block diagonal with all
!    blocks having NCOLS columns using no more storage than it takes to
!    store the interesting part of A.  Such systems occur in the determination
!    of the B-spline coefficients of a spline approximation.
!
!  Block structure of A:  
!
!    The interesting part of A is taken to consist of NBLOKS
!    consecutive blocks, with the I-th block made up of NROWI = INTEGS(1,I)
!    consecutive rows and NCOLS consecutive columns of A, and with
!    the first LASTI = INTEGS(2,I) columns to the left of the next block.
!    These blocks are stored consecutively in the work array W.
!
!    For example, here is an 11th order matrix and its arrangement in
!    the work array W.  (The interesting entries of A are indicated by
!    their row and column index modulo 10.)
!
!                   ---   A   ---                          ---   W   ---
! 
!                      NROW1=3
!           11 12 13 14                                     11 12 13 14
!           21 22 23 24                                     21 22 23 24
!           31 32 33 34      NROW2=2                        31 32 33 34
!    LAST1=2      43 44 45 46                               43 44 45 46
!                 53 54 55 56         NROW3=3               53 54 55 56
!          LAST2=3         66 67 68 69                      66 67 68 69
!                          76 77 78 79                      76 77 78 79
!                          86 87 88 89   NROW4=1            86 87 88 89
!                   LAST3=1   97 98 99 90   NROW5=2         97 98 99 90
!                      LAST4=1   08 09 00 01                08 09 00 01
!                                18 19 10 11                18 19 10 11
!                         LAST5=4
!
!    For this interpretation of A as an almost block diagonal matrix,
!    we have NBLOKS = 5, and the INTEGS array is
!
!                          I = 1   2   3   4   5
!                    K =
!    INTEGS(K,I) =      1      3   2   3   1   2
!                       2      2   3   1   1   4
!
!
!  Method:
!
!    Gauss elimination with scaled partial pivoting is used, but
!    multipliers are not saved in order to save storage.  Rather, the
!    right hand side is operated on during elimination.  The two parameters
!    IPVTEQ and LASTEQ are used to keep track of the action.  IPVTEQ 
!    is the index of the variable to be eliminated next, from equations 
!    IPVTEQ+1,...,LASTEQ, using equation IPVTEQ, possibly after an 
!    interchange, as the pivot equation.  
!
!    The entries in the pivot column are always in column
!    1 of W.  This is accomplished by putting the entries in rows
!    IPVTEQ+1,...,LASTEQ revised by the elimination of the IPVTEQ-th
!    variable one to the left in W.  In this way, the columns of the
!    equations in a given block, as stored in W, will be aligned with
!    those of the next block at the moment when these next equations 
!    become involved in the elimination process.
!
!    Thus, for the above example, the first elimination steps proceed
!    as follows.
!
!    *11 12 13 14    11 12 13 14    11 12 13 14    11 12 13 14
!    *21 22 23 24   *22 23 24       22 23 24       22 23 24
!    *31 32 33 34   *32 33 34      *33 34          33 34
!     43 44 45 46    43 44 45 46   *43 44 45 46   *44 45 46        
!     53 54 55 56    53 54 55 56   *53 54 55 56   *54 55 56
!     66 67 68 69    66 67 68 69    66 67 68 69    66 67 68 69
!
!    In all other respects, the procedure is standard, including the
!    scaled partial pivoting.
!
!  Modified:
!
!    14 February 2007
!
!  Author:
!
!    Carl DeBoor
!
!  Reference:
!
!    Carl DeBoor,
!    A Practical Guide to Splines,
!    Springer, 2001,
!    ISBN: 0387953663,
!    LC: QA1.A647.v27.
!
!    Roger Martin, James Wilkinson,
!    Solution of Symmetric and Unsymmetric Band Equations and
!    the Calculation of Eigenvectors of Band Matrices,
!    Numerische Mathematik,
!    Volume 9, Number 4, December 1976, pages 279-301.
!
!  Parameters:
!
!    Input/output, real ( kind = 8 ) W(NEQU,NCOLS), on input, contains
!    the interesting part of the almost block diagonal coefficient matrix 
!    A.  The array INTEGS describes the storage scheme.  On output, W 
!    contains the upper triangular factor U of the LU factorization of a 
!    possibly permuted version of A.  In particular, the determinant of 
!    A could now be found as
!      IFLAG * W(1,1) * W(2,1) * ... * W(NEQU,1).
!
!    Input/output, real ( kind = 8 ) B(NEQU); on input, the right hand
!    side of the linear system.  On output, B has been overwritten by
!    other information.
!
!    Input, integer ( kind = 4 ) NEQU, the number of equations.
!
!    Input, integer ( kind = 4 ) NCOLS, the block width, that is, the number of 
!    columns in each block.
!
!    Input, integer ( kind = 4 ) INTEGS(2,NEQU), describes the block structure 
!    of A.
!    INTEGS(1,I) = number of rows in block I = NROW.
!    INTEGS(2,I) = number of elimination steps in block I = overhang over 
!    next block = LAST.
!
!    Input, integer ( kind = 4 ) NBOKS, the number of blocks.
!
!    Workspace, real D(NEQU), used to contain row sizes.  If storage is 
!    scarce, the array X could be used in the calling sequence for D.
!
!    Output, real ( kind = 8 ) X(NEQU), the computed solution, if 
!    IFLAG is nonzero.
!
!    Output, integer ( kind = 4 ) IFLAG, error flag.
!    = (-1)**(number of interchanges during elimination) if A is invertible;
!    = 0 if A is singular.
!
  implicit none

  integer ( kind = 4 ) nbloks
  integer ( kind = 4 ) ncols
  integer ( kind = 4 ) nequ

  real ( kind = 8 ) awi1od
  real ( kind = 8 ) b(nequ)
  real ( kind = 8 ) colmax
  real ( kind = 8 ) d(nequ)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) icount
  integer ( kind = 4 ) iflag
  integer ( kind = 4 ) ii
  integer ( kind = 4 ) integs(2,nbloks)
  integer ( kind = 4 ) ipvteq
  integer ( kind = 4 ) ipvtp1
  integer ( kind = 4 ) istar
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jmax
  integer ( kind = 4 ) lastcl
  integer ( kind = 4 ) lasteq
  integer ( kind = 4 ) lasti
  integer ( kind = 4 ) nexteq
  integer ( kind = 4 ) nrowad
  real ( kind = 8 ) ratio
  real ( kind = 8 ) rowmax
  real ( kind = 8 ) sum1
  real ( kind = 8 ) temp
  real ( kind = 8 ) w(nequ,ncols)
  real ( kind = 8 ) x(nequ)

  iflag = 1
  ipvteq = 0
  lasteq = 0
!
!  The I loop runs over the blocks.
!
  do i = 1, nbloks
!
!  The equations for the current block are added to those currently
!  involved in the elimination process, by increasing LASTEQ
!  by INTEGS(1,I) after the row size of these equations has been
!  recorded in the array D.
!
    nrowad = integs(1,i)
    
    do icount = 1, nrowad

      nexteq = lasteq + icount
      rowmax = maxval ( abs ( w(nexteq,1:ncols) ) )
      
      if ( rowmax == 0.0D+00 ) then
        iflag = 0
        return
      end if

      d(nexteq) = rowmax

    end do
   
    lasteq = lasteq + nrowad
!
!  There will be LASTI = INTEGS(2,I) elimination steps before
!  the equations in the next block become involved.
!
!  Further, LASTCL records the number of columns involved in the current
!  elimination step.  It starts equal to NCOLS when a block
!  first becomes involved and then drops by one after each elimination
!  step.
!
    lastcl = ncols
    lasti = integs(2,i)
    
    do icount = 1, lasti
    
      ipvteq = ipvteq + 1

      if ( lasteq <= ipvteq ) then

        if ( d(ipvteq) < abs ( w(ipvteq,1) ) + d(ipvteq) ) then
          exit
        end if

        iflag = 0
        return

      end if
!
!  Determine the smallest ISTAR in (IPVTEQ,LASTEQ) for
!  which abs ( W(ISTAR,1) ) / D(ISTAR) is as large as possible, and
!  interchange equations IPVTEQ and ISTAR in case  IPVTEQ < ISTAR.
!
      colmax = abs ( w(ipvteq,1) ) / d(ipvteq)
      istar = ipvteq
      ipvtp1 = ipvteq + 1
      
      do ii = ipvtp1, lasteq
        awi1od = abs ( w(ii,1) ) / d(ii)
        if ( colmax < awi1od ) then
          colmax = awi1od
          istar = ii
        end if
      end do
      
      if ( abs ( w(istar,1) ) + d(istar) == d(istar) ) then
        iflag = 0
        return
      end if
!
!  Rearrange data because of pivoting.
!
      if ( istar /= ipvteq ) then

        iflag = -iflag

        temp = d(istar)
        d(istar) = d(ipvteq)
        d(ipvteq) = temp
  
        temp = b(istar)
        b(istar) = b(ipvteq)
        b(ipvteq) = temp
      
        do j = 1, lastcl
          temp = w(istar,j)
          w(istar,j) = w(ipvteq,j)
          w(ipvteq,j) = temp
        end do

      end if
!
!  Subtract the appropriate multiple of equation IPVTEQ from
!  equations IPVTEQ+1,...,LASTEQ to make the coefficient of the
!  IPVTEQ-th unknown (presently in column 1 of W) zero, but
!  store the new coefficients in W one to the left from the old.
!  
      do ii = ipvtp1, lasteq
      
        ratio = w(ii,1) / w(ipvteq,1)
        do j = 2, lastcl
          w(ii,j-1) = w(ii,j) - ratio * w(ipvteq,j)
        end do
        w(ii,lastcl) = 0.0D+00
        b(ii) = b(ii) - ratio * b(ipvteq)
        
      end do
   
      lastcl = lastcl - 1
      
    end do

  end do
!
!  At this point, W and B contain an upper triangular linear system
!  equivalent to the original one, with W(I,J) containing entry
!  (I, I-1+J) of the coefficient matrix.  Solve this system by 
!  back substitution, taking into account its block structure.
!
!  I-loop over the blocks, in reverse order.
!
  i = nbloks

  do while ( 0 < i )

    lasti = integs(2,i)
    jmax = ncols - lasti
  
    do icount = 1, lasti
  
      sum1 = dot_product ( x(ipvteq+1:ipvteq+jmax), w(ipvteq,2:jmax+1) )
  
      x(ipvteq) = ( b(ipvteq) - sum1 ) / w(ipvteq,1)
      jmax = jmax + 1
      ipvteq = ipvteq - 1
    
    end do
  
    i = i - 1

  end do

  return
end
subroutine difequ ( mode, xx, v )

!*****************************************************************************80
!
!! DIFEQU returns information about a differential equation.
!
!  Discussion:
!
!    This sample version of DIFEQU is for the example in chapter XV.  It is a
!    nonlinear second order two point boundary value problem.
!
!  Modified:
!
!    16 February 2007
!
!  Author:
!
!    Carl DeBoor
!
!  Reference:
!
!    Carl DeBoor,
!    A Practical Guide to Splines,
!    Springer, 2001,
!    ISBN: 0387953663,
!    LC: QA1.A647.v27.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) MODE, an integer indicating the task to 
!    be performed.
!    1, initialization
!    2, evaluate the differential equation at point XX.
!    3, specify the next side condition
!    4, analyze the approximation
!
!    Input, real ( kind = 8 ) XX, a point at which information is wanted
!
!    Output, real ( kind = 8 ) V, depends on the MODE.
!
  implicit none

  integer ( kind = 4 ), parameter :: npiece = 100
  integer ( kind = 4 ), parameter :: ncoef = 2000

  real ( kind = 8 ) break
  real ( kind = 8 ) coef
  real ( kind = 8 ), save :: eps
  real ( kind = 8 ) ep1
  real ( kind = 8 ) ep2
  real ( kind = 8 ) error
  real ( kind = 8 ), save :: factor
  integer ( kind = 4 ) i
  integer ( kind = 4 ) iside
  integer ( kind = 4 ) itermx
  integer ( kind = 4 ) k
  integer ( kind = 4 ) kpm
  integer ( kind = 4 ) l
  integer ( kind = 4 ) m
  integer ( kind = 4 ) mode
  real ( kind = 8 ) ppvalu
  real ( kind = 8 ) rho
  real ( kind = 8 ), save :: s2ovep
  real ( kind = 8 ) solutn
  real ( kind = 8 ) un
  real ( kind = 8 ) v(20)
  real ( kind = 8 ) value
  real ( kind = 8 ) x
  real ( kind = 8 ) xside
  real ( kind = 8 ) xx

  save / approx /
  save / other /
  save / side /

  common / approx / break(npiece), coef(ncoef), l, kpm
  common / other / itermx, k, rho(19)
  common / side / m, iside, xside(10)
!
!  Initialize everything,  Set the order M of the differential equation, 
!  the nondecreasing sequence XSIDE(1:M), of points at which side 
!  conditions are given and anything else necessary.
!
  if ( mode == 1 ) then

    m = 2
    xside(1) = 0.0D+00
    xside(2) = 1.0D+00
!
!  Print out heading.
!
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Carrier''s nonlinear perturbation problem'
    write ( *, '(a)' ) ' '
  
    eps = 0.005D+00
    write ( *, '(a,g14.6)' ) '  EPS = ', eps
!
!  Set constants used in formula for solution below.
!
    factor = ( sqrt ( 2.0D+00 ) + sqrt ( 3.0D+00 ) )**2
    s2ovep = sqrt ( 2.0D+00 / eps )
!
!  Initial guess for Newton iteration: UN(X) = X*X-1.
!
    l = 1
    break(1) = 0.0D+00
    coef(1:kpm) = 0.0D+00
    coef(1) = -1.0D+00
    coef(3) = 2.0D+00
    itermx = 10
!
!  Provide value of left side coefficients and right hand side at XX.
!  Specifically, at XX the differential equation reads:
!
!    V(M+1) D**M + V(M) D**(M-1) + ... + V(1) D**0 = V(M+2)
!
!  in terms of the quantities V(1:M+2), to be computed here.
!
  else if ( mode == 2 ) then

    v(3) = eps
    v(2) = 0.0D+00

    un = ppvalu ( break, coef, l, kpm, xx, 0 )

    v(1) = 2.0D+00 * un
    v(4) = un**2 + 1.0D+00
!
!  Provide the M side conditions. these conditions are of the form
!
!    V(M+1) D**M + V(M) D**(M-1) + ... + V(1) D**0 = V(M+2)
!
!  in terms of the quantities V(1:M+2), to be specified here.
!  Note that V(M+1) = 0 for customary side conditions.
!
  else if ( mode == 3 ) then

    v(m+1) = 0.0D+00

    if ( iside == 1 ) then
      v(2) = 1.0D+00
      v(1) = 0.0D+00
      v(4) = 0.0D+00
      iside = iside + 1
    else if ( iside == 2 ) then
      v(2) = 0.0D+00
      v(1) = 1.0D+00
      v(4) = 0.0D+00
      iside = iside + 1
    end if
!
!  Calculate the error near the boundary layer at 1.
!
  else if ( mode == 4 ) then

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '      X             G(X)          G(X)-F(X):'
    write ( *, '(a)' ) ' '

    x = 0.75D+00
  
    do i = 1, 9

      ep1 = exp ( s2ovep * ( 1.0D+00 - x ) ) * factor
      ep2 = exp ( s2ovep * ( 1.0D+00 + x ) ) * factor

      solutn = 12.0D+00 / ( 1.0D+00 + ep1 )**2 * ep1 &
             + 12.0D+00 / ( 1.0D+00 + ep2 )**2 * ep2 - 1.0D+00

      value = ppvalu ( break, coef, l, kpm, x, 0 )

      error = solutn - value
      write ( *, '(2x,3g14.6)' ) x, solutn, error
      x = x + 0.03125D+00

    end do

  else

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'DIFEQU - Fatal error!'
    write ( *, '(a)' ) '  Illegal value of MODE:'
    write ( *, '(a,i8)' ) mode
    stop

  end if
  
  return
end
subroutine dtblok ( bloks, integs, nbloks, ipivot, iflag, detsgn, detlog )

!*****************************************************************************80
!
!! DTBLOK gets the determinant of an almost block diagonal matrix.
!
!  Discussion:
!
!    The matrix's PLU factorization must have been obtained 
!    previously by FCBLOK.
!
!    The logarithm of the determinant is computed instead of the
!    determinant itself to avoid the danger of overflow or underflow
!    inherent in this calculation.
!
!  Modified:
!
!    14 February 2007
!
!  Author:
!
!    Carl DeBoor
!
!  Reference:
!
!    Carl DeBoor,
!    A Practical Guide to Splines,
!    Springer, 2001,
!    ISBN: 0387953663,
!    LC: QA1.A647.v27.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) BLOKS(*), the factorization of A computed
!    by FCBLOK.
!
!    Input, integer ( kind = 4 ) INTEGS(3,NBLOKS), describes the block 
!    structure of A.
!
!    Input, integer ( kind = 4 ) NBLOKS, the number of blocks in A.
!
!    Input, integer ( kind = 4 ) IPIVOT(*), pivoting information.
!    The dimension of IPIVOT is the sum ( INTEGS(1,1:NBLOKS) ).
!
!    Input, integer ( kind = 4 ) IFLAG, = (-1)**(number of interchanges during
!    factorization) if successful, otherwise IFLAG = 0.
!
!    Output, real ( kind = 8 ) DETSGN, the sign of the determinant.
!
!    Output, real ( kind = 8 ) DETLOG, the natural logarithm of the 
!    determinant, if the determinant is not zero.  If the determinant
!    is 0, then DETLOG is returned as 0.
!
  implicit none

  integer ( kind = 4 ) nbloks

  real    ( kind = 8 ) bloks(*)
  real    ( kind = 8 ) detlog
  real    ( kind = 8 ) detsgn
  integer ( kind = 4 ) i
  integer ( kind = 4 ) iflag
  integer ( kind = 4 ) index
  integer ( kind = 4 ) indexp
  integer ( kind = 4 ) integs(3,nbloks)
  integer ( kind = 4 ) ip
  integer ( kind = 4 ) ipivot(1)
  integer ( kind = 4 ) k
  integer ( kind = 4 ) last
  integer ( kind = 4 ) nrow

  detsgn = iflag
  detlog = 0.0D+00

  if ( iflag == 0 ) then
    return
  end if

  index = 0
  indexp = 0
  
  do i = 1, nbloks
  
    nrow = integs(1,i)
    last = integs(3,i)
    
    do k = 1, last
      ip = index + nrow * ( k - 1 ) + ipivot(indexp+k)
      detlog = detlog + log ( abs ( bloks(ip) ) )
      detsgn = detsgn * sign ( 1.0D+00, bloks(ip) )
    end do
   
    index = nrow * integs(2,i) + index
    indexp = indexp + nrow
    
  end do
   
  return
end
subroutine eqblok ( t, n, kpm, work1, work2, bloks, lenblk, integs, nbloks, b )

!*****************************************************************************80
!
!! EQBLOK is to be called in COLLOC.
!
!  Method:
!
!    Each breakpoint interval gives rise to a block in the linear system.
!    This block is determined by the K collocation equations in the interval
!    with the side conditions, if any, in the interval interspersed
!    appropriately, and involves the KPM B-splines having the interval in
!    their support.  Correspondingly, such a block has NROW = K + ISIDEL
!    rows, with ISIDEL = number of side conditions in this and the 
!    previous intervals, and NCOL = KPM columns.
!
!    Further, because the interior knots have multiplicity K, we can
!    carry out in SLVBLK K elimination steps in a block before pivoting
!    might involve an equation from the next block.  In the last block,
!    of course, all KPM elimination steps will be carried out in SLVBLK.
!
!    See the detailed comments in SLVBLK for further
!    information about the almost block diagonal form used here.
!
!  Modified:
!
!    14 February 2007
!
!  Author:
!
!    Carl DeBoor
!
!  Reference:
!
!    Carl DeBoor,
!    A Practical Guide to Splines,
!    Springer, 2001,
!    ISBN: 0387953663,
!    LC: QA1.A647.v27.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) T(N+KPM), the knot sequence.
!
!    Input, integer ( kind = 4 ) N, the dimension of the approximating spline 
!    space, that is, the order of the linear system to be constructed.
!
!    Input, integer ( kind = 4 ) KPM, = K + M, the order of the approximating 
!    spline.
!
!    Input, integer ( kind = 4 ) LENBLK, the maximum length of the array BLOKS,
!    as allowed by the dimension statement in COLLOC.
!
!    Workspace, real ( kind = 8 ) WORK1(KPM,KPM), used in PUTIT.
!
!    Workspace, real ( kind = 8 ) WORK2(KPM,M+1), used in PUTIT.
!
!    Output, real ( kind = 8 ) BLOKS(*), the coefficient matrix of the
!    linear system, stored in almost block diagonal form, of size
!    KPM * sum ( INTEGS(1,1:NBLOKS) ).
!
!    Output, integer ( kind = 4 ) INTEGS(3,NBLOKS), describing the block 
!    structure.
!    INTEGS(1,I) = number of rows in block I;
!    INTEGS(2,I) = number of columns in block I;
!    INTEGS(3,I) = number of elimination steps which can be carried out in 
!    block I before pivoting might bring in an equation from the next block.
!
!    Output, integer ( kind = 4 ) NBLOKS, the number of blocks, equals number of 
!    polynomial pieces.
!
!    Output, real ( kind = 8 ) B(*), the right hand side of the linear 
!    system, stored corresponding to the almost block diagonal form, 
!    of size sum ( INTEGS(1,1:NBLOKS) ).
!
  implicit none

  integer ( kind = 4 ) kpm
  integer ( kind = 4 ) n

  real    ( kind = 8 ) b(*)
  real    ( kind = 8 ) bloks(*)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) index
  integer ( kind = 4 ) indexb
  integer ( kind = 4 ) integs(3,*)
  integer ( kind = 4 ) iside
  integer ( kind = 4 ) isidel
  integer ( kind = 4 ) itermx
  integer ( kind = 4 ) k
  integer ( kind = 4 ) left
  integer ( kind = 4 ) lenblk
  integer ( kind = 4 ) m
  integer ( kind = 4 ) nbloks
  integer ( kind = 4 ) nrow
  real    ( kind = 8 ) rho
  real    ( kind = 8 ) t(n+kpm)
  real    ( kind = 8 ) work1(kpm,kpm)
  real    ( kind = 8 ) work2(kpm,*)
  real    ( kind = 8 ) xside

  save / other /
  save / side /

  common / other / itermx, k, rho(19)
  common / side / m, iside, xside(10)

  index = 1
  indexb = 1
  i = 0
  iside = 1
  
  do left = kpm, n, k
  
    i = i + 1
!
!  Determine INTEGS(:,I).
!
    integs(2,i) = kpm
    
    if ( n <= left ) then

      integs(3,i) = kpm
      isidel = m
!
!  At this point, ISIDE - 1 gives the number of side conditions
!  incorporated so far.  Adding to this the side conditions in the
!  current interval gives the number ISIDEL.
!
    else

      integs(3,i) = k

      isidel = iside - 1

      do

        if ( isidel == m ) then
          exit
        end if

        if ( t(left+1) <= xside(isidel+1) ) then
          exit
        end if

        isidel = isidel + 1

      end do

    end if

    nrow = k + isidel
    integs(1,i) = nrow
!
!  The detailed equations for this block are generated and put
!  together in PUTIT.
!
    if ( lenblk < index + nrow * kpm - 1 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'EQBLOK - Fatal error!'
      write ( *, '(a)' ) '  The dimension of BLOKS is too small.'
      write ( *, '(a,i8)' ) '  LENBLK = ', lenblk
      stop
    end if

    call putit ( t, kpm, left, work1, work2, bloks(index), nrow, b(indexb) )

    index = index + nrow * kpm
    indexb = indexb + nrow
    
  end do
 
  nbloks = i

  return
end
subroutine evnnot ( break, coef, l, k, brknew, lnew, coefg )

!*****************************************************************************80
!
!! EVNNOT is a version of NEWNOT returning uniform knots.
!
!  Discussion:
!
!    EVNNOT returns LNEW+1 knots in BRKNEW which are evenly spaced between 
!    BREAK(1) and BREAK(L+1).
!
!  Modified:
!
!    14 February 2007
!
!  Author:
!
!    Carl DeBoor
!
!  Reference:
!
!    Carl DeBoor,
!    A Practical Guide to Splines,
!    Springer, 2001,
!    ISBN: 0387953663,
!    LC: QA1.A647.v27.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) BREAK(L+1), real ( kind = 8 ) COEF(K,L), 
!    integer ( kind = 4 ) L, integer K, the piecewise polynomial representation 
!    of a certain function F of order K.  Specifically,
!      d**(K-1) F(X) = COEF(K,I) for BREAK(I) <= X < BREAK(I+1).
!
!    Input, integer ( kind = 4 ) LNEW, the number of subintervals into which 
!    the interval (A,B) is to be sectioned by the new breakpoint sequence BRKNEW.
!
!    Output, real ( kind = 8 ) BRKNEW(LNEW+1), the new breakpoints.
!
!    Output, real (kind = 8 ) COEFG(2,L), the coefficient part of the 
!    piecewise polynomial representation BREAK, COEFG, L, 2 for the monotone 
!    piecewise linear function G with respect to which BRKNEW will
!    be equidistributed.
!
  implicit none

  integer ( kind = 4 ) k
  integer ( kind = 4 ) l
  integer ( kind = 4 ) lnew

  real    ( kind = 8 ) break(l+1)
  real    ( kind = 8 ) brknew(lnew+1)
  real    ( kind = 8 ) coef(k,l)
  real    ( kind = 8 ) coefg(2,l)
  integer ( kind = 4 ) i

  coefg(2,l) = 0.0D+00

  if ( lnew == 0 ) then

    brknew(1) = 0.5D+00 * ( break(1) + break(l+1) )

  else

    do i = 1, lnew + 1
      brknew(i) = ( real ( lnew - i + 1, kind = 8 ) * break(1) &
                  + real (        i - 1, kind = 8 ) * break(l+1) ) &
                  / real ( lnew,         kind = 8 )
    end do

  end if
 
  return
end
subroutine factrb ( w, ipivot, d, nrow, ncol, last, iflag )

!*****************************************************************************80
!
!! FACTRB constructs a partial PLU factorization.
!
!  Discussion:
!
!    This factorization corresponds to steps 1 through LAST in Gauss 
!    elimination for the matrix W of order ( NROW, NCOL ), using 
!    pivoting of scaled rows.
!
!  Modified:
!
!    14 February 2007
!
!  Author:
!
!    Carl DeBoor
!
!  Reference:
!
!    Carl DeBoor,
!    A Practical Guide to Splines,
!    Springer, 2001,
!    ISBN: 0387953663,
!    LC: QA1.A647.v27.
!
!  Parameters:
!
!    Input/output, real ( kind = 8 ) W(NROW,NCOL); on input, contains the
!    matrix to be partially factored; on output, the partial factorization.
!
!    Output, integer ( kind = 4 ) IPIVOT(NROW), contains a record of the pivoting 
!    strategy used; row IPIVOT(I) is used during the I-th elimination step,
!    for I = 1, ..., LAST.
!
!    Workspace, real ( kind = 8 ) D(NROW), used to store the maximum entry
!    in each row.
!
!    Input, integer ( kind = 4 ) NROW, the number of rows of W.
!
!    Input, integer ( kind = 4 ) NCOL, the number of columns of W.
!
!    Input, integer ( kind = 4 ) LAST, the number of elimination steps to 
!    be carried out.
!
!    Input/output, integer ( kind = 4 ) IFLAG.  On output, equals the input value 
!    times (-1)**(number of row interchanges during the factorization 
!    process), in case no zero pivot was encountered.
!    Otherwise, IFLAG = 0 on output.
!
  implicit none

  integer ( kind = 4 ) ncol
  integer ( kind = 4 ) nrow

  real    ( kind = 8 ) awikdi
  real    ( kind = 8 ) colmax
  real    ( kind = 8 ) d(nrow)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) iflag
  integer ( kind = 4 ) ipivi
  integer ( kind = 4 ) ipivk
  integer ( kind = 4 ) ipivot(nrow)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) kp1
  integer ( kind = 4 ) last
  real    ( kind = 8 ) ratio
  real    ( kind = 8 ) rowmax
  real    ( kind = 8 ) w(nrow,ncol)
!
!  Initialize IPIVOT and D.
!
  do i = 1, nrow
    ipivot(i) = i
  end do

  do i = 1, nrow
    
    rowmax = maxval ( abs ( w(i,1:ncol) ) )
    
    if ( rowmax == 0.0D+00 ) then
      iflag = 0
      return
    end if

    d(i) = rowmax
    
  end do
!
!  Gauss elimination with pivoting of scaled rows, loop over K = 1,..., LAST.
!
  k = 1
!
!  As pivot row for K-th step, pick among the rows not yet used,
!  that is, from rows IPIVOT(K:NROW), the one whose K-th entry, compared 
!  to the row size, is largest. 
!
!  If this row does not turn out to be row IPIVOT(K), redefine IPIVOT(K)
!  appropriately and record this interchange by changing the sign
!  of IFLAG.
!
  do while ( k <= last )

    ipivk = ipivot(k)

    if ( k == nrow ) then
      if ( abs ( w(ipivk,nrow) ) + d(ipivk) <= d(ipivk) ) then
        iflag = 0
      end if
      return
    end if

    j = k
    kp1 = k + 1
    colmax = abs ( w(ipivk,k) ) / d(ipivk)
!
!  Find the largest pivot.
!
    do i = kp1, nrow
      ipivi = ipivot(i)
      awikdi = abs ( w(ipivi,k) ) / d(ipivi)
      if ( colmax < awikdi ) then
        colmax = awikdi
        j = i
      end if
    end do
  
    if ( j /= k ) then
      ipivk = ipivot(j)
      ipivot(j) = ipivot(k)
      ipivot(k) = ipivk
      iflag = - iflag
    end if
!
!  If the pivot element is too small in absolute value, declare
!  the matrix to be noninvertible and quit.
!
    if ( abs ( w(ipivk,k) ) + d(ipivk) <= d(ipivk) ) then
      iflag = 0
      return
    end if
!
!  Otherwise, subtract the appropriate multiple of the pivot
!  row from the remaining rows, that is, the rows IPIVOT(K+1:NROW),
!  to make the K-th entry zero.
!
!  Save the multiplier in its place.
!
    do i = kp1, nrow
  
      ipivi = ipivot(i)
      w(ipivi,k) = w(ipivi,k) / w(ipivk,k)
    
      ratio = - w(ipivi,k)
      w(ipivi,kp1:ncol) = ratio * w(ipivk,kp1:ncol) + w(ipivi,kp1:ncol)
    
    end do
    
    k = kp1

  end do

  return
end
subroutine fcblok ( bloks, integs, nbloks, ipivot, scrtch, iflag )

!*****************************************************************************80
!
!! FCBLOK supervises the PLU factorization of an almost block diagonal matrix.
!
!  Discussion:
!
!    The routine supervises the PLU factorization with pivoting of
!    the scaled rows of an almost block diagonal matrix.
!
!    The almost block diagonal matrix is stored in the arrays
!    BLOKS and INTEGS.
!
!    The FACTRB routine carries out steps 1,..., LAST of Gauss
!    elimination, with pivoting, for an individual block.
!
!    The SHIFTB routine shifts the remaining rows to the top of
!    the next block.
!
!  Modified:
!
!    14 February 2007
!
!  Author:
!
!    Carl DeBoor
!
!  Reference:
!
!    Carl DeBoor,
!    A Practical Guide to Splines,
!    Springer, 2001,
!    ISBN: 0387953663,
!    LC: QA1.A647.v27.
!
!  Parameters:
!
!    Input/output, real ( kind = 8 ) BLOKS(*).  On input, the almost 
!    block diagonal matrix A to be factored.  On output, the
!    factorization of A.
!
!    Input, integer ( kind = 4 ) INTEGS(3,NBLOKS), describes the block 
!    structure of A.
!
!    Input, integer ( kind = 4 ) NBLOKS, the number of blocks in A.
!
!    Output, integer ( kind = 4 ) IPIVOT(*), which will contain pivoting 
!    information.  The dimension of IPIVOT is the sum ( INTEGS(1,1:NBLOKS) ).
!
!    Workspace, real SCRTCH(*), of length maxval ( integs(1,1:NBLOKS) ).
!
!    Output, integer ( kind = 4 ) IFLAG, error flag.
!    = 0,  in case matrix was found to be singular;
!    = (-1)**(number of row interchanges during factorization), otherwise.
!
  implicit none

  integer ( kind = 4 ) nbloks

  real    ( kind = 8 ) bloks(*)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) iflag
  integer ( kind = 4 ) index
  integer ( kind = 4 ) indexb
  integer ( kind = 4 ) indexn
  integer ( kind = 4 ) integs(3,nbloks)
  integer ( kind = 4 ) ipivot(*)
  integer ( kind = 4 ) last
  integer ( kind = 4 ) ncol
  integer ( kind = 4 ) nrow
  real    ( kind = 8 ) scrtch(*)

  iflag = 1
  indexb = 1
  indexn = 1
  i = 1 
!
!  Loop over the blocks.  I is the loop index.
!
  do

    index = indexn
    nrow = integs(1,i)
    ncol = integs(2,i)
    last = integs(3,i)
!
!  Carry out elimination on the I-th block until next block
!  enters, for columns 1 through LAST of I-th block.
!
    call factrb ( bloks(index), ipivot(indexb), scrtch, nrow, ncol, &
      last, iflag )
!
!  Check for having reached a singular block or the last block.
!
    if ( iflag == 0 .or. i == nbloks ) then
      exit
    end if

    i = i + 1
    indexn = nrow * ncol + index
!
!  Put the rest of the I-th block onto the next block.
!
    call shiftb ( bloks(index), ipivot(indexb), nrow, ncol, last, &
      bloks(indexn), integs(1,i), integs(2,i) )

    indexb = indexb + nrow

  end do

  return
end
subroutine interv ( xt, lxt, x, left, mflag )

!*****************************************************************************80
!
!! INTERV brackets a real value in an ascending vector of values.
!
!  Discussion:
!
!    The XT array is a set of increasing values.  The goal of the routine
!    is to determine the largest index I so that XT(I) <= X.
!
!    The routine is designed to be efficient in the common situation
!    that it is called repeatedly, with X taken from an increasing
!    or decreasing sequence.
!
!    This will happen when a piecewise polynomial is to be graphed.
!    The first guess for LEFT is therefore taken to be the value
!    returned at the previous call and stored in the local variable ILO.
!
!    A first check ascertains that ILO < LXT.  This is necessary
!    since the present call may have nothing to do with the previous
!    call.  Then, if 
!
!      XT(ILO) <= X < XT(ILO+1), 
!
!    we set LEFT = ILO and are done after just three comparisons.
!
!    Otherwise, we repeatedly double the difference ISTEP = IHI - ILO
!    while also moving ILO and IHI in the direction of X, until
!
!      XT(ILO) <= X < XT(IHI)
!
!    after which we use bisection to get, in addition, ILO + 1 = IHI.
!    The value LEFT = ILO is then returned.
!
!  Modified:
!
!    14 February 2007
!
!  Author:
!
!    Carl DeBoor
!
!  Reference:
!
!    Carl DeBoor,
!    A Practical Guide to Splines,
!    Springer, 2001,
!    ISBN: 0387953663,
!    LC: QA1.A647.v27.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) XT(LXT), a nondecreasing sequence of values.
!
!    Input, integer ( kind = 4 ) LXT, the dimension of XT.
!
!    Input, real ( kind = 8 ) X, the point whose location with 
!    respect to the sequence XT is to be determined.
!
!    Output, integer ( kind = 4 ) LEFT, the index of the bracketing value:
!      1     if             X  <  XT(1)
!      I     if   XT(I)  <= X  < XT(I+1)
!      LXT   if  XT(LXT) <= X
!
!    Output, integer ( kind = 4 ) MFLAG, indicates whether X lies within the
!    range of the data.
!    -1:            X  <  XT(1)
!     0: XT(I)   <= X  < XT(I+1)
!    +1: XT(LXT) <= X
!
  implicit none

  integer ( kind = 4 ) lxt

  integer ( kind = 4 ) left
  integer ( kind = 4 ) mflag
  integer ( kind = 4 ) ihi
  integer ( kind = 4 ), save :: ilo = 1
  integer ( kind = 4 ) istep
  integer ( kind = 4 ) middle
  real ( kind = 8 ) x
  real ( kind = 8 ) xt(lxt)

  ihi = ilo + 1

  if ( lxt <= ihi ) then

    if ( xt(lxt) <= x ) then
      go to 110
    end if

    if ( lxt <= 1 ) then
      mflag = -1
      left = 1
      return
    end if

    ilo = lxt - 1
    ihi = lxt

  end if

  if ( xt(ihi) <= x ) then
    go to 20
  end if

  if ( xt(ilo) <= x ) then
    mflag = 0
    left = ilo
    return
  end if
!
!  Now X < XT(ILO).  Decrease ILO to capture X.
!
  istep = 1

10 continue

  ihi = ilo
  ilo = ihi - istep

  if ( 1 < ilo ) then
    if ( xt(ilo) <= x ) then
      go to 50
    end if
    istep = istep * 2
    go to 10
  end if

  ilo = 1

  if ( x < xt(1) ) then
    mflag = -1
    left = 1
    return
  end if

  go to 50
!
!  Now XT(IHI) <= X.  Increase IHI to capture X.
!
20 continue

  istep = 1

30 continue

  ilo = ihi
  ihi = ilo + istep

  if ( ihi < lxt ) then

    if ( x < xt(ihi) ) then
      go to 50
    end if

    istep = istep * 2
    go to 30

  end if

  if ( xt(lxt) <= x ) then
    go to 110
  end if
!
!  Now XT(ILO) < = X < XT(IHI).  Narrow the interval.
!
  ihi = lxt

50 continue

  do

    middle = ( ilo + ihi ) / 2

    if ( middle == ilo ) then
      mflag = 0
      left = ilo
      return
    end if
!
!  It is assumed that MIDDLE = ILO in case IHI = ILO+1.
!
    if ( xt(middle) <= x ) then
      ilo = middle
    else
      ihi = middle
    end if

  end do
!
!  Set output and return.
!
110 continue

  mflag = 1

  if ( x == xt(lxt) ) then
    mflag = 0
  end if

  do left = lxt, 1, -1
    if ( xt(left) < xt(lxt) ) then
      return
    end if
  end do

  return
end
subroutine knots ( break, l, kpm, m, t, n )

!*****************************************************************************80
!
!! KNOTS is to be called in COLLOC.
!
!  Discussion:
!
!    Note that the FORTRAN77 calling sequence has been modified, by
!    adding the variable M.
!
!    From the given breakpoint sequence BREAK, this routine constructs the 
!    knot sequence T so that
!
!      SPLINE(K+M,T) = PP(K+M,BREAK) 
!
!    with M-1 continuous derivatives.
!
!    This means that T(1:N+KPM) is equal to BREAK(1) KPM times, then 
!    BREAK(2) through BREAK(L) each K times, then, finally, BREAK(L+1) 
!    KPM times.
!
!  Modified:
!
!    14 February 2007
!
!  Author:
!
!    Carl DeBoor
!
!  Reference:
!
!    Carl DeBoor,
!    A Practical Guide to Splines,
!    Springer, 2001,
!    ISBN: 0387953663,
!    LC: QA1.A647.v27.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) BREAK(L+1), the breakpoint sequence.
!
!    Input, integer ( kind = 4 ) L, the number of intervals or pieces.
!
!    Input, integer ( kind = 4 ) KPM, = K+M, the order of the piecewise polynomial
!    function or spline.
!
!    Input, integer ( kind = 4 ) M, the order of the differential equation.
!
!    Output, real ( kind = 8 ) T(N+KPM), the knot sequence.
!
!    Output, integer ( kind = 4 ) N, = L*K+M = the dimension of SPLINE(K+M,T).
!
  implicit none

  integer ( kind = 4 ) kpm
  integer ( kind = 4 ) l
  integer ( kind = 4 ) n

  real ( kind = 8 ) break(l+1)
  integer ( kind = 4 ) iside
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jj
  integer ( kind = 4 ) jjj
  integer ( kind = 4 ) k
  integer ( kind = 4 ) ll
  integer ( kind = 4 ) m
  real ( kind = 8 ) t(n+kpm)
  real ( kind = 8 ) xside

  k = kpm - m
  n = l * k + m
  jj = n + kpm
  jjj = l + 1
  
  do ll = 1, kpm
    t(jj) = break(jjj)
    jj = jj - 1
  end do
  
  do j = 1, l
    jjj = jjj - 1
    do ll = 1, k
      t(jj) = break(jjj)
      jj = jj - 1
    end do
  end do
   
  t(1:kpm) = break(1)

  return
end
subroutine l2appr ( t, n, k, q, diag, bcoef )

!*****************************************************************************80
!
!! L2APPR constructs a weighted L2 spline approximation to given data.
!
!  Discussion:
!
!    The routine constructs the weighted discrete L2-approximation by 
!    splines of order K with knot sequence T(1:N+K) to 
!    given data points ( TAU(1:NTAU), GTAU(1:NTAU) ).  
!
!    The B-spline coefficients BCOEF of the approximating spline are 
!    determined from the normal equations using Cholesky's method.
!
!  Method:
!
!    The B-spline coefficients of the L2-approximation are determined as the 
!    solution of the normal equations, for 1 <= I <= N:
!
!      sum ( 1 <= J <= N ) ( B(I), B(J) ) * BCOEF(J) = ( B(I), G ).
!
!    Here, B(I) denotes the I-th B-spline, G denotes the function to
!    be approximated, and the inner product of two functions F and G 
!    is given by
!
!      ( F, G ) = sum ( 1 <= I <= NTAU ) WEIGHT(I) * F(TAU(I)) * G(TAU(I)).
!
!    The arrays TAU and WEIGHT are given in common block DATA, as is the 
!    array GTAU(1:NTAU) = G(TAU(1:NTAU)).
!
!    The values of the B-splines B(1:N) are supplied by BSPLVB.
!
!    The coefficient matrix C, with
!
!       C(I,J) = ( B(I), B(J) )
!
!    of the normal equations is symmetric and (2*K-1)-banded, therefore
!    can be specified by giving its K bands at or below the diagonal. 
!
!    For I = 1:N and J = I:min(I+K-1,N), we store
!
!      ( B(I), B(J) ) = C(I,J) 
!
!    in
!
!      Q(I-J+1,J), 
!
!    and the right hand side
!
!      ( B(I), G )  
!
!    in
!
!      BCOEF(I).
!
!    Since B-spline values are most efficiently generated by finding
!    simultaneously the value of every nonzero B-spline at one point,
!    the entries of C (that is, of Q), are generated by computing, for
!    each LL, all the terms involving TAU(LL) simultaneously and adding
!    them to all relevant entries.
!
!  Modified:
!
!    14 February 2007
!
!  Author:
!
!    Carl DeBoor
!
!  Reference:
!
!    Carl DeBoor,
!    A Practical Guide to Splines,
!    Springer, 2001,
!    ISBN: 0387953663,
!    LC: QA1.A647.v27.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) T(N+K), the knot sequence.
!
!    Input, integer ( kind = 4 ) N, the dimension of the space of splines 
!    of order K with knots T.
!
!    Input, integer ( kind = 4 ) K, the order of the splines.
!
!    Workspace, real ( kind = 8 ) Q(K,N), used to store the K lower 
!    diagonals of the Gramian matrix C.
!
!    Workspace, real ( kind = 8 ) DIAG(N), used in BCHFAC.
!
!    Output, real ( kind = 8 ) BCOEF(N), the B-spline coefficients of
!    the L2 approximation to the data.
!
  implicit none

  integer ( kind = 4 ) k
  integer ( kind = 4 ) n
  integer ( kind = 4 ), parameter :: ntmax = 200

  real ( kind = 8 ) bcoef(n)
  real ( kind = 8 ) biatx(k)
  real ( kind = 8 ) diag(n)
  real ( kind = 8 ) dw
  real ( kind = 8 ) gtau
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jj
  integer ( kind = 4 ) left
  integer ( kind = 4 ) leftmk
  integer ( kind = 4 ) ll
  integer ( kind = 4 ) mm
  integer ( kind = 4 ) ntau
  real ( kind = 8 ) q(k,n)
  real ( kind = 8 ) t(n+k)
  real ( kind = 8 ) tau
  real ( kind = 8 ) totalw
  real ( kind = 8 ) weight

  save / i4data /
  save / r8data /

  common / i4data / ntau
  common / r8data / tau(ntmax), gtau(ntmax), weight(ntmax), totalw

  bcoef(1:n) = 0.0D+00
  q(1:k,1:n) = 0.0D+00

  left = k
  leftmk = 0
  
  do ll = 1, ntau
!
!  Locate LEFT such that TAU(LL) is in ( T(LEFT), T(LEFT+1) ).
!
    do

      if ( left == n ) then
        exit
      end if

      if ( tau(ll) < t(left+1) ) then
        exit
      end if

      left = left + 1
      leftmk = leftmk + 1

    end do

    call bsplvb ( t, k, 1, tau(ll), left, biatx )
!
!  BIATX(MM) contains the value of B(LEFT-K+MM) at TAU(LL).
!
!  Hence, with DW = BIATX(MM) * WEIGHT(LL), the number DW * GTAU(LL)
!  is a summand in the inner product
!
!    ( B(LEFT-K+MM), G)
!
!  which goes into  BCOEF(LEFT-K+MM)
!  and the number BIATX(JJ)*DW is a summand in the inner product
!    (B(LEFT-K+JJ), B(LEFT-K+MM)), into  Q(JJ-MM+1,LEFT-K+MM)
!  since  (LEFT-K+JJ)-(LEFT-K+MM)+1 = JJ - MM + 1.
!
    do mm = 1, k
    
      dw = biatx(mm) * weight(ll)
      j = leftmk + mm
      bcoef(j) = dw * gtau(ll) + bcoef(j)
      i = 1
      
      do jj = mm, k
        q(i,j) = biatx(jj) * dw + q(i,j)
        i = i + 1
      end do
      
    end do
    
  end do
!
!  Construct the Cholesky factorization for C in Q, then 
!  use it to solve the normal equations
!
!    C * X = BCOEF
!
!  for X, and store X in BCOEF.
!
  call bchfac ( q, k, n, diag )
  
  call bchslv ( q, k, n, bcoef )
  
  return
end
subroutine l2err ( iprfun, ftau, error )

!*****************************************************************************80
!
!! L2ERR computes the errors of an L2 approximation.
!
!  Discussion:
!
!    This routine computes various errors of the current L2 approximation, 
!    whose piecewise polynomial representation is contained in common 
!    block APPROX, to the given data contained in common block DATA.  
!
!    It prints out the average error ERRL1, the L2 error ERRL2, and the
!    maximum error ERRMAX.
!
!  Modified:
!
!    16 February 2007
!
!  Author:
!
!    Carl DeBoor
!
!  Reference:
!
!    Carl DeBoor,
!    A Practical Guide to Splines,
!    Springer, 2001,
!    ISBN: 0387953663,
!    LC: QA1.A647.v27.
!
!  Parameters: 
!
!    Input, integer ( kind = 4 ) IPRFUN.  If IPRFUN = 1, the routine prints out
!    the value of the approximation as well as its error at
!    every data point.
!
!    Output, real ( kind = 8 ) FTAU(NTAU), contains the value of the computed
!    approximation at each value TAU(1:NTAU).
!
!    Output, real ( kind = 8 ) ERROR(NTAU), with 
!      ERROR(I) = SCALE * ( G - F )(TAU(I)).  Here, SCALE equals 1
!    in case IPRFUN /= 1, or the absolute error is greater than 100 
!    somewhere.  Otherwise, SCALE is such that the maximum of the
!    absolute value of ERROR(1:NTAU) lies between 10 and 100.  This
!    makes the printed output more illustrative.
!
  implicit none

  integer ( kind = 4 ), parameter :: lpkmax = 100
  integer ( kind = 4 ), parameter :: ntmax = 200
  integer ( kind = 4 ), parameter :: ltkmax = 2000

  integer ( kind = 4 ) ntau

  real ( kind = 8 ) break
  real ( kind = 8 ) coef
  real ( kind = 8 ) err
  real ( kind = 8 ) errl1
  real ( kind = 8 ) errl2
  real ( kind = 8 ) errmax
  real ( kind = 8 ) error(ntau)
  real ( kind = 8 ) ftau(ntau)
  real ( kind = 8 ) gtau
  integer ( kind = 4 ) ie
  integer ( kind = 4 ) iprfun
  integer ( kind = 4 ) k
  integer ( kind = 4 ) l
  integer ( kind = 4 ) ll
  real ( kind = 8 ) ppvalu
  real ( kind = 8 ) scale
  real ( kind = 8 ) tau
  real ( kind = 8 ) totalw
  real ( kind = 8 ) weight

  save / approx /
  save / i4data /
  save / r8data /

  common / approx / break(lpkmax), coef(ltkmax), l, k
  common / i4data / ntau
  common / r8data / tau(ntmax), gtau(ntmax), weight(ntmax), totalw

  errl1 = 0.0D+00
  errl2 = 0.0D+00
  errmax = 0.0D+00

  do ll = 1, ntau

    ftau(ll) = ppvalu ( break, coef, l, k, tau(ll), 0 )

    error(ll) = gtau(ll) - ftau(ll)
    err = abs(error(ll))

    if ( errmax < err ) then
      errmax = err
    end if

    errl1 = errl1 + err * weight(ll)
    errl2 = errl2 + err**2 * weight(ll)

  end do

  errl1 = errl1 / totalw
  errl2 = sqrt ( errl2 / totalw )

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  Least square error = ', errl2
  write ( *, '(a,g14.6)' ) '  Average error      = ', errl1
  write ( *, '(a,g14.6)' ) '  Maximum error      = ', errmax
  write ( *, '(a)' ) ' '
  
  if ( iprfun /= 1 ) then
    return
  end if
!
!  Scale error curve and print.
!
  ie = 0
  scale = 1.0D+00

  if ( errmax < 10.0D+00 ) then
  
    do ie = 1, 9
      scale = scale * 10.0D+00
      if ( 10.0D+00 <= errmax * scale ) then
        exit
      end if
    end do

  end if  

  error(1:ntau) = error(1:ntau) * scale
  
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Approximation and scaled error curve'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i1)' ) &
    '       Data point       Approximation   Deviation x 10**', ie
  write ( *, '(a)' ) ' '
  write ( *, '(i4,f16.8,f16.8,f17.6)' ) &
    ( ll, tau(ll), ftau(ll), error(ll), ll = 1, ntau )

  return
end
subroutine l2knts ( break, l, k, t, n )

!*****************************************************************************80
!
!! L2KNTS converts breakpoints to knots.
!
!  Discussion:
!
!    The breakpoint sequence BREAK is converted into a corresponding 
!    knot sequence T to allow the representation of a piecewise
!    polynomial function of order K with K-2 continuous derivatives 
!    as a spline of order K with knot sequence T. 
!
!    This means that T(1:N+K) = BREAK(1) K times, then BREAK(2:L), 
!    then BREAK(L+1) K times.  
!
!    Therefore, N = K - 1 + L.
!
!  Modified:
!
!    14 February 2007
!
!  Author:
!
!    Carl DeBoor
!
!  Reference:
!
!    Carl DeBoor,
!    A Practical Guide to Splines,
!    Springer, 2001,
!    ISBN: 0387953663,
!    LC: QA1.A647.v27.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) K, the order.
!
!    Input, integer ( kind = 4 ) L, the number of polynomial pieces.
!
!    Input, real ( kind = 8 ) BREAK(L+1), the breakpoint sequence.
!
!    Output, real ( kind = 8 ) T(N+K), the knot sequence.
!
!    Output, integer ( kind = 4 ) N, the dimension of the corresponding spline 
!    space of order K.
!
  implicit none

  integer ( kind = 4 ) k
  integer ( kind = 4 ) l
  integer ( kind = 4 ) n

  real ( kind = 8 ) break(l+1)
  real ( kind = 8 ) t(k-1+l+k)

  n = k - 1 + l

  t(1:k-1) = break(1)
  t(k:n) = break(1:l) 
  t(n+1:n+k) = break(l+1)
 
  return
end
subroutine newnot ( break, coef, l, k, brknew, lnew, coefg )

!*****************************************************************************80
!
!! NEWNOT returns LNEW+1 knots which are equidistributed on (A,B).
!
!  Discussion:
!
!    The knots are equidistributed on (A,B) = ( BREAK(1), BREAK(L+1) ) 
!    with respect to a certain monotone function G related to the K-th root of 
!    the K-th derivative of the piecewise polynomial function F whose 
!    piecewise polynomial representation is contained in BREAK, COEF, L, K.
!
!  Method:
!
!    The K-th derivative of the given piecewise polynomial function F does 
!    not exist, except perhaps as a linear combination of delta functions. 
!
!    Nevertheless, we construct a piecewise constant function H with 
!    breakpoint sequence BREAK which is approximately proportional 
!    to abs ( d**K(F) ).
!
!    Specifically, on (BREAK(I), BREAK(I+1)),
!
!          abs(jump at BREAK(I) of PC)     abs(jump at BREAK(I+1) of PC)
!      H = ---------------------------  +  ----------------------------
!          BREAK(I+1) - BREAK(I-1)         BREAK(I+2) - BREAK(I)
!
!    with PC the piecewise constant (K-1)st derivative of F.
!
!    Then, the piecewise linear function G is constructed as
!
!      G(X) = integral ( A <= Y <= X )  H(Y)**(1/K) dY,
!
!    and its piecewise polynomial coefficients are stored in COEFG.
!
!    Then BRKNEW is determined by
!
!      BRKNEW(I) = A + G**(-1)((I-1)*STEP), for I = 1:LNEW+1,
!
!    where STEP = G(B) / LNEW and (A,B) = ( BREAK(1), BREAK(L+1) ).
!
!    In the event that PC = d**(K-1)(F) is constant in ( A, B ) and
!    therefore H = 0 identically, BRKNEW is chosen uniformly spaced.
!
!    If IPRINT is set positive, then the piecewise polynomial coefficients
!    of G will be printed out.
!
!  Modified:
!
!    14 February 2007
!
!  Author:
!
!    Carl DeBoor
!
!  Reference:
!
!    Carl DeBoor,
!    A Practical Guide to Splines,
!    Springer, 2001,
!    ISBN: 0387953663,
!    LC: QA1.A647.v27.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) BREAK(L+1), real ( kind = 8 ) COEF(K,L), 
!    integer ( kind = 4 ) L, integer K, the piecewise polynomial representation
!    of a certain function F of order K.  Specifically,
!      d**(k-1) F(X) = COEF(K,I) for BREAK(I) <= X < BREAK(I+1).
!
!    Input, integer ( kind = 4 ) LNEW, the number of intervals into which the 
!    interval (A,B) is to be divided by the new breakpoint sequence BRKNEW. 
!
!    Output, real ( kind = 8 ) BRKNEW(LNEW+1), the new breakpoint sequence.
!
!    Output, real ( kind = 8 ) COEFG(2,L), the coefficient part of the piecewise
!    polynomial representation BREAK, COEFG, L, 2 for the monotone piecewise 
!    linear function G with respect to which BRKNEW will be equidistributed.
!
  implicit none

  integer ( kind = 4 ) k
  integer ( kind = 4 ) l
  integer ( kind = 4 ) lnew

  real ( kind = 8 ) break(l+1)
  real ( kind = 8 ) brknew(lnew+1)
  real ( kind = 8 ) coef(k,l)
  real ( kind = 8 ) coefg(2,l)
  real ( kind = 8 ) dif
  real ( kind = 8 ) difprv
  integer ( kind = 4 ) i
  integer ( kind = 4 ), save :: iprint = 0
  integer ( kind = 4 ) j
  real ( kind = 8 ) oneovk
  real ( kind = 8 ) step
  real ( kind = 8 ) stepi
!
!  If G is constant, BRKNEW is uniform.
!
  if ( l <= 1 ) then
    call evnnot ( break, coef, l, k, brknew, lnew, coefg )
    return
  end if

  brknew(1) = break(1)
  brknew(lnew+1) = break(l+1)
!
!  Construct the continuous piecewise linear function G.
!
  oneovk = 1.0D+00 / real ( k, kind = 8 )
  coefg(1,1) = 0.0D+00
  difprv = abs ( coef(k,2) - coef(k,1) ) / ( break(3) - break(1) )
  
  do i = 2, l
    dif = abs ( coef(k,i) - coef(k,i-1) ) / ( break(i+1) - break(i-1) )
    coefg(2,i-1) = ( dif + difprv )**oneovk
    coefg(1,i) = coefg(1,i-1) + coefg(2,i-1) * ( break(i) - break(i-1) )
    difprv = dif
  end do
   
  coefg(2,l) = ( 2.0D+00 * difprv )**oneovk
!
!  STEP = G(B) / LNEW.
!
  step = ( coefg(1,l) + coefg(2,l) * ( break(l+1) - break(l) ) ) &
    / real ( lnew, kind = 8 )

  if ( 0 < iprint ) then
    write ( *, '(2x,e16.7)' ) step
    do i = 1, l
      write ( *, '(i5,2e16.5)' ) i, coefg(1:2,i)
    end do
  end if
!
!  If G is constant, BRKNEW is uniform.
!
  if ( step <= 0.0D+00 ) then
    call evnnot ( break, coef, l, k, brknew, lnew, coefg )
    return
  end if
!
!  For I = 2,..., LNEW, construct  BRKNEW(I) = A + G**(-1)(STEPI),
!  with STEPI = ( I - 1 ) * STEP.  
!
!  This requires inversion of the piecewise linear function G.  
!
!  For this, J is found so that
!
!    G(BREAK(J)) <= STEPI <= G(BREAK(J+1))
!
!  and then
!
!    BRKNEW(I) = BREAK(J) + ( STEPI - G(BREAK(J)) ) / DG(BREAK(J) ).
!
!  The midpoint is chosen if DG(BREAK(J)) = 0.
!
  j = 1
  
  do i = 2, lnew
  
    stepi = real ( i - 1, kind = 8 ) * step

    do

      if ( j == l ) then
        exit
      end if

      if ( stepi <= coefg(1,j+1) ) then
        exit
      end if

      j = j + 1

    end do
       
    if ( coefg(2,j) /= 0.0D+00 ) then
      brknew(i) = break(j) + ( stepi - coefg(1,j) ) / coefg(2,j)
    else
      brknew(i) = ( break(j) + break(j+1) ) / 2.0D+00
    end if
    
  end do
  
  return
end
function ppvalu ( break, coef, l, k, x, jderiv )

!*****************************************************************************80
!
!! PPVALU evaluates a piecewise polynomial function or its derivative.
!
!  Discussion:
!
!    PPVALU calculates the value at X of the JDERIV-th derivative of
!    the piecewise polynomial function F from its piecewise
!    polynomial representation.
!
!    The interval index I, appropriate for X, is found through a
!    call to INTERV.  The formula for the JDERIV-th derivative
!    of F is then evaluated by nested multiplication.
!
!    The J-th derivative of F is given by:
!
!      (d**J) F(X) = 
!        COEF(J+1,I) + H * (
!        COEF(J+2,I) + H * (
!        ...
!        COEF(K-1,I) + H * (
!        COEF(K,  I) / (K-J-1) ) / (K-J-2) ... ) / 2 ) / 1
!
!    with
!
!      H = X - BREAK(I)
!
!    and
!
!      I = max ( 1, max ( J, BREAK(J) <= X, 1 <= J <= L ) ).
!
!  Modified:
!
!    16 February 2007
!
!  Author:
!
!    Carl DeBoor
!
!  Reference:
!
!    Carl DeBoor,
!    A Practical Guide to Splines,
!    Springer, 2001,
!    ISBN: 0387953663,
!    LC: QA1.A647.v27.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) BREAK(L+1), real COEF(*), integer L, the
!    piecewise polynomial representation of the function F to be evaluated.
!
!    Input, integer ( kind = 4 ) K, the order of the polynomial pieces that 
!    make up the function F.  The usual value for K is 4, signifying a piecewise 
!    cubic polynomial.
!
!    Input, real ( kind = 8 ) X, the point at which to evaluate F or
!    of its derivatives.
!
!    Input, integer ( kind = 4 ) JDERIV, the order of the derivative to be
!    evaluated.  If JDERIV is 0, then F itself is evaluated,
!    which is actually the most common case.  It is assumed
!    that JDERIV is zero or positive.
!
!    Output, real ( kind = 8 ) PPVALU, the value of the JDERIV-th
!    derivative of F at X.
!
  implicit none

  integer ( kind = 4 ) k
  integer ( kind = 4 ) l

  real ( kind = 8 ) break(l+1)
  real ( kind = 8 ) coef(k,l)
  real ( kind = 8 ) fmmjdr
  real ( kind = 8 ) h
  integer ( kind = 4 ) i
  integer ( kind = 4 ) jderiv
  integer ( kind = 4 ) m
  integer ( kind = 4 ) ndummy
  real ( kind = 8 ) ppvalu
  real ( kind = 8 ) value
  real ( kind = 8 ) x

  value = 0.0D+00

  fmmjdr = k - jderiv
!
!  Derivatives of order K or higher are identically zero.
!
  if ( k <= jderiv ) then
    return
  end if
!
!  Find the index I of the largest breakpoint to the left of X.
!
  call interv ( break, l+1, x, i, ndummy )
!
!  Evaluate the JDERIV-th derivative of the I-th polynomial piece at X.
!
  h = x - break(i)
  m = k
 
  do

    value = ( value / fmmjdr ) * h + coef(m,i)
    m = m - 1
    fmmjdr = fmmjdr - 1.0D+00

    if ( fmmjdr <= 0.0D+00 ) then
      exit
    end if

  end do

  ppvalu = value
 
  return
end
subroutine putit ( t, kpm, left, scrtch, dbiatx, q, nrow, b )

!*****************************************************************************80
!
!! PUTIT puts together one block of the collocation equation system.
!
!  Method:
!
!    The K collocation equations for the interval ( T(LEFT), T(LEFT+1) )
!    are constructed with the aid of the subroutine DIFEQU( 2, ., . ) 
!    and interspersed (in order) with the side conditions, if any, in
!    this interval, using DIFEQU ( 3, ., . )  for the information.
!
!    The block Q has KPM columns, corresponding to the KPM B-splines of order 
!    KPM which have the interval ( T(LEFT), T(LEFT+1) ) in their support. 
!
!    The block's diagonal is part of the diagonal of the total system.
!
!    The first equation in this block not overlapped by the preceding block 
!    is therefore equation LOWROW, with LOWROW = number of side conditions 
!    in preceding intervals (or blocks).
!
!  Modified:
!
!    14 February 2007
!
!  Author:
!
!    Carl DeBoor
!
!  Reference:
!
!    Carl DeBoor,
!    A Practical Guide to Splines,
!    Springer, 2001,
!    ISBN: 0387953663,
!    LC: QA1.A647.v27.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) T(LEFT+KPM), the knot sequence.
!
!    Input, integer ( kind = 4 ) KPM, the order of the spline.
!
!    Input, integer ( kind = 4 ) LEFT, indicates the interval of interest, 
!    that is, the interval ( T(LEFT), T(LEFT+1) ).
!
!    Workspace, real ( kind = 8 ) SCRTCH(KPM,KPM).
!
!    Workspace, real ( kind = 8 ) DBIATX(KPM,M+1), derivatives of B-splines, 
!    with DBIATX(J,I+1) containing the I-th derivative of the J-th B-spline 
!    of interest.
!
!    Output, real ( kind = 8 ) Q(NROW,KPM), the block.
!
!    Input, integer ( kind = 4 ) NROW, number of rows in block to be put together.
!
!    Output, real ( kind = 8 ) B(NROW), the corresponding piece of 
!    the right hand side.
!
  implicit none

  integer ( kind = 4 ) kpm
  integer ( kind = 4 ) left
  integer ( kind = 4 ) nrow

  real ( kind = 8 ) b(nrow)
  real ( kind = 8 ) dbiatx(kpm,*)
  real ( kind = 8 ) dx
  integer ( kind = 4 ) i
  integer ( kind = 4 ) irow
  integer ( kind = 4 ) iside
  integer ( kind = 4 ) itermx
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) ll
  integer ( kind = 4 ) lowrow
  integer ( kind = 4 ) m
  integer ( kind = 4 ) mode
  integer ( kind = 4 ) mp1
  real ( kind = 8 ) q(nrow,kpm)
  real ( kind = 8 ) rho
  real ( kind = 8 ) scrtch(kpm,kpm)
  real ( kind = 8 ) sum1
  real ( kind = 8 ) t(left+kpm)
  real ( kind = 8 ) v(20)
  real ( kind = 8 ) xm
  real ( kind = 8 ) xside
  real ( kind = 8 ) xx

  save / other /
  save / side /

  common / other / itermx, k, rho(19)
  common / side / m, iside, xside(10)

  mp1 = m + 1
  
  q(1:nrow,1:kpm) = 0.0D+00
  
  xm = ( t(left+1) + t(left) ) / 2.0D+00
  dx = ( t(left+1) - t(left) ) / 2.0D+00

  ll = 1
  lowrow = iside

  do irow = lowrow, nrow

    if ( k < ll ) then
      go to 20
    end if

    mode = 2
!
!  Next collocation point:
!
    xx = xm + dx * rho(ll)
    ll = ll + 1
!
!  The corresponding collocation equation is next unless the next side
!  condition occurs at a point at, or to the left of, the next
!  collocation point.
!
    if ( m < iside ) then
      go to 30
    end if

    if ( xx < xside(iside) ) then
      go to 30
    end if

    ll = ll - 1

   20   continue

    mode = 3
    xx = xside(iside)

   30   continue

    call difequ ( mode, xx, v )
!
!  The next equation, a collocation equation (MODE=2) or a side
!  condition (MODE=3), reads
!
!    (*)   (V(M+1)*D**M+V(M)*D**(M-1) +...+ V(1)*D**0)F(XX) = V(M+2)
!
!  in terms of the information supplied by DIFEQU. 
!
!  The corresponding equation for the B-spline coefficients of F therefore
!  has the left side of (*), evaluated at each of the KPM B-splines having
!  XX in their support, as its KPM possibly nonzero coefficients.
!
    call bsplvd ( t, kpm, xx, left, scrtch, dbiatx, mp1 )
    
    do j = 1, kpm
    
      q(irow,j) = dot_product ( dbiatx(j,1:mp1), v(1:mp1) )
            
    end do
   
    b(irow) = v(m+2)
    
  end do
   
  return
end
subroutine r8vec_print ( n, a, title )

!*****************************************************************************80
!
!! R8VEC_PRINT prints an R8VEC.
!
!  Discussion:
!
!    An R8VEC is an array of double precision real values.
!
!  Modified:
!
!    22 August 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of components of the vector.
!
!    Input, real ( kind = 8 ) A(N), the vector to be printed.
!
!    Input, character ( len = * ) TITLE, an optional title.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n)
  integer ( kind = 4 ) i
  character ( len = * ) title

  if ( 0 < len_trim ( title ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) trim ( title )
  end if

  write ( *, '(a)' ) ' '
  do i = 1, n
    write ( *, '(2x,i8,2x,g16.8)' ) i, a(i)
  end do

  return
end
function round ( x, size )

!*****************************************************************************80
!
!! ROUND is called to add some noise to data.
!
!  Discussion:
!
!    This function simply adds plus or minus a perturbation value
!    to the input data.
!
!  Modified:
!
!    14 February 2007
!
!  Author:
!
!    Carl DeBoor
!
!  Reference:
!
!    Carl DeBoor,
!    A Practical Guide to Splines,
!    Springer, 2001,
!    ISBN: 0387953663,
!    LC: QA1.A647.v27.
!
!  Parameters:
!
!    Input. real ( kind = 8 ) X, the value to be perturbed.
!
!    Input, real ( kind = 8 ) SIZE, the size of the perturbation.
!
!    Output, real ( kind = 8 ) ROUND, the perturbed value.
!
  implicit none

  real ( kind = 8 ), save :: flip = -1.0D+00
  real ( kind = 8 ) round
  real ( kind = 8 ) size
  real ( kind = 8 ) x

  flip = -flip
  round = x + flip * size

  return
end
subroutine sbblok ( bloks, integs, nbloks, ipivot, b, x )

!*****************************************************************************80
!
!! SBBLOK solves a linear system that was factored by FCBLOK.
!
!  Discussion:
!
!    The routine supervises the solution, by forward and backward 
!    substitution, of the linear system 
!
!      A * x = b
!
!    for X, with the PLU factorization of A already generated in FCBLOK.
!    Individual blocks of equations are solved via SUBFOR and SUBBAK.
!
!  Modified:
!
!    14 February 2007
!
!  Author:
!
!    Carl DeBoor
!
!  Reference:
!
!    Carl DeBoor,
!    A Practical Guide to Splines,
!    Springer, 2001,
!    ISBN: 0387953663,
!    LC: QA1.A647.v27.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) BLOKS(*), integer INTEGS(3,NBLOKS), integer
!    NBLOKS, integer IPIVOT(*), are as on return from FCBLOK.
!
!    Input, real ( kind = 8 ) B(*), the right hand side, stored corresponding 
!    to the storage of the equations.  See comments in SLVBLK for details.
!
!    Output, real ( kind = 8 ) X(*), the solution vector.
!
  implicit none

  integer ( kind = 4 ) nbloks

  real ( kind = 8 ) b(*)
  real ( kind = 8 ) bloks(*)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) index
  integer ( kind = 4 ) indexb
  integer ( kind = 4 ) indexx
  integer ( kind = 4 ) integs(3,nbloks)
  integer ( kind = 4 ) ipivot(*)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) last
  integer ( kind = 4 ) nbp1
  integer ( kind = 4 ) ncol
  integer ( kind = 4 ) nrow
  real ( kind = 8 ) x(*)
!
!  Forward substitution:
!
  index = 1
  indexb = 1
  indexx = 1

  do i = 1, nbloks

    nrow = integs(1,i)
    last = integs(3,i)

    call subfor ( bloks(index), ipivot(indexb), nrow, last, b(indexb), &
      x(indexx) )

    index = nrow * integs(2,i) + index
    indexb = indexb + nrow
    indexx = indexx + last

  end do
!
!  Back substitution.
!
  nbp1 = nbloks + 1

  do j = 1, nbloks

    i = nbp1 - j
    nrow = integs(1,i)
    ncol = integs(2,i)
    last = integs(3,i)
    index = index - nrow * ncol
    indexb = indexb - nrow
    indexx = indexx - last

    call subbak ( bloks(index), ipivot(indexb), nrow, ncol, last, x(indexx) )

  end do
   
  return
end
subroutine setupq ( x, dx, y, npoint, v, qty )

!*****************************************************************************80
!
!! SETUPQ is to be called in SMOOTH.
!
!  Discussion:
!
!    Put DELX = X(*+1) - X(*) into V(*,4).
!
!    Put the three bands of Q' * D into V(*,1:3).
!
!    Put the three bands of ( D * Q )' * ( D * Q ) at and above the diagonal
!    into V(*,5:7).
!
!    Here, Q is the tridiagonal matrix of order ( NPOINT-2, NPOINT )
!    with general row 
!
!      1/DELX(I), -1/DELX(I)-1/DELX(I+1), 1/DELX(I+1)
!
!    and D is the diagonal matrix with general row DX(I).
!
!  Modified:
!
!    14 February 2007
!
!  Author:
!
!    Carl DeBoor
!
!  Reference:
!
!    Carl DeBoor,
!    A Practical Guide to Splines,
!    Springer, 2001,
!    ISBN: 0387953663,
!    LC: QA1.A647.v27.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X(NPOINT), the abscissas, assumed to be 
!    strictly increasing.
!
!    Input, real ( kind = 8 ) DX(NPOINT), the data uncertainty estimates,
!    which are assumed to be positive.
!
!    Input, real ( kind = 8 ) Y(NPOINT), the corresponding ordinates.
!
!    Input, integer ( kind = 4 ) NPOINT, the number of data points.
!
!    Output, real ( kind = 8 ) V(NPOINT,7), contains data needed for
!    the smoothing computation.
!
!    Output, real ( kind = 8 ) QTY(NPOINT), the value of Q' * Y.
!
  implicit none

  integer ( kind = 4 ) npoint

  real ( kind = 8 ) diff
  real ( kind = 8 ) dx(npoint)
  integer ( kind = 4 ) i
  real ( kind = 8 ) prev
  real ( kind = 8 ) qty(npoint)
  real ( kind = 8 ) v(npoint,7)
  real ( kind = 8 ) x(npoint)
  real ( kind = 8 ) y(npoint)

  v(1:npoint-1,4) = x(2:npoint) - x(1:npoint-1)

  v(2:npoint-1,1) = dx(1:npoint-2) / v(1:npoint-2,4)  
  v(npoint,1) = 0.0D+00

  v(2:npoint-1,2) = - dx(2:npoint-1) / v(2:npoint-1,4) &
                    - dx(2:npoint-1) / v(1:npoint-2,4)

  v(2:npoint-1,3) = dx(3:npoint) / v(2:npoint-1,4)
   
  v(2:npoint-1,5) = v(2:npoint-1,1)**2 &
                  + v(2:npoint-1,2)**2 &
                  + v(2:npoint-1,3)**2
   
  v(2:npoint-2,6) = v(2:npoint-2,2) * v(3:npoint-1,1) &
                  + v(2:npoint-2,3) * v(3:npoint-1,2)
  v(npoint-1,6) = 0.0D+00

  v(2:npoint-3,7) = v(2:npoint-3,3) * v(4:npoint-1,1)   
  v(npoint-2,7) = 0.0D+00
  v(npoint-1,7) = 0.0D+00
!
!  Construct Q' * Y in QTY.
!
  prev = ( y(2) - y(1) ) / v(1,4)
  do i = 2, npoint - 1
    diff = ( y(i+1) - y(i) ) / v(i,4)
    qty(i) = diff - prev
    prev = diff
  end do
  
  return
end
subroutine shiftb ( ai, ipivot, nrowi, ncoli, last, ai1, nrowi1, ncoli1 )

!*****************************************************************************80
!
!! SHIFTB shifts the rows in the current block.
!
!  Discussion:
!
!    This routine shifts rows in the current block, AI, which are not used 
!    as pivot rows, if any, that is, rows IPIVOT(LAST+1) through IPIVOT(NROWI), 
!    onto the first MMAX = NROW - LAST rows of the next block, AI1, 
!    with column LAST + J of AI going to column J, 
!    for J = 1,..., JMAX = NCOLI - LAST.
!
!    The remaining columns of these rows of AI1 are zeroed out.
!
!  Diagram:
!
!       Original situation after         Results in a new block I+1
!       LAST = 2 columns have been       created and ready to be
!       done in FACTRB, assuming no      factored by next FACTRB call.
!       interchanges of rows.
!
!                   1
!              X  X 1X  X  X           X  X  X  X  X
!                   1
!              0  X 1X  X  X           0  X  X  X  X
!  BLOCK I          1                       ---------------
!  NROWI=4     0  0 1X  X  X           0  0 1X  X  X  0  01
!  NCOLI=5          1                       1             1
!  LAST=2      0  0 1X  X  X           0  0 1X  X  X  0  01
!              -------------------          1             1   NEW
!                   1X  X  X  X  X          1X  X  X  X  X1  BLOCK
!                   1                       1             1   I+1
!  BLOCK I+1        1X  X  X  X  X          1X  X  X  X  X1
!  NROWI1= 5        1                       1             1
!  NCOLI1= 5        1X  X  X  X  X          1X  X  X  X  X1
!              -------------------          1-------------1
!                   1
!
!  Modified:
!
!    14 February 2007
!
!  Author:
!
!    Carl DeBoor
!
!  Reference:
!
!    Carl DeBoor,
!    A Practical Guide to Splines,
!    Springer, 2001,
!    ISBN: 0387953663,
!    LC: QA1.A647.v27.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) AI(NROWI,NCOLI), the current block.
!
!    Input, integer ( kind = 4 ) IPIVOT(NROWI), the pivot vector.
!
!    Input, integer ( kind = 4 ) NROWI, NCOLI, the number of rows and columns
!    in block AI.
!
!    Input, integer ( kind = 4 ) LAST, indicates the last row on which pivoting
!    has been carried out.
!
!    Input/output, real ( kind = 8 ) AI1(NROWI1,NCOLI1), the next block.
!
!    Input, integer ( kind = 4 ) NROWI1, NCOLI1, the number of rows and columns
!    in block AI1.
!
  implicit none

  integer ( kind = 4 ) ncoli
  integer ( kind = 4 ) ncoli1
  integer ( kind = 4 ) nrowi1
  integer ( kind = 4 ) nrowi

  real ( kind = 8 ) ai(nrowi,ncoli)
  real ( kind = 8 ) ai1(nrowi1,ncoli1)
  integer ( kind = 4 ) ip
  integer ( kind = 4 ) ipivot(nrowi)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) last
  integer ( kind = 4 ) m

  if ( nrowi - last < 1 ) then
    return
  end if

  if ( ncoli - last < 1 ) then
    return
  end if
!
!  Put the remainder of block I into AI1.
!
  do m = 1, nrowi - last
    ip = ipivot(last+m)
    do j = 1, ncoli - last
      ai1(m,j) = ai(ip,last+j)
    end do
  end do
!
!  Zero out the upper right corner of AI1.
!
  do j = ncoli + 1 - last, ncoli1
    do m = 1, nrowi - last
      ai1(m,j) = 0.0D+00
    end do
  end do
  
  return
end
subroutine slvblk ( bloks, integs, nbloks, b, ipivot, x, iflag )

!*****************************************************************************80
!
!! SLVBLK solves the almost block diagonal linear system A * x = b.  
!
!  Discussion:
!
!    Such almost block diagonal matrices arise naturally in piecewise 
!    polynomial interpolation or approximation and in finite element 
!    methods for two-point boundary value problems.  The PLU factorization 
!    method is implemented here to take advantage of the special structure 
!    of such systems for savings in computing time and storage requirements.
!
!    SLVBLK relies on several auxiliary programs:
!
!    FCBLOK (BLOKS,INTEGS,NBLOKS,IPIVOT,SCRTCH,IFLAG)  
!    factors the matrix A.
!
!    SBBLOK (BLOKS,INTEGS,NBLOKS,IPIVOT,B,X)
!    solves the system A*X=B once A is factored.
!
!    DTBLOK (BLOKS,INTEGS,NBLOKS,IPIVOT,IFLAG,DETSGN,DETLOG) 
!    computes the determinant of A once it has been factored.
!
!  Block structure of A:
!
!    The NBLOKS blocks are stored consecutively in the array BLOKS.
!
!    The first block has its (1,1)-entry at BLOKS(1), and, if the I-th
!    block has its (1,1)-entry at BLOKS(INDEX(I)), then
!
!      INDEX(I+1) = INDEX(I) + NROW(I) * NCOL(I).
!
!    The blocks are pieced together to give the interesting part of A
!    as follows.  For I=1,2,..., NBLOKS-1, the (1,1)-entry of the next
!    block (the (I+1)st block) corresponds to the (LAST+1,LAST+1)-entry
!    of the current I-th block.  Recall LAST = INTEGS(3,I) and note that
!    this means that
!
!    A: every block starts on the diagonal of A.
!
!    B: the blocks overlap (usually). the rows of the (I+1)st block
!       which are overlapped by the I-th block may be arbitrarily 
!       defined initially.  They are overwritten during elimination.
!
!    The right hand side for the equations in the I-th block are stored
!    correspondingly as the last entries of a piece of B of length NROW
!    (= INTEGS(1,I)) and following immediately in B the corresponding
!    piece for the right hand side of the preceding block, with the right 
!    hand side for the first block starting at B(1).  In this, the right 
!    hand side for an equation need only be specified once on input, 
!    in the first block in which the equation appears.
!
!  Example:
!
!    The test driver for this package contains an example, a linear
!    system of order 11, whose nonzero entries are indicated in the
!    following diagram by their row and column index modulo 10.  Next to it
!    are the contents of the INTEGS arrray when the matrix is taken to
!    be almost block diagonal with NBLOKS = 5, and below it are the five
!    blocks.
!
!                        NROW1 = 3, NCOL1 = 4
!             11 12 13 14
!             21 22 23 24   NROW2 = 3, NCOL2 = 3
!             31 32 33 34
!    LAST1 = 2      43 44 45
!                   53 54 55            NROW3 = 3, NCOL3 = 4
!          LAST2 = 3         66 67 68 69   NROW4 = 3, NCOL4 = 4
!                            76 77 78 79      NROW5 = 4, NCOL5 = 4
!                            86 87 88 89
!                   LAST3 = 1   97 98 99 90
!                      LAST4 = 1   08 09 00 01
!                                  18 19 10 11
!                         LAST5 = 4
!
!    Actual input to BLOKS shown by rows of blocks of A.
!    The ** items are arbitrary.
!
!    11 12 13 14  / ** ** **  / 66 67 68 69  / ** ** ** **  / ** ** ** **
!    21 22 23 24 /  43 44 45 /  76 77 78 79 /  ** ** ** ** /  ** ** ** **
!    31 32 33 34/   53 54 55/   86 87 88 89/   97 98 99 90/   08 09 00 01
!                                                             18 19 10 11
!
!    INDEX = 1      INDEX = 13  INDEX = 22     INDEX = 34     INDEX = 46
!
!    Actual right hand side values with ** for arbitrary values:
!
!      B1 B2 B3 ** B4 B5 B6 B7 B8 ** ** B9 ** ** B10 B11
!
!    It would have been more efficient to combine block 3 with block 4.
!
!  Modified:
!
!    14 February 2007
!
!  Author:
!
!    Carl DeBoor
!
!  Reference:
!
!    Carl DeBoor,
!    A Practical Guide to Splines,
!    Springer, 2001,
!    ISBN: 0387953663,
!    LC: QA1.A647.v27.
!
!  Parameters:
!
!    Input/output, real ( kind = 8 ) BLOKS(*), a one-dimenional array, 
!    of length sum ( INTEGS(1,1:NBLOKS) * INTEGS(2,1:NBLOKS) ).
!    On input, contains the blocks of the almost block diagonal matrix A.  
!    The array INTEGS describes the block structure.
!    On output, contains correspondingly the PLU factorization
!    of A, if IFLAG /= 0.  Certain entries in BLOKS are arbitrary, 
!    where the blocks overlap.
!
!    Input, integer ( kind = 4 ) INTEGS(3,NBLOKS), description of the block 
!    structure of A.
!    integs(1,I) = number of rows of block I = nrow;
!    integs(2,I) = number of colums of block I = ncol;
!    integs(3,I) = number of elimination steps in block I = last.
!    The linear system is of order n = sum ( integs(3,i), i=1,...,nbloks ),
!    but the total number of rows in the blocks is
!    nbrows=sum( integs(1,i) ; i = 1,...,nbloks)
!
!    Input, integer ( kind = 4 ) NBLOKS, the number of blocks.
!
!    Input, real ( kind = 8 ) B(NBROWS), the right hand side.  Certain entries 
!    are arbitrary, corresponding to rows of the blocks which overlap.  See 
!    the block structure in the example.
!
!    Output, integer ( kind = 4 ) IPIVOT(NBROWS), the pivoting sequence used.
!
!    Output, real ( kind = 8 ) X(N), the computed solution, if iflag /= 0.
!
!    Output, integer ( kind = 4 ) IFLAG.
!    = (-1)**(number of interchanges during factorization) if A is invertible;
!    = 0 if A is singular.
!
  implicit none

  integer ( kind = 4 ) nbloks

  real ( kind = 8 ) b(*)
  real ( kind = 8 ) bloks(*)
  integer ( kind = 4 ) iflag
  integer ( kind = 4 ) integs(3,nbloks)
  integer ( kind = 4 ) ipivot(*)
  real ( kind = 8 ) x(*)
!
!  In the call to FCBLOK, X is used for temporary storage.
!
  call fcblok ( bloks, integs, nbloks, ipivot, x, iflag )
  
  if ( iflag == 0 ) then
    return
  end if
  
  call sbblok ( bloks, integs, nbloks, ipivot, b, x )
  
  return
end
function smooth ( x, y, dy, npoint, s, v, a )

!*****************************************************************************80
!
!! SMOOTH constructs the cubic smoothing spline to given data.
!
!  Discussion:
!
!    The data is of the form
!
!      ( X(1:NPOINT), Y(1:NPOINT) )
!
!    The cubic smoothing spline has as small a second derivative as
!    possible, while
!
!      S(F) <= S,
!
!    where
!
!      S(F) = sum ( 1 <= I <= NPOINT ) ( ( ( Y(I) - F(X(I)) ) / DY(I) )**2.
!
!  Method:
!
!    The matrices Q' * D and Q' * D**2 * Q are constructed in SETUPQ from
!    X and DY, as is the vector QTY = Q' * Y.
!
!    Then, for given P, the vector U is determined in CHOL1D as
!    the solution of the linear system
!
!      ( 6 * (1-P) * Q' * D**2 * Q + P * R ) * U = QTY.
!
!    From U and this choice of smoothing parameter P, the smoothing spline F
!    is obtained in the sense that:
!
!             F(X(.)) = Y - 6 (1-P) D**2 * Q * U,
!      (d**2) F(X(.)) = 6 * P * U.
!
!    The smoothing parameter P is found, if possible, so that
!
!      SF(P) = S,
!
!    with SF(P) = S(F), where F is the smoothing spline as it depends
!    on P.  If S = 0, then P = 1.  If SF(0) <= S, then P = 0.
!    Otherwise, the secant method is used to locate an appropriate P in
!    the open interval (0,1).
!
!    Specifically,
!
!      P(0) = 0,  P(1) = ( S - SF(0) ) / DSF
!
!    with
!
!      DSF = -24 * U' * R * U
!
!    a good approximation to
!
!      D(SF(0)) = DSF + 60 * (D*Q*U)' * (D*Q*U),
!
!    and U as obtained for P = 0.
!
!    After that, for N = 1, 2,...  until SF(P(N)) <= 1.01 * S, do:
!    determine P(N+1) as the point at which the secant to SF at the
!    points P(N) and P(N-1) takes on the value S.
!
!    If 1 <= P(N+1), choose instead P(N+1) as the point at which
!    the parabola SF(P(N))*((1-.)/(1-P(N)))**2 takes on the value S.
!
!    Note that, in exact arithmetic, it is always the case that
!      P(N+1) < P(N),
!    hence
!      SF(P(N+1)) < SF(P(N)).
!
!    Therefore, also stop the iteration, with final P = 1, in case
!      SF(P(N)) <= SF(P(N+1)).
!
!  Modified:
!
!    16 February 2007
!
!  Author:
!
!    Carl DeBoor
!
!  Reference:
!
!    Carl DeBoor,
!    A Practical Guide to Splines,
!    Springer, 2001,
!    ISBN: 0387953663,
!    LC: QA1.A647.v27.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X(NPOINT), the abscissas, assumed to be
!    strictly increasing.
!
!    Input, real ( kind = 8 ) Y(NPOINT), the corresponding ordinates.
!
!    Input, real ( kind = 8 ) DY(NPOINT), the data uncertainty estimates,
!    which are assumed to be positive.
!
!    Input, integer ( kind = 4 ) NPOINT, the number of data points.
!
!    Input, real ( kind = 8 ) S, an upper bound on the discrete weighted mean
!    square distance of the approximation F from the data.
!
!    Workspace, real ( kind = 8 ) V(NPOINT,7).
!
!    Workspace, real ( kind = 8 ) A(NPOINT,4).
!
!    Output, real ( kind = 8 ) A(NPOINT,4).
!    A(*,1).....contains the sequence of smoothed ordinates.
!    A(I,J) = d**(J-1) F(X(I)), for J = 2:4, I = 1:NPOINT-1.
!    That is, the first three derivatives of the smoothing spline F at the
!    left end of each of the data intervals.  Note that A would have to
!    be transposed before it could be used in PPVALU.
!
!    Output, real ( kind = 8 ) SMOOTH, the value of the smoothing parameter.
!
  implicit none

  integer ( kind = 4 ) npoint

  real ( kind = 8 ) a(npoint,4)
  real ( kind = 8 ) change
  real ( kind = 8 ) dy(npoint)
  integer ( kind = 4 ) i
  real ( kind = 8 ) oosf
  real ( kind = 8 ) ooss
  real ( kind = 8 ) p
  real ( kind = 8 ) prevq
  real ( kind = 8 ) prevsf
  real ( kind = 8 ) q
  real ( kind = 8 ) s
  real ( kind = 8 ) sfq
  real ( kind = 8 ) smooth
  real ( kind = 8 ) utru
  real ( kind = 8 ) v(npoint,7)
  real ( kind = 8 ) x(npoint)
  real ( kind = 8 ) y(npoint)

  call setupq ( x, dy, y, npoint, v, a(1,4) )

  if ( s <= 0.0D+00 ) then

    p = 1.0D+00
    call chol1d ( p, v, a(1,4), npoint, 1, a(1,3), a(1,1) )
    sfq = 0.0D+00

  else

    p = 0.0D+00
    call chol1d ( p, v, a(1,4), npoint, 1, a(1,3), a(1,1) )

    sfq = 36.0D+00 * dot_product ( a(1:npoint,1)**2, dy(1:npoint)**2 )

    if ( s < sfq ) then

      utru = 0.0D+00
      do i = 2, npoint
        utru = utru + v(i-1,4) * ( a(i-1,3) * ( a(i-1,3) + a(i,3) ) &
          + a(i,3)**2 )
      end do

      ooss = 1.0D+00 / sqrt ( s )
      oosf = 1.0D+00 / sqrt ( sfq )
      q = - ( oosf - ooss ) * sfq / ( 6.0D+00 * utru * oosf )
!
!  Secant iteration for the determination of P starts here.
!
      prevq = 0.0D+00
      prevsf = oosf

      do

        call chol1d ( q / ( 1.0D+00 + q ), v, a(1,4), npoint, 1, &
          a(1,3), a(1,1) )

        sfq = 36.0D+00 * dot_product ( a(1:npoint,1)**2, dy(1:npoint)**2 ) &
          / ( 1.0D+00 + q )**2

        if ( abs ( sfq - s ) <= 0.01D+00 * s ) then
          exit
        end if

        oosf = 1.0D+00 / sqrt ( sfq )
        change = ( q - prevq ) / ( oosf - prevsf ) * ( oosf - ooss )
        prevq = q
        q = q - change
        prevsf = oosf

      end do

      p = q / ( 1.0D+00 + q )

    end if

  end if
!
!  Correct value of P has been found.
!  Compute polynomial coefficients from Q * U in A(.,1).
!
  smooth = sfq

  a(1:npoint,1) = y(1:npoint) - 6.0D+00 * ( 1.0D+00 - p ) &
    * dy(1:npoint)**2 * a(1:npoint,1)

  a(1:npoint,3) = a(1:npoint,3) * 6.0D+00 * p

  do i = 1, npoint - 1
    a(i,4) = ( a(i+1,3) - a(i,3) ) / v(i,4)
    a(i,2) = ( a(i+1,1) - a(i,1) ) / v(i,4) &
      - ( a(i,3) + a(i,4) / 3.0D+00 * v(i,4) ) / 2.0D+00 * v(i,4)
  end do

  return
end
subroutine spli2d ( tau, gtau, t, n, k, m, work, q, bcoef, iflag )

!*****************************************************************************80
!
!! SPLI2D produces a interpolatory tensor product spline.
!
!  Discussion:
!
!    SPLI2D is an extended version of SPLINT.
!
!    SPLI2D produces the B-spline coefficients BCOEF(J,.) of the 
!    spline of order K with knots T(1:N+K), which takes on 
!    the value GTAU(I,J) at TAU(I), I=1,..., N, J=1,...,M.
! 
!    The I-th equation of the linear system 
!
!      A * BCOEF = B  
!  
!    for the B-spline coefficients of the interpolant enforces 
!    interpolation at TAU(I), I=1,...,N.  Hence,  B(I) = GTAU(I), 
!    for all I, and A is a band matrix with 2*K-1 bands, if it is 
!    invertible.
! 
!    The matrix A is generated row by row and stored, diagonal by
!    diagonal, in the rows of the array Q, with the main diagonal
!    going into row K.
! 
!    The banded system is then solved by a call to BANFAC, which 
!    constructs the triangular factorization for A and stores it 
!    again in Q, followed by a call to BANSLV, which then obtains 
!    the solution BCOEF by substitution.
!
!     The linear system to be solved is theoretically invertible if
!     and only if 
!       
!       T(I) < TAU(I) < TAU(I+K), for all I.
!         
!     Violation of this condition is certain to lead to IFLAG = 2.
! 
!  Modified:
!
!    14 February 2007
!
!  Author:
!
!    Carl DeBoor
!
!  Reference:
!
!    Carl DeBoor,
!    A Practical Guide to Splines,
!    Springer, 2001,
!    ISBN: 0387953663,
!    LC: QA1.A647.v27.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) TAU(N), contains the data point abscissas.
!    TAU must be strictly increasing
! 
!    Input, real ( kind = 8 ) GTAU(N,M), contains the data point ordinates.
! 
!    Input, real ( kind = 8 ) T(N+K), the knot sequence.
! 
!    Input, integer ( kind = 4 ) N, the number of data points and the 
!    dimension of the spline space SPLINE(K,T)
! 
!    Input, integer ( kind = 4 ) K, the order of the spline.
! 
!    Input, integer ( kind = 4 ) M, the number of data sets.
!
!    Work space, real ( kind = 8 ) WORK(N).
! 
!    Output, real ( kind = 8 ) Q(2*K-1)*N, the triangular 
!    factorization of the coefficient matrix of the linear 
!    system for the B-spline coefficients of the spline interpolant.
!    The B-spline coefficients for the interpolant of an additional 
!    data set ( TAU(I), HTAU(I) ), I=1,...,N  with the same data 
!    abscissae can be obtained without going through all the 
!    calculations in this routine, simply by loading HTAU into 
!    BCOEF and then using the statement
!      CALL BANSLV ( Q, 2*K-1, N, K-1, K-1, BCOEF )
! 
!    Output, real ( kind = 8 ) BCOEF(N), the B-spline coefficients of
!    the interpolant.
! 
!    Output, integer ( kind = 4 ) IFLAG, error indicator.
!    1, no error.   
!    2, an error occurred, which may have been caused by 
!       singularity of the linear system.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) bcoef(m,n)
  real ( kind = 8 ) gtau(n,m)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) iflag
  integer ( kind = 4 ) ilp1mx
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jj
  integer ( kind = 4 ) k
  integer ( kind = 4 ) left
  real ( kind = 8 ) q((2*k-1)*n)
  real ( kind = 8 ) t(n+k)
  real ( kind = 8 ) tau(n)
  real ( kind = 8 ) taui
  real ( kind = 8 ) work(n)

  left = k
  
  q(1:(2*k-1)*n) = 0.0D+00
!
!  Construct the N interpolation equations.
!
  do i = 1, n
  
    taui = tau(i)
    ilp1mx = min ( i + k, n + 1 )
!
!  Find the index LEFT in the closed interval (I,I+K-1) such that:
!
!    T(LEFT) < = TAU(I) < T(LEFT+1)
!
!  The matrix will be singular if this is not possible.
!
    left = max ( left, i )
    
    if ( taui < t(left) ) then
      iflag = 2
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'SPLI2D - Fatal error!'
      write ( *, '(a)' ) '  The TAU array is not strictly increasing.'
      stop
    end if
   
    do while ( t(left+1) <= taui )
   
      left = left + 1

      if ( left < ilp1mx ) then
        cycle
      end if
    
      left = left - 1
    
      if ( t(left+1) < taui ) then
        iflag = 2
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'SPLI2D - Fatal error!'
        write ( *, '(a)' ) '  The TAU array is not strictly increasing.'
        stop
      end if
 
      exit

    end do
!
!  The I-th equation enforces interpolation at TAUI, hence
!
!    A(I,J) = B(J,K,T)(TAUI), for all J. 
!
!  Only the K entries with J = LEFT-K+1, ..., LEFT actually might be 
!  nonzero.  These K numbers are returned, in WORK (used for 
!  temporary storage here), by the following call:
!
    call bsplvb ( t, k, 1, taui, left, work )
!
!  We therefore want  
!
!    WORK(J) = B(LEFT-K+J)(TAUI) 
!
!  to go into
!
!    A(I,LEFT-K+J),
!
!  that is, into  Q(I-(LEFT+J)+2*K,(LEFT+J)-K) since
!  A(I+J,J) is to go into Q(I+K,J), for all I, J, if we consider Q
!  as a two-dimensional array, with  2*K-1 rows.  See comments in
!  BANFAC.
!
!  In the present program, we treat Q as an equivalent one-dimensional 
!  array, because of fortran restrictions on dimension statements.  
!
!  We therefore want WORK(J) to go into the entry of Q with index:
!    I -(LEFT+J)+2*K + ((LEFT+J)-K-1)*(2*K-1)
!    = I-LEFT+1+(LEFT -K)*(2*K-1) + (2*K-2)*J
!
    jj = i - left + 1 + ( left - k ) * ( k + k - 1 )
    
    do j = 1, k
      jj = jj + k + k - 2
      q(jj) = work(j)
    end do
    
  end do
!
!  Factor A, stored again in Q.
!
  call banfac ( q, k+k-1, n, k-1, k-1, iflag )
  
  if ( iflag == 2 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'SPLI2D - Fatal error!'
    write ( *, '(a)' ) '  BANFAC reports that the matrix is singular.'
    stop
  end if
!
!  Solve 
!
!    A * BCOEF = GTAU 
!
!  by back substitution.
!
  do j = 1, m
   
    work(1:n) = gtau(1:n,j)
   
    call banslv ( q, k+k-1, n, k-1, k-1, work )
    
    bcoef(j,1:n) = work(1:n)
    
  end do
   
  return
end
subroutine splint ( tau, gtau, t, n, k, q, bcoef, iflag )

!*****************************************************************************80
!
!! SPLINT produces the B-spline coefficients BCOEF of an interpolating spline.
!
!  Discussion:
!
!    The spline is of order K with knots T(1:N+K), and takes on the 
!    value GTAU(I) at TAU(I), for I = 1 to N.
!
!    The I-th equation of the linear system 
!
!      A * BCOEF = B 
!
!    for the B-spline coefficients of the interpolant enforces interpolation
!    at TAU(1:N).
!
!    Hence, B(I) = GTAU(I), for all I, and A is a band matrix with 2*K-1
!    bands, if it is invertible.
!
!    The matrix A is generated row by row and stored, diagonal by diagonal,
!    in the rows of the array Q, with the main diagonal going
!    into row K.  See comments in the program.
!
!    The banded system is then solved by a call to BANFAC, which 
!    constructs the triangular factorization for A and stores it again in
!    Q, followed by a call to BANSLV, which then obtains the solution
!    BCOEF by substitution.
!
!    BANFAC does no pivoting, since the total positivity of the matrix
!    A makes this unnecessary.
!
!    The linear system to be solved is (theoretically) invertible if
!    and only if
!      T(I) < TAU(I) < TAU(I+K), for all I.
!    Violation of this condition is certain to lead to IFLAG = 2.
!
!  Modified:
!
!    14 February 2007
!
!  Author:
!
!    Carl DeBoor
!
!  Reference:
!
!    Carl DeBoor,
!    A Practical Guide to Splines,
!    Springer, 2001,
!    ISBN: 0387953663,
!    LC: QA1.A647.v27.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) TAU(N), the data point abscissas.  The entries in
!    TAU should be strictly increasing.
!
!    Input, real ( kind = 8 ) GTAU(N), the data ordinates.
!
!    Input, real ( kind = 8 ) T(N+K), the knot sequence.
!
!    Input, integer ( kind = 4 ) N, the number of data points.
!
!    Input, integer ( kind = 4 ) K, the order of the spline.
!
!    Output, real ( kind = 8 ) Q((2*K-1)*N), the triangular factorization
!    of the coefficient matrix of the linear system for the B-coefficients 
!    of the spline interpolant.  The B-coefficients for the interpolant 
!    of an additional data set can be obtained without going through all 
!    the calculations in this routine, simply by loading HTAU into BCOEF 
!    and then executing the call:
!      call banslv ( q, 2*k-1, n, k-1, k-1, bcoef )
!
!    Output, real ( kind = 8 ) BCOEF(N), the B-spline coefficients of 
!    the interpolant.
!
!    Output, integer ( kind = 4 ) IFLAG, error flag.
!    1, = success.
!    2, = failure.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) bcoef(n)
  real ( kind = 8 ) gtau(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) iflag
  integer ( kind = 4 ) ilp1mx
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jj
  integer ( kind = 4 ) k
  integer ( kind = 4 ) kpkm2
  integer ( kind = 4 ) left
  real ( kind = 8 ) q((2*k-1)*n)
  real ( kind = 8 ) t(n+k)
  real ( kind = 8 ) tau(n)
  real ( kind = 8 ) taui

  kpkm2 = 2 * ( k - 1 )
  left = k
  q(1:(2*k-1)*n) = 0.0D+00
!
!  Loop over I to construct the N interpolation equations.
!
  do i = 1, n
  
    taui = tau(i)
    ilp1mx = min ( i + k, n + 1 )
!
!  Find LEFT in the closed interval (I,I+K-1) such that
!
!    T(LEFT) <= TAU(I) < T(LEFT+1)
!
!  The matrix is singular if this is not possible.
!
    left = max ( left, i )

    if ( taui < t(left) ) then
      iflag = 2
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'SPLINT - Fatal Error!'
      write ( *, '(a)' ) '  The linear system is not invertible!'
      return
    end if

    do while ( t(left+1) <= taui )

      left = left + 1

      if ( left < ilp1mx ) then
        cycle
      end if

      left = left - 1

      if ( t(left+1) < taui ) then
        iflag = 2
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'SPLINT - Fatal Error!'
        write ( *, '(a)' ) '  The linear system is not invertible!'
        return
      end if

      exit

    end do
!
!  The I-th equation enforces interpolation at TAUI, hence for all J,
!
!    A(I,J) = B(J,K,T)(TAUI).
!
!  Only the K entries with J = LEFT-K+1,...,LEFT actually might be nonzero.
!
!  These K numbers are returned, in BCOEF (used for temporary storage here),
!  by the following.
!
    call bsplvb ( t, k, 1, taui, left, bcoef )
!
!  We therefore want BCOEF(J) = B(LEFT-K+J)(TAUI) to go into
!  A(I,LEFT-K+J), that is, into Q(I-(LEFT+J)+2*K,(LEFT+J)-K) since
!  A(I+J,J) is to go into Q(I+K,J), for all I, J, if we consider Q
!  as a two-dimensional array, with  2*K-1 rows.  See comments in
!  BANFAC.
!
!  In the present program, we treat Q as an equivalent
!  one-dimensional array, because of fortran restrictions on
!  dimension statements.
!
!  We therefore want  BCOEF(J) to go into the entry of Q with index:
!
!    I -(LEFT+J)+2*K + ((LEFT+J)-K-1)*(2*K-1)
!   = I-LEFT+1+(LEFT -K)*(2*K-1) + (2*K-2)*J
!
    jj = i - left + 1 + ( left - k ) * ( k + k - 1 )

    do j = 1, k
      jj = jj + kpkm2
      q(jj) = bcoef(j)
    end do
    
  end do
!
!  Obtain factorization of A, stored again in Q.
!
  call banfac ( q, k+k-1, n, k-1, k-1, iflag )
  
  if ( iflag == 2 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'SPLINT - Fatal Error!'
    write ( *, '(a)' ) '  The linear system is not invertible!'
    return
  end if
!
!  Solve 
!
!    A * BCOEF = GTAU
!
!  by back substitution.
!
  bcoef(1:n) = gtau(1:n)

  call banslv ( q, k+k-1, n, k-1, k-1, bcoef )

  return
end
subroutine splopt ( tau, n, k, scrtch, t, iflag )

!*****************************************************************************80
!
!! SPLOPT computes the knots for an optimal recovery scheme. 
!
!  Discussion:
!
!    The optimal recovery scheme is of order K for data at TAU(1:N).
!
!    The interior knots T(K+1:N) are determined by Newton's method in 
!    such a way that the signum function which changes sign at
!      T(K+1:N)  and nowhere else in ( TAU(1), TAU(N) ) is 
!    orthogonal to the spline space SPLINE ( K, TAU ) on that interval.
!
!    Let XI(J) be the current guess for T(K+J), J=1,...,N-K.  Then
!    the next Newton iterate is of the form
!
!      XI(J) + (-)**(N-K-J)*X(J),  J=1,...,N-K,
!
!    with X the solution of the linear system
!
!      C * X = D.
!
!    Here, for all J,
!
!      C(I,J) = B(I)(XI(J)), 
!
!    with B(I) the I-th B-spline of order K for the knot sequence TAU, 
!    for all I, and D is the vector given, for each I, by
!
!      D(I) = sum ( -A(J), J=I,...,N ) * ( TAU(I+K) - TAU(I) ) / K,
!
!    with, for I = 1 to N-1:  
!
!      A(I) = sum ( (-)**(N-K-J)*B(I,K+1,TAU)(XI(J)), J=1,...,N-K )
!
!    and  
!
!      A(N) = -0.5.
!
!    See Chapter XIII of text and references there for a derivation.
!
!    The first guess for T(K+J) is sum ( TAU(J+1:J+K-1) ) / ( K - 1 ).
!
!    The iteration terminates if max ( abs ( X(J) ) ) < TOL, with
!
!      TOL = TOLRTE * ( TAU(N) - TAU(1) ) / ( N - K ),
!
!    or else after NEWTMX iterations.
!
!  Modified:
!
!    14 February 2007
!
!  Author:
!
!    Carl DeBoor
!
!  Reference:
!
!    Carl DeBoor,
!    A Practical Guide to Splines,
!    Springer, 2001,
!    ISBN: 0387953663,
!    LC: QA1.A647.v27.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) TAU(N), the interpolation points.
!    assumed to be nondecreasing, with TAU(I) < TAU(I+K), for all I.
!
!    Input, integer ( kind = 4 ) N, the number of data points.
!
!    Input, integer ( kind = 4 ) K, the order of the optimal recovery scheme 
!    to be used.
!
!    Workspace, real ( kind = 8 ) SCRTCH((N-K)*(2*K+3)+5*K+3).  The various
!    contents are specified in the text below.
!
!    Output, real ( kind = 8 ) T(N+K), the optimal knots ready for
!    use in optimal recovery.  Specifically, T(1:K) = TAU(1), 
!    T(N+1:N+K) = TAU(N), while the N - K interior knots T(K+1:N) 
!    are calculated.
!
!    Output, integer ( kind = 4 ) IFLAG, error indicator.
!    = 1, success.  T contains the optimal knots.
!    = 2, failure.  K < 3 or N < K or the linear system was singular.
!
  implicit none

  integer ( kind = 4 ) k
  integer ( kind = 4 ) n

  real ( kind = 8 ) del
  real ( kind = 8 ) delmax
  real ( kind = 8 ) floatk
  integer ( kind = 4 ) i
  integer ( kind = 4 ) id
  integer ( kind = 4 ) iflag
  integer ( kind = 4 ) index
  integer ( kind = 4 ) j
  integer ( kind = 4 ) kp1
  integer ( kind = 4 ) kpkm1
  integer ( kind = 4 ) kpn
  integer ( kind = 4 ) l
  integer ( kind = 4 ) left
  integer ( kind = 4 ) leftmk
  integer ( kind = 4 ) lenw
  integer ( kind = 4 ) ll
  integer ( kind = 4 ) llmax
  integer ( kind = 4 ) llmin
  integer ( kind = 4 ) na
  integer ( kind = 4 ) nb
  integer ( kind = 4 ) nc
  integer ( kind = 4 ) nd
  integer ( kind = 4 ), parameter :: newtmx = 10
  integer ( kind = 4 ) newton
  integer ( kind = 4 ) nmk
  integer ( kind = 4 ) nx
  real ( kind = 8 ) scrtch((n-k)*(2*k+3)+5*k+3)
  real ( kind = 8 ) t(n+k)
  real ( kind = 8 ) tau(n)
  real ( kind = 8 ) sign
  real ( kind = 8 ) signst
  real ( kind = 8 ) sum1
  real ( kind = 8 ) tol
  real ( kind = 8 ), parameter :: tolrte = 0.000001D+00
  real ( kind = 8 ) xij

  nmk = n - k
  
  if ( n < k ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'SPLOPT - Fatal error!'
    write ( *, '(a)' ) '  N < K.'
    iflag = 2
    return
  end if
  
  if ( n == k ) then
    t(1:k) = tau(1)
    t(n+1:n+k) = tau(n)
    return
  end if
  
  if ( k <= 2 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'SPLOPT - Fatal error!'
    write ( *, '(a)' ) '  K < 2.'
    iflag = 2
    stop
  end if
 
  floatk = k
  kp1 = k + 1
  kpkm1 = k + k - 1
  kpn = k + n

  signst = -1.0D+00
  if ( ( nmk / 2 ) * 2 < nmk ) then
    signst = 1.0D+00
  end if
!
!  SCRTCH(I) = TAU-EXTENDED(I), I=1,...,N+K+K
!
  nx = n + k + k + 1
!
!  SCRTCH(I+NX) = XI(I), I=0,...,N-K+1
!
  na = nx + nmk + 1
!
!  SCRTCH(I+NA) = - A(I), I=1,...,N
!
  nd = na + n
!
!  SCRTCH(I+ND) = X(I) or D(I), I=1,...,N-K
!
  nb = nd + nmk
!
!  SCRTCH(I+NB) = BIATX(I), I=1,...,K+1
!
  nc = nb + kp1
!
!  SCRTCH(I+(J-1)*(2K-1)+NC) = W(I,J) = C(I-K+J,J), I=J-K,...,J+K,
!                                                     J=1,...,N-K.
!
  lenw = kpkm1 * nmk
!
!  Extend TAU to a knot sequence and store in SCRTCH.
!
  scrtch(1:k) = tau(1)
  scrtch(k+1:k+n) = tau(1:n)
  scrtch(kpn+1:kpn+k) = tau(n)
!
!  First guess for SCRTCH (.+NX) = XI.
!
  scrtch(nx) = tau(1)
  scrtch(nmk+1+nx) = tau(n)
 
  do j = 1, nmk 
    scrtch(j+nx) = sum ( tau(j+1:j+k-1) ) / real ( k - 1, kind = 8 )
  end do
!
!  Last entry of SCRTCH (.+NA) = -A  is always ...
!
  scrtch(n+na) = 0.5D+00
!
!  Start the Newton iteration.
!
  newton = 1
  tol = tolrte * ( tau(n) - tau(1) ) / real ( nmk, kind = 8 )
!
!  Start the Newton step.
!  Compute the 2*K-1 bands of the matrix C and store in SCRTCH(.+NC),
!  and compute the vector SCRTCH(.+NA) = -A.
!
  do newton = 1, newtmx

    scrtch(nc+1:nc+lenw) = 0.0D+00
    scrtch(na+1:na+n-1) = 0.0D+00
  
    sign = signst
    left = kp1
  
    do j = 1, nmk
  
      xij = scrtch(j+nx)

      do

        if ( xij < scrtch(left+1) ) then
          exit
        end if

        left = left + 1

        if ( kpn <= left ) then
          left = left - 1
          exit
        end if

      end do

      call bsplvb ( scrtch, k, 1, xij, left, scrtch(1+nb) )
!
!  The TAU sequence in SCRTCH is preceded by K additional knots.
!
!  Therefore, SCRTCH(LL+NB) now contains B(LEFT-2K+LL)(XIJ)
!  which is destined for C(LEFT-2K+LL,J), and therefore for
!
!    W(LEFT-K-J+LL,J)= SCRTCH(LEFT-K-J+LL+(J-1)*KPKM1 + NC)
!
!  since we store the 2*K-1 bands of C in the 2*K-1 rows of
!  the work array W, and W in turn is stored in SCRTCH,
!  with W(1,1) = SCRTCH(1+NC).
!
!  Also, C being of order N - K, we would want  
!    1 <= LEFT-2K+LL <= N - K or
!    LLMIN=2K-LEFT  <=  LL <= N-LEFT+K = LLMAX.
!
      leftmk = left - k
      index = leftmk - j + ( j - 1 ) * kpkm1 + nc
      llmin = max ( 1, k - leftmk )
      llmax = min ( k, n - leftmk )
      do ll = llmin, llmax
        scrtch(ll+index) = scrtch(ll+nb)
      end do
    
      call bsplvb ( scrtch, kp1, 2, xij, left, scrtch(1+nb) )

      id = max ( 0, leftmk - kp1 )
      llmin = 1 - min ( 0, leftmk - kp1 )
      do ll = llmin, kp1
        id = id + 1
        scrtch(id+na) = scrtch(id+na) - sign * scrtch(ll+nb)
      end do
    
      sign = - sign
    
    end do
  
    call banfac ( scrtch(1+nc), kpkm1, nmk, k-1, k-1, iflag )
  
    if ( iflag == 2 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'SPLOPT - Fatal error!'
      write ( *, '(a)' ) '  Matrix C is not invertible.'
      stop
    end if
!
!  Compute SCRTCH(.+ND) = D from SCRTCH(.+NA) = -A.
!
    do i = n, 2, -1
      scrtch(i-1+na) = scrtch(i-1+na) + scrtch(i+na)
    end do
  
    do i = 1, nmk
      scrtch(i+nd) = scrtch(i+na) * ( tau(i+k) - tau(i) ) / floatk
    end do
!
!  Compute SCRTCH(.+ND)= X.
!
    call banslv ( scrtch(1+nc), kpkm1, nmk, k-1, k-1, scrtch(1+nd) )
!
!  Compute SCRTCH(.+ND) = change in XI.  Modify, if necessary, to
!  prevent new XI from moving more than 1/3 of the way to its
!  neighbors.  Then add to XI to obtain new XI in SCRTCH(.+NX).
!
    delmax = 0.0D+00
    sign = signst

    do i = 1, nmk

      del = sign * scrtch(i+nd)
      delmax = max ( delmax, abs ( del ) )

      if ( 0.0D+00 < del ) then
        del = min ( del, ( scrtch(i+1+nx) - scrtch(i+nx) ) / 3.0D+00 )
      else
        del = max ( del, ( scrtch(i-1+nx) - scrtch(i+nx) ) / 3.0D+00 )
      end if

      sign = - sign
      scrtch(i+nx) = scrtch(i+nx) + del

    end do
!
!  Call it a day in case change in XI was small enough or too many
!  steps were taken.
!
    if ( delmax < tol ) then
      exit
    end if

  end do

  if ( tol <= delmax ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'SPLOPT - Warning!'
    write ( *, '(a)' ) '  The Newton iteration did not converge.'
  end if

  t(1:k) = tau(1)
  t(k+1:n) = scrtch(nx+1:nx+n-k)  
  t(n+1:n+k) = tau(n)
  
  return
end
subroutine subbak ( w, ipivot, nrow, ncol, last, x )

!*****************************************************************************80
!
!! SUBBAK carries out back substitution for the current block.
!
!  Modified:
!
!    14 February 2007
!
!  Author:
!
!    Carl DeBoor
!
!  Reference:
!
!    Carl DeBoor,
!    A Practical Guide to Splines,
!    Springer, 2001,
!    ISBN: 0387953663,
!    LC: QA1.A647.v27.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) W(NROW,NCOL), integer IPIVOT(NROW), integer
!    NROW, integer NCOL, integer LAST, are as on return from FACTRB.
!
!    Input/output, real ( kind = 8 ) X(NCOL).
!    On input, the right hand side for the equations in this block after 
!    back substitution has been carried out up to, but not including,
!    equation IPIVOT(LAST).  This means that X(1:LAST) contains the right hand
!    sides of equation IPIVOT(1:LAST) as modified during elimination,
!    while X(LAST+1:NCOL) is already a component of the solution vector.
!    On output, the components of the solution corresponding to the present
!    block.
!
  implicit none

  integer ( kind = 4 ) ncol
  integer ( kind = 4 ) nrow

  integer ( kind = 4 ) ip
  integer ( kind = 4 ) ipivot(nrow)
  integer ( kind = 4 ) k
  integer ( kind = 4 ) last
  real ( kind = 8 ) w(nrow,ncol)
  real ( kind = 8 ) x(ncol)

  do k = last, 1, -1

    ip = ipivot(k)
   
    x(k) = ( x(k) - dot_product ( w(ip,k+1:ncol), x(k+1:ncol) ) ) / w(ip,k)

  end do

  return
end
subroutine subfor ( w, ipivot, nrow, last, b, x )

!*****************************************************************************80
!
!! SUBFOR carries out the forward pass of substitution for the current block.
!
!  Discussion:
!
!    The forward pass is the action on the right hand side corresponding to the 
!    elimination carried out in FACTRB for this block.
!
!    At the end, X(1:NROW) contains the right hand side of the transformed
!    IPIVOT(1:NROW)-th equation in this block.  
!
!    Then, since for I=1,...,NROW-LAST, B(NROW+I) is going to be used as 
!    the right hand side of equation I in the next block (shifted over there 
!    from this block during factorization), it is set equal to X(LAST+I) here.
!
!  Modified:
!
!    14 February 2007
!
!  Author:
!
!    Carl DeBoor
!
!  Reference:
!
!    Carl DeBoor,
!    A Practical Guide to Splines,
!    Springer, 2001,
!    ISBN: 0387953663,
!    LC: QA1.A647.v27.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) W(NROW,LAST), integer IPIVOT(NROW), 
!    integer ( kind = 4 ) NROW, integer LAST, are as on return from FACTRB.
!
!    Output, real ( kind = 8 ) B(2*NROW-LAST).  On input, B(1:NROW)
!    contains the right hand sides for this block.  On output,
!    B(NROW+1:2*NROW-LAST) contains the appropriately modified right
!    hand sides for the next block.
!
!    Output, real X(NROW), contains, on output, the appropriately modified
!    right hand sides of equations IPIVOT(1:NROW).
!
  implicit none

  integer ( kind = 4 ) last
  integer ( kind = 4 ) nrow

  real ( kind = 8 ) b(nrow+nrow-last)
  integer ( kind = 4 ) ip
  integer ( kind = 4 ) ipivot(nrow)
  integer ( kind = 4 ) jhi
  integer ( kind = 4 ) k
  real ( kind = 8 ) w(nrow,last)
  real ( kind = 8 ) x(nrow)

  ip = ipivot(1)
  x(1) = b(ip)
  
  do k = 2, nrow
  
    ip = ipivot(k)

    jhi = min ( k - 1, last )
    
    x(k) = b(ip) - dot_product ( w(ip,1:jhi), x(1:jhi) )

  end do
!
!  Transfer modified right hand sides of equations IPIVOT(LAST+1:NROW) 
!  to next block.
!
  b(nrow+1:2*nrow-last) = x(last+1:nrow)
  
  return
end
subroutine tautsp ( tau, gtau, ntau, gamma, s, break, coef, l, k, iflag )

!*****************************************************************************80
!
!! TAUTSP constructs a cubic spline interpolant to given data.
!
!  Discussion:
!
!    If 0 < GAMMA, additional knots are introduced where needed to
!    make the interpolant more flexible locally.  This avoids extraneous
!    inflection points typical of cubic spline interpolation at knots to
!    rapidly changing data.
!
!  Method:  
!
!    On the I-th interval, (TAU(I), TAU(I+1)), the interpolant is of the
!    form:
!
!    (*)  F(U(X)) = A + B * U + C * H(U,Z) + D * H(1-U,1-Z),
!
!    with  
!
!      U = U(X) = ( X - TAU(I) ) / DTAU(I). 
!
!    Here,
!
!      Z(I) = ADDG(I+1) / ( ADDG(I) + ADDG(I+1) )
!
!    but if the denominator vanishes, we set Z(I) = 0.5
!
!    Also, we have
!
!      ADDG(J) = abs ( DDG(J) ), 
!      DDG(J) = DG(J+1) - DG(J),
!      DG(J) = DIVDIF(J) = ( GTAU(J+1) - GTAU(J) ) / DTAU(J)
!
!    and
!
!      H(U,Z) = ALPHA * U**3 
!             + ( 1 - ALPHA ) * ( max ( ( ( U - ZETA ) / ( 1 - ZETA ) ), 0 )**3
!
!    with
!
!      ALPHA(Z) = ( 1 - GAMMA / 3 ) / ZETA
!      ZETA(Z) = 1 - GAMMA * min ( ( 1 - Z ), 1/3 )
!
!    Thus, for 1/3 <= Z <= 2/3, F is just a cubic polynomial on
!    the interval I.  Otherwise, it has one additional knot, at
!
!      TAU(I) + ZETA * DTAU(I).
!
!    As Z approaches 1, H(.,Z) has an increasingly sharp bend near 1,
!    thus allowing F to turn rapidly near the additional knot.
!
!    In terms of F(J) = GTAU(J) and FSECND(J) = second derivative of F 
!    at TAU(J), the coefficients for (*) are given as:
!
!      A = F(I) - D
!      B = ( F(I+1) - F(I) ) - ( C - D )
!      C = FSECND(I+1) * DTAU(I)**2 / HSECND(1,Z)
!      D = FSECND(I) * DTAU(I)**2 / HSECND(1,1-Z)
!
!    Hence these can be computed once FSECND(1:NTAU) is fixed.
!
!    F is automatically continuous and has a continuous second derivative
!    except when Z=0 or 1 for some I.  We determine FSECND from
!    the requirement that the first derivative of F be continuous.
!
!    In addition, we require that the third derivative be continuous
!    across TAU(2) and across TAU(NTAU-1).  This leads to a strictly
!    diagonally dominant tridiagonal linear system for the FSECND(I)
!    which we solve by Gauss elimination without pivoting.
!
!  Modified:
!
!    14 February 2007
!
!  Author:
!
!    Carl DeBoor
!
!  Reference:
!
!    Carl DeBoor,
!    A Practical Guide to Splines,
!    Springer, 2001,
!    ISBN: 0387953663,
!    LC: QA1.A647.v27.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) TAU(NTAU), the sequence of data points.  
!    TAU must be strictly increasing.
!
!    Input, real ( kind = 8 ) GTAU(NTAU), the corresponding sequence of
!    function values.
!
!    Input, integer ( kind = 4 ) NTAU, the number of data points.  
!    NTAU must be at least 4.
!
!    Input, real ( kind = 8 ) GAMMA, indicates whether additional flexibility
!    is desired.
!    GAMMA = 0.0, no additional knots;
!    GAMMA in (0.0,3.0), under certain conditions on the given data at
!    points I-1, I, I+1, and I+2, a knot is added in the I-th interval, 
!    for I = 2,...,NTAU-2.  See description of method.  The interpolant 
!    gets rounded with increasing gamma.  A value of 2.5 for GAMMA is typical.
!    GAMMA in (3.0,6.0), same, except that knots might also be added in
!    intervals in which an inflection point would be permitted.  A value 
!    of 5.5 for GAMMA is typical.
!
!    Output, real ( kind = 8 ) BREAK(L), real ( kind = 8 ) COEF(K,L), 
!    integer ( kind = 4 ) L, integer K, give the piecewise polynomial 
!    representation of the interpolant.  Specifically, 
!    for BREAK(i) <= X <= BREAK(I+1), the interpolant has the form:
!      F(X) = COEF(1,I) +  DX    * ( 
!             COEF(2,I) + (DX/2) * (
!             COEF(3,I) + (DX/3) *
!             COEF(4,I) ) )
!    with  DX = X - BREAK(I) for I = 1,..., L.
!
!    Output, integer ( kind = 4 ) IFLAG, error flag.
!    1, no error.
!    2, input was incorrect.
!
!    Output, real ( kind = 8 ) S(NTAU,6).  The individual columns of this
!    array contain the following quantities mentioned in the write up
!    and below.
!    S(.,1) = DTAU = TAU(.+1)-TAU;
!    S(.,2) = DIAG = diagonal in linear system;
!    S(.,3) = U = upper diagonal in linear system;
!    S(.,4) = R = right hand side for linear system (initially)
!           = FSECND = solution of linear system, namely the second
!             derivatives of interpolant at TAU;
!    S(.,5) = Z = indicator of additional knots;
!    S(.,6) = 1/HSECND(1,X) with X = Z or 1-Z.
!
  implicit none

  integer ( kind = 4 ) ntau

  real ( kind = 8 ) alph
  real ( kind = 8 ) alpha
  real ( kind = 8 ) break(*)
  real ( kind = 8 ) c
  real ( kind = 8 ) coef(4,*)
  real ( kind = 8 ) d
  real ( kind = 8 ) del
  real ( kind = 8 ) denom
  real ( kind = 8 ) divdif
  real ( kind = 8 ) entry
  real ( kind = 8 ) entry3
  real ( kind = 8 ) factor
  real ( kind = 8 ) factr2
  real ( kind = 8 ) gam
  real ( kind = 8 ) gamma
  real ( kind = 8 ) gtau(ntau)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) iflag
  integer ( kind = 4 ) k
  integer ( kind = 4 ) l
  integer ( kind = 4 ) method
  real ( kind = 8 ) onemg3
  real ( kind = 8 ) onemzt
  real ( kind = 8 ) ratio
  real ( kind = 8 ) s(ntau,6)
  real ( kind = 8 ) sixth
  real ( kind = 8 ) tau(ntau)
  real ( kind = 8 ) temp
  real ( kind = 8 ) x
  real ( kind = 8 ) z
  real ( kind = 8 ) zeta
  real ( kind = 8 ) zt2

  alph(x) = min ( 1.0D+00, onemg3 / x )
!
!  There must be at least 4 interpolation points.
!
  if ( ntau < 4 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TAUTSP - Fatal error!'
    write ( *, '(a)' ) '  Input NTAU must be at least 4.'
    write ( *, '(a,i8)' ) '  NTAU = ', ntau
    iflag = 2
    stop
  end if
!
!  Construct delta TAU and first and second (divided) differences of data.
!
  do i = 1, ntau - 1

    s(i,1) = tau(i+1) - tau(i)
    
    if ( s(i,1) <= 0.0D+00 ) then
      write ( *, '(a,i3,a,2e15.6,a)' ) &
        '  Point ', i, ' and the next ', tau(i), tau(i+1), ' are disordered.'
      iflag = 2
      return
    end if
    
    s(i+1,4) = ( gtau(i+1) - gtau(i) ) / s(i,1)

  end do
   
  do i = 2, ntau - 1
    s(i,4) = s(i+1,4) - s(i,4)
  end do
!
!  Construct system of equations for second derivatives at TAU.
!
!  At each interior data point, there is one continuity equation.
!  At the first and last interior data point there is an additional 
!  equation for a total of NTAU equations in NTAU unknowns.
!
  i = 2
  s(2,2) = s(1,1) / 3.0D+00
  sixth = 1.0D+00 / 6.0D+00
  method = 2
  gam = gamma

  if ( gam <= 0.0D+00 ) then
    method = 1
  end if

  if ( 3.0D+00 < gam ) then
    method = 3
    gam = gam - 3.0D+00
  end if

  onemg3 = 1.0D+00 - gam / 3.0D+00
!
!  Loop over I.
!
   70 continue
!
!  Construct Z(I) and ZETA(I).
!
  z = 0.5D+00

  if ( method == 1 ) then
    go to 100
  end if

  if ( method == 3 ) then
    go to 90
  end if

  if ( s(i,4) * s(i+1,4) < 0.0D+00 ) then
    go to 100
  end if

   90 continue

  temp = abs ( s(i+1,4) )
  denom = abs ( s(i,4) ) + temp
  
  if ( denom /= 0.0D+00 ) then
    z = temp / denom
    if ( abs ( z - 0.5D+00 ) <= sixth ) then
      z = 0.5D+00
    end if
  end if
  
  100 continue

  s(i,5) = z
!
!  Set up part of the I-th equation which depends on the I-th interval.
!
  if ( z < 0.5D+00 ) then

    zeta = gam * z
    onemzt = 1.0D+00 - zeta
    zt2 = zeta**2
    alpha = alph(onemzt)
    factor = zeta / ( alpha * ( zt2 - 1.0D+00 ) + 1.0D+00 )
    s(i,6) = zeta * factor / 6.0D+00
    s(i,2) = s(i,2) + s(i,1) &
      * ( ( 1.0D+00 - alpha * onemzt ) * factor / 2.0D+00 - s(i,6) )
!
!  If Z = 0 and the previous Z = 1, then D(I) = 0.  
!  Since then also U(I-1) = L(I+1) = 0, its value does not matter.  
!  Reset D(I) = 1 to insure nonzero pivot in elimination.
!
    if ( s(i,2) <= 0.0D+00 ) then
      s(i,2) = 1.0D+00
    end if

    s(i,3) = s(i,1) / 6.0D+00

  else if ( z == 0.5D+00 ) then

    s(i,2) = s(i,2) + s(i,1) / 3.0D+00
    s(i,3) = s(i,1) / 6.0D+00

  else if ( 0.5D+00 < z ) then

    onemzt = gam * ( 1.0D+00 - z )
    zeta = 1.0D+00 - onemzt
    alpha = alph(zeta)
    factor = onemzt / ( 1.0D+00 - alpha * zeta * ( 1.0D+00 + onemzt ) )
    s(i,6) = onemzt * factor / 6.0D+00
    s(i,2) = s(i,2) + s(i,1) / 3.0D+00
    s(i,3) = s(i,6) * s(i,1)

  end if

  if ( 2 < i ) then
    go to 190
  end if

  s(1,5) = 0.5D+00
!
!  The first two equations enforce continuity of the first and of
!  the third derivative across TAU(2).
!
  s(1,2) = s(1,1) / 6.0D+00
  s(1,3) = s(2,2)
  entry3 = s(2,3)

  if ( z < 0.5D+00 ) then

    factr2 = zeta * ( alpha * ( zt2 - 1.0D+00 ) + 1.0D+00 ) &
      / ( alpha * ( zeta * zt2 - 1.0D+00 ) + 1.0D+00 )

    ratio = factr2 * s(2,1) / s(1,2)
    s(2,2) = factr2 * s(2,1) + s(1,1)
    s(2,3) = - factr2 * s(1,1)
  
  else if ( z == 0.5D+00 ) then

    ratio = s(2,1) / s(1,2)
    s(2,2) = s(2,1) + s(1,1)
    s(2,3) = - s(1,1)
  
  else if ( 0.5D+00 < z ) then

    ratio = s(2,1) / s(1,2)
    s(2,2) = s(2,1) + s(1,1)
    s(2,3) = - s(1,1) * 6.0D+00 * alpha * s(2,6)

  end if
!
!  At this point, the first two equations read:
!
!              DIAG(1)*X1+U(1)*X2 + ENTRY3*X3 = R(2)
!       -RATIO*DIAG(1)*X1+DIAG(2)*X2 + U(2)*X3 = 0.0
!
!  Eliminate first unknown from second equation.
!
  s(2,2) = ratio * s(1,3) + s(2,2)
  s(2,3) = ratio * entry3 + s(2,3)
  s(1,4) = s(2,4)
  s(2,4) = ratio * s(1,4)

  go to 200
  
  190 continue
!
!  The I-th equation enforces continuity of the first derivative
!  across TAU(I).  It now reads:
!
!    - RATIO * DIAG(I-1) * X(I-1) + DIAG(I) * X(I) + U(I) * X(I+1) = R(I).
!
!  Eliminate (I-1)st unknown from this equation
!
  s(i,2) = ratio * s(i-1,3) + s(i,2)
  s(i,4) = ratio * s(i-1,4) + s(i,4)
!
!  Set up the part of the next equation which depends on the I-th interval.
!
  200 continue

  if ( z < 0.5D+00 ) then

    ratio = - s(i,6) * s(i,1) / s(i,2)
    s(i+1,2) = s(i,1) / 3.0D+00

  else if ( z == 0.5D+00 ) then

    ratio = - ( s(i,1) / 6.0D+00 ) / s(i,2)
    s(i+1,2) = s(i,1) / 3.0D+00

  else if ( 0.5D+00 < z ) then

    ratio = - ( s(i,1) / 6.0D+00 ) / s(i,2)
    s(i+1,2) = s(i,1) &
      * ( ( 1.0D+00 - zeta * alpha ) * factor / 2.0D+00 - s(i,6) )

  end if
!
!  End of I loop.
!
  i = i + 1

  if ( i < ntau - 1 ) then
    go to 70
  end if

  s(i,5) = 0.5D+00
!
!  The last two equations enforce continuity of third derivative and
!  of first derivative across TAU(NTAU-1).
!
  entry = ratio * s(i-1,3) + s(i,2) + s(i,1) / 3.0D+00
  s(i+1,2) = s(i,1) / 6.0D+00
  s(i+1,4) = ratio * s(i-1,4) + s(i,4)

  if ( z < 0.5D+00 ) then

    ratio = s(i,1) * 6.0D+00 * s(i-1,6) * alpha / s(i-1,2)
    s(i,2) = ratio * s(i-1,3) + s(i,1) + s(i-1,1)
    s(i,3) = - s(i-1,1)
  
  else if ( z == 0.5D+00 ) then

    ratio = s(i,1) / s(i-1,2)
    s(i,2) = ratio * s(i-1,3) + s(i,1) + s(i-1,1)
    s(i,3) = - s(i-1,1)
  
  else if ( 0.5D+00 < z ) then

    factr2 = onemzt * ( alpha * ( onemzt**2 - 1.0D+00 ) + 1.0D+00 ) &
      / ( alpha * ( onemzt**3 - 1.0D+00 ) + 1.0D+00 )

    ratio = factr2 * s(i,1) / s(i-1,2)
    s(i,2) = ratio * s(i-1,3) + factr2 * s(i-1,1) + s(i,1)
    s(i,3) = - factr2 * s(i-1,1)

  end if
!
!  At this point, the last two equations read:
!
!           DIAG(I)*XI+     U(I)*XI+1 = R(I)
!    -RATIO*DIAG(I)*XI+DIAG(I+1)*XI+1 = R(I+1)
!
!  Eliminate XI from the last equation.
!
  s(i,4) = ratio * s(i-1,4)
  ratio = - entry / s(i,2)
  s(i+1,2) = ratio * s(i,3) + s(i+1,2)
  s(i+1,4) = ratio * s(i,4) + s(i+1,4)
!
!  Back substitution.
!
  s(ntau,4) = s(ntau,4) / s(ntau,2)

  do while ( 1 < i )

    s(i,4) = ( s(i,4) - s(i,3) * s(i+1,4) ) / s(i,2)
    i = i - 1

  end do

  s(1,4) = ( s(1,4) - s(1,3) * s(2,4) - entry3 * s(3,4) ) / s(1,2)
!
!  Construct polynomial pieces. 
!
  break(1) = tau(1)
  l = 1

  do i = 1, ntau - 1

    coef(1,l) = gtau(i)
    coef(3,l) = s(i,4)
    divdif = ( gtau(i+1) - gtau(i) ) / s(i,1)
    z = s(i,5)

    if ( z == 0.0D+00 ) then

      coef(2,l) = divdif
      coef(3,l) = 0D+00
      coef(4,l) = 0.0D+00

    else if ( z < 0.5D+00 ) then

      zeta = gam * z
      onemzt = 1.0D+00 - zeta
      c = s(i+1,4) / 6.0D+00
      d = s(i,4) * s(i,6)
      l = l + 1
      del = zeta * s(i,1)
      break(l) = tau(i) + del
      zt2 = zeta**2
      alpha = alph(onemzt)
      factor = onemzt**2 * alpha
      coef(1,l) = gtau(i) + divdif * del &
        + s(i,1)**2 * ( d * onemzt * ( factor - 1.0D+00 ) &
        + c * zeta * ( zt2 - 1.0D+00 ) )
      coef(2,l) = divdif + s(i,1) * ( d * ( 1.0D+00 - 3.0D+00 * factor ) &
        + c * ( 3.0D+00 * zt2 - 1.0D+00 ) )
      coef(3,l) = 6.0D+00 * ( d * alpha * onemzt + c * zeta )
      coef(4,l) = 6.0D+00 * ( c - d * alpha ) / s(i,1)
      coef(4,l-1) = coef(4,l) &
        - 6.0D+00 * d * ( 1.0D+00 - alpha ) / ( del * zt2 )
      coef(2,l-1) = coef(2,l) - del * ( coef(3,l) &
        - ( del / 2.0D+00 ) * coef(4,l-1))
 
    else if ( z == 0.5D+00 ) then

      coef(2,l) = divdif &
        - s(i,1) * ( 2.0D+00 * s(i,4) + s(i+1,4) ) / 6.0D+00
      coef(4,l) = ( s(i+1,4) - s(i,4) ) / s(i,1)

    else if ( 0.5D+00 <= z ) then

      onemzt = gam * ( 1.0D+00 - z )

      if ( onemzt == 0.0D+00 ) then

        coef(2,l) = divdif
        coef(3,l) = 0D+00
        coef(4,l) = 0.0D+00

      else

        zeta = 1.0D+00 - onemzt
        alpha = alph(zeta)
        c = s(i+1,4) * s(i,6)
        d = s(i,4) / 6.0D+00
        del = zeta * s(i,1)
        break(l+1) = tau(i) + del
        coef(2,l) = divdif - s(i,1) * ( 2.0D+00 * d + c )
        coef(4,l) = 6.0D+00 * ( c * alpha - d ) / s(i,1)
        l = l + 1
        coef(4,l) = coef(4,l-1) + 6.0D+00 * ( 1.0D+00 - alpha ) * c &
          / ( s(i,1) * onemzt**3 )
        coef(3,l) = coef(3,l-1) + del * coef(4,l-1)
        coef(2,l) = coef(2,l-1) + del * ( coef(3,l-1) &
          + ( del / 2.0D+00 ) * coef(4,l-1) )
        coef(1,l) = coef(1,l-1) + del * ( coef(2,l-1) &
          + ( del / 2.0D+00 ) * ( coef(3,l-1) &
          + ( del / 3.0D+00 ) * coef(4,l-1) ) )

      end if

    end if

    l = l + 1
    break(l) = tau(i+1)
    
  end do
  
  l = l - 1
  k = 4
  iflag = 1

  return
end
subroutine titand ( t, g, n )

!*****************************************************************************80
!
!! TITAND represents a temperature dependent property of titanium.
!
!  Discussion:
!
!    The data has been used extensively as an example in spline
!    approximation with variable knots.
!
!  Modified:
!
!    14 February 2007
!
!  Author:
!
!    Carl DeBoor
!
!  Reference:
!
!    Carl DeBoor,
!    A Practical Guide to Splines,
!    Springer, 2001,
!    ISBN: 0387953663,
!    LC: QA1.A647.v27.
!
!  Parameters:
!
!    Output, real ( kind = 8 ) T(N), the location of the data points.
!
!    Output, real ( kind = 8 ) G(N), the value associated with the data points.
!
!    Output, integer ( kind = 4 ) N, the number of data points, which is 49.
!
  implicit none

  real ( kind = 8 ) g(*)
  integer ( kind = 4 ) n
  real ( kind = 8 ) t(*)

  n = 49

  t(1:49) = (/ &
     595.0D+00,  605.0D+00,  615.0D+00,  625.0D+00,  635.0D+00, &
     645.0D+00,  655.0D+00,  665.0D+00,  675.0D+00,  685.0D+00, &
     695.0D+00,  705.0D+00,  715.0D+00,  725.0D+00,  735.0D+00, &
     745.0D+00,  755.0D+00,  765.0D+00,  775.0D+00,  785.0D+00, &
     795.0D+00,  805.0D+00,  815.0D+00,  825.0D+00,  835.0D+00, &
     845.0D+00,  855.0D+00,  865.0D+00,  875.0D+00,  885.0D+00, &
     895.0D+00,  905.0D+00,  915.0D+00,  925.0D+00,  935.0D+00, &
     945.0D+00,  955.0D+00,  965.0D+00,  975.0D+00,  985.0D+00, &
     995.0D+00, 1005.0D+00, 1015.0D+00, 1025.0D+00, 1035.0D+00, &
    1045.0D+00, 1055.0D+00, 1065.0D+00, 1075.0D+00 /)

  g(1:49) = (/ &
    0.644D+00, 0.622D+00, 0.638D+00, 0.649D+00, 0.652D+00, &
    0.639D+00, 0.646D+00, 0.657D+00, 0.652D+00, 0.655D+00, &
    0.644D+00, 0.663D+00, 0.663D+00, 0.668D+00, 0.676D+00, &
    0.676D+00, 0.686D+00, 0.679D+00, 0.678D+00, 0.683D+00, &
    0.694D+00, 0.699D+00, 0.710D+00, 0.730D+00, 0.763D+00, &
    0.812D+00, 0.907D+00, 1.044D+00, 1.336D+00, 1.881D+00, &
    2.169D+00, 2.075D+00, 1.598D+00, 1.211D+00, 0.916D+00, &
    0.746D+00, 0.672D+00, 0.627D+00, 0.615D+00, 0.607D+00, &
    0.606D+00, 0.609D+00, 0.603D+00, 0.601D+00, 0.603D+00, &
    0.601D+00, 0.611D+00, 0.601D+00, 0.608D+00 /)

  return
end
