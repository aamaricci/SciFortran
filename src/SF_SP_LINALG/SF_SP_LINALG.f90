module SF_SP_LINALG
  USE SF_LINALG, only: eye
#ifdef _MPI
  USE SF_MPI
  USE MPI
#endif
  implicit none
  private 

  interface sp_eigh
     module procedure :: lanczos_arpack_d
     module procedure :: lanczos_arpack_c
#ifdef _MPI
     module procedure :: lanczos_parpack_d
     module procedure :: lanczos_parpack_c
#endif
  end interface sp_eigh


  interface sp_lanc_eigh
     module procedure :: lanczos_eigh_d
     module procedure :: lanczos_eigh_c
#ifdef _MPI
     module procedure :: mpi_lanczos_eigh_d
     module procedure :: mpi_lanczos_eigh_c
#endif
  end interface sp_lanc_eigh


  interface sp_lanc_tridiag
     module procedure :: lanczos_tridiag_d
     module procedure :: lanczos_tridiag_c
#ifdef _MPI
     module procedure :: mpi_lanczos_tridiag_d
     module procedure :: mpi_lanczos_tridiag_c
#endif
  end interface sp_lanc_tridiag



  complex(8),parameter              :: zero=(0d0,0d0)
  complex(8),parameter              :: one=(1d0,0d0)
  complex(8),parameter              :: xi=(0d0,1d0)
  integer,allocatable               :: seed_random(:)
  integer                           :: nrandom
  logical                           :: verb=.false.
  real(8)                           :: threshold_=1.d-13
  integer                           :: ncheck_=10



  !****************************************************************************************
  !                                      PUBLIC 
  !****************************************************************************************
  public :: sp_eigh
  public :: sp_lanc_eigh
  public :: sp_lanc_tridiag
  !****************************************************************************************




contains



  !##################################################################
  ! ARPACK METHOD for LOWEST part of the spectrum of a of a SPARSE
  !   MATRIX (defined via H*v)
  ! - DBLE and CMPLX versions included
  ! - SERIAL and PARALLEL-MPI versions included
  !        [COMM    MPI  Communicator for the processor grid.  (INPUT)]
  !        IDO     Integer.  (INPUT/OUTPUT)
  !                Reverse communication flag.  IDO must be zero on the first
  !                call to pdsaupd .  IDO will be set internally to
  !                indicate the type of operation to be performed.  Control is
  !                then given back to the calling routine which has the
  !                responsibility to carry out the requested operation and call
  !                pdsaupd  with the result.  The operand is given in
  !                WORKD(IPNTR(1)), the result must be put in WORKD(IPNTR(2)).
  !                (If Mode = 2 see remark 5 below)
  !                -------------------------------------------------------------
  !                IDO =  0: first call to the reverse communication interface
  !                IDO = -1: compute  Y = OP * X  where
  !                          IPNTR(1) is the pointer into WORKD for X,
  !                          IPNTR(2) is the pointer into WORKD for Y.
  !                          This is for the initialization phase to force the
  !                          starting vector into the range of OP.
  !                IDO =  1: compute  Y = OP * X where
  !                          IPNTR(1) is the pointer into WORKD for X,
  !                          IPNTR(2) is the pointer into WORKD for Y.
  !                          In mode 3,4 and 5, the vector B * X is already
  !                          available in WORKD(ipntr(3)).  It does not
  !                          need to be recomputed in forming OP * X.
  !                IDO =  2: compute  Y = B * X  where
  !                          IPNTR(1) is the pointer into WORKD for X,
  !                          IPNTR(2) is the pointer into WORKD for Y.
  !                IDO =  3: compute the IPARAM(8) shifts where
  !                          IPNTR(11) is the pointer into WORKL for
  !                          placing the shifts. See remark 6 below.
  !                IDO = 99: done
  !                -------------------------------------------------------------
  !        BMAT    Character*1.  (INPUT)
  !                BMAT specifies the type of the matrix B that defines the
  !                semi-inner product for the operator OP.
  !                B = 'I' -> standard eigenvalue problem A*x = lambda*x
  !                B = 'G' -> generalized eigenvalue problem A*x = lambda*B*x
  !        N       Integer.  (INPUT)
  !                Dimension of the eigenproblem.
  !        WHICH   Character*2.  (INPUT)
  !                Specify which of the Ritz values of OP to compute.
  !                'LA' - compute the NEV largest (algebraic) eigenvalues.
  !                'SA' - compute the NEV smallest (algebraic) eigenvalues.
  !                'LM' - compute the NEV largest (in magnitude) eigenvalues.
  !                'SM' - compute the NEV smallest (in magnitude) eigenvalues.
  !                'BE' - compute NEV eigenvalues, half from each end of the
  !                       spectrum.  When NEV is odd, compute one more from the
  !                       high end than from the low end.
  !                 (see remark 1 below)
  !        NEV     Integer.  (INPUT)
  !                Number of eigenvalues of OP to be computed. 0 < NEV < N.
  !        TOL     Double precision  scalar.  (INPUT)
  !                Stopping criterion: the relative accuracy of the Ritz value
  !                is considered acceptable if BOUNDS(I) .LE. TOL*ABS(RITZ(I)).
  !                If TOL .LE. 0. is passed a default is set:
  !                DEFAULT = DLAMCH ('EPS')  (machine precision as computed
  !                          by the LAPACK auxiliary subroutine DLAMCH ).
  !        RESID   Double precision  array of length N.  (INPUT/OUTPUT)
  !                On INPUT:
  !                If INFO .EQ. 0, a random initial residual vector is used.
  !                If INFO .NE. 0, RESID contains the initial residual vector,
  !                                possibly from a previous run.
  !                On OUTPUT:
  !                RESID contains the final residual vector.
  !        NCV     Integer.  (INPUT)
  !                Number of columns of the matrix V (less than or equal to N).
  !                This will indicate how many Lanczos vectors are generated
  !                at each iteration.  After the startup phase in which NEV
  !                Lanczos vectors are generated, the algorithm generates
  !                NCV-NEV Lanczos vectors at each subsequent update iteration.
  !                Most of the cost in generating each Lanczos vector is in the
  !                matrix-vector product OP*x. (See remark 4 below).
  !        V       Double precision  N by NCV array.  (OUTPUT)
  !                The NCV columns of V contain the Lanczos basis vectors.
  !        LDV     Integer.  (INPUT)
  !                Leading dimension of V exactly as declared in the calling
  !                program.
  !        IPARAM  Integer array of length 11.  (INPUT/OUTPUT)
  !                IPARAM(1) = ISHIFT: method for selecting the implicit shifts.
  !                The shifts selected at each iteration are used to restart
  !                the Arnoldi iteration in an implicit fashion.
  !                -------------------------------------------------------------
  !                ISHIFT = 0: the shifts are provided by the user via
  !                            reverse communication.  The NCV eigenvalues of
  !                            the current tridiagonal matrix T are returned in
  !                            the part of WORKL array corresponding to RITZ.
  !                            See remark 6 below.
  !                ISHIFT = 1: exact shifts with respect to the reduced
  !                            tridiagonal matrix T.  This is equivalent to
  !                            restarting the iteration with a starting vector
  !                            that is a linear combination of Ritz vectors
  !                            associated with the "wanted" Ritz values.
  !                -------------------------------------------------------------
  !                IPARAM(2) = LEVEC
  !                No longer referenced. See remark 2 below.
  !                IPARAM(3) = MXITER
  !                On INPUT:  maximum number of Arnoldi update iterations allowed.
  !                On OUTPUT: actual number of Arnoldi update iterations taken.
  !                IPARAM(4) = NB: blocksize to be used in the recurrence.
  !                The code currently works only for NB = 1.
  !                IPARAM(5) = NCONV: number of "converged" Ritz values.
  !                This represents the number of Ritz values that satisfy
  !                the convergence criterion.
  !                IPARAM(6) = IUPD
  !                No longer referenced. Implicit restarting is ALWAYS used.
  !                IPARAM(7) = MODE
  !                On INPUT determines what type of eigenproblem is being solved.
  !                Must be 1,2,3,4,5; See under PURPOSE of pdsaupd  for the
  !                five modes available.
  !                IPARAM(8) = NP
  !                When ido = 3 and the user provides shifts through reverse
  !                communication (IPARAM(1)=0), pdsaupd  returns NP, the number
  !                of shifts the user is to provide. 0 < NP <=NCV-NEV. See Remark
  !                6 below.
  !                IPARAM(9) = NUMOP, IPARAM(10) = NUMOPB, IPARAM(11) = NUMREO,
  !                OUTPUT: NUMOP  = total number of OP*x operations,
  !                        NUMOPB = total number of B*x operations if BMAT='G',
  !                        NUMREO = total number of steps of re-orthogonalization.
  !        IPNTR   Integer array of length 11.  (OUTPUT)
  !                Pointer to mark the starting locations in the WORKD and WORKL
  !                arrays for matrices/vectors used by the Lanczos iteration.
  !                -------------------------------------------------------------
  !                IPNTR(1): pointer to the current operand vector X in WORKD.
  !                IPNTR(2): pointer to the current result vector Y in WORKD.
  !                IPNTR(3): pointer to the vector B * X in WORKD when used in
  !                          the shift-and-invert mode.
  !                IPNTR(4): pointer to the next available location in WORKL
  !                          that is untouched by the program.
  !                IPNTR(5): pointer to the NCV by 2 tridiagonal matrix T in WORKL.
  !                IPNTR(6): pointer to the NCV RITZ values array in WORKL.
  !                IPNTR(7): pointer to the Ritz estimates in array WORKL associated
  !                          with the Ritz values located in RITZ in WORKL.
  !                IPNTR(11): pointer to the NP shifts in WORKL. See Remark 6 below.
  !                Note: IPNTR(8:10) is only referenced by pdseupd . See Remark 2.
  !                IPNTR(8): pointer to the NCV RITZ values of the original system.
  !                IPNTR(9): pointer to the NCV corresponding error bounds.
  !                IPNTR(10): pointer to the NCV by NCV matrix of eigenvectors
  !                           of the tridiagonal matrix T. Only referenced by
  !                           pdseupd  if RVEC = .TRUE. See Remarks.
  !                -------------------------------------------------------------
  !        WORKD   Double precision  work array of length 3*N.  (REVERSE COMMUNICATION)
  !                Distributed array to be used in the basic Arnoldi iteration
  !                for reverse communication.  The user should not use WORKD
  !                as temporary workspace during the iteration. Upon termination
  !                WORKD(1:N) contains B*RESID(1:N). If the Ritz vectors are desired
  !                subroutine pdseupd  uses this output.
  !                See Data Distribution Note below.
  !        WORKL   Double precision  work array of length LWORKL.  (OUTPUT/WORKSPACE)
  !                Private (replicated) array on each PE or array allocated on
  !                the front end.  See Data Distribution Note below.
  !        LWORKL  Integer.  (INPUT)
  !                LWORKL must be at least NCV**2 + 8*NCV .
  !        INFO    Integer.  (INPUT/OUTPUT)
  !                If INFO .EQ. 0, a randomly initial residual vector is used.
  !                If INFO .NE. 0, RESID contains the initial residual vector,
  !                                possibly from a previous run.
  !                Error flag on output.
  !                =  0: Normal exit.
  !                =  1: Maximum number of iterations taken.
  !                      All possible eigenvalues of OP has been found. IPARAM(5)
  !                      returns the number of wanted converged Ritz values.
  !                =  2: No longer an informational error. Deprecated starting
  !                      with release 2 of ARPACK.
  !                =  3: No shifts could be applied during a cycle of the
  !                      Implicitly restarted Arnoldi iteration. One possibility
  !                      is to increase the size of NCV relative to NEV.
  !                      See remark 4 below.
  !                = -1: N must be positive.
  !                = -2: NEV must be positive.
  !                = -3: NCV must be greater than NEV and less than or equal to N.
  !                = -4: The maximum number of Arnoldi update iterations allowed
  !                      must be greater than zero.
  !                = -5: WHICH must be one of 'LM', 'SM', 'LA', 'SA' or 'BE'.
  !                = -6: BMAT must be one of 'I' or 'G'.
  !                = -7: Length of private work array WORKL is not sufficient.
  !                = -8: Error return from trid. eigenvalue calculation;
  !                      Informatinal error from LAPACK routine dsteqr .
  !                = -9: Starting vector is zero.
  !                = -10: IPARAM(7) must be 1,2,3,4,5.
  !                = -11: IPARAM(7) = 1 and BMAT = 'G' are incompatable.
  !                = -12: IPARAM(1) must be equal to 0 or 1.
  !                = -13: NEV and WHICH = 'BE' are incompatable.
  !                = -9999: Could not build an Arnoldi factorization.
  !                         IPARAM(5) returns the size of the current Arnoldi
  !                         factorization. The user is advised to check that
  !                         enough workspace and array storage has been allocated.
  !##################################################################
  include "arpack_serial.f90"
#ifdef _MPI
  include "arpack_mpi.f90"
#endif


  !##################################################################
  ! LANCZOS METHOD for LOWEST EigenSolution OR tri-diagonalization
  !    of a SPARSE MATRIX (defined via H*v)
  ! - DBLE and CMPLX versions included
  ! - SERIAL and PARALLEL-MPI versions included
  !##################################################################
  include "lanczos_serial.f90"
#ifdef _MPI
  include "lanczos_mpi.f90"
#endif











  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !++++++++++++++++++COMPUTATIONAL ROUTINE: TQL2++++++++++++++++++++++++ 
  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !---------------------------------------------------------------------
  ! PURPOSE computes all eigenvalues/vectors, real symmetric tridiagonal matrix.
  !    This subroutine finds the eigenvalues and eigenvectors of a symmetric
  !    tridiagonal matrix by the QL method.  The eigenvectors of a full
  !    symmetric matrix can also be found if TRED2 has been used to reduce this
  !    full matrix to tridiagonal form.
  !  Parameters:
  !    Input, integer ( kind = 4 ) N, the order of the matrix.
  !
  !    Input/output, real ( kind = 8 ) D(N).  On input, the diagonal elements of
  !    the matrix.  On output, the eigenvalues in ascending order.  If an error
  !    exit is made, the eigenvalues are correct but unordered for indices
  !    1,2,...,IERR-1.
  !
  !    Input/output, real ( kind = 8 ) E(N).  On input, E(2:N) contains the
  !    subdiagonal elements of the input matrix, and E(1) is arbitrary.
  !    On output, E has been destroyed.
  !
  !    Input, real ( kind = 8 ) Z(N,N).  On input, the transformation matrix
  !    produced in the reduction by TRED2, if performed.  If the eigenvectors of
  !    the tridiagonal matrix are desired, Z must contain the identity matrix.
  !    On output, Z contains the orthonormal eigenvectors of the symmetric
  !    tridiagonal (or full) matrix.  If an error exit is made, Z contains
  !    the eigenvectors associated with the stored eigenvalues.
  !
  !    Output, integer ( kind = 4 ) IERR, error flag.
  !    0, normal return,
  !    J, if the J-th eigenvalue has not been determined after
  !    30 iterations.
  !
  !---------------------------------------------------------------------
  subroutine tql2 ( n, d, e, z, ierr )
    integer :: n
    real(8) :: c
    real(8) :: c2
    real(8) :: c3
    real(8) :: d(n)
    real(8) :: dl1
    real(8) :: e(n)
    real(8) :: el1
    real(8) :: f
    real(8) :: g
    real(8) :: h
    integer ( kind = 4 ) i
    integer ( kind = 4 ) ierr
    integer ( kind = 4 ) ii
    integer ( kind = 4 ) j
    integer ( kind = 4 ) k
    integer ( kind = 4 ) l
    integer ( kind = 4 ) l1
    integer ( kind = 4 ) l2
    integer ( kind = 4 ) m
    integer ( kind = 4 ) mml
    real(8) :: p
    real(8) :: r
    real(8) :: s
    real(8) :: s2
    real(8) :: tst1
    real(8) :: tst2
    real(8) :: z(n,n)
    ierr = 0
    if ( n == 1 ) then
       return
    end if
    do i = 2, n
       e(i-1) = e(i)
    end do
    f = 0.0D+00
    tst1 = 0.0D+00
    e(n) = 0.0D+00
    do l = 1, n
       j = 0
       h = abs ( d(l) ) + abs ( e(l) )
       tst1 = max ( tst1, h )
       !
       !  Look for a small sub-diagonal element.
       !
       do m = l, n
          tst2 = tst1 + abs ( e(m) )
          if ( tst2 == tst1 ) then
             exit
          end if
       end do
       if ( m == l ) then
          go to 220
       end if
130    continue
       if ( 30 <= j ) then
          ierr = l
          return
       end if
       j = j + 1
       !
       !  Form shift.
       !
       l1 = l + 1
       l2 = l1 + 1
       g = d(l)
       p = ( d(l1) - g ) / ( 2.0D+00 * e(l) )
       r = pythag ( p, 1.0D+00 )
       d(l) = e(l) / ( p + sign ( r, p ) )
       d(l1) = e(l) * ( p + sign ( r, p ) )
       dl1 = d(l1)
       h = g - d(l)
       d(l2:n) = d(l2:n) - h
       f = f + h
       !
       !  QL transformation.
       !
       p = d(m)
       c = 1.0D+00
       c2 = c
       el1 = e(l1)
       s = 0.0D+00
       mml = m - l
       do ii = 1, mml
          c3 = c2
          c2 = c
          s2 = s
          i = m - ii
          g = c * e(i)
          h = c * p
          r = pythag ( p, e(i) )
          e(i+1) = s * r
          s = e(i) / r
          c = p / r
          p = c * d(i) - s * g
          d(i+1) = h + s * ( c * g + s * d(i) )
          !
          !  Form vector.
          !
          do k = 1, n
             h = z(k,i+1)
             z(k,i+1) = s * z(k,i) + c * h
             z(k,i) = c * z(k,i) - s * h
          end do
       end do
       p = - s * s2 * c3 * el1 * e(l) / dl1
       e(l) = s * p
       d(l) = c * p
       tst2 = tst1 + abs ( e(l) )
       if ( tst2 > tst1 ) then
          go to 130
       end if
220    continue
       d(l) = d(l) + f
    end do
    !
    !  Order eigenvalues and eigenvectors.
    !
    do ii = 2, n
       i = ii - 1
       k = i
       p = d(i)
       do j = ii, n
          if ( d(j) < p ) then
             k = j
             p = d(j)
          end if
       end do
       if ( k /= i ) then
          d(k) = d(i)
          d(i) = p
          do j = 1, n
             call r8_swap ( z(j,i), z(j,k) )
          end do
       end if
    end do
    return
  end subroutine tql2


  !---------------------------------------------------------------------
  ! PURPOSE: computes SQRT ( A * A + B * B ) carefully.
  !    The formula
  !    PYTHAG = sqrt ( A * A + B * B )
  !    is reasonably accurate, but can fail if, for example, A**2 is larger
  !    than the machine overflow.  The formula can lose most of its accuracy
  !    if the sum of the squares is very large or very small.
  !  Parameters:
  !    Input, real(8) :: A, B, the two legs of a right triangle.
  !    Output, real(8) :: PYTHAG, the length of the hypotenuse.
  !---------------------------------------------------------------------
  function pythag ( a, b )
    implicit none
    real(8) :: a
    real(8) :: b
    real(8) :: p
    real(8) :: pythag
    real(8) :: r
    real(8) :: s
    real(8) :: t
    real(8) :: u
    p = max ( abs ( a ), abs ( b ) )
    if ( p /= 0.0D+00 ) then
       r = ( min ( abs ( a ), abs ( b ) ) / p )**2
       do
          t = 4.0D+00 + r
          if ( t == 4.0D+00 ) then
             exit
          end if
          s = r / t
          u = 1.0D+00 + 2.0D+00 * s
          p = u * p
          r = ( s / u )**2 * r
       end do
    end if
    pythag = p
    return
  end function pythag

  !---------------------------------------------------------------------
  ! PURPOSE: swaps two R8's.
  !  Parameters:
  !    Input/output, real(8) :: X, Y.  On output, the values of X and
  !    Y have been interchanged.
  !---------------------------------------------------------------------
  subroutine r8_swap ( x, y )
    real(8) :: x
    real(8) :: y
    real(8) :: z
    z = x
    x = y
    y = z
    return
  end subroutine r8_swap




end module SF_SP_LINALG
