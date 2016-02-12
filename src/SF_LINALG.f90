module SF_LINALG
  USE SF_CONSTANTS, only: one,xi,zero
  implicit none
  private


  !>EIGENVALUE PROBLEM:
  !Eigenvalue/-vector problem for real symmetric/complex hermitian matrices using Jacobi method (unsorted out):
  public :: eigh_jacobi
  !Eigenvalue/-vector problem for general matrices:
  public :: eig
  !Eigenvalue/-vector problem for real symmetric/complex hermitian matrices:
  public :: eigh, matrix_diagonalize
  !Eigenvalues for general matrices:
  public :: eigvals
  !Eigenvalues for symmetric/hermitian matrices:
  public :: eigvalsh

  !>MATRIX INVERSION:
  !Matrix inversion for real/complex matrices:
  public :: inv, matrix_inverse
  !Matrix inversion for real/complex symmetric matrices:
  public :: inv_sym, matrix_inverse_sym
  ! matrix inversion for complex hermitian matrices:
  public :: inv_her, matrix_inverse_her
  ! matrix inversion for real/complex triangular matrices:
  public :: inv_triang, matrix_inverse_triang
  ! matrix inversion for real/complex  matrices using Gauss-Jordan elimination:
  public :: inv_gj, matrix_inverse_gj


  !>TRIDIAGONAL MATRICES INVERSION (diagonal elements only):
  !invert a (block) tridiagonal matrix using the iterative algorithm
  public :: inv_tridiag


  !>LINEAR SYSTEM SOLUTION:
  ! solution to linear systems of equation with real/complex coefficients:
  public :: solve, solve_linear_system
  !least square solutions the real/complex systems of equations of possibly non-square shape:
  public :: lstsq


  !>AUXILIARY:
  ! determinants of real/complex square matrices:
  public :: det
  !Returns the real/complex identity matrix of size n x n .
  public :: eye, deye, zeye, deye_tridiag, zeye_tridiag
  !check the matrix is actually (block) tridiagonal
  public :: check_tridiag
  !get the main (block) diagonals from a (block) tridiagonal matrix
  public :: get_tridiag
  !build a (block) tridigonal matrix from the main (block) diagonals
  public :: build_tridiag
  !construction of square matrices from the diagonal elements:
  public :: diag
  !trace of real/complex matrices:
  public :: trace


  !>SVD DECOMPOSITION:
  !singular values of real/complex matrices:
  public :: svdvals
  !singular value decomposition of real/complex matrices:
  public :: svd


  !>EXTERNAL PRODUCTS
  !Kroenecker product of matrices
  public :: kron
  public :: kronecker_product
  public :: kroenecker_product
  !outer product of two 1d arrays to form a matrix
  public :: outerprod
  public :: cross_product
  public :: s3_product


  !>EIGENVALUE PROBLEM:
  interface eig
     module procedure deig
     module procedure zeig
  end interface eig

  interface eigh
     module procedure deigh_generalized
     module procedure zeigh_generalized
     module procedure deigh_simple
     module procedure zeigh_simple
  end interface eigh
  interface matrix_diagonalize
     module procedure deigh_simple
     module procedure zeigh_simple
  end interface matrix_diagonalize

  interface eigh_jacobi
     module procedure d_jacobi
     module procedure c_jacobi
  end interface eigh_jacobi

  interface eigvals
     module procedure deigvals
     module procedure zeigvals
  end interface eigvals

  interface eigvalsh
     module procedure deigvalsh
     module procedure zeigvalsh
  end interface eigvalsh



  !>MATRIX INVERSION:
  interface inv
     module procedure dinv
     module procedure zinv
  end interface inv
  interface matrix_inverse
     module procedure dinv
     module procedure zinv
  end interface matrix_inverse

  interface inv_sym
     module procedure dinv_sym
     module procedure zinv_sym
  end interface inv_sym
  interface matrix_inverse_sym
     module procedure dinv_sym
     module procedure zinv_sym
  end interface matrix_inverse_sym

  interface inv_her
     module procedure zinv_her
  end interface inv_her
  interface matrix_inverse_her
     module procedure zinv_her
  end interface matrix_inverse_her

  interface inv_triang
     module procedure dinv_triang
     module procedure zinv_triang
  end interface inv_triang
  interface matrix_inverse_triang
     module procedure dinv_triang
     module procedure zinv_triang
  end interface matrix_inverse_triang

  interface inv_gj
     module procedure dinv_gj
     module procedure zinv_gj
  end interface inv_gj
  interface matrix_inverse_gj
     module procedure dinv_gj
     module procedure zinv_gj
  end interface matrix_inverse_gj

  !>TRIDIAGONAL MATRICES INVERSION (diagonal elements only):
  interface inv_tridiag
     module procedure d_invert_tridiag_matrix
     module procedure c_invert_tridiag_matrix
     module procedure d_invert_tridiag_block_matrix
     module procedure c_invert_tridiag_block_matrix
     module procedure d_invert_tridiag_matrix_mat
     module procedure c_invert_tridiag_matrix_mat
     module procedure d_invert_tridiag_block_matrix_mat
     module procedure c_invert_tridiag_block_matrix_mat
  end interface inv_tridiag




  !>LINEAR SYSTEM SOLUTION:
  interface solve
     module procedure dsolve_1rhs
     module procedure zsolve_1rhs
     module procedure dsolve_Mrhs
     module procedure zsolve_Mrhs
  end interface solve
  interface solve_linear_system
     module procedure dsolve_1rhs
     module procedure zsolve_1rhs
     module procedure dsolve_Mrhs
     module procedure zsolve_Mrhs
  end interface solve_linear_system

  interface lstsq
     module procedure dlstsq
     module procedure zlstsq
  end interface lstsq




  !>AUXILIARY:
  interface det
     module procedure ddet
     module procedure zdet
  end interface det

  interface eye
     module procedure deye
     module procedure deye_tridiag
  end interface eye

  interface eye_tridiag
     module procedure deye_tridiag
  end interface eye_tridiag

  interface diag
     module procedure ddiag
     module procedure zdiag
  end interface diag

  interface trace
     module procedure dtrace
     module procedure ztrace
  end interface trace

  interface check_tridiag
     module procedure d_check_tridiag
     module procedure c_check_tridiag
     module procedure d_check_tridiag_block
     module procedure c_check_tridiag_block
  end interface check_tridiag

  interface get_tridiag
     module procedure d_get_tridiag
     module procedure c_get_tridiag
     module procedure d_get_tridiag_block
     module procedure c_get_tridiag_block
  end interface get_tridiag

  interface build_tridiag
     module procedure d_build_tridiag
     module procedure c_build_tridiag
     module procedure d_build_tridiag_block
     module procedure c_build_tridiag_block
  end interface build_tridiag


  !>SVD DECOMPOSITION:
  interface svdvals
     module procedure dsvdvals
     module procedure zsvdvals
  end interface svdvals

  interface svd
     module procedure dsvd
     module procedure zsvd
  end interface svd




  !>EXTERNAL PRODUCTS
  interface kron
     module procedure i_kronecker_product,d_kronecker_product,c_kronecker_product
  end interface kron
  interface kronecker_product
     module procedure i_kronecker_product,d_kronecker_product,c_kronecker_product
  end interface kronecker_product
  interface kroenecker_product
     module procedure i_kronecker_product,d_kronecker_product,c_kronecker_product
  end interface kroenecker_product

  interface outerprod
     module procedure outerprod_d,outerprod_c
  end interface outerprod

  interface cross_product
     module procedure cross_3d_d
     module procedure cross_3d_c
  end interface cross_product

  interface s3_product
     module procedure s3_product_d
     module procedure s3_product_c
  end interface s3_product



  !NOT PUBLIC:
  !Assert shape of matrices:
  interface assert_shape
     module procedure dassert_shape
     module procedure zassert_shape
  end interface assert_shape

  !Swap two elements A and B (from nr90)
  interface swap
     module procedure swap_i,swap_r,swap_rv,swap_z,swap_zv,swap_zm
  end interface swap

  !Interface to Lapack function ilaenv
  interface
     integer function ilaenv( ispec, name, opts, n1, n2, n3, n4 )
       character*(*) name, opts
       integer       ispec, n1, n2, n3, n4
     end function ilaenv
  end interface






contains


  !+-----------------------------------------------------------------------------+!
  !PURPOSE:  cross or vector product for 2d and 3d vectors. 
  !+-----------------------------------------------------------------------------+!
  function cross_2d_d(a,b) result(c)
    real(8),dimension(2) :: a,b
    real(8)              :: c
    c = a(1)*b(2) - a(2)*b(1)
  end function cross_2d_d
  function cross_2d_c(a,b) result(c)
    complex(8),dimension(2) :: a,b
    complex(8)              :: c
    c = a(1)*b(2) - a(2)*b(1)
  end function cross_2d_c
  function cross_3d_d(a,b) result(c)
    real(8),dimension(3) :: a,b
    real(8),dimension(3) :: c
    c(1) = a(2)*b(3) - a(3)*b(2)
    c(2) = a(3)*b(1) - a(1)*b(3)
    c(3) = a(1)*b(2) - a(2)*b(1)
  end function cross_3d_d
  function cross_3d_c(a,b) result(c)
    complex(8),dimension(3) :: a,b
    complex(8),dimension(3) :: c
    c(1) = a(2)*b(3) - a(3)*b(2)
    c(2) = a(3)*b(1) - a(1)*b(3)
    c(3) = a(1)*b(2) - a(2)*b(1)
  end function cross_3d_c



  !+-----------------------------------------------------------------------------+!
  !PURPOSE: evaluate the S3 product A.(BxC) for 3d vectors
  !+-----------------------------------------------------------------------------+!
  function s3_product_d(a,b,c) result(s3)
    real(8),dimension(3),intent(in) :: a,b,c
    real(8)                         :: s3
    s3 = dot_product(a,cross_product(b, c))
  end function s3_product_d
  function s3_product_c(a,b,c) result(s3)
    complex(8),dimension(3),intent(in) :: a,b,c
    real(8)                            :: s3
    s3 = dot_product(a,cross_product(b, c))
  end function s3_product_c



  !-------------------------------------------------------------------------------------------
  !PURPOSE:  Eigenvalue/-vector problem for general matrices (in general non-symmetric or 
  ! non-complex-hermitian matrices).
  !-------------------------------------------------------------------------------------------
  !<TODO  add optional switch for left or right eigenvectors in deig() and zeig()?
  !>TODO
  subroutine deig(A,Eval,Evec,jobvl,jobvr)
    !IN matrix - OUT eigenvectors: A.evec = eval.evec ; evec(i,j) = ith component of jth vec.
    real(8),intent(in)                 :: A(:,:)
    complex(8),intent(out)             :: Eval(:)
    complex(8),intent(out)             :: Evec(:,:)
    character(len=1),optional          :: jobvl,jobvr
    character(len=1)                   :: jobvl_,jobvr_
    real(8),dimension(:,:),allocatable :: At,vl,vr
    real(8),dimension(:),allocatable   :: wi,wr
    real(8),dimension(:),allocatable   :: work
    real(8),dimension(1)               :: lwork_guess
    integer                            :: info,lda,ldvl,ldvr,lwork,n,i
    jobvl_='N';if(present(jobvl))jobvl_=jobvl
    jobvr_='V';if(present(jobvr))jobvr_=jobvr
    if(jobvl_=='V')stop "deig error: jobvl = V is not supported yet."
    lda   = size(A(:,1))
    n     = size(A(1,:))
    call assert_shape(A,[n,n],"solve","A")
    call assert_shape(Evec,[n,n],"solve","Evec")
    ldvl  = n
    ldvr  = n
    allocate(At(n,n),wr(n),wi(n),vl(ldvl,n),vr(ldvr,n))
    !Copy the input Matrix
    At    = A
    !1st Call: Query the right size for the working array.
    call dgeev(jobvl_, jobvr_, n, At, lda, wr, wi, vl, ldvl, vr, ldvr, lwork_guess, -1, info)
    if(info /= 0) then
       print*, "dgeev returned info = ",info
       if (info < 0) then
          print*, "the",-info,"-th argument had an illegal value"
       else
          print*, "the QR algorithm failed to compute all the"
          print*, "eigenvalues, and no eigenvectors have been computed;"
          print*, "elements ", info+1, ":", n, "of WR and WI contain eigenvalues which"
          print*, "have converged."
       end if
       stop 'deig error: 1st call dgeev'
    endif
    lwork = int(lwork_guess(1)) ;allocate(work(lwork))
    !2nd Call: Actual solution of the eigenproblem
    call dgeev(jobvl_, jobvr_, n, At, lda, wr, wi, vl, ldvl, vr, ldvr, work, lwork, info)
    if(info /= 0) then
       print *, "dgeev returned info = ",info
       if (info < 0) then
          print*, "the",-info,"-th argument had an illegal value"
       else
          print*, "the QR algorithm failed to compute all the"
          print*, "eigenvalues, and no eigenvectors have been computed;"
          print*, "elements ", info+1, ":", n, "of WR and WI contain eigenvalues which"
          print*, "have converged."
       end if
       stop 'deig error: 2nd call dgeev'
    end if
    Eval = wr + xi*wi
    ! as DGEEV has a rather complicated way of returning the eigenvectors,
    ! it is necessary to build the complex array of eigenvectors from
    ! two real arrays:
    do i = 1,n
       if(wi(i) > 0d0) then  ! first of two conjugate eigenvalues
          Evec(:, i) = vr(:, i) + xi*vr(:, i+1)
       elseif(wi(i) < 0d0) then  ! second of two conjugate eigenvalues
          Evec(:, i) = vr(:, i-1) - xi*vr(:, i)
       else
          Evec(:, i) = vr(:, i)
       end if
    end do
  end subroutine deig

  subroutine zeig(A,Eval,Evec,jobvl,jobvr)
    !IN matrix - OUT eigenvectors: A.evec = eval.evec ; evec(i,j) = ith component of jth vec.
    complex(8),intent(in)                 :: A(:,:)
    complex(8),intent(out)                :: Eval(:)
    complex(8),intent(out)                :: Evec(:,:)
    character(len=1),optional             :: jobvl,jobvr
    character(len=1)                      :: jobvl_,jobvr_
    real(8),dimension(:),allocatable      :: rwork
    complex(8),dimension(:),allocatable   :: work
    complex(8),dimension(1)               :: lwork_guess
    complex(8),dimension(:,:),allocatable :: vl,vr
    integer                               :: info,lda,ldvl,ldvr,lwork,n,i,lrwork

    jobvl_='N';if(present(jobvl))jobvl_=jobvl
    jobvr_='V';if(present(jobvr))jobvr_=jobvr
    if(jobvl_=='V')stop "deig error: jobvl = V is not supported yet."
    lda   = size(A(:,1))
    n     = size(A(1,:))
    call assert_shape(A,[n,n],"solve","A")
    call assert_shape(Evec,[n,n],"solve","Evec")
    ldvl  = n
    ldvr  = n
    lrwork= 2*n
    allocate(vl(ldvl,n),vr(ldvr,n),rwork(lrwork))
    !Copy the input Matrix
    Evec  = A
    !1st Call: Query the right size for the working array.
    call zgeev(jobvl_, jobvr_, n, Evec, lda, Eval, vl, ldvl, vr, ldvr, lwork_guess, -1, rwork, info)
    if(info /= 0) then
       print*, "zgeev returned info = ",info
       if (info < 0) then
          print*, "the",-info,"-th argument had an illegal value"
       else
          print*, "the QR algorithm failed to compute all the"
          print*, "eigenvalues, and no eigenvectors have been computed;"
          print*, "elements ", info+1, ":", n, "of WR and WI contain eigenvalues which"
          print*, "have converged."
       end if
       stop 'deig error: 1st call zgeev'
    endif
    lwork = int(lwork_guess(1)) ;allocate(work(lwork))
    !2nd Call: Actual solution of the eigenproblem
    call zgeev(jobvl_, jobvr_, n, Evec, lda, Eval, vl, ldvl, vr, ldvr, work, lwork, rwork, info)
    if(info /= 0) then
       print *, "zgeev returned info = ",info
       if (info < 0) then
          print*, "the",-info,"-th argument had an illegal value"
       else
          print*, "the QR algorithm failed to compute all the"
          print*, "eigenvalues, and no eigenvectors have been computed;"
          print*, "elements ", info+1, ":", n, "of WR and WI contain eigenvalues which"
          print*, "have converged."
       end if
       stop 'deig error: 2nd call zgeev'
    end if
    Evec = vr
  end subroutine zeig









  !-------------------------------------------------------------------------------------------
  !PURPOSE:  eigenvalue/-vector problem for real symmetric/complex hermitian matrices:
  !-------------------------------------------------------------------------------------------
  subroutine deigh_generalized(Am, Bm, lam, c)
    ! solves generalized eigen value problem for all eigenvalues and eigenvectors
    ! Am must by symmetric, Bm symmetric positive definite. ! Only the lower triangular part of Am and Bm is used.
    real(8), intent(in)  :: Am(:,:)   ! LHS matrix: Am c = lam Bm c
    real(8), intent(in)  :: Bm(:,:)   ! RHS matrix: Am c = lam Bm c
    real(8), intent(out) :: lam(:)   ! eigenvalues: Am c = lam Bm c
    real(8), intent(out) :: c(:,:)   ! eigenvectors: Am c = lam Bm c; c(i,j) = ith component of jth vec.
    integer              :: n
    ! lapack variables
    integer              :: lwork, liwork, info
    integer, allocatable :: iwork(:)
    real(8), allocatable :: Bmt(:,:), work(:)
    ! solve
    n = size(Am,1)
    call assert_shape(Am, [n, n], "eigh", "Am")
    call assert_shape(Bm, [n, n], "eigh", "B")
    call assert_shape(c, [n, n], "eigh", "c")
    lwork = 1 + 6*n + 2*n**2
    liwork = 3 + 5*n
    allocate(Bmt(n,n), work(lwork), iwork(liwork))
    c = Am; Bmt = Bm  ! Bmt temporaries overwritten by dsygvd
    call dsygvd(1,'V','L',n,c,n,Bmt,n,lam,work,lwork,iwork,liwork,info)
    if (info /= 0) then
       print *, "dsygvd returned info =", info
       if (info < 0) then
          print *, "the", -info, "-th argument had an illegal value"
       else if (info <= n) then
          print *, "the algorithm failed to compute an eigenvalue while working"
          print *, "on the submatrix lying in rows and columns", 1d0*info/(n+1)
          print *, "through", mod(info, n+1)
       else
          print *, "The leading minor of order ", info-n, &
               "of B is not positive definite. The factorization of B could ", &
               "not be completed and no eigenvalues or eigenvectors were computed."
       end if
       stop 'deigh_generalized error: dsygvd'
    end if
  end subroutine deigh_generalized
  !
  subroutine zeigh_generalized(Am, Bm, lam, c)
    ! solves generalized eigen value problem for all eigenvalues and eigenvectors
    ! Am must by hermitian, Bm hermitian positive definite.
    ! Only the lower triangular part of Am and Bm is used.
    complex(8), intent(in)  :: Am(:,:)   ! LHS matrix: Am c = lam Bm c
    complex(8), intent(in)  :: Bm(:,:)   ! RHS matrix: Am c = lam Bm c
    real(8), intent(out)    :: lam(:)      ! eigenvalues: Am c = lam Bm c
    complex(8), intent(out) :: c(:,:)   ! eigenvectors: Am c = lam Bm c; c(i,j) = ith component of jth vec.
    ! lapack variables
    integer                 :: info, liwork, lrwork, lwork, n
    integer, allocatable    :: iwork(:)
    real(8), allocatable    :: rwork(:)
    complex(8), allocatable :: Bmt(:,:), work(:)
    n = size(Am,1)
    call assert_shape(Am, [n, n], "eigh", "Am")
    call assert_shape(Bm, [n, n], "eigh", "Bm")
    call assert_shape(c, [n, n], "eigh", "c")
    lwork = 2*n + n**2
    lrwork = 1 + 5*N + 2*n**2
    liwork = 3 + 5*n
    allocate(Bmt(n,n), work(lwork), rwork(lrwork), iwork(liwork))
    c = Am; Bmt = Bm  ! Bmt temporary overwritten by zhegvd
    call zhegvd(1,'V','L',n,c,n,Bmt,n,lam,work,lwork,rwork,lrwork,iwork,liwork,info)
    if (info /= 0) then
       print *, "zhegvd returned info =", info
       if (info < 0) then
          print *, "the", -info, "-th argument had an illegal value"
       else if (info <= n) then
          print *, "the algorithm failed to compute an eigenvalue while working"
          print *, "on the submatrix lying in rows and columns", 1d0*info/(n+1)
          print *, "through", mod(info, n+1)
       else
          print *, "The leading minor of order ", info-n, &
               "of B is not positive definite. The factorization of B could ", &
               "not be completed and no eigenvalues or eigenvectors were computed."
       end if
       stop 'deigh_generalized error: zhegvd'
    end if
  end subroutine zeigh_generalized


  subroutine deigh_simple(M,E,jobz,uplo)
    real(8),dimension(:,:),intent(inout) :: M ! M v = E v/v(i,j) = ith component of jth vec.
    real(8),dimension(:),intent(inout)   :: E ! eigenvalues
    character(len=1),optional            :: jobz,uplo
    character(len=1)                     :: jobz_,uplo_
    integer                              :: i,j,n,lda,info,lwork
    real(8),dimension(:),allocatable     :: work
    real(8),dimension(1)                 :: lwork_guess
    jobz_='V';if(present(jobz))jobz_=jobz
    uplo_='U';if(present(uplo))uplo_=uplo
    lda = max(1,size(M,1))
    n   = size(M,2)
    call assert_shape(M,[n,n],"eigh","M")
    Call dsyev(jobz_,uplo_,n,M,lda,E,lwork_guess,-1,info)
    if (info /= 0) then
       print*, "dsyevd returned info =", info
       if (info < 0) then
          print*, "the", -info, "-th argument had an illegal value"
       else
          print*, "the algorithm failed to compute an eigenvalue while working"
          print*, "on the submatrix lying in rows and columns", 1d0*info/(n+1)
          print*, "through", mod(info, n+1)
       end if
       stop 'error deigh: 1st call dsyev'
    end if
    lwork=lwork_guess(1)
    allocate(work(lwork))
    call dsyev(jobz_,uplo_,n,M,lda,E,work,lwork,info)
    if (info /= 0) then
       print*, "dsyevd returned info =", info
       if (info < 0) then
          print*, "the", -info, "-th argument had an illegal value"
       else
          print*, "the algorithm failed to compute an eigenvalue while working"
          print*, "on the submatrix lying in rows and columns", 1d0*info/(n+1)
          print*, "through", mod(info, n+1)
       end if
       stop 'error deigh: 2nd call dsyev'
    end if
    deallocate(work)
  end subroutine deigh_simple
  !-----------------------------
  subroutine zeigh_simple(M,E,jobz,uplo)
    complex(8),dimension(:,:),intent(inout) :: M! M v = E v/v(i,j) = ith component of jth vec.
    real(8),dimension(:),intent(inout)      :: E! eigenvalues
    character(len=1),optional               :: jobz,uplo
    character(len=1)                        :: jobz_,uplo_
    integer                                 :: i,j,n,lda,info,lwork
    complex(8),dimension(1)                 :: lwork_guess
    complex(8),dimension(:),allocatable     :: work
    real(8),dimension(:),allocatable        :: rwork
    !write(*,*)"matrix_diagonalization called with: jobz="//jobz_//" uplo="//uplo_
    jobz_='V';if(present(jobz))jobz_=jobz
    uplo_='U';if(present(uplo))uplo_=uplo
    lda = max(1,size(M,1))
    n   = size(M,2)
    call assert_shape(M,[n,n],"eigh","M")
    allocate(rwork(max(1,3*N-2)))
    call zheev(jobz_,uplo_,n,M,lda,E,lwork_guess,-1,rwork,info)
    if(info/=0) then
       print*, "zheev returned info =", info
       if (info < 0) then
          print*, "the", -info, "-th argument had an illegal value"
       else
          print*, "the algorithm failed to compute an eigenvalue while working"
          print*, "on the submatrix lying in rows and columns", 1d0*info/(n+1)
          print*, "through", mod(info, n+1)
       end if
       stop 'error zeigh: 1st call zheev'
    end if
    lwork=lwork_guess(1) ; allocate(work(lwork))
    call zheev(jobz_,uplo_,n,M,lda,E,work,lwork,rwork,info)
    if(info/=0) then
       print*, "zheev returned info =", info
       if (info < 0) then
          print*, "the", -info, "-th argument had an illegal value"
       else
          print*, "the algorithm failed to compute an eigenvalue while working"
          print*, "on the submatrix lying in rows and columns", 1d0*info/(n+1)
          print*, "through", mod(info, n+1)
       end if
       stop 'error zeigh: 2nd call zheev'
    end if
    deallocate(work,rwork)
  end subroutine zeigh_simple





  !-------------------------------------------------------------------------------------------
  !PURPOSE:  comment
  !-------------------------------------------------------------------------------------------
  function deigvals(A) result(lam)
    real(8), intent(in)     :: A(:, :) ! matrix for eigenvalue compuation
    complex(8)              :: lam(size(A,1))  ! eigenvalues: A c = lam c
    real(8), allocatable    :: At(:,:),vl(:,: ),vr(:,:),wi(:),work(:),wr(:)
    real(8)                 :: lwork_guess(1)
    integer                 :: info, lda, ldvl, ldvr, lwork, n
    lda   = size(A(:,1))
    n     = size(A(1,:))
    call assert_shape(A,[n,n],"solve","A")
    ldvl  = n
    ldvr  = n
    allocate(At(n,n),wr(n),wi(n),vl(ldvl,n),vr(ldvr,n))
    !Copy the input Matrix
    At    = A
    !1st Call: Query the right size for the working array.
    call dgeev('N', 'N', n, At, lda, wr, wi, vl, ldvl, vr, ldvr,  lwork_guess, -1, info)
    if(info /= 0) then
       print*, "dgeev returned info = ",info
       if (info < 0) then
          print*, "the",-info,"-th argument had an illegal value"
       else
          print*, "the QR algorithm failed to compute all the"
          print*, "eigenvalues, and no eigenvectors have been computed;"
          print*, "elements ", info+1, ":", n, "of WR and WI contain eigenvalues which"
          print*, "have converged."
       end if
       stop 'deig error: 1st call dgeev'
    endif
    lwork = int(lwork_guess(1)) ;allocate(work(lwork))
    !2nd Call: Actual solution of the eigenproblem
    call dgeev('N', 'N', n, At, lda, wr, wi, vl, ldvl, vr, ldvr, work, lwork, info)
    if(info /= 0) then
       print *, "dgeev returned info = ",info
       if (info < 0) then
          print*, "the",-info,"-th argument had an illegal value"
       else
          print*, "the QR algorithm failed to compute all the"
          print*, "eigenvalues, and no eigenvectors have been computed;"
          print*, "elements ", info+1, ":", n, "of WR and WI contain eigenvalues which"
          print*, "have converged."
       end if
       stop 'deig error: 2nd call dgeev'
    end if
    lam = wr + xi*wi
  end function deigvals

  function zeigvals(A) result(lam)
    complex(8), intent(in)  :: A(:, :)  ! matrix to solve eigenproblem for
    complex(8)              :: lam(size(A,1))  ! eigenvalues: A c = lam c
    integer                 :: info, lda, ldvl, ldvr, lwork, n, lrwork
    real(8), allocatable    :: rwork(:)
    complex(8), allocatable :: At(:,:), vl(:,:), vr(:,:), work(:)
    complex(8)              :: lwork_guess(1)
    lda   = size(A(:,1))
    n     = size(A(1,:))
    call assert_shape(A,[n,n],"solve","A")
    ldvl  = n
    ldvr  = n
    lrwork= 2*n
    allocate(vl(ldvl,n),vr(ldvr,n),rwork(lrwork))
    !Copy the input Matrix
    At  = A
    !1st Call: Query the right size for the working array.
    call zgeev('N', 'N', n, At, lda, lam, vl, ldvl, vr, ldvr, lwork_guess, -1, rwork, info)
    if(info /= 0) then
       print*, "zgeev returned info = ",info
       if (info < 0) then
          print*, "the",-info,"-th argument had an illegal value"
       else
          print*, "the QR algorithm failed to compute all the"
          print*, "eigenvalues, and no eigenvectors have been computed;"
          print*, "elements ", info+1, ":", n, "of WR and WI contain eigenvalues which"
          print*, "have converged."
       end if
       stop 'deig error: 1st call zgeev'
    endif
    lwork = int(lwork_guess(1)) ;allocate(work(lwork))
    !2nd Call: Actual solution of the eigenproblem
    call zgeev('N', 'N', n, At, lda, lam, vl, ldvl, vr, ldvr, work, lwork, rwork, info)
    if(info /= 0) then
       print *, "zgeev returned info = ",info
       if (info < 0) then
          print*, "the",-info,"-th argument had an illegal value"
       else
          print*, "the QR algorithm failed to compute all the"
          print*, "eigenvalues, and no eigenvectors have been computed;"
          print*, "elements ", info+1, ":", n, "of WR and WI contain eigenvalues which"
          print*, "have converged."
       end if
       stop 'deig error: 2nd call zgeev'
    end if
  end function zeigvals




  !-------------------------------------------------------------------------------------------
  !PURPOSE:  comment
  !-------------------------------------------------------------------------------------------
  function deigvalsh(A) result(lam)
    real(8), intent(in)  :: A(:, :)         ! matrix for eigenvalue compuation
    real(8)              :: lam(size(A,1))  ! eigenvalues: A c = lam c
    real(8), allocatable :: At(:,:),work(:),iwork(:)
    real(8)              :: lwork_guess(1),liwork_guess(1)
    integer              :: info,lda,lwork,liwork,n
    lda   = size(A(:,1))
    n     = size(A(1,:))
    call assert_shape(A,[n,n],"solve","A")
    allocate(At(n,n))
    !Copy the input Matrix
    At    = A
    !1st Call: Query the right size for the working array.
    call dsyevd('N','U', n, At, lda, lam, lwork_guess, -1, liwork_guess, -1, info )
    if(info /= 0) then
       print*, "dsyevd returned info = ",info
       stop
    endif
    lwork  = int(lwork_guess(1))  ;allocate(work(lwork))
    liwork = int(liwork_guess(1)) ;allocate(iwork(liwork))
    !2nd Call: Actual solution of the eigenproblem
    call dsyevd('N','U', n, At, lda, lam, work, lwork, iwork, liwork, info )
    if(info /= 0) then
       print *, "ssyevd returned info = ",info
       if (info < 0) then
          print*, "the",-info,"-th argument had an illegal value"
       else
          print*,"the algorithm failed to converge"
          print*,"the",info,"off-diagonal elements of an intermediate"
          print*,"tridiagonal form did not converge to zero"
       end if
       stop 'deigvalsh error: 2nd call dsyevd'
    end if
  end function deigvalsh

  function zeigvalsh(A) result(lam)
    complex(8),intent(in) :: A(:, :)         ! matrix for eigenvalue compuation
    real(8)               :: lam(size(A,1))  ! eigenvalues: A c = lam c
    real(8), allocatable  :: At(:,:),work(:),rwork(:),iwork(:)
    real(8)               :: lwork_guess(1),lrwork_guess(1),liwork_guess(1)
    integer               :: info, lda,lwork,lrwork,liwork,n
    lda   = size(A(:,1))
    n     = size(A(1,:))
    call assert_shape(A,[n,n],"solve","A")
    allocate(At(n,n))
    !Copy the input Matrix
    At    = A
    !1st Call: Query the right size for the working array.
    call zheevd('N','U', n, At, lda, lam, lwork_guess,-1, lrwork_guess,-1, liwork_guess,-1, info )
    if(info /= 0) then
       print*, "zsyevd returned info = ",info
       stop
    endif
    lwork  = int(lwork_guess(1))  ;allocate(work(lwork))
    rwork  = int(lrwork_guess(1)) ;allocate(rwork(lrwork))
    liwork = int(liwork_guess(1)) ;allocate(iwork(liwork))
    !2nd Call: Actual solution of the eigenproblem
    call zheevd('N','U', n, At, lda, lam, work,lwork, rwork,lrwork, iwork,liwork, info )
    if(info /= 0) then
       print *, "zheevd returned info = ",info
       if (info < 0) then
          print*, "the",-info,"-th argument had an illegal value"
       else
          print*,"the algorithm failed to converge"
          print*,"the",info,"off-diagonal elements of an intermediate"
          print*,"tridiagonal form did not converge to zero"
       end if
       stop 'zeigvalsh error: 2nd call zheevd'
    end if
  end function zeigvalsh











  !-------------------------------------------------------------------------------------------
  !PURPOSE  : Invert a general m*n matrix using LAPACK library  
  ! M is destroyed and replaces by its inverse M^-1
  !-------------------------------------------------------------------------------------------
  subroutine dinv(Am)
    real(8), intent(inout) :: Am(:,:)            ! matrix to be inverted
    real(8), allocatable   :: Amt(:,:),work(:)  ! temporary work arrays
    integer                :: info,lda,n,lwork,nb
    integer, allocatable   :: ipiv(:)
    n = size(Am(1, :))
    call assert_shape(Am, [n, n], "inv", "Am")
    lda = n
    nb = ilaenv(1, 'DGETRI', "UN", n, -1, -1, -1)  ! TODO: check UN param
    lwork = n*nb
    if (nb < 1) nb = max(1, n)
    allocate(Amt(n,n), work(lwork), ipiv(n))
    Amt = Am
    call dgetrf(n, n, Amt, lda, ipiv, info)
    if(info /= 0) then
       print *, "dgetrf returned info =", info
       if (info < 0) then
          print *, "the", -info, "-th argument had an illegal value"
       else
          print *, "U(", info, ",", info, ") is exactly zero; The factorization"
          print *, "has been completed, but the factor U is exactly"
          print *, "singular, and division by zero will occur if it is used"
          print *, "to solve a system of equations."
       end if
       stop ' D_mat_invert: dgetrf error'
    end if
    call dgetri(n, Amt, n, ipiv, work, lwork, info)
    if (info /= 0) then
       print *, "dgetri returned info =", info
       if (info < 0) then
          print *, "the", -info, "-th argument had an illegal value"
       else
          print *, "U(", info, ",", info, ") is exactly zero; the matrix is"
          print *, "singular and its inverse could not be computed."
       end if
       stop ' D_mat_invert: dgetri error'
    end if
    Am = Amt
  end subroutine Dinv
  !
  subroutine Zinv(Am)
    ! Inverts the general complex matrix Am
    complex(8), intent(inout) :: Am(:,:)   ! Matrix to be inverted
    complex(8), allocatable   :: Amt(:,:), work(:)
    integer                   :: n, nb
    integer                   :: lwork, info
    integer, allocatable      :: ipiv(:)
    n = size(Am, 1)
    call assert_shape(Am, [n, n], "inv", "Am")
    nb = ilaenv(1, 'ZGETRI', "UN", n, -1, -1, -1)  ! TODO: check UN param
    if (nb < 1) nb = max(1, n)
    lwork = n*nb
    allocate(Amt(n,n), ipiv(n), work(lwork))
    Amt = Am
    call zgetrf(n, n, Amt, n, ipiv, info)
    if (info /= 0) then
       print *, "zgetrf returned info =", info
       if (info < 0) then
          print *, "the", -info, "-th argument had an illegal value"
       else
          print *, "U(", info, ",", info, ") is exactly zero; The factorization"
          print *, "has been completed, but the factor U is exactly"
          print *, "singular, and division by zero will occur if it is used"
          print *, "to solve a system of equations."
       end if
       stop 'Z_mat_invert: zgetrf error'
    end if
    call zgetri(n, Amt, n, ipiv, work, lwork, info)
    if (info /= 0) then
       print *, "zgetri returned info =", info
       if (info < 0) then
          print *, "the", -info, "-th argument had an illegal value"
       else
          print *, "U(", info, ",", info, ") is exactly zero; the matrix is"
          print *, "singular and its inverse could not be computed."
       end if
       stop 'Z_mat_invert: zgetri error'
    end if
    Am = Amt
  end subroutine Zinv



  !-------------------------------------------------------------------------------------------
  !PURPOSE  : Invert a general triangular m*n matrix using LAPACK library 
  ! M is destroyed and replaces by its inverse M^-1
  ! M on output is the square matrix n*n
  !-------------------------------------------------------------------------------------------
  subroutine Dinv_triang(A,uplo,diag)
    real(8),dimension(:,:)           :: A
    character(len=1),optional        :: uplo,diag
    character(len=1)                 :: uplo_
    character(len=1)                 :: diag_
    integer                          :: n,lda,info
    uplo_="U";if(present(uplo))uplo_=uplo
    diag_="N";if(present(diag))diag_=diag !not a unit triangular matrix
    lda = max(1,size(A,1))
    n   = size(A,2)
    call dtrtri(uplo_,diag_,n,A,lda,info)
    if(info/=0)stop "Error MATRIX/D_mat_invertTRIANG: dtrtri"
  end subroutine Dinv_triang
  !
  subroutine Zinv_triang(A,uplo,diag)
    complex(8),dimension(:,:)           :: A
    character(len=1),optional        :: uplo,diag
    character(len=1)                 :: uplo_
    character(len=1)                 :: diag_
    integer                          :: ndim1,ndim2
    integer                          :: n,lda,info
    uplo_="U";if(present(uplo))uplo_=uplo
    diag_="N";if(present(diag))diag_=diag !not a unit triangular matrix
    lda = max(1,size(A,1))
    n   = size(A,2)
    call ztrtri(uplo_,diag_,n,A,lda,info)
    if(info/=0)stop "Error MATRIX/D_mat_invertTRIANG: ztrtri"
  end subroutine Zinv_triang





  !-------------------------------------------------------------------------------------------
  !PURPOSE  : Invert a symmetric n*n matrix using LAPACK library 
  ! M is destroyed and replaces by its inverse M^-1
  !-------------------------------------------------------------------------------------------
  subroutine Dinv_sym(A,uplo)
    real(8),dimension(:,:)           :: A
    character(len=*),optional        :: uplo
    character(len=1)                 :: uplo_
    integer                          :: n,lda,info,lwork
    integer,dimension(:),allocatable :: ipvt
    real(8),dimension(:),allocatable :: work
    real(8),dimension(1)             :: lwork_guess
    uplo_="U";if(present(uplo))uplo_=uplo
    lda  = max(1,size(A,1))
    n    = size(A,2)
    allocate(ipvt(n))
    call dsytrf(uplo_,n,A,lda,ipvt,lwork_guess,-1,info)
    if(info/=0)stop "Error MATRIX/D_mat_invertSYM: 1st call dsytrf"
    lwork=lwork_guess(1);allocate(work(lwork))
    call dsytrf(uplo_,n,A,lda,ipvt,work,lwork,info)
    if(info/=0)stop "Error MATRIX/D_mat_invertSYM: 2nd call dsytrf"
    call dsytri(uplo_,n,A,lda,ipvt,work,info)
    if(info/=0)stop "Error MATRIX/D_mat_invertSYM: dsytri"
    deallocate(ipvt,work)
  end subroutine Dinv_SYM
  !-----------------------------
  subroutine Zinv_sym(A,uplo)
    complex(8),dimension(:,:)           :: A
    character(len=*),optional        :: uplo
    character(len=1)                 :: uplo_
    integer                          :: n,lda,info,lwork
    integer,dimension(:),allocatable :: ipvt
    complex(8),dimension(:),allocatable :: work
    complex(8),dimension(1)             :: lwork_guess
    uplo_="U";if(present(uplo))uplo_=uplo
    lda  = max(1,size(A,1))
    n    = size(A,2)
    allocate(ipvt(n))
    call zsytrf(uplo_,n,A,lda,ipvt,lwork_guess,-1,info)
    if(info/=0)stop "Error MATRIX/D_mat_invertSYM: 1st call zsytrf"
    lwork=lwork_guess(1);allocate(work(lwork))
    call zsytrf(uplo_,n,A,lda,ipvt,work,lwork,info)
    if(info/=0)stop "Error MATRIX/D_mat_invertSYM: 2nd call zsytrf"
    call zsytri(uplo_,n,A,lda,ipvt,work,info)
    if(info/=0)stop "Error MATRIX/D_mat_invertSYM: zsytri"
    deallocate(ipvt,work)
  end subroutine Zinv_sym



  !-------------------------------------------------------------------------------------------
  !PURPOSE  : Invert a hermitian n*n matrix using MKL_LAPACK library  
  ! M is destroyed and replaces by its inverse M^-1
  !-------------------------------------------------------------------------------------------
  subroutine zinv_her(A,uplo)
    complex(8),dimension(:,:)                 :: A
    complex(8)                                :: c
    character(len=*),optional                 :: uplo
    character(len=1)                          :: uplo_
    integer                                   :: n,lda,info,lwork,i,j
    integer,dimension(:),allocatable          :: ipvt
    complex(8),dimension(:),allocatable       :: work
    complex(8),dimension(1)                   :: lwork_guess
    logical                                   :: bool
    uplo_="U";if(present(uplo))uplo_=uplo
    lda=max(1,size(A,1))
    n  =size(A,2)
    allocate(ipvt(n))
    !Test hermiticity:
    bool=.false.
    testH: do i=1,size(A,1)
       do j=1,size(A,2)
          c=A(i,j)-conjg(A(j,i))
          if(c/=cmplx(0.d0,0.d0,8))bool=.true.
          exit testH
       enddo
    enddo testH
    if(bool)stop "Error MATRIX/Z_mat_invertHER: A not Hermitian"
    !
    call zhetrf(uplo_,n,A,lda,ipvt,lwork_guess,-1,info)
    if(info/=0)stop "Error MATRIX/Z_mat_invertHER: 1st call zhetrf"
    lwork=lwork_guess(1) ; allocate(work(lwork))
    call zhetrf(uplo_,n,A,lda,ipvt,work,lwork,info)
    if(info/=0)stop "Error MATRIX/Z_mat_invertHERE: 2nd call zhetrf"
    call zhetri(uplo_,n,A,lda,ipvt,work,info)
    if(info/=0)stop "Error MATRIX/Z_mat_invertHERE: zhetri"
    deallocate(ipvt,work)
    if(uplo_=="U")then
       forall(i=1:size(A,1),j=1:size(A,2),i>j)A(i,j)=conjg(A(j,i))
    elseif(uplo_=="L")then
       forall(i=1:size(A,1),j=1:size(A,2),i<j)A(i,j)=conjg(A(j,i))
    endif
  end subroutine Zinv_her







  !+-----------------------------------------------------------------+
  !PURPOSE  : Linear equation solution by Gauss-Jordan elimination
  ! a is an N x N input coefficient matrix. On output, a is replaced 
  !by its matrix inverse.
  !+-----------------------------------------------------------------+
  subroutine Dinv_gj(a)
    real(8), dimension(:,:), intent(inout) :: a
    integer, dimension(size(a,1))      :: ipiv,indxr,indxc
    !these arrays are used for bookkeeping on the pivoting.
    !integer                            :: nn
    logical, dimension(size(a,1))      :: lpiv
    real(8)                            :: pivinv
    real(8), dimension(size(a,1))      :: dumc
    integer, target                    :: irc(2)
    integer                            :: i,l,n
    integer, pointer                   :: irow,icol
    real(8)                            :: zero=0.d0, one=1.d0
    n=size(a,1)
    irow => irc(1)
    icol => irc(2)
    ipiv=0
    do i=1,n
       !main loop over columns to be reduced.
       lpiv = (ipiv == 0)
       !begin search for a pivot element.
       irc=maxloc(abs(a),outerand(lpiv,lpiv))
       ipiv(icol)=ipiv(icol)+1
       if (ipiv(icol) > 1) stop 'gaussj:singular matrix (1)'
       !we now have the pivot element, so we interchange
       !rows, if needed, to put the pivot element on the diagonal. the columns
       !are not physically interchanged, only relabeled:
       !indxc(i),the column of the ith pivot element, is the ith column that is
       !reduced, while indxr(i) is the row in which that pivot element was
       !originally located. if indxr(i) = indxc(i) there is an implied column
       !interchange. with this form of bookkeeping, the inverse matrix will be
       !scrambled by
       !columns.
       if (irow /= icol) call swap(a(irow,:),a(icol,:))
       indxr(i)=irow !we are now ready to divide the pivot row by the pivot element,
       !located at irow and icol.
       indxc(i)=icol
       if (a(icol,icol) == zero) stop 'gaussj:singular matrix (2)'
       pivinv=one/a(icol,icol)
       a(icol,icol)= one !cmplx(one,zero)
       a(icol,:)=a(icol,:)*pivinv
       dumc=a(:,icol)
       !next, we reduce the rows, except for the pivot one, of course.
       a(:,icol)     = zero !cmplx
       a(icol,icol)  = pivinv
       a(1:icol-1,:) = a(1:icol-1,:) - outerprod(dumc(1:icol-1),a(icol,:))
       a(icol+1:,:)  = a(icol+1:,:)  - outerprod(dumc(icol+1:),a(icol,:))
    end do
    !it only remains to unscramble the solution in view of the column
    !interchanges.
    !we do this by interchanging pairs of columns in the reverse order that the
    !permutation
    !was built up.
    do l=n,1,-1
       call swap(a(:,indxr(l)),a(:,indxc(l)))
    end do
  contains
    function outerand(a,b)
      implicit none
      logical, dimension(:), intent(in)   :: a,b
      logical, dimension(size(a),size(b)) :: outerand
      outerand = spread(a,dim=2,ncopies=size(b)).and.spread(b,dim=1,ncopies=size(a))
    end function outerand
    function outerprod(a,b)
      real(8), dimension(:), intent(in) :: a,b
      real(8), dimension(size(a),size(b)) :: outerprod
      outerprod = spread(a,dim=2,ncopies=size(b)) * &
           spread(b,dim=1,ncopies=size(a))
    end function outerprod
  end subroutine Dinv_gj

  subroutine Zinv_gj(a)
    complex(8), dimension(:,:), intent(inout) :: a
    integer, dimension(size(a,1))      :: ipiv,indxr,indxc
    !these arrays are used for bookkeeping on the pivoting.
    !integer                            :: nn
    logical, dimension(size(a,1))      :: lpiv
    complex(8)                         :: pivinv
    complex(8), dimension(size(a,1))   :: dumc
    integer, target                    :: irc(2)
    integer                            :: i,l,n
    integer, pointer                   :: irow,icol
    real(8)                            :: zero=0.d0, one=1.d0
    n=size(a,1)
    irow => irc(1)
    icol => irc(2)
    ipiv=0
    do i=1,n
       !main loop over columns to be reduced.
       lpiv = (ipiv == 0)
       !begin search for a pivot element.
       irc=maxloc(abs(a),outerand(lpiv,lpiv))
       ipiv(icol)=ipiv(icol)+1
       if (ipiv(icol) > 1) stop 'gaussj:singular matrix (1)'
       !we now have the pivot element, so we interchange
       !rows, if needed, to put the pivot element on the diagonal. the columns
       !are not physically interchanged, only relabeled:
       !indxc(i),the column of the ith pivot element, is the ith column that is
       !reduced, while indxr(i) is the row in which that pivot element was
       !originally located. if indxr(i) = indxc(i) there is an implied column
       !interchange. with this form of bookkeeping, the inverse matrix will be
       !scrambled by
       !columns.
       if (irow /= icol) call swap(a(irow,:),a(icol,:))
       indxr(i)=irow !we are now ready to divide the pivot row by the pivot element,
       !located at irow and icol.
       indxc(i)=icol
       if (a(icol,icol) == zero) stop 'gaussj:singular matrix (2)'
       pivinv=one/a(icol,icol)
       a(icol,icol)=cmplx(one,zero,8)
       a(icol,:)=a(icol,:)*pivinv
       dumc=a(:,icol)
       !next, we reduce the rows, except for the pivot one, of course.
       a(:,icol)     = cmplx(zero,zero,8)
       a(icol,icol)  = pivinv
       a(1:icol-1,:) = a(1:icol-1,:) - outerprod(dumc(1:icol-1),a(icol,:))
       a(icol+1:,:)  = a(icol+1:,:)  - outerprod(dumc(icol+1:),a(icol,:))
    end do
    !it only remains to unscramble the solution in view of the column
    !interchanges.
    !we do this by interchanging pairs of columns in the reverse order that the
    !permutation
    !was built up.
    do l=n,1,-1
       call swap(a(:,indxr(l)),a(:,indxc(l)))
    end do
  contains
    function outerand(a,b)
      implicit none
      logical, dimension(:), intent(in)   :: a,b
      logical, dimension(size(a),size(b)) :: outerand
      outerand = spread(a,dim=2,ncopies=size(b)).and.spread(b,dim=1,ncopies=size(a))
    end function outerand
    function outerprod(a,b)
      complex(8), dimension(:), intent(in) :: a,b
      complex(8), dimension(size(a),size(b)) :: outerprod
      outerprod = spread(a,dim=2,ncopies=size(b)) * &
           spread(b,dim=1,ncopies=size(a))
    end function outerprod
  end subroutine Zinv_gj







  !+-----------------------------------------------------------------+
  !PROGRAM  : SOLVE_LINEAR_SYSTEM
  !PURPOSE  : Interface to lapack solution of linear system: A*x=b
  !+-----------------------------------------------------------------+
  subroutine Dsolve_1rhs(A,b,trans)
    real(8),dimension(:,:),intent(in)  :: A
    real(8),dimension(:),intent(inout) :: b
    real(8),dimension(:,:),allocatable :: b_
    character(len=1),optional          :: trans
    character(len=1)                   :: trans_
    integer                            :: m,n,nrhs,lda,ldb
    integer                            :: info
    integer,dimension(:),allocatable   :: ipvt
    trans_="N";if(present(trans))trans_=trans
    lda = max(1,size(A,1))
    ldb = max(1,size(B))
    m   = size(A,1)
    n   = size(A,2)
    nrhs= 1
    allocate(ipvt(min(m,n)))
    call dgetrf(m,n,A,lda,ipvt,info)
    if(info/=0)stop "Error MATRIX/d_mat_solve_linear_system: dgetrf"    
    allocate(b_(ldb,nrhs))
    b_(:,1)=b
    call dgetrs(trans_,n,nrhs,A,lda,ipvt,b_,ldb,info)
    if(info/=0)stop "Error MATRIX/d_mat_solve_linear_system: dgetrs"
    b=b_(:,1)
    deallocate(ipvt,b_)
  end subroutine dsolve_1rhs

  subroutine Zsolve_1rhs(A,b,trans)
    complex(8),dimension(:,:),intent(in)  :: A
    complex(8),dimension(:),intent(inout) :: b
    complex(8),dimension(:,:),allocatable :: b_
    character(len=1),optional          :: trans
    character(len=1)                   :: trans_
    integer                            :: m,n,nrhs,lda,ldb
    integer                            :: info
    integer,dimension(:),allocatable   :: ipvt
    trans_="N";if(present(trans))trans_=trans
    lda = max(1,size(A,1))
    ldb = max(1,size(B))
    m   = size(A,1)
    n   = size(A,2)
    nrhs= 1
    allocate(ipvt(min(m,n)))
    call zgetrf(m,n,A,lda,ipvt,info)
    if(info/=0)stop "Error MATRIX/d_mat_solve_linear_system: dgetrf"    
    allocate(b_(ldb,nrhs))
    b_(:,1)=b
    call zgetrs(trans_,n,nrhs,A,lda,ipvt,b_,ldb,info)
    if(info/=0)stop "Error MATRIX/d_mat_solve_linear_system: dgetrs"
    b=b_(:,1)
    deallocate(ipvt,b_)
  end subroutine Zsolve_1rhs

  subroutine Dsolve_Mrhs(A,b,trans)
    real(8),dimension(:,:),intent(in)    :: A
    real(8),dimension(:,:),intent(inout) :: b
    character(len=1),optional          :: trans
    character(len=1)                   :: trans_
    integer                            :: m,n,nrhs,lda,ldb
    integer                            :: info
    integer,dimension(:),allocatable   :: ipvt
    trans_="N";if(present(trans))trans_=trans
    lda = max(1,size(A,1))
    ldb = max(1,size(B,1))
    m   = size(A,1)
    n   = size(A,2)
    nrhs= size(B,2)
    allocate(ipvt(min(m,n)))
    call dgetrf(m,n,A,lda,ipvt,info)
    if(info/=0)stop "Error MATRIX/d_mat_solve_linear_system: dgetrf"    
    call dgetrs(trans_,n,nrhs,A,lda,ipvt,b,ldb,info)
    if(info/=0)stop "Error MATRIX/d_mat_solve_linear_system: dgetrs"
    deallocate(ipvt)
  end subroutine dsolve_Mrhs

  subroutine Zsolve_Mrhs(A,b,trans)
    complex(8),dimension(:,:),intent(in)    :: A
    complex(8),dimension(:,:),intent(inout) :: b
    character(len=1),optional          :: trans
    character(len=1)                   :: trans_
    integer                            :: m,n,nrhs,lda,ldb
    integer                            :: info
    integer,dimension(:),allocatable   :: ipvt
    trans_="N";if(present(trans))trans_=trans
    lda = max(1,size(A,1))
    ldb = max(1,size(B,1))
    m   = size(A,1)
    n   = size(A,2)
    nrhs= size(B,2)   
    allocate(ipvt(n))
    call zgetrf(m,n,A,lda,ipvt,info)
    if(info/=0)stop "Error MATRIX/d_mat_solve_linear_system: dgetrf"    
    lda=n ; ldb=n ; nrhs=size(b,2)
    call zgetrs(trans_,n,nrhs,A,lda,ipvt,b,ldb,info)
    if(info/=0)stop "Error MATRIX/d_mat_solve_linear_system: dgetrs"
    deallocate(ipvt)
  end subroutine Zsolve_Mrhs








  !-------------------------------------------------------------------------------------------
  !PURPOSE: compute the determinant of a real matrix using an LU factorization
  !-------------------------------------------------------------------------------------------
  function ddet(A) result(x)
    ! compute the determinant of a real matrix using an LU factorization
    real(8), intent(in)  :: A(:, :)
    real(8)              :: x
    integer              :: i
    integer              :: info, n
    integer, allocatable :: ipiv(:)
    real(8), allocatable :: At(:,:)
    n = size(A(1,:))
    call assert_shape(A, [n, n], "det", "A")
    allocate(At(n,n), ipiv(n))
    At = A
    call dgetrf(n, n, At, n, ipiv, info)
    if(info /= 0) then
       print *, "dgetrf returned info =", info
       if (info < 0) then
          print *, "the", -info, "-th argument had an illegal value"
       else
          print *, "U(", info, ",", info, ") is exactly zero; The factorization"
          print *, "has been completed, but the factor U is exactly"
          print *, "singular, and division by zero will occur if it is used"
          print *, "to solve a system of equations."
       end if
       stop 'det: dgetrf error'
    end if
    ! At now contains the LU of the factorization A = PLU
    ! as L has unit diagonal entries, the determinant can be computed
    ! from the product of U's diagonal entries. Additional sign changes
    ! stemming from the permutations P have to be taken into account as well.
    x = 1d0
    do i = 1,n
       if(ipiv(i) /= i) then  ! additional sign change
          x = -x*At(i,i)
       else
          x = x*At(i,i)
       endif
    end do
  end function ddet

  function zdet(A) result(x)
    ! compute the determinant of a real matrix using an LU factorization
    complex(8), intent(in)  :: A(:, :)
    complex(8)              :: x
    integer                 :: i
    integer                 :: info, n
    integer, allocatable    :: ipiv(:)
    complex(8), allocatable :: At(:,:)
    n = size(A(1,:))
    call assert_shape(A, [n, n], "det", "A")
    allocate(At(n,n), ipiv(n))
    At = A
    call zgetrf(n, n, At, n, ipiv, info)
    if(info /= 0) then
       print *, "zgetrf returned info =", info
       if (info < 0) then
          print *, "the", -info, "-th argument had an illegal value"
       else
          print *, "U(", info, ",", info, ") is exactly zero; The factorization"
          print *, "has been completed, but the factor U is exactly"
          print *, "singular, and division by zero will occur if it is used"
          print *, "to solve a system of equations."
       end if
       stop 'zdet error: zgetrf '
    end if
    ! for details on the computation, compare the comment in ddet().
    x = one
    do i = 1,n
       if(ipiv(i) /= i) then  ! additional sign change
          x = -x*At(i,i)
       else
          x = x*At(i,i)
       endif
    end do
  end function zdet






  !-------------------------------------------------------------------------------------------
  !PURPOSE: compute least square solution to A x = b for real A, b
  !-------------------------------------------------------------------------------------------
  function dlstsq(A, b) result(x)
    real(8), intent(in)  :: A(:,:), b(:)
    real(8), allocatable :: x(:)
    integer              :: info, ldb, lwork, m, n, rank
    real(8)              :: rcond
    real(8), allocatable :: work(:), At(:,:), Bt(:,:)
    integer, allocatable :: jpvt(:)
    m = size(A(:,1)) ! = lda
    n = size(A(1,:))
    ldb = size(b)
    allocate(x(n), At(m,n), Bt(ldb,1), jpvt(n), work(1))
    call dgelsy(m, n, 1, At, m, Bt, ldb, jpvt, rcond, rank, work, &
         -1, info)  ! query optimal workspace size
    lwork = int(real(work(1)))
    deallocate(work)
    allocate(work(lwork))  ! allocate with ideal size
    rcond = 0d0
    jpvt(:) = 0
    Bt(:,1) = b(:)  ! only one right-hand side
    At(:,:) = A(:,:)
    call dgelsy(m, n, 1, At, m, Bt, ldb, jpvt, rcond, rank, work, &
         lwork, info)
    if(info /= 0) then
       print *, "dgelsy returned info = ", info
       print *, "the ", -info, "-th argument had an illegal value"
       stop 'dlstsq error: dgelsy'
    endif
    x(:) = Bt(1:n,1)
  end function dlstsq

  function zlstsq(A, b) result(x)
    complex(8), intent(in)  :: A(:,:), b(:)
    complex(8), allocatable :: x(:)
    integer                 :: info, ldb, lwork, m, n, rank
    real(8)                 :: rcond
    complex(8), allocatable :: At(:,:), Bt(:,:), work(:)
    real(8), allocatable    :: rwork(:)
    integer, allocatable    :: jpvt(:)
    m = size(A(:,1)) ! = lda
    n = size(A(1,:))
    ldb = size(b)
    allocate(x(n), At(m,n), Bt(ldb,1), jpvt(n), work(1), rwork(2*n))
    call zgelsy(m, n, 1, At, m, Bt, ldb, jpvt, rcond, rank, work, &
         -1, rwork, info)  ! query optimal workspace size
    lwork = int(real(work(1)))
    deallocate(work)
    allocate(work(lwork))  ! allocate with ideal size
    rcond = 0d0
    jpvt(:) = 0
    Bt(:,1) = b(:)  ! only one right-hand side
    At(:,:) = A(:,:)
    call zgelsy(m, n, 1, At, m, Bt, ldb, jpvt, rcond, rank, work, &
         lwork, rwork, info)
    if(info /= 0) then
       print *, "zgelsy returned info = ", info
       print *, "the ", -info, "-th argument had an illegal value"
       stop 'zlstsq error: zgelsy'
    endif
    x(:) = Bt(1:n,1)
  end function zlstsq




  !-------------------------------------------------------------------------------------------
  !PURPOSE:  compute singular values s_i of a real/complex m x n matrix A
  !-------------------------------------------------------------------------------------------
  function dsvdvals(A) result(s)
    real(8), intent(in)  :: A(:,:)
    real(8), allocatable :: s(:)
    integer              :: info, lwork, m, n
    real(8), allocatable :: work(:), At(:,:)
    real(8)              :: u(1,1), vt(1,1)  ! not used if only s is to be computed
    m = size(A(:,1))  ! = lda
    n = size(A(1,:))
    allocate(At(m,n), s(min(m,n)))
    At(:,:) = A(:, :)  ! A is overwritten in dgesvd
    ! query optimal lwork and allocate workspace:
    allocate(work(1))
    call dgesvd('N', 'N', m, n, At, m, s, u, 1, vt, 1, work, -1, info)
    lwork = int(real(work(1)))
    deallocate(work)
    allocate(work(lwork))
    call dgesvd('N', 'N', m, n, At, m, s, u, 1, vt, 1, work, lwork, info)
    if(info /= 0) then
       print *, "dgesvd returned info = ", info
       if(info < 0) then
          print *, "the ", -info, "-th argument had an illegal value"
       else
          print *, "DBDSQR did not converge, there are ", info
          print *, "superdiagonals of an intermediate bidiagonal form B"
          print *, "did not converge to zero. See the description of WORK"
          print *, "in DGESVD's man page for details."
       endif
       stop 'dsvdvals error: dgesvd '
    endif
  end function dsvdvals

  function zsvdvals(A) result(s)
    complex(8), intent(in)  :: A(:,:)
    real(8), allocatable    :: s(:)
    integer                 :: info, lwork, m, n, lrwork
    complex(8), allocatable :: work(:), At(:,:)
    real(8), allocatable    :: rwork(:)
    complex(8)              :: u(1,1), vt(1,1)  ! not used if only s is to be computed
    m = size(A(:,1))  ! = lda
    n = size(A(1,:))
    lrwork = 5*min(m,n)
    allocate(At(m,n), s(min(m,n)), rwork(lrwork))
    At(:,:) = A(:,:)  ! A is overwritten in zgesvd!
    ! query optimal lwork and allocate workspace:
    allocate(work(1))
    call zgesvd('N', 'N', m, n, At, m, s, u, 1, vt, 1, work, -1, rwork, info)
    lwork = int(real(work(1)))
    deallocate(work)
    allocate(work(lwork))
    call zgesvd('N', 'N', m, n, At, m, s, u, 1, vt, 1, work, lwork, rwork, info)
    if(info /= 0) then
       print *, "zgesvd returned info = ", info
       if(info < 0) then
          print *, "the ", -info, "-th argument had an illegal value"
       else
          print *, "ZBDSQR did not converge, there are ", info
          print *, "superdiagonals of an intermediate bidiagonal form B"
          print *, "did not converge to zero. See the description of RWORK"
          print *, "in ZGESVD's man page for details."
       endif
       stop 'zsvdvals error: zgesvd'
    endif
  end function zsvdvals





  !-------------------------------------------------------------------------------------------
  !PURPOSE: compute the singular value decomposition A = U sigma Vtransp / U sigma V^H of 
  ! real/complex matrix A
  !-------------------------------------------------------------------------------------------
  subroutine dsvd(A, s, U, Vtransp)
    ! real m x n matrix A
    ! U is m x m
    ! Vtransp is n x n
    ! s has size min(m, n) --> sigma matrix is (n x m) with sigma_ii = s_i
    real(8), intent(in)  :: A(:,:)
    real(8), intent(out) :: s(:), U(:,:), Vtransp(:,:)
    integer              :: info, lwork, m, n, ldu
    real(8), allocatable :: work(:), At(:,:)
    m = size(A(:,1))  ! = lda
    n = size(A(1,:))
    ldu = m
    allocate(At(m,n))
    At(:,:) = A(:,:)  ! use a temporary as dgesvd destroys its input
    call assert_shape(U, [m, m], "svd", "U")
    call assert_shape(Vtransp, [n, n], "svd", "Vtransp")
    ! query optimal lwork and allocate workspace:
    allocate(work(1))
    call dgesvd('A', 'A', m, n, At, m, s, U, ldu, Vtransp, n, work, -1, info)
    lwork = int(real(work(1)))
    deallocate(work)
    allocate(work(lwork))
    call dgesvd('A', 'A', m, n, At, m, s, U, ldu, Vtransp, n, work, lwork, info)
    if(info /= 0) then
       print *, "dgesvd returned info = ", info
       if(info < 0) then
          print *, "the ", -info, "-th argument had an illegal value"
       else
          print *, "DBDSQR did not converge, there are ", info
          print *, "superdiagonals of an intermediate bidiagonal form B"
          print *, "did not converge to zero. See the description of WORK"
          print *, "in DGESVD's man page for details."
       endif
       stop 'dsvd errpr: dgesvd'
    endif
  end subroutine dsvd

  subroutine zsvd(A, s, U, Vtransp)
    ! complex m x m matrix A
    ! U is m x min(m, n)
    ! Vtransp is n x n
    ! sigma is m x n with with sigma_ii = s_i
    ! note that this routine returns V^H, not V!
    complex(8), intent(in)  :: A(:,:)
    real(8), intent(out)    :: s(:)
    complex(8), intent(out) :: U(:,:), Vtransp(:,:)
    integer                 :: info, lwork, m, n, ldu, lrwork
    real(8), allocatable    :: rwork(:)
    complex(8), allocatable :: work(:), At(:,:)
    ! TODO: check shapes here and in other routines?
    m = size(A(:,1))  ! = lda
    n = size(A(1,:))
    ldu = m
    lrwork = 5*min(m,n)
    allocate(rwork(lrwork), At(m,n))
    At(:,:) = A(:,:)  ! use a temporary as zgesvd destroys its input
    call assert_shape(U, [m, m], "svd", "U")
    call assert_shape(Vtransp, [n, n], "svd", "Vtransp")
    ! query optimal lwork and allocate workspace:
    allocate(work(1))
    call zgesvd('A', 'A', m, n, At, m, s, U, ldu, Vtransp, n, work, -1,&
         rwork, info)
    lwork = int(real(work(1)))
    deallocate(work)
    allocate(work(lwork))
    call zgesvd('A', 'A', m, n, At, m, s, U, ldu, Vtransp, n, work, &
         lwork, rwork, info)
    if(info /= 0) then
       print *, "zgesvd returned info = ", info
       if(info < 0) then
          print *, "the ", -info, "-th argument had an illegal value"
       else
          print *, "ZBDSQR did not converge, there are ", info
          print *, "superdiagonals of an intermediate bidiagonal form B"
          print *, "did not converge to zero. See the description of WORK"
          print *, "in DGESVD's man page for details."
       endif
       stop 'zsvd error: zgesvd'
    endif
  end subroutine zsvd





  !******************************************************************
  !                TRIDIAGONAL MATRICES ROUTINES
  !******************************************************************
  !+-----------------------------------------------------------------------------+!
  !PURPOSE:  comment
  !+-----------------------------------------------------------------------------+!
  function d_check_tridiag(Amat) result(Mcheck)
    real(8),dimension(:,:)                       :: Amat
    logical,dimension(size(Amat,1),size(Amat,2)) :: Lmat
    logical                                      :: Mcheck
    integer                                      :: i,j,N
    N=size(Amat,1)
    call assert_shape(Amat,[N,N],"d_check_tridiag","Amat")
    Lmat=.true.
    forall(i=1:N-1)
       Lmat(i+1,i)=.false.
       Lmat(i,i)  =.false.
       Lmat(i,i+1)=.false.
    end forall
    Lmat(N,N)=.false.
    Mcheck = .not.(sum(abs(Amat),mask=Lmat)>0d0)
  end function d_check_tridiag
  function c_check_tridiag(Amat) result(Mcheck)
    complex(8),dimension(:,:)                    :: Amat
    logical,dimension(size(Amat,1),size(Amat,2)) :: Lmat
    logical                                      :: Mcheck
    integer                                      :: i,N
    N=size(Amat,1)
    call assert_shape(Amat,[N,N],"c_check_tridiag","Amat")
    Lmat=.true.
    forall(i=1:N-1)
       Lmat(i+1,i)=.false.
       Lmat(i,i)  =.false.
       Lmat(i,i+1)=.false.
    end forall
    Lmat(N,N)=.false.
    Mcheck = .not.(sum(abs(Amat),mask=Lmat)>0d0)
  end function c_check_tridiag
  function d_check_tridiag_block(Nblock,Nsize,Amat) result(Mcheck)
    integer                                          :: Nblock
    integer                                          :: Nsize
    real(8),dimension(Nblock*Nsize,Nblock*Nsize)     :: Amat
    logical,dimension(Nblock*Nsize,Nblock*Nsize)     :: Lmat
    integer                                          :: i,j,iblock,is,js
    logical                                          :: Mcheck
    Lmat=.true.
    do iblock=1,Nblock-1
       do i=1,Nsize
          do j=1,Nsize
             is = i + (iblock-1)*Nsize
             js = j + (iblock-1)*Nsize
             Lmat(Nsize+is,js) =.false.
             Lmat(is,js)       =.false.
             Lmat(is,Nsize+js) =.false.
          enddo
       enddo
    enddo
    do i=1,Nsize
       do j=1,Nsize
          is = i + (Nblock-1)*Nsize
          js = j + (Nblock-1)*Nsize
          Lmat(is,js)=.false.
       enddo
    enddo
    Mcheck = .not.(sum(abs(Amat),mask=Lmat)>0d0)
  end function d_check_tridiag_block
  function c_check_tridiag_block(Nblock,Nsize,Amat) result(Mcheck)
    integer                                         :: Nblock
    integer                                         :: Nsize
    complex(8),dimension(Nblock*Nsize,Nblock*Nsize) :: Amat
    logical,dimension(Nblock*Nsize,Nblock*Nsize)    :: Lmat
    integer                                         :: i,j,iblock,is,js
    logical                                         :: Mcheck
    Lmat=.true.
    do iblock=1,Nblock-1
       do i=1,Nsize
          do j=1,Nsize
             is = i + (iblock-1)*Nsize
             js = j + (iblock-1)*Nsize
             Lmat(Nsize+is,js) =.false.
             Lmat(is,js)       =.false.
             Lmat(is,Nsize+js) =.false.
          enddo
       enddo
    enddo
    do i=1,Nsize
       do j=1,Nsize
          is = i + (Nblock-1)*Nsize
          js = j + (Nblock-1)*Nsize
          Lmat(is,js)=.false.
       enddo
    enddo
    Mcheck = .not.(sum(abs(Amat),mask=Lmat)>0d0)
  end function c_check_tridiag_block



  !+-----------------------------------------------------------------------------+!
  !PURPOSE: Get a the three (sub,main,over) diagonals of a Tridiagonal matrix
  !+-----------------------------------------------------------------------------+!
  subroutine d_get_tridiag(Amat,sub,diag,over)
    real(8),dimension(:,:)                     :: Amat
    real(8),dimension(size(Amat,1))            :: diag
    real(8),dimension(size(Amat,1)-1)          :: sub
    real(8),dimension(size(Amat,1)-1),optional :: over
    real(8),dimension(size(Amat,1)-1)          :: over_
    integer                                    :: i,N
    N=size(Amat,1)
    call assert_shape(Amat,[N,N],"d_get_tridiag","Amat")
    forall(i=1:N-1)
       sub(i)  = Amat(i+1,i)
       diag(i) = Amat(i,i)
       over_(i)= Amat(i,i+1)
    end forall
    diag(N) = Amat(N,N)
    if(present(over))over=over_
  end subroutine d_get_tridiag
  subroutine c_get_tridiag(Amat,sub,diag,over)
    complex(8),dimension(:,:)                     :: Amat
    complex(8),dimension(size(Amat,1))            :: diag
    complex(8),dimension(size(Amat,1)-1)          :: sub
    complex(8),dimension(size(Amat,1)-1),optional :: over
    complex(8),dimension(size(Amat,1)-1)          :: over_
    integer                                       :: i,N
    N=size(Amat,1)
    call assert_shape(Amat,[N,N],"d_get_tridiag","Amat")
    forall(i=1:N-1)
       sub(i)  = Amat(i+1,i)
       diag(i) = Amat(i,i)
       over_(i)= Amat(i,i+1)
    end forall
    diag(N) = Amat(N,N)
    if(present(over))over=over_
  end subroutine c_get_tridiag
  subroutine d_get_tridiag_block(Nblock,Nsize,Amat,sub,diag,over)
    integer                                          :: Nblock
    integer                                          :: Nsize
    real(8),dimension(Nblock*Nsize,Nblock*Nsize)     :: Amat
    real(8),dimension(Nblock-1,Nsize,Nsize)          :: sub
    real(8),dimension(Nblock,Nsize,Nsize)            :: diag
    real(8),dimension(Nblock-1,Nsize,Nsize),optional :: over
    real(8),dimension(Nblock-1,Nsize,Nsize)          :: over_
    integer                                          :: i,j,iblock,is,js
    do iblock=1,Nblock-1
       do i=1,Nsize
          do j=1,Nsize
             is = i + (iblock-1)*Nsize
             js = j + (iblock-1)*Nsize
             Sub(iblock,i,j)   = Amat(Nsize+is,js)
             Diag(iblock,i,j)  = Amat(is,js)
             Over_(iblock,i,j) = Amat(is,Nsize+js)
          enddo
       enddo
    enddo
    do i=1,Nsize
       do j=1,Nsize
          is = i + (Nblock-1)*Nsize
          js = j + (Nblock-1)*Nsize
          Diag(Nblock,i,j) = Amat(is,js)
       enddo
    enddo
    if(present(over))over=over_
  end subroutine d_get_tridiag_block
  subroutine c_get_tridiag_block(Nblock,Nsize,Amat,sub,diag,over)
    integer                                          :: Nblock
    integer                                          :: Nsize
    complex(8),dimension(Nblock*Nsize,Nblock*Nsize)     :: Amat
    complex(8),dimension(Nblock-1,Nsize,Nsize)          :: sub
    complex(8),dimension(Nblock,Nsize,Nsize)            :: diag
    complex(8),dimension(Nblock-1,Nsize,Nsize),optional :: over
    complex(8),dimension(Nblock-1,Nsize,Nsize)          :: over_
    integer                                          :: i,j,iblock,is,js
    do iblock=1,Nblock-1
       do i=1,Nsize
          do j=1,Nsize
             is = i + (iblock-1)*Nsize
             js = j + (iblock-1)*Nsize
             Sub(iblock,i,j)   = Amat(Nsize+is,js)
             Diag(iblock,i,j)  = Amat(is,js)
             Over_(iblock,i,j) = Amat(is,Nsize+js)
          enddo
       enddo
    enddo
    do i=1,Nsize
       do j=1,Nsize
          is = i + (Nblock-1)*Nsize
          js = j + (Nblock-1)*Nsize
          Diag(Nblock,i,j) = Amat(is,js)
       enddo
    enddo
    if(present(over))over=over_
  end subroutine c_get_tridiag_block


  !+-----------------------------------------------------------------------------+!
  !PURPOSE: Build a tridiagonal Matrix Amat from the three (sub,main,over) diagonal.
  ! In this version the over-diagonal is optional
  !+-----------------------------------------------------------------------------+!
  function d_build_tridiag(sub,diag,over) result(Amat)
    real(8),dimension(:)                     :: diag
    real(8),dimension(size(diag)-1)          :: sub
    real(8),dimension(size(diag)-1),optional :: over
    real(8),dimension(size(diag),size(diag)) :: Amat
    real(8),dimension(size(diag)-1)          :: over_
    integer                                  :: i,N
    over_=sub;if(present(over))over_=over
    N=size(diag)
    Amat=0d0
    forall(i=1:N-1)
       Amat(i+1,i) = sub(i)
       Amat(i,i)   = diag(i)
       Amat(i,i+1) = over_(i)
    end forall
    Amat(N,N)=diag(N)
  end function d_build_tridiag
  function c_build_tridiag(sub,diag,over) result(Amat)
    complex(8),dimension(:)                     :: diag
    complex(8),dimension(size(diag)-1)          :: sub
    complex(8),dimension(size(diag)-1),optional :: over
    complex(8),dimension(size(diag),size(diag)) :: Amat
    complex(8),dimension(size(diag)-1)          :: over_
    integer                                     :: i,N
    over_=sub;if(present(over))over_=over
    N=size(diag)
    Amat=dcmplx(0d0,0d0)
    forall(i=1:N-1)
       Amat(i+1,i) = sub(i)
       Amat(i,i)   = diag(i)
       Amat(i,i+1) = over_(i)
    end forall
    Amat(N,N)=diag(N)
  end function c_build_tridiag
  function d_build_tridiag_block(Nblock,Nsize,sub,diag,over) result(Amat)
    integer                                          :: Nblock
    integer                                          :: Nsize
    real(8),dimension(Nblock*Nsize,Nblock*Nsize)     :: Amat
    real(8),dimension(Nblock-1,Nsize,Nsize)          :: sub
    real(8),dimension(Nblock,Nsize,Nsize)            :: diag
    real(8),dimension(Nblock-1,Nsize,Nsize),optional :: over
    real(8),dimension(Nblock-1,Nsize,Nsize)          :: over_
    integer                                          :: i,j,iblock,is,js
    over_=sub;if(present(over))over_=over
    !
    Amat=0d0
    !
    do iblock=1,Nblock-1
       do i=1,Nsize
          do j=1,Nsize
             is = i + (iblock-1)*Nsize
             js = j + (iblock-1)*Nsize
             Amat(Nsize+is,js) = Sub(iblock,i,j)
             Amat(is,js)       = Diag(iblock,i,j)
             Amat(is,Nsize+js) = Over_(iblock,i,j)
          enddo
       enddo
    enddo
    do i=1,Nsize
       do j=1,Nsize
          is = i + (Nblock-1)*Nsize
          js = j + (Nblock-1)*Nsize
          Amat(is,js)       = Diag(Nblock,i,j)
       enddo
    enddo
  end function d_build_tridiag_block
  function c_build_tridiag_block(Nblock,Nsize,sub,diag,over) result(Amat)
    integer                                          :: Nblock
    integer                                          :: Nsize
    complex(8),dimension(Nblock*Nsize,Nblock*Nsize)     :: Amat
    complex(8),dimension(Nblock-1,Nsize,Nsize)          :: sub
    complex(8),dimension(Nblock,Nsize,Nsize)            :: diag
    complex(8),dimension(Nblock-1,Nsize,Nsize),optional :: over
    complex(8),dimension(Nblock-1,Nsize,Nsize)          :: over_
    integer                                          :: i,j,iblock,is,js
    over_=sub;if(present(over))over_=over
    !
    Amat=0d0
    !
    do iblock=1,Nblock-1
       do i=1,Nsize
          do j=1,Nsize
             is = i + (iblock-1)*Nsize
             js = j + (iblock-1)*Nsize
             Amat(Nsize+is,js) = Sub(iblock,i,j)
             Amat(is,js)       = Diag(iblock,i,j)
             Amat(is,Nsize+js) = Over_(iblock,i,j)
          enddo
       enddo
    enddo
    do i=1,Nsize
       do j=1,Nsize
          is = i + (Nblock-1)*Nsize
          js = j + (Nblock-1)*Nsize
          Amat(is,js)       = Diag(Nblock,i,j)
       enddo
    enddo
  end function c_build_tridiag_block









  !+-----------------------------------------------------------------------------+!
  !PURPOSE: return the N diagonal elements of the inverse matrix.
  ! b = sub-diagonal
  ! d = main-diagonal
  ! a = over-diagonal
  !+-----------------------------------------------------------------------------+!
  subroutine d_invert_tridiag_matrix(N,sub_,diag_,over_,Inv)
    integer                         :: i,j,N
    real(8),dimension(N)            :: diag_
    real(8),dimension(N-1)          :: sub_
    real(8),dimension(N-1)          :: over_
    real(8),dimension(N)            :: Inv
    real(8)                         :: Foo,Cleft,Cright
    real(8),dimension(N)            :: dleft
    real(8),dimension(N)            :: dright
    !
    !DOWNWARD:
    dleft(1) = diag_(1)
    if(dleft(1)==0d0)stop "matrix is ill-conditioned: no inverse exists"
    do i=2,N
       foo      = 1d0/dleft(i-1)
       cleft    = sub_(i-1)*foo
       dleft(i) = diag_(i) - cleft*over_(i-1)
       if(dleft(i)==0d0)stop "matrix is ill-conditioned: no inverse exists"
    enddo
    !
    !UPWARD:
    dright(N) = diag_(N)
    if(dright(N)==0d0)stop "matrix is ill-conditioned: no inverse exists"
    do i=N-1,1,-1
       foo       = 1d0/dright(i+1)
       cright    = over_(i)*foo
       dright(i) = diag_(i) - cright*sub_(i)
       if(dright(i)==0d0)stop "matrix is ill-conditioned: no inverse exists"
    enddo
    !
    do i=1,N
       Foo    =  dleft(i) + dright(i) - diag_(i)
       Inv(i) = 1d0/Foo
    end do
  end subroutine d_invert_tridiag_matrix

  subroutine c_invert_tridiag_matrix(N,sub_,diag_,over_,Inv)
    integer                            :: i,j,N
    complex(8),dimension(N)            :: diag_
    complex(8),dimension(N-1)          :: sub_
    complex(8),dimension(N-1)          :: over_
    complex(8),dimension(N)            :: Inv
    complex(8)                         :: Foo,Cleft,Cright
    complex(8),dimension(N)            :: dleft
    complex(8),dimension(N)            :: dright
    !
    !DOWNWARD:
    dleft(1) = diag_(1)
    if(dleft(1)==0d0)stop "matrix is ill-conditioned: no inverse exists"
    do i=2,N
       foo      = 1d0/dleft(i-1)
       cleft    = sub_(i-1)*foo
       dleft(i) = diag_(i) - cleft*over_(i-1) !over_(i-1)/dleft(i-1)*sub_(i)
       if(dleft(i)==0d0)stop "matrix is ill-conditioned: no inverse exists"
    enddo
    !
    !UPWARD:
    dright(N) = diag_(N)
    if(dright(N)==0d0)stop "matrix is ill-conditioned: no inverse exists"
    do i=N-1,1,-1
       foo       = 1d0/dright(i+1)
       cright    = over_(i)*foo
       dright(i) = diag_(i) - cright*sub_(i) !sub_(i+1)/dright(i+1)*over_(i)
       if(dright(i)==0d0)stop "matrix is ill-conditioned: no inverse exists"
    enddo
    !
    do i=1,N
       Foo    =  dleft(i) + dright(i) - diag_(i)
       Inv(i) = 1d0/Foo
    end do
  end subroutine c_invert_tridiag_matrix





  !+-----------------------------------------------------------------------------+!
  !PURPOSE: return the Nb diagonal NxN blocks of the inverse matrix.
  ! b = sub-diagonal
  ! d = main-diagonal
  ! a = over-diagonal = sub-diagonal (symmetric matrix)
  !+-----------------------------------------------------------------------------+!
  subroutine d_invert_tridiag_block_matrix(Nb,N,sub_,diag_,over_,Ainv)
    integer                              :: ib,i,j,Nb,N
    real(8),dimension(Nb,N,N)            :: diag_
    real(8),dimension(Nb-1,N,N)          :: sub_
    real(8),dimension(Nb-1,N,N)          :: over_
    real(8),dimension(Nb,N,N)            :: Ainv
    real(8),dimension(N,N)               :: Foo,Cleft,Cright
    real(8),dimension(Nb,N,N)            :: Dleft
    real(8),dimension(Nb,N,N)            :: Dright
    !
    !DOWNWARD:
    dleft(1,:,:) = diag_(1,:,:)
    do i=2,Nb
       foo  = dleft(i-1,:,:) ; call inv(foo)
       cleft= matmul(sub_(i-1,:,:),foo)
       dleft(i,:,:) = diag_(i,:,:) - matmul(cleft,over_(i-1,:,:))
    enddo
    !
    !BACKWARD:
    dright(Nb,:,:) = diag_(Nb,:,:)
    do i=Nb-1,1,-1
       foo   = dright(i+1,:,:) ; call inv(foo)
       cright= matmul(over_(i,:,:),foo)
       dright(i,:,:) = diag_(i,:,:) - matmul(cright,sub_(i,:,:))
    enddo
    !
    do ib=1,Nb
       Ainv(ib,:,:)    =  dleft(ib,:,:) + dright(ib,:,:) - diag_(ib,:,:)
       call inv(Ainv(ib,:,:))
    end do
  end subroutine d_invert_tridiag_block_matrix

  subroutine c_invert_tridiag_block_matrix(Nb,N,sub_,diag_,over_,Ainv)
    integer                        :: ib,i,j,Nb,N
    complex(8),dimension(Nb,N,N)   :: diag_
    complex(8),dimension(Nb-1,N,N) :: sub_
    complex(8),dimension(Nb-1,N,N) :: over_
    complex(8),dimension(Nb,N,N)   :: Ainv
    complex(8),dimension(N,N)      :: Foo,Cleft,Cright
    complex(8),dimension(Nb,N,N)   :: Dleft
    complex(8),dimension(Nb,N,N)   :: Dright
    !
    !DOWNWARD:
    dleft(1,:,:) = diag_(1,:,:)
    do i=2,Nb
       foo  = dleft(i-1,:,:) ; call inv(foo)
       cleft= matmul(sub_(i-1,:,:),foo)
       dleft(i,:,:) = diag_(i,:,:) - matmul(cleft,over_(i-1,:,:))
    enddo
    !
    !BACKWARD:
    dright(Nb,:,:) = diag_(Nb,:,:)
    do i=Nb-1,1,-1
       foo   = dright(i+1,:,:) ; call inv(foo)
       cright= matmul(over_(i,:,:),foo)
       dright(i,:,:) = diag_(i,:,:) - matmul(cright,sub_(i,:,:))
    enddo
    !
    do ib=1,Nb
       Ainv(ib,:,:)    =  dleft(ib,:,:) + dright(ib,:,:) - diag_(ib,:,:)
       call inv(Ainv(ib,:,:))
    end do
  end subroutine c_invert_tridiag_block_matrix


  !
  !
  !

  subroutine d_invert_tridiag_matrix_mat(Amat)
    real(8),dimension(:,:),intent(inout)         :: Amat
    real(8),dimension(size(Amat,1))              :: diag_
    real(8),dimension(size(Amat,1)-1)            :: sub_
    real(8),dimension(size(Amat,1)-1)            :: over_
    real(8),dimension(size(Amat,1))              :: Inv
    integer                                      :: i,N
    N=size(Amat,1)
    call assert_shape(Amat,[N,N],"d_invert_tridiag_matrix_mat","Amat")
    call get_tridiag(Amat,sub_,diag_,over_)
    call inv_tridiag(N,sub_,diag_,over_,Inv)
    Amat= 0d0
    forall(i=1:N)Amat(i,i)=Inv(i)
  end subroutine d_invert_tridiag_matrix_mat

  subroutine c_invert_tridiag_matrix_mat(Amat)
    complex(8),dimension(:,:),intent(inout)         :: Amat
    complex(8),dimension(size(Amat,1))              :: diag_
    complex(8),dimension(size(Amat,1)-1)            :: sub_
    complex(8),dimension(size(Amat,1)-1)            :: over_
    complex(8),dimension(size(Amat,1))              :: Inv
    integer                                         :: i,N
    N=size(Amat,1)
    call assert_shape(Amat,[N,N],"d_invert_tridiag_matrix_mat","Amat")
    call get_tridiag(Amat,sub_,diag_,over_)
    call inv_tridiag(N,sub_,diag_,over_,Inv)
    Amat= dcmplx(0d0,0d0)
    forall(i=1:N)Amat(i,i)=Inv(i)
  end subroutine c_invert_tridiag_matrix_mat

  subroutine d_invert_tridiag_block_matrix_mat(Nb,N,Amat)
    integer,intent(in)                         :: Nb
    integer,intent(in)                         :: N
    real(8),dimension(Nb*N,Nb*N),intent(inout) :: Amat
    real(8),dimension(Nb-1,N,N)                :: sub_
    real(8),dimension(Nb,N,N)                  :: diag_
    real(8),dimension(Nb-1,N,N)                :: over_
    real(8),dimension(Nb,N,N)                  :: Inv
    integer                                    :: i,j,is,js,iblock
    call get_tridiag(Nb,N,Amat,sub_,diag_,over_)
    call inv_tridiag(Nb,N,sub_,diag_,over_,Inv)
    Amat=0d0
    do iblock=1,Nb
       do i=1,N
          do j=1,N
             is = i + (iblock-1)*N
             js = j + (iblock-1)*N
             Amat(is,js) = Inv(iblock,i,j)
          enddo
       enddo
    enddo
  end subroutine d_invert_tridiag_block_matrix_mat

  subroutine c_invert_tridiag_block_matrix_mat(Nb,N,Amat)
    integer,intent(in)                            :: Nb
    integer,intent(in)                            :: N
    complex(8),dimension(Nb*N,Nb*N),intent(inout) :: Amat
    complex(8),dimension(Nb-1,N,N)                :: sub_
    complex(8),dimension(Nb,N,N)                  :: diag_
    complex(8),dimension(Nb-1,N,N)                :: over_
    complex(8),dimension(Nb,N,N)                  :: Inv
    integer                                       :: i,j,is,js,iblock
    call get_tridiag(Nb,N,Amat,sub_,diag_,over_)
    call inv_tridiag(Nb,N,sub_,diag_,over_,Inv)
    Amat= dcmplx(0d0,0d0)
    do iblock=1,Nb
       do i=1,N
          do j=1,N
             is = i + (iblock-1)*N
             js = j + (iblock-1)*N
             Amat(is,js) = Inv(iblock,i,j)
          enddo
       enddo
    enddo
  end subroutine c_invert_tridiag_block_matrix_mat








  subroutine d_jacobi(a,d,v,nrot)
    real(8),dimension(:,:),intent(inout)               :: a
    real(8),dimension(size(a,1)),intent(out)           :: d
    real(8),dimension(size(a,1),size(a,2)),intent(out) :: v
    integer,intent(out)                                :: nrot
    integer                                            :: i,ip,iq,n
    real(8)                                            :: c,g,h,s,sm,t,tau,theta,tresh
    real(8),dimension(size(d))                         :: b,z
    !
    n=size(a,1);if(size(a,2)/=n)stop "Error in Jacobi: size(a)!=n**2 - a not a square matrix"
    !
    v = eye(n)
    forall(i=1:n)b(i)=a(i,i)
    !
    d=b
    z=0d0
    nrot=0
    !
    do i=1,50
       sm=sum(abs(a),mask=upper_triangle(n,n))
       if (sm == 0.0) return
       tresh=merge(0.2d0*sm/n**2, 0d0, i<4 )
       do ip=1,n-1
          do iq=ip+1,n
             g=100d0*abs(a(ip,iq))
             if ((i > 4) .and. (abs(d(ip))+g == abs(d(ip))) &
                  .and. (abs(d(iq))+g == abs(d(iq)))) then
                a(ip,iq)=0.0
             else if (abs(a(ip,iq)) > tresh) then
                h=d(iq)-d(ip)
                if (abs(h)+g == abs(h)) then
                   t=a(ip,iq)/h
                else
                   theta=0.5d0*h/a(ip,iq)
                   t=1.0d0/(abs(theta)+sqrt(1.0d0+theta**2))
                   if (theta < 0.0) t=-t
                end if
                c=1d0/sqrt(1+t**2)
                s=t*c
                tau=s/(1d0+c)
                h=t*a(ip,iq)
                z(ip)=z(ip)-h
                z(iq)=z(iq)+h
                d(ip)=d(ip)-h
                d(iq)=d(iq)+h
                a(ip,iq)=0d0
                call jrotate(a(1:ip-1,ip),a(1:ip-1,iq))
                call jrotate(a(ip,ip+1:iq-1),a(ip+1:iq-1,iq))
                call jrotate(a(ip,iq+1:n),a(iq,iq+1:n))
                call jrotate(v(:,ip),v(:,iq))
                nrot=nrot+1
             end if
          end do
       end do
       b(:)=b(:)+z(:)
       d(:)=b(:)
       z(:)=0d0
    end do
    stop "Too many iteration in Jacobi"
  contains
    subroutine jrotate(a1,a2)
      real(8),dimension(:), intent(inout) :: a1,a2
      real(8),dimension(size(a1))         :: wk1
      wk1(:)=a1(:)
      a1(:)=a1(:)-s*(a2(:)+a1(:)*tau)
      a2(:)=a2(:)+s*(wk1(:)-a2(:)*tau)
    end subroutine jrotate
  end subroutine d_jacobi

  subroutine c_jacobi(A,D,U,sweep)
    implicit none
    complex(8),dimension(:,:),intent(inout)               :: a
    real(8),dimension(size(a,1)),intent(out)              :: d
    complex(8),dimension(size(a,1),size(a,2)),intent(out) :: U
    integer                                               :: n
    integer                                               :: p, q, j
    real(8)                                               :: red, off, thresh
    real(8)                                               :: t, delta, invc, s
    complex(8)                                            :: x, y, Apq
    real(8)                                               :: ev(2,size(a,1))
    integer                                               :: sweep
    real(8)                                               :: SYM_EPS=tiny(1d0)
    !
    n=size(a,1);if(size(a,2)/=n)stop "Error in Jacobi: size(a)!=n**2 - a not a square matrix"
    !
    do p = 1, n
       ev(1,p) = 0d0
       ev(2,p) = dble(A(p,p))
       d(p) = ev(2,p)
    enddo
    U = eye(n)
    red = 0.04d0/n**4
    do sweep = 1, 50
       off = 0
       do q = 2, n
          do p = 1, q - 1
             off = off + abs(A(p,q))**2
          enddo
       enddo
       if( .not. off .gt. SYM_EPS ) return !goto 10
       thresh = 0
       if( sweep .lt. 4 ) thresh = off*red
       do q = 2, n
          do p = 1, q - 1
             Apq = A(p,q)
             off = abs(Apq)**2
             if( sweep .gt. 4 .and. off .lt. SYM_EPS*(ev(2,p)**2 + ev(2,q)**2) ) then
                A(p,q) = 0
             else if( off .gt. thresh ) then
                t = .5D0*(ev(2,p) - ev(2,q))
                t = 1/(t + sign(sqrt(t**2 + off), t))
                delta = t*off
                ev(1,p) = ev(1,p) + delta
                ev(2,p) = d(p) + ev(1,p)
                ev(1,q) = ev(1,q) - delta
                ev(2,q) = d(q) + ev(1,q)
                invc = sqrt(delta*t + 1)
                s = t/invc
                t = delta/(invc + 1)
                do j = 1, p - 1
                   x = A(j,p)
                   y = A(j,q)
                   A(j,p) = x + s*(conjg(Apq)*y - t*x)
                   A(j,q) = y - s*(Apq*x + t*y)
                enddo
                do j = p + 1, q - 1
                   x = A(p,j)
                   y = A(j,q)
                   A(p,j) = x + s*(Apq*conjg(y) - t*x)
                   A(j,q) = y - s*(Apq*conjg(x) + t*y)
                enddo
                do j = q + 1, n
                   x = A(p,j)
                   y = A(q,j)
                   A(p,j) = x + s*(Apq*y - t*x)
                   A(q,j) = y - s*(conjg(Apq)*x + t*y)
                enddo
                A(p,q) = 0d0
                do j = 1, n
                   x = U(p,j)
                   y = U(q,j)
                   ! #if UCOLS
                   ! U(p,j) = x + s*(conjg(Apq)*y - t*x)
                   ! U(q,j) = y - s*(Apq*x + t*y)
                   ! #else
                   U(p,j) = x + s*(Apq*y - t*x)
                   U(q,j) = y - s*(conjg(Apq)*x + t*y)
                   ! #endif
                enddo
             endif
          enddo
       enddo
       do p = 1, n
          ev(1,p) = 0
          d(p) = ev(2,p)
       enddo
    enddo
    stop "Too many iteration in Jacobi"
    ! 10  continue
    !     ! if( sort .eq. 0 ) return
    !     do p = 1, n - 1
    !        j = p
    !        t = d(p)
    !        do q = p + 1, n
    !           if( 1*(t - d(q)) .gt. 0 ) then
    !              j = q
    !              t = d(q)
    !           endif
    !        enddo
    !        if( j .ne. p ) then
    !           d(j) = d(p)
    !           d(p) = t
    !           do q = 1, n
    !              x = U(p,q)
    !              U(p,q) = U(j,q)
    !              U(j,q) = x
    !           enddo
    !        endif
    !     enddo
  end subroutine C_jacobi

  function upper_triangle(j,k,extra)
    integer,intent(in)           :: j,k
    integer,optional, intent(in) :: extra
    logical,dimension(j,k)       :: upper_triangle
    integer                      :: n
    n=0
    if (present(extra)) n=extra
    upper_triangle=(outerdiff(arth(1,1,j),arth(1,1,k)) < n)
  end function upper_triangle

  function outerdiff(a,b)
    integer,dimension(:), intent(in)   :: a,b
    integer,dimension(size(a),size(b)) :: outerdiff
    outerdiff = spread(a,dim=2,ncopies=size(b)) - &
         spread(b,dim=1,ncopies=size(a))
  end function outerdiff

  function arth(first,increment,n) result(arth_i)
    integer,parameter    :: npar_arth=16,npar2_arth=8
    integer,intent(in)   :: first,increment,n
    integer,dimension(n) :: arth_i
    integer              :: k,k2,temp
    if (n > 0) arth_i(1)=first
    if (n <= npar_arth) then
       do k=2,n
          arth_i(k)=arth_i(k-1)+increment
       end do
    else
       do k=2,npar2_arth
          arth_i(k)=arth_i(k-1)+increment
       end do
       temp=increment*npar2_arth
       k=npar2_arth
       do
          if (k >= n) exit
          k2=k+k
          arth_i(k+1:min(k2,n))=temp+arth_i(1:min(k,n-k))
          temp=temp+temp
          k=k2
       end do
    end if
  end function arth








  !******************************************************************
  !                AUXILIARY ROUTINES
  !******************************************************************

  !-------------------------------------------------------------------------------------------
  !PURPOSE:   ! construct real matrix from diagonal elements
  !-------------------------------------------------------------------------------------------
  function ddiag(x) result(A)
    real(8), intent(in)  :: x(:)
    real(8), allocatable :: A(:,:)
    integer              :: i, n
    n = size(x)
    allocate(A(n,n))
    A(:,:) = 0d0
    forall(i=1:n) A(i,i) = x(i)
  end function ddiag

  function zdiag(x) result(A)
    complex(8), intent(in)  :: x(:)
    complex(8), allocatable :: A(:,:)
    integer                 :: i, n
    n = size(x)
    allocate(A(n,n))
    A(:,:) = zero
    forall(i=1:n) A(i,i) = x(i)
  end function zdiag



  !-------------------------------------------------------------------------------------------
  !PURPOSE:  return trace along the main diagonal
  !-------------------------------------------------------------------------------------------
  function dtrace(A) result(t)
    real(8), intent(in) :: A(:,:)
    real(8)             :: t
    integer             :: i
    t = 0d0
    do i = 1,minval(shape(A))
       t = t + A(i,i)
    end do
  end function dtrace

  function ztrace(A) result(t)
    complex(8), intent(in) :: A(:,:)
    complex(8)             :: t
    integer                :: i
    t = zero
    do i = 1,minval(shape(A))
       t = t + A(i,i)
    end do
  end function ztrace






  !-------------------------------------------------------------------------------------------
  !PURPOSE:  Returns the identity matrix of size n x n and type real.
  !-------------------------------------------------------------------------------------------
  function deye(n) result(A)
    integer, intent(in) :: n
    real(8)             :: A(n, n)
    integer             :: i
    A = 0d0
    do i = 1, n
       A(i,i) = 1d0
    end do
  end function deye

  function zeye(n) result(A)
    integer, intent(in) :: n
    complex(8)          :: A(n, n)
    integer             :: i
    A = zero
    do i = 1, n
       A(i,i) = one
    end do
  end function zeye

  function deye_tridiag(Nblock,N) result(eye_block)
    integer                       :: Nblock
    integer                       :: N
    real(8),dimension(Nblock,N,N) :: eye_block
    integer                       :: iblock
    do iblock=1,Nblock
       eye_block(iblock,:,:) = eye(N)
    enddo
  end function deye_tridiag

  function zeye_tridiag(Nblock,N) result(eye_block)
    integer                          :: Nblock
    integer                          :: N
    complex(8),dimension(Nblock,N,N) :: eye_block
    integer                          :: iblock
    do iblock=1,Nblock
       eye_block(iblock,:,:) = zeye(N)
    enddo
  end function zeye_tridiag


  !---------------------------------------------------------------------
  !PURPOSE: Function to compute the tensor product (M1_kp_M2) of 
  ! two complex matrices M1 and M2. nr1(nr2) and nc1(nc2) are 
  ! the number of rows and columns of the Matrix M1 and M2
  !---------------------------------------------------------------------
  function i_kronecker_product(A,B) result(AxB)
    integer,intent(in) :: A(:,:), B(:,:)
    integer            :: i,j
    integer            :: rowA,colA
    integer            :: rowB,colB
    integer            :: AxB(size(A,1)*size(B,1),size(A,2)*size(B,2))
    AxB = 0
    rowA=size(A,1) ; colA=size(A,2)
    rowB=size(B,1) ; colB=size(B,2)
    forall(i=1:rowA,j=1:colA)
       AxB(1+rowB*(i-1):rowB*i,1+colB*(j-1):colB*j)  =  A(i,j)*B(:,:)
    end forall
  end function i_kronecker_product
  !
  function d_kronecker_product(A,B) result(AxB)
    real(8),intent(in) :: A(:,:), B(:,:)
    integer            :: i,j
    integer            :: rowA,colA
    integer            :: rowB,colB
    real(8)            :: AxB(size(A,1)*size(B,1),size(A,2)*size(B,2))
    AxB = 0
    rowA=size(A,1) ; colA=size(A,2)
    rowB=size(B,1) ; colB=size(B,2)
    forall(i=1:rowA,j=1:colA)
       AxB(1+rowB*(i-1):rowB*i,1+colB*(j-1):colB*j)  =  A(i,j)*B(:,:)
    end forall
  end function d_kronecker_product
  !
  function c_kronecker_product(A,B) result(AxB)
    complex(8),intent(in) :: A(:,:), B(:,:)
    integer               :: i,j
    integer               :: rowA,colA
    integer               :: rowB,colB
    complex(8)            :: AxB(size(A,1)*size(B,1),size(A,2)*size(B,2))
    AxB = 0
    rowA=size(A,1) ; colA=size(A,2)
    rowB=size(B,1) ; colB=size(B,2)
    forall(i=1:rowA,j=1:colA)
       AxB(1+rowB*(i-1):rowB*i,1+colB*(j-1):colB*j)  =  A(i,j)*B(:,:)
    end forall
  end function c_kronecker_product




  !+-----------------------------------------------------------------------------+!
  !PURPOSE: Form a matrix A(:,:) from the outerproduct of two 1d arrays:
  ! A(i,j) = a_i*b_j
  !+-----------------------------------------------------------------------------+!
  function outerprod_d(a,b) result(outerprod)
    real(8), dimension(:), intent(in)   :: a,b
    real(8), dimension(size(a),size(b)) :: outerprod
    outerprod = spread(a,dim=2,ncopies=size(b)) * &
         spread(b,dim=1,ncopies=size(a))
  end function outerprod_d
  function outerprod_c(a,b) result(outerprod)
    complex(8), dimension(:), intent(in)   :: a,b
    complex(8), dimension(size(a),size(b)) :: outerprod
    outerprod = spread(a,dim=2,ncopies=size(b)) * &
         spread(b,dim=1,ncopies=size(a))
  end function outerprod_c



  !-------------------------------------------------------------------------------------------
  !PURPOSE: Asser the correct shape of a matrix
  !-------------------------------------------------------------------------------------------
  subroutine dassert_shape(A, shap, routine, matname)
    real(8),intent(in)  :: A(:,:)
    integer,intent(in)  :: shap(:)
    character(len=*)    :: routine, matname
    if(any(shape(A) /= shap)) then
       print*, "In routine " // routine // " matrix " // matname // " has illegal shape ", shape(A)
       print*, "Shape should be ", shap
       stop "Aborting due to illegal matrix operation"
    end if
  end subroutine dassert_shape
  subroutine zassert_shape(A, shap, routine, matname)
    complex(8),intent(in) :: A(:,:)
    integer,intent(in)    :: shap(:)
    character(len=*)      :: routine, matname
    if(any(shape(A) /= shap)) then
       print*, "In routine " // routine // " matrix " // matname // " has illegal shape ", shape(A)
       print*, "Shape should be ", shap
       stop "Aborting due to illegal matrix operation"
    end if
  end subroutine zassert_shape




  !+-----------------------------------------------------------------+
  !program  : swap
  !+-----------------------------------------------------------------+
  subroutine swap_i(a,b)
    integer, intent(inout) :: a,b
    integer :: dum
    dum=a
    a=b
    b=dum
  end subroutine swap_i
  !-----------------------------
  subroutine swap_r(a,b)
    real(8), intent(inout) :: a,b
    real(8) :: dum
    dum=a
    a=b
    b=dum
  end subroutine swap_r
  !-----------------------------
  subroutine swap_rv(a,b)
    real(8), dimension(:), intent(inout) :: a,b
    real(8), dimension(size(a)) :: dum
    dum=a
    a=b
    b=dum
  end subroutine swap_rv
  !-----------------------------
  subroutine swap_z(a,b)
    complex(8), intent(inout) :: a,b
    complex(8) :: dum
    dum=a
    a=b
    b=dum
  end subroutine swap_z
  !-----------------------------
  subroutine swap_zv(a,b)
    complex(8), dimension(:), intent(inout) :: a,b
    complex(8), dimension(size(a)) :: dum
    dum=a
    a=b
    b=dum
  end subroutine swap_zv
  !-----------------------------
  subroutine swap_zm(a,b)
    complex(8), dimension(:,:), intent(inout) :: a,b
    complex(8), dimension(size(a,1),size(a,2)) :: dum
    dum=a
    a=b
    b=dum
  end subroutine swap_zm







end module SF_LINALG


