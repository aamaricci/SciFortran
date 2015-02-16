module SF_LINALG
  USE SF_CONSTANTS, only: one,xi,zero
  implicit none
  private


  !Eigenvalue/-vector problem for general matrices:
  interface eig
     module procedure deig
     module procedure zeig
  end interface eig
  public :: eig

  ! eigenvalue/-vector problem for real symmetric/complex hermitian matrices:
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
  public :: eigh
  public :: matrix_diagonalize

  ! eigenvalues for general matrices:
  interface eigvals
     module procedure deigvals
     module procedure zeigvals
  end interface eigvals
  public :: eigvals

  ! matrix inversion for real/complex matrices:
  interface inv
     module procedure dinv
     module procedure zinv
  end interface inv
  interface matrix_inverse
     module procedure dinv
     module procedure zinv
  end interface matrix_inverse
  public :: inv
  public :: matrix_inverse

  ! matrix inversion for real/complex symmetric matrices:
  interface inv_sym
     module procedure dinv_sym
     module procedure zinv_sym
  end interface inv_sym
  interface matrix_inverse_sym
     module procedure dinv_sym
     module procedure zinv_sym
  end interface matrix_inverse_sym
  public :: inv_sym
  public :: matrix_inverse_sym


  ! matrix inversion for complex hermitian matrices:
  interface inv_her
     module procedure zinv_her
  end interface inv_her
  interface matrix_inverse_her
     module procedure zinv_her
  end interface matrix_inverse_her
  public :: inv_her
  public :: matrix_inverse_her

  ! matrix inversion for real/complex triangular matrices:
  interface inv_triang
     module procedure dinv_triang
     module procedure zinv_triang
  end interface inv_triang
  interface matrix_inverse_triang
     module procedure dinv_triang
     module procedure zinv_triang
  end interface matrix_inverse_triang
  public :: inv_triang
  public :: matrix_inverse_triang

  ! matrix inversion for real/complex  matrices using Gauss-Jordan elimination:
  interface inv_gj
     module procedure dinv_gj
     module procedure zinv_gj
  end interface inv_gj
  interface matrix_inverse_gj
     module procedure dinv_gj
     module procedure zinv_gj
  end interface matrix_inverse_gj
  public :: inv_gj
  public :: matrix_inverse_gj


  ! solution to linear systems of equation with real/complex coefficients:
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
  public :: solve
  public :: solve_linear_system

  ! determinants of real/complex square matrices:
  interface det
     module procedure ddet
     module procedure zdet
  end interface det
  public :: det

  !Returns the real/complex identity matrix of size n x n .
  interface eye
     module procedure deye
  end interface eye
  public :: eye
  public :: deye
  public :: zeye

  !least square solutions the real/complex systems of equations of possibly non-square shape:
  interface lstsq
     module procedure dlstsq
     module procedure zlstsq
  end interface lstsq
  public :: lstsq

  !construction of square matrices from the diagonal elements:
  interface diag
     module procedure ddiag
     module procedure zdiag
  end interface diag
  public :: diag

  !trace of real/complex matrices:
  interface trace
     module procedure dtrace
     module procedure ztrace
  end interface trace
  public :: trace

  !singular values of real/complex matrices:
  interface svdvals
     module procedure dsvdvals
     module procedure zsvdvals
  end interface svdvals
  public :: svdvals

  !singular value decomposition of real/complex matrices:
  interface svd
     module procedure dsvd
     module procedure zsvd
  end interface svd
  public :: svd

  !Kroenecker product of matrices
  interface kron
     module procedure i_kronecker_product,d_kronecker_product,c_kronecker_product
  end interface kron
  interface kronecker_product
     module procedure i_kronecker_product,d_kronecker_product,c_kronecker_product
  end interface kronecker_product
  interface kroenecker_product
     module procedure i_kronecker_product,d_kronecker_product,c_kronecker_product
  end interface kroenecker_product
  public :: kron
  public :: kronecker_product
  public :: kroenecker_product


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
    complex(8), allocatable :: lam(:)  ! eigenvalues: A c = lam c
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
    complex(8), allocatable :: lam(:)  ! eigenvalues: A c = lam c
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



  !---------------------------------------------------------------------
  !PURPOSE: Function to compute the tensor product (M1_kp_M2) of 
  ! two complex matrices M1 and M2. nr1(nr2) and nc1(nc2) are 
  ! the number of rows and columns of the Matrix M1 and M2
  !---------------------------------------------------------------------
  ! function i_kronecker_product(M1, nr1, nc1, M2, nr2, nc2) result(M1_kp_M2)
  !   integer               :: i, j
  !   integer,intent(in)    :: nr1,nc1,nr2,nc2
  !   integer,intent(in)    :: M1(nr1,nc1), M2(nr2,nc2)
  !   integer               :: M1_kp_M2(nr1*nr2,nc1*nc2)
  !   M1_kp_M2 = zero
  !   forall(i =1:nr1,j=1:nc1)
  !      M1_kp_M2(nr2*(i-1)+1 : nr2*i , nc2*(j-1)+1 : nc2*j)  =  M1(i,j)*M2
  !   end forall
  ! end function i_kronecker_product
  ! !
  ! function d_kronecker_product(M1, nr1, nc1, M2, nr2, nc2) result(M1_kp_M2)
  !   integer               :: i, j
  !   integer,intent(in)    :: nr1,nc1,nr2,nc2
  !   real(8),intent(in)    :: M1(nr1,nc1), M2(nr2,nc2)
  !   real(8)               :: M1_kp_M2(nr1*nr2,nc1*nc2)
  !   M1_kp_M2 = zero
  !   forall(i =1:nr1,j=1:nc1)
  !      M1_kp_M2(nr2*(i-1)+1 : nr2*i , nc2*(j-1)+1 : nc2*j)  =  M1(i,j)*M2
  !   end forall
  ! end function d_kronecker_product
  ! !
  ! function c_kronecker_product(M1, nr1, nc1, M2, nr2, nc2) result(M1_kp_M2)
  !   integer               :: i, j
  !   integer,intent(in)    :: nr1,nc1,nr2,nc2
  !   complex(8),intent(in) :: M1(nr1,nc1), M2(nr2,nc2)
  !   complex(8)            :: M1_kp_M2(nr1*nr2,nc1*nc2)
  !   M1_kp_M2 = zero
  !   forall(i =1:nr1,j=1:nc1)
  !      M1_kp_M2(nr2*(i-1)+1 : nr2*i , nc2*(j-1)+1 : nc2*j)  =  M1(i,j)*M2
  !   end forall
  ! end function c_kronecker_product
  function i_kronecker_product(M1,M2) result(M1_kp_M2)
    integer               :: i, j
    integer               :: nr1,nc1,nr2,nc2
    integer,intent(in)    :: M1(:,:), M2(:,:)
    integer               :: M1_kp_M2(size(M1,1)*size(M2,1),size(M1,2)*size(M2,2))
    M1_kp_M2 = zero
    nr1=size(M1,1)
    nc1=size(M1,2)
    nr2=size(M2,1)
    nc2=size(m2,2)
    forall(i =1:nr1,j=1:nc1)
       M1_kp_M2(nr2*(i-1)+1 : nr2*i , nc2*(j-1)+1 : nc2*j)  =  M1(i,j)*M2
    end forall
  end function i_kronecker_product
  !
  function d_kronecker_product(M1,M2) result(M1_kp_M2)
    integer               :: i, j
    integer               :: nr1,nc1,nr2,nc2
    real(8),intent(in)    :: M1(:,:), M2(:,:)
    real(8)               :: M1_kp_M2(size(M1,1)*size(M2,1),size(M1,2)*size(M2,2))
    M1_kp_M2 = zero
    nr1=size(M1,1)
    nc1=size(M1,2)
    nr2=size(M2,1)
    nc2=size(m2,2)
    forall(i =1:nr1,j=1:nc1)
       M1_kp_M2(nr2*(i-1)+1 : nr2*i , nc2*(j-1)+1 : nc2*j)  =  M1(i,j)*M2
    end forall
  end function d_kronecker_product
  !
  function c_kronecker_product(M1,M2) result(M1_kp_M2)
    integer               :: i, j
    integer               :: nr1,nc1,nr2,nc2
    complex(8),intent(in) :: M1(:,:), M2(:,:)
    complex(8)            :: M1_kp_M2(size(M1,1)*size(M2,1),size(M1,2)*size(M2,2))
    M1_kp_M2 = zero
    nr1=size(M1,1)
    nc1=size(M1,2)
    nr2=size(M2,1)
    nc2=size(m2,2)
    forall(i =1:nr1,j=1:nc1)
       M1_kp_M2(nr2*(i-1)+1 : nr2*i , nc2*(j-1)+1 : nc2*j)  =  M1(i,j)*M2
    end forall
  end function c_kronecker_product








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


