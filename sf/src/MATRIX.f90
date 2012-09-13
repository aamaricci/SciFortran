!###############################################################
! PROGRAM  : MATRIX
! TYPE     : Module
! PURPOSE  : A collection of MATRIX related routines
!###############################################################
module MATRIX
  USE COMMON_VARS
  implicit none
  private
  include "mkl_lapack.fi"

  !----------------------------------------------------------------------------------

  interface matrix_inverse
     module procedure D_mat_invert, Z_mat_invert
  end interface matrix_inverse

  interface matrix_inverse_sym
     module procedure d_mat_invertsym, z_mat_invertsym
  end interface matrix_inverse_sym

  interface matrix_inverse_her
     module procedure z_mat_inverther
  end interface matrix_inverse_her

  interface matrix_inverse_triang
     module procedure d_mat_inverttriang, z_mat_inverttriang
  end interface matrix_inverse_triang

  interface matrix_inverse_gj
     module procedure d_mat_invertgj, z_mat_invertgj
  end interface matrix_inverse_gj

  !----------------------------------------------------------------------------------

  interface m_invert_gj
     module procedure d_invert_matrix_gj,c_invert_matrix_gj
  end interface m_invert_gj

  interface m_invert
     module procedure d_invert_matrix,c_invert_matrix
  end interface m_invert

  !----------------------------------------------------------------------------------

  interface matrix_diagonalize
     module procedure d_mat_diagonalization,z_mat_diagonalization
  end interface matrix_diagonalize

  !----------------------------------------------------------------------------------

  interface solve_linear_system
     module procedure d_mat_solve_linear_system_1rhs,z_mat_solve_linear_system_1rhs,&
          d_mat_solve_linear_system_Mrhs,z_mat_solve_linear_system_Mrhs
  end interface solve_linear_system


  !from nr90
  interface swap
     module procedure swap_i,swap_r,swap_rv,swap_z,swap_zv,swap_zm
  end interface swap

  public :: matrix_diagonalize
  !
  public :: solve_linear_system
  !
  public :: matrix_inverse
  public :: matrix_inverse_sym
  public :: matrix_inverse_her
  public :: matrix_inverse_triang
  public :: matrix_inverse_gj
  !
  public :: m_invert
  public :: m_invert_gj
  ! 


contains



  !+-----------------------------------------------------------------+
  !PROGRAM  : MAT_DIAGONALIZATION
  !PURPOSE  : function interface to matrix inversion subroutine
  !+-----------------------------------------------------------------+
  subroutine d_mat_diagonalization(M,E,jobz,uplo)
    real(8),dimension(:,:),intent(inout) :: M
    real(8),dimension(:),intent(inout)   :: E
    character(len=1),optional            :: jobz,uplo
    character(len=1)                     :: jobz_,uplo_
    integer                              :: i,j,n,lda,info,lwork
    real(8),dimension(1)                 :: lwork_guess
    real(8),dimension(:),allocatable     :: work
    jobz_='N';if(present(jobz))jobz_=jobz
    uplo_='U';if(present(uplo))uplo_=uplo
    lda = max(1,size(M,1))
    n   = size(M,2)
    Call dsyev(jobz_,uplo_,n,M,lda,E,lwork_guess,-1,info)
    if(info /= 0)call error("Error MATRIX/d_mat_diagonalization: 1st call dsyev")
    lwork=lwork_guess(1) ; allocate(work(lwork))
    call dsyev(jobz_,uplo_,n,M,lda,E,work,lwork,info)
    if(info /= 0)call error("Error MATRIX/d_mat_diagonalization: 2ns call dsyev")
    deallocate(work)
  end subroutine d_mat_diagonalization
  !-----------------------------

  !-----------------------------
  subroutine z_mat_diagonalization(M,E,jobz,uplo)
    complex(8),dimension(:,:),intent(inout) :: M
    real(8),dimension(:),intent(inout)      :: E
    character(len=1),optional               :: jobz,uplo
    character(len=1)                        :: jobz_,uplo_
    integer                                 :: i,j,n,lda,info,lwork
    complex(8),dimension(1)                 :: lwork_guess
    complex(8),dimension(:),allocatable     :: work
    real(8),dimension(:),allocatable        :: rwork
    jobz_='N';if(present(jobz))jobz_=jobz
    uplo_='U';if(present(uplo))uplo_=uplo
    lda = max(1,size(M,1))
    n   = size(M,2)
    allocate(rwork(max(1,3*N-2)))
    call zheev(jobz_,uplo_,n,M,lda,E,lwork_guess,-1,rwork,info)
    if(info /= 0)call error("Error MATRIX/d_mat_diagonalization: 1st call zsyev")
    lwork=lwork_guess(1) ; allocate(work(lwork))
    call zheev(jobz_,uplo_,n,M,lda,E,work,lwork,rwork,info)
    if(info /= 0)call error("Error MATRIX/d_mat_diagonalization: 2Nd call zsyev")
    deallocate(work,rwork)
  end subroutine z_mat_diagonalization





  !******************************************************************
  !******************************************************************
  !******************************************************************




  !+-----------------------------------------------------------------+
  !PROGRAM  : SOLVE_LINEAR_SYSTEM
  !PURPOSE  : Interface to lapack solution of linear system: A*x=b
  !+-----------------------------------------------------------------+
  subroutine d_mat_solve_linear_system_1rhs(A,b,trans)
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
    if(info/=0)call error("Error MATRIX/d_mat_solve_linear_system: dgetrf")    
    allocate(b_(ldb,nrhs))
    call dgetrs(trans_,n,nrhs,A,lda,ipvt,b_,ldb,info)
    if(info/=0)call error("Error MATRIX/d_mat_solve_linear_system: dgetrs")
    b=b_(:,1)
    deallocate(ipvt,b_)
  end subroutine d_mat_solve_linear_system_1rhs
  !-----------------------------

  !-----------------------------
  subroutine z_mat_solve_linear_system_1rhs(A,b,trans)
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
    if(info/=0)call error("Error MATRIX/d_mat_solve_linear_system: dgetrf")    
    allocate(b_(ldb,nrhs))
    call zgetrs(trans_,n,nrhs,A,lda,ipvt,b_,ldb,info)
    if(info/=0)call error("Error MATRIX/d_mat_solve_linear_system: dgetrs")
    b=b_(:,1)
    deallocate(ipvt,b_)
  end subroutine z_mat_solve_linear_system_1rhs
  !-----------------------------

  !-----------------------------
  subroutine d_mat_solve_linear_system_Mrhs(A,b,trans)
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
    if(info/=0)call error("Error MATRIX/d_mat_solve_linear_system: dgetrf")    
    call dgetrs(trans_,n,nrhs,A,lda,ipvt,b,ldb,info)
    if(info/=0)call error("Error MATRIX/d_mat_solve_linear_system: dgetrs")
    deallocate(ipvt)
  end subroutine d_mat_solve_linear_system_Mrhs
  !-----------------------------

  !-----------------------------
  subroutine z_mat_solve_linear_system_Mrhs(A,b,trans)
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
    if(info/=0)call error("Error MATRIX/d_mat_solve_linear_system: dgetrf")    
    lda=n ; ldb=n ; nrhs=size(b,2)
    call zgetrs(trans_,n,nrhs,A,lda,ipvt,b,ldb,info)
    if(info/=0)call error("Error MATRIX/d_mat_solve_linear_system: dgetrs")
    deallocate(ipvt)
  end subroutine z_mat_solve_linear_system_Mrhs






  !******************************************************************
  !******************************************************************
  !******************************************************************






  !+-----------------------------------------------------------------+
  !PROGRAM  : invert_matrix_xx
  !PURPOSE  : function interface to matrix inversion subroutine
  !+-----------------------------------------------------------------+
  function d_invert_matrix_gj(M,L) result(M1)
    integer                   :: L
    real(8),dimension(L,L)    :: M,M1
    M1=M
    call matrix_inverse_GJ(M1)
  end function d_invert_matrix_gj
  function c_invert_matrix_gj(M,L) result(M1)
    integer                   :: L
    complex(8),dimension(L,L) :: M,M1
    M1=M
    call matrix_inverse_GJ(M1)
  end function c_invert_matrix_gj
  !
  function d_invert_matrix(M,L) result(M1)
    integer                   :: L
    real(8),dimension(L,L)    :: M,M1
    M1=M
    call matrix_inverse(M1)
  end function d_invert_matrix
  function c_invert_matrix(M,L) result(M1)
    integer                   :: L
    complex(8),dimension(L,L) :: M,M1
    M1=M
    call matrix_inverse(M1)
  end function c_invert_matrix




  !******************************************************************
  !******************************************************************
  !******************************************************************





  !+-----------------------------------------------------------------+
  !PROGRAM  : MAT_INVERSION(M)
  !TYPE     : Subroutine
  !PURPOSE  : Invert a general m*n matrix using LAPACK library 
  !COMMENT  : M is destroyed and replaces by its inverse M^-1
  !+-----------------------------------------------------------------+
  subroutine D_mat_invert(A)
    real(8),dimension(:,:)           :: A
    integer                          :: m,n,lda,info,lwork
    integer,dimension(:),allocatable :: ipvt
    real(8),dimension(:),allocatable :: work
    real(8),dimension(1)             :: lwork_guess
    lda = max(1,size(A,1))
    m = size(A,1)
    n = size(A,2)
    allocate(ipvt(min(m,n)))
    call dgetrf(m,n,A,lda,ipvt,info)
    if(info/=0)call error("Error MATRIX/D_mat_invert: dgetrf")
    call dgetri(n,A,lda,ipvt,lwork_guess,-1,info)
    if(info/=0)call error("Error MATRIX/D_mat_invert: 1st call dgetri")
    lwork=lwork_guess(1) ; allocate(work(lwork))
    call dgetri(n,A,lda,ipvt,work,lwork,info)
    if(info/=0)call error("Error MATRIX/D_mat_invert: 2nd call dgetri")
    deallocate(ipvt,work)
  end subroutine D_mat_invert
  !-----------------------------

  !-----------------------------
  subroutine Z_mat_invert(A)
    complex(8),dimension(:,:)            :: A
    integer                              :: m,n,lda,info,lwork
    integer,dimension(:),allocatable     :: ipvt
    complex(8),dimension(:),allocatable  :: work
    complex(8),dimension(1)              :: lwork_guess
    lda = max(1,size(A,1))
    m   = size(A,1)
    n   = size(A,2)
    allocate(ipvt(min(m,n)))
    call zgetrf(m,n,A,lda,ipvt,info)
    if(info/=0)call error("Error MATRIX/Z_mat_invert: zgetrf")
    call zgetri(n,A,lda,ipvt,lwork_guess,-1,info)
    if(info/=0)call error("Error MATRIX/Z_mat_invert: 1st call dgetri")
    lwork=lwork_guess(1) ; allocate(work(lwork))
    call zgetri(n,A,lda,ipvt,work,lwork,info)
    if(info/=0)call error("Error MATRIX/Z_mat_invert: zgetri")
    deallocate(ipvt,work)
  end subroutine Z_mat_invert




  !******************************************************************
  !******************************************************************
  !******************************************************************











  !+-----------------------------------------------------------------+
  !PROGRAM  : MAT_INVERSION_TRIANG(M)
  !TYPE     : Subroutine
  !PURPOSE  : Invert a general triangular m*n matrix using LAPACK library 
  !COMMENT  : M is destroyed and replaces by its inverse M^-1
  !           M on output is the square matrix n*n
  !+-----------------------------------------------------------------+
  subroutine D_mat_invertTRIANG(A,uplo,diag)
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
    if(info/=0)call error("Error MATRIX/D_mat_invertTRIANG: dtrtri")
  end subroutine D_mat_invertTRIANG
  !-----------------------------

  !-----------------------------
  subroutine Z_mat_invertTRIANG(A,uplo,diag)
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
    if(info/=0)call error("Error MATRIX/D_mat_invertTRIANG: ztrtri")
  end subroutine Z_mat_invertTRIANG





  !******************************************************************
  !******************************************************************
  !******************************************************************











  !+-----------------------------------------------------------------+
  !PROGRAM  : MAT_INVERSION_SYM(M,n)
  !TYPE     : Subroutine
  !PURPOSE  : Invert a symmetric n*n matrix using LAPACK library 
  !COMMENT  : M is destroyed and replaces by its inverse M^-1
  !+-----------------------------------------------------------------+
  subroutine D_mat_invertSYM(A,uplo)
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
    if(info/=0)call error("Error MATRIX/D_mat_invertSYM: 1st call dsytrf")
    lwork=lwork_guess(1);allocate(work(lwork))
    call dsytrf(uplo_,n,A,lda,ipvt,work,lwork,info)
    if(info/=0)call error("Error MATRIX/D_mat_invertSYM: 2nd call dsytrf")
    call dsytri(uplo_,n,A,lda,ipvt,work,info)
    if(info/=0)call error("Error MATRIX/D_mat_invertSYM: dsytri")
    deallocate(ipvt,work)
  end subroutine D_mat_invertSYM
  !-----------------------------

  !-----------------------------
  subroutine Z_mat_invertSYM(A,uplo)
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
    if(info/=0)call error("Error MATRIX/D_mat_invertSYM: 1st call zsytrf")
    lwork=lwork_guess(1);allocate(work(lwork))
    call zsytrf(uplo_,n,A,lda,ipvt,work,lwork,info)
    if(info/=0)call error("Error MATRIX/D_mat_invertSYM: 2nd call zsytrf")
    call zsytri(uplo_,n,A,lda,ipvt,work,info)
    if(info/=0)call error("Error MATRIX/D_mat_invertSYM: zsytri")
    deallocate(ipvt,work)
  end subroutine Z_mat_invertSYM





  !******************************************************************
  !******************************************************************
  !******************************************************************









  !+-----------------------------------------------------------------+
  !PROGRAM  : MAT_INVERSION_HER(M,n)
  !TYPE     : Subroutine
  !PURPOSE  : Invert a hermitian n*n matrix using MKL_LAPACK library 
  !COMMENT  : M is destroyed and replaces by its inverse M^-1
  !+-----------------------------------------------------------------+
  subroutine Z_mat_invertHER(A,uplo)
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
    if(bool)call error("Error MATRIX/Z_mat_invertHER: A not Hermitian")
    !
    call zhetrf(uplo_,n,A,lda,ipvt,lwork_guess,-1,info)
    if(info/=0)call error("Error MATRIX/Z_mat_invertHER: 1st call zhetrf")
    lwork=lwork_guess(1) ; allocate(work(lwork))
    call zhetrf(uplo_,n,A,lda,ipvt,work,lwork,info)
    if(info/=0)call error("Error MATRIX/Z_mat_invertHERE: 2nd call zhetrf")
    call zhetri(uplo_,n,A,lda,ipvt,work,info)
    if(info/=0)call error("Error MATRIX/Z_mat_invertHERE: zhetri")
    deallocate(ipvt,work)
    if(uplo_=="U")then
       forall(i=1:size(A,1),j=1:size(A,2),i>j)A(i,j)=conjg(A(j,i))
    elseif(uplo_=="L")then
       forall(i=1:size(A,1),j=1:size(A,2),i<j)A(i,j)=conjg(A(j,i))
    endif
  end subroutine Z_mat_invertHER





  !******************************************************************
  !******************************************************************
  !******************************************************************








  !+-----------------------------------------------------------------+
  !PURPOSE  : Linear equation solution by Gauss-Jordan elimination
  ! a is an N x N input coefficient matrix. On output, a is replaced 
  !by its matrix inverse.
  !+-----------------------------------------------------------------+
  SUBROUTINE D_mat_invertGJ(a)
    REAL(8), DIMENSION(:,:), INTENT(INOUT) :: a
    INTEGER, DIMENSION(SIZE(a,1))      :: ipiv,indxr,indxc
    !These arrays are used for bookkeeping on the pivoting.
    !INTEGER                            :: nn
    LOGICAL, DIMENSION(SIZE(a,1))      :: lpiv
    REAL(8)                            :: pivinv
    REAL(8), DIMENSION(SIZE(a,1))      :: dumc
    INTEGER, TARGET                    :: irc(2)
    INTEGER                            :: i,l,n
    INTEGER, POINTER                   :: irow,icol
    real(8)                            :: zero=0.d0, one=1.d0
    n=SIZE(a,1)
    irow => irc(1)
    icol => irc(2)
    ipiv=0
    DO i=1,n
       !Main loop over columns to be reduced.
       lpiv = (ipiv == 0)
       !Begin search for a pivot element.
       irc=MAXLOC(ABS(a),outerand(lpiv,lpiv))
       ipiv(icol)=ipiv(icol)+1
       IF (ipiv(icol) > 1) STOP 'gaussj:singular matrix (1)'
       !We now have the pivot element, so we interchange
       !rows, if needed, to put the pivot element on the diagonal. The columns
       !are not physically interchanged, only relabeled:
       !indxc(i),the column of the ith pivot element, is the ith column that is
       !reduced, while indxr(i) is the row in which that pivot element was
       !originally located. If indxr(i) = indxc(i) there is an implied column
       !interchange. With this form of bookkeeping, the inverse matrix will be
       !scrambled by
       !columns.
       IF (irow /= icol) CALL swap(a(irow,:),a(icol,:))
       indxr(i)=irow !We are now ready to divide the pivot row by the pivot element,
       !located at irow and icol.
       indxc(i)=icol
       IF (a(icol,icol) == zero) STOP 'gaussj:singular matrix (2)'
       pivinv=one/a(icol,icol)
       a(icol,icol)= one !CMPLX(one,zero)
       a(icol,:)=a(icol,:)*pivinv
       dumc=a(:,icol)
       !Next, we reduce the rows, except for the pivot one, of course.
       a(:,icol)     = zero !CMPLX
       a(icol,icol)  = pivinv
       a(1:icol-1,:) = a(1:icol-1,:) - outerprod(dumc(1:icol-1),a(icol,:))
       a(icol+1:,:)  = a(icol+1:,:)  - outerprod(dumc(icol+1:),a(icol,:))
    END DO
    !It only remains to unscramble the solution in view of the column
    !interchanges.
    !We do this by interchanging pairs of columns in the reverse order that the
    !permutation
    !was built up.
    DO l=n,1,-1
       CALL swap(a(:,indxr(l)),a(:,indxc(l)))
    END DO
  CONTAINS
    FUNCTION outerand(a,b)
      IMPLICIT NONE
      LOGICAL, DIMENSION(:), INTENT(IN)   :: a,b
      LOGICAL, DIMENSION(SIZE(a),SIZE(b)) :: outerand
      outerand = SPREAD(a,dim=2,ncopies=SIZE(b)).AND.SPREAD(b,dim=1,ncopies=SIZE(a))
    END FUNCTION outerand
    FUNCTION outerprod(a,b)
      real(8), DIMENSION(:), INTENT(IN) :: a,b
      real(8), DIMENSION(size(a),size(b)) :: outerprod
      outerprod = spread(a,dim=2,ncopies=size(b)) * &
           spread(b,dim=1,ncopies=size(a))
    END FUNCTION outerprod
  END SUBROUTINE D_mat_invertGJ
  !-----------------------------

  !-----------------------------
  SUBROUTINE Z_mat_invertGJ(a)
    COMPLEX(8), DIMENSION(:,:), INTENT(INOUT) :: a
    INTEGER, DIMENSION(SIZE(a,1))      :: ipiv,indxr,indxc
    !These arrays are used for bookkeeping on the pivoting.
    !INTEGER                            :: nn
    LOGICAL, DIMENSION(SIZE(a,1))      :: lpiv
    COMPLEX(8)                         :: pivinv
    COMPLEX(8), DIMENSION(SIZE(a,1))   :: dumc
    INTEGER, TARGET                    :: irc(2)
    INTEGER                            :: i,l,n
    INTEGER, POINTER                   :: irow,icol
    real(8)                            :: zero=0.d0, one=1.d0
    n=SIZE(a,1)
    irow => irc(1)
    icol => irc(2)
    ipiv=0
    DO i=1,n
       !Main loop over columns to be reduced.
       lpiv = (ipiv == 0)
       !Begin search for a pivot element.
       irc=MAXLOC(ABS(a),outerand(lpiv,lpiv))
       ipiv(icol)=ipiv(icol)+1
       IF (ipiv(icol) > 1) STOP 'gaussj:singular matrix (1)'
       !We now have the pivot element, so we interchange
       !rows, if needed, to put the pivot element on the diagonal. The columns
       !are not physically interchanged, only relabeled:
       !indxc(i),the column of the ith pivot element, is the ith column that is
       !reduced, while indxr(i) is the row in which that pivot element was
       !originally located. If indxr(i) = indxc(i) there is an implied column
       !interchange. With this form of bookkeeping, the inverse matrix will be
       !scrambled by
       !columns.
       IF (irow /= icol) CALL swap(a(irow,:),a(icol,:))
       indxr(i)=irow !We are now ready to divide the pivot row by the pivot element,
       !located at irow and icol.
       indxc(i)=icol
       IF (a(icol,icol) == zero) STOP 'gaussj:singular matrix (2)'
       pivinv=one/a(icol,icol)
       a(icol,icol)=CMPLX(one,zero,8)
       a(icol,:)=a(icol,:)*pivinv
       dumc=a(:,icol)
       !Next, we reduce the rows, except for the pivot one, of course.
       a(:,icol)     = CMPLX(zero,zero,8)
       a(icol,icol)  = pivinv
       a(1:icol-1,:) = a(1:icol-1,:) - outerprod(dumc(1:icol-1),a(icol,:))
       a(icol+1:,:)  = a(icol+1:,:)  - outerprod(dumc(icol+1:),a(icol,:))
    END DO
    !It only remains to unscramble the solution in view of the column
    !interchanges.
    !We do this by interchanging pairs of columns in the reverse order that the
    !permutation
    !was built up.
    DO l=n,1,-1
       CALL swap(a(:,indxr(l)),a(:,indxc(l)))
    END DO
  CONTAINS
    FUNCTION outerand(a,b)
      IMPLICIT NONE
      LOGICAL, DIMENSION(:), INTENT(IN)   :: a,b
      LOGICAL, DIMENSION(SIZE(a),SIZE(b)) :: outerand
      outerand = SPREAD(a,dim=2,ncopies=SIZE(b)).AND.SPREAD(b,dim=1,ncopies=SIZE(a))
    END FUNCTION outerand
    FUNCTION outerprod(a,b)
      complex(8), DIMENSION(:), INTENT(IN) :: a,b
      complex(8), DIMENSION(size(a),size(b)) :: outerprod
      outerprod = spread(a,dim=2,ncopies=size(b)) * &
           spread(b,dim=1,ncopies=size(a))
    END FUNCTION outerprod
  END SUBROUTINE Z_mat_invertGJ





  !******************************************************************
  !******************************************************************
  !******************************************************************









  !+-----------------------------------------------------------------+
  !PROGRAM  : SWAP
  !TYPE     : Subroutine
  !PURPOSE  : 
  !COMMENT  : 
  !+-----------------------------------------------------------------+
  SUBROUTINE swap_i(a,b)
    INTEGER, INTENT(INOUT) :: a,b
    INTEGER :: dum
    dum=a
    a=b
    b=dum
  END SUBROUTINE swap_i
  !-----------------------------
  !-----------------------------
  !-----------------------------
  SUBROUTINE swap_r(a,b)
    REAL(8), INTENT(INOUT) :: a,b
    REAL(8) :: dum
    dum=a
    a=b
    b=dum
  END SUBROUTINE swap_r
  !-----------------------------
  !-----------------------------
  !-----------------------------
  SUBROUTINE swap_rv(a,b)
    REAL(8), DIMENSION(:), INTENT(INOUT) :: a,b
    REAL(8), DIMENSION(SIZE(a)) :: dum
    dum=a
    a=b
    b=dum
  END SUBROUTINE swap_rv
  !-----------------------------
  !-----------------------------
  !-----------------------------
  SUBROUTINE swap_z(a,b)
    COMPLEX(8), INTENT(INOUT) :: a,b
    COMPLEX(8) :: dum
    dum=a
    a=b
    b=dum
  END SUBROUTINE swap_z
  !-----------------------------
  !-----------------------------
  !-----------------------------
  SUBROUTINE swap_zv(a,b)
    COMPLEX(8), DIMENSION(:), INTENT(INOUT) :: a,b
    COMPLEX(8), DIMENSION(SIZE(a)) :: dum
    dum=a
    a=b
    b=dum
  END SUBROUTINE swap_zv
  !-----------------------------
  !-----------------------------
  !-----------------------------
  SUBROUTINE swap_zm(a,b)
    COMPLEX(8), DIMENSION(:,:), INTENT(INOUT) :: a,b
    COMPLEX(8), DIMENSION(size(a,1),size(a,2)) :: dum
    dum=a
    a=b
    b=dum
  END SUBROUTINE swap_zm




  !******************************************************************
  !******************************************************************
  !******************************************************************




end module MATRIX
