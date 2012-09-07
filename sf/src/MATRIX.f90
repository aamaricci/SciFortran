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

  interface mat_inversion
     module procedure D_mat_invert, Z_mat_invert
  end interface mat_inversion

  interface mat_inversion_sym
     module procedure d_mat_invertsym, z_mat_invertsym
  end interface mat_inversion_sym

  interface mat_inversion_her
     module procedure z_mat_inverther
  end interface mat_inversion_her

  interface mat_inversion_triang
     module procedure d_mat_inverttriang, z_mat_inverttriang
  end interface mat_inversion_triang

  interface mat_inversion_gj
     module procedure d_mat_invertgj, z_mat_invertgj
  end interface mat_inversion_gj

  !----------------------------------------------------------------------------------

  interface invert_matrix_gj
     module procedure d_invert_matrix_gj,c_invert_matrix_gj
  end interface invert_matrix_gj

  interface invert_matrix
     module procedure d_invert_matrix,c_invert_matrix
  end interface invert_matrix

  !----------------------------------------------------------------------------------

  interface mat_diagonalization
     module procedure d_mat_diagonalization,z_mat_diagonalization
  end interface mat_diagonalization

  !----------------------------------------------------------------------------------

  interface solve_linear_system
     module procedure d_mat_solve_linear_system_1rhs,z_mat_solve_linear_system_1rhs,&
          d_mat_solve_linear_system_Mrhs,z_mat_solve_linear_system_Mrhs
  end interface solve_linear_system


  !from nr90
  interface swap
     module procedure swap_i,swap_r,swap_rv,swap_z,swap_zv,swap_zm
  end interface swap

  public :: mat_diagonalization
  !
  public :: solve_linear_system
  !
  public :: mat_inversion
  public :: mat_inversion_sym
  public :: mat_inversion_her
  public :: mat_inversion_triang
  public :: mat_inversion_gj
  !
  public :: invert_matrix
  public :: invert_matrix_gj
  ! 


contains



  !+-----------------------------------------------------------------+
  !PROGRAM  : MAT_DIAGONALIZATION
  !PURPOSE  : function interface to matrix inversion subroutine
  !+-----------------------------------------------------------------+
  subroutine d_mat_diagonalization(M,E,jobz,uplo)
    real(8),dimension(:,:),intent(inout)       :: M
    real(8),dimension(size(M,1)),intent(inout) :: E
    character(len=1),optional                  :: jobz,uplo
    character(len=1)                           :: jobz_,uplo_
    integer                                    :: i,j,Nsize,info
    integer                                    :: lwork
    real(8),dimension(1)                       :: lwork_guess
    real(8),dimension(:),allocatable           :: work
    Nsize=size(M,1)
    if(Nsize/=size(M,2))call error("Error MATRIX/d_mat_diagonalization: size(M,1) /= size(M,2)")
    jobz_='N';if(present(jobz))jobz_=jobz
    uplo_='U';if(present(uplo))uplo_=uplo
    call dsyev(jobz_,uplo_,Nsize,M,Nsize,E,lwork_guess,-1,info)
    if(info /= 0)call error("Error MATRIX/d_mat_diagonalization: 1st call dsyev")
    lwork=lwork_guess(1) ; allocate(work(lwork))
    call dsyev(jobz_,uplo_,Nsize,M,Nsize,E,work,lwork,info)
    if(info /= 0)call error("Error MATRIX/d_mat_diagonalization: 2ns call dsyev")
    deallocate(work)
  end subroutine d_mat_diagonalization
  !-----------------------------

  !-----------------------------
  subroutine z_mat_diagonalization(M,E,jobz,uplo)
    complex(8),dimension(:,:),intent(inout)    :: M
    real(8),dimension(size(M,1)),intent(inout) :: E
    character(len=1),optional                  :: jobz,uplo
    character(len=1)                           :: jobz_,uplo_
    integer                                    :: i,j,Nsize,info
    integer                                    :: lwork
    complex(8),dimension(1)                    :: lwork_guess
    complex(8),dimension(:),allocatable        :: work
    real(8),dimension(:),allocatable           :: rwork
    Nsize=size(M,1)
    if(Nsize/=size(M,2))call error("Error MATRIX/z_mat_diagonalization: size(M,1) /= size(M,2)")
    jobz_='N';if(present(jobz))jobz_=jobz
    uplo_='U';if(present(uplo))uplo_=uplo
    allocate(rwork(3*Nsize))
    call zheev(jobz_,uplo_,Nsize,M,Nsize,E,lwork_guess,-1,rwork,info)
    if(info /= 0)call error("Error MATRIX/d_mat_diagonalization: 1st call zsyev")
    lwork=lwork_guess(1) ; allocate(work(lwork))
    call zheev(jobz_,uplo_,Nsize,M,Nsize,E,work,lwork,rwork,info)
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
    integer                            :: n,nrhs,lda,ldb
    integer                            :: info
    integer,dimension(:),allocatable   :: ipvt
    trans_="N";if(present(trans))trans_=trans
    n=size(A,1) 
    if(n/=size(A,2))call error("Error in MATRIX/d_mat_solve_linear_system:not a square matrix")
    if(n/=size(b))  call error("Error in MATRIX/d_mat_solve_linear_system:b has wrong dimension")    
    allocate(ipvt(n))
    call dgetrf(n,n,A,n,ipvt,info)
    if(info/=0)call error("Error MATRIX/d_mat_solve_linear_system: dgetrf")    
    lda=n ; ldb=n ; nrhs=1
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
    integer                            :: n,nrhs,lda,ldb
    integer                            :: info
    integer,dimension(:),allocatable   :: ipvt
    trans_="N";if(present(trans))trans_=trans
    n=size(A,1) 
    if(n/=size(A,2))call error("Error in MATRIX/d_mat_solve_linear_system:not a square matrix")
    if(n/=size(b))  call error("Error in MATRIX/d_mat_solve_linear_system:b has wrong dimension")    
    allocate(ipvt(n))
    call zgetrf(n,n,A,n,ipvt,info)
    if(info/=0)call error("Error MATRIX/d_mat_solve_linear_system: dgetrf")    
    lda=n ; ldb=n ; nrhs=1
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
    integer                            :: n,nrhs,lda,ldb
    integer                            :: info
    integer,dimension(:),allocatable   :: ipvt
    trans_="N";if(present(trans))trans_=trans
    n=size(A,1) 
    if(n/=size(A,2))call error("Error in MATRIX/d_mat_solve_linear_system:not a square matrix")
    if(n/=size(b,1))call error("Error in MATRIX/d_mat_solve_linear_system:b has wrong dimension")    
    allocate(ipvt(n))
    call dgetrf(n,n,A,n,ipvt,info)
    if(info/=0)call error("Error MATRIX/d_mat_solve_linear_system: dgetrf")    
    lda=n ; ldb=n ; nrhs=size(b,2)
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
    integer                            :: n,nrhs,lda,ldb
    integer                            :: info
    integer,dimension(:),allocatable   :: ipvt
    trans_="N";if(present(trans))trans_=trans
    n=size(A,1) 
    if(n/=size(A,2))call error("Error in MATRIX/d_mat_solve_linear_system:not a square matrix")
    if(n/=size(b,1))call error("Error in MATRIX/d_mat_solve_linear_system:b has wrong dimension")    
    allocate(ipvt(n))
    call zgetrf(n,n,A,n,ipvt,info)
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
    call mat_inversion_GJ(M1)
  end function d_invert_matrix_gj
  function c_invert_matrix_gj(M,L) result(M1)
    integer                   :: L
    complex(8),dimension(L,L) :: M,M1
    M1=M
    call mat_inversion_GJ(M1)
  end function c_invert_matrix_gj
  !
  function d_invert_matrix(M,L) result(M1)
    integer                   :: L
    real(8),dimension(L,L)    :: M,M1
    M1=M
    call mat_inversion(M1)
  end function d_invert_matrix
  function c_invert_matrix(M,L) result(M1)
    integer                   :: L
    complex(8),dimension(L,L) :: M,M1
    M1=M
    call mat_inversion(M1)
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
  subroutine D_mat_invert(M)
    integer                          :: ndim,ndim1,ndim2
    integer                          :: info
    integer                          :: lwork
    integer,dimension(:),allocatable :: ipvt
    real(8),dimension(:,:)           :: M
    real(8),dimension(:),allocatable :: work
    real(8),dimension(1)             :: lwork_guess
    ndim1=size(M,1) ; ndim2=size(M,2) 
    if(ndim1<max(1,ndim2))call error("Error in MATRIX/D_mat_inversion: ndim1<ndim2")
    ndim=max(ndim1,ndim2)
    allocate(ipvt(ndim))
    call dgetrf(ndim1,ndim2,M,ndim1,ipvt,info)
    if(info/=0)call error("Error MATRIX/D_mat_invert: dgetrf")
    call dgetri(ndim2,M,ndim1,ipvt,lwork_guess,-1,info)
    if(info/=0)call error("Error MATRIX/D_mat_invert: 1st call dgetri")
    lwork=lwork_guess(1) ; allocate(work(lwork))
    call dgetri(ndim2,M,ndim1,ipvt,work,lwork,info)
    if(info/=0)call error("Error MATRIX/D_mat_invert: 2nd call dgetri")
    deallocate(ipvt,work)
  end subroutine D_mat_invert
  !-----------------------------

  !-----------------------------
  subroutine Z_mat_invert(M)
    integer                              :: ndim,ndim1,ndim2
    integer                              :: info
    integer                              :: lwork
    integer,dimension(:),allocatable     :: ipvt
    complex(8),dimension(:,:)            :: M
    complex(8),dimension(:),allocatable  :: work
    complex(8),dimension(1)              :: lwork_guess
    ndim1=size(M,1) ; ndim2=size(M,2) 
    if(ndim1<max(1,ndim2))call error("Error in MATRIX/Z_mat_inversion: ndim1<ndim2")
    ndim=max(ndim1,ndim2)
    lwork=ndim*8
    allocate(ipvt(ndim))
    call zgetrf(ndim1,ndim2,M,ndim1,ipvt,info)
    if(info/=0)call error("Error MATRIX/Z_mat_invert: zgetrf")
    call zgetri(ndim2,M,ndim1,ipvt,lwork_guess,-1,info)
    if(info/=0)call error("Error MATRIX/Z_mat_invert: 1st call dgetri")
    lwork=lwork_guess(1) ; allocate(work(lwork))
    call zgetri(ndim2,M,ndim1,ipvt,work,lwork,info)
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
  subroutine D_mat_invertTRIANG(M,char)
    character(len=*),optional        :: char
    character(len=1)                 :: uplo
    character(len=1)                 :: diag
    integer                          :: ndim1,ndim2
    integer                          :: info
    real(8),dimension(:,:)           :: M
    ndim1=size(M,1);ndim2=size(M,2)
    if(ndim1<max(1,ndim2))stop "Error in mat_inversion: ndim1<ndim2"
    uplo="U";if(present(char))uplo=trim(char)
    diag="N" !not a unit triangular matrix
    call dtrtri(uplo,diag,ndim2,M,ndim1,info)
    if(info/=0)stop"Error in dtrtri"
  end subroutine D_mat_invertTRIANG
  !-----------------------------

  !-----------------------------
  subroutine Z_mat_invertTRIANG(M,char)
    character(len=*),optional        :: char
    character(len=1)                 :: uplo
    character(len=1)                 :: diag
    integer                          :: ndim1,ndim2
    integer                          :: info
    complex(8),dimension(:,:)        :: M
    ndim1=size(M,1);ndim2=size(M,2)
    if(ndim1<max(1,ndim2))stop "Error in mat_inversion: ndim1<ndim2"
    uplo="U";if(present(char))uplo=trim(char)
    diag="N" !not a unit triangular matrix
    call ztrtri(uplo,diag,ndim2,M,ndim1,info)
    if(info/=0)stop "Error in ztrtri"
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
  subroutine D_mat_invertSYM(M,ndim,char)
    character(len=*),optional        :: char
    character(len=10)                :: uplo
    integer                          :: ndim
    integer                          :: info
    integer                          :: lwork
    integer,dimension(ndim)          :: ipvt
    real(8),dimension(ndim,ndim)     :: M
    real(8),dimension(:),allocatable :: work
    real(8),dimension(1)             :: lwork_guess
    uplo="U";if(present(char))uplo=trim(char)
    call dsytrf(trim(uplo),ndim,M,ndim,ipvt,lwork_guess,-1,info)
    if(info/=0)call error("Error MATRIX/D_mat_invertSYM: 1st call dsytrf")
    lwork=lwork_guess(1);allocate(work(lwork))
    call dsytrf(trim(uplo),ndim,M,ndim,ipvt,work,lwork,info)
    if(info/=0)call error("Error MATRIX/D_mat_invertSYM: 2nd call dsytrf")
    call dsytri(trim(uplo),ndim,M,ndim,ipvt,work,info)
    if(info/=0)call error("Error MATRIX/D_mat_invertSYM: dsytri")
    deallocate(work)
  end subroutine D_mat_invertSYM
  !-----------------------------

  !-----------------------------
  subroutine Z_mat_invertSYM(M,ndim,char)
    character(len=*),optional    :: char
    character(len=10)            :: uplo
    integer                      :: ndim
    integer                      :: info
    integer                      :: lwork
    integer,dimension(ndim)      :: ipvt
    complex(8),dimension(ndim,ndim) :: M
    complex(8),dimension(:),allocatable :: work
    complex(8),dimension(1)             :: lwork_guess
    uplo="U";if(present(char))uplo=trim(char)
    call zsytrf(trim(uplo),ndim,M,ndim,ipvt,lwork_guess,-1,info)
    if(info/=0)call error("Error MATRIX/Z_mat_invertSYM: 1st call zsytrf")
    lwork=lwork_guess(1);allocate(work(lwork))
    call zsytrf(trim(uplo),ndim,M,ndim,ipvt,work,lwork,info)
    if(info/=0)call error("Error MATRIX/Z_mat_invertSYM: 2nd call zsytrf")
    call zsytri(trim(uplo),ndim,M,ndim,ipvt,work,info)
    if(info/=0)call error("Error MATRIX/Z_mat_invertSYM: zsytri")
    deallocate(work)
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
  subroutine Z_mat_invertHER(M,ndim,char)
    character(len=*),optional    :: char
    character(len=10)            :: uplo
    integer                      :: ndim
    integer                      :: info
    integer                      :: lwork
    integer,dimension(ndim)      :: ipvt
    complex(8),dimension(ndim,ndim) :: M
    complex(8),dimension(:),allocatable :: work
    complex(8),dimension(1)             :: lwork_guess
    uplo="U";if(present(char))uplo=trim(char)
    call zhetrf(trim(uplo),ndim,M,ndim,ipvt,lwork_guess,-1,info)
    if(info/=0)call error("Error MATRIX/Z_mat_invertSYM: 1st call zhetrf")
    lwork=lwork_guess(1) ; allocate(work(lwork))
    call zhetrf(trim(uplo),ndim,M,ndim,ipvt,work,lwork,info)
    if(info/=0)call error("Error MATRIX/Z_mat_invertSYM: 2nd call zhetrf")
    call zhetri(trim(uplo),ndim,M,ndim,ipvt,work,info)
    if(info/=0)call error("Error MATRIX/Z_mat_invertSYM: zhetri")
    deallocate(work)
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
