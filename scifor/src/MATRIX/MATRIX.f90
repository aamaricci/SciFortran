
module MATRIX
  !###############################################################
  ! PROGRAM  : MATRIX
  ! TYPE     : Module
  ! PURPOSE  : A collection of MATRIX related routines
  !###############################################################
  implicit none
  private
  include "mkl_lapack.fi"

  interface mat_inversion
     module procedure D_mat_invert, Z_mat_invert
  end interface mat_inversion

  interface mat_inversion_SYM
     module procedure D_mat_invertSYM, Z_mat_invertSYM
  end interface mat_inversion_SYM

  interface mat_inversion_HER
     module procedure Z_mat_invertHER
  end interface mat_inversion_HER

  interface mat_inversion_TRIANG
     module procedure D_mat_invertTRIANG, Z_mat_invertTRIANG
  end interface mat_inversion_TRIANG

  interface mat_inversion_GJ
     module procedure D_mat_invertGJ, Z_mat_invertGJ
  end interface mat_inversion_GJ

  interface invert_matrix_gj
     module procedure d_invert_matrix_gj,c_invert_matrix_gj
  end interface invert_matrix_gj

  interface invert_matrix
     module procedure d_invert_matrix,c_invert_matrix
  end interface invert_matrix

  !From NR90
  interface swap
     module procedure swap_i,swap_r,swap_rv,swap_z,swap_zv,swap_zm
  end interface swap

  public :: mat_inversion
  public :: mat_inversion_SYM
  public :: mat_inversion_HER
  public :: mat_inversion_TRIANG
  public :: mat_inversion_GJ
  public :: invert_matrix,invert_matrix_gj
contains

  function d_invert_matrix_gj(M,L) result(M1)
    integer :: L
    real(8),dimension(L,L) :: M,M1
    M1=M
    call mat_inversion_GJ(M1)
  end function d_invert_matrix_gj
  function c_invert_matrix_gj(M,L) result(M1)
    integer :: L
    complex(8),dimension(L,L) :: M,M1
    M1=M
    call mat_inversion_GJ(M1)
  end function c_invert_matrix_gj

  function d_invert_matrix(M,L) result(M1)
    integer :: L
    real(8),dimension(L,L) :: M,M1
    M1=M
    call mat_inversion(M1)
  end function d_invert_matrix
  function c_invert_matrix(M,L) result(M1)
    integer :: L
    complex(8),dimension(L,L) :: M,M1
    M1=M
    call mat_inversion(M1)
  end function c_invert_matrix


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
    ndim1=size(M,1) ; ndim2=size(M,2) 
    if(ndim1<max(1,ndim2))stop "Error in mat_inversion: ndim1<ndim2"
    ndim=max(ndim1,ndim2)
    lwork=ndim*8
    allocate(ipvt(ndim),work(lwork))
    call dgetrf(ndim1,ndim2,M,ndim1,ipvt,info)
    if(info/=0)stop "Error in dgetrf"
    call dgetri(ndim2,M,ndim1,ipvt,work,lwork,info)
    if(info/=0)stop "Error in dgetri"
  end subroutine D_mat_invert
  !-----------------------------
  !-----------------------------
  !-----------------------------
  subroutine Z_mat_invert(M)
    integer                              :: ndim,ndim1,ndim2
    integer                              :: info
    integer                              :: lwork
    integer,dimension(:),allocatable     :: ipvt
    complex(8),dimension(:,:)            :: M
    complex(8),dimension(:),allocatable  :: work
    ndim1=size(M,1) ; ndim2=size(M,2) 
    if(ndim1<max(1,ndim2))stop "Error in mat_inversion: ndim1<ndim2"
    ndim=max(ndim1,ndim2)
    lwork=ndim*8
    allocate(ipvt(ndim),work(lwork))
    call zgetrf(ndim1,ndim2,M,ndim1,ipvt,info)
    if(info/=0)stop "Error in zgetrf"
    call zgetri(ndim2,M,ndim1,ipvt,work,lwork,info)
    if(info/=0)stop "Error in zgetri"
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
    integer                          :: ndim,ndim1,ndim2
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
  !-----------------------------
  subroutine Z_mat_invertTRIANG(M,char)
    character(len=*),optional        :: char
    character(len=1)                 :: uplo
    character(len=1)                 :: diag
    integer                          :: ndim,ndim1,ndim2
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
    character(len=*),optional    :: char
    character(len=10)            :: uplo
    integer                      :: ndim
    integer                      :: info
    integer                      :: lwork
    integer,dimension(ndim)      :: ipvt
    real(8),dimension(ndim,ndim) :: M
    real(8),dimension(ndim*8)    :: work
    lwork=ndim*8
    uplo="U";if(present(char))uplo=trim(char)
    call dsytrf(trim(uplo),ndim,M,ndim,ipvt,work,lwork,info)
    if(info/=0)stop "Error in dsytrf"
    call dsytri(trim(uplo),ndim,M,ndim,ipvt,work,info)
    if(info/=0)stop "Error in dsytri"
  end subroutine D_mat_invertSYM
  !-----------------------------
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
    complex(8),dimension(ndim*8)    :: work
    lwork=ndim*8
    uplo="U";if(present(char))uplo=trim(char)
    call zsytrf(trim(uplo),ndim,M,ndim,ipvt,work,lwork,info)
    if(info/=0)stop "Error in zsytrf"
    call zsytri(trim(uplo),ndim,M,ndim,ipvt,work,info)
    if(info/=0)stop "Error in zsytri"
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
    complex(8),dimension(ndim*8)    :: work
    lwork=ndim*8
    uplo="U";if(present(char))uplo=trim(char)
    call zhetrf(trim(uplo),ndim,M,ndim,ipvt,work,lwork,info)
    if(info/=0)stop "Error in zsytrf"
    call zhetri(trim(uplo),ndim,M,ndim,ipvt,work,info)
    if(info/=0)stop "Error in zsytri"
  end subroutine Z_mat_invertHER
  !******************************************************************
  !******************************************************************
  !******************************************************************








  !+-----------------------------------------------------------------+
  !PROGRAM  : 
  !TYPE     : Subroutine
  !PURPOSE  : Linear equation solution by Gauss-Jordan elimination
  ! a is an N x N input coefficient matrix. On output, a is replaced 
  !by its matrix inverse.
  !+-----------------------------------------------------------------+
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
       a(icol,icol)=CMPLX(one,zero)
       a(icol,:)=a(icol,:)*pivinv
       dumc=a(:,icol)
       !Next, we reduce the rows, except for the pivot one, of course.
       a(:,icol)     = CMPLX(zero,zero)
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
  !PROGRAM  : 
  !TYPE     : Subroutine
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
       a(:,icol)     = zero !CMPLX(zero,zero)
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
