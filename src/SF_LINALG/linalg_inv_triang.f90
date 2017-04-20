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
