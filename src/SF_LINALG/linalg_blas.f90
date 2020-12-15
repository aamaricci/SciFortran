subroutine d_matmul(A,B,C,alfa,beta)
  real(8),dimension(:,:),intent(inout) :: A ![N,K]
  real(8),dimension(:,:),intent(inout) :: B ![K,M]
  real(8),dimension(:,:),intent(inout) :: C ![N,M]
  real(8),optional                     :: alfa,beta
  real(8)                              :: alfa_,beta_
  integer                              :: N,K,M
  !
  ! C = alfa*A*B + beta*C
  !
  alfa_     = 1d0 ;if(present(alfa))alfa_=alfa
  beta_     = 0d0 ;if(present(beta))beta_=beta
  !
  !
  N = size(A,1)
  K = size(A,2) !==size(B,1)
  M = size(B,2)
  if(any(shape(B)/=[K,M]))stop "d_matmul error: B has illegal shape"
  if(any(shape(C)/=[N,M]))stop "d_matmul error: C has illegal shape"
  !
  call DGEMM('N', 'N', N, M, K, alfa_, A, N, B, K, beta_, C, M)
  !
  return
end subroutine d_matmul

subroutine z_matmul(A,B,C,alfa,beta)
  complex(8),dimension(:,:),intent(inout) :: A ![N,K]
  complex(8),dimension(:,:),intent(inout) :: B ![K,M]
  complex(8),dimension(:,:),intent(inout) :: C ![N,M]
  complex(8),optional                     :: alfa,beta
  complex(8)                              :: alfa_,beta_
  integer                                 :: N,K,M
  !
  ! C = alfa*A*B + beta*C
  !
  alfa_     = dcmplx(1d0,0d0) ;if(present(alfa))alfa_=alfa
  beta_     = dcmplx(0d0,0d0) ;if(present(beta))beta_=beta
  !
  !
  N = size(A,1)
  K = size(A,2) !==size(B,1)
  M = size(B,2)
  if(any(shape(B)/=[K,M]))stop "z_matmul error: B has illegal shape"
  if(any(shape(C)/=[N,M]))stop "z_matmul error: C has illegal shape"
  !
  call ZGEMM('N', 'N', N, M, K, alfa_, A, N, B, K, beta_, C, M)
  !
  return
end subroutine z_matmul




!############## OVERLOAD MATMUL OPERATOR --> .x. #################


function d_matmul_(A,B) result(C)
  real(8),dimension(:,:),intent(in)      :: A ![N,K]
  real(8),dimension(:,:),intent(in)      :: B ![K,M]
  real(8),dimension(size(A,1),size(B,2)) :: C ![N,M]
  integer                                :: N,K,M
  !
  ! C = alfa*A*B + beta*C
  !
  N = size(A,1)
  K = size(A,2) !==size(B,1)
  M = size(B,2)
  if(any(shape(B)/=[K,M]))stop "d_matmul error: B has illegal shape"
  if(any(shape(C)/=[N,M]))stop "d_matmul error: C has illegal shape"
  !
  call DGEMM('N', 'N', N, M, K, 1d0, A, N, B, K, 0d0, C, M)
  !
  return
end function d_matmul_

function z_matmul_(A,B) result(C)
  complex(8),dimension(:,:),intent(in)      :: A ![N,K]
  complex(8),dimension(:,:),intent(in)      :: B ![K,M]
  complex(8),dimension(size(A,1),size(B,2)) :: C ![N,M]
  integer                                :: N,K,M
  !
  ! C = alfa*A*B + beta*C
  !
  N = size(A,1)
  K = size(A,2) !==size(B,1)
  M = size(B,2)
  if(any(shape(B)/=[K,M]))stop "d_matmul error: B has illegal shape"
  if(any(shape(C)/=[N,M]))stop "d_matmul error: C has illegal shape"
  !
  call ZGEMM('N', 'N', N, M, K, dcmplx(1d0,0d0), A, N, B, K, dcmplx(0d0,0d0), C, M)
  !
  return
end function z_matmul_




