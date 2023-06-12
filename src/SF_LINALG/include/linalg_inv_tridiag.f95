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
