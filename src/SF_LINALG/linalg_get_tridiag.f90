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
