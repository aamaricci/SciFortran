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
