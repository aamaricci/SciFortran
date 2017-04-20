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
