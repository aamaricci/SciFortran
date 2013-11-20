subroutine sreadV_II(pname,X,Y1,Y2,Y3,Y4,Y5,Y6,Y7,Y8)
  integer                             :: i,Np
  character(len=*)                    :: pname
  integer,dimension(:)                :: X
  integer,dimension(size(X))          :: Y1
  integer,dimension(size(X)),optional :: Y2,Y3,Y4,Y5,Y6,Y7,Y8
  Np=size(X)
  open(719,file=adjustl(trim(pname)))
  if(present(Y8))then
     do i=1,Np
        read(719,"(9(I15))")X(i),Y1(i),Y2(i),Y3(i),Y4(i),Y5(i),Y6(i),Y7(i),Y8(i)
     enddo
  elseif(present(Y7))then
     do i=1,Np
        read(719,"(8(I15))")X(i),Y1(i),Y2(i),Y3(i),Y4(i),Y5(i),Y6(i),Y7(i)
     enddo
  elseif(present(Y6))then
     do i=1,Np
        read(719,"(7(I15))")X(i),Y1(i),Y2(i),Y3(i),Y4(i),Y5(i),Y6(i)
     enddo
  elseif(present(Y5))then
     do i=1,Np
        read(719,"(6(I15))")X(i),Y1(i),Y2(i),Y3(i),Y4(i),Y5(i)
     enddo
  elseif(present(Y4))then
     do i=1,Np
        read(719,"(5(I15))")X(i),Y1(i),Y2(i),Y3(i),Y4(i)
     enddo
  elseif(present(Y3))then
     do i=1,Np
        read(719,"(4(I15))")X(i),Y1(i),Y2(i),Y3(i)
     enddo
  elseif(present(Y2))then
     do i=1,Np
        read(719,*)X(i),Y1(i),Y2(i)
     enddo
  else
     do i=1,Np
        read(719,*)X(i),Y1(i)
     enddo
  endif
  close(719)
end subroutine sreadV_II

!----------------------------

subroutine sreadV_IR(pname,X,Y1,Y2,Y3,Y4,Y5,Y6,Y7,Y8)
  integer                             :: i,Np
  character(len=*)                    :: pname
  integer,dimension(:)                :: X
  real(8),dimension(size(X))          :: Y1
  real(8),dimension(size(X)),optional :: Y2,Y3,Y4,Y5,Y6,Y7,Y8
  Np=size(X)
  open(719,file=adjustl(trim(pname)))

  if(present(Y8))then
     do i=1,Np
        read(719,"(I15,8(F21.12))")X(i),Y1(i),Y2(i),Y3(i),Y4(i),Y5(i),Y6(i),Y7(i),Y8(i)
     enddo
  elseif(present(Y7))then
     do i=1,Np
        read(719,"(I15,7(F21.12))")X(i),Y1(i),Y2(i),Y3(i),Y4(i),Y5(i),Y6(i),Y7(i)
     enddo
  elseif(present(Y6))then
     do i=1,Np
        read(719,"(I15,6(F21.12))")X(i),Y1(i),Y2(i),Y3(i),Y4(i),Y5(i),Y6(i)
     enddo
  elseif(present(Y5))then
     do i=1,Np
        read(719,"(I15,5(F21.12))")X(i),Y1(i),Y2(i),Y3(i),Y4(i),Y5(i)
     enddo
  elseif(present(Y4))then
     do i=1,Np
        read(719,"(I15,4(F21.12))")X(i),Y1(i),Y2(i),Y3(i),Y4(i)
     enddo
  elseif(present(Y3))then
     do i=1,Np
        read(719,"(I15,3(F21.12))")X(i),Y1(i),Y2(i),Y3(i)
     enddo
  elseif(present(Y2))then
     do i=1,Np
        read(719,*)X(i),Y1(i),Y2(i)
     enddo
  else
     do i=1,Np
        read(719,*)X(i),Y1(i)
     enddo
  endif

  close(719)
end subroutine sreadV_IR

!----------------------------

subroutine sreadV_IC(pname,X,Y1,Y2,Y3,Y4)
  integer                                :: i,Np
  character(len=*)                       :: pname
  integer,dimension(:)                   :: X
  complex(8),dimension(size(X))          :: Y1
  complex(8),dimension(size(X)),optional :: Y2,Y3,Y4
  real(8),allocatable,dimension(:,:)     :: reY,imY

  Np=size(X)
  open(719,file=adjustl(trim(pname)))
  if(present(Y4))then
     allocate(reY(4,size(X)),imY(4,size(X)))
     do i=1,Np
        read(719,"(I15,8(F21.12))")X(i),reY(1,i),imY(1,i),reY(2,i),imY(2,i),reY(3,i),imY(3,i),reY(4,i),imY(4,i)
     enddo
     Y1=cmplx(reY(1,:),imY(1,:),8)
     Y2=cmplx(reY(2,:),imY(2,:),8)
     Y3=cmplx(reY(3,:),imY(3,:),8)
     Y4=cmplx(reY(4,:),imY(4,:),8)
  elseif(present(Y3))then
     allocate(reY(3,size(X)),imY(3,size(X)))
     do i=1,Np
        read(719,"(I15,6(F21.12))")X(i),reY(1,i),imY(1,i),reY(2,i),imY(2,i),reY(3,i),imY(3,i)
     enddo
     Y1=cmplx(reY(1,:),imY(1,:),8)
     Y2=cmplx(reY(2,:),imY(2,:),8)
     Y3=cmplx(reY(3,:),imY(3,:),8)
  elseif(present(Y2))then
     allocate(reY(2,size(X)),imY(2,size(X)))
     do i=1,Np
        read(719,"(I15,4(F21.12))")X(i),reY(1,i),imY(1,i),reY(2,i),imY(2,i)
     enddo
     Y1=cmplx(reY(1,:),imY(1,:),8)
     Y2=cmplx(reY(2,:),imY(2,:),8)
  else
     allocate(reY(1,size(X)),imY(1,size(X)))
     do i=1,Np
        read(719,*)X(i),reY(1,i),imY(1,i)
     enddo
     Y1=cmplx(reY(1,:),imY(1,:),8)
  endif
  close(719)
end subroutine sreadV_IC

!----------------------------

subroutine sreadV_RI(pname,X,Y1,Y2,Y3,Y4,Y5,Y6,Y7,Y8)
  integer                             :: i,Np
  character(len=*)                    :: pname
  real(8),dimension(:)                :: X
  integer,dimension(size(X))          :: Y1
  integer,dimension(size(X)),optional :: Y2,Y3,Y4,Y5,Y6,Y7,Y8
  Np=size(X)
  open(719,file=adjustl(trim(pname)))
  if(present(Y8))then
     do i=1,Np
        read(719,"(F21.12,8(I15))")X(i),Y1(i),Y2(i),Y3(i),Y4(i),Y5(i),Y6(i),Y7(i),Y8(i)
     enddo
  elseif(present(Y7))then
     do i=1,Np
        read(719,"(F21.12,7(I15))")X(i),Y1(i),Y2(i),Y3(i),Y4(i),Y5(i),Y6(i),Y7(i)
     enddo
  elseif(present(Y6))then
     do i=1,Np
        read(719,"(F21.12,6(I15))")X(i),Y1(i),Y2(i),Y3(i),Y4(i),Y5(i),Y6(i)
     enddo
  elseif(present(Y5))then
     do i=1,Np
        read(719,"(F21.12,5(I15))")X(i),Y1(i),Y2(i),Y3(i),Y4(i),Y5(i)
     enddo
  elseif(present(Y4))then
     do i=1,Np
        read(719,"(F21.12,4(I15))")X(i),Y1(i),Y2(i),Y3(i),Y4(i)
     enddo
  elseif(present(Y3))then
     do i=1,Np
        read(719,"(F21.12,3(I15))")X(i),Y1(i),Y2(i),Y3(i)
     enddo
  elseif(present(Y2))then
     do i=1,Np
        read(719,*)X(i),Y1(i),Y2(i)
     enddo
  else
     do i=1,Np
        read(719,*)X(i),Y1(i)
     enddo
  endif
  close(719)
end subroutine sreadV_RI

!----------------------------

subroutine sreadV_RR(pname,X,Y1,Y2,Y3,Y4,Y5,Y6,Y7,Y8)
  integer                             :: i,Np
  character(len=*)                    :: pname
  real(8),dimension(:)                :: X
  real(8),dimension(size(X))          :: Y1
  real(8),dimension(size(X)),optional :: Y2,Y3,Y4,Y5,Y6,Y7,Y8
  Np=size(X)
  open(719,file=adjustl(trim(pname)))
  if(present(Y8))then
     do i=1,Np
        read(719,"(9(F21.12))")X(i),Y1(i),Y2(i),Y3(i),Y4(i),Y5(i),Y6(i),Y7(i),Y8(i)
     enddo
  elseif(present(Y7))then
     do i=1,Np
        read(719,"(8(F21.12))")X(i),Y1(i),Y2(i),Y3(i),Y4(i),Y5(i),Y6(i),Y7(i)
     enddo
  elseif(present(Y6))then
     do i=1,Np
        read(719,"(7(F21.12))")X(i),Y1(i),Y2(i),Y3(i),Y4(i),Y5(i),Y6(i)
     enddo
  elseif(present(Y5))then
     do i=1,Np
        read(719,"(6(F21.12))")X(i),Y1(i),Y2(i),Y3(i),Y4(i),Y5(i)
     enddo
  elseif(present(Y4))then
     do i=1,Np
        read(719,"(5(F21.12))")X(i),Y1(i),Y2(i),Y3(i),Y4(i)
     enddo
  elseif(present(Y3))then
     do i=1,Np
        read(719,"(4(F21.12))")X(i),Y1(i),Y2(i),Y3(i)
     enddo
  elseif(present(Y2))then
     do i=1,Np
        read(719,*)X(i),Y1(i),Y2(i)
     enddo
  else
     do i=1,Np
        read(719,*)X(i),Y1(i)
     enddo
  endif
  close(719)
end subroutine sreadV_RR

!----------------------------

subroutine sreadV_RC(pname,X,Y1,Y2,Y3,Y4)
  integer                                :: i,Np
  character(len=*)                       :: pname
  real(8),dimension(:)                   :: X
  complex(8),dimension(size(X))          :: Y1
  complex(8),dimension(size(X)),optional :: Y2,Y3,Y4
  real(8),allocatable,dimension(:,:)     :: reY,imY

  Np=size(X)
  open(719,file=trim(adjustl(trim(pname))))
  if(present(Y4))then
     allocate(reY(4,size(X)),imY(4,size(X)))
     do i=1,Np
        read(719,"(F21.12,8(F21.12))")X(i),imY(1,i),reY(1,i),imY(2,i),reY(2,i),imY(3,i),reY(3,i),imY(4,i),reY(4,i)
     enddo
     Y1=cmplx(reY(1,:),imY(1,:),8)
     Y2=cmplx(reY(2,:),imY(2,:),8)
     Y3=cmplx(reY(3,:),imY(3,:),8)
     Y4=cmplx(reY(4,:),imY(4,:),8)
  elseif(present(Y3))then
     allocate(reY(3,size(X)),imY(3,size(X)))
     do i=1,Np
        read(719,"(F21.12,6(F21.12))")X(i),imY(1,i),reY(1,i),imY(2,i),reY(2,i),imY(3,i),reY(3,i)
     enddo
     Y1=cmplx(reY(1,:),imY(1,:),8)
     Y2=cmplx(reY(2,:),imY(2,:),8)
     Y3=cmplx(reY(3,:),imY(3,:),8)
  elseif(present(Y2))then
     allocate(reY(2,size(X)),imY(2,size(X)))
     do i=1,Np
        read(719,"(F21.12,4(F21.12))")X(i),imY(1,i),reY(1,i),imY(2,i),reY(2,i)
     enddo
     Y1=cmplx(reY(1,:),imY(1,:),8)
     Y2=cmplx(reY(2,:),imY(2,:),8)
  else
     allocate(reY(1,size(X)),imY(1,size(X)))
     do i=1,Np
        read(719,*)X(i),imY(1,i),reY(1,i)
     enddo
     Y1=cmplx(reY(1,:),imY(1,:),8)
  endif
  close(719)
end subroutine sreadV_RC

!----------------------------


