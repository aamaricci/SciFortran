subroutine sreadM_II(pname,X,Y1,Y2,Y3,Y4)
  integer                                  :: i,j,Np,Ny1,Ny2
  character(len=*)                         :: pname
  integer,dimension(:,:)                   :: Y1
  integer,dimension(size(Y1,1),size(Y1,2)),optional :: Y2,Y3,Y4
  integer,dimension(size(Y1,2))            :: X
  Ny1=size(Y1,1) ; Ny2=size(Y1,2)
  open(719,file=adjustl(trim(pname)))
  if(present(Y4))then
     do i=1,Ny1
        do j=1,Ny2
           write(719,'(5(I15))')X(j),Y1(i,j),Y2(i,j),Y3(i,j),Y4(i,j)
        enddo
     enddo
  elseif(present(Y3))then
     do i=1,Ny1
        do j=1,Ny2
           write(719,'(4(I15))')X(j),Y1(i,j),Y2(i,j),Y3(i,j)
        enddo
     enddo
  elseif(present(Y2))then
     do i=1,Ny1
        do j=1,Ny2
           write(719,'(3(I15))')X(j),Y1(i,j),Y2(i,j)
        enddo
     enddo
  else
     do i=1,Ny1
        do j=1,Ny2
           read(719,*)X(j),Y1(i,j)
        enddo
     enddo
  endif
  close(719)
end subroutine sreadM_II

subroutine sreadM_IR(pname,X,Y1,Y2,Y3,Y4)
  integer                                :: i,j,Ny1,Ny2
  character(len=*)                       :: pname
  real(8),dimension(:,:)                   :: Y1
  real(8),dimension(size(Y1,1),size(Y1,2)),optional :: Y2,Y3,Y4
  integer,dimension(size(Y1,2))            :: X
  Ny1=size(Y1,1) ; Ny2=size(Y1,2)
  open(719,file=adjustl(trim(pname)))
  if(present(Y4))then
     do i=1,Ny1
        do j=1,Ny2
           write(719,"(I15,4(F21.12))")X(j),Y1(i,j),Y2(i,j),Y3(i,j),Y4(i,j)
        enddo
     enddo
  elseif(present(Y3))then
     do i=1,Ny1
        do j=1,Ny2
           write(719,"(I15,3(F21.12))")X(j),Y1(i,j),Y2(i,j),Y3(i,j)
        enddo
     enddo
  elseif(present(Y2))then
     do i=1,Ny1
        do j=1,Ny2
           write(719,"(I15,2(F21.12))")X(j),Y1(i,j),Y2(i,j)
        enddo
     enddo
  else
     do i=1,Ny1
        do j=1,Ny2
           read(719,*)X(j),Y1(i,j)
        enddo
     enddo
  endif
  close(719)
end subroutine sreadM_IR

subroutine sreadM_IC(pname,X,Y1,Y2)
  integer                                              :: i,j,Ny1,Ny2
  character(len=*)                                     :: pname
  complex(8),dimension(:,:)                            :: Y1
  complex(8),dimension(size(Y1,1),size(Y1,2)),optional :: Y2
  integer,dimension(size(Y1,2))                        :: X
  real(8),allocatable,dimension(:,:,:)     :: reY,imY
  Ny1=size(Y1,1) ; Ny2=size(Y1,2)
  open(719,file=adjustl(trim(pname)))
  if(present(Y2))then
     allocate(reY(2,size(Y1,1),size(Y1,2)) ,&
          imY(2,size(Y1,1),size(Y1,2)))
     do i=1,Ny1
        do j=1,Ny2
           read(719,"(I15,4(F21.12))")X(j),imY(1,i,j),reY(1,i,j),imY(2,i,j),reY(2,i,j)
        enddo
     enddo
     Y1=cmplx(reY(1,:,:),imY(1,:,:),8)
     Y2=cmplx(reY(2,:,:),imY(2,:,:),8)
  else
     allocate(reY(2,size(Y1,1),size(Y1,2)))
     do i=1,Ny1
        do j=1,Ny2
           read(719,*)X(j),imY(1,i,j),reY(1,i,j)
        enddo
     enddo
     Y1=cmplx(reY(1,:,:),imY(1,:,:),8)
  endif
  close(719)
end subroutine sreadM_IC


!----------------------------



subroutine sreadM_RI(pname,X,Y1,Y2,Y3,Y4)
  integer                                :: i,j,Ny1,Ny2
  character(len=*)                       :: pname
  integer,dimension(:,:)                   :: Y1
  integer,dimension(size(Y1,1),size(Y1,2)),optional :: Y2,Y3,Y4
  real(8),dimension(size(Y1,2))            :: X
  Ny1=size(Y1,1) ; Ny2=size(Y1,2)
  open(719,file=adjustl(trim(pname)))
  if(present(Y4))then
     do i=1,Ny1
        do j=1,Ny2
           write(719,"(F21.12,4(I15))")X(j),Y1(i,j),Y2(i,j),Y3(i,j),Y4(i,j)
        enddo
     enddo
  elseif(present(Y3))then
     do i=1,Ny1
        do j=1,Ny2
           write(719,"(F21.12,3(I15))")X(j),Y1(i,j),Y2(i,j),Y3(i,j)
        enddo
     enddo
  elseif(present(Y2))then
     do i=1,Ny1
        do j=1,Ny2
           read(719,*)X(j),Y1(i,j),Y2(i,j)
        enddo
     enddo
  else
     do i=1,Ny1
        do j=1,Ny2
           read(719,*)X(j),Y1(i,j)
        enddo
     enddo
  endif
  close(719)
end subroutine sreadM_RI

subroutine sreadM_RR(pname,X,Y1,Y2,Y3,Y4)
  integer                                :: i,j,Ny1,Ny2
  character(len=*)                       :: pname
  real(8),dimension(:,:)                   :: Y1
  real(8),dimension(size(Y1,1),size(Y1,2)),optional :: Y2,Y3,Y4
  real(8),dimension(size(Y1,2))            :: X
  Ny1=size(Y1,1) ; Ny2=size(Y1,2)
  open(719,file=adjustl(trim(pname)))
  if(present(Y4))then
     do i=1,Ny1
        do j=1,Ny2
           write(719,"(5(F21.12))")X(j),Y1(i,j),Y2(i,j),Y3(i,j),Y4(i,j)
        enddo
     enddo
  elseif(present(Y3))then
     do i=1,Ny1
        do j=1,Ny2
           write(719,"(4(F21.12))")X(j),Y1(i,j),Y2(i,j),Y3(i,j)
        enddo
     enddo
  elseif(present(Y2))then
     do i=1,Ny1
        do j=1,Ny2
           read(719,*)X(j),Y1(i,j),Y2(i,j)
        enddo
     enddo
  else
     do i=1,Ny1
        do j=1,Ny2
           read(719,*)X(j),Y1(i,j)
        enddo
     enddo
  endif
  close(719)
end subroutine sreadM_RR

subroutine sreadM_RC(pname,X,Y1,Y2)
  integer                                              :: i,j,Ny1,Ny2
  character(len=*)                                     :: pname
  complex(8),dimension(:,:)                            :: Y1
  complex(8),dimension(size(Y1,1),size(Y1,2)),optional :: Y2
  real(8),dimension(size(Y1,2))                        :: X
  real(8),allocatable,dimension(:,:,:)                 :: reY,imY
  Ny1=size(Y1,1) ; Ny2=size(Y1,2)
  open(719,file=adjustl(trim(pname)))
  if(present(Y2))then
     allocate(reY(2,size(Y1,1),size(Y1,2)) ,&
          imY(2,size(Y1,1),size(Y1,2)))
     do i=1,Ny1
        do j=1,Ny2
           read(719,"(F21.12,4(F21.12))")X(j),imY(1,i,j),reY(1,i,j),imY(2,i,j),reY(2,i,j)
        enddo
     enddo
     Y1=cmplx(reY(1,:,:),imY(1,:,:),8)
     Y2=cmplx(reY(2,:,:),imY(2,:,:),8)
  else
     allocate(reY(2,size(Y1,1),size(Y1,2)))
     do i=1,Ny1
        do j=1,Ny2
           read(719,*)X(j),imY(1,i,j),reY(1,i,j)
        enddo
     enddo
     Y1=cmplx(reY(1,:,:),imY(1,:,:),8)
  endif
  close(719)
end subroutine sreadM_RC


!----------------------------
!----------------------------
!----------------------------



subroutine sreadA3_II(pname,X,Y1,Y2)
  integer                                                      :: i,j,k,Ny1,Ny2,Ny3
  character(len=*)                                             :: pname
  integer,dimension(:,:,:)                                     :: Y1
  integer,dimension(size(Y1,1),size(Y1,2),size(Y1,3)),optional :: Y2
  integer,dimension(size(Y1,3))                                :: X

  Ny1=size(Y1,1) ; Ny2=size(Y1,2); Ny3=size(Y1,3)
  open(719,file=adjustl(trim(pname)))
  if(present(Y2))then
     do i=1,Ny1
        do j=1,Ny2
           do k=1,Ny3
              read(719,*)X(k),Y1(i,j,k),Y2(i,j,k)
           enddo

        enddo
     enddo
  else
     do i=1,Ny1
        do j=1,Ny2
           do k=1,Ny3
              read(719,*)X(k),Y1(i,j,k)
           enddo

        enddo
     enddo
  endif
  close(719)
end subroutine sreadA3_II

subroutine sreadA3_IR(pname,X,Y1,Y2)
  integer                                                      :: i,j,k,Ny1,Ny2,Ny3
  character(len=*)                                             :: pname
  real(8),dimension(:,:,:)                                     :: Y1
  real(8),dimension(size(Y1,1),size(Y1,2),size(Y1,3)),optional :: Y2
  integer,dimension(size(Y1,3))                                :: X

  Ny1=size(Y1,1) ; Ny2=size(Y1,2); Ny3=size(Y1,3)
  open(719,file=adjustl(trim(pname)))
  if(present(Y2))then
     do i=1,Ny1
        do j=1,Ny2
           do k=1,Ny3
              read(719,*)X(k),Y1(i,j,k),Y2(i,j,k)
           enddo

        enddo
     enddo
  else
     do i=1,Ny1
        do j=1,Ny2
           do k=1,Ny3
              read(719,*)X(k),Y1(i,j,k)
           enddo

        enddo
     enddo
  endif
  close(719)
end subroutine sreadA3_IR

subroutine sreadA3_IC(pname,X,Y1,Y2)
  integer                                                         :: i,j,k,Ny1,Ny2,Ny3
  character(len=*)                                                :: pname
  complex(8),dimension(:,:,:)                                     :: Y1
  complex(8),dimension(size(Y1,1),size(Y1,2),size(Y1,3)),optional :: Y2
  integer,dimension(size(Y1,3))                                   :: X
  real(8),allocatable,dimension(:,:,:,:)                          :: reY,imY
  Ny1=size(Y1,1) ; Ny2=size(Y1,2); Ny3=size(Y1,3)
  open(719,file=adjustl(trim(pname)))
  if(present(Y2))then
     allocate(reY(2,size(Y1,1),size(Y1,2),size(Y1,3)),&
          imY(2,size(Y1,1),size(Y1,2),size(Y1,3)))
     do i=1,Ny1
        do j=1,Ny2
           do k=1,Ny3
              read(719,"(I15,4(F21.12))")X(j),imY(1,i,j,k),reY(1,i,j,k),imY(2,i,j,k),reY(2,i,j,k)
           enddo
        enddo
     enddo
     Y1=cmplx(reY(1,:,:,:),imY(1,:,:,:),8)
     Y2=cmplx(reY(2,:,:,:),imY(2,:,:,:),8)
  else
     allocate(reY(2,size(Y1,1),size(Y1,2),size(Y1,3)))
     do i=1,Ny1
        do j=1,Ny2
           do k=1,Ny3
              read(719,*)X(j),imY(1,i,j,k),reY(1,i,j,k)
           enddo
        enddo
     enddo
     Y1=cmplx(reY(1,:,:,:),imY(1,:,:,:),8)
  endif
  close(719)
end subroutine sreadA3_IC

!----------------------------

subroutine sreadA3_RI(pname,X,Y1,Y2)
  integer                                                      :: i,j,k,Ny1,Ny2,Ny3
  character(len=*)                                             :: pname
  integer,dimension(:,:,:)                                     :: Y1
  integer,dimension(size(Y1,1),size(Y1,2),size(Y1,3)),optional :: Y2
  real(8),dimension(size(Y1,3))                                :: X
  Ny1=size(Y1,1) ; Ny2=size(Y1,2); Ny3=size(Y1,3)
  open(719,file=adjustl(trim(pname)))
  if(present(Y2))then
     do i=1,Ny1
        do j=1,Ny2
           do k=1,Ny3
              read(719,*)X(k),Y1(i,j,k),Y2(i,j,k)
           enddo
        enddo
     enddo
  else
     do i=1,Ny1
        do j=1,Ny2
           do k=1,Ny3
              read(719,*)X(k),Y1(i,j,k)
           enddo
        enddo
     enddo
  endif
  close(719)
end subroutine sreadA3_RI

subroutine sreadA3_RR(pname,X,Y1,Y2)
  integer                                                      :: i,j,k,Ny1,Ny2,Ny3
  character(len=*)                                             :: pname
  real(8),dimension(:,:,:)                                     :: Y1
  real(8),dimension(size(Y1,1),size(Y1,2),size(Y1,3)),optional :: Y2
  real(8),dimension(size(Y1,3))                                :: X
  Ny1=size(Y1,1) ; Ny2=size(Y1,2); Ny3=size(Y1,3)
  open(719,file=adjustl(trim(pname)))
  if(present(Y2))then
     do i=1,Ny1
        do j=1,Ny2
           do k=1,Ny3
              read(719,*)X(k),Y1(i,j,k),Y2(i,j,k)
           enddo
        enddo
     enddo
  else
     do i=1,Ny1
        do j=1,Ny2
           do k=1,Ny3
              read(719,*)X(k),Y1(i,j,k)
           enddo
        enddo
     enddo
  endif
  close(719)
end subroutine sreadA3_RR

subroutine sreadA3_RC(pname,X,Y1,Y2)
  integer                                                         :: i,j,k,Ny1,Ny2,Ny3
  character(len=*)                                                :: pname
  complex(8),dimension(:,:,:)                                     :: Y1
  complex(8),dimension(size(Y1,1),size(Y1,2),size(Y1,3)),optional :: Y2
  real(8),dimension(size(Y1,3))                                   :: X
  real(8),allocatable,dimension(:,:,:,:)                          :: reY,imY
  Ny1=size(Y1,1) ; Ny2=size(Y1,2); Ny3=size(Y1,3)
  open(719,file=adjustl(trim(pname)))
  if(present(Y2))then
     allocate(reY(2,size(Y1,1),size(Y1,2),size(Y1,3)),&
          imY(2,size(Y1,1),size(Y1,2),size(Y1,3)))
     do i=1,Ny1
        do j=1,Ny2
           do k=1,Ny3
              read(719,"(F21.12,4(F21.12))")X(j),imY(1,i,j,k),reY(1,i,j,k),imY(2,i,j,k),reY(2,i,j,k)
           enddo
        enddo
     enddo
     Y1=cmplx(reY(1,:,:,:),imY(1,:,:,:),8)
     Y2=cmplx(reY(2,:,:,:),imY(2,:,:,:),8)
  else
     allocate(reY(2,size(Y1,1),size(Y1,2),size(Y1,3)))
     do i=1,Ny1
        do j=1,Ny2
           do k=1,Ny3
              read(719,*)X(j),imY(1,i,j,k),reY(1,i,j,k)
           enddo
        enddo
     enddo
     Y1=cmplx(reY(1,:,:,:),imY(1,:,:,:),8)
  endif
  close(719)
end subroutine sreadA3_RC
!----------------------------
!----------------------------
!----------------------------
