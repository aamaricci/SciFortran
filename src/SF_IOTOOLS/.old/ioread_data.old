subroutine data_readV_I(pname,Y1)
  integer               :: i,Np
  character(len=*)      :: pname
  integer,dimension(:)  :: Y1
  Np=size(Y1)  
  call data_open(trim(pname))
  include "ioread_control.f90"
  open(719,file=str(pname))
  do i=1,Np
     read(719,*)Y1(i)
  enddo
  close(719)
  call data_store(pname)
end subroutine data_readV_I
!------------------------------------------------------------------!
subroutine data_readV_R(pname,Y1)
  integer               :: i,Np
  character(len=*)      :: pname
  real(8),dimension(:)  :: Y1
  Np=size(Y1)
  call data_open(trim(pname))
  include "ioread_control.f90"
  open(719,file=str(pname))
  do i=1,Np
     read(719,*)Y1(i)
  enddo
  close(719)
  call data_store(pname)
end subroutine data_readV_R
!------------------------------------------------------------------!
subroutine data_readV_C(pname,Y1)
  integer                :: i,Np
  character(len=*)       :: pname
  complex(8),dimension(:):: Y1
  real(8)                :: reY1,imY1
  Np=size(Y1)  
  call data_open(trim(pname))
  include "ioread_control.f90"
  open(719,file=str(pname))
  do i=1,Np
     read(719,*)imY1,reY1
     Y1(i)=cmplx(reY1,imY1,8)
  enddo
  close(719)
  call data_store(pname)
end subroutine data_readV_C

!------------------------------------------------------------------!
!------------------------------------------------------------------!
!------------------------------------------------------------------!

subroutine data_readM_I(pname,Y1,X)
  integer                       :: i,j,Ny1,Ny2
  character(len=*)              :: pname
  integer,dimension(:,:)        :: Y1
  real(8),optional,dimension(:) :: X(size(Y1,2))
  Ny1=size(Y1,1) ; Ny2=size(Y1,2)
  call data_open(trim(pname))
  include "ioread_control.f90"
  open(719,file=str(pname))
  do i=1,Ny1
     do j=1,Ny2
        if(present(X))then
           read(719,*)X(j),Y1(i,j)
        else
           read(719,*)Y1(i,j)
        endif
     enddo
  enddo
  close(719)
  call data_store(pname)
end subroutine data_readM_I
!------------------------------------------------------------------!
subroutine data_readM_R(pname,Y1,X)
  integer                       :: i,j,Ny1,Ny2
  character(len=*)              :: pname
  real(8),dimension(:,:)        :: Y1
  real(8),optional,dimension(:) :: X(size(Y1,2))
  Ny1=size(Y1,1) ; Ny2=size(Y1,2)
  call data_open(trim(pname))
  include "ioread_control.f90"
  open(719,file=str(pname))
  do i=1,Ny1
     do j=1,Ny2
        if(present(X))then
           read(719,*)X(j),Y1(i,j)
        else
           read(719,*)Y1(i,j)
        endif
     enddo
  enddo
  close(719)
  call data_store(pname)
end subroutine data_readM_R
!------------------------------------------------------------------!
subroutine data_readM_C(pname,Y1,X)
  integer                       :: i,j,Ny1,Ny2
  character(len=*)              :: pname
  complex(8),dimension(:,:)     :: Y1
  real(8),optional,dimension(:) :: X(size(Y1,2))
  real(8)                       :: reY,imY
  Ny1=size(Y1,1) ; Ny2=size(Y1,2)
  call data_open(trim(pname))
  include "ioread_control.f90"
  open(719,file=str(pname))
  do i=1,Ny1
     do j=1,Ny2
        if(present(X))then
           read(719,*)X(j),imY,reY!Y1(i,j)
        else
           read(719,*)imY,reY!Y1(i,j)
        endif
        Y1(i,j)=dcmplx(reY,imY)
     enddo
  enddo
  close(719)
  call data_store(pname)
end subroutine data_readM_C

!------------------------------------------------------------------!
!------------------------------------------------------------------!
!------------------------------------------------------------------!


subroutine data_readA3_I(pname,Y1,X)
  integer                       :: i,j,k,Ny1,Ny2,Ny3
  character(len=*)              :: pname
  integer,dimension(:,:,:)      :: Y1
  real(8),optional,dimension(:) :: X(size(Y1,3))
  Ny1=size(Y1,1) ; Ny2=size(Y1,2) ; Ny3=size(Y1,3) 
  call data_open(trim(pname))
  include "ioread_control.f90"
  open(719,file=str(pname))
  do i=1,Ny1
     do j=1,Ny2
        do k=1,Ny3
           if(present(X))then
              read(719,*)X(k),Y1(i,j,k)
           else
              read(719,*)Y1(i,j,k)
           endif
        enddo
     enddo
  enddo
  close(719)
  call data_store(pname)
end subroutine data_readA3_I
!------------------------------------------------------------------!
subroutine data_readA3_R(pname,Y1,X)
  integer                       :: i,j,k,Ny1,Ny2,Ny3
  character(len=*)              :: pname
  real(8),dimension(:,:,:)      :: Y1
  real(8),optional,dimension(:) :: X(size(Y1,3))
  Ny1=size(Y1,1) ; Ny2=size(Y1,2) ; Ny3=size(Y1,3) 
  call data_open(trim(pname))
  include "ioread_control.f90"
  open(719,file=str(pname))
  do i=1,Ny1
     do j=1,Ny2
        do k=1,Ny3
           if(present(X))then
              read(719,*)X(k),Y1(i,j,k)
           else
              read(719,*)Y1(i,j,k)
           endif
        enddo
     enddo
  enddo
  close(719)
  call data_store(pname)
end subroutine data_readA3_R
!------------------------------------------------------------------!
subroutine data_readA3_C(pname,Y1,X)
  integer                       :: i,j,k,Ny1,Ny2,Ny3
  character(len=*)              :: pname
  complex(8),dimension(:,:,:)   :: Y1
  real(8),optional,dimension(:) :: X(size(Y1,3))
  real(8)                       :: reY,imY
  Ny1=size(Y1,1) ; Ny2=size(Y1,2) ; Ny3=size(Y1,3) 
  call data_open(trim(pname))
  include "ioread_control.f90"

  open(719,file=str(pname))
  do i=1,Ny1
     do j=1,Ny2
        do k=1,Ny3
           if(present(X))then
              read(719,*)X(k),imY,reY
           else
              read(719,*)imY,reY
           endif
           Y1(i,j,k)=dcmplx(reY,imY)
        enddo
     enddo
  enddo
  close(719)
  call data_store(pname)
end subroutine data_readA3_C

!------------------------------------------------------------------!
!------------------------------------------------------------------!
!------------------------------------------------------------------!
