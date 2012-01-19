!----------------------------

subroutine sreadV_II(pname,X,Y1,Y2)
  integer                       :: i,j,Np,Ny1,Ny2,Nx
  character(len=*)              :: pname
  integer,dimension(:)          :: X
  integer,dimension(:)          :: Y1
  integer,dimension(:),optional :: Y2
  Nx=size(X) ; Ny1=size(Y1);   if(present(Y2))Ny2=size(Y2)
  Np=min(Nx,Ny1);if(present(Y2))Np=min(Nx,Ny1,Ny2)
  if(Nx/=Ny1.OR.(present(Y2).AND.Nx/=Ny2))write(*,"(a,1x,I6,I6,I6)")"problem while printing "//trim(pname)//" Nx,Ny1,Ny2",Nx,Ny1,Ny2
  open(719,file=adjustl(trim(pname)))
  if(present(Y2))then
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

subroutine sreadV_IR(pname,X,Y1,Y2)
  integer                       :: i,j,Np,Ny1,Ny2,Nx
  character(len=*)              :: pname
  integer,dimension(:)          :: X
  real(8),dimension(:)          :: Y1
  real(8),dimension(:),optional :: Y2
  Nx=size(X) ; Ny1=size(Y1);   if(present(Y2))Ny2=size(Y2)
  Np=min(Nx,Ny1);if(present(Y2))Np=min(Nx,Ny1,Ny2)
  if(Nx/=Ny1.OR.(present(Y2).AND.Nx/=Ny2))&
       write(*,"(a,1x,I6,I6,I6)")"problem while printing "//trim(pname)//" Nx,Ny1,Ny2",Nx,Ny1,Ny2
  open(719,file=adjustl(trim(pname)))
  if(present(Y2))then
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

subroutine sreadV_IC(pname,X,Y1)
  integer                       :: i,j,Np,Ny1,Ny2,Nx
  character(len=*)              :: pname
  integer,dimension(:)          :: X
  complex(8),dimension(:)       :: Y1
  real(8),dimension(size(Y1))   :: reY1,imY1
  Nx=size(X) ; Ny1=size(Y1)
  if(Nx/=Ny1)write(*,"(a,1x,I6,I6,I6)")"problem while printing "//trim(pname)//" Nx,Ny1",Nx,Ny1
  Np=min(Nx,Ny1)
  open(719,file=adjustl(trim(pname)))
  do i=1,Np
     read(719,*)X(i),reY1(i),imY1(i)
  enddo
  Y1=cmplx(reY1,imY1,8)
  close(719)
end subroutine sreadV_IC

!----------------------------

subroutine sreadV_RI(pname,X,Y1,Y2)
  integer                       :: i,j,Np,Ny1,Ny2,Nx
  character(len=*)              :: pname
  real(8),dimension(:)          :: X
  integer,dimension(:)          :: Y1
  integer,dimension(:),optional :: Y2
  Nx=size(X) ; Ny1=size(Y1);   if(present(Y2))Ny2=size(Y2)
  Np=min(Nx,Ny1);if(present(Y2))Np=min(Nx,Ny1,Ny2)
  if(Nx/=Ny1.OR.(present(Y2).AND.Nx/=Ny2))write(*,"(a,1x,I6,I6,I6)")"problem while printing "//trim(pname)//" Nx,Ny1,Ny2",Nx,Ny1,Ny2
  open(719,file=adjustl(trim(pname)))
  if(present(Y2))then
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

subroutine sreadV_RR(pname,X,Y1,Y2)
  integer                       :: i,j,Np,Ny1,Ny2,Nx
  character(len=*)              :: pname
  real(8),dimension(:)          :: X
  real(8),dimension(:)          :: Y1
  real(8),dimension(:),optional :: Y2
  Nx=size(X) ; Ny1=size(Y1);   if(present(Y2))Ny2=size(Y2)
  Np=min(Nx,Ny1);if(present(Y2))Np=min(Nx,Ny1,Ny2)
  if(Nx/=Ny1.OR.(present(Y2).AND.Nx/=Ny2))write(*,"(a,1x,I6,I6,I6)")"problem while printing "//trim(pname)//" Nx,Ny1,Ny2",Nx,Ny1,Ny2
  open(719,file=adjustl(trim(pname)))
  if(present(Y2))then
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

subroutine sreadV_RC(pname,X,Y1)
  integer                       :: i,j,Np,Ny1,Ny2,Nx
  character(len=*)              :: pname
  real(8),dimension(:)          :: X
  complex(8),dimension(:)       :: Y1
  real(8),dimension(size(Y1))   :: reY1,imY1
  Nx=size(X) ; Ny1=size(Y1)
  Np=min(Nx,Ny1)
  if(Nx/=Ny1)write(*,"(a,1x,I6,I6,I6)")"problem while printing "//trim(pname)//" Nx,Ny",Nx,Ny1
  open(719,file=trim(adjustl(trim(pname))))
  do i=1,Np
     read(719,*)X(i),imY1(i),reY1(i)
  enddo
  Y1=cmplx(reY1,imY1,8)
  close(719)
end subroutine sreadV_RC

!----------------------------


