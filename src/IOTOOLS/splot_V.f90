!----------------------------

subroutine splotV_II(pname,X,Y1,Y2,append)
  integer                       :: i,j,Np,Ny1,Ny2,Nx
  character(len=*)              :: pname
  integer,dimension(:)          :: X
  integer,dimension(:)          :: Y1
  integer,dimension(:),optional :: Y2
  logical,optional              :: append
  logical                       :: check
  Nx=size(X) ; Ny1=size(Y1);   if(present(Y2))Ny2=size(Y2)
  Np=min(Nx,Ny1);if(present(Y2))Np=min(Nx,Ny1,Ny2)
  if(Nx/=Ny1.OR.(present(Y2).AND.Nx/=Ny2))write(*,"(a,1x,I6,I6,I6)")"problem while printing "//trim(pname)//" Nx,Ny1,Ny2",Nx,Ny1,Ny2

  if(present(append).AND. append==.true.)then
     inquire(file=trim(adjustl(trim(pname))),exist=check)
     open(719,file=adjustl(trim(pname)),access="append")
     if(check)write(719,*)
  else
     open(719,file=adjustl(trim(pname)))
  endif
  if(present(Y2))then
     do i=1,Np
        write(719,*)X(i),Y1(i),Y2(i)
     enddo
  else
     do i=1,Np
        write(719,*)X(i),Y1(i)
     enddo
  endif
  !if(present(append))  write(719,*)
  close(719)
end subroutine splotV_II

!----------------------------

subroutine splotV_IR(pname,X,Y1,Y2,append)
  integer                       :: i,j,Np,Ny1,Ny2,Nx
  character(len=*)              :: pname
  integer,dimension(:)          :: X
  real(8),dimension(:)          :: Y1
  real(8),dimension(:),optional :: Y2
  logical,optional              :: append
  logical                       :: check
  Nx=size(X) ; Ny1=size(Y1);   if(present(Y2))Ny2=size(Y2)
  Np=min(Nx,Ny1);if(present(Y2))Np=min(Nx,Ny1,Ny2)
  if(Nx/=Ny1.OR.(present(Y2).AND.Nx/=Ny2))&
       write(*,"(a,1x,I6,I6,I6)")"problem while printing "//trim(pname)//" Nx,Ny1,Ny2",Nx,Ny1,Ny2
  if(present(append).AND. append==.true.)then
     inquire(file=trim(adjustl(trim(pname))),exist=check)
     open(719,file=adjustl(trim(pname)),access="append")
     if(check)write(719,*)
  else
     open(719,file=adjustl(trim(pname)))
  endif
  if(present(Y2))then
     do i=1,Np
        write(719,*)X(i),Y1(i),Y2(i)
     enddo
  else
     do i=1,Np
        write(719,*)X(i),Y1(i)
     enddo
  endif
  !if(present(append))  write(719,*)
  close(719)
end subroutine splotV_IR

!----------------------------

subroutine splotV_IC(pname,X,Y1,append)
  integer                       :: i,j,Np,Ny1,Ny2,Nx
  character(len=*)              :: pname
  integer,dimension(:)          :: X
  complex(8),dimension(:)          :: Y1
logical,optional              :: append
  logical                       :: check
  Nx=size(X) ; Ny1=size(Y1) 
  Np=min(Nx,Ny1)
  if(Nx/=Ny1)write(*,"(a,1x,I6,I6,I6)")"problem while printing "//trim(pname)//" Nx,Ny1",Nx,Ny1
  if(present(append).AND. append==.true.)then
     inquire(file=trim(adjustl(trim(pname))),exist=check)
     open(719,file=adjustl(trim(pname)),access="append")
     if(check)write(719,*)
  else
     open(719,file=adjustl(trim(pname)))
  endif
  do i=1,Np
     write(719,*)X(i),real(Y1(i)),aimag(Y1(i))
  enddo
  !if(present(append))  write(719,*)
  close(719)
end subroutine splotV_IC

!----------------------------

subroutine splotV_RI(pname,X,Y1,Y2,append)
  integer                       :: i,j,Np,Ny1,Ny2,Nx
  character(len=*)              :: pname
  real(8),dimension(:)          :: X
  integer,dimension(:)          :: Y1
  integer,dimension(:),optional :: Y2
logical,optional              :: append
  logical                       :: check
  Nx=size(X) ; Ny1=size(Y1);   if(present(Y2))Ny2=size(Y2)
  Np=min(Nx,Ny1);if(present(Y2))Np=min(Nx,Ny1,Ny2)
  if(Nx/=Ny1.OR.(present(Y2).AND.Nx/=Ny2))write(*,"(a,1x,I6,I6,I6)")"problem while printing "//trim(pname)//" Nx,Ny1,Ny2",Nx,Ny1,Ny2
  if(present(append).AND. append==.true.)then
     inquire(file=trim(adjustl(trim(pname))),exist=check)
     open(719,file=adjustl(trim(pname)),access="append")
     if(check)write(719,*)
  else
     open(719,file=adjustl(trim(pname)))
  endif
  if(present(Y2))then
     do i=1,Np
        write(719,*)X(i),Y1(i),Y2(i)
     enddo
  else
     do i=1,Np
        write(719,*)X(i),Y1(i)
     enddo
  endif
  !if(present(append))  write(719,*)
  close(719)
end subroutine splotV_RI

!----------------------------

subroutine splotV_RR(pname,X,Y1,Y2,append)
  integer                       :: i,j,Np,Ny1,Ny2,Nx
  character(len=*)              :: pname
  real(8),dimension(:)          :: X
  real(8),dimension(:)          :: Y1
  real(8),dimension(:),optional :: Y2
logical,optional              :: append
  logical                       :: check
  Nx=size(X) ; Ny1=size(Y1);   if(present(Y2))Ny2=size(Y2)
  Np=min(Nx,Ny1);if(present(Y2))Np=min(Nx,Ny1,Ny2)
  if(Nx/=Ny1.OR.(present(Y2).AND.Nx/=Ny2))write(*,"(a,1x,I6,I6,I6)")"problem while printing "//trim(pname)//" Nx,Ny1,Ny2",Nx,Ny1,Ny2
  if(present(append).AND. append==.true.)then
     inquire(file=trim(adjustl(trim(pname))),exist=check)
     open(719,file=adjustl(trim(pname)),access="append")
     if(check)write(719,*)
  else
     open(719,file=adjustl(trim(pname)))
  endif
  if(present(Y2))then
     do i=1,Np
        write(719,*)X(i),Y1(i),Y2(i)
     enddo
  else
     do i=1,Np
        write(719,*)X(i),Y1(i)
     enddo
  endif
  !if(present(append))  write(719,*)
  close(719)
end subroutine splotV_RR

!----------------------------

subroutine splotV_RC(pname,X,Y1,append)
  integer                       :: i,j,Np,Ny1,Ny2,Nx
  character(len=*)              :: pname
  real(8),dimension(:)          :: X
  complex(8),dimension(:)       :: Y1
  logical,optional              :: append
  logical                       :: check
  Nx=size(X) ; Ny1=size(Y1)
  Np=min(Nx,Ny1)
  if(Nx/=Ny1)write(*,"(a,1x,I6,I6,I6)")"problem while printing "//trim(pname)//" Nx,Ny",Nx,Ny1
  if(present(append).AND. append==.true.)then
     inquire(file=trim(adjustl(trim(pname))),exist=check)
     open(719,file=adjustl(trim(pname)),access="append")
     if(check)write(719,*)
  else
     open(719,file=adjustl(trim(pname)))
  endif
  do i=1,Np
     write(719,*)X(i),aimag(Y1(i)),real(Y1(i))
  enddo
  !if(present(append))  write(719,*)
  close(719)
end subroutine splotV_RC

!----------------------------


