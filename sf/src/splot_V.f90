subroutine splotV_II(pname,X,Y1,Y2,Y3,Y4,Y5,Y6,Y7,Y8,append)
  integer                             :: i,Np!,Ny1,Ny2,Nx
  character(len=*)                    :: pname
  integer,dimension(:)                :: X,Y1
  integer,dimension(size(X)),optional :: Y2,Y3,Y4,Y5,Y6,Y7,Y8
  logical,optional                    :: append
  logical                             :: check
  Np=size(X)!Nx
  !here as a reminder of the older version:
  ! Ny1=size(Y1)
  ! if(present(Y2))Ny2=size(Y2)
  ! Np=min(Nx,Ny1);if(present(Y2))Np=min(Nx,Ny1,Ny2)
  ! if(Nx/=Ny1.OR.(present(Y2).AND.Nx/=Ny2))&
  !      write(*,"(a,1x,I6,I6,I6)")"problem while printing "//trim(pname)//" Nx,Ny1,Ny2",Nx,Ny1,Ny2
  if(present(append).AND. append.eqv..true.)then
     inquire(file=trim(adjustl(trim(pname))),exist=check)
     open(719,file=adjustl(trim(pname)),access="append")
     if(check)write(719,*)
  else
     open(719,file=adjustl(trim(pname)))
  endif
  if(present(Y8))then
     do i=1,Np
        write(719,"(9(I15))")X(i),Y1(i),Y2(i),Y3(i),Y4(i),Y5(i),Y6(i),Y7(i),Y8(i)
     enddo
  elseif(present(Y7))then
     do i=1,Np
        write(719,"(8(I15))")X(i),Y1(i),Y2(i),Y3(i),Y4(i),Y5(i),Y6(i),Y7(i)
     enddo
  elseif(present(Y6))then
     do i=1,Np
        write(719,"(7(I15))")X(i),Y1(i),Y2(i),Y3(i),Y4(i),Y5(i),Y6(i)
     enddo
  elseif(present(Y5))then
     do i=1,Np
        write(719,"(6(I15))")X(i),Y1(i),Y2(i),Y3(i),Y4(i),Y5(i)
     enddo
  elseif(present(Y4))then
     do i=1,Np
        write(719,"(5(I15))")X(i),Y1(i),Y2(i),Y3(i),Y4(i)
     enddo
  elseif(present(Y3))then
     do i=1,Np
        write(719,"(4(I15))")X(i),Y1(i),Y2(i),Y3(i)
     enddo
  elseif(present(Y2))then
     do i=1,Np
        write(719,*)X(i),Y1(i),Y2(i)
     enddo
  else
     do i=1,Np
        write(719,*)X(i),Y1(i)
     enddo
  endif
  close(719)
end subroutine splotV_II

!----------------------------

subroutine splotV_IR(pname,X,Y1,Y2,Y3,Y4,Y5,Y6,Y7,Y8,append)
  integer                             :: i,Np
  character(len=*)                    :: pname
  integer,dimension(:)                :: X
  real(8),dimension(size(X))          :: Y1
  real(8),dimension(size(X)),optional :: Y2,Y3,Y4,Y5,Y6,Y7,Y8
  logical,optional                    :: append
  logical                             :: check
  Np=size(X)
  if(present(append).AND. append.eqv..true.)then
     inquire(file=trim(adjustl(trim(pname))),exist=check)
     open(719,file=adjustl(trim(pname)),access="append")
     if(check)write(719,*)
  else
     open(719,file=adjustl(trim(pname)))
  endif
  if(present(Y8))then
     do i=1,Np
        write(719,"(I15,8(F18.10))")X(i),Y1(i),Y2(i),Y3(i),Y4(i),Y5(i),Y6(i),Y7(i),Y8(i)
     enddo
  elseif(present(Y7))then
     do i=1,Np
        write(719,"(I15,7(F18.10))")X(i),Y1(i),Y2(i),Y3(i),Y4(i),Y5(i),Y6(i),Y7(i)
     enddo
  elseif(present(Y6))then
     do i=1,Np
        write(719,"(I15,6(F18.10))")X(i),Y1(i),Y2(i),Y3(i),Y4(i),Y5(i),Y6(i)
     enddo
  elseif(present(Y5))then
     do i=1,Np
        write(719,"(I15,5(F18.10))")X(i),Y1(i),Y2(i),Y3(i),Y4(i),Y5(i)
     enddo
  elseif(present(Y4))then
     do i=1,Np
        write(719,"(I15,4(F18.10))")X(i),Y1(i),Y2(i),Y3(i),Y4(i)
     enddo
  elseif(present(Y3))then
     do i=1,Np
        write(719,"(I15,3(F18.10))")X(i),Y1(i),Y2(i),Y3(i)
     enddo
  elseif(present(Y2))then
     do i=1,Np
        write(719,*)X(i),Y1(i),Y2(i)
     enddo
  else
     do i=1,Np
        write(719,*)X(i),Y1(i)
     enddo
  endif
  close(719)
end subroutine splotV_IR

!----------------------------

subroutine splotV_IC(pname,X,Y1,Y2,Y3,Y4,append)
  integer                                :: i,Np
  character(len=*)                       :: pname
  integer,dimension(:)                   :: X
  complex(8),dimension(:)                :: Y1
  complex(8),dimension(size(X)),optional :: Y2,Y3,Y4
  logical,optional                       :: append
  logical                                :: check
  Np=size(X)
  if(present(append).AND. append.eqv..true.)then
     inquire(file=trim(adjustl(trim(pname))),exist=check)
     open(719,file=adjustl(trim(pname)),access="append")
     if(check)write(719,*)
  else
     open(719,file=adjustl(trim(pname)))
  endif
  if(present(Y4))then
     do i=1,Np
        write(719,"(I15,8(F18.10))")X(i),real(Y1(i),8),dimag(Y1(i)),real(Y2(i),8),dimag(Y2(i)),real(Y3(i),8),dimag(Y3(i)),&
             real(Y4(i),8),dimag(Y4(i))
     enddo
  elseif(present(Y3))then
     do i=1,Np
        write(719,"(I15,6(F18.10))")X(i),real(Y1(i),8),dimag(Y1(i)),real(Y2(i),8),dimag(Y2(i)),real(Y3(i),8),dimag(Y3(i))
     enddo
  elseif(present(Y2))then
     do i=1,Np
        write(719,"(I15,4(F18.10))")X(i),real(Y1(i),8),dimag(Y1(i)),real(Y2(i),8),dimag(Y2(i))
     enddo
  else
     do i=1,Np
        write(719,*)X(i),real(Y1(i),8),dimag(Y1(i))
     enddo
  endif
  close(719)
end subroutine splotV_IC

!----------------------------
!----------------------------
!----------------------------

subroutine splotV_RI(pname,X,Y1,Y2,Y3,Y4,Y5,Y6,Y7,Y8,append)
  integer                             :: i,Np
  character(len=*)                    :: pname
  real(8),dimension(:)                :: X
  integer,dimension(size(X))          :: Y1
  integer,dimension(size(X)),optional :: Y2,Y3,Y4,Y5,Y6,Y7,Y8
  logical,optional                    :: append
  logical                             :: check
  Np=size(X)
  if(present(append).AND. append.eqv..true.)then
     inquire(file=trim(adjustl(trim(pname))),exist=check)
     open(719,file=adjustl(trim(pname)),access="append")
     if(check)write(719,*)
  else
     open(719,file=adjustl(trim(pname)))
  endif
  if(present(Y8))then
     do i=1,Np
        write(719,"(F18.10,8(I15))")X(i),Y1(i),Y2(i),Y3(i),Y4(i),Y5(i),Y6(i),Y7(i),Y8(i)
     enddo
  elseif(present(Y7))then
     do i=1,Np
        write(719,"(F18.10,7(I15))")X(i),Y1(i),Y2(i),Y3(i),Y4(i),Y5(i),Y6(i),Y7(i)
     enddo
  elseif(present(Y6))then
     do i=1,Np
        write(719,"(F18.10,6(I15))")X(i),Y1(i),Y2(i),Y3(i),Y4(i),Y5(i),Y6(i)
     enddo
  elseif(present(Y5))then
     do i=1,Np
        write(719,"(F18.10,5(I15))")X(i),Y1(i),Y2(i),Y3(i),Y4(i),Y5(i)
     enddo
  elseif(present(Y4))then
     do i=1,Np
        write(719,"(F18.10,4(I15))")X(i),Y1(i),Y2(i),Y3(i),Y4(i)
     enddo
  elseif(present(Y3))then
     do i=1,Np
        write(719,"(F18.10,3(I15))")X(i),Y1(i),Y2(i),Y3(i)
     enddo
  elseif(present(Y2))then
     do i=1,Np
        write(719,*)X(i),Y1(i),Y2(i)
     enddo
  else
     do i=1,Np
        write(719,*)X(i),Y1(i)
     enddo
  endif
  close(719)
end subroutine splotV_RI

!----------------------------

subroutine splotV_RR(pname,X,Y1,Y2,Y3,Y4,Y5,Y6,Y7,Y8,append)
  integer                       :: i,Np
  character(len=*)              :: pname
  real(8),dimension(:)          :: X
  real(8),dimension(:)          :: Y1
  real(8),dimension(:),optional :: Y2,Y3,Y4,Y5,Y6,Y7,Y8
  logical,optional              :: append
  logical                       :: check
  Np=size(X)
  if(present(append).AND. append.eqv..true.)then
     inquire(file=trim(adjustl(trim(pname))),exist=check)
     open(719,file=adjustl(trim(pname)),access="append")
     if(check)write(719,*)
  else
     open(719,file=adjustl(trim(pname)))
  endif
  if(present(Y8))then
     do i=1,Np
        write(719,"(9(F18.10))")X(i),Y1(i),Y2(i),Y3(i),Y4(i),Y5(i),Y6(i),Y7(i),Y8(i)
     enddo
  elseif(present(Y7))then
     do i=1,Np
        write(719,"(8(F18.10))")X(i),Y1(i),Y2(i),Y3(i),Y4(i),Y5(i),Y6(i),Y7(i)
     enddo
  elseif(present(Y6))then
     do i=1,Np
        write(719,"(7(F18.10))")X(i),Y1(i),Y2(i),Y3(i),Y4(i),Y5(i),Y6(i)
     enddo
  elseif(present(Y5))then
     do i=1,Np
        write(719,"(6(F18.10))")X(i),Y1(i),Y2(i),Y3(i),Y4(i),Y5(i)
     enddo
  elseif(present(Y4))then
     do i=1,Np
        write(719,"(5(F18.10))")X(i),Y1(i),Y2(i),Y3(i),Y4(i)
     enddo
  elseif(present(Y3))then
     do i=1,Np
        write(719,"(4(F18.10))")X(i),Y1(i),Y2(i),Y3(i)
     enddo
  elseif(present(Y2))then
     do i=1,Np
        write(719,*)X(i),Y1(i),Y2(i)
     enddo
  else
     do i=1,Np
        write(719,*)X(i),Y1(i)
     enddo
  endif
  close(719)
end subroutine splotV_RR

!----------------------------

subroutine splotV_RC(pname,X,Y1,Y2,Y3,Y4,append)
  integer                          :: i,Np
  character(len=*)                 :: pname
  real(8),dimension(:)             :: X
  complex(8),dimension(:)          :: Y1
  complex(8),dimension(:),optional :: Y2,Y3,Y4
  logical,optional                 :: append
  logical                          :: check
  Np=size(X)
  if(present(append).AND. append.eqv..true.)then
     inquire(file=trim(adjustl(trim(pname))),exist=check)
     open(719,file=adjustl(trim(pname)),access="append")
     if(check)write(719,*)
  else
     open(719,file=adjustl(trim(pname)))
  endif
  if(present(Y4))then
     do i=1,Np
        write(719,"(F18.10,8(F18.10))")X(i),aimag(Y1(i)),real(Y1(i)),aimag(Y2(i)),real(Y2(i)),aimag(Y3(i)),real(Y3(i)),&
             aimag(Y4(i)),real(Y4(i))
     enddo
  elseif(present(Y3))then
     do i=1,Np
        write(719,"(F18.10,6(F18.10))")X(i),aimag(Y1(i)),real(Y1(i)),aimag(Y2(i)),real(Y2(i)),aimag(Y3(i)),real(Y3(i))
     enddo
  elseif(present(Y2))then
     do i=1,Np
        write(719,"(F18.10,4(F18.10))")X(i),aimag(Y1(i)),real(Y1(i)),aimag(Y2(i)),real(Y2(i))
     enddo
  else
     do i=1,Np
        write(719,*)X(i),aimag(Y1(i)),real(Y1(i))
     enddo
  endif
  close(719)
end subroutine splotV_RC

!----------------------------


