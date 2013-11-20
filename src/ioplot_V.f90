subroutine splotV_II(pname,X,Y1,Y2,Y3,Y4,Y5,Y6,Y7,Y8,append)
  integer                             :: i,Np
  character(len=*)                    :: pname
  integer,dimension(:)                :: X
  integer,dimension(size(X))          :: Y1
  integer,dimension(size(X)),optional :: Y2,Y3,Y4,Y5,Y6,Y7,Y8
  logical,optional                    :: append
  logical                             :: check,rw
  rw=.false.;if(present(append))rw=append
  Np=size(X)
  if(rw)then
     inquire(file=trim(adjustl(trim(pname))),exist=check)
     open(719,file=adjustl(trim(pname)),position="append")
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
  logical                             :: check,rw
  rw=.false.;if(present(append))rw=append
  Np=size(X)
  if(rw)then
     inquire(file=trim(adjustl(trim(pname))),exist=check)
     open(719,file=adjustl(trim(pname)),position="append")
     if(check)write(719,*)
  else
     open(719,file=adjustl(trim(pname)))
  endif
  if(present(Y8))then
     do i=1,Np
        write(719,"(I15,8(F21.12))")X(i),Y1(i),Y2(i),Y3(i),Y4(i),Y5(i),Y6(i),Y7(i),Y8(i)
     enddo
  elseif(present(Y7))then
     do i=1,Np
        write(719,"(I15,7(F21.12))")X(i),Y1(i),Y2(i),Y3(i),Y4(i),Y5(i),Y6(i),Y7(i)
     enddo
  elseif(present(Y6))then
     do i=1,Np
        write(719,"(I15,6(F21.12))")X(i),Y1(i),Y2(i),Y3(i),Y4(i),Y5(i),Y6(i)
     enddo
  elseif(present(Y5))then
     do i=1,Np
        write(719,"(I15,5(F21.12))")X(i),Y1(i),Y2(i),Y3(i),Y4(i),Y5(i)
     enddo
  elseif(present(Y4))then
     do i=1,Np
        write(719,"(I15,4(F21.12))")X(i),Y1(i),Y2(i),Y3(i),Y4(i)
     enddo
  elseif(present(Y3))then
     do i=1,Np
        write(719,"(I15,3(F21.12))")X(i),Y1(i),Y2(i),Y3(i)
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
  complex(8),dimension(size(X))          :: Y1
  complex(8),dimension(size(X)),optional :: Y2,Y3,Y4
  logical,optional                       :: append
  logical                             :: check,rw
  rw=.false.;if(present(append))rw=append
  Np=size(X)
  if(rw)then
     inquire(file=trim(adjustl(trim(pname))),exist=check)
     open(719,file=adjustl(trim(pname)),position="append")
     if(check)write(719,*)
  else
     open(719,file=adjustl(trim(pname)))
  endif
  if(present(Y4))then
     do i=1,Np
        write(719,"(I15,8(F21.12))")X(i),dreal(Y1(i)),dimag(Y1(i)),&
             dreal(Y2(i)),dimag(Y2(i)),dreal(Y3(i)),dimag(Y3(i)),&
             dreal(Y4(i)),dimag(Y4(i))
     enddo
  elseif(present(Y3))then
     do i=1,Np
        write(719,"(I15,6(F21.12))")X(i),dreal(Y1(i)),dimag(Y1(i)),&
             dreal(Y2(i)),dimag(Y2(i)),dreal(Y3(i)),dimag(Y3(i))
     enddo
  elseif(present(Y2))then
     do i=1,Np
        write(719,"(I15,4(F21.12))")X(i),dreal(Y1(i)),dimag(Y1(i)),&
             dreal(Y2(i)),dimag(Y2(i))
     enddo
  else
     do i=1,Np
        write(719,*)X(i),dreal(Y1(i)),dimag(Y1(i))
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
  logical                             :: check,rw
  rw=.false.;if(present(append))rw=append
  Np=size(X)
  if(rw)then
     inquire(file=trim(adjustl(trim(pname))),exist=check)
     open(719,file=adjustl(trim(pname)),position="append")
     if(check)write(719,*)
  else
     open(719,file=adjustl(trim(pname)))
  endif
  if(present(Y8))then
     do i=1,Np
        write(719,"(F21.12,8(I15))")X(i),Y1(i),Y2(i),Y3(i),Y4(i),Y5(i),Y6(i),Y7(i),Y8(i)
     enddo
  elseif(present(Y7))then
     do i=1,Np
        write(719,"(F21.12,7(I15))")X(i),Y1(i),Y2(i),Y3(i),Y4(i),Y5(i),Y6(i),Y7(i)
     enddo
  elseif(present(Y6))then
     do i=1,Np
        write(719,"(F21.12,6(I15))")X(i),Y1(i),Y2(i),Y3(i),Y4(i),Y5(i),Y6(i)
     enddo
  elseif(present(Y5))then
     do i=1,Np
        write(719,"(F21.12,5(I15))")X(i),Y1(i),Y2(i),Y3(i),Y4(i),Y5(i)
     enddo
  elseif(present(Y4))then
     do i=1,Np
        write(719,"(F21.12,4(I15))")X(i),Y1(i),Y2(i),Y3(i),Y4(i)
     enddo
  elseif(present(Y3))then
     do i=1,Np
        write(719,"(F21.12,3(I15))")X(i),Y1(i),Y2(i),Y3(i)
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
  integer                             :: i,Np
  character(len=*)                    :: pname
  real(8),dimension(:)                :: X
  real(8),dimension(size(X))          :: Y1
  real(8),dimension(size(X)),optional :: Y2,Y3,Y4,Y5,Y6,Y7,Y8
  logical,optional                    :: append
  logical                             :: check,rw
  rw=.false.;if(present(append))rw=append
  Np=size(X)
  if(rw)then
     inquire(file=trim(adjustl(trim(pname))),exist=check)
     open(719,file=adjustl(trim(pname)),position="append")
     if(check)write(719,*)
  else
     open(719,file=adjustl(trim(pname)))
  endif
  if(present(Y8))then
     do i=1,Np
        write(719,"(9(F21.12))")X(i),Y1(i),Y2(i),Y3(i),Y4(i),Y5(i),Y6(i),Y7(i),Y8(i)
     enddo
  elseif(present(Y7))then
     do i=1,Np
        write(719,"(8(F21.12))")X(i),Y1(i),Y2(i),Y3(i),Y4(i),Y5(i),Y6(i),Y7(i)
     enddo
  elseif(present(Y6))then
     do i=1,Np
        write(719,"(7(F21.12))")X(i),Y1(i),Y2(i),Y3(i),Y4(i),Y5(i),Y6(i)
     enddo
  elseif(present(Y5))then
     do i=1,Np
        write(719,"(6(F21.12))")X(i),Y1(i),Y2(i),Y3(i),Y4(i),Y5(i)
     enddo
  elseif(present(Y4))then
     do i=1,Np
        write(719,"(5(F21.12))")X(i),Y1(i),Y2(i),Y3(i),Y4(i)
     enddo
  elseif(present(Y3))then
     do i=1,Np
        write(719,"(4(F21.12))")X(i),Y1(i),Y2(i),Y3(i)
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
  integer                                :: i,Np
  character(len=*)                       :: pname
  real(8),dimension(:)                   :: X
  complex(8),dimension(size(X))          :: Y1
  complex(8),dimension(size(X)),optional :: Y2,Y3,Y4
  logical,optional                       :: append
  logical                                :: check,rw
  rw=.false.;if(present(append))rw=append
  Np=size(X)
  if(rw)then
     inquire(file=trim(adjustl(trim(pname))),exist=check)
     open(719,file=adjustl(trim(pname)),position="append")
     if(check)write(719,*)
  else
     open(719,file=adjustl(trim(pname)))
  endif
  if(present(Y4))then
     do i=1,Np
        write(719,"(F21.12,8(F21.12))")X(i),dimag(Y1(i)),dreal(Y1(i)),&
             dimag(Y2(i)),dreal(Y2(i)),dimag(Y3(i)),dreal(Y3(i)),&
             dimag(Y4(i)),dreal(Y4(i))
     enddo
  elseif(present(Y3))then
     do i=1,Np
        write(719,"(F21.12,6(F21.12))")X(i),dimag(Y1(i)),dreal(Y1(i)),&
             dimag(Y2(i)),dreal(Y2(i)),dimag(Y3(i)),dreal(Y3(i))
     enddo
  elseif(present(Y2))then
     do i=1,Np
        write(719,"(F21.12,4(F21.12))")X(i),dimag(Y1(i)),dreal(Y1(i)),&
             dimag(Y2(i)),dreal(Y2(i))
     enddo
  else
     do i=1,Np
        write(719,*)X(i),dimag(Y1(i)),dreal(Y1(i))
     enddo
  endif
  close(719)
end subroutine splotV_RC

!----------------------------


