

!----------------------------

subroutine splotP_II(pname,X,Y1,Y2,Y3,Y4,Y5,Y6,Y7,Y8,append)
  character(len=*) :: pname
  integer          :: X
  integer          :: Y1
  integer,optional :: Y2,Y3,Y4,Y5,Y6,Y7,Y8
  logical,optional :: append
<<<<<<< HEAD
  if(present(append).AND. append.eqv..true.)then
=======
  logical          :: rw
  rw=.false.;if(present(append))rw=append
  if(rw)then
>>>>>>> devel
     open(719,file=adjustl(trim(pname)),position="append")
  else
     open(719,file=adjustl(trim(pname)))
  endif
  if(present(Y8))then
     write(719,"(9(I15))")X,Y1,Y2,Y3,Y4,Y5,Y6,Y7,Y8
  elseif(present(Y7))then
     write(719,"(8(I15))")X,Y1,Y2,Y3,Y4,Y5,Y6,Y7
  elseif(present(Y6))then
     write(719,"(7(I15))")X,Y1,Y2,Y3,Y4,Y5,Y6
  elseif(present(Y5))then
     write(719,"(6(I15))")X,Y1,Y2,Y3,Y4,Y5
  elseif(present(Y4))then
     write(719,"(5(I15))")X,Y1,Y2,Y3,Y4
  elseif(present(Y3))then
     write(719,"(4(I15))")X,Y1,Y2,Y3
  elseif(present(Y2))then
     write(719,*)X,Y1,Y2
  else
     write(719,*)X,Y1
  endif
  close(719)
end subroutine splotP_II

!----------------------------

subroutine splotP_IR(pname,X,Y1,Y2,Y3,Y4,Y5,Y6,Y7,Y8,append)
  character(len=*) :: pname
  integer          :: X
  real(8)          :: Y1
  real(8),optional :: Y2,Y3,Y4,Y5,Y6,Y7,Y8
  logical,optional :: append
<<<<<<< HEAD
  if(present(append).AND. append.eqv..true.)then
=======
  logical          :: rw
  rw=.false.;if(present(append))rw=append
  if(rw)then
>>>>>>> devel
     open(719,file=adjustl(trim(pname)),position="append")
  else
     open(719,file=adjustl(trim(pname)))
  endif
  if(present(Y8))then
     write(719,"(I15,8(F18.10))")X,Y1,Y2,Y3,Y4,Y5,Y6,Y7,Y8
  elseif(present(Y7))then
     write(719,"(I15,7(F18.10))")X,Y1,Y2,Y3,Y4,Y5,Y6,Y7
  elseif(present(Y6))then
     write(719,"(I15,6(F18.10))")X,Y1,Y2,Y3,Y4,Y5,Y6
  elseif(present(Y5))then
     write(719,"(I15,5(F18.10))")X,Y1,Y2,Y3,Y4,Y5
  elseif(present(Y4))then
     write(719,"(I15,4(F18.10))")X,Y1,Y2,Y3,Y4
  elseif(present(Y3))then
     write(719,"(I15,3(F18.10))")X,Y1,Y2,Y3
  elseif(present(Y2))then
     write(719,*)X,Y1,Y2
  else
     write(719,*)X,Y1
  endif
  close(719)
end subroutine splotP_IR

!----------------------------

subroutine splotP_IC(pname,X,Y1,Y2,Y3,Y4,append)
  character(len=*) :: pname
  integer          :: X
  complex(8)       :: Y1
  complex(8),optional:: Y2,Y3,Y4
  logical,optional :: append
<<<<<<< HEAD
  if(present(append).AND. append.eqv..true.)then
=======
  logical          :: rw
  rw=.false.;if(present(append))rw=append
  if(rw)then
>>>>>>> devel
     open(719,file=adjustl(trim(pname)),position="append")
  else
     open(719,file=adjustl(trim(pname)))
  endif
  if(present(Y4))then
     write(719,"(I15,8(F18.10))")X,dreal(Y1),dimag(Y1),&
          dreal(Y2),dimag(Y2),dreal(Y3),dimag(Y3),dreal(Y4),dimag(Y4)
  elseif(present(Y3))then
     write(719,"(I15,6(F18.10))")X,dreal(Y1),dimag(Y1),&
          dreal(Y2),dimag(Y2),dreal(Y3),dimag(Y3)
  elseif(present(Y2))then
     write(719,"(I15,4(F18.10))")X,dreal(Y1),dimag(Y1),dreal(Y2),dimag(Y2)
  else
     write(719,*)X,dreal(Y1),aimag(Y1)
  endif
  close(719)
end subroutine splotP_IC

!----------------------------

subroutine splotP_RI(pname,X,Y1,Y2,Y3,Y4,Y5,Y6,Y7,Y8,append)
  character(len=*) :: pname
  real(8)          :: X
  integer          :: Y1
  integer,optional :: Y2,Y3,Y4,Y5,Y6,Y7,Y8
  logical,optional :: append
<<<<<<< HEAD
  if(present(append).AND. append.eqv..true.)then
=======
  logical          :: rw
  rw=.false.;if(present(append))rw=append
  if(rw)then
>>>>>>> devel
     open(719,file=adjustl(trim(pname)),position="append")
  else
     open(719,file=adjustl(trim(pname)))
  endif
  if(present(Y8))then
     write(719,"(F18.10,8I15)")X,Y1,Y2,Y3,Y4,Y5,Y6,Y7,Y8
  elseif(present(Y7))then
     write(719,"(F18.10,7I15)")X,Y1,Y2,Y3,Y4,Y5,Y6,Y7
  elseif(present(Y6))then
     write(719,"(F18.10,6I15)")X,Y1,Y2,Y3,Y4,Y5,Y6
  elseif(present(Y5))then
     write(719,"(F18.10,5I15)")X,Y1,Y2,Y3,Y4,Y5
  elseif(present(Y4))then
     write(719,"(F18.10,4I15)")X,Y1,Y2,Y3,Y4
  elseif(present(Y3))then
     write(719,"(F18.10,3I15)")X,Y1,Y2,Y3
  elseif(present(Y2))then
     write(719,*)X,Y1,Y2
  else
     write(719,*)X,Y1
  endif
  close(719)
end subroutine splotP_RI

!----------------------------

subroutine splotP_RR(pname,X,Y1,Y2,Y3,Y4,Y5,Y6,Y7,Y8,append)
  character(len=*) :: pname
  real(8)          :: X
  real(8)          :: Y1
  real(8),optional :: Y2,Y3,Y4,Y5,Y6,Y7,Y8
  logical,optional :: append
<<<<<<< HEAD
  if(present(append).AND. append.eqv..true.)then
=======
  logical          :: rw
  rw=.false.;if(present(append))rw=append
  if(rw)then
>>>>>>> devel
     open(719,file=adjustl(trim(pname)),position="append")
  else
     open(719,file=adjustl(trim(pname)))
  endif
  if(present(Y8))then
     write(719,"(9F18.10)")X,Y1,Y2,Y3,Y4,Y5,Y6,Y7,Y8
  elseif(present(Y7))then
     write(719,"(8F18.10)")X,Y1,Y2,Y3,Y4,Y5,Y6,Y7
  elseif(present(Y6))then
     write(719,"(7F18.10)")X,Y1,Y2,Y3,Y4,Y5,Y6
  elseif(present(Y5))then
     write(719,"(6F18.10)")X,Y1,Y2,Y3,Y4,Y5
  elseif(present(Y4))then
     write(719,"(5F18.10)")X,Y1,Y2,Y3,Y4
  elseif(present(Y3))then
     write(719,"(4F18.10)")X,Y1,Y2,Y3
  elseif(present(Y2))then
     write(719,*)X,Y1,Y2
  else
     write(719,*)X,Y1
  endif
  close(719)
end subroutine splotP_RR

!----------------------------

subroutine splotP_RC(pname,X,Y1,Y2,Y3,Y4,append)
  character(len=*)    :: pname
  real(8)             :: X
  complex(8)          :: Y1
  complex(8),optional :: Y2,Y3,Y4
  logical,optional    :: append
  if(present(append).AND. append.eqv..true.)then
     open(719,file=adjustl(trim(pname)),position="append")
  else
     open(719,file=adjustl(trim(pname)))
  endif
  if(present(Y4))then
     write(719,"(F18.10,8(F18.10))")X,dreal(Y1),dimag(Y1),&
          dreal(Y2),dimag(Y2),dreal(Y3),dimag(Y3),dreal(Y4),dimag(Y4)
  elseif(present(Y3))then
     write(719,"(F18.10,6(F18.10))")X,dreal(Y1),dimag(Y1),&
          dreal(Y2),dimag(Y2),dreal(Y3),dimag(Y3)
  elseif(present(Y2))then
     write(719,"(F18.10,4(F18.10))")X,dreal(Y1),dimag(Y1),dreal(Y2),dimag(Y2)
  else
     write(719,*)X,dreal(Y1),dimag(Y1)
  endif
  close(719)
end subroutine splotP_RC
