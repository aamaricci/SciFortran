

!----------------------------

subroutine splotP_II(pname,X,Y1,Y2,Y3,Y4,Y5,Y6,Y7,Y8,append)
  character(len=*) :: pname
  integer          :: X
  integer          :: Y1
  integer,optional :: Y2,Y3,Y4,Y5,Y6,Y7,Y8
  logical,optional :: append
  logical          :: rw
  rw=.false.;if(present(append))rw=append
  if(rw)then
     open(719,file=reg((pname)),position="append")
  else
     open(719,file=reg((pname)))
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
  logical          :: rw
  rw=.false.;if(present(append))rw=append
  if(rw)then

     open(719,file=reg((pname)),position="append")
  else
     open(719,file=reg((pname)))
  endif
  if(present(Y8))then
     write(719,"(I15,8(F21.12))")X,Y1,Y2,Y3,Y4,Y5,Y6,Y7,Y8
  elseif(present(Y7))then
     write(719,"(I15,7(F21.12))")X,Y1,Y2,Y3,Y4,Y5,Y6,Y7
  elseif(present(Y6))then
     write(719,"(I15,6(F21.12))")X,Y1,Y2,Y3,Y4,Y5,Y6
  elseif(present(Y5))then
     write(719,"(I15,5(F21.12))")X,Y1,Y2,Y3,Y4,Y5
  elseif(present(Y4))then
     write(719,"(I15,4(F21.12))")X,Y1,Y2,Y3,Y4
  elseif(present(Y3))then
     write(719,"(I15,3(F21.12))")X,Y1,Y2,Y3
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
  logical          :: rw
  rw=.false.;if(present(append))rw=append
  if(rw)then

     open(719,file=reg((pname)),position="append")
  else
     open(719,file=reg((pname)))
  endif
  if(present(Y4))then
     write(719,"(I15,8(F21.12))")X,dreal(Y1),dimag(Y1),&
          dreal(Y2),dimag(Y2),dreal(Y3),dimag(Y3),dreal(Y4),dimag(Y4)
  elseif(present(Y3))then
     write(719,"(I15,6(F21.12))")X,dreal(Y1),dimag(Y1),&
          dreal(Y2),dimag(Y2),dreal(Y3),dimag(Y3)
  elseif(present(Y2))then
     write(719,"(I15,4(F21.12))")X,dreal(Y1),dimag(Y1),dreal(Y2),dimag(Y2)
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
  logical          :: rw
  rw=.false.;if(present(append))rw=append
  if(rw)then

     open(719,file=reg((pname)),position="append")
  else
     open(719,file=reg((pname)))
  endif
  if(present(Y8))then
     write(719,"(F21.12,8I15)")X,Y1,Y2,Y3,Y4,Y5,Y6,Y7,Y8
  elseif(present(Y7))then
     write(719,"(F21.12,7I15)")X,Y1,Y2,Y3,Y4,Y5,Y6,Y7
  elseif(present(Y6))then
     write(719,"(F21.12,6I15)")X,Y1,Y2,Y3,Y4,Y5,Y6
  elseif(present(Y5))then
     write(719,"(F21.12,5I15)")X,Y1,Y2,Y3,Y4,Y5
  elseif(present(Y4))then
     write(719,"(F21.12,4I15)")X,Y1,Y2,Y3,Y4
  elseif(present(Y3))then
     write(719,"(F21.12,3I15)")X,Y1,Y2,Y3
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
  logical          :: rw
  rw=.false.;if(present(append))rw=append
  if(rw)then

     open(719,file=reg((pname)),position="append")
  else
     open(719,file=reg((pname)))
  endif
  if(present(Y8))then
     write(719,"(9F21.12)")X,Y1,Y2,Y3,Y4,Y5,Y6,Y7,Y8
  elseif(present(Y7))then
     write(719,"(8F21.12)")X,Y1,Y2,Y3,Y4,Y5,Y6,Y7
  elseif(present(Y6))then
     write(719,"(7F21.12)")X,Y1,Y2,Y3,Y4,Y5,Y6
  elseif(present(Y5))then
     write(719,"(6F21.12)")X,Y1,Y2,Y3,Y4,Y5
  elseif(present(Y4))then
     write(719,"(5F21.12)")X,Y1,Y2,Y3,Y4
  elseif(present(Y3))then
     write(719,"(4F21.12)")X,Y1,Y2,Y3
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
     open(719,file=reg((pname)),position="append")
  else
     open(719,file=reg((pname)))
  endif
  if(present(Y4))then
     write(719,"(F21.12,8(F21.12))")X,dreal(Y1),dimag(Y1),&
          dreal(Y2),dimag(Y2),dreal(Y3),dimag(Y3),dreal(Y4),dimag(Y4)
  elseif(present(Y3))then
     write(719,"(F21.12,6(F21.12))")X,dreal(Y1),dimag(Y1),&
          dreal(Y2),dimag(Y2),dreal(Y3),dimag(Y3)
  elseif(present(Y2))then
     write(719,"(F21.12,4(F21.12))")X,dreal(Y1),dimag(Y1),dreal(Y2),dimag(Y2)
  else
     write(719,*)X,dreal(Y1),dimag(Y1)
  endif
  close(719)
end subroutine splotP_RC
