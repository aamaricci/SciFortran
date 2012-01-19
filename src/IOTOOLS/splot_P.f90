
!----------------------------

subroutine splotP_II(pname,X,Y1,Y2,Y3,Y4,Y5,Y6,Y7,Y8,append)
  character(len=*) :: pname
  integer          :: X
  integer          :: Y1
  integer,optional :: Y2,Y3,Y4,Y5,Y6,Y7,Y8
  logical,optional :: append
  if(present(append).AND. append==.true.)then
     open(719,file=adjustl(trim(pname)),access="append")
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
  if(present(append).AND. append==.true.)then
     open(719,file=adjustl(trim(pname)),access="append")
  else
     open(719,file=adjustl(trim(pname)))
  endif
  if(present(Y8))then
     write(719,"(I15,8(F14.9))")X,Y1,Y2,Y3,Y4,Y5,Y6,Y7,Y8
  elseif(present(Y7))then
     write(719,"(I15,7(F14.9))")X,Y1,Y2,Y3,Y4,Y5,Y6,Y7
  elseif(present(Y6))then
     write(719,"(I15,6(F14.9))")X,Y1,Y2,Y3,Y4,Y5,Y6
  elseif(present(Y5))then
     write(719,"(I15,5(F14.9))")X,Y1,Y2,Y3,Y4,Y5
  elseif(present(Y4))then
     write(719,"(I15,4(F14.9))")X,Y1,Y2,Y3,Y4
  elseif(present(Y3))then
     write(719,"(I15,3(F14.9))")X,Y1,Y2,Y3
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
  if(present(append).AND. append==.true.)then
     open(719,file=adjustl(trim(pname)),access="append")
  else
     open(719,file=adjustl(trim(pname)))
  endif
  if(present(Y4))then
     write(719,"(I15,8(F14.9))")X,real(Y1,8),dimag(Y1),&
          real(Y2,8),dimag(Y2),real(Y3,8),dimag(Y3),real(Y4,8),dimag(Y4)
  elseif(present(Y3))then
     write(719,"(I15,6(F14.9))")X,real(Y1,8),dimag(Y1),&
          real(Y2,8),dimag(Y2),real(Y3,8),dimag(Y3)
  elseif(present(Y2))then
     write(719,"(I15,4(F14.9))")X,real(Y1,8),dimag(Y1),real(Y2,8),dimag(Y2)
  else
     write(719,*)X,real(Y1,8),aimag(Y1)
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
  if(present(append).AND. append==.true.)then
     open(719,file=adjustl(trim(pname)),access="append")
  else
     open(719,file=adjustl(trim(pname)))
  endif
  if(present(Y8))then
     write(719,"(F14.9,8I15)")X,Y1,Y2,Y3,Y4,Y5,Y6,Y7,Y8
  elseif(present(Y7))then
     write(719,"(F14.9,7I15)")X,Y1,Y2,Y3,Y4,Y5,Y6,Y7
  elseif(present(Y6))then
     write(719,"(F14.9,6I15)")X,Y1,Y2,Y3,Y4,Y5,Y6
  elseif(present(Y5))then
     write(719,"(F14.9,5I15)")X,Y1,Y2,Y3,Y4,Y5
  elseif(present(Y4))then
     write(719,"(F14.9,4I15)")X,Y1,Y2,Y3,Y4
  elseif(present(Y3))then
     write(719,"(F14.9,3I15)")X,Y1,Y2,Y3
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
  if(present(append).AND. append==.true.)then
     open(719,file=adjustl(trim(pname)),access="append")
  else
     open(719,file=adjustl(trim(pname)))
  endif
  if(present(Y8))then
     write(719,"(9F14.9)")X,Y1,Y2,Y3,Y4,Y5,Y6,Y7,Y8
  elseif(present(Y7))then
     write(719,"(8F14.9)")X,Y1,Y2,Y3,Y4,Y5,Y6,Y7
  elseif(present(Y6))then
     write(719,"(7F14.9)")X,Y1,Y2,Y3,Y4,Y5,Y6
  elseif(present(Y5))then
     write(719,"(6F14.9)")X,Y1,Y2,Y3,Y4,Y5
  elseif(present(Y4))then
     write(719,"(5F14.9)")X,Y1,Y2,Y3,Y4
  elseif(present(Y3))then
     write(719,"(4F14.9)")X,Y1,Y2,Y3
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
  if(present(append).AND. append==.true.)then
     open(719,file=adjustl(trim(pname)),access="append")
  else
     open(719,file=adjustl(trim(pname)))
  endif
  if(present(Y4))then
     write(719,"(F14.9,8(F14.9))")X,real(Y1,8),dimag(Y1),&
          real(Y2,8),dimag(Y2),real(Y3,8),dimag(Y3),real(Y4,8),dimag(Y4)
  elseif(present(Y3))then
     write(719,"(F14.9,6(F14.9))")X,real(Y1,8),dimag(Y1),&
          real(Y2,8),dimag(Y2),real(Y3,8),dimag(Y3)
  elseif(present(Y2))then
     write(719,"(F14.9,4(F14.9))")X,real(Y1,8),dimag(Y1),real(Y2,8),dimag(Y2)
  else
     write(719,*)X,real(Y1,8),aimag(Y1)
  endif
  close(719)
end subroutine splotP_RC
