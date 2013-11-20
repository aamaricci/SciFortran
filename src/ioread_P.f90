
!----------------------------

subroutine sreadP_II(pname,X,Y1,Y2,Y3,Y4,Y5,Y6,Y7,Y8)
  character(len=*) :: pname
  integer          :: X
  integer          :: Y1
  integer,optional :: Y2,Y3,Y4,Y5,Y6,Y7,Y8
  open(719,file=trim(adjustl(trim(pname))))

  if(present(Y8))then
     read(719,"(9(I15))")X,Y1,Y2,Y3,Y4,Y5,Y6,Y7,Y8
  elseif(present(Y7))then
     read(719,"(8(I15))")X,Y1,Y2,Y3,Y4,Y5,Y6,Y7
  elseif(present(Y6))then
     read(719,"(7(I15))")X,Y1,Y2,Y3,Y4,Y5,Y6
  elseif(present(Y5))then
     read(719,"(6(I15))")X,Y1,Y2,Y3,Y4,Y5
  elseif(present(Y4))then
     read(719,"(5(I15))")X,Y1,Y2,Y3,Y4
  elseif(present(Y3))then
     read(719,"(4(I15))")X,Y1,Y2,Y3
  elseif(present(Y2))then
     read(719,*)X,Y1,Y2
  else
     read(719,*)X,Y1
  endif

  close(719)
end subroutine sreadP_II

!----------------------------

subroutine sreadP_IR(pname,X,Y1,Y2,Y3,Y4,Y5,Y6,Y7,Y8)
  character(len=*) :: pname
  integer          :: X
  real(8)          :: Y1
  real(8),optional :: Y2,Y3,Y4,Y5,Y6,Y7,Y8
  open(719,file=trim(adjustl(trim(pname))))
  if(present(Y8))then
     read(719,"(I15,8(F21.12))")X,Y1,Y2,Y3,Y4,Y5,Y6,Y7,Y8
  elseif(present(Y7))then
     read(719,"(I15,7(F21.12))")X,Y1,Y2,Y3,Y4,Y5,Y6,Y7
  elseif(present(Y6))then
     read(719,"(I15,6(F21.12))")X,Y1,Y2,Y3,Y4,Y5,Y6
  elseif(present(Y5))then
     read(719,"(I15,5(F21.12))")X,Y1,Y2,Y3,Y4,Y5
  elseif(present(Y4))then
     read(719,"(I15,4(F21.12))")X,Y1,Y2,Y3,Y4
  elseif(present(Y3))then
     read(719,"(I15,3(F21.12))")X,Y1,Y2,Y3
  elseif(present(Y2))then
     read(719,*)X,Y1,Y2
  else
     read(719,*)X,Y1
  endif
  close(719)
end subroutine sreadP_IR

!----------------------------

subroutine sreadP_IC(pname,X,Y1,Y2,Y3,Y4)
  character(len=*)    :: pname
  integer             :: X
  complex(8)          :: Y1
  complex(8),optional :: Y2,Y3,Y4
  real(8)             :: reY1,imY1
  real(8)             :: reY2,imY2,reY3,imY3,reY4,imY4
  open(719,file=trim(adjustl(trim(pname))))
  if(present(Y4))then
     read(719,"(I15,8(F21.12))")X,reY1,imY1,reY2,imY2,reY3,imY3,reY4,imY4
     Y1=cmplx(reY1,imY1);Y2=cmplx(reY2,imY2);Y3=cmplx(reY3,imY3);Y4=cmplx(reY4,imY4)
  elseif(present(Y3))then
     read(719,"(I15,8(F21.12))")X,reY1,imY1,reY2,imY2,reY3,imY3
     Y1=cmplx(reY1,imY1);Y2=cmplx(reY2,imY2);Y3=cmplx(reY3,imY3)
  elseif(present(Y2))then
     read(719,"(I15,8(F21.12))")X,reY1,imY1,reY2,imY2
     Y1=cmplx(reY1,imY1);Y2=cmplx(reY2,imY2)
  else
     read(719,*)X,reY1,imY1 ; Y1=cmplx(reY1,imY1)
  endif
  close(719)
end subroutine sreadP_IC

!----------------------------

subroutine sreadP_RI(pname,X,Y1,Y2,Y3,Y4,Y5,Y6,Y7,Y8)
  character(len=*) :: pname
  real(8)          :: X
  integer          :: Y1
  integer,optional :: Y2,Y3,Y4,Y5,Y6,Y7,Y8
  open(719,file=trim(adjustl(trim(pname))))
  if(present(Y8))then
     read(719,"(F21.12,8I15)")X,Y1,Y2,Y3,Y4,Y5,Y6,Y7,Y8
  elseif(present(Y7))then
     read(719,"(F21.12,7I15)")X,Y1,Y2,Y3,Y4,Y5,Y6,Y7
  elseif(present(Y6))then
     read(719,"(F21.12,6I15)")X,Y1,Y2,Y3,Y4,Y5,Y6
  elseif(present(Y5))then
     read(719,"(F21.12,5I15)")X,Y1,Y2,Y3,Y4,Y5
  elseif(present(Y4))then
     read(719,"(F21.12,4I15)")X,Y1,Y2,Y3,Y4
  elseif(present(Y3))then
     read(719,"(F21.12,3I15)")X,Y1,Y2,Y3
  elseif(present(Y2))then
     read(719,*)X,Y1,Y2
  else
     read(719,*)X,Y1
  endif
  close(719)
end subroutine sreadP_RI

!----------------------------

subroutine sreadP_RR(pname,X,Y1,Y2,Y3,Y4,Y5,Y6,Y7,Y8)
  character(len=*) :: pname
  real(8)          :: X
  real(8)          :: Y1
  real(8),optional :: Y2,Y3,Y4,Y5,Y6,Y7,Y8
  open(719,file=trim(adjustl(trim(pname))))
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
end subroutine sreadP_RR

!----------------------------

subroutine sreadP_RC(pname,X,Y1,Y2,Y3,Y4)
  character(len=*)    :: pname
  complex(8)          :: Y1
  complex(8),optional :: Y2,Y3,Y4
  real(8)             :: X,reY1,imY1
  real(8)             :: reY2,imY2,reY3,imY3,reY4,imY4
  open(719,file=trim(adjustl(trim(pname))))
  if(present(Y4))then
     read(719,"(I15,8(F21.12))")X,reY1,imY1,reY2,imY2,reY3,imY3,reY4,imY4
     Y1=cmplx(reY1,imY1);Y2=cmplx(reY2,imY2);Y3=cmplx(reY3,imY3);Y4=cmplx(reY4,imY4)
  elseif(present(Y3))then
     read(719,"(I15,8(F21.12))")X,reY1,imY1,reY2,imY2,reY3,imY3
     Y1=cmplx(reY1,imY1);Y2=cmplx(reY2,imY2);Y3=cmplx(reY3,imY3)
  elseif(present(Y2))then
     read(719,"(I15,8(F21.12))")X,reY1,imY1,reY2,imY2
     Y1=cmplx(reY1,imY1);Y2=cmplx(reY2,imY2)
  else
     read(719,*)X,reY1,imY1 ; Y1=cmplx(reY1,imY1)
  endif
  close(719)
end subroutine sreadP_RC
