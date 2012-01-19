!+-----------------------------------------------------------------+
!PROGRAM  : SPLOT
!TYPE     : Subroutine
!PURPOSE  : 
!+-----------------------------------------------------------------+
subroutine splotP_II(pname,X,Y1,Y2,append)
  character(len=*) :: pname
  integer          :: X
  integer          :: Y1
  integer,optional :: Y2
  logical,optional :: append
  if(.not.present(append))then
     open(719,file=trim(adjustl(trim(pname))))
  else
     open(719,file=trim(adjustl(trim(pname))),access="append")
     !write(719,*)""
  endif
  if(present(Y2))then
     write(719,*)X,Y1,Y2
  else
     write(719,*)X,Y1
  endif
  close(719)
end subroutine splotP_II
!----------------------------
!----------------------------
!----------------------------
subroutine splotP_IR(pname,X,Y1,Y2,append)
  character(len=*) :: pname
  integer          :: X
  real(8)          :: Y1
  real(8),optional :: Y2
  logical,optional :: append
  if(.not.present(append))then
     open(719,file=adjustl(trim(pname)))
  else
     open(719,file=adjustl(trim(pname)),access="append")
     !write(719,*)""
  endif
  if(present(Y2))then
     write(719,*)X,Y1,Y2
  else
     write(719,*)X,Y1
  endif
  close(719)
end subroutine splotP_IR
!----------------------------
!----------------------------
!----------------------------
subroutine splotP_RR(pname,X,Y1,Y2,append)
  character(len=*) :: pname
  real(8)          :: X
  real(8)          :: Y1
  real(8),optional :: Y2
  logical,optional :: append
  if(.not.present(append))then
     open(719,file=adjustl(trim(pname)))
  else
     open(719,file=adjustl(trim(pname)),access="append")
     !write(719,*)""
  endif
  if(present(Y2))then
     write(719,*)X,Y1,Y2
  else
     write(719,*)X,Y1
  endif
  close(719)
end subroutine splotP_RR
!----------------------------
!----------------------------
!----------------------------
subroutine splotP_RC(pname,X,Y1,append)
  character(len=*) :: pname
  real(8)          :: X
  complex(8)       :: Y1
  logical,optional :: append
  if(.not.present(append))then
     open(719,file=adjustl(trim(pname)))
  else
     open(719,file=adjustl(trim(pname)),access="append")
     !write(719,*)""
  endif
  write(719,*)X,real(Y1),aimag(Y1)
  close(719)
end subroutine splotP_RC
!----------------------------
!----------------------------
!----------------------------

subroutine splotV_C(pname,X,Y1,append)
  integer                       :: i,j,Np,Ny1,Ny2,Nx
  character(len=*)              :: pname
  real(8),dimension(:)          :: X
  complex(8),dimension(:)       :: Y1
  logical,optional              :: append
  Nx=size(X) ; Ny1=size(Y1)
  Np=min(Nx,Ny1)
  if(Nx/=Ny1)write(*,"(a,1x,I6,I6,I6)")"problem while printing "//trim(pname)//" Nx,Ny",Nx,Ny1
  if(.not.present(append))then
     open(719,file=trim(adjustl(trim(pname))))
  else
     open(719,file=trim(adjustl(trim(pname))),access="append")
     write(719,*)""
  endif
  do i=1,Np
     write(719,*)X(i),aimag(Y1(i)),real(Y1(i))
  enddo
  if(present(append))write(719,*)""
  close(719)
end subroutine splotV_C
!----------------------------
!----------------------------
!----------------------------
subroutine splotV_R(pname,X,Y1,Y2,append)
  integer                       :: i,j,Np,Ny1,Ny2,Nx
  character(len=*)              :: pname
  real(8),dimension(:)          :: X
  real(8),dimension(:)          :: Y1
  real(8),dimension(:),optional :: Y2
  logical,optional              :: append
  Nx=size(X) ; Ny1=size(Y1);   if(present(Y2))Ny2=size(Y2)
  Np=min(Nx,Ny1);if(present(Y2))Np=min(Nx,Ny1,Ny2)
  if(Nx/=Ny1.OR.(present(Y2).AND.Nx/=Ny2))write(*,"(a,1x,I6,I6,I6)")"problem while printing "//trim(pname)//" Nx,Ny1,Ny2",Nx,Ny1,Ny2
  if(.not.present(append))then
     open(719,file=adjustl(trim(pname)))
  else
     open(719,file=adjustl(trim(pname)),access="append")
     write(719,*)""
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
  if(present(append))write(719,*)""
  close(719)
end subroutine splotV_R
!----------------------------
!----------------------------
!----------------------------

subroutine splotM_C(pname,Y1,X,append)
  integer                       :: i,j,Np,Ny1,Ny2,Nx
  character(len=*)              :: pname
  real(8),optional,dimension(:) :: X
  complex(8),dimension(:,:)     :: Y1
  logical,optional              :: append
  Ny1=size(Y1,1) ; Ny2=size(Y1,2)
  if(present(X))then
     Nx=size(X)
     if(Nx < Ny1) then
        write(*,*)"big problem while printing "//trim(pname)//": skipping"
        return
     endif
     if(Nx/=Ny1 .OR. Nx/=Ny2) write(*,"(a,1x,I6,I6,I6)")"problem while printing "//trim(pname)//" Nx,Ny1,Ny2",Nx,Ny1,Ny2
  endif

  if(.not.present(append))then
     open(719,file=adjustl(trim(pname)))
  else
     open(719,file=adjustl(trim(pname)),access="append")
     write(719,*)""
  endif

  do i=1,Ny1
     do j=1,Ny2
        if(present(X))then
           write(719,*)X(i),Y1(i,j)
        else
           write(719,*)Y1(i,j)
        endif
     enddo
     write(719,*)""
  enddo
  if(present(append))write(719,*)""
  close(719)
end subroutine splotM_C
!----------------------------
!----------------------------
!----------------------------
subroutine splotM_R(pname,Y1,X,append)
  integer                       :: i,j,Np,Ny1,Ny2,Nx
  character(len=*)              :: pname
  real(8),optional,dimension(:) :: X
  real(8),dimension(:,:)        :: Y1
  logical,optional              :: append
  Ny1=size(Y1,1) ; Ny2=size(Y1,2)
  if(present(X))then
     Nx=size(X)
     if(Nx < Ny1) then
        write(*,*)"big problem while printing "//trim(pname)//": skipping"
        return
     endif
     if(Nx/=Ny1 .OR. Nx/=Ny2) write(*,"(a,1x,I6,I6,I6)")"problem while printing "//trim(pname)//" Nx,Ny1,Ny2",Nx,Ny1,Ny2
  endif

  if(.not.present(append))then
     open(719,file=adjustl(trim(pname)))
  else
     open(719,file=adjustl(trim(pname)),access="append")
     write(719,*)""
  endif

  do i=1,Ny1
     do j=1,Ny2
        if(present(X))then
           write(719,*)X(i),Y1(i,j)
        else
           write(719,*)Y1(i,j)
        endif
     enddo
     write(719,*)""
  enddo
  if(present(append))write(719,*)""
  close(719)
end subroutine splotM_R
!----------------------------
!----------------------------
!----------------------------

subroutine splot3D(pname,X1,X2,Y,wlines,nlines)
  integer                              :: i,j,Ny,Nx1,Nx2,count,Nl
  character(len=*)                     :: pname
  real(8),dimension(:)                 :: X1
  real(8),dimension(:)                 :: X2
  real(8),dimension(size(X1),size(X2)) :: Y
  integer,optional                     :: nlines
  logical,optional                     :: wlines
  Nx1=size(X1) ; Nx2=size(X2)
  Nl=5; if(present(nlines))Nl=nlines
  open(719,file=adjustl(trim(pname)))
  if(present(wlines))open(720,file=adjustl(trim(pname))//"_withlines")
  do i=1,Nx1
     count=mod(i,Nl)
     do j=1,Nx2
        write(719,*)X1(i),X2(j),Y(i,j)
        if(present(wlines).AND.count==0)write(720,*)X1(i),X2(j),Y(i,j)
     enddo
     write(719,*)""
     if(present(wlines))write(720,*)""
     if(present(wlines))write(720,*)""
  enddo
  close(719)
  if(present(wlines))close(720)
  call system("echo set nokey > plot_"//adjustl(trim(pname)) )
  call system("echo set grid >> plot_"//adjustl(trim(pname)) )
  call system("echo set view 50,170,1,1 >> plot_"//adjustl(trim(pname)) )
  !call system("echo set zlabel \'"//adjustl(trim(pname))//"\' >> plot_"//adjustl(trim(pname)) )
  call system("echo splot \'"//trim(pname)//"\' with pm3d >> plot_"//adjustl(trim(pname)) )
  if(present(wlines))call system("echo rep \'"//trim(pname)//"_withlines\' with lines >> plot_"//adjustl(trim(pname)) )
  call system("echo '#'set term png size 1024,768 >> plot_"//adjustl(trim(pname)) )
  call system("echo '#'set out\'"//adjustl(trim(pname))//".png\' >> plot_"//adjustl(trim(pname)) )
  call system("echo '#'rep >> plot_"//adjustl(trim(pname)) )
end subroutine splot3D

subroutine splot3D_(pname,X1,X2,Y,wlines,nlines)
  integer                              :: i,j,Ny,Nx1,Nx2,count,Nl,Nk
  character(len=*)                     :: pname
  real(8),dimension(:)                 :: X1 !(0:Nt)
  real(8),dimension(:,:)               :: X2 !(0:Nt,Lk)
  real(8),dimension(size(X2,1),size(X2,2)) :: Y  !(0:Nt,Lk)
  integer,optional                     :: nlines
  logical,optional                     :: wlines
  Nx1=size(X1) ; Nk=size(X2,2)
  Nl=5; if(present(nlines))Nl=nlines
  open(719,file=adjustl(trim(pname)))
  if(present(wlines))open(720,file=adjustl(trim(pname))//"_withlines")

  do i=1,Nx1
     count=mod(i,Nl)
     do j=1,Nk
        write(719,*)X1(i),X2(i,j),Y(i,j)
        if(present(wlines).AND.count==0)write(720,*)X1(i),X2(i,j),Y(i,j)
     enddo
     write(719,*)""
     if(present(wlines))write(720,*)""
     if(present(wlines))write(720,*)""
  enddo
  close(719)
  if(present(wlines))close(720)
  call system("echo set nokey > plot_"//adjustl(trim(pname)) )
  call system("echo set grid >> plot_"//adjustl(trim(pname)) )
  call system("echo set view 50,170,1,1 >> plot_"//adjustl(trim(pname)) )
  !call system("echo set zlabel \'"//adjustl(trim(pname))//"\' >> plot_"//adjustl(trim(pname)) )
  call system("echo splot \'"//trim(pname)//"\' with pm3d >> plot_"//adjustl(trim(pname)) )
  if(present(wlines))call system("echo rep \'"//trim(pname)//"_withlines\' with lines >> plot_"//adjustl(trim(pname)) )
  call system("echo '#'set term png size 1024,768 >> plot_"//adjustl(trim(pname)) )
  call system("echo '#'set out\'"//adjustl(trim(pname))//".png\' >> plot_"//adjustl(trim(pname)) )
  call system("echo '#'rep >> plot_"//adjustl(trim(pname)) )
end subroutine splot3D_
!********************************************************************
!********************************************************************
!********************************************************************

