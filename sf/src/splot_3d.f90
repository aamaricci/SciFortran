subroutine splot3D(pname,X1,X2,Y,wlines,nlines)
  integer                              :: i,j,Nx1,Nx2,count,Nl
  character(len=*)                     :: pname
  real(8),dimension(:)                 :: X1
  real(8),dimension(:)                 :: X2
  real(8),dimension(size(X1),size(X2)) :: Y
  integer,optional                     :: nlines
  logical,optional                     :: wlines
  real(8)                              :: X1min,X1max
  real(8)                              :: X2min,X2max
  character(len=12)                    :: minx,miny,maxx,maxy
  character(len=256)                   :: fname,dname
  fname=get_filename(pname)
  dname=get_filepath(pname)
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
  X1min=minval(X1);X1max=maxval(X1)
  X2min=minval(X2);X2max=maxval(X2)
  write(minx,"(f12.4)")X1min
  write(maxx,"(f12.4)")X1max
  write(miny,"(f12.4)")X2min
  write(maxy,"(f12.4)")X2max
  if(present(wlines))close(720)

  open(10,file=adjustl(trim(pname))//".gp")
  write(10,*)"gnuplot -persist << EOF"
  write(10,*)"set term wxt"
  write(10,*)"set title '"//trim(fname)//"'"
  write(10,*)"set pm3d map"
  write(10,*)"set size square"
  write(10,*)"set xrange ["//trim(adjustl(trim(minx)))//":"//trim(adjustl(trim(maxx)))//"]"
  write(10,*)"set yrange ["//trim(adjustl(trim(miny)))//":"//trim(adjustl(trim(maxy)))//"]"
  write(10,*)"splot '"//trim(fname)//"'"
  write(10,*)"#set term png size 1920,1280"
  write(10,*)"#set out '"//adjustl(trim(fname))//".png'"
  write(10,*)"#rep"
  write(10,*)"#"
  write(10,"(A)")"EOF"
  !
  write(10,*)"gnuplot -persist << EOF"
  write(10,*)"set term wxt"
  write(10,*)"set title '"//trim(fname)//"'"
  write(10,*)"set nokey"
  write(10,*)"set grid"
  write(10,*)"set view 50,10,1,1"
  write(10,*)"splot '"//trim(fname)//"' with pm3d"
  if(present(wlines))write(10,*)"rep '"//trim(pname)//"_withlines' with lines"
  write(10,*)"#set term png size 1920,1280"
  write(10,*)"#set out '"//adjustl(trim(fname))//".png'"
  write(10,*)"#rep"
  write(10,"(A)")"EOF"
  close(10)
  call system("chmod +x "//adjustl(trim(pname))//".gp")

end subroutine splot3D


! subroutine splot3D_(pname,X1,X2,Y,wlines,nlines)
!   integer                              :: i
!   character(len=*)                     :: pname
!   real(8),dimension(:)                 :: X1
!   real(8),dimension(:)                 :: X2
!   complex(8),dimension(size(X1),size(X2)) :: Y
!   integer,optional                     :: nlines
!   logical,optional                     :: wlines
!   character(len=256)                   :: fname,dname
!   fname=get_filename(pname)
!   dname=get_filepath(pname)
!   if(present(nlines).and.present(wlines))then
!      call splot3d(reg_filename(dname)//"re_"//reg_filename(fname),X1,X2,real(Y,8),wlines,nlines)
!      call splot3d(reg_filename(dname)//"im_"//reg_filename(fname),X1,X2,dimag(Y),wlines,nlines)
!   elseif(present(wlines))then
!      call splot3d(reg_filename(dname)//"re_"//reg_filename(fname),X1,X2,real(Y,8),wlines)
!      call splot3d(reg_filename(dname)//"im_"//reg_filename(fname),X1,X2,dimag(Y),wlines)
!   elseif(present(nlines))then
!      call splot3d(reg_filename(dname)//"re_"//reg_filename(fname),X1,X2,real(Y,8),nlines=nlines)
!      call splot3d(reg_filename(dname)//"im_"//reg_filename(fname),X1,X2,dimag(Y),nlines=nlines)
!   else
!      call splot3d(reg_filename(dname)//"re_"//reg_filename(fname),X1,X2,real(Y,8))
!      call splot3d(reg_filename(dname)//"im_"//reg_filename(fname),X1,X2,dimag(Y))
!   endif
! end subroutine splot3D_

subroutine splot3D_(pname,X1,X2,Y,wlines,nlines)
  integer                              :: i,j,Nx1,Nx2,count,Nl
  character(len=*)                     :: pname
  real(8),dimension(:)                 :: X1
  real(8),dimension(:)                 :: X2
  complex(8),dimension(size(X1),size(X2)) :: Y
  integer,optional                     :: nlines
  logical,optional                     :: wlines
  real(8)                              :: X1min,X1max
  real(8)                              :: X2min,X2max
  character(len=12)                    :: minx,miny,maxx,maxy
  character(len=256)                   :: fname,dname
  fname=get_filename(pname)
  dname=get_filepath(pname)
  Nx1=size(X1) ; Nx2=size(X2)
  Nl=5; if(present(nlines))Nl=nlines

  open(619,file=adjustl(trim(dname))//"re_"//adjustl(trim(fname)))
  open(719,file=adjustl(trim(dname))//"im_"//adjustl(trim(fname)))
  if(present(wlines))then
     open(620,file=adjustl(trim(dname))//"re_"//adjustl(trim(fname))//"_withlines")
     open(720,file=adjustl(trim(dname))//"im_"//adjustl(trim(fname))//"_withlines")
  endif
  do i=1,Nx1
     count=mod(i,Nl)
     do j=1,Nx2
        write(619,*)X1(i),X2(j),real(Y(i,j),8)
        write(719,*)X1(i),X2(j),aimag(Y(i,j))
        if(present(wlines).AND.count==0)write(620,*)X1(i),X2(j),real(Y(i,j),8)
        if(present(wlines).AND.count==0)write(720,*)X1(i),X2(j),aimag(Y(i,j))
     enddo
     write(619,*)""
     write(719,*)""
     if(present(wlines))then
        write(620,*)""
        write(620,*)""
        write(720,*)""
        write(720,*)""
     endif
  enddo
  close(619)
  close(719)

  X1min=minval(X1);X1max=maxval(X1)
  X2min=minval(X2);X2max=maxval(X2)
  write(minx,"(f12.4)")X1min
  write(maxx,"(f12.4)")X1max
  write(miny,"(f12.4)")X2min
  write(maxy,"(f12.4)")X2max
  if(present(wlines))then
     close(620)
     close(720)
  endif

  open(10,file=adjustl(trim(pname))//".gp")
  !Re:
  write(10,*)"gnuplot -persist << EOF"
  write(10,*)"set term wxt"
  write(10,*)"set title 'Re_"//trim(fname)//"'"
  write(10,*)"set pm3d map"
  write(10,*)"set size square"
  write(10,*)"set xrange ["//trim(adjustl(trim(minx)))//":"//trim(adjustl(trim(maxx)))//"]"
  write(10,*)"set yrange ["//trim(adjustl(trim(miny)))//":"//trim(adjustl(trim(maxy)))//"]"
  write(10,*)"splot 're_"//trim(fname)//"'"
  write(10,*)"#set term png size 1920,1280"
  write(10,*)"#set out 're_"//adjustl(trim(fname))//".png'"
  write(10,*)"#rep"
  write(10,"(A)")"EOF"
  !Im
  write(10,*)"gnuplot -persist << EOF"
  write(10,*)"set term wxt"
  write(10,*)"set title 'Im_"//trim(fname)//"'"
  write(10,*)"set pm3d map"
  write(10,*)"set size square"
  write(10,*)"set xrange ["//trim(adjustl(trim(minx)))//":"//trim(adjustl(trim(maxx)))//"]"
  write(10,*)"set yrange ["//trim(adjustl(trim(miny)))//":"//trim(adjustl(trim(maxy)))//"]"
  write(10,*)"splot 'im_"//trim(fname)//"'"
  write(10,*)"#set term png size 1920,1280"
  write(10,*)"#set out 'im_"//adjustl(trim(fname))//".png'"
  write(10,*)"#rep"
  write(10,"(A)")"EOF"
  !
  write(10,*)"gnuplot -persist << EOF"
  write(10,*)"set term wxt"
  write(10,*)"set title 'Re_"//trim(fname)//"'"
  write(10,*)"set nokey"
  write(10,*)"set grid"
  write(10,*)"set view 50,10,1,1"
  write(10,*)"splot 're_"//trim(fname)//"' with pm3d"
  if(present(wlines))write(10,*)"rep 're_"//trim(pname)//"_withlines' with lines"
  write(10,*)"#set term png size 1920,1280"
  write(10,*)"#set out 're_"//adjustl(trim(fname))//".png'"
  write(10,*)"#rep"
  write(10,"(A)")"EOF"
  !  
  write(10,*)"gnuplot -persist << EOF"
  write(10,*)"set term wxt"
  write(10,*)"set title 'Im_"//trim(fname)//"'"
  write(10,*)"set nokey"
  write(10,*)"set grid"
  write(10,*)"set view 50,10,1,1"
  write(10,*)"splot 'im_"//trim(fname)//"' with pm3d"
  if(present(wlines))write(10,*)"rep 'im_"//trim(pname)//"_withlines' with lines"
  write(10,*)"#set term png size 1920,1280"
  write(10,*)"#set out 'im_"//adjustl(trim(fname))//".png'"
  write(10,*)"#rep"
  write(10,"(A)")"EOF"
  close(10)
  call system("chmod +x "//adjustl(trim(pname))//".gp")
end subroutine splot3D_






subroutine splot3D__(pname,X1,X2,Y,wlines,nlines)
  integer                              :: i,j,Nx1,count,Nl,Nk
  character(len=*)                     :: pname
  real(8),dimension(:)                 :: X1 !(0:Nt)
  real(8),dimension(:,:)               :: X2 !(0:Nt,Lk)
  real(8),dimension(size(X2,1),size(X2,2)) :: Y  !(0:Nt,Lk)
  integer,optional                     :: nlines
  logical,optional                     :: wlines
  real(8)                              :: X1min,X1max
  real(8)                              :: X2min,X2max
  character(len=12)                   :: minx,miny,maxx,maxy
  character(len=256)                   :: fname,dname
  fname=get_filename(pname)
  dname=get_filepath(pname)
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
  X1min=minval(X1);X1max=maxval(X1)
  X2min=minval(X2);X2max=maxval(X2)
  write(minx,"(f12.4)")X1min
  write(maxx,"(f12.4)")X1max
  write(miny,"(f12.4)")X2min
  write(maxy,"(f12.4)")X2max
  call system("echo set term wxt 0 persist > "//adjustl(trim(pname))//".gp" )
  call system("echo set title \'"//trim(fname)//"\'  >> "//adjustl(trim(pname))//".gp" )
  call system("echo set pm3d map >> "//adjustl(trim(pname))//".gp" )
  call system("echo set size square >> "//adjustl(trim(pname))//".gp" )
  call system("echo set xrange ["//trim(adjustl(trim(minx)))//":"//trim(adjustl(trim(maxx)))//"]>> "//adjustl(trim(pname))//".gp" )
  call system("echo set yrange ["//trim(adjustl(trim(miny)))//":"//trim(adjustl(trim(maxy)))//"]>> "//adjustl(trim(pname))//".gp" )
  call system("echo splot \'"//trim(fname)//"\'  >> "//adjustl(trim(pname))//".gp" )
  call system("echo '#'set term png size 1920,1280 >> "//adjustl(trim(pname))//".gp" )
  call system("echo '#'set out \'"//adjustl(trim(fname))//".png\' >> "//adjustl(trim(pname))//".gp" )
  call system("echo '#'rep >> "//adjustl(trim(pname))//".gp" )
  !
  call system("echo '#' >> "//adjustl(trim(pname))//".gp" )
  call system("echo reset >> "//adjustl(trim(pname))//".gp" )
  !
  call system("echo set term wxt 1 persist >> "//adjustl(trim(pname))//".gp" )
  call system("echo set title \'"//trim(fname)//"\'  >> "//adjustl(trim(pname))//".gp" )
  call system("echo set nokey >> "//adjustl(trim(pname))//".gp" )
  call system("echo set grid >> "//adjustl(trim(pname))//".gp" )
  call system("echo set view 50,10,1,1 >> "//adjustl(trim(pname))//".gp" )
  call system("echo splot \'"//trim(fname)//"\' with pm3d >> "//adjustl(trim(pname))//".gp" )
  if(present(wlines))call system("echo rep \'"//trim(pname)//"_withlines\' with lines >> "//adjustl(trim(pname))//".gp" )
  call system("echo '#'set term png size 1920,1280 >> "//adjustl(trim(pname))//".gp" )
  call system("echo '#'set out \'"//adjustl(trim(fname))//".png\' >> "//adjustl(trim(pname))//".gp" )
  call system("echo '#'rep >> "//adjustl(trim(pname))//".gp" )
end subroutine splot3D__
!********************************************************************
!********************************************************************
!********************************************************************
