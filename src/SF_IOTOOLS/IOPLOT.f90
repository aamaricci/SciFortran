module IOPLOT
  USE IOFILE
  implicit none
  private

  interface splot
  ! SPLOT arrays (1--7)
     module procedure :: splotA1_RR
     module procedure :: splotA1_RC
     module procedure :: splotA2_RR
     module procedure :: splotA2_RC
     module procedure :: splotA3_RR
     module procedure :: splotA3_RC
     module procedure :: splotA4_RR
     module procedure :: splotA4_RC
     module procedure :: splotA5_RR
     module procedure :: splotA5_RC
     module procedure :: splotA6_RR
     module procedure :: splotA6_RC
     module procedure :: splotA7_RR
     module procedure :: splotA7_RC
  end interface splot


  interface splot3d
  ! SPLOT 3D:
     module procedure :: d_splot3d
     module procedure :: c_splot3d
     module procedure :: d_splot3d_animate
     module procedure :: c_splot3d_animate
  end interface splot3d

  interface save_array
  ! SAVE_ARRAY arrays (1--7)
     module procedure :: data_saveA1_R
     module procedure :: data_saveA1_C
     module procedure :: data_saveA2_R
     module procedure :: data_saveA2_C
     module procedure :: data_saveA3_R
     module procedure :: data_saveA3_C
     module procedure :: data_saveA4_R
     module procedure :: data_saveA4_C
     module procedure :: data_saveA5_R
     module procedure :: data_saveA5_C
     module procedure :: data_saveA6_R
     module procedure :: data_saveA6_C
     module procedure :: data_saveA7_R
     module procedure :: data_saveA7_C
  end interface save_array




  public :: splot
  public :: splot3d
  public :: save_array

  integer            :: unit
  character(len=128) :: fmt


contains


  ! SPLOT arrays (1--7)
subroutine splotA1_RR(pname,X,Y1,append)
   integer                       :: i,Np
   character(len=*)              :: pname
   real(8),dimension(:)          :: X
   real(8),dimension(size(X))    :: Y1
   logical,optional              :: append
   logical                       :: check,append_
   append_=.false.;if(present(append))append_=append
   if(append_)then
      inquire(file=reg(pname),exist=check)
      open(free_unit(unit),file=reg(pname),position="append")
      if(check)write(unit,*)
   else
      open(free_unit(unit),file=reg(pname))
   endif
   Np=size(X)
   do i=1,Np
      write(unit,*)X(i),Y1(i)
   enddo
   close(unit)
 end subroutine splotA1_RR
 subroutine splotA1_RC(pname,X,Y1,append)
   integer                       :: i,j,Np
   character(len=*)              :: pname
   real(8),dimension(:)          :: X
   complex(8),dimension(size(X)) :: Y1
   logical,optional              :: append
   logical                       :: check,append_
   append_=.false.;if(present(append))append_=append
   if(append_)then
      inquire(file=reg(pname),exist=check)
      open(free_unit(unit),file=reg(pname),position="append")
      if(check)write(unit,*)
   else
      open(free_unit(unit),file=reg(pname))
   endif
   Np=size(X)
   do i=1,Np
      write(unit,*)X(i),dimag(Y1(i)),dreal(Y1(i))
   enddo
   close(unit)
 end subroutine splotA1_RC
 !
 
 
 
 
 
 
 
 subroutine splotA2_RR(pname,X,Y1,append)
   integer                       :: i,j,Ny1,Ny2
   character(len=*)              :: pname
   real(8),dimension(:,:)        :: Y1
   real(8),dimension(size(Y1,2)) :: X
   logical,optional              :: append
   logical                       :: check,append_
   append_=.false.;if(present(append))append_=append
   if(append_)then
      inquire(file=reg(pname),exist=check)
      open(free_unit(unit),file=reg(pname),position="append")
      if(check)write(unit,*)
   else
      open(free_unit(unit),file=reg(pname))
   endif
   !
   Ny1=size(Y1,1)
   Ny2=size(Y1,2)
   !
   do i=1,Ny1
      do j=1,Ny2
         write(unit,*)X(j),Y1(i,j)
      enddo
      write(unit,*)""
   enddo
   close(unit)
 end subroutine splotA2_RR
 !
 subroutine splotA2_RC(pname,X,Y1,append)
   integer                       :: i,j,Ny1,Ny2
   character(len=*)              :: pname
   complex(8),dimension(:,:)     :: Y1
   real(8),dimension(size(Y1,2)) :: X
   logical,optional              :: append
   logical                       :: check,append_
   append_=.false.;if(present(append))append_=append
   if(append_)then
      inquire(file=reg(pname),exist=check)
      open(free_unit(unit),file=reg(pname),position="append")
      if(check)write(unit,*)
   else
      open(free_unit(unit),file=reg(pname))
   endif
   !
   Ny1=size(Y1,1)
   Ny2=size(Y1,2)
   !
   do i=1,Ny1
      do j=1,Ny2
         write(unit,*)X(j),dimag(Y1(i,j)),dreal(Y1(i,j))
      enddo
      write(unit,*)""
   enddo
   close(unit)
 end subroutine splotA2_RC
 
 
 
 
 
 subroutine splotA3_RR(pname,X,Y1,append)
   integer                       :: i,j,k,Ny1,Ny2,Ny3
   character(len=*)              :: pname
   real(8),dimension(:,:,:)      :: Y1
   real(8),dimension(size(Y1,3)) :: X
   logical,optional              :: append
   logical                       :: check,append_
   append_=.false.;if(present(append))append_=append
   if(append_)then
      inquire(file=reg(pname),exist=check)
      open(free_unit(unit),file=reg(pname),position="append")
      if(check)write(unit,*)
   else
      open(free_unit(unit),file=reg(pname))
   endif
   !
   Ny1=size(Y1,1)
   Ny2=size(Y1,2)
   Ny3=size(Y1,3)
   !
   do i=1,Ny1
      do j=1,Ny2
         do k=1,Ny3
            write(unit,*)X(k),Y1(i,j,k)
         enddo
         write(unit,*)""
      enddo
   enddo
   close(unit)
 end subroutine splotA3_RR
 !
 subroutine splotA3_RC(pname,X,Y1,append)
   integer                       :: i,j,k,Ny1,Ny2,Ny3
   character(len=*)              :: pname
   complex(8),dimension(:,:,:)   :: Y1
   real(8),dimension(size(Y1,3)) :: X
   logical,optional              :: append
   logical                       :: check,append_
   append_=.false.;if(present(append))append_=append
   if(append_)then
      inquire(file=reg(pname),exist=check)
      open(free_unit(unit),file=reg(pname),position="append")
      if(check)write(unit,*)
   else
      open(free_unit(unit),file=reg(pname))
   endif
   !
   Ny1=size(Y1,1)
   Ny2=size(Y1,2)
   Ny3=size(Y1,3)
   !
   do i=1,Ny1
      do j=1,Ny2
         do k=1,Ny3
            write(unit,*)X(k),dimag(Y1(i,j,k)),dreal(Y1(i,j,k))
         enddo
         write(unit,*)""
      enddo
   enddo
   close(unit)
 end subroutine splotA3_RC
 
 
 
 
 
 subroutine splotA4_RR(pname,X,Y1,append)
   integer                       :: Ny1,Ny2,Ny3,Ny4
   integer                       :: i1,i2,i3,i4
   character(len=*)              :: pname
   real(8),dimension(:,:,:,:)    :: Y1
   real(8),dimension(size(Y1,4)) :: X
   logical,optional              :: append
   logical                       :: check,append_
   append_=.false.;if(present(append))append_=append
   if(append_)then
      inquire(file=reg(pname),exist=check)
      open(free_unit(unit),file=reg(pname),position="append")
      if(check)write(unit,*)
   else
      open(free_unit(unit),file=reg(pname))
   endif
   !
   Ny1=size(Y1,1)
   Ny2=size(Y1,2)
   Ny3=size(Y1,3)
   Ny4=size(Y1,4)
   !
   do i1=1,Ny1
      do i2=1,Ny2
         do i3=1,Ny3
            do i4=1,Ny4
               write(unit,*)X(i4),Y1(i1,i2,i3,i4)
            enddo
            write(unit,*)""
         enddo
      enddo
   enddo
   close(unit)
 end subroutine splotA4_RR
 !
 subroutine splotA4_RC(pname,X,Y1,append)
   integer                       :: Ny1,Ny2,Ny3,Ny4
   integer                       :: i1,i2,i3,i4
   character(len=*)              :: pname
   complex(8),dimension(:,:,:,:) :: Y1
   real(8),dimension(size(Y1,4)) :: X
   logical,optional              :: append
   logical                       :: check,append_
   append_=.false.;if(present(append))append_=append
   if(append_)then
      inquire(file=reg(pname),exist=check)
      open(free_unit(unit),file=reg(pname),position="append")
      if(check)write(unit,*)
   else
      open(free_unit(unit),file=reg(pname))
   endif
   !
   !
   Ny1=size(Y1,1)
   Ny2=size(Y1,2)
   Ny3=size(Y1,3)
   Ny4=size(Y1,4)
   !
   do i1=1,Ny1
      do i2=1,Ny2
         do i3=1,Ny3
            do i4=1,Ny4
               write(unit,*)X(i4),dimag(Y1(i1,i2,i3,i4)),dreal(Y1(i1,i2,i3,i4))
            enddo
            write(unit,*)""
         enddo
      enddo
   enddo
   close(unit)
 end subroutine splotA4_RC
 
 
 
 
 
 
 
 subroutine splotA5_RR(pname,X,Y1,append)
   integer                         :: Ny1,Ny2,Ny3,Ny4,Ny5
   integer                         :: i1,i2,i3,i4,i5
   character(len=*)                :: pname
   real(8),dimension(:,:,:,:,:)    :: Y1
   real(8),dimension(size(Y1,5))   :: X
   logical,optional                :: append
   logical                         :: check,append_
   append_=.false.;if(present(append))append_=append
   if(append_)then
      inquire(file=reg(pname),exist=check)
      open(free_unit(unit),file=reg(pname),position="append")
      if(check)write(unit,*)
   else
      open(free_unit(unit),file=reg(pname))
   endif
   !
   Ny1=size(Y1,1)
   Ny2=size(Y1,2)
   Ny3=size(Y1,3)
   Ny4=size(Y1,4)
   Ny5=size(Y1,5)
   !
   do i1=1,Ny1
      do i2=1,Ny2
         do i3=1,Ny3
            do i4=1,Ny4
               do i5=1,Ny5
                  write(unit,*)X(i5),Y1(i1,i2,i3,i4,i5)
               enddo
               write(unit,*)""
            enddo
         enddo
      enddo
   enddo
   close(unit)
 end subroutine splotA5_RR
 !
 subroutine splotA5_RC(pname,X,Y1,append)
   integer                         :: Ny1,Ny2,Ny3,Ny4,Ny5
   integer                         :: i1,i2,i3,i4,i5
   character(len=*)                :: pname
   complex(8),dimension(:,:,:,:,:) :: Y1
   real(8),dimension(size(Y1,5))   :: X
   logical,optional                :: append
   logical                         :: check,append_
   append_=.false.;if(present(append))append_=append
   if(append_)then
      inquire(file=reg(pname),exist=check)
      open(free_unit(unit),file=reg(pname),position="append")
      if(check)write(unit,*)
   else
      open(free_unit(unit),file=reg(pname))
   endif
   !
   !
   Ny1=size(Y1,1)
   Ny2=size(Y1,2)
   Ny3=size(Y1,3)
   Ny4=size(Y1,4)
   Ny5=size(Y1,5)
   !
   do i1=1,Ny1
      do i2=1,Ny2
         do i3=1,Ny3
            do i4=1,Ny4
               do i5=1,Ny5
                  write(unit,*)X(i5),dimag(Y1(i1,i2,i3,i4,i5)),dreal(Y1(i1,i2,i3,i4,i5))
               enddo
               write(unit,*)""
            enddo
         enddo
      enddo
   enddo
   close(unit)
 end subroutine splotA5_RC
 
 
 
 
 
 
 
 subroutine splotA6_RR(pname,X,Y1,append)
   integer                           :: Ny1,Ny2,Ny3,Ny4,Ny5,Ny6
   integer                           :: i1,i2,i3,i4,i5,i6
   character(len=*)                  :: pname
   real(8),dimension(:,:,:,:,:,:)    :: Y1
   real(8),dimension(size(Y1,6))     :: X
   logical,optional                  :: append
   logical                           :: check,append_
   append_=.false.;if(present(append))append_=append
   if(append_)then
      inquire(file=reg(pname),exist=check)
      open(free_unit(unit),file=reg(pname),position="append")
      if(check)write(unit,*)
   else
      open(free_unit(unit),file=reg(pname))
   endif
   !
   Ny1=size(Y1,1)
   Ny2=size(Y1,2)
   Ny3=size(Y1,3)
   Ny4=size(Y1,4)
   Ny5=size(Y1,5)
   Ny6=size(Y1,6)
   !
   do i1=1,Ny1
      do i2=1,Ny2
         do i3=1,Ny3
            do i4=1,Ny4
               do i5=1,Ny5
                  do i6=1,Ny6
                     write(unit,*)X(i6),Y1(i1,i2,i3,i4,i5,i6)
                  enddo
                  write(unit,*)""
               enddo
            enddo
         enddo
      enddo
   enddo
   close(unit)
 end subroutine splotA6_RR
 !
 subroutine splotA6_RC(pname,X,Y1,append)
   integer                           :: Ny1,Ny2,Ny3,Ny4,Ny5,Ny6
   integer                           :: i1,i2,i3,i4,i5,i6
   character(len=*)                  :: pname
   complex(8),dimension(:,:,:,:,:,:) :: Y1
   real(8),dimension(size(Y1,6))     :: X
   logical,optional                  :: append
   logical                           :: check,append_
   append_=.false.;if(present(append))append_=append
   if(append_)then
      inquire(file=reg(pname),exist=check)
      open(free_unit(unit),file=reg(pname),position="append")
      if(check)write(unit,*)
   else
      open(free_unit(unit),file=reg(pname))
   endif
   !
   !
   Ny1=size(Y1,1)
   Ny2=size(Y1,2)
   Ny3=size(Y1,3)
   Ny4=size(Y1,4)
   Ny5=size(Y1,5)
   Ny6=size(Y1,6)
   !
   do i1=1,Ny1
      do i2=1,Ny2
         do i3=1,Ny3
            do i4=1,Ny4
               do i5=1,Ny5
                  do i6=1,Ny6
                     write(unit,*)X(i6),dimag(Y1(i1,i2,i3,i4,i5,i6)),dreal(Y1(i1,i2,i3,i4,i5,i6))
                  enddo
                  write(unit,*)""
               enddo
            enddo
         enddo
      enddo
   enddo
   close(unit)
 end subroutine splotA6_RC
 
 
 
 
 
 
 
 subroutine splotA7_RR(pname,X,Y1,append)
   integer                             :: Ny1,Ny2,Ny3,Ny4,Ny5,Ny6,Ny7
   integer                             :: i1,i2,i3,i4,i5,i6,i7
   character(len=*)                    :: pname
   real(8),dimension(:,:,:,:,:,:,:)    :: Y1
   real(8),dimension(size(Y1,7))       :: X
   logical,optional                    :: append
   logical                             :: check,append_
   append_=.false.;if(present(append))append_=append
   if(append_)then
      inquire(file=reg(pname),exist=check)
      open(free_unit(unit),file=reg(pname),position="append")
      if(check)write(unit,*)
   else
      open(free_unit(unit),file=reg(pname))
   endif
   !
   Ny1=size(Y1,1)
   Ny2=size(Y1,2)
   Ny3=size(Y1,3)
   Ny4=size(Y1,4)
   Ny5=size(Y1,5)
   Ny6=size(Y1,6)
   Ny7=size(Y1,7)
   !
   do i1=1,Ny1
      do i2=1,Ny2
         do i3=1,Ny3
            do i4=1,Ny4
               do i5=1,Ny5
                  do i6=1,Ny6
                     do i7=1,Ny7
                        write(unit,*)X(i7),Y1(i1,i2,i3,i4,i5,i6,i7)
                     enddo
                     write(unit,*)""
                  enddo
               enddo
            enddo
         enddo
      enddo
   enddo
   close(unit)
 end subroutine splotA7_RR
 !
 subroutine splotA7_RC(pname,X,Y1,append)
   integer                             :: Ny1,Ny2,Ny3,Ny4,Ny5,Ny6,Ny7
   integer                             :: i1,i2,i3,i4,i5,i6,i7
   character(len=*)                    :: pname
   complex(8),dimension(:,:,:,:,:,:,:) :: Y1
   real(8),dimension(size(Y1,7))       :: X
   logical,optional                    :: append
   logical                             :: check,append_
   append_=.false.;if(present(append))append_=append
   if(append_)then
      inquire(file=reg(pname),exist=check)
      open(free_unit(unit),file=reg(pname),position="append")
      if(check)write(unit,*)
   else
      open(free_unit(unit),file=reg(pname))
   endif
   !
   !
   Ny1=size(Y1,1)
   Ny2=size(Y1,2)
   Ny3=size(Y1,3)
   Ny4=size(Y1,4)
   Ny5=size(Y1,5)
   Ny6=size(Y1,6)
   Ny7=size(Y1,7)
   !
   do i1=1,Ny1
      do i2=1,Ny2
         do i3=1,Ny3
            do i4=1,Ny4
               do i5=1,Ny5
                  do i6=1,Ny6
                     do i7=1,Ny7
                        write(unit,*)X(i7),dimag(Y1(i1,i2,i3,i4,i5,i6,i7)),dreal(Y1(i1,i2,i3,i4,i5,i6,i7))
                     enddo
                     write(unit,*)""
                  enddo
               enddo
            enddo
         enddo
      enddo
   enddo
   close(unit)
 end subroutine splotA7_RC
 


  ! SPLOT 3D:
 subroutine d_splot3D(pname,X1,X2,Y,xmin,xmax,ymin,ymax,nosurface,wlines,nlines)
   integer                              :: i,j,Nx1,Nx2,count,Nl
   character(len=*)                     :: pname
   real(8),dimension(:)                 :: X1
   real(8),dimension(:)                 :: X2
   real(8),dimension(size(X1),size(X2)) :: Y
   real(8),optional                     :: xmin,xmax,ymin,ymax
   integer,optional                     :: nlines
   logical,optional                     :: wlines,nosurface
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
   X1min=minval(X1);if(present(xmin))X1min=xmin
   X1max=maxval(X1);if(present(xmax))X1max=xmax
   X2min=minval(X2);if(present(ymin))X2min=ymin
   X2max=maxval(X2);if(present(ymax))X2max=ymax
   write(minx,"(f12.4)")X1min
   write(maxx,"(f12.4)")X1max
   write(miny,"(f12.4)")X2min
   write(maxy,"(f12.4)")X2max
   if(present(wlines))close(720)
 
   open(10,file=adjustl(trim(pname))//"_map.gp")
   
   write(10,*)"set title '"//trim(fname)//"'"
   write(10,*)"set pm3d map"
   write(10,*)"set size square"
   write(10,*)"set xrange ["//trim(adjustl(trim(minx)))//":"//trim(adjustl(trim(maxx)))//"]"
   write(10,*)"set yrange ["//trim(adjustl(trim(miny)))//":"//trim(adjustl(trim(maxy)))//"]"
   write(10,*)"splot '"//trim(fname)//"'"
   write(10,*)"#set term png size 1920,1280"
   write(10,*)"#set out '"//adjustl(trim(fname))//".png'"
   write(10,*)"#rep"
   !
   close(10)
   if(present(nosurface).AND.nosurface)return
   open(10,file=adjustl(trim(pname))//"_surface.gp")
   
   write(10,*)"set title '"//trim(fname)//"'"
   write(10,*)"unset key"
   write(10,*)"unset grid"
   write(10,*)"set view 50,10,1,1"
   write(10,*)"splot '"//trim(fname)//"' with pm3d"
   if(present(wlines))write(10,*)"rep '"//trim(pname)//"_withlines' with lines"
   write(10,*)"#set term png size 1920,1280"
   write(10,*)"#set out '"//adjustl(trim(fname))//".png'"
   write(10,*)"#rep"
   close(10)
 end subroutine d_splot3D
 
 
 subroutine c_splot3D(pname,X1,X2,Y,xmin,xmax,ymin,ymax,nosurface,wlines,nlines)
   integer                                 :: i,j,Nx1,Nx2,count,Nl
   character(len=*)                        :: pname
   real(8),dimension(:)                    :: X1
   real(8),dimension(:)                    :: X2
   complex(8),dimension(size(X1),size(X2)) :: Y
   real(8),optional                        :: xmin,xmax,ymin,ymax
   integer,optional                        :: nlines
   logical,optional                        :: wlines,nosurface
   real(8)                                 :: X1min,X1max
   real(8)                                 :: X2min,X2max
   character(len=12)                       :: minx,miny,maxx,maxy
   character(len=256)                      :: fname,dname
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
         write(619,*)X1(i),X2(j),dreal(Y(i,j))
         write(719,*)X1(i),X2(j),dimag(Y(i,j))
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
   X1min=minval(X1);if(present(xmin))X1min=xmin
   X1max=maxval(X1);if(present(xmax))X1max=xmax
   X2min=minval(X2);if(present(ymin))X2min=ymin
   X2max=maxval(X2);if(present(ymax))X2max=ymax
   write(minx,"(f12.4)")X1min
   write(maxx,"(f12.4)")X1max
   write(miny,"(f12.4)")X2min
   write(maxy,"(f12.4)")X2max
   if(present(wlines))then
      close(620)
      close(720)
   endif
 
 
   !Re:
   open(10,file=adjustl(trim(pname))//"_re_map.gp")
   
   write(10,*)"set title 'Re_"//trim(fname)//"'"
   write(10,*)"set pm3d map"
   write(10,*)"set size square"
   write(10,*)"set xrange ["//trim(adjustl(trim(minx)))//":"//trim(adjustl(trim(maxx)))//"]"
   write(10,*)"set yrange ["//trim(adjustl(trim(miny)))//":"//trim(adjustl(trim(maxy)))//"]"
   write(10,*)"splot 're_"//trim(fname)//"'"
   write(10,*)"#set term png size 1920,1280"
   write(10,*)"#set out 're_"//adjustl(trim(fname))//".png'"
   write(10,*)"#rep"
   close(10)
   !Im
   open(10,file=adjustl(trim(pname))//"_im_map.gp")
   
   write(10,*)"set title 'Im_"//trim(fname)//"'"
   write(10,*)"set pm3d map"
   write(10,*)"set size square"
   write(10,*)"set xrange ["//trim(adjustl(trim(minx)))//":"//trim(adjustl(trim(maxx)))//"]"
   write(10,*)"set yrange ["//trim(adjustl(trim(miny)))//":"//trim(adjustl(trim(maxy)))//"]"
   write(10,*)"splot 'im_"//trim(fname)//"'"
   write(10,*)"#set term png size 1920,1280"
   write(10,*)"#set out 'im_"//adjustl(trim(fname))//".png'"
   write(10,*)"#rep"
   close(10)
   if(present(nosurface).AND.nosurface)return
   !Re
   open(10,file=adjustl(trim(pname))//"_re_surface.gp")
   
   write(10,*)"set title 'Re_"//trim(fname)//"'"
   write(10,*)"unset key"
   write(10,*)"unset grid"
   write(10,*)"set view 50,10,1,1"
   write(10,*)"splot 're_"//trim(fname)//"' with pm3d"
   if(present(wlines))write(10,*)"rep 're_"//trim(pname)//"_withlines' with lines"
   write(10,*)"#set term png size 1920,1280"
   write(10,*)"#set out 're_"//adjustl(trim(fname))//".png'"
   write(10,*)"#rep"
   close(10)
   !Im
   open(10,file=adjustl(trim(pname))//"_im_surface.gp")
   
   write(10,*)"set title 'Im_"//trim(fname)//"'"
   write(10,*)"set nokey"
   write(10,*)"set grid"
   write(10,*)"set view 50,10,1,1"
   write(10,*)"splot 'im_"//trim(fname)//"' with pm3d"
   if(present(wlines))write(10,*)"rep 'im_"//trim(pname)//"_withlines' with lines"
   write(10,*)"#set term png size 1920,1280"
   write(10,*)"#set out 'im_"//adjustl(trim(fname))//".png'"
   write(10,*)"#rep"
   close(10)
 end subroutine c_splot3D
 
 
 
 subroutine d_splot3d_animate(pname,X1,X2,Y,xmin,xmax,ymin,ymax)
   integer                  :: i,j,m,Nx1,Nx2,Nt
   character(len=*)         :: pname
   real(8),dimension(:)     :: X1
   real(8),dimension(:)     :: X2
   real(8),dimension(:,:,:) :: Y
   real(8)                  :: X1min,X1max
   real(8)                  :: X2min,X2max
   real(8),optional         :: xmin,xmax,ymin,ymax
   real(8)                  :: Rmin(3),Rmax(3),Zmin,Zmax
   Character(len=256)       :: fname,dname
   fname=get_filename(pname)
   dname=get_filepath(pname)
   Nx1=size(X1) ; Nx2=size(X2)
   if(size(Y,1)/=Nx1) stop "Error Nx1"
   if(size(Y,2)/=Nx2) stop "Error Nx2"
   Nt=size(Y,3)
   write(*,*) "splot3d_animate "//str(fname)//" ("//str(Nx1)//","//str(Nx2)//","//str(Nt)//")"
   open(719,file=adjustl(trim(pname)))
   do m=1,Nt
      do i=1,Nx1
         do j=1,Nx2
            write(719,*)X1(i),X2(j),Y(i,j,m)
         enddo
         write(719,*)""
      enddo
   enddo
   close(719)
   X1min=minval(X1);if(present(xmin))X1min=xmin
   X1max=maxval(X1);if(present(xmax))X1max=xmax
   X2min=minval(X2);if(present(ymin))X2min=ymin
   X2max=maxval(X2);if(present(ymax))X2max=ymax
   Rmin=minval(Y);Rmax=maxval(Y)
   Zmin=minval(Rmin);Zmax=maxval(Rmax)
   open(10,file=adjustl(trim(pname))//"_map.gp")
   ! write(10,*)"gnuplot -persist << EOF"
   write(10,*)"reset"
   write(10,*)"#set term gif animate"
   write(10,*)"#set output '"//trim(fname)//".gif'"
   
   write(10,*)"set pm3d map"
   write(10,*)"set size square"
   write(10,*)"set xrange ["//trim(adjustl(trim(txtfy(X1min))))//":"//trim(adjustl(trim(txtfy(X1max))))//"]"
   write(10,*)"set yrange ["//trim(adjustl(trim(txtfy(X2min))))//":"//trim(adjustl(trim(txtfy(X2max))))//"]"
   write(10,*)"set cbrange ["//trim(adjustl(trim(txtfy(Zmin))))//":"//trim(adjustl(trim(txtfy(Zmax))))//"]"
   write(10,*)"n="//trim(txtfy(Nt))
   write(10,*)"do for [i=0:n-1]{"
   write(10,*)"set title sprintf('"//trim(fname)//"; i=%i',i+1)"
   write(10,*)"splot '"//trim(fname)//"' every :::("//trim(adjustl(trim(txtfy(Nx2))))//&
        "*i)::("//trim(adjustl(trim(txtfy(Nx2))))//"*i+"//trim(adjustl(trim(txtfy(Nx1-1))))//") title '' "
   write(10,*)"}"
   write(10,*)"#set output"
   ! write(10,"(A)")"EOF"
   close(10)
   ! call system("chmod +x "//adjustl(trim(pname))//"_map.gp")
 end subroutine d_splot3D_animate
 
 
 
 subroutine c_splot3d_animate(pname,X1,X2,Y,xmin,xmax,ymin,ymax)
   integer                     :: i,j,m,Nx1,Nx2,Nt
   character(len=*)            :: pname
   real(8),dimension(:)        :: X1
   real(8),dimension(:)        :: X2
   complex(8),dimension(:,:,:) :: Y
   real(8)                     :: X1min,X1max
   real(8)                     :: X2min,X2max
   real(8),optional            :: xmin,xmax,ymin,ymax
   real(8)                     :: Rmin(3),Rmax(3),Zmin,Zmax  
   character(len=256)          :: fname,dname
   fname=get_filename(pname)
   dname=get_filepath(pname)
   Nx1=size(X1) ; Nx2=size(X2)
   if(size(Y,1)/=Nx1) stop "Error Nx1"
   if(size(Y,2)/=Nx2) stop "Error Nx2"
   Nt=size(Y,3)
   write(*,*) "splot3d_animate "//str(fname)//" ("//str(Nx1)//","//str(Nx2)//","//str(Nt)//")"
   open(619,file=adjustl(trim(dname))//"re_"//adjustl(trim(fname)))
   open(719,file=adjustl(trim(dname))//"im_"//adjustl(trim(fname)))
   do m=1,Nt
      do i=1,Nx1
         do j=1,Nx2
            write(619,*)X1(i),X2(j),real(Y(i,j,m),8)
            write(719,*)X1(i),X2(j),aimag(Y(i,j,m))
         enddo
         write(619,*)""
         write(719,*)""
      enddo
   enddo
   close(619)
   close(719)
   X1min=minval(X1);if(present(xmin))X1min=xmin
   X1max=maxval(X1);if(present(xmax))X1max=xmax
   X2min=minval(X2);if(present(ymin))X2min=ymin
   X2max=maxval(X2);if(present(ymax))X2max=ymax
   Rmin=minval(dreal(Y));Rmax=maxval(dreal(Y))
   Zmin=minval(Rmin);Zmax=maxval(Rmax)
   open(10,file=adjustl(trim(pname))//"_re_map.gp")
   ! write(10,*)"gnuplot -persist << EOF"
   write(10,*)"reset"
   write(10,*)"#set term gif animate"
   write(10,*)"#set output 're_"//trim(fname)//".gif'"
 
   write(10,*)"set pm3d map"
   write(10,*)"set size square"
   write(10,*)"set xrange ["//trim(adjustl(trim(txtfy(X1min))))//":"//trim(adjustl(trim(txtfy(X1max))))//"]"
   write(10,*)"set yrange ["//trim(adjustl(trim(txtfy(X2min))))//":"//trim(adjustl(trim(txtfy(X2max))))//"]"
   write(10,*)"set cbrange ["//trim(adjustl(trim(txtfy(Zmin))))//":"//trim(adjustl(trim(txtfy(Zmax))))//"]"
   write(10,*)"n="//trim(txtfy(Nt))
   write(10,*)"do for [i=0:n-1]{"
   write(10,*)"set title sprintf('Re("//trim(fname)//"); i=%i',i+1)"
   write(10,*)"splot 're_"//trim(fname)//"' every :::("//trim(adjustl(trim(txtfy(Nx2))))//&
        "*i)::("//trim(adjustl(trim(txtfy(Nx2))))//"*i+"//trim(adjustl(trim(txtfy(Nx1-1))))//") title '' "
   write(10,*)"}"
   write(10,*)"#set output"
   ! write(10,"(A)")"EOF"
   close(10)
 
   Rmin=minval(dimag(Y));Rmax=maxval(dimag(Y))
   Zmin=minval(Rmin);Zmax=maxval(Rmax)
   open(10,file=adjustl(trim(pname))//"_im_map.gp")
   ! write(10,*)"gnuplot -persist << EOF"
   write(10,*)"reset"
   write(10,*)"#set term gif animate"
   write(10,*)"#set output 'im_"//trim(fname)//".gif'"
 
   write(10,*)"set pm3d map"
   write(10,*)"set size square"
   write(10,*)"set xrange ["//trim(adjustl(trim(txtfy(X1min))))//":"//trim(adjustl(trim(txtfy(X1max))))//"]"
   write(10,*)"set yrange ["//trim(adjustl(trim(txtfy(X2min))))//":"//trim(adjustl(trim(txtfy(X2max))))//"]"
   write(10,*)"set cbrange ["//trim(adjustl(trim(txtfy(Zmin))))//":"//trim(adjustl(trim(txtfy(Zmax))))//"]"
   write(10,*)"n="//trim(txtfy(Nt))
   write(10,*)"do for [i=0:n-1]{"
   write(10,*)"set title sprintf('Im("//trim(fname)//"); i=%i',i+1)"
   write(10,*)"splot 'im_"//trim(fname)//"' every :::("//trim(adjustl(trim(txtfy(Nx2))))//&
        "*i)::("//trim(adjustl(trim(txtfy(Nx2))))//"*i+"//trim(adjustl(trim(txtfy(Nx1-1))))//") title '' "
   write(10,*)"}"
   write(10,*)"#set output"
   ! write(10,"(A)")"EOF"
   close(10)
   ! call system("chmod +x "//adjustl(trim(pname))//".gp")
 end subroutine c_splot3d_animate
 
 
 ! subroutine splot3D__(pname,X1,X2,Y,wlines,nlines)
 !   integer                              :: i,j,Nx1,count,Nl,Nk
 !   character(len=*)                     :: pname
 !   real(8),dimension(:)                 :: X1 !(0:Nt)
 !   real(8),dimension(:,:)               :: X2 !(0:Nt,Lk)
 !   real(8),dimension(size(X2,1),size(X2,2)) :: Y  !(0:Nt,Lk)
 !   integer,optional                     :: nlines
 !   logical,optional                     :: wlines
 !   real(8)                              :: X1min,X1max
 !   real(8)                              :: X2min,X2max
 !   character(len=12)                   :: minx,miny,maxx,maxy
 !   character(len=256)                   :: fname,dname
 !   fname=get_filename(pname)
 !   dname=get_filepath(pname)
 !   Nx1=size(X1) ; Nk=size(X2,2)
 !   Nl=5; if(present(nlines))Nl=nlines
 !   open(719,file=adjustl(trim(pname)))
 !   if(present(wlines))open(720,file=adjustl(trim(pname))//"_withlines")
 !   do i=1,Nx1
 !      count=mod(i,Nl)
 !      do j=1,Nk
 !         write(719,*)X1(i),X2(i,j),Y(i,j)
 !         if(present(wlines).AND.count==0)write(720,*)X1(i),X2(i,j),Y(i,j)
 !      enddo
 !      write(719,*)""
 !      if(present(wlines))write(720,*)""
 !      if(present(wlines))write(720,*)""
 !   enddo
 !   close(719)
 !   if(present(wlines))close(720)
 !   X1min=minval(X1);X1max=maxval(X1)
 !   X2min=minval(X2);X2max=maxval(X2)
 !   write(minx,"(f12.4)")X1min
 !   write(maxx,"(f12.4)")X1max
 !   write(miny,"(f12.4)")X2min
 !   write(maxy,"(f12.4)")X2max
 !   call system("echo set term wxt 0 persist > "//adjustl(trim(pname))//".gp" )
 !   call system("echo set title \'"//trim(fname)//"\'  >> "//adjustl(trim(pname))//".gp" )
 !   call system("echo set pm3d map >> "//adjustl(trim(pname))//".gp" )
 !   call system("echo set size square >> "//adjustl(trim(pname))//".gp" )
 !   call system("echo set xrange ["//trim(adjustl(trim(minx)))//":"//trim(adjustl(trim(maxx)))//"]>> "//adjustl(trim(pname))//".gp" )
 !   call system("echo set yrange ["//trim(adjustl(trim(miny)))//":"//trim(adjustl(trim(maxy)))//"]>> "//adjustl(trim(pname))//".gp" )
 !   call system("echo splot \'"//trim(fname)//"\'  >> "//adjustl(trim(pname))//".gp" )
 !   call system("echo '#'set term png size 1920,1280 >> "//adjustl(trim(pname))//".gp" )
 !   call system("echo '#'set out \'"//adjustl(trim(fname))//".png\' >> "//adjustl(trim(pname))//".gp" )
 !   call system("echo '#'rep >> "//adjustl(trim(pname))//".gp" )
 !   !
 !   call system("echo '#' >> "//adjustl(trim(pname))//".gp" )
 !   call system("echo reset >> "//adjustl(trim(pname))//".gp" )
 !   !
 !   call system("echo set term wxt 1 persist >> "//adjustl(trim(pname))//".gp" )
 !   call system("echo set title \'"//trim(fname)//"\'  >> "//adjustl(trim(pname))//".gp" )
 !   call system("echo set nokey >> "//adjustl(trim(pname))//".gp" )
 !   call system("echo set grid >> "//adjustl(trim(pname))//".gp" )
 !   call system("echo set view 50,10,1,1 >> "//adjustl(trim(pname))//".gp" )
 !   call system("echo splot \'"//trim(fname)//"\' with pm3d >> "//adjustl(trim(pname))//".gp" )
 !   if(present(wlines))call system("echo rep \'"//trim(pname)//"_withlines\' with lines >> "//adjustl(trim(pname))//".gp" )
 !   call system("echo '#'set term png size 1920,1280 >> "//adjustl(trim(pname))//".gp" )
 !   call system("echo '#'set out \'"//adjustl(trim(fname))//".png\' >> "//adjustl(trim(pname))//".gp" )
 !   call system("echo '#'rep >> "//adjustl(trim(pname))//".gp" )
 ! end subroutine splot3D__
 ! !********************************************************************
 ! !********************************************************************
 ! !********************************************************************
 

  
  ! SAVE_ARRAY arrays (1--7)
 subroutine data_saveA1_R(pname,Y1)
   integer                              :: i,Np
   character(len=*)                     :: pname
   real(8),dimension(:)                 :: Y1
   Np=size(Y1)
   open(free_unit(unit),file=reg(pname))
   do i=1,Np
      write(unit,*)Y1(i)
   enddo
   close(unit)
   call file_gzip(reg(pname))
 end subroutine data_saveA1_R
 
 subroutine data_saveA1_C(pname,Y1)
   integer                              :: i,Np
   character(len=*)                     :: pname
   complex(8),dimension(:)              :: Y1
   Np=size(Y1)  
   open(free_unit(unit),file=reg(pname))
   do i=1,Np
      write(unit,*)Y1(i)
   enddo
   close(unit)
   call file_gzip(reg(pname))
 end subroutine data_saveA1_C
 
 
 
 
 
 subroutine data_saveA2_R(pname,Y1,order,wspace)
   integer                   :: i,j,Ny1,Ny2
   character(len=*)          :: pname
   real(8),dimension(:,:)    :: Y1
   character(len=*),optional :: order
   logical,optional          :: wspace
   character(len=1)          :: order_
   logical                   :: wspace_
   order_ = "R"   ; if(present(order))order_=trim(order(1:1))
   wspace_= .true.; if(present(wspace))wspace_=wspace
   Ny1=size(Y1,1) ; Ny2=size(Y1,2)
   open(free_unit(unit),file=reg(pname))
   select case(order_)
   case default
      stop "save_array: order != Row-major, Col-major"
   case ("R")
      do i=1,Ny1
         do j=1,Ny2
            write(unit,*)Y1(i,j)
         enddo
         if(wspace_)write(unit,*)
      enddo
   case ("C")
      do j=1,Ny2
         do i=1,Ny1
            write(unit,*)Y1(i,j)
         enddo
         if(wspace_)write(unit,*)
      enddo
   end select
   close(unit)
   call file_gzip(reg(pname))
 end subroutine data_saveA2_R
 
 subroutine data_saveA2_C(pname,Y1,order,wspace)
   integer                   :: i,j,Ny1,Ny2
   character(len=*)          :: pname
   complex(8),dimension(:,:) :: Y1
   character(len=*),optional :: order
   logical,optional          :: wspace
   character(len=1)          :: order_
   logical                   :: wspace_
   order_ = "R"   ; if(present(order))order_=trim(order(1:1))
   wspace_= .true.; if(present(wspace))wspace_=wspace
   Ny1=size(Y1,1) ; Ny2=size(Y1,2)
   open(free_unit(unit),file=reg(pname))
   select case(order_)
   case default
      stop "save_array: order != Row-major, Col-major"
   case ("R")
      do i=1,Ny1
         do j=1,Ny2
            write(unit,*)Y1(i,j)
         enddo
         if(wspace_)write(unit,*)
      enddo
   case ("C")
      do j=1,Ny2
         do i=1,Ny1
            write(unit,*)Y1(i,j)
         enddo
         if(wspace_)write(unit,*)
      enddo
   end select
   close(unit)
   call file_gzip(reg(pname))
 end subroutine data_saveA2_C
 
 
 
 
 
 
 subroutine data_saveA3_R(pname,Y1,order,wspace)
   integer                     :: Ny1,Ny2,Ny3
   integer                     :: i1,i2,i3
   character(len=*)            :: pname
   real(8),dimension(:,:,:)    :: Y1
   character(len=*),optional   :: order
   logical,optional            :: wspace
   character(len=1)            :: order_
   logical                     :: wspace_
   order_ = "R"   ; if(present(order))order_=trim(order(1:1))
   wspace_= .true.; if(present(wspace))wspace_=wspace
   !
   open(free_unit(unit),file=reg(pname))
   !
   Ny1=size(Y1,1)
   Ny2=size(Y1,2)
   Ny3=size(Y1,3)
   !
   select case(order_)
   case default
      stop "save_array: order != Row-major, Col-major"
   case ("R")
      do i1=1,Ny1
         do i2=1,Ny2
            do i3=1,Ny3
               write(unit,*)Y1(i1,i2,i3)
            enddo
            if(wspace_)write(unit,*)
         enddo
      enddo
   case ("C")
      do i3=1,Ny3
         do i2=1,Ny2
            do i1=1,Ny1
               write(unit,*)Y1(i1,i2,i3)
            enddo
            if(wspace_)write(unit,*)
         enddo
      enddo
   end select
   close(unit)
   call file_gzip(reg(pname))
 end subroutine data_saveA3_R
 
 subroutine data_saveA3_C(pname,Y1,order,wspace)
   integer                     :: Ny1,Ny2,Ny3
   integer                     :: i1,i2,i3
   character(len=*)            :: pname
   complex(8),dimension(:,:,:) :: Y1
   character(len=*),optional   :: order
   logical,optional            :: wspace
   character(len=1)            :: order_
   logical                     :: wspace_
   order_ = "R"   ; if(present(order))order_=trim(order(1:1))
   wspace_= .true.; if(present(wspace))wspace_=wspace
   !
   open(free_unit(unit),file=reg(pname))
   !
   Ny1=size(Y1,1)
   Ny2=size(Y1,2)
   Ny3=size(Y1,3)
   !
   select case(order_)
   case default
      stop "save_array: order != Row-major, Col-major"
   case ("R")
      do i1=1,Ny1
         do i2=1,Ny2
            do i3=1,Ny3
               write(unit,*)Y1(i1,i2,i3)
            enddo
            if(wspace_)write(unit,*)
         enddo
      enddo
   case ("C")
      do i3=1,Ny3
         do i2=1,Ny2
            do i1=1,Ny1
               write(unit,*)Y1(i1,i2,i3)
            enddo
            if(wspace_)write(unit,*)
         enddo
      enddo
   end select
   close(unit)
   call file_gzip(reg(pname))
 end subroutine data_saveA3_C
 
 
 
 
 
 
 
 
 
 subroutine data_saveA4_R(pname,Y1,order,wspace)
   integer                    :: Ny1,Ny2,Ny3,Ny4
   integer                    :: i1,i2,i3,i4
   character(len=*)           :: pname
   real(8),dimension(:,:,:,:) :: Y1
   character(len=*),optional  :: order
   logical,optional           :: wspace
   character(len=1)           :: order_
   logical                    :: wspace_
   order_ = "R"   ; if(present(order))order_=trim(order(1:1))
   wspace_= .true.; if(present(wspace))wspace_=wspace
   open(free_unit(unit),file=reg(pname))
   Ny1=size(Y1,1)
   Ny2=size(Y1,2)
   Ny3=size(Y1,3)
   Ny4=size(Y1,4)
   select case(order_)
   case default
      stop "save_array: order != Row-major, Col-major"
   case ("R")
      do i1=1,Ny1
         do i2=1,Ny2
            do i3=1,Ny3
               do i4=1,Ny4
                  write(unit,*)Y1(i1,i2,i3,i4)
               enddo
               if(wspace_)write(unit,*)
            enddo
         enddo
      enddo
   case ("C")
      do i4=1,Ny4
         do i3=1,Ny3
            do i2=1,Ny2
               do i1=1,Ny1
                  write(unit,*)Y1(i1,i2,i3,i4)
               enddo
               if(wspace_)write(unit,*)
            enddo
         enddo
      enddo
   end select
   close(unit)
   call file_gzip(reg(pname))
 end subroutine data_saveA4_R
 
 subroutine data_saveA4_C(pname,Y1,order,wspace)
   integer                       :: Ny1,Ny2,Ny3,Ny4
   integer                       :: i1,i2,i3,i4
   character(len=*)              :: pname
   complex(8),dimension(:,:,:,:) :: Y1
   character(len=*),optional     :: order
   logical,optional              :: wspace
   character(len=1)              :: order_
   logical                       :: wspace_
   order_ = "R"   ; if(present(order))order_=trim(order(1:1))
   wspace_= .true.; if(present(wspace))wspace_=wspace
 
 
   open(free_unit(unit),file=reg(pname))
   Ny1=size(Y1,1)
   Ny2=size(Y1,2)
   Ny3=size(Y1,3)
   Ny4=size(Y1,4)
   select case(order_)
   case default
      stop "save_array: order != Row-major, Col-major"
   case ("R")
      do i1=1,Ny1
         do i2=1,Ny2
            do i3=1,Ny3
               do i4=1,Ny4
                  write(unit,*)Y1(i1,i2,i3,i4)              
               enddo
               if(wspace_)write(unit,*)
            enddo
         enddo
      enddo
   case ("C")
      do i4=1,Ny4
         do i3=1,Ny3
            do i2=1,Ny2
               do i1=1,Ny1
                  write(unit,*)Y1(i1,i2,i3,i4)
               enddo
               if(wspace_)write(unit,*)
            enddo
         enddo
      enddo
   end select
   close(unit)
   call file_gzip(reg(pname))
 end subroutine data_saveA4_C
 
 
 
 
 
 
 
 
 
 subroutine data_saveA5_R(pname,Y1,order,wspace)
   integer                      :: Ny1,Ny2,Ny3,Ny4,Ny5
   integer                      :: i1,i2,i3,i4,i5
   character(len=*)             :: pname
   real(8),dimension(:,:,:,:,:) :: Y1
   character(len=*),optional    :: order
   logical,optional             :: wspace
   character(len=1)             :: order_
   logical                      :: wspace_
   order_ = "R"   ; if(present(order))order_=trim(order(1:1))
   wspace_= .true.; if(present(wspace))wspace_=wspace
   open(free_unit(unit),file=reg(pname))
   Ny1=size(Y1,1)
   Ny2=size(Y1,2)
   Ny3=size(Y1,3)
   Ny4=size(Y1,4)
   Ny5=size(Y1,5)
   select case(order_)
   case default
      stop "save_array: order != Row-major, Col-major"
   case ("R")
      do i1=1,Ny1
         do i2=1,Ny2
            do i3=1,Ny3
               do i4=1,Ny4
                  do i5=1,Ny5
                     write(unit,*)Y1(i1,i2,i3,i4,i5)
                  enddo
                  if(wspace_)write(unit,*)
               enddo
            enddo
         enddo
      enddo
   case ("C")
      do i5=1,Ny5
         do i4=1,Ny4
            do i3=1,Ny3
               do i2=1,Ny2
                  do i1=1,Ny1
                     write(unit,*)Y1(i1,i2,i3,i4,i5)
                  enddo
                  if(wspace_)write(unit,*)
               enddo
            enddo
         enddo
      enddo
   end select
   close(unit)
   call file_gzip(reg(pname))
 end subroutine data_saveA5_R
 
 subroutine data_saveA5_C(pname,Y1,order,wspace)
   integer                         :: Ny1,Ny2,Ny3,Ny4,Ny5
   integer                         :: i1,i2,i3,i4,i5
   character(len=*)                :: pname
   complex(8),dimension(:,:,:,:,:) :: Y1
   character(len=*),optional       :: order
   logical,optional                :: wspace
   character(len=1)                :: order_
   logical                         :: wspace_
   order_ = "R"   ; if(present(order))order_=trim(order(1:1))
   wspace_= .true.; if(present(wspace))wspace_=wspace
   open(free_unit(unit),file=reg(pname))
   Ny1=size(Y1,1)
   Ny2=size(Y1,2)
   Ny3=size(Y1,3)
   Ny4=size(Y1,4)
   Ny5=size(Y1,5)
   select case(order_)
   case default
      stop "save_array: order != Row-major, Col-major"
   case ("R")
      do i1=1,Ny1
         do i2=1,Ny2
            do i3=1,Ny3
               do i4=1,Ny4
                  do i5=1,Ny5
                     write(unit,*)Y1(i1,i2,i3,i4,i5)
                  enddo
                  if(wspace_)write(unit,*)
               enddo
            enddo
         enddo
      enddo
   case ("C")
      do i5=1,Ny5
         do i4=1,Ny4
            do i3=1,Ny3
               do i2=1,Ny2
                  do i1=1,Ny1
                     write(unit,*)Y1(i1,i2,i3,i4,i5)
                  enddo
                  if(wspace_)write(unit,*)
               enddo
            enddo
         enddo
      enddo
   end select
   close(unit)
   call file_gzip(reg(pname))
 end subroutine data_saveA5_C
 
 
 
 
 
 
 
 
 
 
 
 subroutine data_saveA6_R(pname,Y1,order,wspace)
   integer                        :: Ny1,Ny2,Ny3,Ny4,Ny5,Ny6
   integer                        :: i1,i2,i3,i4,i5,i6
   character(len=*)               :: pname
   real(8),dimension(:,:,:,:,:,:) :: Y1
   character(len=*),optional      :: order
   logical,optional               :: wspace
   character(len=1)               :: order_
   logical                        :: wspace_
   order_ = "R"   ; if(present(order))order_=trim(order(1:1))
   wspace_= .true.; if(present(wspace))wspace_=wspace
   open(free_unit(unit),file=reg(pname))
   Ny1=size(Y1,1)
   Ny2=size(Y1,2)
   Ny3=size(Y1,3)
   Ny4=size(Y1,4)
   Ny5=size(Y1,5)
   Ny6=size(Y1,6)
   !
   select case(order_)
   case default
      stop "save_array: order != Row-major, Col-major"
   case ("R")
      do i1=1,Ny1
         do i2=1,Ny2
            do i3=1,Ny3
               do i4=1,Ny4
                  do i5=1,Ny5
                     do i6=1,Ny6
                        write(unit,*)Y1(i1,i2,i3,i4,i5,i6)
                     enddo
                     if(wspace_)write(unit,*)
                  enddo
               enddo
            enddo
         enddo
      enddo
   case ("C")
      do i6=1,Ny6
         do i5=1,Ny5
            do i4=1,Ny4
               do i3=1,Ny3
                  do i2=1,Ny2
                     do i1=1,Ny1
                        write(unit,*)Y1(i1,i2,i3,i4,i5,i6)
                     enddo
                     if(wspace_)write(unit,*)
                  enddo
               enddo
            enddo
         enddo
      enddo
   end select
   close(unit)
   call file_gzip(reg(pname))
 end subroutine data_saveA6_R
 
 subroutine data_saveA6_C(pname,Y1,order,wspace)
   integer                           :: Ny1,Ny2,Ny3,Ny4,Ny5,Ny6
   integer                           :: i1,i2,i3,i4,i5,i6
   character(len=*)                  :: pname
   complex(8),dimension(:,:,:,:,:,:) :: Y1
   character(len=*),optional         :: order
   logical,optional                  :: wspace
   character(len=1)                  :: order_
   logical                           :: wspace_
   order_ = "R"   ; if(present(order))order_=trim(order(1:1))
   wspace_= .true.; if(present(wspace))wspace_=wspace
   open(free_unit(unit),file=reg(pname))
   Ny1=size(Y1,1)
   Ny2=size(Y1,2)
   Ny3=size(Y1,3)
   Ny4=size(Y1,4)
   Ny5=size(Y1,5)
   Ny6=size(Y1,6)
   !
   select case(order_)
   case default
      stop "save_array: order != Row-major, Col-major"
   case ("R")
      do i1=1,Ny1
         do i2=1,Ny2
            do i3=1,Ny3
               do i4=1,Ny4
                  do i5=1,Ny5
                     do i6=1,Ny6
                        write(unit,*)Y1(i1,i2,i3,i4,i5,i6)
                     enddo
                     if(wspace_)write(unit,*)
                  enddo
               enddo
            enddo
         enddo
      enddo
   case ("C")
      do i6=1,Ny6
         do i5=1,Ny5
            do i4=1,Ny4
               do i3=1,Ny3
                  do i2=1,Ny2
                     do i1=1,Ny1
                        write(unit,*)Y1(i1,i2,i3,i4,i5,i6)
                     enddo
                     if(wspace_)write(unit,*)
                  enddo
               enddo
            enddo
         enddo
      enddo
   end select
   close(unit)
   call file_gzip(reg(pname))
 end subroutine data_saveA6_C
 
 
 
 
 
 
 
 
 
 
 
 
 subroutine data_saveA7_R(pname,Y1,order,wspace)
   integer                          :: Ny1,Ny2,Ny3,Ny4,Ny5,Ny6,Ny7
   integer                          :: i1,i2,i3,i4,i5,i6,i7
   character(len=*)                 :: pname
   real(8),dimension(:,:,:,:,:,:,:) :: Y1
   character(len=*),optional        :: order
   logical,optional                 :: wspace
   character(len=1)                 :: order_
   logical                          :: wspace_
   order_ = "R"   ; if(present(order))order_=trim(order(1:1))
   wspace_= .true.; if(present(wspace))wspace_=wspace
   open(free_unit(unit),file=reg(pname))
   Ny1=size(Y1,1)
   Ny2=size(Y1,2)
   Ny3=size(Y1,3)
   Ny4=size(Y1,4)
   Ny5=size(Y1,5)
   Ny6=size(Y1,6)
   Ny7=size(Y1,7)
   !
   select case(order_)
   case default
      stop "save_array: order != Row-major, Col-major"
   case ("R")
      do i1=1,Ny1
         do i2=1,Ny2
            do i3=1,Ny3
               do i4=1,Ny4
                  do i5=1,Ny5
                     do i6=1,Ny6
                        do i7=1,Ny7
                           write(unit,*)Y1(i1,i2,i3,i4,i5,i6,i7)
                        enddo
                        if(wspace_)write(unit,*)
                     enddo
                  enddo
               enddo
            enddo
         enddo
      enddo
   case ("C")
      do i7=1,Ny7
         do i6=1,Ny6
            do i5=1,Ny5
               do i4=1,Ny4
                  do i3=1,Ny3
                     do i2=1,Ny2
                        do i1=1,Ny1
                           write(unit,*)Y1(i1,i2,i3,i4,i5,i6,i7)
                        enddo
                        if(wspace_)write(unit,*)
                     enddo
                  enddo
               enddo
            enddo
         enddo
      enddo
   end select
   close(unit)
   call file_gzip(reg(pname))
 end subroutine data_saveA7_R
 
 subroutine data_saveA7_C(pname,Y1,order,wspace)
   integer                             :: Ny1,Ny2,Ny3,Ny4,Ny5,Ny6,Ny7
   integer                             :: i1,i2,i3,i4,i5,i6,i7
   character(len=*)                    :: pname
   complex(8),dimension(:,:,:,:,:,:,:) :: Y1
   character(len=*),optional           :: order
   logical,optional                    :: wspace
   character(len=1)                    :: order_
   logical                             :: wspace_
   order_ = "R"   ; if(present(order))order_=trim(order(1:1))
   wspace_= .true.; if(present(wspace))wspace_=wspace
   open(free_unit(unit),file=reg(pname))
   Ny1=size(Y1,1)
   Ny2=size(Y1,2)
   Ny3=size(Y1,3)
   Ny4=size(Y1,4)
   Ny5=size(Y1,5)
   Ny6=size(Y1,6)
   Ny7=size(Y1,7)
   !
   select case(order_)
   case default
      stop "save_array: order != Row-major, Col-major"
   case ("R")
      do i1=1,Ny1
         do i2=1,Ny2
            do i3=1,Ny3
               do i4=1,Ny4
                  do i5=1,Ny5
                     do i6=1,Ny6
                        do i7=1,Ny7
                           write(unit,*)Y1(i1,i2,i3,i4,i5,i6,i7)
                        enddo
                        if(wspace_)write(unit,*)
                     enddo
                  enddo
               enddo
            enddo
         enddo
      enddo
   case ("C")
      do i7=1,Ny7
         do i6=1,Ny6
            do i5=1,Ny5
               do i4=1,Ny4
                  do i3=1,Ny3
                     do i2=1,Ny2
                        do i1=1,Ny1
                           write(unit,*)Y1(i1,i2,i3,i4,i5,i6,i7)
                        enddo
                        if(wspace_)write(unit,*)
                     enddo
                  enddo
               enddo
            enddo
         enddo
      enddo
   end select
   close(unit)
   call file_gzip(reg(pname))
 end subroutine data_saveA7_C
 


end module IOPLOT
