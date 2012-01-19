!###############################################################
!     PROGRAM  : SLPLOT
!     TYPE     : Module
!     PURPOSE  : SIMPLE PLOTTING LIBRARY FOR FORTRAN 90/95
!     AUTHORS  : Adriano Amaricci (SISSA)
!###############################################################
module SLPLOT
  USE IOTOOLS
  implicit none
  private

  interface splot
     module procedure &
          splotP_II,splotP_IR,splotP_IC, &
          splotP_RI,splotP_RR,splotP_RC, &
          splotV_II, splotV_IR,splotV_IC,&
          splotV_RI, splotV_RR,splotV_RC,&
          data_saveV_I,&
          data_saveV_R,&
          data_saveV_C,&
          data_saveM_I,&
          data_saveM_R,&
          data_saveM_C,&
          data_saveA3_I,&
          data_saveA3_R,&
          data_saveA3_C,&
          splot3D,splot3D_,splot3D__
  end interface splot

  public :: splot

contains

  ! 0-dim array
  ! X=int  ; Y=int,dble,cmplx
  ! X=dble ; Y=int,dble,cmplx
  include "splot_P.f90"

  ! 1-dim array
  ! X=int  ; Y=int,dble,cmplx
  ! X=dble ; Y=int,dble,cmplx
  include "splot_V.f90"

  ! 1,2-dim arrays
  ! Y=int,dble,cmplx
  ! X=int,dble [only for 2-dim, optional]
  include "data_save.f90"


  subroutine splot3D(pname,X1,X2,Y,wlines,nlines)
    integer                              :: i,j,Ny,Nx1,Nx2,count,Nl
    character(len=*)                     :: pname
    real(8),dimension(:)                 :: X1
    real(8),dimension(:)                 :: X2
    real(8),dimension(size(X1),size(X2)) :: Y
    integer,optional                     :: nlines
    logical,optional                     :: wlines
    real(8)                              :: X1min,X1max
    real(8)                              :: X2min,X2max
    character(len=12)                    :: minx,miny,maxx,maxy
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
    call system("echo set term wxt 0 persist > plot_"//adjustl(trim(pname)) )
    call system("echo set title \'"//trim(pname)//"\'  >> plot_"//adjustl(trim(pname)) )
    call system("echo set pm3d map >> plot_"//adjustl(trim(pname)) )
    call system("echo set size square >> plot_"//adjustl(trim(pname)) )
    call system("echo set xrange ["//trim(adjustl(trim(minx)))//":"//trim(adjustl(trim(maxx)))//"]>> plot_"//adjustl(trim(pname)) )
    call system("echo set yrange ["//trim(adjustl(trim(miny)))//":"//trim(adjustl(trim(maxy)))//"]>> plot_"//adjustl(trim(pname)) )
    call system("echo splot \'"//trim(pname)//"\'  >> plot_"//adjustl(trim(pname)) )
    call system("echo '#'set term png size 1920,1280 >> plot_"//adjustl(trim(pname)) )
    call system("echo '#'set out\'"//adjustl(trim(pname))//".png\' >> plot_"//adjustl(trim(pname)) )
    call system("echo '#'rep >> plot_"//adjustl(trim(pname)) )
    !
    call system("echo '#' >> plot_"//adjustl(trim(pname)) )
    call system("echo reset >> plot_"//adjustl(trim(pname)) )
    !
    call system("echo set term wxt 1 persist >> plot_"//adjustl(trim(pname)) )
    call system("echo set title \'"//trim(pname)//"\'  >> plot_"//adjustl(trim(pname)) )
    call system("echo set nokey >> plot_"//adjustl(trim(pname)) )
    call system("echo set grid >> plot_"//adjustl(trim(pname)) )
    call system("echo set view 50,10,1,1 >> plot_"//adjustl(trim(pname)) )
    call system("echo splot \'"//trim(pname)//"\' with pm3d >> plot_"//adjustl(trim(pname)) )
    if(present(wlines))call system("echo rep \'"//trim(pname)//"_withlines\' with lines >> plot_"//adjustl(trim(pname)) )
    call system("echo '#'set term png size 1920,1280 >> plot_"//adjustl(trim(pname)) )
    call system("echo '#'set out\'"//adjustl(trim(pname))//".png\' >> plot_"//adjustl(trim(pname)) )
    call system("echo '#'rep >> plot_"//adjustl(trim(pname)) )
  end subroutine splot3D


  subroutine splot3D_(pname,X1,X2,Y,wlines,nlines)
    integer                              :: i,j,Ny,Nx1,Nx2,count,Nl
    character(len=*)                     :: pname
    real(8),dimension(:)                 :: X1
    real(8),dimension(:)                 :: X2
    complex(8),dimension(size(X1),size(X2)) :: Y
    integer,optional                     :: nlines
    logical,optional                     :: wlines
    real(8)                              :: X1min,X1max
    real(8)                              :: X2min,X2max
    character(len=6)                     :: minx,miny,maxx,maxy
    if(present(nlines).and.present(wlines))then
       call splot3d("re"//pname,X1,X2,real(Y,8),wlines,nlines)
       call splot3d("im"//pname,X1,X2,dimag(Y),wlines,nlines)
    elseif(present(wlines))then
       call splot3d("re"//pname,X1,X2,real(Y,8),wlines)
       call splot3d("im"//pname,X1,X2,dimag(Y),wlines)
    elseif(present(nlines))then
       call splot3d("re"//pname,X1,X2,real(Y,8),nlines=nlines)
       call splot3d("im"//pname,X1,X2,dimag(Y),nlines=nlines)
    else
       call splot3d("re"//pname,X1,X2,real(Y,8))
       call splot3d("im"//pname,X1,X2,dimag(Y))
    endif
  end subroutine splot3D_


  subroutine splot3D__(pname,X1,X2,Y,wlines,nlines)
    integer                              :: i,j,Ny,Nx1,Nx2,count,Nl,Nk
    character(len=*)                     :: pname
    real(8),dimension(:)                 :: X1 !(0:Nt)
    real(8),dimension(:,:)               :: X2 !(0:Nt,Lk)
    real(8),dimension(size(X2,1),size(X2,2)) :: Y  !(0:Nt,Lk)
    integer,optional                     :: nlines
    logical,optional                     :: wlines
    real(8)                              :: X1min,X1max
    real(8)                              :: X2min,X2max
    character(len=12)                   :: minx,miny,maxx,maxy
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
    call system("echo set term wxt 0 persist > plot_"//adjustl(trim(pname)) )
    call system("echo set title \'"//trim(pname)//"\'  >> plot_"//adjustl(trim(pname)) )
    call system("echo set pm3d map >> plot_"//adjustl(trim(pname)) )
    call system("echo set size square >> plot_"//adjustl(trim(pname)) )
    call system("echo set xrange ["//trim(adjustl(trim(minx)))//":"//trim(adjustl(trim(maxx)))//"]>> plot_"//adjustl(trim(pname)) )
    call system("echo set yrange ["//trim(adjustl(trim(miny)))//":"//trim(adjustl(trim(maxy)))//"]>> plot_"//adjustl(trim(pname)) )
    call system("echo splot \'"//trim(pname)//"\'  >> plot_"//adjustl(trim(pname)) )
    call system("echo '#'set term png size 1920,1280 >> plot_"//adjustl(trim(pname)) )
    call system("echo '#'set out\'"//adjustl(trim(pname))//".png\' >> plot_"//adjustl(trim(pname)) )
    call system("echo '#'rep >> plot_"//adjustl(trim(pname)) )
    !
    call system("echo '#' >> plot_"//adjustl(trim(pname)) )
    call system("echo reset >> plot_"//adjustl(trim(pname)) )
    !
    call system("echo set term wxt 1 persist >> plot_"//adjustl(trim(pname)) )
    call system("echo set title \'"//trim(pname)//"\'  >> plot_"//adjustl(trim(pname)) )
    call system("echo set nokey >> plot_"//adjustl(trim(pname)) )
    call system("echo set grid >> plot_"//adjustl(trim(pname)) )
    call system("echo set view 50,10,1,1 >> plot_"//adjustl(trim(pname)) )
    call system("echo splot \'"//trim(pname)//"\' with pm3d >> plot_"//adjustl(trim(pname)) )
    if(present(wlines))call system("echo rep \'"//trim(pname)//"_withlines\' with lines >> plot_"//adjustl(trim(pname)) )
    call system("echo '#'set term png size 1920,1280 >> plot_"//adjustl(trim(pname)) )
    call system("echo '#'set out\'"//adjustl(trim(pname))//".png\' >> plot_"//adjustl(trim(pname)) )
    call system("echo '#'rep >> plot_"//adjustl(trim(pname)) )
  end subroutine splot3D__
  !********************************************************************
  !********************************************************************
  !********************************************************************

end module SLPLOT
