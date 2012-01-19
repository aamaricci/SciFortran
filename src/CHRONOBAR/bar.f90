!+-------------------------------------------------------------------+
!PROGRAM  : 
!TYPE     : Subroutine
!PURPOSE  : 
!+-------------------------------------------------------------------+
subroutine print_bar(i,imax,leta)
  integer :: i,imax
  logical,optional :: leta
  if(mpiID==0)then
     if(present(leta))then
        call delete_bar_eta(i,imax)
        old_time=time
        time=timer()
        dtime        = time-old_time     
        elapsed_time = elapsed_time + dtime
        dtime        = elapsed_time/real(i,4)
        eta_time     = dtime*imax - elapsed_time
        call plot_bar_eta(i,imax,eta_time)
     else
        call delete_bar(i,imax)
        call plot_bar(i,imax)
     endif
     if(i==imax)call close_bar
  endif
end subroutine print_bar
!*********************************************************************
!*********************************************************************
!*********************************************************************


!+-------------------------------------------------------------------+
!PROGRAM  : 
!TYPE     : Subroutine
!PURPOSE  : 
!+-------------------------------------------------------------------+
subroutine plot_bar_eta(i,imax,time)
  integer :: i,imax,k
  real    :: time
  character(len=1) :: bar
  bar = "="
  ms=fraction(time)*1000
  h=time/secs_in_one_hour
  m=(time - h*secs_in_one_hour)/secs_in_one_min
  s=time - h*secs_in_one_hour - m*secs_in_one_min
  write(*,'(2x,1i9,2x,1i3,1a7,i2,a1,i2.2,a1,i2.2,a1,i3.3,1a1,2x,1a2)', advance='no') &
       i,100*i/imax,"% |ETA:",h,":",m,":",s,".",ms,' |'
end subroutine plot_bar_eta
!*********************************************************************
!*********************************************************************
!*********************************************************************



!+-------------------------------------------------------------------+
!PROGRAM  : 
!TYPE     : Subroutine
!PURPOSE  : 
!+-------------------------------------------------------------------+
subroutine delete_bar_eta(i,imax)
  integer :: i,imax,k
  character(len=1) :: back
  back = char(8)
  write(*,'(40a1)', advance='no') (back, k =1,2+9+2+3+7+2+1+2+1+2+1+3+1+2+2)
end subroutine delete_bar_eta
!*********************************************************************
!*********************************************************************
!*********************************************************************




!+-------------------------------------------------------------------+
!PROGRAM  : 
!TYPE     : Subroutine
!PURPOSE  : 
!+-------------------------------------------------------------------+
subroutine plot_bar(i,imax)
  integer :: i,imax,k
  character(len=1) :: bar
  bar = "="
  write(*,'(2x,1i9,2x,1i3,1a1,2x,1a1,256a1)', advance='no') &
       i,100*i/imax,'%','|', (bar, k =1,50*i/imax)
end subroutine plot_bar
!*********************************************************************
!*********************************************************************
!*********************************************************************



!+-------------------------------------------------------------------+
!PROGRAM  : 
!TYPE     : Subroutine
!PURPOSE  : 
!+-------------------------------------------------------------------+
subroutine delete_bar(i,imax)
  integer :: i,imax,k
  character(len=1) :: back
  back = char(8)
  write(*,'(256a1)', advance='no') (back, k =1,(50*i/imax)+9+11)
end subroutine delete_bar
!*********************************************************************
!*********************************************************************
!*********************************************************************



!+-------------------------------------------------------------------+
!PROGRAM  : 
!TYPE     : Subroutine
!PURPOSE  : 
!+-------------------------------------------------------------------+
subroutine close_bar
  write(*,'(a)') '| done.'
  write(*,*)""
end subroutine close_bar
!*********************************************************************
!*********************************************************************
!*********************************************************************
