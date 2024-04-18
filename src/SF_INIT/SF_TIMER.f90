module SF_TIMER
  implicit none
  private

  integer(4),dimension(8),save             :: data
  integer(4),save                          :: year
  integer(4),save                          :: mese
  integer(4),save                          :: day
  integer(4),save                          :: h
  integer(4),save                          :: m
  integer(4),save                          :: s
  integer(4),save                          :: ms
  character(len=9),parameter,dimension(12) :: month = (/ &
       'January  ', 'February ', 'March    ', 'April    ', &
       'May      ', 'June     ', 'July     ', 'August   ', &
       'September', 'October  ', 'November ', 'December ' /)

  integer,save                             :: timer_index=0
  integer(4),dimension(100,8),save         :: timer_start,timer_stop,timer0,timer1

  real,save                                :: time,old_time,dtime,elapsed_time,eta_time

  integer,parameter :: secs_in_one_day=86400
  integer,parameter :: secs_in_one_hour=3600
  integer,parameter :: secs_in_one_min=60

  integer  :: Funit

  public :: start_timer,stop_timer
  public :: start_progress,stop_progress
  public :: eta
  public :: progress
  public :: progress_bar
  public :: progress_bar_eta

contains

  !+-------------------------------------------------------------------+
  !PURPOSE  : start a timer to measure elapsed time between two call
  !+-------------------------------------------------------------------+
  subroutine start_timer(title)
    character(len=*),optional :: title
    if(present(title))write(*,"(A)")trim(title)//":"
    timer_index=timer_index+1
    if(timer_index>size(timer_start,1))then
       stop "Error in cronograph: too many timers started"
    endif
    call date_and_time(values=timer_start(timer_index,:))
    timer0(timer_index,:)=timer_start(timer_index,:)
    !init variables for ETA:
    elapsed_time =0.0
    old_time     =0.0
    time         =0.0
    !
  end subroutine start_timer
  !
  subroutine start_progress(unit)
    integer,optional :: unit
    funit=6;if(present(unit))funit=unit
    !open(funit,carriagecontrol='fortran')
    open(funit)
    call start_timer
  end subroutine start_progress

  !+-------------------------------------------------------------------+
  !PURPOSE  : stop the timer and get the partial time
  !+-------------------------------------------------------------------+
  subroutine stop_timer(title,unit)
    character(len=*),optional :: title
    integer,optional          :: unit
    integer                   :: unit_
    integer,dimension(8)      :: itimer
    unit_=6;if(present(unit))unit_=unit
    call date_and_time(values=timer_stop(timer_index,:))
    itimer=time_difference(timer_stop(timer_index,:),timer_start(timer_index,:))
    if(present(title))then
       call print_total_time(itimer,unit,title)
    else
       call print_total_time(itimer,unit)
    endif
    timer_start(timer_index,:)=0
    timer_stop(timer_index,:)=0
    if(timer_index>1)then
       timer_index=timer_index-1
    else
       timer_index=0
    endif
  end subroutine stop_timer
  !
  subroutine stop_progress
    close(funit)
    call stop_timer
  end subroutine stop_progress


  subroutine print_total_time(dummy,unit,title)
    integer,optional          :: unit
    character(len=*),optional :: title
    integer                   :: unit_
    integer(4),dimension(8)   :: dummy
    unit_=6;if(present(unit))unit_=unit
    year = dummy(1)
    mese = dummy(2)
    day  = dummy(3)
    h    = dummy(5)
    m    = dummy(6)
    s    = dummy(7)
    ms   = dummy(8)
    if(present(title))then
       write(unit_,"(a,i3,a1,i2.2,a1,i2.2,a1,i3.3,A2,A)")"Total time [h:m:s.ms]: ",h,":",m,":",s,".",ms," :",trim(title)
    else
       write(unit_,"(a,i3,a1,i2.2,a1,i2.2,a1,i3.3)")"Total time [h:m:s.ms]: ",h,":",m,":",s,".",ms
    endif
    !write(unit_,*)""
  end subroutine print_total_time





  !+-------------------------------------------------------------------+
  !PURPOSE  : get time difference between two events
  !+-------------------------------------------------------------------+
  function time_difference(data1,data0)
    integer(4),dimension(8) :: data1,data0,dummy,time_difference
    dummy =data1-data0
    year = dummy(1)
    mese = dummy(2)
    day  = dummy(3)
    h    = dummy(5)
    m    = dummy(6)
    s    = dummy(7)
    ms   = dummy(8)
    if(h<0)then
       day=day-1
       h=24+h
    endif
    if(m<0)then
       h=h-1
       m=60+m
    endif
    if(s<0)then
       m=m-1
       s=60+s
    endif
    if(ms<0)then
       s=s-1
       ms=1000+ms
    endif
    time_difference(1)=year
    time_difference(2)=mese
    time_difference(3)=day
    time_difference(5)=h
    time_difference(6)=m
    time_difference(7)=s
    time_difference(8)=ms
  end function time_difference



  !+-------------------------------------------------------------------+
  !PURPOSE  : get Expected Time of Arrival
  !+-------------------------------------------------------------------+
  subroutine eta(i,L,unit,file,step)
    integer                 :: i,L!,k
    integer,optional        :: step
    integer,save            :: mod_print
    integer                 :: percent,iprint
    integer,save            :: older=0,oldiprint=0
    logical                 :: esc,fullprint
    integer(4),dimension(8) :: dummy
    integer,optional        :: unit
    integer,save            :: unit_
    character(len=*),optional :: file
    character(len=16)      :: string
    character(len=80)      :: message
    logical,save            :: lentry=.true.
    !Preambolo:
    if(lentry)then
       unit_=6 ; if(present(unit))unit_=unit
       mod_print=10;if(present(step))mod_print=step
       if(unit_==5)unit_=6     !do not write to stdin
       if(present(file))then
          unit_=719
          open(unit_,file=trim(adjustl(trim(file))))
          write(unit_,*)
          write(*,"(2x,A)")"+ETA --> "//trim(adjustl(trim(file)))
       else
          write(string,"(I4)")unit_
          write(unit_,*)""
          write(*,"(2x,A,I3)")"+ETA --> fort."//trim(adjustl(trim(string)))
       endif
       lentry=.false.
    endif
    !
    if(i==L)lentry=.true.
    !
    !avoid repetition of percentage (within the error)
    percent=100*i/L
    if(percent==0)return
    if(percent==older)return
    if(percent<mod_print)return
    older=percent
    !
    !set step for printing:
    esc=.true.
    iprint=percent/mod_print
    !if(percent<=mod_print .OR. iprint/=oldiprint)esc=.false.
    if(iprint/=oldiprint)esc=.false.
    if(esc)return
    oldiprint=iprint
    !
    !check if fullprint (date) is needed
    fullprint=.false.;if(percent<=1 .OR. percent==10 .OR. percent==50)fullprint=.true.
    !
    old_time=time
    time=total_time()
    dtime        = time-old_time     
    elapsed_time = elapsed_time + dtime
    dtime        = elapsed_time/real(i,4)
    eta_time     = dtime*L - elapsed_time
    ms=int(fraction(eta_time)*1000.0)
    h =int(eta_time/secs_in_one_hour)
    m =int((eta_time - h*secs_in_one_hour)/secs_in_one_min)
    s =int(eta_time - h*secs_in_one_hour - m*secs_in_one_min)
    !
    call date_and_time (values=dummy)
    if(fullprint)then
       write(message,"(1i3,1a7,i2,a1,i2.2,a1,i2.2,a1,i3.3,a2,i2,1x,a,1x,i4,2x,i2,a1,i2.2,a1,i2.2)")&
            percent,"% |ETA: ",h,":",m,":",s,".",ms," @",dummy(3),trim(month(dummy(2))),dummy(1), &
            dummy(5),':',dummy(6),':',dummy(7)
    else
       write(message,"(1i3,1a7,i2,a1,i2.2,a1,i2.2,a1,i3.3)")percent,"% |ETA: ",h,":",m,":",s,".",ms
    endif
    write(unit_,*)trim(message)
    !if(lentry.AND.present(file))close(unit_)
  end subroutine eta



  !+-------------------------------------------------------------------+
  !PURPOSE  : a total_time
  !+-------------------------------------------------------------------+
  function total_time()
    real(4)                 :: total_time
    integer(4),dimension(8) :: dummy
    call date_and_time(values=timer1(timer_index,:))
    dummy= time_difference(timer1(timer_index,:),timer0(timer_index,:))
    year = dummy(1)
    mese = dummy(2)
    day  = dummy(3)
    h    = dummy(5)
    m    = dummy(6)
    s    = dummy(7)
    ms   = dummy(8)
    total_time= real(dble(ms)/1000.d0 &
         + dble(s)               &
         + dble(m)*60.d0         &
         + dble(h)*60.d0**2      &
         + dble(day)*24.d0*60.d0**2,4)
  end function total_time



  !+-------------------------------------------------------------------+
  !TYPE     : Subroutine
  !+-------------------------------------------------------------------+
  subroutine progress(i,imax)
    integer            :: i,imax,k,jmax
    character(len=7)  :: bar="> ???% "
    write(unit=bar(3:5),fmt="(I3,$)")100*i/imax
    write(unit=funit,fmt="(A1,A1,A7,$)")'+',char(13), bar
    if(i==imax)write(funit,*)
  end subroutine progress



  !+-------------------------------------------------------------------+
  !TYPE     : Subroutine
  !+-------------------------------------------------------------------+
  subroutine progress_bar(i,imax)
    integer            :: i,imax,k,jmax
    character(len=57)  :: bar="???% |                                                  |"
    write(unit=bar(1:3),fmt="(I3,$)")100*i/imax
    jmax=50*i/imax
    do k=1,jmax
       bar(6+k:6+k)="*"
    enddo
    write(unit=funit,fmt="(A1,A1,A57,$)")'+',char(13), bar
    if(i==imax)write(funit,*)
  end subroutine progress_bar



  !+-------------------------------------------------------------------+
  !TYPE     : Subroutine
  !+-------------------------------------------------------------------+
  subroutine progress_bar_eta(i,imax)
    integer           :: i,imax
    integer           :: k,jmax
    !                   !bar="100% |-------------50char-------------------------------| ETA "
    character(len=62) :: bar="???% |                                                  | ETA "
    old_time=time
    time=total_time()
    dtime        = time-old_time     
    elapsed_time = elapsed_time + dtime
    dtime        = elapsed_time/real(i,4)
    eta_time     = dtime*imax - elapsed_time
    ms=int(fraction(eta_time)*1000.0)
    h =int(eta_time/secs_in_one_hour)
    m =int((eta_time - h*secs_in_one_hour)/secs_in_one_min)
    s =int(eta_time - h*secs_in_one_hour - m*secs_in_one_min)
    write(unit=bar(1:3),fmt="(I3,$)")100*i/imax
    jmax=50*i/imax
    do k=1,jmax
       bar(6+k:6+k)="*"
    enddo
    write(unit=funit,fmt="(A1,A1,A62,I2,A1,I2.2,A1,I2.2,A1,I3.3,$)")'+',char(13), bar,h,":",m,":",s,".",ms
    if(i==imax)write(funit,*)
  end subroutine progress_bar_eta


end module SF_TIMER

