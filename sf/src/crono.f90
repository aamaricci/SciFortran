
!+-------------------------------------------------------------------+
!PURPOSE  : start a timer to measure elapsed time between two call
!+-------------------------------------------------------------------+
subroutine start_timer
  if(mpiID==0)then
     timer_index=timer_index+1
     if(timer_index>size(timer_start,1))then
        call error("Error in cronograph: too many timer started")
     endif
     call date_and_time(values=timer_start(timer_index,:))
     timer0(timer_index,:)=timer_start(timer_index,:)
     !init variables for ETA:
     elapsed_time =0.0
     old_time     =0.0
     time         =0.0
  endif
end subroutine start_timer



!*********************************************************************
!*********************************************************************
!*********************************************************************



!+-------------------------------------------------------------------+
!PURPOSE  : stop the timer and get the partial time
!+-------------------------------------------------------------------+
subroutine stop_timer
  integer(4),dimension(8) :: itimer
  if(mpiID==0)then
     call date_and_time(values=timer_stop(timer_index,:))
     itimer=time_difference(timer_stop(timer_index,:),timer_start(timer_index,:))
     call print_timer(itimer)
     timer_start(timer_index,:)=0
     timer_stop(timer_index,:)=0
     if(timer_index>1)then
        timer_index=timer_index-1
     else
        timer_index=0
     endif
  endif
end subroutine stop_timer


!*********************************************************************
!*********************************************************************
!*********************************************************************



!+-------------------------------------------------------------------+
!PURPOSE  : a timer
!+-------------------------------------------------------------------+
function timer()
  real(4)                 :: timer
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
  timer= real(dble(ms)/1000.d0 &
       + dble(s)               &
       + dble(m)*60.d0         &
       + dble(h)*60.d0**2      &
       + dble(day)*24.d0*60.d0**2,4)
end function timer
!*********************************************************************
!*********************************************************************
!*********************************************************************




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



!*********************************************************************
!*********************************************************************
!*********************************************************************




subroutine print_timer(dummy)
  integer(4),dimension(8) :: dummy
  year = dummy(1)
  mese = dummy(2)
  day  = dummy(3)
  h    = dummy(5)
  m    = dummy(6)
  s    = dummy(7)
  ms   = dummy(8)
  write(*,"(a,i3,a1,i2.2,a1,i2.2,a1,i3.3)")"Total time [h:m:s.ms]: ",h,":",m,":",s,".",ms
  write(*,*)""
end subroutine print_timer



!*********************************************************************
!*********************************************************************
!*********************************************************************



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
  if(mpiID==0)then
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

     if(i==L)lentry=.true.

     !avoid repetition of percentage (within the error)
     percent=100*i/L ; if(percent==older)return
     older=percent

     !set step for printing:
     esc=.true.
     iprint=percent/mod_print
     if(percent<=mod_print .OR. iprint/=oldiprint)esc=.false.
     if(esc)return
     oldiprint=iprint

     !check if fullprint (date) is needed
     fullprint=.false.;if(percent<=1 .OR. percent==10 .OR. percent==50)fullprint=.true.

     old_time=time
     time=timer()
     dtime        = time-old_time     
     elapsed_time = elapsed_time + dtime
     dtime        = elapsed_time/real(i,4)
     eta_time     = dtime*L - elapsed_time
     ms=int(fraction(eta_time)*1000.0)
     h =int(eta_time/secs_in_one_hour)
     m =int((eta_time - h*secs_in_one_hour)/secs_in_one_min)
     s =int(eta_time - h*secs_in_one_hour - m*secs_in_one_min)

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
  endif
end subroutine eta
!*********************************************************************
!*********************************************************************
!*********************************************************************



