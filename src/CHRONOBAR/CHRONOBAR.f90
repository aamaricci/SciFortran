!########################################################################
!PROGRAM  : CHRONOBAR
!TYPE     : module
!PURPOSE  : set of routines for timing procedures
!AUTHORS  : Adriano Amaricci
!########################################################################
module CHRONOBAR
  USE COMMON_VARS
  USE MPI
  implicit none
  private

  integer,save                             :: timer_index=0
  integer(4),dimension(100,8),save         :: timer_start,timer_stop,timer0,timer1
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
  real,save                                :: time,old_time,dtime,elapsed_time,eta_time

  integer,parameter :: secs_in_one_day=86400
  integer,parameter :: secs_in_one_hour=3600
  integer,parameter :: secs_in_one_min=60


  public :: print_bar
  public :: timestamp
  public :: start_timer,stop_timer
  public :: eta
  public :: timer

contains
  !+-------------------------------------------------------------------+
  !PURPOSE  : prints the current YMDHMS date as a time stamp.
  ! Example: 31 May 2001   9:45:54.872 AM
  !+-------------------------------------------------------------------+
  subroutine timestamp
    if(mpiID==0)then
       call date_and_time(values=data)
       call print_date(data)
    endif
    return
  end subroutine timestamp



  !*********************************************************************
  !*********************************************************************
  !*********************************************************************


  !+-------------------------------------------------------------------+
  !PURPOSE  : print actual date
  !+-------------------------------------------------------------------+
  subroutine print_date(dummy)
    integer(4),dimension(8) :: dummy
    year = dummy(1)
    mese = dummy(2)
    day  = dummy(3)
    h    = dummy(5)
    m    = dummy(6)
    s    = dummy(7)
    ms   = dummy(8)
    write(*,"(i2,1x,a,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3)")&
         day,trim(month(mese)),year, h,':',m,':',s,'.',ms
    write(*,*)""
  end subroutine print_date

  include "crono.f90"

  include "bar.f90"

end module CHRONOBAR

