!########################################################################
!PROGRAM  : CHRONOBAR
!TYPE     : module
!PURPOSE  : set of routines for timing procedures
!AUTHORS  : Adriano Amaricci
!########################################################################
module TIMER
  USE COMMON_VARS
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


  public :: print_bar
  public :: start_timer,stop_timer
  public :: eta

contains

  include "chrono.f90"

  include "bar.f90"

end module TIMER

