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
 
  real,save                                :: time,old_time,dtime,elapsed_time,eta_time

  integer,parameter :: secs_in_one_day=86400
  integer,parameter :: secs_in_one_hour=3600
  integer,parameter :: secs_in_one_min=60


  public :: print_bar
  public :: start_timer,stop_timer
  public :: eta
  public :: timer

contains

  include "crono.f90"

  include "bar.f90"

end module CHRONOBAR

