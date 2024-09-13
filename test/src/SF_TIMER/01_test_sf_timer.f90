program test_SF_TIMER
  USE SF_CONSTANTS
  USE SF_TIMER
  USE ASSERTING
  implicit none

  integer :: i
  integer(8) :: rate,t0,t1

  call start_timer()

  call system_clock(count_rate=rate);print*,rate
  call system_clock(t0)
  call start_timer("Loop 100x30ms")
  do i=1,100
     call wait(1)
     call eta(i,100)
  enddo
  call stop_timer("timer + ETA example")
  print*,""
  call system_clock(t1)
  print*,real(t1-t0)/rate

  call start_progress("Loop 100x25ms")
  do i=1,100
     call wait(25)
     call progress(i,100)
  enddo
  call stop_progress("progress + ETA example")
  print*,""

  call stop_timer("total time CFR w/ time")
  
end program test_SF_TIMER
