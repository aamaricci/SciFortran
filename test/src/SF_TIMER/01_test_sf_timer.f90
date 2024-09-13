program test_SF_TIMER
  USE SF_CONSTANTS
  USE SF_TIMER
  USE ASSERTING
  implicit none

  integer :: i


  call start_timer("wait 100x30ms")
  do i=1,100
     call wait(30)
     call eta(i,100)
  enddo
  call stop_timer("timer example")
  write(*,*) ""


  call start_progress("progress ETA 100x40ms")
  do i=1,100
     call wait(40)
     call progress_bar_eta(i,100)
  enddo
  call stop_progress()
  write(*,*) ""

end program test_SF_TIMER
