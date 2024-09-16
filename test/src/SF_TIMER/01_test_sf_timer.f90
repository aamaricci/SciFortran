program test_SF_TIMER
  USE SF_CONSTANTS
  USE SF_TIMER
  USE ASSERTING
  USE SF_MPI
  implicit none

  integer :: i
  integer :: rank,size
  
  call init_MPI()
  rank = get_rank_MPI()
  size = get_size_MPI()
  
  call start_timer()

  call start_timer("Loop 100x30ms = 3s",unit=6+2*rank)
  do i=1+rank,100,size
     call wait(30)
     call eta(i,100)
  enddo
  call stop_timer("timer + ETA example")

  
  call start_progress("Loop 100x25ms")
  do i=1,100
     call wait(25)
     call progress(i,100)
  enddo
  call stop_progress("progress + ETA example")


  call stop_timer("Total time")
  
  call finalize_MPI()
end program test_SF_TIMER
