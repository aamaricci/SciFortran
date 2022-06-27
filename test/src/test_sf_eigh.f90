program test_eigh
  use scifor
  use asserting
  implicit none

  real(8) :: d(2),u(1)
  real(8), dimension(2,2) :: Ev

  real(8) :: e(2)
  logical :: info
  
  e(1) = -1.d0  !+4.d0 !Induced Error
  e(2) = 3.d0
  d=1.d0
  u=2.d0

  call eigh(d,u,Ev)
  call assert_arr_d(d,e,info)
  
  if(info) then  
     write(*,"(A)") "test_eigh exit status is "//achar(27)//"[32m POSITIVE "//achar(27)//"[0m."
  else
     write(*,"(A)") "test_eigh exit status is "//achar(27)//"[31m NEGATIVE "//achar(27)//"[0m."
     error stop 2
  endif  
  
  
end program test_eigh
