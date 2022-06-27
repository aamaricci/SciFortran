program test_eigh
  use scifor
  use asserting
  implicit none

  real(8),dimension(:),allocatable :: f
  real(8) :: aa,bb,int,res=0.d0,dx
  integer :: N,j
  logical :: info
  N=1000000
  aa=0.d0
  bb=pi2
  dx=(bb-aa)/dble(N)
  allocate(f(N))
  do j=1,N
     f(j)=cos(aa+dx*j)
  enddo
  int = trapz(f,aa,bb)
  call assert_d(int,res,info)
  
  if(info) then  
     write(*,"(A)") "test_trapz exit status is "//achar(27)//"[32m POSITIVE "//achar(27)//"[0m."
  else
     write(*,"(A)") "test_trapz exit status is "//achar(27)//"[31m NEGATIVE "//achar(27)//"[0m."
     error stop 2
  endif  
  
  
end program test_eigh

