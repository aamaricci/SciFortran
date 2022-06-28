program test_SF_CONSTANTS
  USE SF_CONSTANTS
  USE ASSERTING
  implicit none

  real(8) :: my_pi,my_e,my_sqrt6,my_c

  my_pi = acos(-1d0)
  call assert(my_pi,pi,"PI",tol=1d-12)

  my_e = exp(1d0)
  call assert(my_e,euler,"E",tol=1d-12)

  my_sqrt6=sqrt(6d0)
  call assert(my_sqrt6,sqrt6,"E",tol=1d-12)


  my_c = 2.99792458000D+08
  call assert(my_c,speed_of_light_in_vacuum,"SPEED OF LIGHT IN VACUUM",tol=1d-8)
  
  call assert(zero,dcmplx(0d0,0d0),"Zero")
  call assert(one,dcmplx(1d0,0d0),"One")
  call assert(xi,dcmplx(0d0,1d0),"Xi")

  my_pi=0d0
  call assert(isinfty(1d0/my_pi),.true.,"IS INFTY")
  call assert(isnan(my_pi/my_pi),.true.,"IS NaN")

  call timestamp()
  call wait(1.8437d0)
  call timestamp()

end program test_SF_CONSTANTS
