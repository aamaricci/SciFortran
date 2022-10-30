program testDERIVE
  USE SCIFOR
  USE ASSERTING
  implicit none

  integer,parameter    :: L=20
  real(8),dimension(L) :: x,f,df,df2,df3,df4,df5,df6
  real(8),dimension(L) :: test
  real(8)              :: dx
  integer :: i

  x=linspace(0d0,pi,L,mesh=dx)
  f   =  cos(x)
  df  = -sin(x)
  df2 = -cos(x)
  df3 =  sin(x)
  df4 =  cos(x)                  != f
  df5 = -sin(x)
  df6 = -cos(x)



  test  = derivative(f,dx,1)
  call assert(test(2:L-1),df(2:L-1),"derivative rule 121 df",tol=1d-2)
  !
  test = derivative(f,dx,2)
  call assert(test,df,"derivative rule 222 df",tol=1d-2)
  !
  test = derivative(f,dx)       !default 4
  call assert(test,df,"derivative rule 444 df",tol=1d-4)
  !
  test = derivative(f,dx,6)
  call assert(test,df,"derivative rule 666 df",tol=1d-5)




  test = derivative2(f,dx,2)
  call assert(test,df2,"derivative rule 222 d^2f",tol=1d-1)
  !
  test = derivative2(f,dx)
  call assert(test,df2,"derivative rule 444 d^2f",tol=1d-2)
  !
  test = derivative2(f,dx,6)
  call assert(test,df2,"derivative rule 666 d^2f",tol=1d-4)




  test = derivative3(f,dx,2)
  call assert(test,df3,"derivative rule 222 d^3f",tol=1d-1)
  !
  test = derivative3(f,dx)
  call assert(test,df3,"derivative rule 444 d^3f",tol=1d-2)
  !
  test = derivative3(f,dx,6)
  call assert(test,df3,"derivative rule 666 d^3f",tol=1d-4)




  test = derivative4(f,dx,2)
  call assert(test,df4,"derivative rule 222 d^4f",tol=1d-1)
  !
  test = derivative4(f,dx)
  call assert(test,df4,"derivative rule 444 d^4f",tol=1d-2)
  !
  test = derivative4(f,dx,6)
  call assert(test,df4,"derivative rule 666 d^4f",tol=1d-3)


  test = derivativeN(f,dx,3)
  call assert(test,df3,"derivative N d^3f",tol=1d-2)

end program testDERIVE
