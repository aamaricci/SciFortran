program test_fsolve
  USE SCIFOR
  USE ASSERTING
  implicit none
  integer, parameter :: n = 2
  integer            :: iflag
  integer            :: info
  real(8)            :: tol
  real(8)            :: x(n),x0(n),y(n),a(3)
  real(8),dimension(50) :: xdata,ydata,ynoise,xwork,ywork,fvec

  x = [4d0,0d0]
  y = hybrd_func(x)
  tol = 1.d-15
  call fsolve(hybrd_func,x,tol,info)
  y = hybrd_func(x)
  call assert([x,y],[1d0,1d0,0d0,0d0],"Fsolve no Grad with Function interface:")


  x = [4d0,0d0]
  call hybrd_sub(x,y)
  tol = 1.d-15
  call fsolve(hybrd_sub,x,tol,info)
  call hybrd_sub(x,y)
  call assert([x,y],[1d0,1d0,0d0,0d0],"Fsolve no Grad with Subroutine interface:")


  x = [4d0,0d0]
  y = hybrj_func(x)
  tol = 1.d-15
  call fsolve(hybrj_func,hybrj_dfunc,x,tol,info)
  y = hybrj_func(x)
  call assert([x,y],[1d0,1d0,0d0,0d0],"Fsolve + Grad with function interface:")


  x = [4d0,0d0]
  call hybrj_sub(x,y)
  tol = 1.d-15
  call fsolve(hybrj_sub,hybrj_dsub,x,tol,info)
  call hybrj_sub(x,y)
  call assert([x,y],[1d0,1d0,0d0,0d0],"Fsolve + Grad with Subroutine interface:")




contains




  function hybrd_func(x) result(ffunc)
    implicit none
    real(8),dimension(:),intent(in) :: x
    real(8),dimension(size(x))      :: ffunc
    ffunc(1) = x(1)**2 - 10.d0*x(1) + x(2)**2 + 8.d0
    ffunc(2) = x(1)*x(2)**2 + x(1) - 10.d0*x(2) + 8.d0
    write(*,"(4F16.9)")x(1),x(2),ffunc(1),ffunc(2)
  end function hybrd_func


  subroutine hybrd_sub(x,ff)
    implicit none
    real(8),dimension(:),intent(in) :: x
    real(8),dimension(size(x))      :: ff
    ff(1) = x(1)**2 - 10.d0*x(1) + x(2)**2 + 8.d0
    ff(2) = x(1)*x(2)**2 + x(1) - 10.d0*x(2) + 8.d0
    write(*,"(4F16.9)")x(1),x(2),ff(1),ff(2)
  end subroutine hybrd_sub



  function hybrj_func(x) result(ffunc)
    real(8),dimension(:),intent(in) :: x
    real(8),dimension(size(x))      :: ffunc
    ffunc(1) = x(1)**2 - 10.d0*x(1) + x(2)**2 + 8.d0
    ffunc(2) = x(1)*x(2)**2 + x(1) - 10.d0*x(2) + 8.d0
    write(*,"(4F16.9)")x(1),x(2),ffunc(1),ffunc(2)
  end function hybrj_func
  function hybrj_dfunc(x) result(dfunc)
    real(8),dimension(:),intent(in)    :: x
    real(8),dimension(size(x),size(x)) :: dfunc
    dfunc(1,1) = 2d0*x(1) - 10d0
    dfunc(1,2) = 2d0*x(2)
    dfunc(2,1) = x(2)*x(2) + 1d0
    dfunc(2,2) = 2d0*x(1)*x(2) - 10d0
  end function hybrj_dfunc


  subroutine hybrj_sub(x,ffunc)
    real(8),dimension(:),intent(in) :: x
    real(8),dimension(size(x))      :: ffunc
    ffunc(1) = x(1)**2 - 10.d0*x(1) + x(2)**2 + 8.d0
    ffunc(2) = x(1)*x(2)**2 + x(1) - 10.d0*x(2) + 8.d0
    write(*,"(4F16.9)")x(1),x(2),ffunc(1),ffunc(2)
  end subroutine hybrj_sub
  subroutine hybrj_dsub(x,dfunc)
    real(8),dimension(:),intent(in)    :: x
    real(8),dimension(size(x),size(x)) :: dfunc
    dfunc(1,1) = 2d0*x(1) - 10d0
    dfunc(1,2) = 2d0*x(2)
    dfunc(2,1) = x(2)*x(2) + 1d0
    dfunc(2,2) = 2d0*x(1)*x(2) - 10d0
  end subroutine hybrj_dsub



end program test_fsolve
