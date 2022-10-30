program testMINPACK
  USE SCIFOR
  USE ASSERTING
  implicit none
  integer, parameter :: n = 2
  integer            :: iflag
  integer            :: info
  real(8)            :: tol
  real(8)            :: x(n),x0(n),y(n),a(3)
  real(8),dimension(50) :: xdata,ydata,ynoise,xwork,ywork,fvec

  xdata = linspace(0d0,3d0,size(xdata))
  ynoise = 0.1d0*ran_normal(size(xdata))
  ydata = exp_func(xdata,[2d0,1d0,0.2d0])

  call splot("curve.dat",xdata,ydata)
  call splot("curve_noise.dat",xdata,ydata+ynoise)


  xwork = xdata
  ywork = ydata+ynoise
  fvec  = ywork
  a = [2d0,1d0,0.2d0]
  call curvefit(exp_func,a,xwork,ywork)
  ywork = exp_func(xwork,a)
  call assert(sum(ywork-fvec),0d0,"CurveFit LMDIF with Function interface:",tol=1d-6)


  xwork = xdata
  ywork = ydata+ynoise
  fvec  = ywork
  a = [2d0,1d0,0.2d0]
  call curvefit(exp_sub,a,xwork,ywork)
  call exp_sub(xwork,a,ywork)
  call assert(sum(ywork-fvec),0d0,"CurveFit LMDIF with Subroutine interface:",tol=1d-6)

  xwork = xdata
  ywork = ydata+ynoise
  fvec  = ywork
  a = [2d0,1d0,0.2d0]
  call curvefit(exp_func,exp_dfunc,a,xwork,ywork)
  ywork = exp_func(xwork,a)
  call assert(sum(ywork-fvec),0d0,"CurveFit + Grad LMDER with Function interface:",tol=1d-6)


  xwork = xdata
  ywork = ydata+ynoise
  fvec  = ywork
  a = [2d0,1d0,0.2d0]
  call curvefit(exp_sub,exp_dsub,a,xwork,ywork)
  call exp_sub(xwork,a,ywork)
  call assert(sum(ywork-fvec),0d0,"CurveFit + Grad LMDER with Subroutine interface:",tol=1d-6)



contains




  function exp_func(x,a)
    real(8),dimension(:)       :: x
    real(8)                    :: a(:)
    real(8),dimension(size(x)) :: exp_func
    exp_func = a(1)*exp(-a(2)*x) + a(3)
  end function exp_func

  subroutine exp_sub(x,a,func)
    real(8),dimension(:)       :: x
    real(8)                    :: a(:)
    real(8),dimension(size(x)) :: func
    func = a(1)*exp(-a(2)*x) + a(3)
  end subroutine exp_sub

  function exp_dfunc(x,a)
    real(8),dimension(:)               :: x
    real(8),dimension(:)               :: a
    real(8),dimension(size(x),size(a)) :: exp_dfunc
    exp_dfunc(1:,1) = exp(-a(2)*x(:))
    exp_dfunc(1:,2) = -a(1)*x(:)*exp(-a(2)*x(:))
    exp_dfunc(1:,3) = 1d0
  end function exp_dfunc

  subroutine exp_dsub(x,a,exp_dfunc)
    real(8),dimension(:)               :: x
    real(8),dimension(:)               :: a
    real(8),dimension(size(x),size(a)) :: exp_dfunc
    exp_dfunc(1:,1) = exp(-a(2)*x(:))
    exp_dfunc(1:,2) = -a(1)*x(:)*exp(-a(2)*x(:))
    exp_dfunc(1:,3) = 1d0
  end subroutine exp_dsub


  function ran_normal(n)
    integer              :: n
    real(8),dimension(n) :: ran_normal
    integer              :: i
    do i=1,n
       ran_normal(i) = random_normal()
    end do
  end function ran_normal




end program testMINPACK
