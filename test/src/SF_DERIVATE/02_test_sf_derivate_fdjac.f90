program testFDJAC
  USE SCIFOR
  USE ASSERTING
  implicit none
  integer,parameter            :: Ndim=2,Mdim=3
  real(8),dimension(Ndim)      :: xc,fc
  real(8),dimension(Ndim,Ndim) :: fjac,exact_fjac
  real(8),dimension(Ndim)      :: fgrad,exact_fgrad
  real(8),dimension(Mdim,Ndim) :: gjac,exact_gjac
  integer                      :: i,j

  xc=[1.d0,1.d0]
  fc=fcn(xc)
  !
  exact_fjac(1,1)=1.d0
  exact_fjac(1,2)=1.d0
  exact_fjac(2,1)=0.d0
  exact_fjac(2,2)=-1.d0

  exact_fgrad(1)=2.d0
  exact_fgrad(2)=3.d0

  exact_gjac(1,1)=1.d0
  exact_gjac(1,2)=1.d0
  exact_gjac(2,1)=0.d0
  exact_gjac(2,2)=-1.d0
  exact_gjac(3,1)=2.d0
  exact_gjac(3,2)=3.d0


  !#################################
  ! TEST JACOBIAN:
  !#################################
  call djacobian(fcn,xc,fjac)
  call assert(fjac,exact_fjac,"df(x_0) NxN subroutine w/ function interface",tol=1d-8)

  call djacobian(funcn,xc,fjac)
  call assert(fjac,exact_fjac,"df(x_0) NxN subroutine w/ subroutine interface",tol=1d-8)

  fjac = f_djacobian(fcn,xc)
  call assert(fjac,exact_fjac,"df(x_0) NxN function w/ function interface",tol=1d-8)

  fjac = f_djacobian(funcn,xc)
  call assert(fjac,exact_fjac,"df(x_0) NxN function w/ subroutine interface",tol=1d-8)




  !#################################
  ! TEST JACOBIAN:
  !#################################
  call djacobian(fcm,xc,Mdim,gjac)
  call assert(gjac,exact_gjac,"dg(x_0) NxN subroutine w/ function interface",tol=1d-8)

  call djacobian(funcm,xc,Mdim,gjac)
  call assert(gjac,exact_gjac,"dg(x_0) NxN subroutine w/ subroutine interface",tol=1d-8)

  gjac = f_djacobian(fcm,xc,mdim)
  call assert(gjac,exact_gjac,"dg(x_0) NxN function w/ function interface",tol=1d-8)

  gjac = f_djacobian(funcm,xc,mdim)
  call assert(gjac,exact_gjac,"dg(x_0) NxN function w/ subroutine interface",tol=1d-8)




  !#################################
  ! TEST GRADIENT:
  !#################################
  call dgradient(fc1,xc,fgrad)
  call assert(fgrad,exact_fgrad,"\gradf(x_0) subroutine w/ function interface",tol=1d-8)

  call dgradient(func1,xc,fgrad)
  call assert(fgrad,exact_fgrad,"\gradf(x_0) subroutine w/ subroutine interface",tol=1d-8)

  fgrad = f_dgradient(fc1,xc)
  call assert(fgrad,exact_fgrad,"\gradf(x_0) function w/ function interface",tol=1d-8)

  fgrad = f_dgradient(func1,xc)
  call assert(fgrad,exact_fgrad,"\gradf(x_0) function w/ subroutine interface",tol=1d-8)



contains



  function fcn(x) 
    real(8),dimension(:),intent(in) :: x
    real(8),dimension(size(x))      :: fcn
    fcn(1) = x(1)*x(2) - 2.d0
    fcn(2) = x(1) - x(1)*x(2) + 1.d0
  end function fcn

  subroutine funcn(x,fc) 
    real(8),dimension(:),intent(in) :: x
    real(8),dimension(size(x))      :: fc
    fc = fcn(x)
  end subroutine funcn


  function fcm(x,m)
    real(8),dimension(:),intent(in) :: x
    integer                         :: m
    real(8),dimension(m)            :: fcm
    fcm(1) = x(1)*x(2) - 2.d0
    fcm(2) = x(1) - x(1)*x(2) + 1.d0
    fcm(3) = x(1)**2 +3.d0*x(2)
  end function fcm

  subroutine funcm(x,m,fc)
    real(8),dimension(:),intent(in) :: x
    integer                         :: m
    real(8),dimension(m)            :: fc
    fc = fcm(x,m)
  end subroutine funcm

  function fc1(x) 
    real(8),dimension(:),intent(in) :: x
    real(8)                         :: fc1
    fc1 = x(1)**2 + 3.d0*x(2)
  end function fc1

  subroutine func1(x,fc) 
    real(8),dimension(:),intent(in) :: x
    real(8)                         :: fc
    fc = fc1(x)
  end subroutine func1


end program testFDJAC
