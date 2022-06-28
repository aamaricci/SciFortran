program testSF_GAUSS
  USE SCIFOR
  USE ASSERTING
  implicit none

  integer             :: L,Lf,N
  real(8)             :: epsabs,epsrel
  integer             :: ier,i,j,dim
  integer             :: neval,key,inf,integr
  real(8)             :: result,abserr,true,resabs,resasc,tol
  real(8),allocatable :: fx(:),x(:),a(:),b(:),y(:),fxy(:,:)

  key=10
  tol=1d-6
  N=1000
  epsrel=tol

  !1D:
  dim = 1
  allocate(a(dim),b(dim))
  a = 0d0
  b = pi
  true=0.06278740d0
  !
  neval=0
  call gauss_quad(f_cos100sin,a(1),b(1),result,method=key,tol=epsrel,ierr=ier,err=abserr)
  call assert(result,true,"GAUSS QUADRATURE cos(100sin(x)) Dim=1 key="//str(key),tol=1d-6)  
  deallocate(a,b)


  !2D:
  dim = 2
  allocate(a(dim),b(dim))
  a = [0.3827812d0,22.91281d0]
  b = [1.02928d0,-11.928d0]
  true=test_2d_integral_xy(a,b)
  neval=0
  call gauss_quad(test_2d_func_xy,a,b,result,&
       method=key,tol=epsrel,ierr=ier,err=abserr)
  call assert(result,true,"GAUSS QUADRATURE x*y Dim=2 key="//str(key),tol=1d-6)  
  !

  a = [0d0,0d0] ![xl,yl]
  b = [pi,2*pi]   ![xu,yu]
  true=test_2d_integral_sx_cy(a,b)
  neval=0
  call gauss_quad(test_2d_func_sx_cy,a,b,result,&
       method=key,tol=epsrel,ierr=ier,err=abserr)
  call assert(result,true,"GAUSS QUADRATURE sin(x)+cos(y) Dim=2 key="//str(key),tol=1d-6)  
  !
  deallocate(a,b)


  !3D:
  dim = 3
  allocate(a(dim),b(dim))
  a = [2d0,1d0,0d0]
  b = [3d0,2d0,1d0]
  true=test_3d_integral_xyz()
  neval=0
  call gauss_quad(test_3d_func_xyz,a,b,result,&
       method=key,tol=epsrel,ierr=ier,err=abserr)
  call assert(result,true,"GAUSS QUADRATURE x*y*z Dim=3 key="//str(key),tol=1d-6)  
  deallocate(a,b)


  !4D:
  dim = 4
  allocate(a(dim),b(dim))
  a = 0d0
  b = 1d0
  true=test_nd_integral()
  neval=0
  call gauss_quad(test_nd_func,a,b,result,&
       method=key,tol=epsrel,ierr=ier,err=abserr)
  call assert(result,true,"GAUSS QUADRATURE 1  Dim=4 key="//str(key),tol=1d-6)  
  deallocate(a,b)



contains

  subroutine write_out(dim,key)
    integer :: dim,key
    real(8) :: pow,error
    pow = 1d0/dim
    error = abs(true - result)
    write(*,*)'Numerical integral                   =', result
    write(*,*)'Error                                =', error
    write(*,*)'Estimated error                      =', abserr
    write(*,*)'Number of eval                       =', neval,"==",int(neval**pow),"^",dim
    if(error<epsrel*1d-2)then
       write ( *, '(a,i8,a)' ) 'adaptive automatic integrator using a Gauss-Kronrod rules. METHOD=',key," PASSED"
    else
       write ( *, '(a,i8,a)' ) 'adaptive automatic integrator using a Gauss-Kronrod rules. METHOD=',key," FAILED"
    endif
    write(*,*)""
  end subroutine write_out


  !! F02 is the integrand function COS(100*SIN(X)).
  function f_cos100sin ( x ) result(f02)
    real(8) ::  x
    real(8) ::  f02
    f02 = cos ( 100d0 * sin ( x ) )
    neval = neval + 1
  end function f_cos100sin


  function test_2d_func_sx_cy(x) result(f)
    real(8),dimension(:) :: x
    real(8)                         :: f
    f  = sin(x(1)) + cos(x(2))
    neval = neval + 1
  end function test_2d_func_sx_cy

  function test_2d_func_xy(x) result(f)
    real(8),dimension(:) :: x
    real(8)                         :: f
    f = x(1) * x(2)
    neval = neval + 1
  end function test_2d_func_xy

  function test_3d_func_xyz(x) result(f)
    real(8),dimension(:) :: x
    real(8)                         :: f
    f = 8d0 * x(1) * x(2) * x(3)
    neval = neval + 1
  end function test_3d_func_xyz


  function test_nd_func(x) result(f)
    real(8),dimension(:) :: x
    real(8)              :: f
    f = 1d0
    neval = neval + 1
  end function test_nd_func




  function test_2d_integral_sx_cy(xlo,xup) result(f)
    real(8),dimension(2) :: xlo,xup
    real(8)                         :: xu,xl,yu,yl
    real(8)                         :: f
    xl = xlo(1);yl=xlo(2)
    xu = xup(1);yu=xup(2)
    f  = (-cos(xu)+cos(xl))*(yu-yl) + (xu-xl)*(sin(yu)-sin(yl))
  end function test_2d_integral_sx_cy

  function test_2d_integral_xy(xlo,xup) result(f)
    real(8),dimension(2) :: xlo,xup
    real(8)                         :: xu,xl,yu,yl
    real(8)                         :: f
    xl = xlo(1);yl=xlo(2)
    xu = xup(1);yu=xup(2)
    f = ( xu**2/2d0 - xl**2/2d0 ) * ( yu**2/2d0 - yl**2/2d0 )
  end function test_2d_integral_xy

  function test_3d_integral_xyz() result(f)
    real(8)              :: f
    f = 15d0
  end function test_3d_integral_xyz

  function test_nd_integral() result(f)
    real(8)              :: f
    f = 1d0
  end function test_nd_integral




end program testSF_GAUSS
