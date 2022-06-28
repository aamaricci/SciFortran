program test_SF_QUAD
  USE SCIFOR
  USE ASSERTING
  implicit none


  integer             :: L,Lf,N
  real(8)             :: dh,a,b,c,omega
  real(8)             :: epsabs,epsrel,gamma,ci,alfa,beta,si
  integer             :: ier,i
  integer             :: neval,key,inf,integr
  real(8)             :: result,abserr,true,resabs,resasc,integral
  real(8),allocatable :: func(:),x(:),singular_points(:)


  L=1000
  allocate(func(L),x(L))
  a=0d0
  b=1d0
  x=linspace(a,b,L,mesh=dh)
  func=sqrt(x)  


  integral = trapz(fsqrt,a,b,L)
  call assert(integral,2d0/3d0,"TRAPZ",tol=1d-2)
  integral = trapz(func,dh)
  call assert(integral,2d0/3d0,"TRAPZ",tol=1d-2)
  integral = trapz(func,a,b)
  call assert(integral,2d0/3d0,"TRAPZ",tol=1d-2)

  integral = simps(fsqrt,a,b,L)
  call assert(integral,2d0/3d0,"SIMPS",tol=1d-4)
  integral = simps(func,dh)
  call assert(integral,2d0/3d0,"SIMPS",tol=1d-4)
  integral = simps(func,a,b)
  call assert(integral,2d0/3d0,"SIMPS",tol=1d-4)



  deallocate(func,x)




  do key=1,6
     a = 0d0
     b = acos(-1d0)
     epsabs=0d0
     epsrel=1d-5
     true=0.06278740D+00
     call quad(f_cos100sin,a,b,key=key,epsabs=epsabs,epsrel=epsrel,verbose=.false.,result=result)
     call assert(result,true,"QAG cos(100sin(x)) with key="//str(key),tol=10*epsrel)
  enddo




  a = 0d0
  epsabs=0d0
  epsrel=1d-4
  true=-acos(-1d0)*log(10d0)/20d0
  inf=1                         !inf=-1 (-infty:a); inf=1 (a:infty);  inf=2 (-infty:infty)
  call quad(f_log_over_px,a,inf=inf,epsabs=epsabs,epsrel=epsrel,verbose=.false.,result=result)
  call assert(result,true,"QAGI LOG(X)/(1+100X^2) infinite interval",tol=epsrel)


  a = 0d0
  b = 1d0
  epsabs=0d0
  epsrel=1d-5
  true=-4d0
  call quad(f_log_over_sqrt,a,b,singular_endpoint=.true.,epsabs=epsabs,epsrel=epsrel,verbose=.false.,result=result)
  call assert(result,true,"QAGS LOG(X)/SQRT(X) endpoint singularities",tol=10*epsrel)





  write ( 100, '(a)' ) 'Test QAGP adaptive integrator that can handle singularities of the integrand at user specified points.'
  a = 0d0
  b = 3d0
  epsabs=0d0
  epsrel=1d-6
  true = 61d0*log(2d0) + 77d0*log(7d0)/4d0 - 27d0
  !
  allocate(singular_points(2))
  singular_points(1) = 1d0
  singular_points(2) = sqrt(2d0)
  !
  call quad(f_x3_log,a,b,singular_points=singular_points,epsabs=epsabs,epsrel=epsrel,verbose=.false.,result=result)
  call assert(result,true,"QAGP x^3.log(abs((x^2-1)*(x^2-2))) specified singularities",tol=10*epsrel)



  a = -1d0
  b =  5d0
  epsabs=0d0
  epsrel=1d-5
  true = log(125d0/631d0)/18d0
  c = 0d0
  call quad(f09,a,b,cpole=c,epsabs=epsabs,epsrel=epsrel,verbose=.false.,result=result)
  call assert(result,true,"QAWC  1/(x*(5x^3+6)) Cauchy principal value of f(x).w(x)",tol=10*epsrel)




  a = 0d0
  b = 1d0
  epsabs=0d0
  epsrel=1d-12
  !
  !Gamma is Euler's constant.
  gamma =  0.5772156649d0
  !ci is the cosine integral:
  !ci(x) = integral ( x <= v < +oo ) - cos ( v ) / v dv.
  ci    = -0.001007d+00
  true = -0.1281316d0
  !
  integr= 2                  !1: w(x)=cos(omega*x);2: w(x)=sin(omega*x)
  omega = pi*10d0
  !
  call quad(f_log,a,b,omega=omega,weight_func=integr,epsabs=epsabs,epsrel=epsrel,verbose=.false.,result=result)
  call assert(result,true,"QAWO  log(x)*cos(w.x) from a to b",tol=1d-5)







  a = 0d0
  b = 1d0
  epsabs=0d0
  epsrel=1d-6
  !
  ci =  0.33740392290096813466d0
  si =  0.94608307036718301494d0
  true =  0.5d0*(ci*sin(1d0)+(pi/2d0 - si)*cos(1d0) - 1d0 ) 
  !
  !    = 1  (x-a)**alfa*(b-x)**beta
  !    = 2  (x-a)**alfa*(b-x)**beta*log(x-a)
  !    = 3  (x-a)**alfa*(b-x)**beta*log(b-x)
  !    = 4  (x-a)**alfa*(b-x)**beta*log(x-a)*log(b-x)
  integr= 2
  alfa = 0d0
  beta = 0d0
  call quad(f08,a,b,alfa=alfa,beta=beta,weight_func=integr,epsabs=epsabs,epsrel=epsrel,verbose=.false.,result=result)
  call assert(result,true,"QAWS 1/(1+LOG(X)**2)**2*w(x)  singularities at endpoint",tol=1d-4)





  a=0d0
  b=1d0
  epsabs=0d0
  epsrel=1d-4
  true=-4d0/9d0
  call quad(f_sqrt_log,a,b,epsabs=epsabs,epsrel=epsrel,verbose=.false.,result=result)
  call assert(result,true,"QNG  log(x)*sqrt(x) non-adaptive automatic",tol=1d-5)



contains

  function fsqrt(x)
    real(8) :: x
    real(8) :: fsqrt
    fsqrt = sqrt(x)
  end function fsqrt


  !! F01 is the integrand function SQRT(X) * LOG(X).
  function f_sqrt_log ( x ) result(f01)
    real(8) :: f01
    real(8) :: x
    if ( x <= 0d0 )then
       f01 = 0d0
    else
       f01 = sqrt ( x ) * log ( x )
    end if

    return
  end function f_sqrt_log


  !! F02 is the integrand function COS(100*SIN(X)).
  function f_cos100sin ( x ) result(f02)
    real(8) ::  f02
    real(8) ::  x
    f02 = cos ( 100d0 * sin ( x ) )
    return
  end function f_cos100sin


  !! F03 is the integrand function LOG(X)/SQRT(X).
  function f_log_over_sqrt ( x ) result(f03)
    real(8) ::  f03
    real(8) ::  x
    if ( x <= 0d0 ) then
       f03 = 0d0
    else
       f03 = log ( x ) / sqrt ( x )
    end if
    return
  end function f_log_over_sqrt


  !! F04 is the integrand function X^3 LOG((X^2-1)*(X^2-2))  
  function f_x3_log ( x ) result(f04)
    real(8) ::  f04
    real(8) ::  x
    f04 = x**3 * log ( abs ( ( x**2 - 1d0 ) * ( x**2 - 2d0 ) ) )
    return
  end function f_x3_log

  !! F05 is the integrand function LOG(X)/(1+100X^2).
  function f_log_over_px ( x ) result(f05)
    real(8) ::  f05
    real(8) ::  x
    f05 = log ( x ) / ( 1d0 + 100d0 * x**2 )
    return
  end function f_log_over_px

  !! F06 is the integrand function LOG(X).  
  function f_log ( x ) result(f06)
    real(8) ::  f06
    real(8) ::  x
    if ( x <= 0d0 ) then
       f06 = 0d0
    else
       f06 = log ( x )
    end if
    return
  end function f_log

  !! F07 is the integrand function 1/SQRT(X).  
  function f_1_over_sqrt ( x ) result(f07)
    real(8) ::  f07
    real(8) ::  x
    if ( x <= 0d0 ) then
       f07 = 0d0
    else
       f07 = 1d0 / sqrt ( x )
    end if
    return
  end function f_1_over_sqrt

  !! F08 is the integrand function 1/(1+LOG(X)**2)**2
  function f08 ( x ) 
    real(8) ::  f08
    real(8) ::  x
    if ( 0d0 < x ) then
       f08 = 1d0 / ( 1d0 + log ( x )**2 )**2
    else
       f08 = 0d0
    end if
    return
  end function f08

  !! F09 is the integrand function 1 / ( 5 X^3 + 6 ).
  function f09 ( x ) 
    real(8) ::  f09
    real(8) ::  x
    f09 = 1d0 / ( 5d0 * x**3 + 6d0 )
    return
  end function f09




end program test_SF_QUAD
