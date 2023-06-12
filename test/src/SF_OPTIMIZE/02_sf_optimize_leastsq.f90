program testMINPACK
  USE SCIFOR
  USE ASSERTING
  implicit none
  integer, parameter   :: n = 2, npar=3, m=50
  integer              :: iflag
  integer              :: info
  real(8)              :: tol
  real(8)              :: x(n),x0(n),y(n),a(npar)
  real(8),dimension(m) :: xdata,ydata,ynoise,xwork,ywork,fvec



  xdata = linspace(0d0,4d0,size(xdata))
  ynoise = 0.1d0*ran_normal(size(xdata)) 
  ydata = exp_func(xdata,[2.5d0,1.3d0,0.5d0])    


  xwork = xdata
  ywork = ydata+ynoise
  a = [2.5d0,1.3d0,0.5d0]
  call leastsq(lmdif_func,a,m)
  ywork = exp_func(xwork,a)
  call assert(ywork,ydata,"LeastSq  LMDIF with Function interface:",tol=0.2d0)


  xwork = xdata
  ywork = ydata+ynoise
  a = [2.5d0,1.3d0,0.5d0]
  call leastsq(lmdif_sub,a,m)
  ywork = exp_func(xwork,a)
  call assert(ywork,ydata,"LeastSq  LMDIF with Subroutine interface:",tol=0.2d0)




  xwork = xdata
  ywork = ydata+ynoise
  a = [2.5d0,1.3d0,0.5d0]
  call leastsq(lmder_func,lmder_dfunc,a,m)
  ywork = exp_func(xwork,a)
  call assert(ywork,ydata,"LeastSq + Grad LMDER with Function interface:",tol=0.2d0)




  xwork = xdata
  ywork = ydata+ynoise
  a = [2.5d0,1.3d0,0.5d0]
  call leastsq(lmder_sub,lmder_dsub,a,m)
  ywork = exp_func(xwork,a)
  call assert(ywork,ydata,"LeastSq + Grad LMDER with Function interface:",tol=0.2d0)









contains




  function lmdif_func(a,m) result(ffunc)
    real(8),dimension(:) :: a
    integer              :: m
    real(8),dimension(m) :: ffunc
    ffunc = exp_func(xwork,a) - ywork
    ! ffunc = a(1)*exp(-a(2)*xdata(:)) + a(3) - ydata(:)
  end function lmdif_func


  subroutine lmdif_sub(a,m,ffunc)
    real(8),dimension(:) :: a
    integer              :: m
    real(8),dimension(m) :: ffunc
    ffunc = exp_func(xwork,a) - ywork
  end subroutine lmdif_sub





  function lmder_func(a,m) result(ffunc)
    real(8),dimension(:)           :: a
    integer                        :: m
    real(8),dimension(m) :: ffunc
    ffunc = exp_func(xwork,a) - ywork
  end function lmder_func
  function lmder_dfunc(a,m) result(dfunc)
    real(8),dimension(:)         :: a
    integer                      :: m
    real(8),dimension(m,size(a)) :: dfunc
    dfunc(1:,1) = exp(-a(2)*xwork(:))
    dfunc(1:,2) = -a(1)*xwork(:)*exp(-a(2)*xwork(:))
    dfunc(1:,3) = 1d0
  end function lmder_dfunc




  subroutine lmder_sub(a,m,ffunc)
    real(8),dimension(:)           :: a
    integer                        :: m
    real(8),dimension(m) :: ffunc
    ffunc = exp_func(xwork,a) - ywork
  end subroutine lmder_sub
  subroutine lmder_dsub(a,m,dfunc)
    real(8),dimension(:)         :: a
    integer                      :: m
    real(8),dimension(m,size(a)) :: dfunc
    dfunc(1:,1) = exp(-a(2)*xwork(:))
    dfunc(1:,2) = -a(1)*xwork(:)*exp(-a(2)*xwork(:))
    dfunc(1:,3) = 1d0
  end subroutine lmder_dsub





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



  function ran_normal(n)
    integer              :: n
    real(8),dimension(n) :: ran_normal
    integer              :: i
    do i=1,n
       ran_normal(i) = random_normal()
    end do
  end function ran_normal




end program testMINPACK
