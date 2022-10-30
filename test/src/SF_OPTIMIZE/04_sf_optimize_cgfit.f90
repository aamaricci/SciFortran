program testCGFIT
  USE SCIFOR
  USE ASSERTING
  implicit none

  integer,parameter    :: npar=2
  integer              :: i
  real(8),dimension(2) :: x
  integer              :: niter
  real(8)              :: fresult

  x=-2.d0
  call fmin_cg(x,fcn,niter,fresult,ftol=1d-10)
  call assert([x,fresult],[1d0,1d0,0d0],"Fmin_CG no Grad Rosenbrack function:",tol=1d-5)


  !Simple test:
  x=-2.d0
  call fmin_cg(x,fcn,dfcn,niter,fresult,ftol=1d-10)
  call assert([x,fresult],[1d0,1d0,0d0],"Fmin_CG no Grad Rosenbrack function:",tol=1d-5)


contains


  function fcn(x)
    real(8),dimension(:) :: x
    real(8)              :: fcn
    fcn = 100d0*((x(2) - x(1)**2)**2) + (1d0 - x(1))**2
  end function fcn

  function dfcn(x)
    real(8),dimension(:)       :: x
    real(8),dimension(size(x)) :: dfcn
    dfcn(1) = 200d0*(x(2) - x(1)**2)*(-2*x(1)) - 2*(1 - x(1))
    dfcn(2) = 200d0*(x(2) - x(1)**2)
  end function dfcn


end program testCGFIT



