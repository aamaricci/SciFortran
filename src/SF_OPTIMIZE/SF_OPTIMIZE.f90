MODULE SF_OPTIMIZE
  USE OPTIMIZE_ROOT_FINDING
  USE OPTIMIZE_MINIMIZE
  private

  !OPTIMIZATION:
  !
  !General purpose
  public :: fmin                !Nelder-Mead
  public :: fmin_cg             !Conjugate-Gradient 1
  public :: fmin_cgplus         !Conjugate-Gradient 2
  public :: fmin_cgminimize     !Conjugate-Gradient 3 (very old f77)
  !
  public :: leastsq             
  public :: curvefit            
  !
  !Constrained (multivariate)
  public :: fmin_bfgs           !BFGS (constrained and not)
  !Global
  ! 
  !Scalar function minimizers
  public :: brent
  public :: dbrent
  public :: bracket


  !ROOT FINDING:
  !Scalar functions
  public :: brentq
  public :: bisect
  public :: newton
  public :: fzero

  !Multidimensional
  !General nonlinear solvers:
  public :: fsolve              !
  public :: broyden1
  !Large-scale nonlinear solvers:


END MODULE SF_OPTIMIZE
