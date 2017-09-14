MODULE SF_OPTIMIZE
  USE OPTIMIZE_ROOT_FINDING
  USE OPTIMIZE_MINIMIZE
  private

  !OPTIMIZATION:
  !
  !General purpose
  public :: fmin_cg
  public :: fmin_cgplus
  public :: fmin_cgminimize
  !
  !Constrained (multivariate)
  !
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
  public :: fsolve
  public :: broyden1
  !Large-scale nonlinear solvers:


END MODULE SF_OPTIMIZE
