MODULE SF_OPTIMIZE
  USE OPTIMIZE_ROOT_FINDING
  USE OPTIMIZE_MINIMIZE
  private

  !ROOT FINDING:
  public :: fzero_brentq
  public :: zbrent

  public :: fsolve
  public :: fzero_hybrd

  public :: fzero_broyden
  public :: broydn            !backward compatibility

  !MINIMIZE
  public :: fmin_cg
  public :: fmin_cgplus
  public :: fmin_cgminimize

END MODULE SF_OPTIMIZE
