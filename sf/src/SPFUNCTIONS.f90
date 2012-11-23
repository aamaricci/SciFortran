  !###############################################################
  ! PROGRAM  : SPFUNCTIONS
  ! PURPOSE  : Special functions
  !###############################################################  
module SPFUNCTIONS
  USE COMMON_VARS
  implicit none
  private

  !FUNCTIONS:
  public :: heaviside
  public :: step
  public :: fermi
  public :: sgn
  public :: dens_hyperc

  !ERROR FUNCS
  public :: wfun         !complex error function (Faddeeva function)
  public :: zerf   

  !BETHE:
  public :: gfbethe
  public :: gfbether
  public :: bethe_lattice
  public :: dens_bethe

end module SPFUNCTIONS
