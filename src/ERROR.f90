!###############################################################
! PURPOSE :
!###############################################################  
module ERROR
  USE COMMON_VARS
  USE IOTOOLS, only: reg
  implicit none
  private


  !CONVERGENCE
  public :: check_convergence,check_convergence_global,check_convergence_local


  interface check_convergence
     module procedure &
          i0_check_convergence_relative,&
          i1_check_convergence_relative,&
          i2_check_convergence_relative,&
          d0_check_convergence_relative_,&
          d0_check_convergence_relative,&
          d1_check_convergence_relative,&
          d2_check_convergence_relative,&
          z0_check_convergence_relative,&
          z1_check_convergence_relative,&
          z2_check_convergence_relative
  end interface check_convergence


  interface check_convergence_local
     module procedure &
          i0_check_convergence_local,&
          i1_check_convergence_local,&
          i2_check_convergence_local,&
          d0_check_convergence_local,&
          d1_check_convergence_local,&
          d2_check_convergence_local,&
          z0_check_convergence_local,&
          z1_check_convergence_local,&
          z2_check_convergence_local
  end interface check_convergence_local

  interface check_convergence_global
     module procedure &
          i0_check_convergence_global,&
          i1_check_convergence_global,&
          i2_check_convergence_global,&
          d0_check_convergence_global,&
          d1_check_convergence_global,&
          d2_check_convergence_global,&
          z0_check_convergence_global,&
          z1_check_convergence_global,&
          z2_check_convergence_global
  end interface check_convergence_global



contains

  include "error_convergence_relative.f90"

  include "error_convergence_local.f90"

  include "error_convergence_global.f90"

end module ERROR
