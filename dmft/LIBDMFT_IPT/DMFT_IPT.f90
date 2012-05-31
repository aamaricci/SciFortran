!###############################################################
!     PROGRAM  : IPT_LIB
!     TYPE     : Module
!     PURPOSE  : Contains all the modules used to solve DMFT with IPT
!     AUTHORS  : Adriano Amaricci
!###############################################################
module DMFT_IPT
  !GLOBAL VARIABLES:
  USE IPT_VARS_GLOBAL
  !IPT:
  USE IPT_SOPT
  USE IPT_MATS
  USE IPT_KELDYSH
  !IPT+SC:
  USE IPT_SC_SOPT
  USE IPT_SC_MATS
  !IPT+AF:
  USE IPT_AF_MATS
  implicit none
end module DMFT_IPT
