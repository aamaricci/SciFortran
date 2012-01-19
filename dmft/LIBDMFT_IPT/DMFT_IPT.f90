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
  !USE IPT_ZEROT
  !IPT+SC:
  USE IPT_SC_SOPT
  USE IPT_SC_MATS
  !MPT:
  !USE MPT_SOPT
  !USE MPT_MATS
  !MPT+SC:
  !USE MPT_SC_SOPT
  !USE MPT_SC_MATS
  implicit none
end module DMFT_IPT
