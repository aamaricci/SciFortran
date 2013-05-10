!###############################################################
!     PROGRAM  : SLREAD
!     PURPOSE  : VERY SIMPLE READING LIBRARY FOR FORTRAN 90/95
!     AUTHORS  : A.Amaricci (SISSA)
!###############################################################
module SLREAD
  USE IOFILE
  USE COMMON_VARS
  implicit none
  private
  logical           :: control

  interface sread
     module procedure &
          sreadP_II,sreadP_IR,sreadP_IC, &
          sreadP_RI,sreadP_RR,sreadP_RC, &
          sreadV_II, sreadV_IR,sreadV_IC,&
          sreadV_RI, sreadV_RR,sreadV_RC,&
          sreadM_II,sreadM_IR,sreadM_IC,&
          sreadM_RI,sreadM_RR,sreadM_RC,&
          sreadA3_II,sreadA3_IR,sreadA3_IC,&
          sreadA3_RI,sreadA3_RR,sreadA3_RC!,&
          ! data_readV_I,&
          ! data_readV_R,&
          ! data_readV_C,&
          ! data_readM_I,&
          ! data_readM_R,&
          ! data_readM_C
  end interface sread

  interface read_data
     module procedure &
          data_readV_I,&
          data_readV_R,&
          data_readV_C,&
          data_readM_I,&
          data_readM_R,&
          data_readM_C
  end interface read_data

  public :: sread
  public :: read_data

contains
  ! 0-dim array
  include "slread_sread_P.f90"

  ! 1-dim array
  include "slread_sread_V.f90"

  ! N=2,3-dim array
  include "slread_sread_M.f90"

  ! 1,2-dim arrays
  include "slread_data_read.f90"

end module SLREAD
