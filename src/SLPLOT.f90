!###############################################################
!     PROGRAM  : SLPLOT
!     PURPOSE  : VERY SIMPLE PLOTTING LIBRARY
!     AUTHORS  : A.Amaricci (SISSA)
!###############################################################
module SLPLOT
  USE IOFILE
  USE COMMON_VARS
  implicit none
  private

  interface splot
     module procedure &
          splotP_II,splotP_IR,splotP_IC, &
          splotP_RI,splotP_RR,splotP_RC, &
          splotV_II, splotV_IR,splotV_IC,&
          splotV_RI, splotV_RR,splotV_RC,&
          splotM_II,splotM_IR,splotM_IC,&
          splotM_RI,splotM_RR,splotM_RC,&
          splotA3_II,splotA3_IR,splotA3_IC,&
          splotA3_RI,splotA3_RR,splotA3_RC!,&
          ! data_saveV_I,data_saveV_R,data_saveV_C,&
          ! data_saveM_I,data_saveM_R,data_saveM_C,&
          ! data_saveA3_I,data_saveA3_R,data_saveA3_C
  end interface splot


  interface splot3d
     module procedure &
          d_splot3d,c_splot3d,d_splot3d_animate,c_splot3d_animate
  end interface splot3d


  interface store_data
     module procedure &
          data_saveV_I,data_saveV_R,data_saveV_C,&
          data_saveM_I,data_saveM_R,data_saveM_C,&
          data_saveA3_I,data_saveA3_R,data_saveA3_C
  end interface store_data

  public :: splot
  public :: splot3d
  public :: store_data


contains

  ! 0-dim array
  include "slplot_splot_P.f90"

  ! 1-dim array
  include "slplot_splot_V.f90"

  ! N=2,3-dim array
  include "slplot_splot_M.f90"

  ! 3Dplot:
  include "slplot_splot_3d.f90"

  ! STORE arrays
  include "slplot_data_save.f90"

end module SLPLOT
