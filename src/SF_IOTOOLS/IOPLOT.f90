module IOPLOT
  USE IOFILE
  implicit none
  private

  interface splot
     module procedure :: splotA1_RR
     module procedure :: splotA1_RC
     module procedure :: splotA2_RR
     module procedure :: splotA2_RC
     module procedure :: splotA3_RR
     module procedure :: splotA3_RC
     module procedure :: splotA4_RR
     module procedure :: splotA4_RC
     module procedure :: splotA5_RR
     module procedure :: splotA5_RC
     module procedure :: splotA6_RR
     module procedure :: splotA6_RC
     module procedure :: splotA7_RR
     module procedure :: splotA7_RC
  end interface splot



  interface store_data
     module procedure :: data_saveA1_R
     module procedure :: data_saveA1_C
     module procedure :: data_saveA2_R
     module procedure :: data_saveA2_C
     module procedure :: data_saveA3_R
     module procedure :: data_saveA3_C
     module procedure :: data_saveA4_R
     module procedure :: data_saveA4_C
     module procedure :: data_saveA5_R
     module procedure :: data_saveA5_C
     module procedure :: data_saveA6_R
     module procedure :: data_saveA6_C
     module procedure :: data_saveA7_R
     module procedure :: data_saveA7_C
  end interface store_data


  interface splot3d
     module procedure :: d_splot3d
     module procedure :: c_splot3d
     module procedure :: d_splot3d_animate
     module procedure :: c_splot3d_animate
  end interface splot3d


  public :: splot
  public :: store_data
  public :: splot3d


  integer            :: unit
  character(len=128) :: fmt


contains


  ! SPLOT ararys (1--7)
  include "ioplot_splot.f90"


  ! STORE_DATA arrays (1--7)
  include "ioplot_store_data.f90"


  ! 3Dplot:
  include "ioplot_splot3d.f90"

end module IOPLOT
