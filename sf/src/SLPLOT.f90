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
          data_saveV_I,data_saveV_R,data_saveV_C,&
          data_saveM_I,data_saveM_R,data_saveM_C,&
          data_saveA3_I,data_saveA3_R,data_saveA3_C,&
          d_splot3d,c_splot3d,d_splot3d_animate,c_splot3d_animate,splot3D__
  end interface splot

  public :: splot

contains

  ! 0-dim array
  ! X=int  ; Y=int,dble,cmplx
  ! X=dble ; Y=int,dble,cmplx
  include "splot_P.f90"

  ! 1-dim array
  ! X=int  ; Y=int,dble,cmplx
  ! X=dble ; Y=int,dble,cmplx
  include "splot_V.f90"

  ! 1,2-dim arrays
  ! Y=int,dble,cmplx
  ! X=int,dble [only for 2-dim, optional]
  include "data_save.f90"


  ! 3Dplot:
  include "splot_3d.f90"


end module SLPLOT
