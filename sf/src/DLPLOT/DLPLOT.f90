!###############################################################
!     PROGRAM  : DLPLOT
!     TYPE     : Module
!     PURPOSE  : PLOTTING LIBRARY FOR FORTRAN 90/95
!     AUTHORS  : Adriano Amaricci (Rutgers/SISSA) .AND. &
!                Cedric Weber (Rutgers)
!     DOC      : cf. the related man page (under development)
!###############################################################
module DLPLOT
  implicit none
  private

  public :: dplot,plot_dislin,plot_xmgrace
  public :: plot_dislin_movie
  public :: plot_dislin_VF
  public :: plot_dislin_3D,plot3Dsurface,plot3Dintensity
  public :: plot_dislin_3D_movie,plot_3D_surface_movie

  interface dplot
     module procedure D_plot, Z_plot, X_plot
  end interface dplot

  interface plot_dislin_movie
     module procedure plot_dislin2Dmovie_, plot_dislin2Dmovie__
  end interface plot_dislin_movie

  interface plot_dislin_VF
     module procedure plot_dislinVF_,plot_dislinVF__
  end interface plot_dislin_VF

  interface plot_dislin_3D
     module procedure R_plot_dislin3D, C_plot_dislin3D !
  end interface plot_dislin_3D

  interface dumpxmgrace
     module procedure dumpxmgrace_,dumpxmgrace__
  end interface dumpxmgrace
contains


  !*******************************************************************
  ! 2D - ELABORATED PLOT ROUTINES
  !*******************************************************************
  include "PLOT.f90"
  !*******************************************************************



  !*******************************************************************
  ! 2D - PLOT 2D MOVIES
  !*******************************************************************
  include "PLOT_MOVIE.f90"
  !*******************************************************************



  !********************************************************************
  ! 2D VECTOR FIELD - PLOT ROUTINES
  !********************************************************************
  include "PLOT_VF.f90"
  !********************************************************************



  !********************************************************************
  ! 3D - PLOT ROUTINES
  !********************************************************************
  include "PLOT_3D.f90"
  !********************************************************************



  !********************************************************************
  ! 3D - PLOT MOVIE ROUTINES
  !********************************************************************
  include "PLOT_3D_MOVIE.f90"
  !********************************************************************



  !********************************************************************
  ! 2N - XMGRACE WRAPPER ROUTINES
  !********************************************************************
  include "XMGRACE.f90"
  !********************************************************************




  !+------------------------------------------------------------------+
  !PROGRAM  : DRANRD
  !TYPE     : function
  !PURPOSE  : Get a random number
  !COMMENTS : you have call rand_init(iseed) to initialize the 
  !random seed, not mandatory if you're not interestedin getting 
  !different sequences.
  !+------------------------------------------------------------------+
  function drand()
    implicit none
    real(8) :: drand
    real(4) :: r
    call random_number(r)
    drand=dble(r)
  end function drand
  !*******************************************************************
  !*******************************************************************
  !*******************************************************************

end module DLPLOT
