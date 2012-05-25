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

  public :: plot_dislin
  public :: plot_xmgrace
  public :: dplot
  public :: dplot_3d
  public :: dplot_3d_surface
  public :: dplot_3d_intensity
  public :: dplot_3d_surface_rotating
  public :: dplot_vector_field



  public :: dplot_movie
  public :: dplot_3d_intensity_animated
  public :: dplot_3d_surface_animated



  !I am lazy and I don't want to fix names of routines so I add 1000 interfaces here:
  interface dumpxmgrace
     module procedure dumpxmgrace_,dumpxmgrace__
  end interface dumpxmgrace

  interface dplot
     module procedure D_plot, Z_plot, X_plot
  end interface dplot

  interface dplot_vector_field
     module procedure plot_dislinVF_,plot_dislinVF__
  end interface dplot_vector_field

  interface dplot_3d
     module procedure R_plot_dislin3D, C_plot_dislin3D !
  end interface dplot_3D

  interface dplot_3d_surface
     module procedure plot3Dsurface
  end interface dplot_3d_surface

  interface dplot_3d_intensity
     module procedure plot3Dintensity
  end interface dplot_3d_intensity

  interface dplot_3d_surface_rotating
     module procedure plot3Dsurface_rotating
  end interface dplot_3d_surface_rotating

  interface dplot_movie
     module procedure plot_dislin2Dmovie_, plot_dislin2Dmovie__
  end interface dplot_movie

  interface dplot_3d_intensity_animated
     module procedure plot_dislin_3D_movie
  end interface dplot_3d_intensity_animated

  interface dplot_3d_surface_animated
     module procedure plot_3D_surface_movie
  end interface dplot_3d_surface_animated



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
