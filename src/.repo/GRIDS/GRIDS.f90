module GRIDS
  !###############################################################
  !     PROGRAM  : GRIDS
  !     TYPE     : Module
  !     PURPOSE  : Setup some grids useful in calculations
  !###############################################################
  implicit none
  private

  !Grid arrays:
  !=========================================================
  REAL(8),DIMENSION(:),ALLOCATABLE,PUBLIC :: wr,wm,t,tau

  public :: init_wgrid
  public :: init_wmgrid
  public :: init_tgrid
  public :: init_taugrid
contains
  !+----------------------------------------------------------------+
  !PROGRAM  : init_wgrid
  !TYPE     : subroutine
  !PURPOSE  : 
  !+----------------------------------------------------------------+
  subroutine init_wgrid(wr_,dw,M)
    integer             :: i,M
    real(8)             :: dw
    real(8),allocatable :: wr_(:)
    allocate(wr_(-M:M))
    forall(i=-M:M)wr_(i)=dble(i)*dw
  end subroutine init_wgrid
  !******************************************************************
  !******************************************************************
  !******************************************************************



  !+----------------------------------------------------------------+
  !PROGRAM  : init_wmgrid
  !TYPE     : subroutine
  !PURPOSE  : 
  !+----------------------------------------------------------------+
  subroutine init_wmgrid(wm_,beta,M)
    integer             :: i,M
    real(8)             :: beta,pi
    real(8),allocatable :: wm_(:)
    pi=3.14159265358979d0
    allocate(wm_(M))!;M=size(wm_)
    forall(i=1:M)wm_(i)=pi*dble(2*i-1)/beta
  end subroutine init_wmgrid
  !******************************************************************
  !******************************************************************
  !******************************************************************



  !+-----------------------------------------------------------------+
  !PROGRAM  : init_tgrid
  !TYPE     : Subroutine
  !PURPOSE  : initialize the real time axis grid of length given
  !+-----------------------------------------------------------------+
  subroutine init_tgrid(t_,dt,M)
    integer             :: i,M
    real(8)             :: dt
    real(8),allocatable :: t_(:)   
    allocate(t_(-M:M))
    forall(i=-M:M)t_(i)=dble(i)*dt
  end subroutine init_tgrid
  !******************************************************************
  !******************************************************************
  !******************************************************************



  !+-----------------------------------------------------------------+
  !PROGRAM  : init_taugrid
  !TYPE     : Subroutine
  !PURPOSE  : initialize the real time axis grid of length given
  !+-----------------------------------------------------------------+
  subroutine init_taugrid(tau_,dtau,M)
    integer              :: i,M,N
    real(8)              :: dtau
    real(8),allocatable  :: tau_(:)
    allocate(tau_(0:M))!;M=size(tau_)-1
    forall(i=0:M)tau_(i)=dble(i)*dtau
  end subroutine init_taugrid
  !******************************************************************
  !******************************************************************
  !******************************************************************



END module GRIDS
