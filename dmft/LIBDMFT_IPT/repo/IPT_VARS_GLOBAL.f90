!###############################################################
!     PROGRAM  : VARS_GLOBAL
!     TYPE     : Module
!     PURPOSE  : Contains global variables
!     AUTHORS  : Adriano Amaricci
!###############################################################
module IPT_VARS_GLOBAL
  !LIBRARY: put here all the LIBS in common
  USE COMMON_VARS
  USE GREENFUNX
  USE CHRONOBAR
  USE SLPLOT
  USE FFT_MKL
  USE INTEGRATE
  USE SPLINE
  USE TOOLS
  USE GRIDS
  implicit none

  !Size of the problem: some parameters:
  !=========================================================
  integer     :: L !once it was protected


  !local variables, others are defined in COMMON_VARS 
  !=========================================================
  integer               :: i,j,Nx
  logical               :: printf
  real(8) 	        :: wmin,wmax
  real(8)               :: eps_error
  integer               :: Nsuccess
  real(8)               :: weigth
  real(8)               :: deltasc
  character(len=64)     :: label

  !Namelists:
  !=========================================================
  namelist/variables/&
       L,        &
       beta,     &
       U,        &
       ts,tsp,   &
       xmu,      &
       wmax,     &
       Nx,       &
       nloop,    &
       eps,      &
       printf,   &
       weigth,   &
       eps_error,&
       Nsuccess, &
       deltasc,  &
       label

contains

  !+----------------------------------------------------------------+
  !PROGRAM  : READinput
  !TYPE     : subroutine
  !PURPOSE  : Read input file
  !+----------------------------------------------------------------+
  subroutine read_input(inputFILE)
    character(len=*) :: inputFILE
    integer :: i
    !variables: default values
    U     = 2.d0
    beta  = 100.d0
    ts    = 0.5d0
    tsp   = 0.d0
    xmu   = 0.d0
    Nx    = 20
    nloop = 10
    eps   = 0.01d0
    wmax  = 5.d0
    L     = 2048
    printf= .true.
    weigth= 0.9d0
    eps_error= 1.d-4
    Nsuccess = 2
    deltasc  = 0.1d0
    label    = ""

    n_command_arg=command_argument_count()
    open(10,file=adjustl(trim(inputFILE)))
    read(10,nml=variables)
    close(10)

    !Process command line variable change:
    if(n_command_arg/=0)then
       do i=1,n_command_arg
          call get_command_argument(i,arg_buffer)
          pos      = scan(arg_buffer,"=")
          nml_name = arg_buffer(1:pos-1);call s_cap(nml_name)
          nml_value= arg_buffer(pos+1:)
          select case(nml_name)
          case("U")     ;read(nml_value,*)U
          case("BETA")  ;read(nml_value,*)beta
          case("TS")    ;read(nml_value,*)ts
          case("TSP")   ;read(nml_value,*)tsp
          case("XMU")   ;read(nml_value,*)xmu
          case("NX")    ;read(nml_value,*)Nx
          case("NLOOP") ;read(nml_value,*)nloop
          case("L")     ;read(nml_value,*)L
          case("EPS")   ;read(nml_value,*)eps
          case("WMAX")  ;read(nml_value,*)wmax
          case("PRINTF");read(nml_value,*)printf
          case("EPS_ERROR");read(nml_value,*)eps_error
          case("NSUCCESS");read(nml_value,*)Nsuccess
          case("DELTASC");read(nml_value,*)deltasc
          case("LABEL");read(nml_value,*)label
          case default
             print*,"No corresponging variable in NML"
          end select
       enddo
    endif

    write(*,nml=variables)
    open(10,file="used."//adjustl(trim(inputFILE)))
    write(10,nml=variables)
    close(10)

    D     = 2.d0*ts
    wmin  =-wmax
    fmesh = abs(wmax-wmin)/dble(2*L) 
    dt    = pi/fmesh/dble(L)
    dtau  = beta/dble(L)
    call init_wmgrid(wm,beta,L)
    call init_taugrid(tau,dtau,L)
    call init_wgrid(wr,fmesh,L)
    call init_tgrid(t,dt,L)
    return
  end subroutine read_input
  !******************************************************************
  !******************************************************************
  !******************************************************************

end module IPT_VARS_GLOBAL
