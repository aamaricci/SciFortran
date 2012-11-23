!###############################################################
!     PROGRAM  : VARS_GLOBAL
!     TYPE     : Module
!     PURPOSE  : Contains global variables
!     AUTHORS  : Adriano Amaricci
!###############################################################
module IPT_VARS_GLOBAL
  !LIBRARY: put here all the LIBS in common
  USE COMMON_VARS
  USE PARSE_CMD
  USE GREENFUNX
  USE FFTGF
  USE INTEGRATE
  USE TOOLS
  implicit none

  !Size of the problem: some parameters:
  !=========================================================
  integer :: L              !a large number fix the size of the problem
  integer :: nloop    !dmft loop variables
  real(8) :: d              !bandwidth
  real(8) :: ts,tsp,tpp,tdd !n.n./n.n.n. hopping amplitude
  real(8) :: u,v            !local,non-local interaction
  real(8) :: tpd,vpd        !hybridization,band-band coupling
  real(8) :: ed0,ep0        !orbital energies
  real(8) :: xmu            !chemical potential
  real(8) :: dt,dtau        !time step
  real(8) :: fmesh          !freq. step
  real(8) :: beta           !inverse temperature
  real(8) :: eps            !broadening
  integer :: Nx
  logical :: printf
  real(8) :: wmin,wmax
  real(8) :: eps_error
  integer :: Nsuccess
  real(8) :: weight
  real(8) :: deltasc
  real(8) :: nread,nerror,ndelta

  !Namelists:
  !=========================================================
  namelist/variables/&
       L,        &
       beta,     &
       U,        &
       ts,tsp,   &
       tdd,tpp,  &
       xmu,      &
       wmax,     &
       Nx,       &
       nloop,    &
       eps,      &
       printf,   &
       weight,   &
       eps_error,&
       Nsuccess, &
       deltasc,&
       nread,nerror,ndelta


contains

  !+----------------------------------------------------------------+
  !PROGRAM  : READinput
  !TYPE     : subroutine
  !PURPOSE  : Read input file
  !+----------------------------------------------------------------+
  subroutine read_input(inputFILE)
    character(len=*) :: inputFILE
    integer :: i
    logical :: back,control
    character(len=256),allocatable,dimension(:) :: help_buffer
    !variables: default values
    U     = 2.d0
    beta  = 100.d0
    ts    = 0.5d0
    tsp   = 0.d0
    tdd   = 0.d0
    tpp   = 0.d0
    vpd=0.4d0
    ed0=0.d0
    ep0=0.d0
    xmu   = 0.d0
    Nx    = 20
    nloop = 10
    eps   = 0.01d0
    wmax  = 5.d0
    L     = 2048
    printf= .true.
    weight= 0.9d0
    eps_error= 1.d-4
    Nsuccess = 2
    deltasc  = 0.1d0
    nread=0.d0
    nerror=1.d-4
    ndelta=0.1d0


    allocate(help_buffer(39))
    help_buffer=([character(len=512)::&
         'NAME',&
         '  library:dmft_ipt',&
         '  ',&
         'DESCRIPTION',&
         '  provide access to a set of solvers for the dmft based on IPT',&
         '  ',&
         'OPTIONS',&
         '  u=[2]      -- interaction',&
         '  beta=[100] -- inverse temperature',&
         '  xmu=[0]    -- chemical potential',&
         '  nloop=[10] -- max number of iterations',&
         '  L=[2048]   -- number of frequencies',&
         '  printf=[T] -- print flag',&
         '  ts=[0.5]   -- n.n. hopping parameter',&
         '  tsp=[0]    -- n.n.n. hopping parameter',&
         '  nx=[20]    -- number of points in energy/k-grid',&
         '  wmax=[5]   -- max frequency on real axis',&
         '  eps=[0.01] -- broadening parameter',&
         '  deltasc=[0.1]     -- breaking symmetry parameter',&
         '  nread=[0.0]       -- required density look for mu',&
         '  nerror=[1.d-4]    -- error treshold for mu-loop',&
         '  ndelta=[0.1]      -- mu-step in mu-loop for fixed density',&
         '  eps_error=[1.D-4] -- error treshold',&
         '  success=[2]       -- number of converged events',&
         '  ',&
         'ROUTINES',&
         '  normal:',&
         '  solve_ipt_keldysh(fg0_,wr_,t_) result(sigma_) dim(-L:L) ',&
         '  solve_ipt_matsubara(fg0_)      result(sigma_) dim(L)',&
         '  solve_ipt_sopt(fg0_,wr_)       result(sigma_) dim(-L:L)',&
         '  solve_ipt_zeroT(fg0_,dw_,dt_)  result(sigma_) dim(-L:L)',&
         '  solve_mpt_sopt(fg0_,wr_,n_,n0_,xmu0_)  result(sigma_) dim(-L:L)',&
         '  solve_mpt_matsubara(fg0_,n_,n0_,xmu0_) result(sigma_) dim(L)',&
         '  ',&
         '  supercond:',&
         '  solve_ipt_sc_sopt(fg0_,wr_,delta_)  result(sigma_) dim(2,-L:L)',&
         '  solve_ipt_sc_matsubara(fg0_,delta_) result(sigma_) dim(2,L)',&
         '  solve_mpt_sc_sopt(fg0_,wr_,n_,n0_,delta_,delta0_)  result(sigma_) dim(2,-L:L)',&
         '  solve_mpt_sc_matsubara(fg0_,n_,n0_,delta_,delta0_) result(sigma_) dim(2,L)',&
         ' ',&
         '  '])

    call parse_cmd_help(help_buffer,status=back)
    deallocate(help_buffer) ; 
    if(back)return

    inquire(file=adjustl(trim(inputFILE)),exist=control)
    if(control)then
       open(10,file=adjustl(trim(inputFILE)))
       read(10,nml=variables)
       close(10)
    else
       open(10,file="default."//adjustl(trim(inputFILE)))
       write(10,nml=variables)
       close(10)
       call error("can not open INPUT file, dumping a default version in default."//adjustl(trim(inputFILE)))
    endif

    call parse_cmd_variable(u,"U")
    call parse_cmd_variable(beta,"BETA")
    call parse_cmd_variable(ts,"TS")
    call parse_cmd_variable(tsp,"TSP")
    call parse_cmd_variable(xmu,"XMU")
    call parse_cmd_variable(vpd,"VPD")
    call parse_cmd_variable(ed0,"ED0")
    call parse_cmd_variable(ep0,"EP0")
    call parse_cmd_variable(nx,"NX")
    call parse_cmd_variable(nloop,"NLOOP")
    call parse_cmd_variable(L,"L")
    call parse_cmd_variable(eps,"EPS")
    call parse_cmd_variable(wmax,"WMAX")
    call parse_cmd_variable(weight,"WEIGTH","WEIGHT")
    call parse_cmd_variable(printf,"PRINTF")
    call parse_cmd_variable(eps_error,"EPS_ERROR")
    call parse_cmd_variable(nsuccess,"NSUCCESS")
    call parse_cmd_variable(deltasc,"DELTASC")
    call parse_cmd_variable(nread,"NREAD")
    call parse_cmd_variable(nerror,"NERROR")
    call parse_cmd_variable(ndelta,"NDELTA")

    if(mpiID==0)then
       write(*,nml=variables)
       open(10,file="used."//adjustl(trim(inputFILE)))
       write(10,nml=variables)
       close(10)
    endif
    return
  end subroutine read_input
  !******************************************************************
  !******************************************************************
  !******************************************************************




end module IPT_VARS_GLOBAL
