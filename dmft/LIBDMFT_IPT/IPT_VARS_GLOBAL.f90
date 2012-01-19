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
  USE FFTGF
  USE INTEGRATE
  USE TOOLS
  implicit none

  !Size of the problem: some parameters:
  !=========================================================
  integer     :: L !once it was protected


  !local variables, others are defined in COMMON_VARS 
  !=========================================================
  integer               :: Nx
  logical               :: printf
  real(8) 	        :: wmin,wmax
  real(8)               :: eps_error
  integer               :: Nsuccess
  real(8)               :: weigth
  real(8)               :: deltasc
  character(len=64)     :: label

  !real(8),allocatable,dimension(:) :: wr,wm,t,tau

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
    logical :: back

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

    open(10,file=adjustl(trim(inputFILE)))
    read(10,nml=variables)
    close(10)

    allocate(help_buffer(39))
    help_buffer=([&
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
         '  label=[""] -- a useful label',&
         '  eps=[0.01] -- broadening parameter',&
         '  deltasc=[0.1]     -- breaking symmetry parameter',&
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
         '  solve_mpt_sc_sopt(fg0_,wr_,n_,n0_,delta_,delta0_)  result(sigma_) dim(2,-L:L) *not working*',&
         '  solve_mpt_sc_matsubara(fg0_,n_,n0_,delta_,delta0_) result(sigma_) dim(2,L)',&
         ' ',&
         '  '])

    call parse_cmd_help(help_buffer,status=back)
    deallocate(help_buffer)
    if(back)return

    call parse_cmd_variable(u,"U")
    call parse_cmd_variable(beta,"BETA")
    call parse_cmd_variable(ts,"TS")
    call parse_cmd_variable(tsp,"TSP")
    call parse_cmd_variable(xmu,"XMU")
    call parse_cmd_variable(nx,"NX")
    call parse_cmd_variable(nloop,"NLOOP")
    call parse_cmd_variable(L,"L")
    call parse_cmd_variable(eps,"EPS")
    call parse_cmd_variable(wmax,"WMAX")
    call parse_cmd_variable(weigth,"WEIGTH","WEIGHT")
    call parse_cmd_variable(printf,"PRINTF")
    call parse_cmd_variable(eps_error,"EPS_ERROR")
    call parse_cmd_variable(nsuccess,"NSUCCESS")
    call parse_cmd_variable(deltasc,"DELTASC")
    call parse_cmd_variable(label,"LABEL")


    ! !Process command line variable change:
    ! do i=1,command_argument_count()
    !    nml_var=get_cmd_variable(i)
    !    select case(nml_var%name)
    !    case("U")     ;read(nml_var%value,*)U
    !    case("BETA")  ;read(nml_var%value,*)beta
    !    case("TS")    ;read(nml_var%value,*)ts
    !    case("TSP")   ;read(nml_var%value,*)tsp
    !    case("XMU")   ;read(nml_var%value,*)xmu
    !    case("NX")    ;read(nml_var%value,*)Nx
    !    case("NLOOP") ;read(nml_var%value,*)nloop
    !    case("L")     ;read(nml_var%value,*)L
    !    case("EPS")   ;read(nml_var%value,*)eps
    !    case("WMAX")  ;read(nml_var%value,*)wmax
    !    case("PRINTF");read(nml_var%value,*)printf
    !    case("EPS_ERROR");read(nml_var%value,*)eps_error
    !    case("NSUCCESS");read(nml_var%value,*)Nsuccess
    !    case("DELTASC");read(nml_var%value,*)deltasc
    !    case("LABEL");read(nml_var%value,*)label
    !    case default
    !       print*,"No corresponging variable in NML"
    !    end select
    ! enddo

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
