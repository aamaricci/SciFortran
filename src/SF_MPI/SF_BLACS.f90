MODULE SF_BLACS
  implicit none
  private

  integer :: p_rank=0
  integer :: p_size=1
  integer :: p_Nx=0
  integer :: p_Ny=0
  integer :: p_context=0
  logical :: blacs_status=.false.

  public :: p_rank
  public :: p_size
  public :: p_Nx,p_Ny
  public :: p_context
  public :: blacs_status


#ifdef _SCALAPACK
  public :: init_BLACS
  public :: finalize_BLACS
  public :: get_master_BLACS
  public :: get_rank_BLACS
  public :: get_size_BLACS
  !


contains

  subroutine init_BLACS()
    integer :: i
    !< Initialize BLACS processor grid (like MPI)
    call blacs_setup(p_rank,p_size)  ![id, size]
    !
    do i=1,int( sqrt( dble(p_size) ) + 1 )
       if(mod(p_size,i)==0) p_Nx = i
    end do
    p_Ny = p_size/p_Nx
    !
    !< Init context with p_Nx,p_Ny procs
    call sl_init(p_context,p_Nx,p_Ny)
    !
    blacs_status=.true.
  end subroutine init_BLACS


  subroutine finalize_BLACS(blacs_end)
    integer,optional :: blacs_end
    integer :: blacs_end_
    blacs_end_=0;if(present(blacs_end))blacs_end_=blacs_end
    call blacs_gridexit(p_context)
    call blacs_exit(blacs_end_)
    blacs_status=.false.
  end subroutine finalize_BLACS


  function get_master_BLACS() result(master)
    logical :: master
    if(.not.blacs_status)stop "get_master_BLACS error: blacs_status=F. call init_BLACS"
    master = p_rank==0
  end function get_master_BLACS

  function get_rank_BLACS() result(id)
    integer :: id
    if(.not.blacs_status)stop "get_master_BLACS error: blacs_status=F. call init_BLACS"
    id = p_rank
  end function get_rank_BLACS


  function get_size_BLACS() result(numproc)
    integer :: numproc
    if(.not.blacs_status)stop "get_master_BLACS error: blacs_status=F. call init_BLACS"
    numproc = p_size
  end function get_size_BLACS

#else


  public :: init_BLACS
  public :: finalize_BLACS
  public :: get_master_BLACS
  public :: get_rank_BLACS
  public :: get_size_BLACS
  !


contains

  !****************************************
  !              BLACS START/STOP
  !****************************************
  subroutine init_BLACS()
    return
  end subroutine init_BLACS

  subroutine finalize_BLACS(blacs_end)
    integer,optional :: blacs_end
    return
  end subroutine finalize_BLACS



  !****************************************
  !            BLACS BASIC FUNCTIONS
  !****************************************
  function get_master_BLACS() result(master)
    logical :: master
    master = .true.
  end function get_master_BLACS

  function get_rank_BLACS() result(id)
    integer :: id
    id = 0
  end function get_rank_BLACS

  function get_size_BLACS() result(numproc)
    integer :: numproc
    numproc = 1
  end function get_size_BLACS

#endif

END MODULE SF_BLACS
