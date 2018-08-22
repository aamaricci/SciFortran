MODULE SF_MPI
#ifdef _MPI
  USE MPI
#endif

  implicit none
  private


#ifdef _MPI
  interface Bcast_MPI
     module procedure :: MPI_Bcast_Bool_0
     module procedure :: MPI_Bcast_Bool_1
     module procedure :: MPI_Bcast_Bool_2
     module procedure :: MPI_Bcast_Bool_3
     module procedure :: MPI_Bcast_Bool_4
     module procedure :: MPI_Bcast_Bool_5
     module procedure :: MPI_Bcast_Bool_6
     module procedure :: MPI_Bcast_Bool_7
     !
     module procedure :: MPI_Bcast_Int_0
     module procedure :: MPI_Bcast_Int_1
     module procedure :: MPI_Bcast_Int_2
     module procedure :: MPI_Bcast_Int_3
     module procedure :: MPI_Bcast_Int_4
     module procedure :: MPI_Bcast_Int_5
     module procedure :: MPI_Bcast_Int_6
     module procedure :: MPI_Bcast_Int_7
     !
     module procedure :: MPI_Bcast_Dble_0
     module procedure :: MPI_Bcast_Dble_1
     module procedure :: MPI_Bcast_Dble_2
     module procedure :: MPI_Bcast_Dble_3
     module procedure :: MPI_Bcast_Dble_4
     module procedure :: MPI_Bcast_Dble_5
     module procedure :: MPI_Bcast_Dble_6
     module procedure :: MPI_Bcast_Dble_7
     !
     module procedure :: MPI_Bcast_Cmplx_0
     module procedure :: MPI_Bcast_Cmplx_1
     module procedure :: MPI_Bcast_Cmplx_2
     module procedure :: MPI_Bcast_Cmplx_3
     module procedure :: MPI_Bcast_Cmplx_4
     module procedure :: MPI_Bcast_Cmplx_5
     module procedure :: MPI_Bcast_Cmplx_6
     module procedure :: MPI_Bcast_Cmplx_7     
  end interface Bcast_MPI



  interface AllGather_MPI
     module procedure :: MPI_AllGather_Bool_0
     module procedure :: MPI_AllGather_Bool_1
     module procedure :: MPI_AllGather_Bool_2
     module procedure :: MPI_AllGather_Bool_3
     module procedure :: MPI_AllGather_Bool_4
     module procedure :: MPI_AllGather_Bool_5
     module procedure :: MPI_AllGather_Bool_6
     module procedure :: MPI_AllGather_Bool_7
     !
     module procedure :: MPI_AllGather_Int_0
     module procedure :: MPI_AllGather_Int_1
     module procedure :: MPI_AllGather_Int_2
     module procedure :: MPI_AllGather_Int_3
     module procedure :: MPI_AllGather_Int_4
     module procedure :: MPI_AllGather_Int_5
     module procedure :: MPI_AllGather_Int_6
     module procedure :: MPI_AllGather_Int_7
     !
     module procedure :: MPI_AllGather_Dble_0
     module procedure :: MPI_AllGather_Dble_1
     module procedure :: MPI_AllGather_Dble_2
     module procedure :: MPI_AllGather_Dble_3
     module procedure :: MPI_AllGather_Dble_4
     module procedure :: MPI_AllGather_Dble_5
     module procedure :: MPI_AllGather_Dble_6
     module procedure :: MPI_AllGather_Dble_7
     !
     module procedure :: MPI_AllGather_Cmplx_0
     module procedure :: MPI_AllGather_Cmplx_1
     module procedure :: MPI_AllGather_Cmplx_2
     module procedure :: MPI_AllGather_Cmplx_3
     module procedure :: MPI_AllGather_Cmplx_4
     module procedure :: MPI_AllGather_Cmplx_5
     module procedure :: MPI_AllGather_Cmplx_6
     module procedure :: MPI_AllGather_Cmplx_7     
  end interface AllGather_MPI




  interface AllReduce_MPI
     module procedure :: MPI_AllReduce_Bool_0
     module procedure :: MPI_AllReduce_Bool_1
     module procedure :: MPI_AllReduce_Bool_2
     module procedure :: MPI_AllReduce_Bool_3
     module procedure :: MPI_AllReduce_Bool_4
     module procedure :: MPI_AllReduce_Bool_5
     module procedure :: MPI_AllReduce_Bool_6
     module procedure :: MPI_AllReduce_Bool_7
     !
     module procedure :: MPI_AllReduce_Int_0
     module procedure :: MPI_AllReduce_Int_1
     module procedure :: MPI_AllReduce_Int_2
     module procedure :: MPI_AllReduce_Int_3
     module procedure :: MPI_AllReduce_Int_4
     module procedure :: MPI_AllReduce_Int_5
     module procedure :: MPI_AllReduce_Int_6
     module procedure :: MPI_AllReduce_Int_7
     !
     module procedure :: MPI_AllReduce_Dble_0
     module procedure :: MPI_AllReduce_Dble_1
     module procedure :: MPI_AllReduce_Dble_2
     module procedure :: MPI_AllReduce_Dble_3
     module procedure :: MPI_AllReduce_Dble_4
     module procedure :: MPI_AllReduce_Dble_5
     module procedure :: MPI_AllReduce_Dble_6
     module procedure :: MPI_AllReduce_Dble_7
     !
     module procedure :: MPI_AllReduce_Cmplx_0
     module procedure :: MPI_AllReduce_Cmplx_1
     module procedure :: MPI_AllReduce_Cmplx_2
     module procedure :: MPI_AllReduce_Cmplx_3
     module procedure :: MPI_AllReduce_Cmplx_4
     module procedure :: MPI_AllReduce_Cmplx_5
     module procedure :: MPI_AllReduce_Cmplx_6
     module procedure :: MPI_AllReduce_Cmplx_7     
  end interface AllReduce_MPI


  public :: Init_MPI
  public :: Finalize_MPI
  public :: StartMsg_MPI
  public :: Barrier_MPI
  !
  public :: Get_Size_MPI
  public :: Get_Rank_MPI
  public :: Get_Master_MPI
  public :: Get_Last_MPI
  public :: Get_Q_MPI
  public :: Get_R_MPI
  public :: Get_Chunk_MPI
  public :: Get_Wtime_MPI
  public :: Get_Wtick_MPI
  public :: Error_MPI
  !
  public :: Bcast_MPI
  public :: AllGather_MPI
  public :: AllReduce_MPI


  integer :: ierr
  integer :: rank

contains


  !****************************************
  !              MPI START/STOP
  !****************************************
  subroutine Init_MPI(comm,msg)
    integer,optional :: comm
    logical,optional :: msg
    call MPI_Init(ierr)
    call Error_MPI(ierr,"MPI_Start")
    if(present(comm))comm=MPI_COMM_WORLD
    if(present(msg).AND.(msg.eqv..true.))call StartMsg_MPI(MPI_COMM_WORLD)
  end subroutine Init_MPI

  subroutine Finalize_MPI(comm)
    integer,optional :: comm
    call MPI_Finalize(ierr)
    call Error_MPI(ierr,"MPI_Stop")
    if(present(comm))comm=MPI_COMM_NULL
  end subroutine Finalize_MPI

  subroutine StartMsg_MPI(comm)
    integer,optional :: comm
    integer          :: comm_,size
    integer          :: i
    comm_=MPI_COMM_WORLD;if(present(comm))comm_=comm
    rank = Get_Rank_MPI(comm_)
    size = Get_Size_MPI(comm_)
    if(rank==0)write(*,'(a)')"---------------MPI----------------"
    do i=0,size-1
       call MPI_Barrier(comm_,ierr)
       if(rank==i)write(*,"(A,I6,A,I6,A)")"rank:",rank," of ",size," alive"          
    enddo
    call MPI_Barrier(comm_,ierr)
    if(rank==0)write(*,'(a)')"----------------------------------"
    if(rank==0)write(*,'(a)')""
  end subroutine StartMsg_MPI

  subroutine Barrier_MPI(comm)
    integer,optional :: comm
    integer          :: comm_
    comm_=MPI_COMM_WORLD;if(present(comm))comm_=comm
    call MPI_Barrier(comm_,ierr)
    call Error_MPI(ierr,"Barrier_MPI")
  end subroutine Barrier_MPI


  !****************************************
  !              MPI TOOLS
  !****************************************
  function Get_size_MPI(comm) result(size)
    integer :: comm
    integer :: size
    call MPI_Comm_size(comm,size,ierr)
    call Error_MPI(ierr,"Get_Size_MPI")
  end function Get_size_MPI

  function Get_rank_MPI(comm) result(rank)
    integer :: comm
    integer :: rank
    call MPI_Comm_rank(comm,rank,ierr)
    call Error_MPI(ierr,"Get_Rank_MPI")
  end function Get_rank_MPI

  function Get_Processor_MPI() result(workstation)
    integer                               :: istat
    character(len=MPI_MAX_PROCESSOR_NAME) :: workstation
    call MPI_GET_PROCESSOR_NAME(workstation,istat,ierr)
    call Error_MPI(ierr,"Get_Processor_MPI")
  end function Get_Processor_MPI

  function Get_master_MPI(comm) result(master)
    integer :: comm
    integer :: size
    integer :: rank
    logical :: master
    call MPI_Comm_rank(comm,rank,ierr)
    call Error_MPI(ierr,"Get_Master_MPI")
    master=.false.
    if(rank==0)master=.true.
  end function Get_master_MPI

  function Get_last_MPI(comm) result(last)
    integer :: comm
    integer :: size
    integer :: rank
    logical :: last
    call MPI_Comm_rank(comm,rank,ierr)
    call MPI_Comm_size(comm,size,ierr)    
    last=.false.
    if(rank==size-1)last=.true.
  end function Get_last_MPI

  function Get_Q_MPI(comm,N) result(mpiQ)
    integer :: comm
    integer :: N
    integer :: size
    integer :: rank
    integer :: mpiQ
    size = Get_size_MPI(comm)
    mpiQ = N/size
  end function Get_Q_MPI

  function Get_R_MPI(comm,N) result(mpiR)
    integer :: comm
    integer :: N
    integer :: size
    integer :: rank
    integer :: mpiR
    logical :: last
    size = Get_size_MPI(comm)
    last = Get_last_MPI(comm)
    mpiR=0
    if(last)mpiR = mod(N,size)
  end function Get_R_MPI

  function Get_Chunk_MPI(comm,N) result(Nchunk)
    integer :: comm
    integer :: N
    integer :: Nchunk
    Nchunk = Get_Q_MPI(comm,N)+Get_R_MPI(comm,N)
  end function Get_Chunk_MPI

  
  !returns an elapsed time on the calling processor
  function Get_Wtime_MPI() result(time)
    real(8) :: time
    time = MPI_WTIME()
  end function Get_Wtime_MPI


  !returns the resolution of MPI_Wtime
  function Get_Wtick_MPI() result(tick)
    real(8) :: tick
    tick = MPI_WTICK()
  end function Get_Wtick_MPI







  !****************************************
  !              MPI BROADCAST
  !****************************************
  include "mpi_bcast.f90"



  !****************************************
  !              MPI ALLGATHER
  !****************************************
  include "mpi_allgather.f90"



  !****************************************
  !              MPI ALLREDUCE
  !****************************************
  include "mpi_allreduce.f90"




# define STR_ERR_COMM      'invalid communicator in mpi call.'
# define STR_ERR_COUNT     'invalid count in mpi call.'
# define STR_ERR_TYPE      'invalid datatype in mpi call.'
# define STR_ERR_BUFFER    'invalid buffer in mpi call.'
# define STR_ERR_ROOT      'invalid root in mpi call.'
# define STR_ERR_ARG       'invalid argument in mpi call.'
# define STR_ERR_TAG       'invalid tag in mpi call.'
# define STR_ERR_RANK      'invalid rank in mpi call.'
# define STR_ERR_GROUP     'null group passed to mpi call.'
# define STR_ERR_OP        'invalid operation in mpi call.'
# define STR_ERR_TOPOLOGY  'invalid topology in mpi call.'
# define STR_ERR_DIMS      'illegal dimension argument in mpi call.'
# define STR_ERR_UNKNOWN   'unknown error in mpi call.'
# define STR_ERR_TRUNCATE  'message truncated on receive in mpi call.'
# define STR_ERR_OTHER     'other error in mpi call.'
# define STR_ERR_INTERN    'internal error code in mpi call.'
# define STR_ERR_IN_STATUS 'look in status for error value.'
# define STR_ERR_PENDING   'pending request in mpi call.'
# define STR_ERR_REQUEST   'illegal mpi_request handle in mpi call.'
# define STR_ERR_LASTCODE  'last error code in mpi call.'
  subroutine Error_MPI(err,sub)
    integer,optional,intent(in)          :: err
    character(len=*),optional,intent(in) :: sub
    integer                              :: err_
    character(len=128)                   :: sub_
    err_=ierr            ; if(present(err))err_=err
    sub_="MPI_Get_Error:"; if(present(sub))sub_=sub
    select case (err_)
    case (MPI_SUCCESS)
       return
       !
    case (MPI_ERR_COMM)
       write(*,'(2A)')  trim(sub_),STR_ERR_COMM
       !
    case (MPI_ERR_COUNT)
       write(*,'(2A)')  trim(sub_),STR_ERR_COUNT
       !
    case (MPI_ERR_TYPE)
       write(*,'(2A)')  trim(sub_),STR_ERR_TYPE
       !
    case (MPI_ERR_BUFFER)
       write(*,'(2A)')  trim(sub_),STR_ERR_BUFFER
       !
    case (MPI_ERR_ROOT)
       write(*,'(2A)')  trim(sub_),STR_ERR_ROOT
       !
    case (MPI_ERR_ARG)
       write(*,'(2A)')  trim(sub_),STR_ERR_ARG
       !
    case (MPI_ERR_TAG)
       write(*,'(2A)')  trim(sub_),STR_ERR_TAG
       !
    case (MPI_ERR_RANK)
       write(*,'(2A)')  trim(sub_),STR_ERR_RANK
       !
    case (MPI_ERR_GROUP)
       write(*,'(2A)')  trim(sub_),STR_ERR_GROUP
       !
    case (MPI_ERR_OP)
       write(*,'(2A)')  trim(sub_),STR_ERR_OP
       !
    case (MPI_ERR_TOPOLOGY)
       write(*,'(2A)')  trim(sub_),STR_ERR_TOPOLOGY
       !
    case (MPI_ERR_DIMS)
       write(*,'(2A)')  trim(sub_),STR_ERR_DIMS
       !
    case (MPI_ERR_UNKNOWN)
       write(*,'(2A)')  trim(sub_),STR_ERR_UNKNOWN
       !
    case (MPI_ERR_TRUNCATE)
       write(*,'(2A)')  trim(sub_),STR_ERR_TRUNCATE
       !
    case (MPI_ERR_OTHER)
       write(*,'(2A)')  trim(sub_),STR_ERR_OTHER
       !
    case (MPI_ERR_INTERN)
       write(*,'(2A)')  trim(sub_),STR_ERR_INTERN
       !
    case (MPI_ERR_IN_STATUS)
       write(*,'(2A)')  trim(sub_),STR_ERR_IN_STATUS
       !
    case (MPI_ERR_PENDING)
       write(*,'(2A)')  trim(sub_),STR_ERR_PENDING
       !
    case (MPI_ERR_REQUEST)
       write(*,'(2A)')  trim(sub_),STR_ERR_REQUEST
       !
    case (MPI_ERR_LASTCODE)
       write(*,'(2A)')  trim(sub_),STR_ERR_LASTCODE
       !
    case default
       return
       !
    end select
  end subroutine Error_MPI



#else



  public :: Init_MPI
  public :: Finalize_MPI
  public :: StartMsg_MPI
  !
  public :: Get_Size_MPI
  public :: Get_Rank_MPI
  public :: Get_Master_MPI
  public :: Get_Last_MPI
  public :: Get_Q_MPI
  public :: Get_R_MPI
  public :: Get_Chunk_MPI
  !

  integer :: size
  integer :: rank
  integer :: ierr

contains


  !****************************************
  !              MPI START/STOP
  !****************************************
  subroutine Init_MPI()
    return
  end subroutine Init_MPI

  subroutine Finalize_MPI()
    return
  end subroutine Finalize_MPI

  subroutine StartMsg_MPI(comm)
    integer :: comm
    return
  end subroutine StartMsg_MPI



  !****************************************
  !              MPI TOOLS
  !****************************************
  function Get_size_MPI(comm) result(size)
    integer :: comm
    integer :: size
    size=1
  end function Get_size_MPI

  function Get_rank_MPI(comm) result(rank)
    integer :: comm
    integer :: rank
    rank=0
  end function Get_rank_MPI

  function Get_master_MPI(comm) result(master)
    integer :: comm
    logical :: master
    master=.true.
  end function Get_master_MPI

  function Get_last_MPI(comm) result(last)
    integer :: comm
    logical :: last
    last=.true.
  end function Get_last_MPI

  function Get_Q_MPI(comm,N) result(mpiQ)
    integer :: comm
    integer :: N
    integer :: mpiQ
    mpiQ = N
  end function Get_Q_MPI

  function Get_R_MPI(comm,N) result(mpiR)
    integer :: comm
    integer :: N
    integer :: mpiR
    mpiR=0
  end function Get_R_MPI

  function Get_Chunk_MPI(comm,N) result(Nchunk)
    integer :: comm
    integer :: N
    integer :: Nchunk
    Nchunk = N
  end function Get_Chunk_MPI

#endif

END MODULE SF_MPI


