!########################################################################
!PROGRAM  : COMVARS
!TYPE     : Module
!PURPOSE  : Declare all the common variables usually in use within codes
!########################################################################
module COMMON_VARS
  USE MPI
  USE OMP_LIB
  implicit none
  private

  include "scifor_version.inc"

  !PARAMETERS
  !===============================================================
  COMPLEX(8),PARAMETER,public :: ZERO=(0.D0,0.D0)
  COMPLEX(8),PARAMETER,public :: XI=(0.D0,1.D0)
  COMPLEX(8),PARAMETER,public :: ONE=(1.D0,0.D0)
  REAL(8),PARAMETER,public    :: SQRT2 = 1.41421356237309504880169D0
  REAL(8),PARAMETER,public    :: SQRT3 = 1.73205080756887729352745D0
  REAL(8),PARAMETER,public    :: SQRT6 = 2.44948974278317809819728D0
  REAL(8),PARAMETER,public    :: PI    = 3.14159265358979323846264338327950288419716939937510D0
  REAL(8),PARAMETER,public    :: PI2   = 6.28318530717959D0
  REAL(8),PARAMETER,public    :: GAMMA_EULER = 0.57721566490153286060D0  !EULER S CONSTANT
  REAL(8),PARAMETER,public    :: EULER= 2.7182818284590452353602874713526624977572470936999595749669676277240766303535D0
  INTEGER,PARAMETER,public    :: MAX_INT  = HUGE(1) 
  REAL(8),PARAMETER,public    :: MAX_REAL = HUGE(1.D0)
  REAL(8),PARAMETER,public    :: MAX_EXP  =  700.D0
  REAL(8),PARAMETER,public    :: MIN_EXP  = -700.D0
  REAL(8),PARAMETER,public    :: MAX_EXP_QUAD = 11000.D0
  REAL(8),PARAMETER,public    :: MIN_EXP_QUAD =-11000.D0
  REAL(8),PARAMETER,public    :: MAX_EXP_R =  81.D0
  REAL(8),PARAMETER,public    :: MIN_EXP_R = -81.D0
  !REAL(8),PARAMETER,public    :: ERROR=EPSILON(1.D0)*100.D0
  REAL(8),PARAMETER,public    :: EPSILONR=EPSILON(1.D0),EPSILONQ=1.D-30
  INTEGER,PARAMETER,public    :: DBL=8,DP=8        ! "DOUBLE" PRECISION
  INTEGER,PARAMETER,public    :: DDP=16            ! "QUAD"   PRECISION
  INTEGER,PARAMETER,public    :: SP = KIND(1.0)    ! "SINGLE" PRECISION
  REAL(DBL),PARAMETER,public  :: TINY_   = 1.D-12
  REAL(DBL),PARAMETER,public  :: HUGE_   = 1.D+12
  REAL(DBL),PARAMETER,public  :: QUARTER = 0.25_DBL
  REAL(DBL),PARAMETER,public  :: THIRD   = 0.3333333333333_DBL
  REAL(DBL),PARAMETER,public  :: HALF    = 0.5_DBL
  LOGICAL,PARAMETER,public    :: TT=.TRUE., FF=.FALSE.


  !GLOABL  VARIABLES
  !=========================================================
  INTEGER,public  :: ILOOP,NLOOP    !DMFT LOOP VARIABLES
  REAL(8),public  :: D              !BANDWIDTH
  REAL(8),public  :: TS,TSP,TPP     !N.N./N.N.N. HOPPING AMPLITUDE
  REAL(8),public  :: U,V            !LOCAL,NON-LOCAL INTERACTION
  REAL(8),public  :: TPD,VPD        !HYBRIDIZATION,BAND-BAND COUPLING
  REAL(8),public  :: ED0,EP0        !ORBITAL ENERGIES
  REAL(8),public  :: XMU            !CHEMICAL POTENTIAL
  REAL(8),public  :: DT,DTAU        !TIME STEP
  REAL(8),public  :: FMESH          !FREQ. STEP
  REAL(8),public  :: BETA           !INVERSE TEMPERATURE
  REAL(8),public  :: TEMP           !TEMPERATURE
  REAL(8),public  :: EPS            !BROADENING


  !CMD LINE VARIABLES:
  !=========================================================
  TYPE,public:: CMD_VARIABLE
     CHARACTER(LEN=64)           :: NAME
     CHARACTER(LEN=64)           :: VALUE
  END TYPE CMD_VARIABLE
  CHARACTER(LEN=512),ALLOCATABLE,public :: HELP_BUFFER(:)
  TYPE(CMD_VARIABLE),public             :: CMD_VAR,NML_VAR

  INTERFACE PARSE_CMD_VARIABLE
     MODULE PROCEDURE D_PARSE_VARIABLE,&
          I_PARSE_VARIABLE,CH_PARSE_VARIABLE,L_PARSE_VARIABLE
  END INTERFACE PARSE_CMD_VARIABLE


  !MPI VARS:
  !=========================================================
  INTEGER,public          :: MPIERR,MPISIZE,MPIID
  CHARACTER(LEN=3),public :: MPICHAR


  !OMP VARS:
  !=========================================================
  INTEGER,public :: OMP_NUM_THREADS,OMP_ID,OMP_SIZE

  !Date variables:

  integer(4)                          :: year
  integer(4)                          :: mese
  integer(4)                          :: day
  integer(4)                          :: h
  integer(4)                          :: m
  integer(4)                          :: s
  integer(4)                          :: ms
  character(len=9),parameter,dimension(12) :: month = (/ &
       'January  ', 'February ', 'March    ', 'April    ', &
       'May      ', 'June     ', 'July     ', 'August   ', &
       'September', 'October  ', 'November ', 'December ' /)


  !===============================================
  public :: parse_cmd_variable,parse_cmd_help,get_cmd_variable
  public :: version
  public :: timestamp
  public :: abort
  public :: msg
  public :: start_mpi,close_mpi
  public :: bold
  public :: underline
  public :: highlight
  public :: erased
  public :: red,green,yellow,blue,purple,cyan
  public :: bold_red,bold_green,bold_yellow,bold_blue,bold_purple,bold_cyan
  public :: bg_red,bg_green,bg_yellow,bg_blue,bg_purple,bg_cyan


contains


  !+-------------------------------------------------------------------+
  !PURPOSE  : print actual version of the software (if any)
  !+-------------------------------------------------------------------+
  subroutine version(revision)
    character(len=*) :: revision
    if(mpiID==0)then
       write(*,"(A)")bg_green("SCIFOR VERSION (GIT): "//trim(adjustl(trim(sf_version))))
       write(*,"(A)")bg_green("CODE VERSION (GIT): "//trim(adjustl(trim(revision))))
       call timestamp

       open(10,file="version.inc")
       write(10,"(A)")"SCIFOR VERSION (GIT): "//trim(adjustl(trim(sf_version)))
       write(10,"(A)")"CODE VERSION (GIT): "//trim(adjustl(trim(revision)))
       call timestamp(unit=10)
       close(10)
    endif
  end subroutine version


  !******************************************************************
  !******************************************************************
  !******************************************************************


  !+-------------------------------------------------------------------+
  !PURPOSE  : prints the current YMDHMS date as a time stamp.
  ! Example: 31 May 2001   9:45:54.872 AM
  !+-------------------------------------------------------------------+
  subroutine timestamp(unit)
    integer,optional        :: unit
    integer                 :: unit_
    integer(4),dimension(8) :: data
    unit_=6;if(present(unit))unit_=unit
    if(mpiID==0)then
       call date_and_time(values=data)
       call print_date(data,unit_)
    endif
  end subroutine timestamp

  !******************************************************************
  !******************************************************************
  !******************************************************************


  !+-------------------------------------------------------------------+
  !PURPOSE  : print actual date
  !+-------------------------------------------------------------------+
  subroutine print_date(dummy,unit)
    integer(4),dimension(8) :: dummy
    integer                 :: unit
    integer(4)                          :: year
    integer(4)                          :: mese
    integer(4)                          :: day
    integer(4)                          :: h
    integer(4)                          :: m
    integer(4)                          :: s
    integer(4)                          :: ms
    character(len=9),parameter,dimension(12) :: month = (/ &
         'January  ', 'February ', 'March    ', 'April    ', &
         'May      ', 'June     ', 'July     ', 'August   ', &
         'September', 'October  ', 'November ', 'December ' /)
    year = dummy(1)
    mese = dummy(2)
    day  = dummy(3)
    h    = dummy(5)
    m    = dummy(6)
    s    = dummy(7)
    ms   = dummy(8)
    write(unit,"(A,i2,1x,a,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3)")&
         "Timestamp: +",day,trim(month(mese)),year, h,':',m,':',s,'.',ms
    write(unit,*)""
  end subroutine print_date



  !******************************************************************
  !******************************************************************
  !******************************************************************



  !+-------------------------------------------------------------------+
  !PURPOSE  : send abort message to std.out and exit 
  !+-------------------------------------------------------------------+
  subroutine abort(text,id,stop)
    character(len=*) :: text
    character(len=4) :: char_id
    integer,optional :: id
    integer          :: i,id_
    logical,optional :: stop
    logical          :: stop_
    id_=0;if(present(id))id_=id
    stop_=.true.;if(present(stop))stop_=stop
    if(id_ > mpiSIZE)id_=0
    if(mpiID==id_)then
       if(mpiID==0)then
          write(*,'(A)',advance="no")bold_red("error:")
       else
          write(char_id,"(I4)")id_
          write(*,'(A)',advance="no")bold_red("error from cpu"//char_id//":")
       endif
       write(*,'(A)')bg_red(text)
    endif
    if(stop_)stop
  end subroutine abort


  !******************************************************************
  !******************************************************************
  !******************************************************************


  subroutine msg(message,lines,id)
    character(len=*) :: message
    integer,optional :: lines,id
    integer          :: i,id_
    id_=0;if(present(id))id_=id
    if(id_ > mpiSIZE)id_=0
    if(mpiID==id_)then
       if(mpiID==0)then
          write(*,'(A)',advance="no")"msg:"
       else
          write(*,'(A,I3,A)',advance="no")"msg from cpu",id_,": "
       endif
       write(*,'(A)') message
       if(present(lines))then
          do i=1,lines
             write(*,'(A)')""
          enddo
       endif
    endif
  end subroutine msg


  !******************************************************************
  !******************************************************************
  !******************************************************************


  subroutine start_mpi
    call MPI_INIT(mpiERR)
    call MPI_COMM_RANK(MPI_COMM_WORLD,mpiID,mpiERR)
    call MPI_COMM_SIZE(MPI_COMM_WORLD,mpiSIZE,mpiERR)
    write(*,"(A,I4,A,I4,A)")'Processor ',mpiID,' of ',mpiSIZE,' is alive'
    call MPI_BARRIER(MPI_COMM_WORLD,mpiERR)
  end subroutine start_mpi


  !******************************************************************
  !******************************************************************
  !******************************************************************



  subroutine close_mpi
    call MPI_BARRIER(MPI_COMM_WORLD,mpiERR)
    call MPI_FINALIZE(mpiERR)
  end subroutine close_mpi



  !******************************************************************
  !******************************************************************
  !******************************************************************


  subroutine parse_cmd_help(buffer,status)
    character(len=256)            :: cmd_line
    character(len=*),dimension(:) :: buffer
    integer                       :: i,lines,N
    logical,optional              :: status
    if(present(status))status=.false.
    N=size(buffer)
    do i=1,command_argument_count()
       call get_command_argument(i,cmd_line)
       if(cmd_line=="--help" .OR. &
            cmd_line=="-h"   .OR. &
            cmd_line=="info" .OR. &
            cmd_line=="--h"  .OR. &
            cmd_line=="help")then
          do lines=1,N
             write(*,"(256A)")trim(adjustl(trim(buffer(lines))))
          enddo
          if(present(status))then
             status=.true.
             return
          else
             stop
          endif
       endif
    enddo
  end subroutine parse_cmd_help


  function get_cmd_variable(i)  result(var)
    integer            :: i,ncount,pos
    type(cmd_variable) :: var
    character(len=512) :: buffer
    ncount=command_argument_count()
    if(i>ncount)then
       write(*,*)"get_cmd_variable: i > cmdline_argument_count!"
       return
    endif
    call get_command_argument(i,buffer)
    pos      = scan(buffer,"=")
    var%name = buffer(1:pos-1);call s_cap(var%name)
    var%value= buffer(pos+1:)
  end function get_cmd_variable


  subroutine i_parse_variable(variable,name,name2,default)
    integer                   :: variable
    integer,optional          :: default
    character(len=*)          :: name
    character(len=*),optional :: name2
    character(len=512)        :: buffer
    type(cmd_variable)        :: var
    integer                   :: i
    if(present(default))variable=default
    do i=1,command_argument_count()
       var = get_cmd_variable(i)
       if(var%name==name)read(var%value,*)variable
       if(present(name2))then
          if(var%name==name2)read(var%value,*)variable
       endif
    enddo
  end subroutine i_parse_variable

  subroutine d_parse_variable(variable,name,name2,default)
    real(8)                   :: variable
    real(8),optional          :: default
    character(len=*)          :: name
    character(len=*),optional :: name2
    character(len=512)        :: buffer
    type(cmd_variable)        :: var
    integer                   :: i
    if(present(default))variable=default
    do i=1,command_argument_count()
       var = get_cmd_variable(i)
       if(var%name==name)read(var%value,*)variable
       if(present(name2))then
          if(var%name==name2)read(var%value,*)variable
       endif
    enddo
  end subroutine d_parse_variable

  subroutine ch_parse_variable(variable,name,name2,default)
    character(len=*)          :: variable
    character(len=*),optional :: default
    character(len=*)          :: name
    character(len=*),optional :: name2
    character(len=512)        :: buffer
    type(cmd_variable)        :: var
    integer                   :: i
    if(present(default))variable=default
    do i=1,command_argument_count()
       var = get_cmd_variable(i)
       if(var%name==name)read(var%value,*)variable
       if(present(name2))then
          if(var%name==name2)read(var%value,*)variable
       endif
    enddo
  end subroutine ch_parse_variable

  subroutine l_parse_variable(variable,name,name2,default)
    logical                   :: variable
    logical,optional          :: default
    character(len=*)          :: name
    character(len=*),optional :: name2
    character(len=512)        :: buffer
    type(cmd_variable)        :: var
    integer                   :: i
    if(present(default))variable=default
    do i=1,command_argument_count()
       var = get_cmd_variable(i)
       if(var%name==name)read(var%value,*)variable
       if(present(name2))then
          if(var%name==name2)read(var%value,*)variable
       endif
    enddo
  end subroutine l_parse_variable


  !******************************************************************
  !******************************************************************
  !******************************************************************

  function bold(text) result(textout)
    character(len=*) :: text
    character(len=8+len(text)) :: textout
    textout=achar(27)//"[1m"//text//achar(27)//"[0m"
  end function bold

  function underline(text) result(textout)
    character(len=*) :: text
    character(len=8+len(text)) :: textout
    textout=achar(27)//"[4m"//text//achar(27)//"[0m"
  end function underline

  function highlight(text) result(textout)
    character(len=*) :: text
    character(len=8+len(text)) :: textout
    textout=achar(27)//"[7m"//text//achar(27)//"[0m"
  end function highlight

  function erased(text) result(textout)
    character(len=*) :: text
    character(len=8+len(text)) :: textout
    textout=achar(27)//"[9m"//text//achar(27)//"[0m"
  end function erased

  function red(text) result(textout)
    character(len=*) :: text
    character(len=9+len(text)) :: textout
    textout=achar(27)//"[91m"//text//achar(27)//"[0m"
  end function red

  function green(text) result(textout)
    character(len=*) :: text
    character(len=9+len(text)) :: textout
    textout=achar(27)//"[92m"//text//achar(27)//"[0m"
  end function green

  function yellow(text) result(textout)
    character(len=*) :: text
    character(len=9+len(text)) :: textout
    textout=achar(27)//"[93m"//text//achar(27)//"[0m"
  end function yellow

  function blue(text) result(textout)
    character(len=*) :: text
    character(len=9+len(text)) :: textout
    textout=achar(27)//"[94m"//text//achar(27)//"[0m"
  end function blue

  function purple(text) result(textout)
    character(len=*) :: text
    character(len=9+len(text)) :: textout
    textout=achar(27)//"[95m"//text//achar(27)//"[0m"
  end function purple

  function cyan(text) result(textout)
    character(len=*) :: text
    character(len=9+len(text)) :: textout
    textout=achar(27)//"[96m"//text//achar(27)//"[0m"
  end function cyan


  function bold_red(text) result(textout)
    character(len=*) :: text
    character(len=11+len(text)) :: textout
    textout=achar(27)//"[1;91m"//text//achar(27)//"[0m"
  end function bold_red

  function bold_green(text) result(textout)
    character(len=*) :: text
    character(len=11+len(text)) :: textout
    textout=achar(27)//"[1;92m"//text//achar(27)//"[0m"
  end function bold_green

  function bold_yellow(text) result(textout)
    character(len=*) :: text
    character(len=11+len(text)) :: textout
    textout=achar(27)//"[1;93m"//text//achar(27)//"[0m"
  end function bold_yellow

  function bold_blue(text) result(textout)
    character(len=*) :: text
    character(len=11+len(text)) :: textout
    textout=achar(27)//"[1;94m"//text//achar(27)//"[0m"
  end function bold_blue

  function bold_purple(text) result(textout)
    character(len=*) :: text
    character(len=11+len(text)) :: textout
    textout=achar(27)//"[1;95m"//text//achar(27)//"[0m"
  end function bold_purple

  function bold_cyan(text) result(textout)
    character(len=*) :: text
    character(len=11+len(text)) :: textout
    textout=achar(27)//"[1;96m"//text//achar(27)//"[0m"
  end function bold_cyan



  function bg_red(text) result(textout)
    character(len=*) :: text
    character(len=9+len(text)) :: textout
    textout=achar(27)//"[41m"//text//achar(27)//"[0m"
  end function bg_red

  function bg_green(text) result(textout)
    character(len=*) :: text
    character(len=9+len(text)) :: textout
    textout=achar(27)//"[42m"//text//achar(27)//"[0m"
  end function bg_green

  function bg_yellow(text) result(textout)
    character(len=*) :: text
    character(len=9+len(text)) :: textout
    textout=achar(27)//"[43m"//text//achar(27)//"[0m"
  end function bg_yellow

  function bg_blue(text) result(textout)
    character(len=*) :: text
    character(len=9+len(text)) :: textout
    textout=achar(27)//"[44m"//text//achar(27)//"[0m"
  end function bg_blue

  function bg_purple(text) result(textout)
    character(len=*) :: text
    character(len=9+len(text)) :: textout
    textout=achar(27)//"[45m"//text//achar(27)//"[0m"
  end function bg_purple

  function bg_cyan(text) result(textout)
    character(len=*) :: text
    character(len=9+len(text)) :: textout
    textout=achar(27)//"[46m"//text//achar(27)//"[0m"
  end function bg_cyan



  !******************************************************************
  !******************************************************************
  !******************************************************************


  subroutine s_cap ( s )
    character              ch
    integer   ( kind = 4 ) i
    character ( len = * )  s
    integer   ( kind = 4 ) s_length
    s_length = len_trim ( s )
    do i = 1, s_length
       ch = s(i:i)
       call ch_cap ( ch )
       s(i:i) = ch
    end do
  end subroutine s_cap

  subroutine ch_cap ( ch )
    character              ch
    integer   ( kind = 4 ) itemp
    itemp = iachar ( ch )
    if ( 97 <= itemp .and. itemp <= 122 ) then
       ch = achar ( itemp - 32 )
    end if
  end subroutine ch_cap

END MODULE COMMON_VARS
