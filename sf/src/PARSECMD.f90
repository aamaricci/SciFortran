!########################################################################
!PROGRAM  : PARSECMD
!PURPOSE  : Declare all the common variables usually in use within codes
!########################################################################
module PARSE_CMD
  USE COMMON_VARS, only: msg
  implicit none
  private 

  !cmd line variables:
  !=========================================================
  type,public:: cmd_variable
     character(len=64)           :: name
     character(len=64)           :: value
  end type cmd_variable
  type(cmd_variable),public      :: cmd_var,nml_var

  interface parse_cmd_variable
     module procedure d_parse_variable,&
          i_parse_variable,ch_parse_variable,l_parse_variable
  end interface parse_cmd_variable

  public :: parse_cmd_variable
  public :: get_cmd_variable
  public :: parse_cmd_help
  public :: print_cmd_help

contains


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


  !******************************************************************
  !******************************************************************
  !******************************************************************



  subroutine print_cmd_help(buffer)
    character(len=*),dimension(:) :: buffer
    integer                       :: lines,N
    N=size(buffer)
    do lines=1,N
       write(*,"(256A)")trim(adjustl(trim(buffer(lines))))
    enddo
    stop
  end subroutine print_cmd_help


  !******************************************************************
  !******************************************************************
  !******************************************************************



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
    var%name = buffer(1:pos-1);call upper_case(var%name)
    var%value= buffer(pos+1:)
  end function get_cmd_variable


  !******************************************************************
  !******************************************************************
  !******************************************************************



  subroutine i_parse_variable(variable,name,name2,default)
    integer                   :: variable
    integer,optional          :: default
    character(len=*)          :: name
    character(len=*),optional :: name2
    type(cmd_variable)        :: var
    integer                   :: i
    If(present(default))variable=default
    do i=1,command_argument_count()
       var = get_cmd_variable(i)
       if(var%name==name)then
          read(var%value,*)variable
          call msg(("Variable "//trim(var%name)//" update to "//trim(var%value)))
       endif
       if(present(name2) .AND. var%name==name2)then
          read(var%value,*)variable
          call msg(("Variable "//trim(var%name)//" update to "//trim(var%value)))
       endif
    enddo
  end subroutine i_parse_variable

  subroutine d_parse_variable(variable,name,name2,default)
    real(8)                   :: variable
    real(8),optional          :: default
    character(len=*)          :: name
    character(len=*),optional :: name2
    type(cmd_variable)        :: var
    integer                   :: i
    if(present(default))variable=default
    do i=1,command_argument_count()
       var = get_cmd_variable(i)
       if(var%name==name)then
          read(var%value,*)variable
          call msg(("Variable "//trim(var%name)//" update to "//trim(var%value)))
       endif
       if(present(name2) .AND. var%name==name2)then
          read(var%value,*)variable
          call msg(("Variable "//trim(var%name)//" update to "//trim(var%value)))
       endif
    enddo
  end subroutine d_parse_variable

  subroutine ch_parse_variable(variable,name,name2,default)
    character(len=*)          :: variable
    character(len=*),optional :: default
    character(len=*)          :: name
    character(len=*),optional :: name2
    type(cmd_variable)        :: var
    integer                   :: i
    if(present(default))variable=default
    do i=1,command_argument_count()
       var = get_cmd_variable(i)
       if(var%name==name)then
          read(var%value,*)variable
          call msg(("Variable "//trim(var%name)//" update to "//trim(var%value)))
       endif
       if(present(name2) .AND. var%name==name2)then
          read(var%value,*)variable
          call msg(("Variable "//trim(var%name)//" update to "//trim(var%value)))
       endif
    enddo
  end subroutine ch_parse_variable

  subroutine l_parse_variable(variable,name,name2,default)
    logical                   :: variable
    logical,optional          :: default
    character(len=*)          :: name
    character(len=*),optional :: name2
    type(cmd_variable)        :: var
    integer                   :: i
    if(present(default))variable=default
    do i=1,command_argument_count()
       var = get_cmd_variable(i)
       if(var%name==name)then
          read(var%value,*)variable
          call msg(("Variable "//trim(var%name)//" update to "//trim(var%value)))
       endif
       if(present(name2) .AND. var%name==name2)then
          read(var%value,*)variable
          call msg(("Variable "//trim(var%name)//" update to "//trim(var%value)))
       endif
    enddo
  end subroutine l_parse_variable



  !******************************************************************
  !******************************************************************
  !******************************************************************



  subroutine upper_case(s)
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
  end subroutine upper_case

  subroutine ch_cap(ch)
    character              ch
    integer   ( kind = 4 ) itemp
    itemp = iachar ( ch )
    if ( 97 <= itemp .and. itemp <= 122 ) then
       ch = achar ( itemp - 32 )
    end if
  end subroutine ch_cap


end module PARSE_CMD
