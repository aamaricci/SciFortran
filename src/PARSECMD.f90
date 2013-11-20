!########################################################################
!Program  : PARSECMD
!PURPOSE  : Declare all the common variables usually in use within codes
!########################################################################
module PARSE_CMD
  USE COMMON_VARS, only: msg
  implicit none
  private 

  !cmd line variables:
  !=========================================================
  type,public:: cmd_variable
     character(len=64)                      :: name
     character(len=256)                     :: value
     character(len=20),dimension(:),allocatable :: args
  end type cmd_variable
  type(cmd_variable),public      :: cmd_var,nml_var

  interface parse_cmd_variable
     module procedure &
          i_parse_variable, iv_parse_variable, im_parse_variable, &
          d_parse_variable, dv_parse_variable, dm_parse_variable,&
          l_parse_variable, lv_parse_variable, lm_parse_variable,&
          ch_parse_variable,chv_parse_variable,chm_parse_variable
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

  function check_cmd_vector_size(ndim,var) result(nargs)
    integer            :: ndim,ncount,j,nargs
    type(cmd_variable) :: var
    logical            :: iscalar
    iscalar=(scan(var%value,",")==0)
    if(iscalar)then
       print*,"error in parse_cmd array:",trim(var%name)
       print*,"expecting a comma separated list of: ",ndim
       stop
    endif
    ncount=0
    do j=1,len(var%value)
       if(var%value(j:j)==",")ncount=ncount+1
    enddo
    nargs=ncount+1
    if(nargs/=ndim)then
       print*,"parse_variable wrong dimensions: ",trim(var%name)
       print*,"expecting a comma separated list of: ",ndim
       stop
    endif
  end function check_cmd_vector_size


  function check_cmd_matrix_size(ndim,var) result(nargs)
    integer            :: ndim,ncount,j,nargs
    type(cmd_variable) :: var
    logical            :: iscalar
    iscalar=(scan(var%value,",")==0)
    if(iscalar)then
       print*,"error in parse_cmd array:",trim(var%name)
       print*,"expecting a comma separated list of: ",ndim
       stop
    endif
    ncount=0
    do j=1,len(var%value)
       if(var%value(j:j)==",")ncount=ncount+1
    enddo
    nargs=ncount+1
    if(nargs/=ndim)then
       print*,"parse_variable wrong dimensions: ",trim(var%name)
       print*,"expecting a comma separated list of: ",ndim
       stop
    endif
  end function check_cmd_matrix_size







  !******************************************************************
  !******************************************************************
  !******************************************************************



  subroutine i_parse_variable(variable,name,default)
    integer                   :: variable
    integer,optional          :: default
    character(len=*)          :: name
    character(len=len(name)) :: name_
    type(cmd_variable)        :: var
    integer                   :: i
    If(present(default))variable=default
    name_=name;call upper_case(name_)
    do i=1,command_argument_count()
       var = get_cmd_variable(i)
       if(var%name==name_)then
          read(var%value,*)variable
          write(*,*)"Variable "//trim(var%name)//" update to "//trim(var%value)
       endif
    enddo
  end subroutine i_parse_variable

  subroutine iv_parse_variable(variable,name,default)
    integer,dimension(:)                       :: variable
    integer,dimension(size(variable)),optional :: default
    character(len=*)                           :: name
    character(len=len(name))                  :: name_
    type(cmd_variable)                         :: var
    integer                                    :: i,j,ndim,ncount,nargs,pos0,iarg
    logical                                    :: iscalar
    If(present(default))variable=default
    ndim=size(variable)
    name_=name;call upper_case(name_)
    do i=1,command_argument_count()
       var = get_cmd_variable(i)
       if(var%name==name_)then
          iscalar=(scan(var%value,",")==0)
          nargs=check_cmd_vector_size(ndim,var)
          allocate(var%args(nargs))
          iarg=0
          pos0=0
          do j=1,len(var%value)
             if(var%value(j:j)==",")then
                iarg=iarg+1
                var%args(iarg)=var%value(pos0+1:j-1)
                pos0=j
             endif
          enddo
          var%args(nargs)=var%value(pos0+1:)
          do iarg=1,nargs
             read(var%args(iarg),*)variable(iarg)
          enddo
          write(*,"(A,100I6)")"Variable "//trim(var%name)//" update to ",(variable(iarg),iarg=1,ndim)
       endif
    enddo
  end subroutine iv_parse_variable

  subroutine im_parse_variable(variable,name,default)
    integer,dimension(:,:)                     :: variable
    integer,dimension(size(variable))          :: dummy_var
    integer,dimension(size(variable)),optional :: default
    character(len=*)                           :: name
    character(len=len(name))                   :: name_
    type(cmd_variable)                         :: var
    integer                                    :: i,j,ndim,nargs,pos0,iarg
    If(present(default))variable=transpose(reshape(default,shape(variable)))
    ndim=size(variable)
    name_=name;call upper_case(name_)
    do i=1,command_argument_count()
       var = get_cmd_variable(i)
       if(var%name==name_)then
          nargs=check_cmd_vector_size(ndim,var)
          allocate(var%args(nargs))
          iarg=0
          pos0=0
          do j=1,len(var%value)
             if(var%value(j:j)==",")then
                iarg=iarg+1
                var%args(iarg)=var%value(pos0+1:j-1)
                pos0=j
             endif
          enddo
          var%args(nargs)=var%value(pos0+1:)
          do iarg=1,nargs
             read(var%args(iarg),*)dummy_var(iarg)
          enddo
          variable=transpose(reshape(dummy_var,shape(variable)))
          write(*,"(A,100I6)")"Variable "//trim(var%name)//" update to ",(dummy_var(iarg),iarg=1,ndim)
       endif
    enddo
  end subroutine im_parse_variable



  !******************************************************************
  !******************************************************************
  !******************************************************************



  subroutine d_parse_variable(variable,name,default)
    real(8)                   :: variable
    real(8),optional          :: default
    character(len=*)          :: name
    character(len=len(name)) :: name_
    type(cmd_variable)        :: var
    integer                   :: i
    if(present(default))variable=default
    name_=name;call upper_case(name_)
    do i=1,command_argument_count()
       var = get_cmd_variable(i)
       if(var%name==name_)then
          read(var%value,*)variable
          write(*,*)"Variable "//trim(var%name)//" update to "//trim(var%value)
       endif
    enddo
  end subroutine d_parse_variable

  subroutine dv_parse_variable(variable,name,default)
    real(8),dimension(:)                       :: variable
    real(8),dimension(size(variable)),optional :: default
    character(len=*)                           :: name
    character(len=len(name))                  :: name_
    type(cmd_variable)                         :: var
    integer                                    :: i,j,ndim,ncount,nargs,pos0,iarg
    logical                                    :: iscalar
    If(present(default))variable=default
    ndim=size(variable)
    name_=name;call upper_case(name_)
    do i=1,command_argument_count()
       var = get_cmd_variable(i)
       if(var%name==name_)then
          nargs=check_cmd_vector_size(ndim,var)
          allocate(var%args(nargs))
          iarg=0
          pos0=0
          do j=1,len(var%value)
             if(var%value(j:j)==",")then
                iarg=iarg+1
                var%args(iarg)=var%value(pos0+1:j-1)
                pos0=j
             endif
          enddo
          var%args(nargs)=var%value(pos0+1:)
          do iarg=1,nargs
             read(var%args(iarg),*)variable(iarg)
          enddo
          write(*,"(A,100F18.9)")"Variable "//trim(var%name)//" update to ",(variable(iarg),iarg=1,ndim)
       endif
    enddo
  end subroutine dv_parse_variable

  subroutine dm_parse_variable(variable,name,default)
    real(8),dimension(:,:)                     :: variable
    real(8),dimension(size(variable))          :: dummy_var
    real(8),dimension(size(variable)),optional :: default
    character(len=*)                           :: name
    character(len=len(name))                   :: name_
    type(cmd_variable)                         :: var
    integer                                    :: i,j,ndim,nargs,pos0,iarg
    If(present(default))variable=transpose(reshape(default,shape(variable)))
    ndim=size(variable)
    name_=name;call upper_case(name_)
    do i=1,command_argument_count()
       var = get_cmd_variable(i)
       if(var%name==name_)then
          nargs=check_cmd_vector_size(ndim,var)
          allocate(var%args(nargs))
          iarg=0
          pos0=0
          do j=1,len(var%value)
             if(var%value(j:j)==",")then
                iarg=iarg+1
                var%args(iarg)=var%value(pos0+1:j-1)
                pos0=j
             endif
          enddo
          var%args(nargs)=var%value(pos0+1:)
          do iarg=1,nargs
             read(var%args(iarg),*)dummy_var(iarg)
          enddo
          variable=transpose(reshape(dummy_var,shape(variable)))
          write(*,"(A,100I6)")"Variable "//trim(var%name)//" update to ",(dummy_var(iarg),iarg=1,ndim)
       endif
    enddo
  end subroutine dm_parse_variable


  !******************************************************************
  !******************************************************************
  !******************************************************************


  subroutine ch_parse_variable(variable,name,default)
    character(len=*)          :: variable
    character(len=*),optional :: default
    character(len=*)          :: name
    character(len=len(name)) :: name_
    type(cmd_variable)        :: var
    integer                   :: i
    if(present(default))variable=default
    name_=name;call upper_case(name_)
    do i=1,command_argument_count()
       var = get_cmd_variable(i)
       if(var%name==name_)then
          read(var%value,*)variable
          write(*,*)"Variable "//trim(var%name)//" update to "//trim(var%value)
       endif
    enddo
  end subroutine ch_parse_variable

  subroutine chv_parse_variable(variable,name,default)
    character(len=*),dimension(:)                       :: variable
    character(len=*),dimension(size(variable)),optional :: default
    character(len=*)                           :: name
    character(len=len(name))                  :: name_
    type(cmd_variable)                         :: var
    integer                                    :: i,j,ndim,ncount,nargs,pos0,iarg
    logical                                    :: iscalar
    If(present(default))variable=default
    ndim=size(variable)
    name_=name;call upper_case(name_)
    do i=1,command_argument_count()
       var = get_cmd_variable(i)
       if(var%name==name_)then
          nargs=check_cmd_vector_size(ndim,var)
          allocate(var%args(nargs))
          iarg=0
          pos0=0
          do j=1,len(var%value)
             if(var%value(j:j)==",")then
                iarg=iarg+1
                var%args(iarg)=var%value(pos0+1:j-1)
                pos0=j
             endif
          enddo
          var%args(nargs)=var%value(pos0+1:)
          do iarg=1,nargs
             read(var%args(iarg),*)variable(iarg)
          enddo
          write(*,"(A,100A20)")"Variable "//trim(var%name)//" update to ",(variable(iarg),iarg=1,ndim)
       endif
    enddo
  end subroutine chv_parse_variable

  subroutine chm_parse_variable(variable,name,default)
    character(len=*),dimension(:,:)                     :: variable
    character(len=20),dimension(size(variable))         :: dummy_var
    character(len=*),dimension(size(variable)),optional :: default
    character(len=*)                                    :: name
    character(len=len(name))                            :: name_
    type(cmd_variable)                                  :: var
    integer                                             :: i,j,ndim,nargs,pos0,iarg
    If(present(default))variable=transpose(reshape(default,shape(variable)))
    ndim=size(variable)
    name_=name;call upper_case(name_)
    do i=1,command_argument_count()
       var = get_cmd_variable(i)
       if(var%name==name_)then
          nargs=check_cmd_vector_size(ndim,var)
          allocate(var%args(nargs))
          iarg=0
          pos0=0
          do j=1,len(var%value)
             if(var%value(j:j)==",")then
                iarg=iarg+1
                var%args(iarg)=var%value(pos0+1:j-1)
                pos0=j
             endif
          enddo
          var%args(nargs)=var%value(pos0+1:)
          do iarg=1,nargs
             read(var%args(iarg),*)dummy_var(iarg)
          enddo
          variable=transpose(reshape(dummy_var,shape(variable)))
          write(*,"(A,100I6)")"Variable "//trim(var%name)//" update to ",(dummy_var(iarg),iarg=1,ndim)
       endif
    enddo
  end subroutine chm_parse_variable



  !******************************************************************
  !******************************************************************
  !******************************************************************


  subroutine l_parse_variable(variable,name,default)
    logical                   :: variable
    logical,optional          :: default
    character(len=*)          :: name
    character(len=len(name)) :: name_
    type(cmd_variable)        :: var
    integer                   :: i
    if(present(default))variable=default
    name_=name;call upper_case(name_)
    do i=1,command_argument_count()
       var = get_cmd_variable(i)
       if(var%name==name_)then
          read(var%value,*)variable
          write(*,*)"Variable "//trim(var%name)//" update to "//trim(var%value)
       endif
    enddo
  end subroutine l_parse_variable

  subroutine lv_parse_variable(variable,name,default)
    logical,dimension(:)                       :: variable
    logical,dimension(size(variable)),optional :: default
    character(len=*)                           :: name
    character(len=len(name))                  :: name_
    type(cmd_variable)                         :: var
    integer                                    :: i,j,ndim,ncount,nargs,pos0,iarg
    logical                                    :: iscalar
    If(present(default))variable=default
    ndim=size(variable)
    name_=name;call upper_case(name_)
    do i=1,command_argument_count()
       var = get_cmd_variable(i)
       if(var%name==name_)then
          nargs=check_cmd_vector_size(ndim,var)
          allocate(var%args(nargs))
          iarg=0
          pos0=0
          do j=1,len(var%value)
             if(var%value(j:j)==",")then
                iarg=iarg+1
                var%args(iarg)=var%value(pos0+1:j-1)
                pos0=j
             endif
          enddo
          var%args(nargs)=var%value(pos0+1:)
          do iarg=1,nargs
             read(var%args(iarg),*)variable(iarg)
          enddo
          write(*,"(A,100L3)")"Variable "//trim(var%name)//" update to ",(variable(iarg),iarg=1,ndim)
       endif
    enddo
  end subroutine lv_parse_variable

  subroutine lm_parse_variable(variable,name,default)
    logical,dimension(:,:)                     :: variable
    logical,dimension(size(variable))          :: dummy_var
    logical,dimension(size(variable)),optional :: default
    character(len=*)                           :: name
    character(len=len(name))                   :: name_
    type(cmd_variable)                         :: var
    integer                                    :: i,j,ndim,nargs,pos0,iarg
    If(present(default))variable=transpose(reshape(default,shape(variable)))
    ndim=size(variable)
    name_=name;call upper_case(name_)
    do i=1,command_argument_count()
       var = get_cmd_variable(i)
       if(var%name==name_)then
          nargs=check_cmd_vector_size(ndim,var)
          allocate(var%args(nargs))
          iarg=0
          pos0=0
          do j=1,len(var%value)
             if(var%value(j:j)==",")then
                iarg=iarg+1
                var%args(iarg)=var%value(pos0+1:j-1)
                pos0=j
             endif
          enddo
          var%args(nargs)=var%value(pos0+1:)
          do iarg=1,nargs
             read(var%args(iarg),*)dummy_var(iarg)
          enddo
          variable=transpose(reshape(dummy_var,shape(variable)))
          write(*,"(A,100I6)")"Variable "//trim(var%name)//" update to ",(dummy_var(iarg),iarg=1,ndim)
       endif
    enddo
  end subroutine lm_parse_variable




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

  subroutine lower_case(s)
    integer   ( kind = 4 ) i
    character ( len = * )  s
    integer   ( kind = 4 ) s_length
    s_length = len_trim ( s )
    do i = 1, s_length
       call ch_low ( s(i:i) )
    end do
  end subroutine lower_case

  subroutine ch_cap(ch)
    character              ch
    integer   ( kind = 4 ) itemp
    itemp = iachar ( ch )
    if ( 97 <= itemp .and. itemp <= 122 ) then
       ch = achar ( itemp - 32 )
    end if
  end subroutine ch_cap

  subroutine ch_low ( ch )
    character ch
    integer ( kind = 4 ) i
    i = iachar ( ch )
    if ( 65 <= i .and. i <= 90 ) then
       ch = achar ( i + 32 )
    end if
  end subroutine ch_low


end module PARSE_CMD
