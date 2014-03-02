!########################################################################
!Program  : PARSECMD
!PURPOSE  : Declare all the common variables usually in use within codes
!########################################################################
module PARSE_INPUT
  implicit none
  private

  !cmd line variables:
  !=========================================================
  type,public :: input_variable
     character(len=64)                          :: name
     character(len=256)                         :: value
     character(len=20),dimension(:),allocatable :: args
  end type input_variable


  public :: parse_cmd_variable
  interface parse_cmd_variable
     module procedure &
          i_parse_variable, iv_parse_variable, im_parse_variable, &
          d_parse_variable, dv_parse_variable, dm_parse_variable,&
          l_parse_variable, lv_parse_variable, lm_parse_variable,&
          ch_parse_variable,chv_parse_variable,chm_parse_variable
  end interface parse_cmd_variable


  public :: parse_input_variable
  interface parse_input_variable
     module procedure &
          i_parse_input, iv_parse_input, im_parse_input, &
          d_parse_input, dv_parse_input, dm_parse_input,&
          l_parse_input, lv_parse_input, lm_parse_input,&
          ch_parse_input,chv_parse_input,chm_parse_input
  end interface parse_input_variable

  interface append_to_used_input
     module procedure &
          i_append_to_used_file,iv_append_to_used_file,im_append_to_used_file,&
          d_append_to_used_file,dv_append_to_used_file,dm_append_to_used_file,&
          ch_append_to_used_file,chv_append_to_used_file,chm_append_to_used_file,&
          l_append_to_used_file,lv_append_to_used_file,lm_append_to_used_file
  end interface append_to_used_input

  interface txtfy
     module procedure i_to_ch,r_to_ch,c_to_ch,l_to_ch
  end interface txtfy


  public :: get_cmd_variable
  !
  public :: parse_cmd_help

  integer            :: unit
  character(len=7)   :: file_status='replace'
  character(len=255) :: p_buffer

contains




  !---------------------------------------------------------------------
  !PURPOSE:
  !---------------------------------------------------------------------
  function get_cmd_variable(i)  result(var)
    integer            :: i,ncount,pos
    type(input_variable) :: var
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

  function get_input_variable(buffer)  result(var)
    character(len=*)   :: buffer
    integer            :: i,pos,iloc,len
    type(input_variable) :: var
    call s_blank_delete(buffer)
    pos      = scan(buffer,"=")
    var%name = buffer(1:pos-1)
    call upper_case(var%name)
    len=len_trim(buffer)
    if(buffer(len:len)==',')then
       var%value= buffer(pos+1:len-1)
    else
       var%value= buffer(pos+1:)
    endif
  end function get_input_variable





















  !---------------------------------------------------------------------
  !PURPOSE:
  !---------------------------------------------------------------------
  subroutine i_parse_variable(variable,name,default)
    integer                   :: variable
    integer,optional          :: default
    character(len=*)          :: name
    character(len=len(name)) :: name_
    type(input_variable)        :: var
    integer                   :: i
    If(present(default))variable=default
    name_=name;call upper_case(name_)
    do i=1,command_argument_count()
       var = get_cmd_variable(i)
       if(var%name==name_)then
          read(var%value,*)variable
          write(0,*)"Variable "//trim(var%name)//" updated to "//trim(var%value)
       endif
    enddo
  end subroutine i_parse_variable

  subroutine iv_parse_variable(variable,name,default)
    integer,dimension(:)                       :: variable
    integer,dimension(size(variable)),optional :: default
    character(len=*)                           :: name
    character(len=len(name))                  :: name_
    type(input_variable)                         :: var
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
          write(0,"(A,100I6)")" Variable "//trim(var%name)//" updated to ",(variable(iarg),iarg=1,ndim)
       endif
    enddo
  end subroutine iv_parse_variable

  subroutine im_parse_variable(variable,name,default)
    integer,dimension(:,:)                     :: variable
    integer,dimension(size(variable))          :: dummy_var
    integer,dimension(size(variable)),optional :: default
    character(len=*)                           :: name
    character(len=len(name))                   :: name_
    type(input_variable)                         :: var
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
          write(0,"(A,100I6)")" Variable "//trim(var%name)//" updated to ",(dummy_var(iarg),iarg=1,ndim)
       endif
    enddo
  end subroutine im_parse_variable

  subroutine i_parse_input(variable,name,file,default)
    integer                  :: variable
    integer,optional         :: default
    character(len=*)         :: name
    character(len=*)         :: file
    character(len=len(name)) :: name_
    type(input_variable)     :: var
    integer                  :: i
    integer                  :: status
    logical                  :: bool,idefault=.true.
    character(len=255)       :: buffer
    if(present(default))variable=default
    name_=name;call upper_case(name_)
    inquire(file=file,exist=bool)
    if(.not.bool)stop "PARSE_INPUT: input file does not exist."
    unit=free_unit()
    open(unit,file=file)
    status=0
    var_search: do while(status>=0)
       read(unit,"(A)",iostat=status)buffer
       var = get_input_variable(trim(buffer))
       if(var%name==name_)then
          read(var%value,*)variable
          exit var_search
       endif
    enddo var_search
    call parse_cmd_variable(variable,name_)
    close(unit)
    call append_to_used_input(variable,name_,file)
    return
  end subroutine i_parse_input

  subroutine iv_parse_input(variable,name,file,default)
    integer,dimension(:)                       :: variable
    integer,dimension(size(variable)),optional :: default
    character(len=*)                           :: name
    character(len=*)                           :: file
    character(len=len(name))                   :: name_
    type(input_variable)                         :: var
    integer                                    :: i,j,ndim,ncount,nargs,pos0,iarg
    logical                                    :: iscalar
    integer                                    :: status
    logical                                    :: bool
    character(len=255)                         :: buffer
    If(present(default))variable=default
    ndim=size(variable)
    name_=name;call upper_case(name_)
    inquire(file=file,exist=bool)
    if(.not.bool)stop "PARSE_INPUT: input file does not exist."
    unit=free_unit()
    open(unit,file=file)
    status=0
    var_search: do while(status>=0)
       read(unit,"(A)",iostat=status)buffer
       var = get_input_variable(trim(buffer))
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
          exit var_search
       endif
    enddo var_search
    call parse_cmd_variable(variable,name_)
    close(unit)
    call append_to_used_input(variable,name_,file)
    return
  end subroutine iv_parse_input

  subroutine im_parse_input(variable,name,file,default)
    integer,dimension(:,:)                     :: variable
    integer,dimension(size(variable))          :: dummy_var
    integer,dimension(size(variable)),optional :: default
    character(len=*)                           :: name
    character(len=*)                           :: file
    character(len=len(name))                   :: name_
    type(input_variable)                         :: var
    integer                                    :: i,j,ndim,nargs,pos0,iarg
    integer                                    :: status
    logical                                    :: bool
    character(len=255)                         :: buffer
    If(present(default))variable=transpose(reshape(default,shape(variable)))
    ndim=size(variable)
    name_=name;call upper_case(name_)
    inquire(file=file,exist=bool)
    if(.not.bool)stop "PARSE_INPUT: input file does not exist."
    unit=free_unit()
    open(unit,file=file)
    status=0
    var_search: do while(status>=0)
       read(unit,"(A)",iostat=status)buffer
       var = get_input_variable(trim(buffer))
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
          exit var_search
       endif
    enddo var_search
    call parse_cmd_variable(variable,name_)
    close(unit)
    call append_to_used_input(variable,name_,file)
    return
  end subroutine im_parse_input














  !---------------------------------------------------------------------
  !PURPOSE:
  !---------------------------------------------------------------------
  subroutine d_parse_variable(variable,name,default)
    real(8)                   :: variable
    real(8),optional          :: default
    character(len=*)          :: name
    character(len=len(name)) :: name_
    type(input_variable)        :: var
    integer                   :: i
    if(present(default))variable=default
    name_=name;call upper_case(name_)
    do i=1,command_argument_count()
       var = get_cmd_variable(i)
       if(var%name==name_)then
          read(var%value,*)variable
          write(0,*)"Variable "//trim(var%name)//" updated to "//trim(var%value)
       endif
    enddo
  end subroutine d_parse_variable

  subroutine dv_parse_variable(variable,name,default)
    real(8),dimension(:)                       :: variable
    real(8),dimension(size(variable)),optional :: default
    character(len=*)                           :: name
    character(len=len(name))                  :: name_
    type(input_variable)                         :: var
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
          write(0,"(A,100F18.9)")" Variable "//trim(var%name)//" updated to ",(variable(iarg),iarg=1,ndim)
       endif
    enddo
  end subroutine dv_parse_variable

  subroutine dm_parse_variable(variable,name,default)
    real(8),dimension(:,:)                     :: variable
    real(8),dimension(size(variable))          :: dummy_var
    real(8),dimension(size(variable)),optional :: default
    character(len=*)                           :: name
    character(len=len(name))                   :: name_
    type(input_variable)                         :: var
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
          write(0,"(A,100I6)")" Variable "//trim(var%name)//" updated to ",(dummy_var(iarg),iarg=1,ndim)
       endif
    enddo
  end subroutine dm_parse_variable


  subroutine d_parse_input(variable,name,file,default)
    real(8)                  :: variable
    real(8),optional         :: default
    character(len=*)         :: name
    character(len=*)         :: file
    character(len=len(name)) :: name_
    type(input_variable)       :: var
    integer                  :: i
    integer                  :: status
    logical                  :: bool
    character(len=255)       :: buffer
    if(present(default))variable=default
    name_=name;call upper_case(name_)
    inquire(file=file,exist=bool)
    if(.not.bool)stop "PARSE_INPUT: input file does not exist."
    unit=free_unit()
    open(unit,file=file)
    status=0
    var_search: do while(status>=0)
       read(unit,"(A)",iostat=status)buffer
       var = get_input_variable(trim(buffer))
       if(var%name==name_)then
          read(var%value,*)variable
          exit var_search
       endif
    enddo var_search
    call parse_cmd_variable(variable,name_)
    close(unit)
    call append_to_used_input(variable,name_,file)
    return
  end subroutine d_parse_input

  subroutine dv_parse_input(variable,name,file,default)
    real(8),dimension(:)                       :: variable
    real(8),dimension(size(variable)),optional :: default
    character(len=*)                           :: name
    character(len=*)                           :: file
    character(len=len(name))                   :: name_
    type(input_variable)                        :: var
    integer                                    :: i,j,ndim,ncount,nargs,pos0,iarg
    logical                                    :: iscalar
    integer                                    :: status
    logical                                    :: bool
    character(len=255)                         :: buffer
    If(present(default))variable=default
    ndim=size(variable)
    name_=name;call upper_case(name_)
    inquire(file=file,exist=bool)
    if(.not.bool)stop "PARSE_INPUT: input file does not exist."
    unit=free_unit()
    open(unit,file=file)
    status=0
    var_search: do while(status>=0)
       read(unit,"(A)",iostat=status)buffer
       var = get_input_variable(trim(buffer))
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
          exit var_search
       endif
    enddo var_search
    call parse_cmd_variable(variable,name_)
    close(unit)
    call append_to_used_input(variable,name_,file)
    return
  end subroutine dv_parse_input

  subroutine dm_parse_input(variable,name,file,default)
    real(8),dimension(:,:)                     :: variable
    real(8),dimension(size(variable))          :: dummy_var
    real(8),dimension(size(variable)),optional :: default
    character(len=*)                           :: name
    character(len=*)                           :: file
    character(len=len(name))                   :: name_
    type(input_variable)                         :: var
    integer                                    :: i,j,ndim,nargs,pos0,iarg
    integer                                    :: status
    logical                                    :: bool
    character(len=255)                         :: buffer
    If(present(default))variable=transpose(reshape(default,shape(variable)))
    ndim=size(variable)
    name_=name;call upper_case(name_)
    inquire(file=file,exist=bool)
    if(.not.bool)stop "PARSE_INPUT: input file does not exist."
    unit=free_unit()
    open(unit,file=file)
    status=0
    var_search: do while(status>=0)
       read(unit,"(A)",iostat=status)buffer
       var = get_input_variable(trim(buffer))
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
          exit var_search
       endif
    enddo var_search
    call parse_cmd_variable(variable,name_)
    close(unit)
    call append_to_used_input(variable,name_,file)
    return
  end subroutine dm_parse_input














  !---------------------------------------------------------------------
  !PURPOSE:
  !---------------------------------------------------------------------
  subroutine ch_parse_variable(variable,name,default)
    character(len=*)          :: variable
    character(len=*),optional :: default
    character(len=*)          :: name
    character(len=len(name)) :: name_
    type(input_variable)        :: var
    integer                   :: i
    if(present(default))variable=default
    name_=name;call upper_case(name_)
    do i=1,command_argument_count()
       var = get_cmd_variable(i)
       if(var%name==name_)then
          read(var%value,*)variable
          write(0,*)"Variable "//trim(var%name)//" updated to "//trim(var%value)
       endif
    enddo
  end subroutine ch_parse_variable

  subroutine chv_parse_variable(variable,name,default)
    character(len=*),dimension(:)                       :: variable
    character(len=*),dimension(size(variable)),optional :: default
    character(len=*)                           :: name
    character(len=len(name))                  :: name_
    type(input_variable)                         :: var
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
          write(0,"(A,100A20)")"Variable "//trim(var%name)//" updated to ",(variable(iarg),iarg=1,ndim)
       endif
    enddo
  end subroutine chv_parse_variable

  subroutine chm_parse_variable(variable,name,default)
    character(len=*),dimension(:,:)                     :: variable
    character(len=20),dimension(size(variable))         :: dummy_var
    character(len=*),dimension(size(variable)),optional :: default
    character(len=*)                                    :: name
    character(len=len(name))                            :: name_
    type(input_variable)                                  :: var
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
          write(0,"(A,100A20)")"Variable "//trim(var%name)//" updated to ",(dummy_var(iarg),iarg=1,ndim)
       endif
    enddo
  end subroutine chm_parse_variable


  subroutine ch_parse_input(variable,name,file,default)
    character(len=*)                  :: variable
    character(len=*),optional         :: default
    character(len=*)         :: name
    character(len=*)         :: file
    character(len=len(name)) :: name_
    type(input_variable)       :: var
    integer                  :: i
    integer                  :: status
    logical                  :: bool
    character(len=255)       :: buffer
    if(present(default))variable=default
    name_=name;call upper_case(name_)
    inquire(file=file,exist=bool)
    if(.not.bool)stop "PARSE_INPUT: input file does not exist."
    unit=free_unit()
    open(unit,file=file)
    status=0
    var_search: do while(status>=0)
       read(unit,"(A)",iostat=status)buffer
       var = get_input_variable(trim(buffer))
       if(var%name==name_)then
          read(var%value,*)variable
          exit var_search
       endif
    enddo var_search
    call parse_cmd_variable(variable,name_)
    close(unit)
    call append_to_used_input(variable,name_,file)
    return
  end subroutine ch_parse_input

  subroutine chv_parse_input(variable,name,file,default)
    character(len=*),dimension(:)                       :: variable
    character(len=*),dimension(size(variable)),optional :: default
    character(len=*)                           :: name
    character(len=*)                           :: file
    character(len=len(name))                   :: name_
    type(input_variable)                         :: var
    integer                                    :: i,j,ndim,ncount,nargs,pos0,iarg
    logical                                    :: iscalar
    integer                                    :: status
    logical                                    :: bool
    character(len=255)                         :: buffer
    If(present(default))variable=default
    ndim=size(variable)
    name_=name;call upper_case(name_)
    inquire(file=file,exist=bool)
    if(.not.bool)stop "PARSE_INPUT: input file does not exist."
    unit=free_unit()
    open(unit,file=file)
    status=0
    var_search: do while(status>=0)
       read(unit,"(A)",iostat=status)buffer
       var = get_input_variable(trim(buffer))
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
          exit var_search
       endif
    enddo var_search
    call parse_cmd_variable(variable,name_)
    close(unit)
    call append_to_used_input(variable,name_,file)
    return
  end subroutine chv_parse_input

  subroutine chm_parse_input(variable,name,file,default)
    character(len=*),dimension(:,:)                     :: variable
    character(len=20),dimension(size(variable))          :: dummy_var
    character(len=*),dimension(size(variable)),optional :: default
    character(len=*)                           :: name
    character(len=*)                           :: file
    character(len=len(name))                   :: name_
    type(input_variable)                         :: var
    integer                                    :: i,j,ndim,nargs,pos0,iarg
    integer                                    :: status
    logical                                    :: bool
    character(len=255)                         :: buffer
    If(present(default))variable=transpose(reshape(default,shape(variable)))
    ndim=size(variable)
    name_=name;call upper_case(name_)
    inquire(file=file,exist=bool)
    if(.not.bool)stop "PARSE_INPUT: input file does not exist."
    unit=free_unit()
    open(unit,file=file)
    status=0
    var_search: do while(status>=0)
       read(unit,"(A)",iostat=status)buffer
       var = get_input_variable(trim(buffer))
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
          exit var_search
       endif
    enddo var_search
    call parse_cmd_variable(variable,name_)
    close(unit)
    call append_to_used_input(variable,name_,file)
    return
  end subroutine chm_parse_input












  !---------------------------------------------------------------------
  !PURPOSE:
  !---------------------------------------------------------------------
  subroutine l_parse_variable(variable,name,default)
    logical                   :: variable
    logical,optional          :: default
    character(len=*)          :: name
    character(len=len(name)) :: name_
    type(input_variable)        :: var
    integer                   :: i
    if(present(default))variable=default
    name_=name;call upper_case(name_)
    do i=1,command_argument_count()
       var = get_cmd_variable(i)
       if(var%name==name_)then
          read(var%value,*)variable
          write(0,*)"Variable "//trim(var%name)//" updated to "//trim(var%value)
       endif
    enddo
  end subroutine l_parse_variable

  subroutine lv_parse_variable(variable,name,default)
    logical,dimension(:)                       :: variable
    logical,dimension(size(variable)),optional :: default
    character(len=*)                           :: name
    character(len=len(name))                  :: name_
    type(input_variable)                         :: var
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
          write(0,"(A,100L3)")"Variable "//trim(var%name)//" updated to ",(variable(iarg),iarg=1,ndim)
       endif
    enddo
  end subroutine lv_parse_variable

  subroutine lm_parse_variable(variable,name,default)
    logical,dimension(:,:)                     :: variable
    logical,dimension(size(variable))          :: dummy_var
    logical,dimension(size(variable)),optional :: default
    character(len=*)                           :: name
    character(len=len(name))                   :: name_
    type(input_variable)                         :: var
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
          write(0,"(A,100L3)")"Variable "//trim(var%name)//" updated to ",(dummy_var(iarg),iarg=1,ndim)
       endif
    enddo
  end subroutine lm_parse_variable



  subroutine l_parse_input(variable,name,file,default)
    logical               :: variable
    logical,optional         :: default
    character(len=*)         :: name
    character(len=*)         :: file
    character(len=len(name)) :: name_
    type(input_variable)       :: var
    integer                  :: i
    integer                  :: status
    logical                  :: bool
    character(len=255)       :: buffer
    if(present(default))variable=default
    name_=name;call upper_case(name_)
    inquire(file=file,exist=bool)
    if(.not.bool)stop "PARSE_INPUT: input file does not exist."
    unit=free_unit()
    open(unit,file=file)
    status=0
    var_search: do while(status>=0)
       read(unit,"(A)",iostat=status)buffer
       var = get_input_variable(trim(buffer))
       if(var%name==name_)then
          read(var%value,*)variable
          exit var_search
       endif
    enddo var_search
    call parse_cmd_variable(variable,name_)
    close(unit)
    call append_to_used_input(variable,name_,file)
    return
  end subroutine l_parse_input

  subroutine lv_parse_input(variable,name,file,default)
    logical,dimension(:)                       :: variable
    logical,dimension(size(variable)),optional :: default
    character(len=*)                           :: name
    character(len=*)                           :: file
    character(len=len(name))                   :: name_
    type(input_variable)                         :: var
    integer                                    :: i,j,ndim,ncount,nargs,pos0,iarg
    logical                                    :: iscalar
    integer                                    :: status
    logical                                    :: bool
    character(len=255)                         :: buffer
    If(present(default))variable=default
    ndim=size(variable)
    name_=name;call upper_case(name_)
    inquire(file=file,exist=bool)
    if(.not.bool)stop "PARSE_INPUT: input file does not exist."
    unit=free_unit()
    open(unit,file=file)
    status=0
    var_search: do while(status>=0)
       read(unit,"(A)",iostat=status)buffer
       var = get_input_variable(trim(buffer))
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
          exit var_search
       endif
    enddo var_search
    call parse_cmd_variable(variable,name_)
    close(unit)
    call append_to_used_input(variable,name_,file)
    return
  end subroutine lv_parse_input

  subroutine lm_parse_input(variable,name,file,default)
    logical,dimension(:,:)                     :: variable
    logical,dimension(size(variable))          :: dummy_var
    logical,dimension(size(variable)),optional :: default
    character(len=*)                           :: name
    character(len=*)                           :: file
    character(len=len(name))                   :: name_
    type(input_variable)                         :: var
    integer                                    :: i,j,ndim,nargs,pos0,iarg
    integer                                    :: status
    logical                                    :: bool
    character(len=255)                         :: buffer
    If(present(default))variable=transpose(reshape(default,shape(variable)))
    ndim=size(variable)
    name_=name;call upper_case(name_)
    inquire(file=file,exist=bool)
    if(.not.bool)stop "PARSE_INPUT: input file does not exist."
    unit=free_unit()
    open(unit,file=file)
    status=0
    var_search: do while(status>=0)
       read(unit,"(A)",iostat=status)buffer
       var = get_input_variable(trim(buffer))
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
          exit var_search
       endif
    enddo var_search
    call parse_cmd_variable(variable,name_)
    close(unit)
    call append_to_used_input(variable,name_,file)
    return
  end subroutine lm_parse_input











  !---------------------------------------------------------------------
  !PURPOSE:
  !---------------------------------------------------------------------
  subroutine i_append_to_used_file(variable,name,file)
    integer          :: variable
    character(len=*) :: name
    character(len=*) :: file
    integer          :: unit
    call s_blank_delete(name)
    p_buffer=trim(name)//"="//txtfy(variable)//","
    call s_blank_delete(p_buffer)
    unit=free_unit()
    open(unit,file="used."//file,position='append',status=trim(file_status));file_status='old'
    write(unit,"(1x,A)")trim(p_buffer)
    write(*,"(1x,A)")trim(p_buffer)
    close(unit)
  end subroutine i_append_to_used_file

  subroutine iv_append_to_used_file(variable,name,file)
    integer,dimension(:) :: variable
    character(len=*)     :: name
    character(len=*)     :: file
    integer              :: unit,i
    call s_blank_delete(name)
    p_buffer=trim(name)//"="
    do i=1,size(variable)
       p_buffer=trim(p_buffer)//txtfy(variable(i))//","
    end do
    call s_blank_delete(p_buffer)
    unit=free_unit()
    open(unit,file="used."//file,position='append',status=trim(file_status));file_status='old'
    write(unit,"(1x,A)")trim(p_buffer)
    write(*,"(1x,A)")trim(p_buffer)
    close(unit)
  end subroutine iv_append_to_used_file

  subroutine im_append_to_used_file(variable,name,file)
    integer,dimension(:,:) :: variable
    character(len=*)     :: name
    character(len=*)     :: file
    integer              :: unit,i,j
    call s_blank_delete(name)
    p_buffer=trim(name)//"="
    do i=1,size(variable,1)
       do j=1,size(variable,2)
          p_buffer=trim(p_buffer)//txtfy(variable(i,j))//","
       end do
    end do
    call s_blank_delete(p_buffer)
    unit=free_unit()
    open(unit,file="used."//file,position='append',status=trim(file_status));file_status='old'
    write(unit,"(1x,A)")trim(p_buffer)
    write(*,"(1x,A)")trim(p_buffer)
    close(unit)
  end subroutine im_append_to_used_file


  subroutine d_append_to_used_file(variable,name,file)
    real(8)          :: variable
    character(len=*) :: name
    character(len=*) :: file
    integer          :: unit
    call s_blank_delete(name)
    p_buffer=trim(name)//"="//txtfy(variable)
    call s_blank_delete(p_buffer)
    unit=free_unit()
    open(unit,file="used."//file,position='append',status=trim(file_status));file_status='old'
    write(unit,"(1x,A)")trim(p_buffer)
    write(*,"(1x,A)")trim(p_buffer)
    close(unit)
  end subroutine d_append_to_used_file

  subroutine dv_append_to_used_file(variable,name,file)
    real(8),dimension(:) :: variable
    character(len=*)     :: name
    character(len=*)     :: file
    integer              :: unit,i
    call s_blank_delete(name)
    p_buffer=trim(name)//"="
    do i=1,size(variable)
       p_buffer=trim(p_buffer)//txtfy(variable(i))//","
    end do
    call s_blank_delete(p_buffer)
    unit=free_unit()
    open(unit,file="used."//file,position='append',status=trim(file_status));file_status='old'
    write(unit,"(1x,A)")trim(p_buffer)
    write(*,"(1x,A)")trim(p_buffer)
    close(unit)
  end subroutine dv_append_to_used_file

  subroutine dm_append_to_used_file(variable,name,file)
    real(8),dimension(:,:) :: variable
    character(len=*)     :: name
    character(len=*)     :: file
    integer              :: unit,i,j
    call s_blank_delete(name)
    p_buffer=trim(name)//"="
    do i=1,size(variable,1)
       do j=1,size(variable,2)
          p_buffer=trim(p_buffer)//txtfy(variable(i,j))//","
       end do
    end do
    call s_blank_delete(p_buffer)
    unit=free_unit()
    open(unit,file="used."//file,position='append',status=trim(file_status));file_status='old'
    write(unit,"(1x,A)")trim(p_buffer)
    write(*,"(1x,A)")trim(p_buffer)
    close(unit)
  end subroutine dm_append_to_used_file



  subroutine ch_append_to_used_file(variable,name,file)
    character(len=*) :: variable
    character(len=*) :: name
    character(len=*) :: file
    integer          :: unit
    call s_blank_delete(name)
    p_buffer=trim(name)//"="//trim(variable)
    call s_blank_delete(p_buffer)
    unit=free_unit()
    open(unit,file="used."//file,position='append',status=trim(file_status));file_status='old'
    write(unit,"(1x,A)")trim(p_buffer)
    write(*,"(1x,A)")trim(p_buffer)
    close(unit)
  end subroutine ch_append_to_used_file

  subroutine chv_append_to_used_file(variable,name,file)
    character(len=*),dimension(:) :: variable
    character(len=*)     :: name
    character(len=*)     :: file
    integer              :: unit,i
    call s_blank_delete(name)
    p_buffer=trim(name)//"="
    do i=1,size(variable)
       p_buffer=trim(p_buffer)//trim(variable(i))//","
    end do
    call s_blank_delete(p_buffer)
    unit=free_unit()
    open(unit,file="used."//file,position='append',status=trim(file_status));file_status='old'
    write(unit,"(1x,A)")trim(p_buffer)
    write(*,"(1x,A)")trim(p_buffer)
    close(unit)
  end subroutine chv_append_to_used_file

  subroutine chm_append_to_used_file(variable,name,file)
    character(len=*),dimension(:,:) :: variable
    character(len=*)     :: name
    character(len=*)     :: file
    integer              :: unit,i,j
    call s_blank_delete(name)
    p_buffer=trim(name)//"="
    do i=1,size(variable,1)
       do j=1,size(variable,2)
          p_buffer=trim(p_buffer)//trim(variable(i,j))//","
       end do
    end do
    call s_blank_delete(p_buffer)
    unit=free_unit()
    open(unit,file="used."//file,position='append',status=trim(file_status));file_status='old'
    write(unit,"(1x,A)")trim(p_buffer)
    write(*,"(1x,A)")trim(p_buffer)
    close(unit)
  end subroutine chm_append_to_used_file



  subroutine l_append_to_used_file(variable,name,file)
    logical          :: variable
    character(len=*) :: name
    character(len=*) :: file
    integer          :: unit
    call s_blank_delete(name)
    p_buffer=trim(name)//"="//txtfy(variable)
    call s_blank_delete(p_buffer)
    unit=free_unit()
    open(unit,file="used."//file,position='append',status=trim(file_status));file_status='old'
    write(unit,"(1x,A)")trim(p_buffer)
    write(*,"(1x,A)")trim(p_buffer)
    close(unit)
  end subroutine l_append_to_used_file

  subroutine lv_append_to_used_file(variable,name,file)
    logical,dimension(:) :: variable
    character(len=*)     :: name
    character(len=*)     :: file
    integer              :: unit,i
    call s_blank_delete(name)
    p_buffer=trim(name)//"="
    do i=1,size(variable)
       p_buffer=trim(p_buffer)//txtfy(variable(i))//","
    end do
    call s_blank_delete(p_buffer)
    unit=free_unit()
    open(unit,file="used."//file,position='append',status=trim(file_status));file_status='old'
    write(unit,"(1x,A)")trim(p_buffer)
    write(*,"(1x,A)")trim(p_buffer)
    close(unit)
  end subroutine lv_append_to_used_file

  subroutine lm_append_to_used_file(variable,name,file)
    logical,dimension(:,:) :: variable
    character(len=*)     :: name
    character(len=*)     :: file
    integer              :: unit,i,j
    call s_blank_delete(name)
    p_buffer=trim(name)//"="
    do i=1,size(variable,1)
       do j=1,size(variable,2)
          p_buffer=trim(p_buffer)//txtfy(variable(i,j))//","
       end do
    end do
    call s_blank_delete(p_buffer)
    unit=free_unit()
    open(unit,file="used."//file,position='append',status=trim(file_status));file_status='old'
    write(unit,"(1x,A)")trim(p_buffer)
    write(*,"(1x,A)")trim(p_buffer)
    close(unit)
  end subroutine lm_append_to_used_file





  !---------------------------------------------------------------------
  !PURPOSE:
  !---------------------------------------------------------------------
  function check_cmd_vector_size(ndim,var) result(nargs)
    integer            :: ndim,ncount,j,nargs
    type(input_variable) :: var
    logical            :: iscalar
    if(ndim==1)then
       nargs=ndim
       return
    endif
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
    type(input_variable) :: var
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







  !---------------------------------------------------------------------
  !PURPOSE:
  !---------------------------------------------------------------------
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


  subroutine s_blank_delete ( s )
    !! S_BLANK_DELETE removes blanks from a string, left justifying the remainder.
    !    All TAB characters are also removed.
    !    Input/output, character ( len = * ) S, the string to be transformed.
    implicit none
    character              ch
    integer   ( kind = 4 ) get
    integer   ( kind = 4 ) put
    character ( len = * )  s
    integer   ( kind = 4 ) s_length
    character, parameter :: tab = achar ( 9 )
    put = 0
    s_length = len_trim ( s )
    do get = 1, s_length
       ch = s(get:get)
       if ( ch /= ' ' .and. ch /= tab ) then
          put = put + 1
          s(put:put) = ch
       end if
    end do
    s(put+1:s_length) = ' '
    return
  end subroutine s_blank_delete


  function i_to_ch(i4) result(string)
    character(len=32) :: string
    integer           :: i4
    call i4_to_s_left(i4,string)
  end function i_to_ch

  function r_to_ch(r8) result(string)
    character(len=32) :: string
    character(len=16) :: string_
    real(8)           :: r8
    call r8_to_s_left(r8,string_)
    string=adjustl(string_)
  end function r_to_ch

  function c_to_ch(c) result(string)
    character(len=32+3) :: string
    character(len=16) :: sre,sim
    complex(8)        :: c
    real(8)           :: re,im
    re=real(c,8);im=aimag(c)
    call r8_to_s_left(re,sre)
    call r8_to_s_left(im,sim)
    string="("//trim(sre)//","//trim(sim)//")"
  end function c_to_ch

  function l_to_ch(bool) result(string)
    logical :: bool
    character(len=1) :: string
    string='F'
    if(bool)string='T'
  end function l_to_ch

  subroutine i4_to_s_left ( i4, s )
    character :: c
    integer   :: i
    integer   :: i4
    integer   :: idig
    integer   :: ihi
    integer   :: ilo
    integer   :: ipos
    integer   :: ival
    character(len=*) ::  s
    s = ' '
    ilo = 1
    ihi = len ( s )
    if ( ihi <= 0 ) then
       return
    end if
    !  Make a copy of the integer.
    ival = i4
    !  Handle the negative sign.
    if ( ival < 0 ) then
       if ( ihi <= 1 ) then
          s(1:1) = '*'
          return
       end if
       ival = -ival
       s(1:1) = '-'
       ilo = 2
    end if
    !  The absolute value of the integer goes into S(ILO:IHI).
    ipos = ihi
    !  Find the last digit of IVAL, strip it off, and stick it into the string.
    do
       idig = mod ( ival, 10 )
       ival = ival / 10
       if ( ipos < ilo ) then
          do i = 1, ihi
             s(i:i) = '*'
          end do
          return
       end if
       call digit_to_ch ( idig, c )
       s(ipos:ipos) = c
       ipos = ipos - 1
       if ( ival == 0 ) then
          exit
       end if
    end do
    !  Shift the string to the left.
    s(ilo:ilo+ihi-ipos-1) = s(ipos+1:ihi)
    s(ilo+ihi-ipos:ihi) = ' '
  end subroutine i4_to_s_left

  subroutine r8_to_s_left ( r8, s )
    integer :: i
    real(8) :: r8
    integer :: s_length
    character(len=*) ::  s
    character(len=16) :: s2
    s_length = len ( s )
    if ( s_length < 16 ) then
       do i = 1, s_length
          s(i:i) = '*'
       end do
    else if ( r8 == 0.0D+00 ) then
       s(1:16) = '     0.0      '
    else
       write ( s2, '(g16.9)' ) r8
       s(1:16) = s2
    end if
    !  Shift the string left.
    s = adjustl ( s )
  end subroutine r8_to_s_left


  subroutine digit_to_ch(digit,ch)
    character :: ch
    integer   :: digit
    if ( 0 <= digit .and. digit <= 9 ) then
       ch = achar ( digit + 48 )
    else
       ch = '*'
    end if
  end subroutine digit_to_ch


  function free_unit() result(unit_)
    integer :: unit_,ios
    logical :: is_it_opened
    unit_=100
    do 
       unit_=unit_+1
       INQUIRE(unit=unit_,OPENED=is_it_opened,iostat=ios)
       if(.not.is_it_opened.AND.ios==0)return 
       if(unit_>900) stop "ERROR free_unit: no unit free smaller than 900. Possible BUG"
    enddo
  end function free_unit


end module PARSE_INPUT



! program test_parse_vars
!   USE PARSE_INPUT
!   integer              :: L,M(2),A(2,2)
!   real(8)        :: x,xx(2)
!   L=0
!   M=0
!   A=0
!   x=0.d0
!   xx=0.d0
!   call parse_input_variable(L,"l","input.in")
!   call parse_input_variable(M,"m","input.in",default=[2,2])
!   call parse_input_variable(A,"A","input.in",default=[2,2,2,2])
!   call parse_input_variable(x,"x","input.in",default=1.d0)
!   call parse_input_variable(xx,"xx","input.in",default=[1.d0,0.d0])
! end program test_parse_vars
