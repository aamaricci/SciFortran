!=====================0-dimension=====================================
subroutine i_parse_input(variable,name,file,default,comment)
  integer                  :: variable
  integer,optional         :: default
  character(len=*)         :: name
  character(len=*),optional:: comment
  character(len=*)         :: file
  character(len=len(name)) :: name_
  type(input_variable)     :: var
  integer                  :: i,unit,pos
  integer                  :: status
  logical                  :: bool,idefault=.true.
  character(len=255)       :: buffer
  if(present(default))variable=default
  name_=name;call upper_case(name_)
  inquire(file=file,exist=bool)
  if(.not.bool)IOinput=.false.
  if(bool)then
     unit=free_unit()
     open(unit,file=file)
     status=0
     var_search: do while(status>=0)
        read(unit,"(A255)",iostat=status)buffer
        pos=scan(buffer,"!");if(pos/=0)buffer=buffer(1:pos-1)
        var = get_input_variable(trim(buffer))
        if(var%name==name_)then
           read(var%value,*)variable
           exit var_search
        endif
     enddo var_search
     close(unit)
  endif
  call parse_cmd_variable(variable,name_)
  if(present(comment))then
     call append_to_input_list(variable,name_,comment)
  else
     call append_to_input_list(variable,name_)
  endif
  return
end subroutine i_parse_input

subroutine d_parse_input(variable,name,file,default,comment)
  real(8)                  :: variable
  real(8),optional         :: default
  character(len=*)         :: name
  character(len=*),optional:: comment
  character(len=*)         :: file
  character(len=len(name)) :: name_
  type(input_variable)       :: var
  integer                  ::i,unit,pos
  integer                  :: status
  logical                  :: bool
  character(len=255)       :: buffer
  if(present(default))variable=default
  name_=name;call upper_case(name_)
  inquire(file=file,exist=bool)
  if(.not.bool)IOinput=.false.
  if(bool)then
     unit=free_unit()
     open(unit,file=file)
     status=0
     var_search: do while(status>=0)
        read(unit,"(A255)",iostat=status)buffer
        pos=scan(buffer,"!");if(pos/=0)buffer=buffer(1:pos-1)
        var = get_input_variable(trim(buffer))
        if(var%name==name_)then
           read(var%value,*)variable
           exit var_search
        endif
     enddo var_search
     close(unit)
  endif
  call parse_cmd_variable(variable,name_)
  if(present(comment))then
     call append_to_input_list(variable,name_,comment)
  else
     call append_to_input_list(variable,name_)
  endif
  return
end subroutine d_parse_input

subroutine ch_parse_input(variable,name,file,default,comment)
  character(len=*)                  :: variable
  character(len=*),optional         :: default
  character(len=*)         :: name
  character(len=*),optional:: comment
  character(len=*)         :: file
  character(len=len(name)) :: name_
  type(input_variable)       :: var
  integer                  ::i,unit,pos
  integer                  :: status
  logical                  :: bool
  character(len=255)       :: buffer
  if(present(default))variable=default
  name_=name;call upper_case(name_)
  inquire(file=file,exist=bool)
  if(.not.bool)IOinput=.false.
  if(bool)then
     unit=free_unit()
     open(unit,file=file)
     status=0
     var_search: do while(status>=0)
        read(unit,"(A255)",iostat=status)buffer
        pos=scan(buffer,"!");if(pos/=0)buffer=buffer(1:pos-1)

        var = get_input_variable(trim(buffer))
        if(var%name==name_)then
           read(var%value,*)variable
           exit var_search
        endif
     enddo var_search
     close(unit)
  endif
  call parse_cmd_variable(variable,name_)
  if(present(comment))then
     call append_to_input_list(variable,name_,comment)
  else
     call append_to_input_list(variable,name_)
  endif
  return
end subroutine ch_parse_input

subroutine l_parse_input(variable,name,file,default,comment)
  logical               :: variable
  logical,optional         :: default
  character(len=*)         :: name
  character(len=*),optional:: comment
  character(len=*)         :: file
  character(len=len(name)) :: name_
  type(input_variable)       :: var
  integer                  ::i,unit,pos
  integer                  :: status
  logical                  :: bool
  character(len=255)       :: buffer
  if(present(default))variable=default
  name_=name;call upper_case(name_)
  inquire(file=file,exist=bool)
  if(.not.bool)IOinput=.false.
  if(bool)then
     unit=free_unit()
     open(unit,file=file)
     status=0
     var_search: do while(status>=0)
        read(unit,"(A255)",iostat=status)buffer
        pos=scan(buffer,"!");if(pos/=0)buffer=buffer(1:pos-1)

        var = get_input_variable(trim(buffer))
        if(var%name==name_)then
           read(var%value,*)variable
           exit var_search
        endif
     enddo var_search
     close(unit)
  endif
  call parse_cmd_variable(variable,name_)
  if(present(comment))then
     call append_to_input_list(variable,name_,comment)
  else
     call append_to_input_list(variable,name_)
  endif
  return
end subroutine l_parse_input

!=====================1-dimension=====================================

subroutine iv_parse_input(variable,name,file,default,comment)
  integer,dimension(:)                       :: variable
  integer,dimension(size(variable)),optional :: default
  character(len=*)                           :: name
  character(len=*),optional:: comment
  character(len=*)                           :: file
  character(len=len(name))                   :: name_
  type(input_variable)                         :: var
  integer                                    ::i,unit,pos,j,ndim,ncount,nargs,pos0,iarg
  integer                                    :: status
  logical                                    :: bool
  character(len=255)                         :: buffer
  If(present(default))variable=default
  ndim=size(variable)
  name_=name;call upper_case(name_)
  inquire(file=file,exist=bool)
  if(.not.bool)IOinput=.false.
  if(bool)then
     unit=free_unit()
     open(unit,file=file)
     status=0
     var_search: do while(status>=0)
        read(unit,"(A255)",iostat=status)buffer
        pos=scan(buffer,"!");if(pos/=0)buffer=buffer(1:pos-1)

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
              read(var%args(iarg),*)variable(iarg)
           enddo
           exit var_search
        endif
     enddo var_search
     close(unit)
  endif
  call parse_cmd_variable(variable,name_)
  if(present(comment))then
     call append_to_input_list(variable,name_,comment)
  else
     call append_to_input_list(variable,name_)
  endif
  return
end subroutine iv_parse_input


subroutine dv_parse_input(variable,name,file,default,comment)
  real(8),dimension(:)                       :: variable
  real(8),dimension(size(variable)),optional :: default
  character(len=*)                           :: name
  character(len=*),optional:: comment
  character(len=*)                           :: file
  character(len=len(name))                   :: name_
  type(input_variable)                        :: var
  integer                                    ::i,unit,pos,j,ndim,ncount,nargs,pos0,iarg
  integer                                    :: status
  logical                                    :: bool
  character(len=255)                         :: buffer
  If(present(default))variable=default
  ndim=size(variable)
  name_=name;call upper_case(name_)
  inquire(file=file,exist=bool)
  if(.not.bool)IOinput=.false.
  if(bool)then
     unit=free_unit()
     open(unit,file=file)
     status=0
     var_search: do while(status>=0)
        read(unit,"(A255)",iostat=status)buffer
        pos=scan(buffer,"!");if(pos/=0)buffer=buffer(1:pos-1)
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
              read(var%args(iarg),*)variable(iarg)
           enddo
           exit var_search
        endif
     enddo var_search
     close(unit)
  endif
  call parse_cmd_variable(variable,name_)
  if(present(comment))then
     call append_to_input_list(variable,name_,comment)
  else
     call append_to_input_list(variable,name_)
  endif
  return
end subroutine dv_parse_input

subroutine chv_parse_input(variable,name,file,default,comment)
  character(len=*),dimension(:)                       :: variable
  character(len=*),dimension(size(variable)),optional :: default
  character(len=*)                           :: name
  character(len=*),optional:: comment
  character(len=*)                           :: file
  character(len=len(name))                   :: name_
  type(input_variable)                         :: var
  integer                                    ::i,unit,pos,j,ndim,ncount,nargs,pos0,iarg
  integer                                    :: status
  logical                                    :: bool
  character(len=255)                         :: buffer
  If(present(default))variable=default
  ndim=size(variable)
  name_=name;call upper_case(name_)
  inquire(file=file,exist=bool)
  if(.not.bool)IOinput=.false.
  if(bool)then
     unit=free_unit()
     open(unit,file=file)
     status=0
     var_search: do while(status>=0)
        read(unit,"(A255)",iostat=status)buffer
        pos=scan(buffer,"!");if(pos/=0)buffer=buffer(1:pos-1)
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
              read(var%args(iarg),*)variable(iarg)
           enddo
           exit var_search
        endif
     enddo var_search
     close(unit)
  endif
  call parse_cmd_variable(variable,name_)
  if(present(comment))then
     call append_to_input_list(variable,name_,comment)
  else
     call append_to_input_list(variable,name_)
  endif
  return
end subroutine chv_parse_input

subroutine lv_parse_input(variable,name,file,default,comment)
  logical,dimension(:)                       :: variable
  logical,dimension(size(variable)),optional :: default
  character(len=*)                           :: name
  character(len=*),optional:: comment
  character(len=*)                           :: file
  character(len=len(name))                   :: name_
  type(input_variable)                         :: var
  integer                                    ::i,unit,pos,j,ndim,ncount,nargs,pos0,iarg
  integer                                    :: status
  logical                                    :: bool
  character(len=255)                         :: buffer
  If(present(default))variable=default
  ndim=size(variable)
  name_=name;call upper_case(name_)
  inquire(file=file,exist=bool)
  if(.not.bool)IOinput=.false.
  if(bool)then
     unit=free_unit()
     open(unit,file=file)
     status=0
     var_search: do while(status>=0)
        read(unit,"(A255)",iostat=status)buffer
        pos=scan(buffer,"!");if(pos/=0)buffer=buffer(1:pos-1)
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
              read(var%args(iarg),*)variable(iarg)
           enddo
           exit var_search
        endif
     enddo var_search
     close(unit)
  endif
  call parse_cmd_variable(variable,name_)
  if(present(comment))then
     call append_to_input_list(variable,name_,comment)
  else
     call append_to_input_list(variable,name_)
  endif
  return
end subroutine lv_parse_input

! !=====================2-dimension=====================================

! subroutine im_parse_input(variable,name,file,default,comment)
!   integer,dimension(:,:)                     :: variable
!   integer,dimension(size(variable))          :: dummy_var
!   integer,dimension(size(variable)),optional :: default
!   character(len=*)                           :: name
!   character(len=*),optional:: comment
!   character(len=*)                           :: file
!   character(len=len(name))                   :: name_
!   type(input_variable)                         :: var
!   integer                                    ::i,unit,pos,j,ndim,nargs,pos0,iarg
!   integer                                    :: status
!   logical                                    :: bool
!   character(len=255)                         :: buffer
!   If(present(default))variable=transpose(reshape(default,shape(variable)))
!   ndim=size(variable)
!   name_=name;call upper_case(name_)
!   inquire(file=file,exist=bool)
!   if(.not.bool)IOinput=.false.
!   if(bool)then
!      unit=free_unit()
!      open(unit,file=file)
!      status=0
!      var_search: do while(status>=0)
!         read(unit,"(A255)",iostat=status)buffer
!         var = get_input_variable(trim(buffer))
!         if(var%name==name_)then
!            nargs=check_cmd_vector_size(ndim,var)
!            allocate(var%args(nargs))
!            iarg=0
!            pos0=0
!            do j=1,len(var%value)
!               if(var%value(j:j)==",")then
!                  iarg=iarg+1
!                  var%args(iarg)=var%value(pos0+1:j-1)
!                  pos0=j
!               endif
!            enddo
!            var%args(nargs)=var%value(pos0+1:)
!            do iarg=1,nargs
!               read(var%args(iarg),*)dummy_var(iarg)
!            enddo
!            variable=transpose(reshape(dummy_var,shape(variable)))
!            exit var_search
!         endif
!      enddo var_search
!   endif
!   call parse_cmd_variable(variable,name_)
!   close(unit)
!   if(present(comment))then
! call append_to_input_list(variable,name_,comment)
! else
! call append_to_input_list(variable,name_)
! endif
!   return
! end subroutine im_parse_input


! subroutine dm_parse_input(variable,name,file,default,comment)
!   real(8),dimension(:,:)                     :: variable
!   real(8),dimension(size(variable))          :: dummy_var
!   real(8),dimension(size(variable)),optional :: default
!   character(len=*)                           :: name
!   character(len=*),optional:: comment
!   character(len=*)                           :: file
!   character(len=len(name))                   :: name_
!   type(input_variable)                         :: var
!   integer                                    ::i,unit,pos,j,ndim,nargs,pos0,iarg
!   integer                                    :: status
!   logical                                    :: bool
!   character(len=255)                         :: buffer
!   If(present(default))variable=transpose(reshape(default,shape(variable)))
!   ndim=size(variable)
!   name_=name;call upper_case(name_)
!   inquire(file=file,exist=bool)
!   if(.not.bool)IOinput=.false.
!   if(bool)then
!      unit=free_unit()
!      open(unit,file=file)
!      status=0
!      var_search: do while(status>=0)
!         read(unit,"(A255)",iostat=status)buffer
!         var = get_input_variable(trim(buffer))
!         if(var%name==name_)then
!            nargs=check_cmd_vector_size(ndim,var)
!            allocate(var%args(nargs))
!            iarg=0
!            pos0=0
!            do j=1,len(var%value)
!               if(var%value(j:j)==",")then
!                  iarg=iarg+1
!                  var%args(iarg)=var%value(pos0+1:j-1)
!                  pos0=j
!               endif
!            enddo
!            var%args(nargs)=var%value(pos0+1:)
!            do iarg=1,nargs
!               read(var%args(iarg),*)dummy_var(iarg)
!            enddo
!            variable=transpose(reshape(dummy_var,shape(variable)))
!            exit var_search
!         endif
!      enddo var_search
!   endif
!   call parse_cmd_variable(variable,name_)
!   close(unit)
!   if(present(comment))then
! call append_to_input_list(variable,name_,comment)
! else
! call append_to_input_list(variable,name_)
! endif
!   return
! end subroutine dm_parse_input
