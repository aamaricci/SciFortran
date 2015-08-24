!=====================0-dimension=====================================
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

!=====================1-dimension=====================================

subroutine iv_parse_variable(variable,name,default)
  integer,dimension(:)                       :: variable
  integer,dimension(size(variable)),optional :: default
  character(len=*)                           :: name
  character(len=len(name))                  :: name_
  type(input_variable)                         :: var
  integer                                    :: i,j,ndim,ncount,nargs,pos0,iarg
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
        write(0,"(A,100I6)")"Variable "//trim(var%name)//" updated to ",(variable(iarg),iarg=1,ndim)
     endif
  enddo
end subroutine iv_parse_variable

subroutine dv_parse_variable(variable,name,default)
  real(8),dimension(:)                       :: variable
  real(8),dimension(size(variable)),optional :: default
  character(len=*)                           :: name
  character(len=len(name))                  :: name_
  type(input_variable)                         :: var
  integer                                    :: i,j,ndim,ncount,nargs,pos0,iarg
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
        write(0,"(A,100F18.9)")"Variable "//trim(var%name)//" updated to ",(variable(iarg),iarg=1,ndim)
     endif
  enddo
end subroutine dv_parse_variable

subroutine chv_parse_variable(variable,name,default)
  character(len=*),dimension(:)                       :: variable
  character(len=*),dimension(size(variable)),optional :: default
  character(len=*)                           :: name
  character(len=len(name))                  :: name_
  type(input_variable)                         :: var
  integer                                    :: i,j,ndim,ncount,nargs,pos0,iarg
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

subroutine lv_parse_variable(variable,name,default)
  logical,dimension(:)                       :: variable
  logical,dimension(size(variable)),optional :: default
  character(len=*)                           :: name
  character(len=len(name))                  :: name_
  type(input_variable)                         :: var
  integer                                    :: i,j,ndim,ncount,nargs,pos0,iarg
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

! !=====================2-dimension=====================================

! subroutine im_parse_variable(variable,name,default)
!   integer,dimension(:,:)                     :: variable
!   integer,dimension(size(variable))          :: dummy_var
!   integer,dimension(size(variable)),optional :: default
!   character(len=*)                           :: name
!   character(len=len(name))                   :: name_
!   type(input_variable)                         :: var
!   integer                                    :: i,j,ndim,nargs,pos0,iarg
!   If(present(default))variable=transpose(reshape(default,shape(variable)))
!   ndim=size(variable)
!   name_=name;call upper_case(name_)
!   do i=1,command_argument_count()
!      var = get_cmd_variable(i)
!      if(var%name==name_)then
!         nargs=check_cmd_vector_size(ndim,var)
!         allocate(var%args(nargs))
!         iarg=0
!         pos0=0
!         do j=1,len(var%value)
!            if(var%value(j:j)==",")then
!               iarg=iarg+1
!               var%args(iarg)=var%value(pos0+1:j-1)
!               pos0=j
!            endif
!         enddo
!         var%args(nargs)=var%value(pos0+1:)
!         do iarg=1,nargs
!            read(var%args(iarg),*)dummy_var(iarg)
!         enddo
!         variable=transpose(reshape(dummy_var,shape(variable)))
!         write(*,"(A,100I6)")" Variable "//trim(var%name)//" updated to ",(dummy_var(iarg),iarg=1,ndim)
!      endif
!   enddo
! end subroutine im_parse_variable

! subroutine dm_parse_variable(variable,name,default)
!   real(8),dimension(:,:)                     :: variable
!   real(8),dimension(size(variable))          :: dummy_var
!   real(8),dimension(size(variable)),optional :: default
!   character(len=*)                           :: name
!   character(len=len(name))                   :: name_
!   type(input_variable)                         :: var
!   integer                                    :: i,j,ndim,nargs,pos0,iarg
!   If(present(default))variable=transpose(reshape(default,shape(variable)))
!   ndim=size(variable)
!   name_=name;call upper_case(name_)
!   do i=1,command_argument_count()
!      var = get_cmd_variable(i)
!      if(var%name==name_)then
!         nargs=check_cmd_vector_size(ndim,var)
!         allocate(var%args(nargs))
!         iarg=0
!         pos0=0
!         do j=1,len(var%value)
!            if(var%value(j:j)==",")then
!               iarg=iarg+1
!               var%args(iarg)=var%value(pos0+1:j-1)
!               pos0=j
!            endif
!         enddo
!         var%args(nargs)=var%value(pos0+1:)
!         do iarg=1,nargs
!            read(var%args(iarg),*)dummy_var(iarg)
!         enddo
!         variable=transpose(reshape(dummy_var,shape(variable)))
!         write(*,"(A,100I6)")" Variable "//trim(var%name)//" updated to ",(dummy_var(iarg),iarg=1,ndim)
!      endif
!   enddo
! end subroutine dm_parse_variable
