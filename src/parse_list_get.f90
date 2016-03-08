!========================0-dimension==================================
subroutine i_get_input_variable(variable,name,list)
  integer                   :: variable
  character(len=*)          :: name
  type(input_list),optional :: list
  integer                   :: i,counter,unit,size_
  type(input_node),pointer  :: c
  logical                   :: bool
  character(len=len(name)) :: name_
  name_=name;call upper_case(name_)
  if(present(list))then
     c => list%root%next
  else
     c => default_list%root%next
  endif
  counter = 0
  unit=free_unit()
  size_=default_list%size
  if(present(list))size_=list%size
  if(size_>0)then
     do
        if(.not.associated(c))exit
        counter=counter+1
        if(trim(c%name)==trim(name_))then
           variable=c%i(1)
           c=>null()
           return
        endif
        c => c%next
     enddo
     write(*,"(A)")"Can not find variable "//trim(name_)//" in the default input list" ; stop "exiting"
  else
     write(*,"(A)")"input list: empty"
     return
  endif
  c => null()
end subroutine i_get_input_variable

subroutine d_get_input_variable(variable,name,list)
  real(8)                   :: variable
  character(len=*)          :: name
  type(input_list),optional :: list
  integer                   :: i,counter,unit,size_
  type(input_node),pointer  :: c
  logical                   :: bool
  character(len=len(name)) :: name_
  name_=name;call upper_case(name_)
  if(present(list))then
     c => list%root%next
  else
     c => default_list%root%next
  endif
  counter = 0 
  unit=free_unit()
  size_=default_list%size
  if(present(list))size_=list%size
  if(size_>0)then
     do
        if(.not.associated(c))exit
        counter=counter+1
        if(trim(c%name)==trim(name_))then
           variable=c%d(1)
           c=>null()
           return
        endif
        c => c%next
     enddo
     write(*,"(A)")"Can not find variable "//trim(name_)//" in the default input list" ; stop "exiting"
  else
     write(*,"(A)")"input list: empty"
     return
  endif
  c => null()
end subroutine d_get_input_variable

subroutine l_get_input_variable(variable,name,list)
  logical                   :: variable
  character(len=*)          :: name
  type(input_list),optional :: list
  integer                   :: i,counter,unit,size_
  type(input_node),pointer  :: c
logical                   :: bool
  character(len=len(name)) :: name_
  name_=name;call upper_case(name_)
  if(present(list))then
     c => list%root%next
  else
     c => default_list%root%next
  endif
  counter = 0 
  unit=free_unit()
  size_=default_list%size
  if(present(list))size_=list%size
  if(size_>0)then
     do
        if(.not.associated(c))exit
        counter=counter+1
        if(trim(c%name)==trim(name_))then
           variable=c%l(1)
           c=>null()
           return
        endif
        c => c%next
     enddo
     write(*,"(A)")"Can not find variable "//trim(name_)//" in the default input list" ; stop "exiting"
  else
     write(*,"(A)")"input list: empty"
     return
  endif
  c => null()
end subroutine l_get_input_variable

subroutine ch_get_input_variable(variable,name,list)
  character(len=*)          :: variable
  character(len=*)          :: name
  type(input_list),optional :: list
  integer                   :: i,counter,unit,size_
  type(input_node),pointer  :: c
logical                   :: bool
  character(len=len(name)) :: name_
  name_=name;call upper_case(name_)
  if(present(list))then
     c => list%root%next
  else
     c => default_list%root%next
  endif
  counter = 0 
  unit=free_unit()
  size_=default_list%size
  if(present(list))size_=list%size
  if(size_>0)then
     do
        if(.not.associated(c))exit
        counter=counter+1
        if(trim(c%name)==trim(name_))then
           variable=c%ch(1)
           c=>null()
           return
        endif
        c => c%next
     enddo
     write(*,"(A)")"Can not find variable "//trim(name_)//" in the default input list" ; stop "exiting"
  else
     write(*,"(A)")"input list: empty"
     return
  endif
  c => null()
end subroutine ch_get_input_variable






!========================1-dimension==================================
subroutine iv_get_input_variable(variable,name,list)
  integer,dimension(:)      :: variable
  character(len=*)          :: name
  type(input_list),optional :: list
  integer                   :: i,counter,unit,size_
  type(input_node),pointer  :: c
logical                   :: bool
  character(len=len(name)) :: name_
  name_=name;call upper_case(name_)
  if(present(list))then
     c => list%root%next
  else
     c => default_list%root%next
  endif
  counter = 0 
  unit=free_unit()
  size_=default_list%size
  if(present(list))size_=list%size
  if(size_>0)then
     do
        if(.not.associated(c))exit
        counter=counter+1
        if(trim(c%name)==trim(name_))then
           variable=c%i(1:size(variable))
           if(size(variable)/=size(c%i))write(*,"(A)")"get_input_variable warning: variable has wrong dimensions"
           c=>null()
           return
        endif
        c => c%next
     enddo
     write(*,"(A)")"Can not find variable "//trim(name_)//" in the default input list" ; stop "exiting"
  else
     write(*,"(A)")"input list: empty"
     return
  endif
  c => null()
end subroutine iv_get_input_variable

subroutine dv_get_input_variable(variable,name,list)
  real(8),dimension(:)      :: variable
  character(len=*)          :: name
  type(input_list),optional :: list
  integer                   :: i,counter,unit,size_
  type(input_node),pointer  :: c
  logical                   :: bool
  character(len=len(name)) :: name_
  name_=name;call upper_case(name_)
  if(present(list))then
     c => list%root%next
  else
     c => default_list%root%next
  endif
  counter = 0 
  unit=free_unit()
  size_=default_list%size
  if(present(list))size_=list%size
  if(size_>0)then
     do
        if(.not.associated(c))exit
        counter=counter+1
        if(trim(c%name)==trim(name_))then
           variable=c%d(1:size(variable))
           if(size(variable)/=size(c%d))write(*,"(A)")"get_input_variable warning: variable has wrong dimensions"
           c=>null()
           return
        endif
        c => c%next
     enddo
     write(*,"(A)")"Can not find variable "//trim(name_)//" in the default input list" ; stop "exiting"
  else
     write(*,"(A)")"input list: empty"
     return
  endif
  c => null()
end subroutine dv_get_input_variable

subroutine lv_get_input_variable(variable,name,list)
  logical,dimension(:)      :: variable
  character(len=*)          :: name
  type(input_list),optional :: list
  integer                   :: i,counter,unit,size_
  type(input_node),pointer  :: c
  logical                   :: bool
  character(len=len(name)) :: name_
  name_=name;call upper_case(name_)
  if(present(list))then
     c => list%root%next
  else
     c => default_list%root%next
  endif
  counter = 0 
  unit=free_unit()
  size_=default_list%size
  if(present(list))size_=list%size
  if(size_>0)then
     do
        if(.not.associated(c))exit
        counter=counter+1
        if(trim(c%name)==trim(name_))then
           variable=c%l(1:size(variable))
           if(size(variable)/=size(c%l))write(*,"(A)")"get_input_variable warning: variable has wrong dimensions"
           c=>null()
           return
        endif
        c => c%next
     enddo
     write(*,"(A)")"Can not find variable "//trim(name_)//" in the default input list" ; stop "exiting"
  else
     write(*,"(A)")"input list: empty"
     return
  endif
  c => null()
end subroutine lv_get_input_variable

subroutine chv_get_input_variable(variable,name,list)
  character(len=*),dimension(:) :: variable
  character(len=*)              :: name
  type(input_list),optional     :: list
  integer                       :: i,counter,unit,size_
  type(input_node),pointer      :: c
  logical                       :: bool
  character(len=len(name)) :: name_
  name_=name;call upper_case(name_)
  if(present(list))then
     c => list%root%next
  else
     c => default_list%root%next
  endif
  counter = 0 
  unit=free_unit()
  size_=default_list%size
  if(present(list))size_=list%size
  if(size_>0)then
     do
        if(.not.associated(c))exit
        counter=counter+1
        if(trim(c%name)==trim(name_))then
           do i=1,size(variable)
              variable(i)=trim(c%ch(i))
           enddo
           if(size(variable)/=size(c%ch))write(*,"(A)")"get_input_variable warning: variable has wrong dimensions"
           c=>null()
           return
        endif
        c => c%next
     enddo
     write(*,"(A)")"Can not find variable "//trim(name_)//" in the default input list" ; stop "exiting"
  else
     write(*,"(A)")"input list: empty"
     return
  endif
  c => null()
end subroutine chv_get_input_variable
