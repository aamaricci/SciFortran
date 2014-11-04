subroutine i_append_to_input_list(variable,name,comment)
  integer                  :: variable
  character(len=*)         :: name
  character(len=*),optional:: comment
  type(input_node),pointer :: p,c
  if(.not.default_list%status)call init_input_list()
  p => default_list%root
  c => p%next
  do                            !traverse the list until obj < value (ordered list)
     if(.not.associated(c))exit !empty list or beginning of the list
     p => c
     c => c%next
  end do
  allocate(p%next)                !Create a new element in the list
  !
  allocate(p%next%i(1))
  p%next%i(1)  = variable
  p%next%name= name
  p%next%type='i'
  p%next%comment=""
  if(present(comment))p%next%comment=trim(comment)
  !
  default_list%size=default_list%size+1
  if(.not.associated(c))then !end of the list special case (current=>current%next)
     p%next%next  => null()
  else
     p%next%next  => c      !the %next of the new node come to current
  end if
  p=>null()
  c=>null()
end subroutine i_append_to_input_list

subroutine d_append_to_input_list(variable,name,comment)
  real(8)                :: variable
  character(len=*)         :: name
  character(len=*),optional:: comment
  type(input_node),pointer :: p,c
  if(.not.default_list%status)call init_input_list()
  p => default_list%root
  c => p%next
  do                            !traverse the list until obj < value (ordered list)
     if(.not.associated(c))exit !empty list or beginning of the list
     p => c
     c => c%next
  end do
  allocate(p%next)                !Create a new element in the list
  !
  allocate(p%next%d(1))
  p%next%d(1) = variable
  p%next%name= name
  p%next%type='d'
  p%next%comment=""
  if(present(comment))p%next%comment=trim(comment)
  !
  default_list%size=default_list%size+1
  if(.not.associated(c))then !end of the list special case (current=>current%next)
     p%next%next  => null()
  else
     p%next%next  => c      !the %next of the new node come to current
  end if
  p=>null()
  c=>null()
end subroutine d_append_to_input_list

subroutine l_append_to_input_list(variable,name,comment)
  logical                  :: variable
  character(len=*)         :: name
  character(len=*),optional:: comment
  type(input_node),pointer :: p,c
  if(.not.default_list%status)call init_input_list()
  p => default_list%root
  c => p%next
  do                            !traverse the list until obj < value (ordered list)
     if(.not.associated(c))exit !empty list or beginning of the list
     p => c
     c => c%next
  end do
  allocate(p%next)                !Create a new element in the list
  !
  allocate(p%next%l(1))
  p%next%l(1) = variable
  p%next%name= name
  p%next%type='l'
  p%next%comment=""
  if(present(comment))p%next%comment=trim(comment)
  !
  default_list%size=default_list%size+1
  if(.not.associated(c))then !end of the list special case (current=>current%next)
     p%next%next  => null()
  else
     p%next%next  => c      !the %next of the new node come to current
  end if
  p=>null()
  c=>null()
end subroutine l_append_to_input_list

subroutine ch_append_to_input_list(variable,name,comment)
  character(len=*)         :: variable
  character(len=*)         :: name
  character(len=*),optional:: comment
  type(input_node),pointer :: p,c
  if(.not.default_list%status)call init_input_list()
  p => default_list%root
  c => p%next
  do                            !traverse the list until obj < value (ordered list)
     if(.not.associated(c))exit !empty list or beginning of the list
     p => c
     c => c%next
  end do
  allocate(p%next)                !Create a new element in the list
  !
  allocate(p%next%ch(1))
  p%next%ch(1) = variable
  p%next%name= name
  p%next%type='ch'
  p%next%comment=""
  if(present(comment))p%next%comment=trim(comment)
  !
  default_list%size=default_list%size+1
  if(.not.associated(c))then !end of the list special case (current=>current%next)
     p%next%next  => null()
  else
     p%next%next  => c      !the %next of the new node come to current
  end if
  p=>null()
  c=>null()
end subroutine ch_append_to_input_list

!========================1-dimension==================================

subroutine iv_append_to_input_list(variable,name,comment)
  integer,dimension(:)     :: variable
  character(len=*)         :: name
  character(len=*),optional:: comment
  type(input_node),pointer :: p,c
  if(.not.default_list%status)call init_input_list()
  p => default_list%root
  c => p%next
  do                            !traverse the list until obj < value (ordered list)
     if(.not.associated(c))exit !empty list or beginning of the list
     p => c
     c => c%next
  end do
  allocate(p%next)                !Create a new element in the list
  !
  allocate(p%next%i(size(variable)))
  p%next%i   = variable
  p%next%name= name
  p%next%type='i'
  p%next%comment=""
  if(present(comment))p%next%comment=trim(comment)
  !
  default_list%size=default_list%size+1
  if(.not.associated(c))then !end of the list special case (current=>current%next)
     p%next%next  => null()
  else
     p%next%next  => c      !the %next of the new node come to current
  end if
  p=>null()
  c=>null()
end subroutine iv_append_to_input_list

subroutine dv_append_to_input_list(variable,name,comment)
  real(8),dimension(:)     :: variable
  character(len=*)         :: name
  character(len=*),optional:: comment
  type(input_node),pointer :: p,c
  if(.not.default_list%status)call init_input_list()
  p => default_list%root
  c => p%next
  do                            !traverse the list until obj < value (ordered list)
     if(.not.associated(c))exit !empty list or beginning of the list
     p => c
     c => c%next
  end do
  allocate(p%next)                !Create a new element in the list
  !
  allocate(p%next%d(size(variable)))
  p%next%d   = variable
  p%next%name= name
  p%next%type='d'
  p%next%comment=""
  if(present(comment))p%next%comment=trim(comment)
  !
  default_list%size=default_list%size+1
  if(.not.associated(c))then !end of the list special case (current=>current%next)
     p%next%next  => null()
  else
     p%next%next  => c      !the %next of the new node come to current
  end if
  p=>null()
  c=>null()
end subroutine dv_append_to_input_list

subroutine lv_append_to_input_list(variable,name,comment)
  logical,dimension(:)     :: variable
  character(len=*)         :: name
  character(len=*),optional:: comment
  type(input_node),pointer :: p,c
  if(.not.default_list%status)call init_input_list()
  p => default_list%root
  c => p%next
  do                            !traverse the list until obj < value (ordered list)
     if(.not.associated(c))exit !empty list or beginning of the list
     p => c
     c => c%next
  end do
  allocate(p%next)                !Create a new element in the list
  !
  allocate(p%next%l(size(variable)))
  p%next%l   = variable
  p%next%name= name
  p%next%type='l'
  p%next%comment=""
  if(present(comment))p%next%comment=trim(comment)
  !
  default_list%size=default_list%size+1
  if(.not.associated(c))then !end of the list special case (current=>current%next)
     p%next%next  => null()
  else
     p%next%next  => c      !the %next of the new node come to current
  end if
  p=>null()
  c=>null()
end subroutine lv_append_to_input_list

subroutine chv_append_to_input_list(variable,name,comment)
  character(len=*),dimension(:) :: variable
  character(len=*)              :: name
  character(len=*),optional:: comment
  type(input_node),pointer      :: p,c
  integer :: i
  if(.not.default_list%status)call init_input_list()
  p => default_list%root
  c => p%next
  do                            !traverse the list until obj < value (ordered list)
     if(.not.associated(c))exit !empty list or beginning of the list
     p => c
     c => c%next
  end do
  allocate(p%next)                !Create a new element in the list
  !
  allocate(p%next%ch(size(variable)))
  do i=1,size(variable)
     p%next%ch(i) = trim(variable(i))
  enddo
  p%next%name= name
  p%next%type='ch'
  p%next%comment=""
  if(present(comment))p%next%comment=trim(comment)
  !
  default_list%size=default_list%size+1
  if(.not.associated(c))then !end of the list special case (current=>current%next)
     p%next%next  => null()
  else
     p%next%next  => c      !the %next of the new node come to current
  end if
  p=>null()
  c=>null()
end subroutine chv_append_to_input_list

! !========================2-dimension==================================

! subroutine im_append_to_input_list(variable,name,comment)
!   integer,dimension(:,:)                               :: variable
!   character(len=*)                                     :: name
! character(len=*),optional:: comment
!   type(input_node),pointer                             :: p,c
!   if(.not.default_list%status)call init_input_list()
!   p => default_list%root
!   c => p%next
!   do                            !traverse the list until obj < value (ordered list)
!      if(.not.associated(c))exit !empty list or beginning of the list
!      p => c
!      c => c%next
!   end do
!   allocate(p%next)                !Create a new element in the list
!   !
!   allocate(p%next%i(size(variable)))
!   p%next%i   = pack(variable,.true.)
!   p%next%name= name
!   p%next%type='i'
! p%next%comment=""
! if(present(comment))p%next%comment=trim(comment)
!   !
!   default_list%size=default_list%size+1
!   if(.not.associated(c))then !end of the list special case (current=>current%next)
!      p%next%next  => null()
!   else
!      p%next%next  => c      !the %next of the new node come to current
!   end if
!   p=>null()
!   c=>null()
! end subroutine im_append_to_input_list

! subroutine dm_append_to_input_list(variable,name,comment)
!   real(8),dimension(:,:)                               :: variable
!   character(len=*)                                     :: name
!character(len=*),optional:: comment
!   type(input_node),pointer                             :: p,c
!   if(.not.default_list%status)call init_input_list()
!   p => default_list%root
!   c => p%next
!   do                            !traverse the list until obj < value (ordered list)
!      if(.not.associated(c))exit !empty list or beginning of the list
!      p => c
!      c => c%next
!   end do
!   allocate(p%next)                !Create a new element in the list
!   !
!   allocate(p%next%d(size(variable)))
!   p%next%d   = pack(variable,.true.)
!   p%next%name= name
!   p%next%type='d'
! p%next%comment=""
! if(present(comment))p%next%comment=trim(comment)
!   !
!   default_list%size=default_list%size+1
!   if(.not.associated(c))then !end of the list special case (current=>current%next)
!      p%next%next  => null()
!   else
!      p%next%next  => c      !the %next of the new node come to current
!   end if
!   p=>null()
!   c=>null()
! end subroutine dm_append_to_input_list
