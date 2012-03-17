module D_UNORDERED_LIST
  implicit none
  private

  type,public :: d_list_node
     private
     real(8)          :: obj  !object value of the node (content of the box)
     type(d_list_node),pointer :: prev !link to prev box (chain)
     type(d_list_node),pointer :: next !link to next box (chain)
  end type d_list_node

  type,public :: d_linked_list
     integer                   :: size !size of the list
     type(d_list_node),pointer :: root !head/root of the list\== list itself
  end type d_linked_list

  public :: init_list
  public :: destroy_list
  public :: add_element
  public :: remove_element
  public :: get_value
  public :: get_node
  public :: print_list,dump_list

contains        !some routine to perform simple operation on the lists

  function init_list() result(new_list)
    type(d_linked_list) :: new_list
    allocate(new_list%root)
    new_list%root%prev => null()
    new_list%root%next => null()
    new_list%size=0
  end function init_list


  subroutine destroy_list(list)
    type(d_linked_list),intent(inout) :: list
    type(d_list_node),pointer     :: current
    do
       current => list%root%next         !current is the first node (root's next)
       if(.not.associated(current))exit  !empty list
       ! here we delete the current node by:
       !1 - pushing the previous node %next arrow to the next %next arrow (skip current)
       !2 - pushing the next node %prev arrow to the previous %prev arriw (skip current)
       !3 - nullify current
       current%prev%next => current%next !
       if(associated(current%next))current%next%prev => current%prev
       current%prev=>null()
       current%next=>null()
       deallocate(current)
    end do
    deallocate(list%root)
  end subroutine destroy_list


  subroutine add_element(list,obj)
    type(d_linked_list),intent(inout) :: list
    real(8) ,intent(in)      :: obj
    type(d_list_node),pointer         :: previous,current
    integer                           :: i
    previous => list%root
    current  => previous%next
    do i=1,list%size                    !traverse the list
       if(.not.associated(current))exit !beginning of the list
       previous => current
       current  => current%next
    end do
    allocate(previous%next)                !Create a new element in the list
    previous%next%obj = obj
    list%size=list%size+1
    if(.not.associated(current))then !end of the list special case (current=>current%next)
       previous%next%next  => null()
       previous%next%prev  => previous     !the %prev of the new node go to previous node 
    else
       previous%next%next  => current      !the %next of the new node come to current
       previous%next%prev  => previous
       current%prev  => previous%next      !our %prev must go to the new node (previous%next)
    end if
  end subroutine add_element


  subroutine remove_element(list,obj,n)
    type(d_linked_list),intent(inout) :: list
    real(8) ,intent(in)      :: obj
    integer,optional                  :: n
    integer                           :: i,pos
    type(d_list_node),pointer         :: previous,current
    previous => list%root
    current  => previous%next
    pos=list%size ; if(present(n))pos=n
    do i=1,pos
       if(.not.associated(current))return
       previous => current   !previous%next
       current  => current%next
    end do
    previous%next     => current%next !reallocate skipping the deleted link
    current%next%prev => previous
    deallocate(current)           !free link
    list%size=list%size-1
  end subroutine remove_element


  function get_value(node) result(value)
    type(d_list_node),intent(in) :: node
    real(8)           :: value
    value = node%obj
  end function get_value


  function get_node(list,n) result(node)
    type(d_linked_list),intent(in) :: list
    integer                        :: i,n
    type(d_list_node),pointer      :: node
    type(d_list_node),pointer      :: current
    current => list%root%next
    if(n>list%size)then
       print*,"error in get_node: n > list.size!"
       return
    endif
    do i=1,n
       if(.not.associated(current))then
          print*,"there are no node with this value"
          node=>list%root
          return
       endif
       current => current%next
    end do
    node => current
  end function get_node


  subroutine print_list(list)
    type(d_linked_list),intent(in) :: list
    type(d_list_node),pointer      :: current
    integer                        :: counter
    current => list%root%next   !assume is associated,ie list exists
    counter = 0
    do
       if(.not.associated(current))exit
       counter=counter+1
       write(*,"(A,I5,A)",advance='no')"element: ",counter," |"
       write(*,*)" value=",current%obj
       current => current%next  !traverse list
    end do
  end subroutine print_list



  subroutine dump_list(list,vector)
    type(d_linked_list),intent(in) :: list
    type(d_list_node),pointer      :: current
    integer                        :: i,N
    real(8),dimension(:)        :: vector
    current => list%root%next   !assume is associated,ie list exists
    if(size(vector) > list%size)then
       print*,"error in dump_list : vector.size > list.size"
       return
    endif
    do i=1,size(vector)
       if(.not.associated(current))exit
       vector(i) = current%obj
       current => current%next  !traverse list
    end do
  end subroutine dump_list

end module D_UNORDERED_LIST
