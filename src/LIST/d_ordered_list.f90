module D_ORDERED_LIST
  implicit none
  private

  type,public :: node_object
     real(8) :: t
  end type node_object

  interface operator(<)
     module procedure less_than_obj
  end interface operator(<)

  interface operator(<=)
     module procedure less_or_equal_than_obj
  end interface operator(<=)

  interface operator(>)
     module procedure greater_than_obj
  end interface operator(>)

  interface operator(>=)
     module procedure greater_or_equal_than_obj
  end interface operator(>=)

  interface operator(==)
     module procedure equal_to_obj
  end interface operator(==)

  type,public :: d_list_node
     private
     type(node_object)       :: obj  !object value of the node (content of the box)
     type(d_list_node),pointer :: prev !link to prev box (chain)
     type(d_list_node),pointer :: next !link to next box (chain)
  end type d_list_node

  type,public :: d_linked_list
     private
     type(d_list_node),pointer :: root !head/root of the list\== list itself
  end type d_linked_list

  interface insert_element
     module procedure insert_element_before
  end interface insert_element

  public :: init_list
  public :: destroy_list
  public :: insert_element,insert_element_before,insert_element_after
  public :: remove_element
  public :: get_value
  public :: get_node
  public :: print_list,dump_list

  ! public :: operator(<)
  ! public :: operator(<=)
  ! public :: operator(>)
  ! public :: operator(>=)
  ! public :: operator(==)

contains        !some routine to perform simple operation on the lists

  function init_list result(new_list)
    type(d_linked_list) :: new_list
    allocate(new_list%root)
    new_list%root%prev => null()
    new_list%root%next => null()
  end function init_list


  subroutine destroy_list(list)
    type(d_linked_list),intent(in) :: list
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


  subroutine insert_element_before(list,obj)
    type(d_linked_list),intent(inout) :: list
    type(node_object),intent(in)    :: obj
    type(d_list_node),pointer         :: previous,current
    previous => list%root
    current  => previous%next
    do                               !traverse the list until obj < value (ordered list)
       if(.not.associated(current))exit !beginning of the list
       if(obj <= current%obj) exit
       previous => current
       current  => current%next
    end do
    allocate(previous%next)                !Create a new element in the list
    previous%next%obj = obj
    if(.not.associated(current))then !end of the list special case (current=>current%next)
       previous%next%next  => null()
       previous%next%prev  => previous     !the %prev of the new node go to previous node 
    else
       previous%next%next  => current      !the %next of the new node come to current
       previous%next%prev  => previous
       current%prev  => previous%next      !our %prev must go to the new node (previous%next)
    end if
  end subroutine insert_element_before

  subroutine insert_element_after(list,obj)
    type(d_linked_list),intent(inout) :: list
    type(node_object),intent(in)    :: obj
    type(d_list_node),pointer         :: previous,current
    previous => list%root
    current  => previous%next
    do                               !traverse the list until obj < value (ordered list)
       if(.not.associated(current))exit !beginning of the list
       if(obj < current%obj) exit
       previous => current
       current  => current%next
    end do
    allocate(previous%next)                !Create a new element in the list
    previous%next%obj = obj
    if(.not.associated(current))then !end of the list special case (current=>current%next)
       previous%next%next  => null()
       previous%next%prev  => previous     !the %prev of the new node go to previous node 
    else
       previous%next%next  => current      !the %next of the new node come to current
       previous%next%prev  => previous
       current%prev  => previous%next      !our %prev must go to the new node (previous%next)
    end if
  end subroutine insert_element_after



  subroutine remove_element(list,obj,found)
    type(d_linked_list),intent(inout) :: list
    type(node_object),intent(in)    :: obj
    logical,intent(out)             :: found
    type(d_list_node),pointer         :: previous,current
    previous => list%root
    current  => previous%next
    found    =  .false.
    do 
       if(found .OR. (.not.associated(current)))return
       if(obj == current%obj)then
          found=.true.
          exit       
       else
          previous => current   !previous%next
          current  => current%next
       endif
    end do
    !if found then delete this node and relocate the links:
    if(found)then
       previous%next     => current%next !reallocate skipping the deleted link
       current%next%prev => previous
       deallocate(current)           !free link
    end if
  end subroutine remove_element


  function get_value(node) result(value)
    type(d_list_node),intent(in) :: node
    type(node_object)          :: value
    value = node%obj
  end function get_value


  function get_node(list,value) result(node)
    type(d_linked_list),intent(in) :: list
    type(node_object)              :: value
    type(d_list_node),pointer      :: node
    type(d_list_node),pointer      :: current
    current => list%root%next
    do 
       if(.not.associated(current))then
          print*,"there are no node with this value"
          node=>list%root
          return
       endif
       if(current%obj == value)exit
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
       call print_obj(current%obj)
       current => current%next  !traverse list
    end do
  end subroutine print_list



  subroutine dump_list(list,vector)
    type(d_linked_list),intent(in) :: list
    type(d_list_node),pointer      :: current
    integer                        :: counter
    real(8),dimension(:)           :: vector
    current => list%root%next   !assume is associated,ie list exists
    counter = 0
    do
       if(.not.associated(current))exit
       counter=counter+1
       if(counter > size(vector))then
          print*,"error in data dump:return"
          return
       endif
       vector(counter) = current%obj%t
       current => current%next  !traverse list
    end do
  end subroutine dump_list

  function less_than_obj(O1,O2) result(boolean)
    type(node_object),intent(in) :: O1,O2
    logical                      :: boolean
    boolean = O1%t < O2%t
  end function less_than_obj
  !
  function less_or_equal_than_obj(O1,O2) result(boolean)
    type(node_object),intent(in) :: O1,O2
    logical                      :: boolean
    boolean = O1%t <= O2%t
  end function less_or_equal_than_obj
  !
  function greater_than_obj(O1,O2) result(boolean)
    type(node_object),intent(in) :: O1,O2
    logical                      :: boolean
    boolean = O1%t > O2%t
  end function greater_than_obj
  !
  function greater_or_equal_than_obj(O1,O2) result(boolean)
    type(node_object),intent(in) :: O1,O2
    logical                      :: boolean
    boolean = O1%t >= O2%t
  end function greater_or_equal_than_obj
  !
  function equal_to_obj(O1,O2) result(boolean)
    type(node_object),intent(in) :: O1,O2
    logical                      :: boolean
    boolean = O1%t == O2%t
  end function equal_to_obj
  !
  subroutine print_obj(obj)
    type(node_object) :: obj
    write(*,"(A,f12.7)")"value: ",obj%t
  end subroutine print_obj
end module D_ORDERED_LIST
