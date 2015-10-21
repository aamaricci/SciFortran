module PARSE_LIST_INPUT
  implicit none
  private

  type input_node
     integer,dimension(:),pointer            :: i
     real(8),dimension(:),pointer            :: d
     logical,dimension(:),pointer            :: l
     character(len=100),dimension(:),pointer :: ch
     !
     character(len=3)                        :: type
     character(len=100)                      :: name
     character(len=512)                      :: comment
     type(input_node),pointer                :: next !link to next box
  end type input_node

  type input_list
     logical                   :: status=.false.
     integer                   :: size 
     type(input_node),pointer  :: root
  end type input_list


  interface append_to_input_list
     module procedure &
          i_append_to_input_list,d_append_to_input_list,ch_append_to_input_list,l_append_to_input_list,&
          iv_append_to_input_list,dv_append_to_input_list,chv_append_to_input_list,lv_append_to_input_list
  end interface append_to_input_list


  interface txtfy
     module procedure i_to_ch,r_to_ch,c_to_ch,l_to_ch
  end interface txtfy

  public :: input_list
  public :: init_input_list
  public :: destroy_input_list
  public :: size_input_list
  public :: append_to_input_list
  public :: print_input_list
  public :: print_help_list

  type(input_list)   :: default_list
  character(len=255) :: p_buffer
  character(len=7)   :: file_status
  integer,parameter  :: pos_comment=46 !72

contains  



  !+------------------------------------------------------------------+
  !PURPOSE: init the input list
  !+------------------------------------------------------------------+
  subroutine init_input_list(list)
    type(input_list),optional :: list
    if(present(list))then
       allocate(list%root)    
       list%size=0
       list%status=.true.
       list%root%next=>null()
    else
       allocate(default_list%root)    
       default_list%size=0
       default_list%status=.true.
       default_list%root%next=>null()
    endif
  end subroutine init_input_list

  !+------------------------------------------------------------------+
  !PURPOSE: delete the list
  !+------------------------------------------------------------------+
  subroutine destroy_input_list(list)
    type(input_list),optional :: list
    type(input_node),pointer  :: p,c
    if(present(list))then
       do
          p => list%root
          c => p%next
          if(.not.associated(c))exit  !empty list
          p%next => c%next !
          c%next=>null()
          deallocate(c)
       end do
       list%status=.false.
    else
       do
          p => default_list%root
          c => p%next
          if(.not.associated(c))exit  !empty list
          p%next => c%next !
          c%next=>null()
          deallocate(c)
       end do
       default_list%status=.false.
    endif
  end subroutine destroy_input_list


  !+------------------------------------------------------------------+
  !PURPOSE: get list size
  !+------------------------------------------------------------------+
  function size_input_list(list) result(size)
    type(input_list),optional :: list
    integer                   :: size
    size=default_list%size
    if(present(list))size=list%size
  end function size_input_list



  !+------------------------------------------------------------------+
  !PURPOSE: print the list to file
  !+------------------------------------------------------------------+
  subroutine print_input_list(file,list)
    type(input_list),optional :: list
    character(len=*)          :: file
    integer                   :: i,counter,unit
    type(input_node),pointer  :: c
    logical                   :: bool
    if(present(list))then
       c => list%root%next
    else
       c => default_list%root%next
    endif
    counter = 0 
    unit=free_unit()
    file_status='replace'
    if(default_list%size>0)then
       do
          if(.not.associated(c))exit
          counter=counter+1
          call print_input_node(c,file)
          c => c%next
       enddo
    else
       write(*,*)"input list: empty"
       return
    endif
    c => null()
  end subroutine print_input_list



  !+------------------------------------------------------------------+
  !PURPOSE: print the list to file
  !+------------------------------------------------------------------+
  subroutine print_help_list(list)
    type(input_list),optional :: list
    integer                   :: i,counter,unit
    type(input_node),pointer  :: c
    logical                   :: bool
    if(present(list))then
       c => list%root%next
    else
       c => default_list%root%next
    endif
    counter = 0 
    if(default_list%size>0)then
       do
          if(.not.associated(c))exit
          counter=counter+1
          call help_input_node(c)
          c => c%next
       enddo
    else
       write(*,*)"input list empty: no help."
       return
    endif
    c => null()
  end subroutine print_help_list



  !+------------------------------------------------------------------+
  !PURPOSE: add input to the list, print to file, aux routines
  !+------------------------------------------------------------------+
  include 'parse_list_append.f90'
  include 'parse_list_print.f90'
  include 'parse_list_aux.f90'

end module PARSE_LIST_INPUT
