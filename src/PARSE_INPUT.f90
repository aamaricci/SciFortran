  !########################################################################
  !Program  : PARSECMD
  !PURPOSE  : Declare all the common variables usually in use within codes
  !########################################################################
  include "PARSE_LIST_INPUT.f90"
  module PARSE_INPUT
    USE PARSE_LIST_INPUT
    implicit none
    private

    !cmd line variables:
    !=========================================================
    type input_variable
       character(len=64)                          :: name
       character(len=256)                         :: value
       character(len=20),dimension(:),allocatable :: args
    end type input_variable

    interface parse_cmd_variable
       module procedure &
            i_parse_variable,  d_parse_variable,  l_parse_variable,  ch_parse_variable,&
            iv_parse_variable, dv_parse_variable, lv_parse_variable, chv_parse_variable
    end interface parse_cmd_variable

    interface parse_input_variable
       module procedure &
            i_parse_input,  d_parse_input,  l_parse_input,  ch_parse_input, &
            iv_parse_input, dv_parse_input, lv_parse_input, chv_parse_input
    end interface parse_input_variable


    public  :: input_variable
    public  :: parse_cmd_variable
    public  :: parse_input_variable
    public  :: get_cmd_variable
    public  :: save_input_file
    public  :: parse_cmd_help
    logical :: IOinput=.true.


  contains

    !---------------------------------------------------------------------
    !PURPOSE: PARSE CMD LINE VARIABLES:
    !---------------------------------------------------------------------
    include 'parse_cmd_variable.f90'



    !---------------------------------------------------------------------
    !PURPOSE: PARSE INPUT FILE VARIABLES
    !---------------------------------------------------------------------
    include 'parse_input_variable.f90'




    !---------------------------------------------------------------------
    !PURPOSE:
    !---------------------------------------------------------------------
    subroutine save_input_file(file)
      character(len=*) :: file
      if(.not.IOinput)then
         write(*,*)"input file can not be found. dumped a default version in used."//trim(file)
         call print_input_list(trim(file))
         stop 
      else
         call print_input_list(trim(file))
      endif
    end subroutine save_input_file



    !---------------------------------------------------------------------
    !PURPOSE: ACCESSORY ROUTINES
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



    !---------------------------------------------------------------------
    !PURPOSE: AUXILIARY ROUTINES
    !---------------------------------------------------------------------
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
