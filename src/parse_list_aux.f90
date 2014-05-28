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
     s(1:16) = '     0.d0     '
  else
     if(abs(r8) < 1.d0)then
        write ( s2, '(ES16.9)' ) r8
     else
        write ( s2, '(F16.9)' ) r8
     endif
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
