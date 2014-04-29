!###############################################################
!     PROGRAM  : IOTOOLS
!     TYPE     : Module
!     PURPOSE  : SIMPLE PLOTTING/READING LIBRARY FOR FORTRAN 90/95
!     AUTHORS  : Adriano Amaricci (SISSA)
!###############################################################
module IOFILE
  USE MPI_VARS
  implicit none
  private

  !file size to be stored automagically (in Kb)
  integer,public :: store_size=2048

  interface txtfy
     module procedure i_to_ch,r_to_ch,c_to_ch,l_to_ch
  end interface txtfy

  interface reg
     module procedure reg_filename
  end interface reg

  interface txtfit
     module procedure reg_filename
  end interface txtfit

  interface txtcut
     module procedure reg_filename
  end interface txtcut

  interface create_dir
     module procedure create_data_dir
  end interface create_dir

  public :: txtfy
  public :: file_size
  public :: file_length
  public :: file_info
  public :: free_unit,free_units
  public :: data_open
  public :: data_store
  public :: set_store_size
  public :: reg_filename,reg,txtfit,txtcut
  public :: create_data_dir,create_dir
  public :: close_file
  public :: get_filename
  public :: get_filepath

contains


  !+-----------------------------------------------------------------+
  !PURPOSE  : 
  !+-----------------------------------------------------------------+
  function reg_filename(file) result(reg)
    character(len=*)                                   :: file    
    character(len=len_trim(trim(adjustl(trim(file))))) :: reg
    reg=trim(adjustl(trim(file)))
  end function reg_filename

  function get_filename(string) result(fname)
    character(len=*) :: string
    character(len=len_trim(string)) :: fname
    integer :: i,slen
    slen=len_trim(string)
    do i=slen,1,-1
       if(string(i:i)== '/')exit
    enddo
    fname=string(i+1:slen)
  end function get_filename

  function get_filepath(string) result(pname)
    character(len=*) :: string
    character(len=len_trim(string)) :: pname
    integer :: i,slen
    slen=len_trim(string)
    do i=slen,1,-1
       if(string(i:i)== '/')exit
    enddo
    pname=string(1:i)
  end function get_filepath

  subroutine close_file(pname)
    character(len=*) :: pname
    open(10,file=reg(pname),position="APPEND")
    write(10,*)""
    close(10)
  end subroutine close_file

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

  function free_units(n) result(unit)
    integer :: n
    integer :: unit(n)
    integer :: i,unit_,ios
    logical :: is_it_opened
    unit_=100
    do i=1,n
       search: do 
          unit_=unit_+1
          inquire(unit=unit_,OPENED=is_it_opened,iostat=ios)
          if(.not.is_it_opened.AND.ios==0)exit search
          if(unit_>900) stop "ERROR free_unit: no unit free smaller than 900. Possible BUG"
       enddo search
       unit(i)=unit_
    enddo
  end function free_units




  !+-----------------------------------------------------------------+
  !PURPOSE  : 
  !+-----------------------------------------------------------------+
  function file_size(file,printf) result(size)
    integer               :: size,status
    character(len=*)      :: file
    integer,dimension(13) :: buff
    logical,optional      :: printf
    logical               :: control
    inquire(file=reg(file),exist=control)
    if(.not.control)then
       write(*,*) 'Cannot read '//reg(file)//'. Skip file_size'
       return
    endif
    open(10,file=reg(file))
    call fstat(10,buff,status)
    size=nint(dble(buff(8))/dble(1024))
    if(present(printf).AND.printf.eqv..true.)&
         write(*,"(A,A,A,f9.6,A)")"file: **",reg(file),"** is ",size," Kb"
  end function file_size



  !******************************************************************
  !******************************************************************
  !******************************************************************



  !+-----------------------------------------------------------------+
  !PURPOSE  : 
  !+-----------------------------------------------------------------+
  function file_info(file)
    integer               :: file_info
    character(len=*)      :: file
    integer,dimension(13) :: buff
    integer               :: status,ifile
    logical               :: IOfile
    inquire(file=reg(file),exist=IOfile)
    if(.not.IOfile)then
       print*,'Cannot read ',reg(file),': skip file_size'
       file_info=0
       return
    endif
    open(10,file=reg(file))
    call fstat(10,buff,status)
    if(status == 0)then
       WRITE (*, FMT="('Device ID:',               T30, I19)") buff(1)
       WRITE (*, FMT="('Inode number:',            T30, I19)") buff(2)
       WRITE (*, FMT="('File mode (octal):',       T30, O19)") buff(3)
       WRITE (*, FMT="('Number of links:',         T30, I19)") buff(4)
       WRITE (*, FMT="('Owner''s uid:',            T30, I19)") buff(5)
       WRITE (*, FMT="('Owner''s gid:',            T30, I19)") buff(6)
       WRITE (*, FMT="('Device where located:',    T30, I19)") buff(7)
       WRITE (*, FMT="('File size:',               T30, I19)") buff(8)
       WRITE (*, FMT="('Last access time:',        T30, A19)") buff(9)
       WRITE (*, FMT="('Last modification time',   T30, A19)") buff(10)
       WRITE (*, FMT="('Last status change time:', T30, A19)") buff(11)
       WRITE (*, FMT="('Preferred block size:',    T30, I19)") buff(12)
       WRITE (*, FMT="('No. of blocks allocated:', T30, I19)") buff(13)
    endif
    close(10)
  end function file_info



  !******************************************************************
  !******************************************************************
  !******************************************************************



  !+-----------------------------------------------------------------+
  !PURPOSE  : 
  !+-----------------------------------------------------------------+
  function file_length(file) result(lines)
    integer           :: lines
    character(len=*)  :: file
    integer           :: ifile,ierr,pos
    logical           :: IOfile,bool,bool1,bool2
    character(len=256)::buffer
    inquire(file=reg(file),exist=IOfile)
    if(.not.IOfile)then
       inquire(file=reg(file)//".gz",exist=IOfile)
       if(IOfile)call data_open(reg(file))
    endif
    lines=0
    if(.not.IOfile)then
       write(*,*) 'Cannot read +'//reg(file)//'. Skip file_size'
       return
    endif
    open(99,file=reg(file))
    ierr=0
    do while(ierr==0)
       lines=lines+1
       read(99,*,iostat=ierr)buffer
       bool1=scan(buffer,"#").ne.0
       bool2=len_trim(buffer).eq.0       
       if(bool1 .OR. bool2)lines=lines-1
    enddo
    lines=lines-1
    write(*,'(A,I9,A)') 'there are', lines,' lines in +'//reg(file)
    rewind(99)
    close(99)
  end function file_length



  !******************************************************************
  !******************************************************************
  !******************************************************************






  !+-----------------------------------------------------------------+
  !PURPOSE  : 
  !+-----------------------------------------------------------------+
  subroutine data_store(file,size)
    character(len=*)  :: file
    integer,optional  :: size
    logical           :: control
    character(len=9)  :: csize 
    integer           :: cstatus,fsize
    write(*,*) "storing "//file
    inquire(file=reg(file),exist=control)
    if(control)then
       fsize=store_size;if(present(size))fsize=size
       if(file_size(reg(file))>fsize)&
            call system("gzip -fv "//reg(file))
    endif
  end subroutine data_store

  subroutine set_store_size(size)
    integer :: size
    store_size=size
    write(*,*)"store size ="//trim(txtfy(size))//"Kb"
  end subroutine set_store_size

  !******************************************************************
  !******************************************************************
  !******************************************************************



  !+-----------------------------------------------------------------+
  !PURPOSE  : 
  !+-----------------------------------------------------------------+
  subroutine data_open(filename,tar)
    character(len=*)               :: filename
    character(len=10)              :: type
    logical                        :: compressed,control
    logical,optional,intent(out)   :: tar
    type=".gz"
    inquire(file=reg(filename),exist=control)
    if(control)then             !If exist return (no untar)
       if(present(tar))tar  = .false.
       return
    else                        !else search the correct compress format
       inquire(file=reg(filename)//reg(type),exist=compressed)
       if(present(tar))tar =compressed
       if(.not.compressed)return
    endif

    write(*,*) "deflate "//reg(filename)//reg(type)
    call system("gunzip "//reg(filename)//(type))
    return
  end subroutine data_open


  !******************************************************************
  !******************************************************************
  !******************************************************************




  !+-----------------------------------------------------------------+
  !PURPOSE  : 
  !+-----------------------------------------------------------------+
  subroutine create_data_dir(dir_name,id)
    character(len=*),optional :: dir_name
    character(len=256)        :: name
    logical                   :: control
    integer,optional :: id
    integer          :: id_
    id_=0         ;if(present(id))id_=id
    name="DATAsrc";if(present(dir_name))name=dir_name
    if(mpiID==id_)then
       call system("mkdir -v "//reg(name))
    endif
  end subroutine create_data_dir


  !******************************************************************
  !******************************************************************
  !******************************************************************


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





  !******************************************************************
  !******************************************************************
  !******************************************************************

  subroutine i4_to_s_left ( i4, s )
    !! I4_TO_S_LEFT converts an I4 to a left-justified string.
    !  Example:
    !    Assume that S is 6 characters long:
    !        I4  S
    !         1  1
    !        -1  -1
    !         0  0
    !      1952  1952
    !    123456  123456
    !   1234567  ******  <-- Not enough room!
    !  Parameters:
    !    Input, integer ( kind = 4 ) I4, an integer to be converted.
    !    Output, character ( len = * ) S, the representation of the integer.
    !    The integer will be left-justified.  If there is not enough space,
    !    the string will be filled with stars.
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
    !! R8_TO_S_LEFT writes an R8 into a left justified string.
    !    An R8 is a real ( kind = 8 ) value.
    !    A 'G14.6' format is used with a WRITE statement.
    !  Parameters:
    !    Input, real ( kind = 8 ) R8, the number to be written into the string.
    !    Output, character ( len = * ) S, the string into which
    !    the real number is to be written.  If the string is less than 14
    !    characters long, it will will be returned as a series of asterisks.
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
    !! DIGIT_TO_CH returns the character representation of a decimal digit.
    !    Instead of CHAR, we now use the ACHAR function, which
    !    guarantees the ASCII collating sequence.
    !  Example:
    !    DIGIT   CH
    !    -----  ---
    !      0    '0'
    !      1    '1'
    !    ...    ...
    !      9    '9'
    !     17    '*'
    !  Parameters:
    !    Input, integer ( kind = 4 ) DIGIT, the digit value between 0 and 9.
    !    Output, character CH, the corresponding character.
    character :: ch
    integer   :: digit
    if ( 0 <= digit .and. digit <= 9 ) then
       ch = achar ( digit + 48 )
    else
       ch = '*'
    end if
  end subroutine digit_to_ch

end module IOFILE




