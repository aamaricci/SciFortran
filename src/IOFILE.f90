!###############################################################
!     PROGRAM  : IOTOOLS
!     TYPE     : Module
!     PURPOSE  : SIMPLE PLOTTING/READING LIBRARY FOR FORTRAN 90/95
!     AUTHORS  : Adriano Amaricci (SISSA)
!###############################################################
module IOFILE
  USE COMMON_VARS
  implicit none
  private

  !file size to be stored automagically (in Kb)
  integer,public :: store_size=2048

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

  public :: file_size
  public :: file_length
  public :: file_info
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
       if(unit_>900) call error("ERROR free_unit: no unit free smaller than 900. Possible BUG")
    enddo
  end function free_unit




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
       call msg('Cannot read '//reg(file)//'. Skip file_size')
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
       call msg('Cannot read +'//reg(file)//'. Skip file_size')
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
    call msg("storing "//file)
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
    call warning("store size ="//trim(txtfy(size))//"Kb")
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

    call msg("deflate "//reg(filename)//reg(type))
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


end module IOFILE




