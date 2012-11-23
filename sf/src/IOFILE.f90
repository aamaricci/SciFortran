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
    open(10,file=trim(adjustl(trim(pname))),access="APPEND")
    write(10,*)""
    close(10)
  end subroutine close_file


  !+-----------------------------------------------------------------+
  !PURPOSE  : 
  !+-----------------------------------------------------------------+
  function file_size(file,printf) result(size)
    integer               :: size,status
    character(len=*)      :: file
    integer,dimension(13) :: buff
    logical,optional      :: printf
    logical               :: control
    inquire(file=trim(adjustl(trim(file))),exist=control)
    if(.not.control)then
       call msg('Cannot read '//trim(adjustl(trim(file)))//'. Skip file_size')
       return
    endif
    open(10,file=trim(adjustl(trim(file))))
    call fstat(10,buff,status)
    size=nint(dble(buff(8))/dble(1024))
    if(present(printf).AND.printf.eqv..true.)&
         write(*,"(A,A,A,f9.6,A)")"file: **",trim(adjustl(trim(file))),"** is ",size," Kb"
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
    inquire(file=trim(adjustl(trim(file))),exist=IOfile)
    if(.not.IOfile)then
       print*,'Cannot read ',trim(adjustl(trim(file))),': skip file_size'
       file_info=0
       return
    endif
    open(10,file=trim(adjustl(trim(file))))
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
    inquire(file=trim(adjustl(trim(file))),exist=IOfile)
    if(.not.IOfile)then
       inquire(file=trim(adjustl(trim(file)))//".gz",exist=IOfile)
       if(IOfile)call data_open(trim(adjustl(trim(file))))
    endif
    lines=0
    if(.not.IOfile)then
       call msg('Cannot read +'//trim(adjustl(trim(file)))//'. Skip file_size')
       return
    endif
    open(99,file=trim(adjustl(trim(file))))
    ierr=0
    do while(ierr==0)
       lines=lines+1
       read(99,*,iostat=ierr)buffer
       bool1=scan(buffer,"#").ne.0
       bool2=len_trim(buffer).eq.0       
       if(bool1 .OR. bool2)lines=lines-1
    enddo
    lines=lines-1
    write(*,'(A,I9,A)') 'there are', lines,' lines in +'//trim(adjustl(trim(file)))
    rewind(99)
    close(99)
  end function file_length



  !******************************************************************
  !******************************************************************
  !******************************************************************



  !+-----------------------------------------------------------------+
  !PURPOSE  : 
  !+-----------------------------------------------------------------+
  function reg_filename(file) result(reg)
    character(len=*)                                   :: file    
    character(len=len_trim(trim(adjustl(trim(file))))) :: reg
    reg=trim(adjustl(trim(file)))
  end function reg_filename


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
    !Check file exists:
    inquire(file=reg_filename(file),exist=control)
    if(control)then
       !Check if store_size is environment variable:
       call get_environment_variable("STORE_SIZE",csize,STATUS=cstatus)
       if(cstatus/=1)read(csize,"(I9)")store_size
       fsize=store_size;if(present(size))fsize=size
       if(file_size(reg_filename(file))>fsize)&
            call system("gzip -fv "//trim(adjustl(trim(file))))
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
    inquire(file=reg_filename(filename),exist=control)
    if(control)then             !If exist return (no untar)
       if(present(tar))tar  = .false.
       return
    else                        !else search the correct compress format
       inquire(file=reg_filename(filename)//reg_filename(type),exist=compressed)
       if(present(tar))tar =compressed
       if(.not.compressed)return
    endif

    call msg("deflate "//reg_filename(filename)//reg_filename(type))
    call system("gunzip "//reg_filename(filename)//reg_filename(type))
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
       ! control = check_data_dir(name)
       ! if(control)then
       !    call warning("directory +"//trim(adjustl(trim(name)))//" exists")
       !    return
       ! else
       call system("mkdir -v "//trim(adjustl(trim(name))))
       ! endif
    endif
  end subroutine create_data_dir


  !******************************************************************
  !******************************************************************
  !******************************************************************



  !+-----------------------------------------------------------------+
  !PURPOSE  : 
  !+-----------------------------------------------------------------+
  function check_data_dir(dir_name) result(logic)
    character(len=*),optional :: dir_name
    logical                   :: logic
    call system("rm -f dir_exist")
    call system("if [ -d "//trim(adjustl(trim(dir_name)))//" ]; then echo > dir_exist; fi")
    inquire(file="dir_exist",EXIST=logic)
    call system("rm -f dir_exist")
  end function check_data_dir



  !******************************************************************
  !******************************************************************
  !******************************************************************


end module IOFILE




