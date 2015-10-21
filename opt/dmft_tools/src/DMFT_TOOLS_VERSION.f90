MODULE DMFT_TOOLS_VERSION
  implicit none
  include "dmft_tools_version.inc"

contains

  !+-------------------------------------------------------------------+
  !PURPOSE  : print actual version of the software (if any)
  !+-------------------------------------------------------------------+
  subroutine dmftt_version(revision)
    character(len=*)                         :: revision
    integer(4),dimension(8)                  :: dummy
    integer(4)                               :: year
    integer(4)                               :: mese
    integer(4)                               :: day
    integer(4)                               :: h
    integer(4)                               :: m
    integer(4)                               :: s
    integer(4)                               :: ms
    character(len=9),parameter,dimension(12) :: month = (/ &
         'January  ', 'February ', 'March    ', 'April    ', &
         'May      ', 'June     ', 'July     ', 'August   ', &
         'September', 'October  ', 'November ', 'December ' /)
    write(*,"(A)")("SCIFOR VERSION (GIT): "//trim(adjustl(trim(sf_version))))
    write(*,"(A)")("CODE VERSION (GIT): "//trim(adjustl(trim(revision))))
    call date_and_time(values=dummy)
    year = dummy(1)
    mese = dummy(2)
    day  = dummy(3)
    h    = dummy(5)
    m    = dummy(6)
    s    = dummy(7)
    ms   = dummy(8)
    write(*,"(A,i2,1x,a,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3)")&
         "Timestamp: +",day,trim(month(mese)),year, h,':',m,':',s,'.',ms
    write(*,*)""
    open(10,file="version.inc")
    write(10,"(A)")"SCIFOR VERSION (GIT): "//trim(adjustl(trim(sf_version)))
    write(10,"(A)")"CODE VERSION (GIT): "//trim(adjustl(trim(revision)))
    write(10,"(A,i2,1x,a,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3)")&
         "Timestamp: +",day,trim(month(mese)),year, h,':',m,':',s,'.',ms
    write(10,*)""
    close(10)
  end subroutine dmftt_version

END MODULE DMFT_TOOLS_VERSION
