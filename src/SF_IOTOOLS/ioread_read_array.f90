subroutine data_readA1_R(pname,Y1)
  integer               :: i,Np
  character(len=*)      :: pname
  real(8),dimension(:)  :: Y1
  !
  call file_gunzip(reg(pname))
  include "ioread_control.f90"
  !
  open(free_unit(unit),file=reg(pname))
  !
  Np=size(Y1)
  !
  do i=1,Np
     read(unit,*)Y1(i)
  enddo
  close(unit)
  call file_gzip(reg(pname))
end subroutine data_readA1_R

subroutine data_readA1_C(pname,Y1)
  integer                :: i,Np
  character(len=*)       :: pname
  complex(8),dimension(:):: Y1
  !
  call file_gunzip(reg(pname))
  include "ioread_control.f90"
  !
  open(free_unit(unit),file=reg(pname))
  !
  Np=size(Y1)
  !
  do i=1,Np
     read(unit,*)Y1(i)
  enddo
  close(unit)
  call file_gzip(reg(pname))
end subroutine data_readA1_C



!------------------------------------------------------------------!
!------------------------------------------------------------------!
!------------------------------------------------------------------!



subroutine data_readA2_R(pname,Y1)
  integer                       :: i,j,Ny1,Ny2
  character(len=*)              :: pname
  real(8),dimension(:,:)        :: Y1
  !
  call file_gunzip(reg(pname))
  include "ioread_control.f90"
  !
  open(free_unit(unit),file=reg(pname))
  !
  Ny1=size(Y1,1) ; Ny2=size(Y1,2)
  !
  do i=1,Ny1
     do j=1,Ny2
        read(unit,*)Y1(i,j)
     enddo
  enddo
  close(unit)
  call file_gzip(reg(pname))
end subroutine data_readA2_R

subroutine data_readA2_C(pname,Y1)
  integer                       :: i,j,Ny1,Ny2
  character(len=*)              :: pname
  complex(8),dimension(:,:)     :: Y1
  !
  call file_gunzip(reg(pname))
  include "ioread_control.f90"
  !
  open(free_unit(unit),file=reg(pname))
  !
  Ny1=size(Y1,1) ; Ny2=size(Y1,2)
  !
  do i=1,Ny1
     do j=1,Ny2
        read(unit,*)Y1(i,j)
     enddo
  enddo
  close(unit)
  call file_gzip(reg(pname))
end subroutine data_readA2_C

!------------------------------------------------------------------!
!------------------------------------------------------------------!
!------------------------------------------------------------------!

subroutine data_readA3_R(pname,Y1)
  integer                  :: Ny1,Ny2,Ny3
  integer                  :: i1,i2,i3
  character(len=*)         :: pname
  real(8),dimension(:,:,:) :: Y1
  !
  call file_gunzip(reg(pname))
  include "ioread_control.f90"
  !
  open(free_unit(unit),file=reg(pname))
  !
  Ny1=size(Y1,1)
  Ny2=size(Y1,2)
  Ny3=size(Y1,3)
  !
  do i1=1,Ny1
     do i2=1,Ny2
        do i3=1,Ny3
           read(unit,*)Y1(i1,i2,i3)
        enddo
     enddo
  enddo
  close(unit)
  call file_gzip(reg(pname))
end subroutine data_readA3_R

subroutine data_readA3_C(pname,Y1)
  integer                     :: Ny1,Ny2,Ny3
  integer                     :: i1,i2,i3
  character(len=*)            :: pname
  complex(8),dimension(:,:,:) :: Y1
  !
  call file_gunzip(reg(pname))
  include "ioread_control.f90"
  !
  open(free_unit(unit),file=reg(pname))
  !
  Ny1=size(Y1,1)
  Ny2=size(Y1,2)
  Ny3=size(Y1,3)
  !
  do i1=1,Ny1
     do i2=1,Ny2
        do i3=1,Ny3
           read(unit,*)Y1(i1,i2,i3)
        enddo
     enddo
  enddo
  close(unit)
  call file_gzip(reg(pname))
end subroutine data_readA3_C


!------------------------------------------------------------------!
!------------------------------------------------------------------!
!------------------------------------------------------------------!





subroutine data_readA4_R(pname,Y1)
  integer                                :: Ny1,Ny2,Ny3,Ny4
  integer                                :: i1,i2,i3,i4
  character(len=*)                       :: pname
  real(8),dimension(:,:,:,:)             :: Y1
  !
  call file_gunzip(reg(pname))
  include "ioread_control.f90"
  !
  open(free_unit(unit),file=reg(pname))
  !
  Ny1=size(Y1,1)
  Ny2=size(Y1,2)
  Ny3=size(Y1,3)
  Ny4=size(Y1,4)
  do i1=1,Ny1
     do i2=1,Ny2
        do i3=1,Ny3
           do i4=1,Ny4
              read(unit,*)Y1(i1,i2,i3,i4)
           enddo
        enddo
     enddo
  enddo
  close(unit)
  call file_gzip(reg(pname))
end subroutine data_readA4_R

subroutine data_readA4_C(pname,Y1)
  integer                       :: Ny1,Ny2,Ny3,Ny4
  integer                       :: i1,i2,i3,i4
  character(len=*)              :: pname
  complex(8),dimension(:,:,:,:) :: Y1
  !
  call file_gunzip(reg(pname))
  include "ioread_control.f90"
  !
  open(free_unit(unit),file=reg(pname))
  !
  Ny1=size(Y1,1)
  Ny2=size(Y1,2)
  Ny3=size(Y1,3)
  Ny4=size(Y1,4)
  do i1=1,Ny1
     do i2=1,Ny2
        do i3=1,Ny3
           do i4=1,Ny4
              read(unit,*)Y1(i1,i2,i3,i4)              
           enddo
        enddo
     enddo
  enddo
  close(unit)
  call file_gzip(reg(pname))
end subroutine data_readA4_C



!----------------------------
!----------------------------
!----------------------------



subroutine data_readA5_R(pname,Y1)
  integer                      :: Ny1,Ny2,Ny3,Ny4,Ny5
  integer                      :: i1,i2,i3,i4,i5
  character(len=*)             :: pname
  real(8),dimension(:,:,:,:,:) :: Y1
  !
  call file_gunzip(reg(pname))
  include "ioread_control.f90"
  !
  open(free_unit(unit),file=reg(pname))
  !
  Ny1=size(Y1,1)
  Ny2=size(Y1,2)
  Ny3=size(Y1,3)
  Ny4=size(Y1,4)
  Ny5=size(Y1,5)
  do i1=1,Ny1
     do i2=1,Ny2
        do i3=1,Ny3
           do i4=1,Ny4
              do i5=1,Ny5
                 read(unit,*)Y1(i1,i2,i3,i4,i5)
              enddo
           enddo
        enddo
     enddo
  enddo
  close(unit)
  call file_gzip(reg(pname))
end subroutine data_readA5_R

subroutine data_readA5_C(pname,Y1)
  integer                         :: Ny1,Ny2,Ny3,Ny4,Ny5
  integer                         :: i1,i2,i3,i4,i5
  character(len=*)                :: pname
  complex(8),dimension(:,:,:,:,:) :: Y1
  !
  call file_gunzip(reg(pname))
  include "ioread_control.f90"
  !
  open(free_unit(unit),file=reg(pname))
  !
  Ny1=size(Y1,1)
  Ny2=size(Y1,2)
  Ny3=size(Y1,3)
  Ny4=size(Y1,4)
  Ny5=size(Y1,5)
  do i1=1,Ny1
     do i2=1,Ny2
        do i3=1,Ny3
           do i4=1,Ny4
              do i5=1,Ny5
                 read(unit,*)Y1(i1,i2,i3,i4,i5)
              enddo
           enddo
        enddo
     enddo
  enddo
  close(unit)
  call file_gzip(reg(pname))
end subroutine data_readA5_C






!----------------------------
!----------------------------
!----------------------------


subroutine data_readA6_R(pname,Y1)
  integer                                :: Ny1,Ny2,Ny3,Ny4,Ny5,Ny6
  integer                                :: i1,i2,i3,i4,i5,i6
  character(len=*)                       :: pname
  real(8),dimension(:,:,:,:,:,:)         :: Y1
  !
  call file_gunzip(reg(pname))
  include "ioread_control.f90"
  !
  open(free_unit(unit),file=reg(pname))
  !
  Ny1=size(Y1,1)
  Ny2=size(Y1,2)
  Ny3=size(Y1,3)
  Ny4=size(Y1,4)
  Ny5=size(Y1,5)
  Ny6=size(Y1,6)
  !
  do i1=1,Ny1
     do i2=1,Ny2
        do i3=1,Ny3
           do i4=1,Ny4
              do i5=1,Ny5
                 do i6=1,Ny6
                    read(unit,*)Y1(i1,i2,i3,i4,i5,i6)
                 enddo
              enddo
           enddo
        enddo
     enddo
  enddo
  close(unit)
  call file_gzip(reg(pname))
end subroutine data_readA6_R

subroutine data_readA6_C(pname,Y1)
  integer                           :: Ny1,Ny2,Ny3,Ny4,Ny5,Ny6
  integer                           :: i1,i2,i3,i4,i5,i6
  character(len=*)                  :: pname
  complex(8),dimension(:,:,:,:,:,:) :: Y1
  !
  call file_gunzip(reg(pname))
  include "ioread_control.f90"
  !
  open(free_unit(unit),file=reg(pname))
  !
  Ny1=size(Y1,1)
  Ny2=size(Y1,2)
  Ny3=size(Y1,3)
  Ny4=size(Y1,4)
  Ny5=size(Y1,5)
  Ny6=size(Y1,6)
  !
  do i1=1,Ny1
     do i2=1,Ny2
        do i3=1,Ny3
           do i4=1,Ny4
              do i5=1,Ny5
                 do i6=1,Ny6
                    read(unit,*)Y1(i1,i2,i3,i4,i5,i6)
                 enddo
              enddo
           enddo
        enddo
     enddo
  enddo
  close(unit)
  call file_gzip(reg(pname))
end subroutine data_readA6_C






!----------------------------
!----------------------------
!----------------------------



subroutine data_readA7_R(pname,Y1)
  integer                          :: Ny1,Ny2,Ny3,Ny4,Ny5,Ny6,Ny7
  integer                          :: i1,i2,i3,i4,i5,i6,i7
  character(len=*)                 :: pname
  real(8),dimension(:,:,:,:,:,:,:) :: Y1
  !
  call file_gunzip(reg(pname))
  include "ioread_control.f90"
  !
  open(free_unit(unit),file=reg(pname))
  !
  Ny1=size(Y1,1)
  Ny2=size(Y1,2)
  Ny3=size(Y1,3)
  Ny4=size(Y1,4)
  Ny5=size(Y1,5)
  Ny6=size(Y1,6)
  Ny7=size(Y1,7)
  !
  do i1=1,Ny1
     do i2=1,Ny2
        do i3=1,Ny3
           do i4=1,Ny4
              do i5=1,Ny5
                 do i6=1,Ny6
                    do i7=1,Ny7
                       read(unit,*)Y1(i1,i2,i3,i4,i5,i6,i7)
                    enddo
                 enddo
              enddo
           enddo
        enddo
     enddo
  enddo
  close(unit)
  call file_gzip(reg(pname))
end subroutine data_readA7_R

subroutine data_readA7_C(pname,Y1)
  integer                             :: Ny1,Ny2,Ny3,Ny4,Ny5,Ny6,Ny7
  integer                             :: i1,i2,i3,i4,i5,i6,i7
  character(len=*)                    :: pname
  complex(8),dimension(:,:,:,:,:,:,:) :: Y1
  !
  call file_gunzip(reg(pname))
  include "ioread_control.f90"
  !
  open(free_unit(unit),file=reg(pname))
  !
  Ny1=size(Y1,1)
  Ny2=size(Y1,2)
  Ny3=size(Y1,3)
  Ny4=size(Y1,4)
  Ny5=size(Y1,5)
  Ny6=size(Y1,6)
  Ny7=size(Y1,7)
  !
  do i1=1,Ny1
     do i2=1,Ny2
        do i3=1,Ny3
           do i4=1,Ny4
              do i5=1,Ny5
                 do i6=1,Ny6
                    do i7=1,Ny7
                       read(unit,*)Y1(i1,i2,i3,i4,i5,i6,i7)
                    enddo
                 enddo
              enddo
           enddo
        enddo
     enddo
  enddo
  close(unit)
  call file_gzip(reg(pname))
end subroutine data_readA7_C
