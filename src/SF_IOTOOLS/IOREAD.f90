module IOREAD

   USE IOFILE
   implicit none
   private

   integer            :: unit
   character(len=128) :: fmt
   logical            :: control

   interface sread
      ! SPLOT arrays (1--7)
      module procedure :: sreadA1_RR
      module procedure :: sreadA1_RC
      module procedure :: sreadA2_RR
      module procedure :: sreadA2_RC
      module procedure :: sreadA3_RR
      module procedure :: sreadA3_RC
      module procedure :: sreadA4_RR
      module procedure :: sreadA4_RC
      module procedure :: sreadA5_RR
      module procedure :: sreadA5_RC
      module procedure :: sreadA6_RR
      module procedure :: sreadA6_RC
      module procedure :: sreadA7_RR
      module procedure :: sreadA7_RC
   end interface sread


   interface read_array
      ! READ_ARRAY arrays (1--7)
      module procedure :: data_readA1_R
      module procedure :: data_readA1_C
      module procedure :: data_readA2_R
      module procedure :: data_readA2_C
      module procedure :: data_readA3_R
      module procedure :: data_readA3_C
      module procedure :: data_readA4_R
      module procedure :: data_readA4_C
      module procedure :: data_readA5_R
      module procedure :: data_readA5_C
      module procedure :: data_readA6_R
      module procedure :: data_readA6_C
      module procedure :: data_readA7_R
      module procedure :: data_readA7_C
   end interface read_array


   public :: sread
   public :: read_array

contains


   ! SPLOT arrays (1--7)
   subroutine sreadA1_RR(pname,X,Y1)
      integer                             :: i,Np
      character(len=*)                    :: pname
      real(8),dimension(:)                :: X
      real(8),dimension(size(X))          :: Y1
      call ioread_control(pname)
      open(free_unit(unit),file=reg(pname))
      Np=size(X)
      do i=1,Np
         read(unit,*)X(i),Y1(i)
      enddo
      close(unit)
   end subroutine sreadA1_RR

   subroutine sreadA1_RC(pname,X,Y1)
      integer                       :: i,Np
      character(len=*)              :: pname
      real(8),dimension(:)          :: X
      complex(8),dimension(size(X)) :: Y1
      real(8),dimension(size(X))    :: reY,imY
      call ioread_control(pname)
      open(free_unit(unit),file=reg(pname))
      Np=size(X)
      do i=1,Np
         read(unit,*)X(i),imY(i),reY(i)
      enddo
      Y1=dcmplx(reY,imY)
      close(unit)
   end subroutine sreadA1_RC









   subroutine sreadA2_RR(pname,X,Y1)
      integer                       :: i,j,Ny1,Ny2
      character(len=*)              :: pname
      real(8),dimension(:,:)        :: Y1
      real(8),dimension(size(Y1,2)) :: X
      call ioread_control(pname)
      open(free_unit(unit),file=reg(pname))
      !
      Ny1=size(Y1,1)
      Ny2=size(Y1,2)
      !
      do i=1,Ny1
         do j=1,Ny2
            read(unit,*)X(j),Y1(i,j)
         enddo
      enddo
      !
      close(unit)
   end subroutine sreadA2_RR

   subroutine sreadA2_RC(pname,X,Y1)
      integer                                  :: i,j,Ny1,Ny2
      character(len=*)                         :: pname
      complex(8),dimension(:,:)                :: Y1
      real(8),dimension(size(Y1,2))            :: X
      real(8),dimension(size(Y1,1),size(Y1,2)) :: reY,imY
      call ioread_control(pname)
      open(free_unit(unit),file=reg(pname))
      !
      Ny1=size(Y1,1)
      Ny2=size(Y1,2)
      do i=1,Ny1
         do j=1,Ny2
            read(unit,*)X(j),imY(i,j),reY(i,j)
         enddo
      enddo
      close(unit)
      Y1=dcmplx(reY,imY)
   end subroutine sreadA2_RC









   subroutine sreadA3_RR(pname,X,Y1)
      integer                       :: i,j,k,Ny1,Ny2,Ny3
      character(len=*)              :: pname
      real(8),dimension(:,:,:)      :: Y1
      real(8),dimension(size(Y1,3)) :: X
      call ioread_control(pname)
      open(free_unit(unit),file=reg(pname))
      !
      Ny1=size(Y1,1)
      Ny2=size(Y1,2)
      Ny3=size(Y1,3)
      !
      do i=1,Ny1
         do j=1,Ny2
            do k=1,Ny3
               read(unit,*)X(k),Y1(i,j,k)
            enddo
         enddo
      enddo
      close(unit)
   end subroutine sreadA3_RR

   subroutine sreadA3_RC(pname,X,Y1)
      integer                                             :: i,j,k,Ny1,Ny2,Ny3
      character(len=*)                                    :: pname
      complex(8),dimension(:,:,:)                         :: Y1
      real(8),dimension(size(Y1,3))                       :: X
      real(8),dimension(size(Y1,1),size(Y1,2),size(Y1,3)) :: reY,imY
      call ioread_control(pname)
      open(free_unit(unit),file=reg(pname))
      !
      Ny1=size(Y1,1)
      Ny2=size(Y1,2)
      Ny3=size(Y1,3)
      !
      do i=1,Ny1
         do j=1,Ny2
            do k=1,Ny3
               read(unit,*)X(k),imY(i,j,k),reY(i,j,k)
            enddo
         enddo
      enddo
      close(unit)
      Y1=dcmplx(reY,imY)
   end subroutine sreadA3_RC







   subroutine sreadA4_RR(pname,X,Y1)
      integer                       :: Ny1,Ny2,Ny3,Ny4
      integer                       :: i1,i2,i3,i4
      character(len=*)              :: pname
      real(8),dimension(:,:,:,:)    :: Y1
      real(8),dimension(size(Y1,4)) :: X
      call ioread_control(pname)
      open(free_unit(unit),file=reg(pname))
      !
      Ny1=size(Y1,1)
      Ny2=size(Y1,2)
      Ny3=size(Y1,3)
      Ny4=size(Y1,4)
      !
      do i1=1,Ny1
         do i2=1,Ny2
            do i3=1,Ny3
               do i4=1,Ny4
                  read(unit,*)X(i4),Y1(i1,i2,i3,i4)
               enddo
            enddo
         enddo
      enddo
      close(unit)
   end subroutine sreadA4_RR

   subroutine sreadA4_RC(pname,X,Y1)
      integer                       :: Ny1,Ny2,Ny3,Ny4
      integer                       :: i1,i2,i3,i4
      character(len=*)              :: pname
      complex(8),dimension(:,:,:,:) :: Y1
      real(8),dimension(size(Y1,4)) :: X
      real(8),dimension(&
         size(Y1,1),&
         size(Y1,2),&
         size(Y1,3),&
         size(Y1,4))              :: reY,imY
      !
      call ioread_control(pname)
      open(free_unit(unit),file=reg(pname))
      !
      !
      Ny1=size(Y1,1)
      Ny2=size(Y1,2)
      Ny3=size(Y1,3)
      Ny4=size(Y1,4)
      !
      do i1=1,Ny1
         do i2=1,Ny2
            do i3=1,Ny3
               do i4=1,Ny4
                  read(unit,*)X(i4),imY(i1,i2,i3,i4),reY(i1,i2,i3,i4)
               enddo
            enddo
         enddo
      enddo
      close(unit)
      Y1=dcmplx(reY,imY)
   end subroutine sreadA4_RC








   subroutine sreadA5_RR(pname,X,Y1)
      integer                         :: Ny1,Ny2,Ny3,Ny4,Ny5
      integer                         :: i1,i2,i3,i4,i5
      character(len=*)                :: pname
      real(8),dimension(:,:,:,:,:)    :: Y1
      real(8),dimension(size(Y1,5))   :: X
      !
      call ioread_control(pname)
      open(free_unit(unit),file=reg(pname))
      !
      Ny1=size(Y1,1)
      Ny2=size(Y1,2)
      Ny3=size(Y1,3)
      Ny4=size(Y1,4)
      Ny5=size(Y1,5)
      !
      do i1=1,Ny1
         do i2=1,Ny2
            do i3=1,Ny3
               do i4=1,Ny4
                  do i5=1,Ny5
                     read(unit,*)X(i5),Y1(i1,i2,i3,i4,i5)
                  enddo
               enddo
            enddo
         enddo
      enddo
      close(unit)
   end subroutine sreadA5_RR

   subroutine sreadA5_RC(pname,X,Y1)
      integer                         :: Ny1,Ny2,Ny3,Ny4,Ny5
      integer                         :: i1,i2,i3,i4,i5
      character(len=*)                :: pname
      complex(8),dimension(:,:,:,:,:) :: Y1
      real(8),dimension(size(Y1,5))   :: X
      real(8),dimension(&
         size(Y1,1),&
         size(Y1,2),&
         size(Y1,3),&
         size(Y1,4),&
         size(Y1,5))                :: reY,imY
      !
      call ioread_control(pname)
      open(free_unit(unit),file=reg(pname))
      !
      !
      Ny1=size(Y1,1)
      Ny2=size(Y1,2)
      Ny3=size(Y1,3)
      Ny4=size(Y1,4)
      Ny5=size(Y1,5)
      !
      do i1=1,Ny1
         do i2=1,Ny2
            do i3=1,Ny3
               do i4=1,Ny4
                  do i5=1,Ny5
                     read(unit,*)X(i5),imY(i1,i2,i3,i4,i5),reY(i1,i2,i3,i4,i5)
                  enddo
               enddo
            enddo
         enddo
      enddo
      close(unit)
      Y1=dcmplx(reY,imY)
   end subroutine sreadA5_RC







   subroutine sreadA6_RR(pname,X,Y1)
      integer                           :: Ny1,Ny2,Ny3,Ny4,Ny5,Ny6
      integer                           :: i1,i2,i3,i4,i5,i6
      character(len=*)                  :: pname
      real(8),dimension(:,:,:,:,:,:)    :: Y1
      real(8),dimension(size(Y1,6))     :: X
      !
      call ioread_control(pname)
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
                        read(unit,*)X(i6),Y1(i1,i2,i3,i4,i5,i6)
                     enddo
                  enddo
               enddo
            enddo
         enddo
      enddo
      close(unit)
   end subroutine sreadA6_RR

   subroutine sreadA6_RC(pname,X,Y1)
      integer                           :: Ny1,Ny2,Ny3,Ny4,Ny5,Ny6
      integer                           :: i1,i2,i3,i4,i5,i6
      character(len=*)                  :: pname
      complex(8),dimension(:,:,:,:,:,:) :: Y1
      real(8),dimension(size(Y1,6))     :: X
      real(8),dimension(&
         size(Y1,1),&
         size(Y1,2),&
         size(Y1,3),&
         size(Y1,4),&
         size(Y1,5),&
         size(Y1,6))                  :: reY,imY
      !
      call ioread_control(pname)
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
                        read(unit,*)X(i6),imY(i1,i2,i3,i4,i5,i6),reY(i1,i2,i3,i4,i5,i6)
                     enddo
                  enddo
               enddo
            enddo
         enddo
      enddo
      close(unit)
      Y1=dcmplx(reY,imY)
   end subroutine sreadA6_RC








   subroutine sreadA7_RR(pname,X,Y1)
      integer                             :: Ny1,Ny2,Ny3,Ny4,Ny5,Ny6,Ny7
      integer                             :: i1,i2,i3,i4,i5,i6,i7
      character(len=*)                    :: pname
      real(8),dimension(:,:,:,:,:,:,:)    :: Y1
      real(8),dimension(size(Y1,7))       :: X
      !
      call ioread_control(pname)
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
                           read(unit,*)X(i7),Y1(i1,i2,i3,i4,i5,i6,i7)
                        enddo
                     enddo
                  enddo
               enddo
            enddo
         enddo
      enddo
      close(unit)
   end subroutine sreadA7_RR

   subroutine sreadA7_RC(pname,X,Y1)
      integer                             :: Ny1,Ny2,Ny3,Ny4,Ny5,Ny6,Ny7
      integer                             :: i1,i2,i3,i4,i5,i6,i7
      character(len=*)                    :: pname
      complex(8),dimension(:,:,:,:,:,:,:) :: Y1
      real(8),dimension(size(Y1,7))       :: X
      real(8),dimension(&
         size(Y1,1),&
         size(Y1,2),&
         size(Y1,3),&
         size(Y1,4),&
         size(Y1,5),&
         size(Y1,6),&
         size(Y1,7))                    :: reY,imY
      !
      call ioread_control(pname)
      open(free_unit(unit),file=reg(pname))
      !
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
                           read(unit,*)X(i7),imY(i1,i2,i3,i4,i5,i6,i7),reY(i1,i2,i3,i4,i5,i6,i7)
                        enddo
                     enddo
                  enddo
               enddo
            enddo
         enddo
      enddo
      close(unit)
      Y1=dcmplx(reY,imY)
   end subroutine sreadA7_RC


   ! READ_ARRAY arrays (1--7)
   subroutine data_readA1_R(pname,Y1)
      integer               :: i,Np
      character(len=*)      :: pname
      real(8),dimension(:)  :: Y1
      !
      call file_gunzip(reg(pname))
      call ioread_control(pname)
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
      call ioread_control(pname)
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









   subroutine data_readA2_R(pname,Y1,order,wspace)
      integer                   :: i,j,Ny1,Ny2
      character(len=*)          :: pname
      real(8),dimension(:,:)    :: Y1
      character(len=*),optional :: order
      logical,optional          :: wspace
      character(len=1)          :: order_
      logical                   :: wspace_
      order_ = "R"   ; if(present(order))order_=trim(order(1:1))
      wspace_= .true.; if(present(wspace))wspace_=wspace
      !
      call file_gunzip(reg(pname))
      call ioread_control(pname)
      !
      open(free_unit(unit),file=reg(pname))
      !
      Ny1=size(Y1,1) ; Ny2=size(Y1,2)
      !
      select case(order_)
       case default
         stop "read_array: order != Row-major, Col-major"
       case ("R")
         do i=1,Ny1
            do j=1,Ny2
               read(unit,*)Y1(i,j)
            enddo
         enddo
       case ("C")
         do j=1,Ny2
            do i=1,Ny1
               read(unit,*)Y1(i,j)
            enddo
         enddo
      end select
      close(unit)
      call file_gzip(reg(pname))
   end subroutine data_readA2_R

   subroutine data_readA2_C(pname,Y1,order,wspace)
      integer                   :: i,j,Ny1,Ny2
      character(len=*)          :: pname
      complex(8),dimension(:,:) :: Y1
      character(len=*),optional :: order
      logical,optional          :: wspace
      character(len=1)          :: order_
      logical                   :: wspace_
      order_ = "R"   ; if(present(order))order_=trim(order(1:1))
      wspace_= .true.; if(present(wspace))wspace_=wspace
      !
      call file_gunzip(reg(pname))
      call ioread_control(pname)
      !
      open(free_unit(unit),file=reg(pname))
      !
      Ny1=size(Y1,1) ; Ny2=size(Y1,2)
      !
      select case(order_)
       case default
         stop "read_array: order != Row-major, Col-major"
       case ("R")
         do i=1,Ny1
            do j=1,Ny2
               read(unit,*)Y1(i,j)
            enddo
         enddo
       case ("C")
         do j=1,Ny2
            do i=1,Ny1
               read(unit,*)Y1(i,j)
            enddo
         enddo
      end select
      close(unit)
      call file_gzip(reg(pname))
   end subroutine data_readA2_C





   subroutine data_readA3_R(pname,Y1,order,wspace)
      integer                   :: Ny1,Ny2,Ny3
      integer                   :: i1,i2,i3
      character(len=*)          :: pname
      real(8),dimension(:,:,:)  :: Y1
      character(len=*),optional :: order
      logical,optional          :: wspace
      character(len=1)          :: order_
      logical                   :: wspace_
      order_ = "R"   ; if(present(order))order_=trim(order(1:1))
      wspace_= .true.; if(present(wspace))wspace_=wspace
      !
      call file_gunzip(reg(pname))
      call ioread_control(pname)
      !
      open(free_unit(unit),file=reg(pname))
      !
      Ny1=size(Y1,1)
      Ny2=size(Y1,2)
      Ny3=size(Y1,3)
      !
      select case(order_)
       case default
         stop "read_array: order != Row-major, Col-major"
       case ("R")
         do i1=1,Ny1
            do i2=1,Ny2
               do i3=1,Ny3
                  read(unit,*)Y1(i1,i2,i3)
               enddo
            enddo
         enddo
       case ("C")
         do i3=1,Ny3
            do i2=1,Ny2
               do i1=1,Ny1
                  read(unit,*)Y1(i1,i2,i3)
               enddo
            enddo
         enddo
      end select
      close(unit)
      call file_gzip(reg(pname))
   end subroutine data_readA3_R

   subroutine data_readA3_C(pname,Y1,order,wspace)
      integer                     :: Ny1,Ny2,Ny3
      integer                     :: i1,i2,i3
      character(len=*)            :: pname
      complex(8),dimension(:,:,:) :: Y1
      character(len=*),optional   :: order
      logical,optional            :: wspace
      character(len=1)            :: order_
      logical                     :: wspace_
      order_ = "R"   ; if(present(order))order_=trim(order(1:1))
      wspace_= .true.; if(present(wspace))wspace_=wspace
      !
      call file_gunzip(reg(pname))
      call ioread_control(pname)
      !
      open(free_unit(unit),file=reg(pname))
      !
      Ny1=size(Y1,1)
      Ny2=size(Y1,2)
      Ny3=size(Y1,3)
      !
      select case(order_)
       case default
         stop "read_array: order != Row-major, Col-major"
       case ("R")
         do i1=1,Ny1
            do i2=1,Ny2
               do i3=1,Ny3
                  read(unit,*)Y1(i1,i2,i3)
               enddo
            enddo
         enddo
       case ("C")
         do i3=1,Ny3
            do i2=1,Ny2
               do i1=1,Ny1
                  read(unit,*)Y1(i1,i2,i3)
               enddo
            enddo
         enddo
      end select
      close(unit)
      call file_gzip(reg(pname))
   end subroutine data_readA3_C










   subroutine data_readA4_R(pname,Y1,order,wspace)
      integer                    :: Ny1,Ny2,Ny3,Ny4
      integer                    :: i1,i2,i3,i4
      character(len=*)           :: pname
      real(8),dimension(:,:,:,:) :: Y1
      character(len=*),optional  :: order
      logical,optional           :: wspace
      character(len=1)           :: order_
      logical                    :: wspace_
      order_ = "R"   ; if(present(order))order_=trim(order(1:1))
      wspace_= .true.; if(present(wspace))wspace_=wspace
      !
      call file_gunzip(reg(pname))
      call ioread_control(pname)
      !
      open(free_unit(unit),file=reg(pname))
      !
      Ny1=size(Y1,1)
      Ny2=size(Y1,2)
      Ny3=size(Y1,3)
      Ny4=size(Y1,4)
      !
      select case(order_)
       case default
         stop "read_array: order != Row-major, Col-major"
       case ("R")
         do i1=1,Ny1
            do i2=1,Ny2
               do i3=1,Ny3
                  do i4=1,Ny4
                     read(unit,*)Y1(i1,i2,i3,i4)
                  enddo
               enddo
            enddo
         enddo
       case ("C")
         do i4=1,Ny4
            do i3=1,Ny3
               do i2=1,Ny2
                  do i1=1,Ny1
                     read(unit,*)Y1(i1,i2,i3,i4)
                  enddo
               enddo
            enddo
         enddo
      end select
      close(unit)
      call file_gzip(reg(pname))
   end subroutine data_readA4_R

   subroutine data_readA4_C(pname,Y1,order,wspace)
      integer                       :: Ny1,Ny2,Ny3,Ny4
      integer                       :: i1,i2,i3,i4
      character(len=*)              :: pname
      complex(8),dimension(:,:,:,:) :: Y1
      character(len=*),optional     :: order
      logical,optional              :: wspace
      character(len=1)              :: order_
      logical                       :: wspace_
      order_ = "R"   ; if(present(order))order_=trim(order(1:1))
      wspace_= .true.; if(present(wspace))wspace_=wspace
      !
      call file_gunzip(reg(pname))
      call ioread_control(pname)
      !
      open(free_unit(unit),file=reg(pname))
      !
      Ny1=size(Y1,1)
      Ny2=size(Y1,2)
      Ny3=size(Y1,3)
      Ny4=size(Y1,4)
      !
      select case(order_)
       case default
         stop "read_array: order != Row-major, Col-major"
       case ("R")
         do i1=1,Ny1
            do i2=1,Ny2
               do i3=1,Ny3
                  do i4=1,Ny4
                     read(unit,*)Y1(i1,i2,i3,i4)
                  enddo
               enddo
            enddo
         enddo
       case ("C")
         do i4=1,Ny4
            do i3=1,Ny3
               do i2=1,Ny2
                  do i1=1,Ny1
                     read(unit,*)Y1(i1,i2,i3,i4)
                  enddo
               enddo
            enddo
         enddo
      end select
      !
      close(unit)
      call file_gzip(reg(pname))
   end subroutine data_readA4_C









   subroutine data_readA5_R(pname,Y1,order,wspace)
      integer                      :: Ny1,Ny2,Ny3,Ny4,Ny5
      integer                      :: i1,i2,i3,i4,i5
      character(len=*)             :: pname
      real(8),dimension(:,:,:,:,:) :: Y1
      character(len=*),optional    :: order
      logical,optional             :: wspace
      character(len=1)             :: order_
      logical                      :: wspace_
      order_ = "R"   ; if(present(order))order_=trim(order(1:1))
      wspace_= .true.; if(present(wspace))wspace_=wspace
      !
      call file_gunzip(reg(pname))
      call ioread_control(pname)
      !
      open(free_unit(unit),file=reg(pname))
      !
      Ny1=size(Y1,1)
      Ny2=size(Y1,2)
      Ny3=size(Y1,3)
      Ny4=size(Y1,4)
      Ny5=size(Y1,5)
      !
      select case(order_)
       case default
         stop "read_array: order != Row-major, Col-major"
       case ("R")
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
       case ("C")
         do i5=1,Ny5
            do i4=1,Ny4
               do i3=1,Ny3
                  do i2=1,Ny2
                     do i1=1,Ny1
                        read(unit,*)Y1(i1,i2,i3,i4,i5)
                     enddo
                  enddo
               enddo
            enddo
         enddo
      end select
      !
      close(unit)
      call file_gzip(reg(pname))
   end subroutine data_readA5_R

   subroutine data_readA5_C(pname,Y1,order,wspace)
      integer                         :: Ny1,Ny2,Ny3,Ny4,Ny5
      integer                         :: i1,i2,i3,i4,i5
      character(len=*)                :: pname
      complex(8),dimension(:,:,:,:,:) :: Y1
      character(len=*),optional       :: order
      logical,optional                :: wspace
      character(len=1)                :: order_
      logical                         :: wspace_
      order_ = "R"   ; if(present(order))order_=trim(order(1:1))
      wspace_= .true.; if(present(wspace))wspace_=wspace
      !
      call file_gunzip(reg(pname))
      call ioread_control(pname)
      !
      open(free_unit(unit),file=reg(pname))
      !
      Ny1=size(Y1,1)
      Ny2=size(Y1,2)
      Ny3=size(Y1,3)
      Ny4=size(Y1,4)
      Ny5=size(Y1,5)
      !
      select case(order_)
       case default
         stop "read_array: order != Row-major, Col-major"
       case ("R")
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
       case ("C")
         do i5=1,Ny5
            do i4=1,Ny4
               do i3=1,Ny3
                  do i2=1,Ny2
                     do i1=1,Ny1
                        read(unit,*)Y1(i1,i2,i3,i4,i5)
                     enddo
                  enddo
               enddo
            enddo
         enddo
      end select
      !
      close(unit)
      call file_gzip(reg(pname))
   end subroutine data_readA5_C











   subroutine data_readA6_R(pname,Y1,order,wspace)
      integer                        :: Ny1,Ny2,Ny3,Ny4,Ny5,Ny6
      integer                        :: i1,i2,i3,i4,i5,i6
      character(len=*)               :: pname
      real(8),dimension(:,:,:,:,:,:) :: Y1
      character(len=*),optional      :: order
      logical,optional               :: wspace
      character(len=1)               :: order_
      logical                        :: wspace_
      order_ = "R"   ; if(present(order))order_=trim(order(1:1))
      wspace_= .true.; if(present(wspace))wspace_=wspace
      !
      call file_gunzip(reg(pname))
      call ioread_control(pname)
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
      select case(order_)
       case default
         stop "read_array: order != Row-major, Col-major"
       case ("R")
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
       case ("C")
         do i6=1,Ny6
            do i5=1,Ny5
               do i4=1,Ny4
                  do i3=1,Ny3
                     do i2=1,Ny2
                        do i1=1,Ny1
                           read(unit,*)Y1(i1,i2,i3,i4,i5,i6)
                        enddo
                     enddo
                  enddo
               enddo
            enddo
         enddo
      end select
      close(unit)
      call file_gzip(reg(pname))
   end subroutine data_readA6_R

   subroutine data_readA6_C(pname,Y1,order,wspace)
      integer                           :: Ny1,Ny2,Ny3,Ny4,Ny5,Ny6
      integer                           :: i1,i2,i3,i4,i5,i6
      character(len=*)                  :: pname
      complex(8),dimension(:,:,:,:,:,:) :: Y1
      character(len=*),optional         :: order
      logical,optional                  :: wspace
      character(len=1)                  :: order_
      logical                           :: wspace_
      order_ = "R"   ; if(present(order))order_=trim(order(1:1))
      wspace_= .true.; if(present(wspace))wspace_=wspace
      !
      call file_gunzip(reg(pname))
      call ioread_control(pname)
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
      select case(order_)
       case default
         stop "read_array: order != Row-major, Col-major"
       case ("R")
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
       case ("C")
         do i6=1,Ny6
            do i5=1,Ny5
               do i4=1,Ny4
                  do i3=1,Ny3
                     do i2=1,Ny2
                        do i1=1,Ny1
                           read(unit,*)Y1(i1,i2,i3,i4,i5,i6)
                        enddo
                     enddo
                  enddo
               enddo
            enddo
         enddo
      end select
      close(unit)
      call file_gzip(reg(pname))
   end subroutine data_readA6_C












   subroutine data_readA7_R(pname,Y1,order,wspace)
      integer                          :: Ny1,Ny2,Ny3,Ny4,Ny5,Ny6,Ny7
      integer                          :: i1,i2,i3,i4,i5,i6,i7
      character(len=*)                 :: pname
      real(8),dimension(:,:,:,:,:,:,:) :: Y1
      character(len=*),optional        :: order
      logical,optional                 :: wspace
      character(len=1)                 :: order_
      logical                          :: wspace_
      order_ = "R"   ; if(present(order))order_=trim(order(1:1))
      wspace_= .true.; if(present(wspace))wspace_=wspace
      !
      call file_gunzip(reg(pname))
      call ioread_control(pname)
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
      select case(order_)
       case default
         stop "read_array: order != Row-major, Col-major"
       case ("R")
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
       case ("C")
         do i7=1,Ny7
            do i6=1,Ny6
               do i5=1,Ny5
                  do i4=1,Ny4
                     do i3=1,Ny3
                        do i2=1,Ny2
                           do i1=1,Ny1
                              read(unit,*)Y1(i1,i2,i3,i4,i5,i6,i7)
                           enddo
                        enddo
                     enddo
                  enddo
               enddo
            enddo
         enddo
      end select
      close(unit)
      call file_gzip(reg(pname))
   end subroutine data_readA7_R

   subroutine data_readA7_C(pname,Y1,order,wspace)
      integer                             :: Ny1,Ny2,Ny3,Ny4,Ny5,Ny6,Ny7
      integer                             :: i1,i2,i3,i4,i5,i6,i7
      character(len=*)                    :: pname
      complex(8),dimension(:,:,:,:,:,:,:) :: Y1
      character(len=*),optional           :: order
      logical,optional                    :: wspace
      character(len=1)                    :: order_
      logical                             :: wspace_
      order_ = "R"   ; if(present(order))order_=trim(order(1:1))
      wspace_= .true.; if(present(wspace))wspace_=wspace
      !
      call file_gunzip(reg(pname))
      call ioread_control(pname)
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
      select case(order_)
       case default
         stop "read_array: order != Row-major, Col-major"
       case ("R")
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
       case ("C")
         do i7=1,Ny7
            do i6=1,Ny6
               do i5=1,Ny5
                  do i4=1,Ny4
                     do i3=1,Ny3
                        do i2=1,Ny2
                           do i1=1,Ny1
                              read(unit,*)Y1(i1,i2,i3,i4,i5,i6,i7)
                           enddo
                        enddo
                     enddo
                  enddo
               enddo
            enddo
         enddo
      end select
      close(unit)
      call file_gzip(reg(pname))
   end subroutine data_readA7_C

   subroutine ioread_control(pname)
      character(len=*)    :: pname
      inquire(file=trim(pname),exist=control)
      if(.not.control)inquire(file=trim(pname)//".gz",exist=control)
      if(.not.control)then
         write(*,"(A)")"I can not read : +"//trim(pname)//" SKIP"
         call sleep(5)
         return
      else
         write(*,"(A,A)")"read:     "//trim(pname)
      endif
   end subroutine ioread_control

end module IOREAD
