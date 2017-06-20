program testSLPLOTeREAD
  USE IOFILE
  USE IOPLOT
  USE IOREAD
  implicit none

  integer,parameter :: L=5
  integer :: unit

  real(8),dimension(L) :: X,X_
  real(8),dimension(L) :: A1,A1_
  real(8),dimension(L,L) :: A2,A2_
  real(8),dimension(L,L,L) :: A3,A3_
  real(8),dimension(L,L,L,L) :: A4,A4_
  real(8),dimension(L,L,L,L,L) :: A5,A5_
  real(8),dimension(L,L,L,L,L,L) :: A6,A6_
  real(8),dimension(L,L,L,L,L,L,L) :: A7,A7_


  complex(8),dimension(L) :: B1,B1_
  complex(8),dimension(L,L) :: B2,B2_
  complex(8),dimension(L,L,L) :: B3,B3_
  complex(8),dimension(L,L,L,L) :: B4,B4_
  complex(8),dimension(L,L,L,L,L) :: B5,B5_
  complex(8),dimension(L,L,L,L,L,L) :: B6,B6_
  complex(8),dimension(L,L,L,L,L,L,L) :: B7,B7_

  real(8) :: const

  logical :: bool

  ! inquire(file="kk.dat",opened=bool,number=unit)
  ! print*,bool,unit
  ! open(100,file="kk.dat")
  ! inquire(file="kk.dat",opened=bool,number=unit)
  ! print*,bool,unit
  ! if(bool)close(100)
  ! inquire(file="kk.dat",opened=bool,number=unit)
  ! print*,bool,unit


  ! !>DEBUG
  ! stop
  ! !<DEBUG

  const=acos(-1.d0)

  x=const

  A1=drand()
  A2=drand()
  A3=drand()
  A4=drand()
  A5=drand()
  A6=drand()
  A7=drand()

  B1=crand()
  B2=crand()
  B3=crand()
  B4=crand()
  B5=crand()
  B6=crand()
  B7=crand()


  print*,"             SPLOT <--> SREAD                   "
  print*,"         Sum(x-x_)   ","               Sum(A-A_)"
  print*,""


  call splot("A1.dat",x,A1)
  call sread("A1.dat",x_,A1_)
  print*,"A1:",sum(x-x_),sum(A1-A1_)
  call splot("B1.dat",x,B1)
  call sread("B1.dat",x_,B1_)
  print*,"B1:",sum(x-x_),sum(B1-B1_)

  call splot("A2.dat",x,A2)
  call sread("A2.dat",x_,A2_)
  print*,"A2:",sum(x-x_),sum(A2-A2_)
  call splot("B2.dat",x,B2)
  call sread("B2.dat",x_,B2_)
  print*,"B2:",sum(x-x_),sum(B2-B2_)

  call splot("A3.dat",x,A3)
  call sread("A3.dat",x_,A3_)
  print*,"A3:",sum(x-x_),sum(A3-A3_)
  call splot("B3.dat",x,B3)
  call sread("B3.dat",x_,B3_)
  print*,"B3:",sum(x-x_),sum(B3-B3_)

  call splot("A4.dat",x,A4)
  call sread("A4.dat",x_,A4_)
  print*,"A4:",sum(x-x_),sum(A4-A4_)
  call splot("B4.dat",x,B4)
  call sread("B4.dat",x_,B4_)
  print*,"B4:",sum(x-x_),sum(B4-B4_)

  call splot("A5.dat",x,A5)
  call sread("A5.dat",x_,A5_)
  print*,"A5:",sum(x-x_),sum(A5-A5_)
  call splot("B5.dat",x,B5)
  call sread("B5.dat",x_,B5_)
  print*,"B5:",sum(x-x_),sum(B5-B5_)

  call splot("A6.dat",x,A6)
  call sread("A6.dat",x_,A6_)
  print*,"A6:",sum(x-x_),sum(A6-A6_)
  call splot("B6.dat",x,B6)
  call sread("B6.dat",x_,B6_)
  print*,"B6:",sum(x-x_),sum(B6-B6_)

  call splot("A7.dat",x,A7)
  call sread("A7.dat",x_,A7_)
  print*,"A7:",sum(x-x_),sum(A7-A7_)
  call splot("B7.dat",x,B7)
  call sread("B7.dat",x_,B7_)
  print*,"B7:",sum(x-x_),sum(B7-B7_)



  call set_store_size(10)

  print*,""
  print*,""
  print*,"             SAVE <--> READ                   "
  print*,""
  call save_array("A1.dat",A1)
  call read_array("A1.dat",A1_)
  print*,"A1:",sum(A1-A1_)
  call save_array("B1.dat",B1)
  call read_array("B1.dat",B1_)
  print*,"B1:",sum(B1-B1_)

  call save_array("A2.dat",A2)
  call read_array("A2.dat",A2_)
  print*,"A2:",sum(A2-A2_)
  call save_array("B2.dat",B2)
  call read_array("B2.dat",B2_)
  print*,"B2:",sum(B2-B2_)

  call save_array("A3.dat",A3)
  call read_array("A3.dat",A3_)
  print*,"A3:",sum(A3-A3_)
  call save_array("B3.dat",B3)
  call read_array("B3.dat",B3_)
  print*,"B3:",sum(B3-B3_)

  call save_array("A4.dat",A4)
  call read_array("A4.dat",A4_)
  print*,"A4:",sum(A4-A4_)
  call save_array("B4.dat",B4)
  call read_array("B4.dat",B4_)
  print*,"B4:",sum(B4-B4_)

  call save_array("A5.dat",A5)
  call read_array("A5.dat",A5_)
  print*,"A5:",sum(A5-A5_)
  call save_array("B5.dat",B5)
  call read_array("B5.dat",B5_)
  print*,"B5:",sum(B5-B5_)

  call save_array("A6.dat",A6)
  call read_array("A6.dat",A6_)
  print*,"A6:",sum(A6-A6_)
  call save_array("B6.dat",B6)
  call read_array("B6.dat",B6_)
  print*,"B6:",sum(B6-B6_)

  call save_array("A7.dat",A7)
  call read_array("A7.dat",A7_)
  print*,"A7:",sum(A7-A7_)
  call save_array("B7.dat",B7)
  call read_array("B7.dat",B7_)
  print*,"B7:",sum(B7-B7_)



contains

  function drand() result(r)
    real(8) :: r
    call random_number(r)
  end function drand

  function crand() result(r)
    real(8)    :: re,im
    complex(8) :: r
    call random_number(re)
    call random_number(im)
    r=dcmplx(re,im)
  end function crand

end program testSLPLOTeREAD
