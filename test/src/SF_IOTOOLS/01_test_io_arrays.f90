program testIOTOOLS
  USE SF_IOTOOLS
  USE ASSERTING
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

  const=acos(-1.d0)

  x=const

  call fill_with_drand(A1)
  call fill_with_drand(A2)
  call fill_with_drand(A3)
  call fill_with_drand(A4)
  call fill_with_drand(A5)
  call fill_with_drand(A6)
  call fill_with_drand(A7)

  call fill_with_crand(B1)
  call fill_with_crand(B2)
  call fill_with_crand(B3)
  call fill_with_crand(B4)
  call fill_with_crand(B5)
  call fill_with_crand(B6)
  call fill_with_crand(B7)



  print*,""
  print*,"             SPLOT <--> SREAD                   "
  print*,""

  call splot("A1.dat",x,A1)
  call sread("A1.dat",x_,A1_)
  call assert(x,x_,"COORDINATES")
  call assert(A1,A1_,"1D REAL ARRAY")
  call splot("B1.dat",x,B1)
  call sread("B1.dat",x_,B1_)
  call assert(x,x_,"COORDINATES")
  call assert(B1,B1_,"1D COMPLEX ARRAY")

  call splot("A2.dat",x,A2)
  call sread("A2.dat",x_,A2_)
  call assert(x,x_,"COORDINATES")
  call assert(A2,A2_,"2D REAL ARRAY")
  call splot("B2.dat",x,B2)
  call sread("B2.dat",x_,B2_)
  call assert(x,x_,"COORDINATES")
  call assert(B2,B2_,"2D COMPLEX ARRAY")

  call splot("A3.dat",x,A3)
  call sread("A3.dat",x_,A3_)
  call assert(x,x_,"COORDINATES")
  call assert(A3,A3_,"3D REAL ARRAY")
  call splot("B3.dat",x,B3)
  call sread("B3.dat",x_,B3_)
  call assert(x,x_,"COORDINATES")
  call assert(B3,B3_,"3D COMPLEX ARRAY")

  call splot("A4.dat",x,A4)
  call sread("A4.dat",x_,A4_)
  call assert(x,x_,"COORDINATES")
  call assert(A4,A4_,"4D REAL ARRAY")
  call splot("B4.dat",x,B4)
  call sread("B4.dat",x_,B4_)
  call assert(x,x_,"COORDINATES")
  call assert(B4,B4_,"4D COMPLEX ARRAY")

  call splot("A5.dat",x,A5)
  call sread("A5.dat",x_,A5_)
  call assert(x,x_,"COORDINATES")
  call assert(A5,A5_,"5D REAL ARRAY")
  call splot("B5.dat",x,B5)
  call sread("B5.dat",x_,B5_)
  call assert(x,x_,"COORDINATES")
  call assert(B5,B5_,"5D COMPLEX ARRAY")

  call splot("A6.dat",x,A6)
  call sread("A6.dat",x_,A6_)
  call assert(x,x_,"COORDINATES")
  call assert(A6,A6_,"6D REAL ARRAY")
  call splot("B6.dat",x,B6)
  call sread("B6.dat",x_,B6_)
  call assert(x,x_,"COORDINATES")
  call assert(B6,B6_,"6D COMPLEX ARRAY")

  call splot("A7.dat",x,A7)
  call sread("A7.dat",x_,A7_)
  call assert(x,x_,"COORDINATES")
  call assert(A7,A7_,"7D REAL ARRAY")
  call splot("B7.dat",x,B7)
  call sread("B7.dat",x_,B7_)
  call assert(x,x_,"COORDINATES")
  call assert(B7,B7_,"7D COMPLEX ARRAY")

  

  print*,""
  print*,"             SAVE <--> READ                   "
  print*,""

  call set_store_size(10)

  call save_array("A1.dat",A1)
  call read_array("A1.dat",A1_)
  call assert(A1,A1_,"1D REAL ARRAY")
  call save_array("B1.dat",B1)
  call read_array("B1.dat",B1_)
  call assert(B1,B1_,"1D COMPLEX ARRAY")

  call save_array("A2.dat",A2)
  call read_array("A2.dat",A2_)
  call assert(A2,A2_,"2D REAL ARRAY")
  call save_array("B2.dat",B2)
  call read_array("B2.dat",B2_)
  call assert(B2,B2_,"2D COMPLEX ARRAY")

  call save_array("A3.dat",A3)
  call read_array("A3.dat",A3_)
  call assert(A3,A3_,"3D REAL ARRAY")
  call save_array("B3.dat",B3)
  call read_array("B3.dat",B3_)
  call assert(B3,B3_,"3D COMPLEX ARRAY")

  call save_array("A4.dat",A4)
  call read_array("A4.dat",A4_)
  call assert(A4,A4_,"4D REAL ARRAY")
  call save_array("B4.dat",B4)
  call read_array("B4.dat",B4_)
  call assert(B4,B4_,"4D COMPLEX ARRAY")

  call save_array("A5.dat",A5)
  call read_array("A5.dat",A5_)
  call assert(A5,A5_,"5D REAL ARRAY")
  call save_array("B5.dat",B5)
  call read_array("B5.dat",B5_)
  call assert(B5,B5_,"5D COMPLEX ARRAY")

  call save_array("A6.dat",A6)
  call read_array("A6.dat",A6_)
  call assert(A6,A6_,"6D REAL ARRAY")
  call save_array("B6.dat",B6)
  call read_array("B6.dat",B6_)
  call assert(B6,B6_,"6D COMPLEX ARRAY")

  call save_array("A7.dat",A7)
  call read_array("A7.dat",A7_)
  call assert(A7,A7_,"7D REAL ARRAY")
  call save_array("B7.dat",B7)
  call read_array("B7.dat",B7_)
  call assert(B7,B7_,"7D COMPLEX ARRAY")



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

  impure elemental subroutine fill_with_drand(A)
    real(8),intent(inout) :: A
    A = drand()
  end subroutine

  impure elemental subroutine fill_with_crand(B)
    complex(8),intent(inout) :: B
    B = crand()
  end subroutine
  
end program testIOTOOLS