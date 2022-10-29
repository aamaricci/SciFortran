MODULE ASSERTING
  USE SCIFOR, only: bold_red,bold_green,str
  implicit none
  private

  interface assert
     module procedure :: assert_i0
     module procedure :: assert_d0
     module procedure :: assert_z0
     module procedure :: assert_ch0
     module procedure :: assert_b0
     !
     module procedure :: assert_i1
     module procedure :: assert_d1
     module procedure :: assert_z1
     module procedure :: assert_ch1
     module procedure :: assert_b1
     !
     module procedure :: assert_i2
     module procedure :: assert_d2
     module procedure :: assert_z2
     module procedure :: assert_ch2
     module procedure :: assert_b2
     !
     module procedure :: assert_i3
     module procedure :: assert_d3
     module procedure :: assert_z3
     module procedure :: assert_ch3
     module procedure :: assert_b3
     !
     module procedure :: assert_i4
     module procedure :: assert_d4
     module procedure :: assert_z4
     module procedure :: assert_ch4
     module procedure :: assert_b4
     !
     module procedure :: assert_i5
     module procedure :: assert_d5
     module procedure :: assert_z5
     module procedure :: assert_ch5
     module procedure :: assert_b5
     !
     module procedure :: assert_i6
     module procedure :: assert_d6
     module procedure :: assert_z6
     module procedure :: assert_ch6
     module procedure :: assert_b6
     !
     module procedure :: assert_i7
     module procedure :: assert_d7
     module procedure :: assert_z7
     module procedure :: assert_ch7
     module procedure :: assert_b7
  end interface assert

  public :: assert


  logical :: test
  real(8) :: tol_

contains



  subroutine assert_i0(a,b,fname)
    integer, intent(in)         :: a,b
    character(len=*),intent(in) :: fname
    test=.false.
    if(a==b)test=.true.
    call assert_msg(fname,test)
  end subroutine assert_i0

  subroutine assert_d0(a,b,fname,tol)
    real(8), intent(in)           :: a,b
    character(len=*),intent(in)   :: fname
    real(8), intent(in), optional :: tol
    tol_=1d-12 ; if(present(tol)) tol_=tol
    test=.false.
    if(abs(a-b)<tol_)test=.true.
    call assert_msg(fname,test)
  end subroutine assert_d0

  subroutine assert_z0(a,b,fname,tol)
    complex(8), intent(in)        :: a,b
    character(len=*),intent(in)   :: fname
    real(8), intent(in), optional :: tol
    tol_=1d-12 ; if(present(tol)) tol_=tol
    test=.false.
    if(abs(a-b)<tol_)test=.true.
    call assert_msg(fname,test)
  end subroutine assert_z0

  subroutine assert_ch0(a,b,fname)
    character(len=*), intent(in)  :: a,b
    character(len=*),intent(in)   :: fname
    test=.false.
    if(a==b)test=.true.
    call assert_msg(fname,test)
  end subroutine assert_ch0

  subroutine assert_b0(a,b,fname)
    logical, intent(in)         :: a,b
    character(len=*),intent(in) :: fname
    test=.false.
    if(a.eqv.b)test=.true.
    call assert_msg(fname,test)
  end subroutine assert_b0


  !#############################################################

  subroutine assert_i1(a,b,fname)
    integer,intent(in),dimension(:) :: a,b
    character(len=*),intent(in)     :: fname
    test=.false.
    if(all(a==b)) test=.true.
    call assert_msg(fname,test)
  end subroutine assert_i1

  subroutine assert_d1(a,b,fname,tol)
    real(8),intent(in),dimension(:) :: a,b
    character(len=*),intent(in)     :: fname
    real(8),intent(in),optional     :: tol
    tol_=1d-12 ; if(present(tol)) tol_=tol
    test=.false.
    if(all(abs(a-b)<tol_)) test=.true.
    call assert_msg(fname,test)
  end subroutine assert_d1

  subroutine assert_z1(a,b,fname,tol)
    complex(8),intent(in),dimension(:) :: a,b
    character(len=*),intent(in)        :: fname
    real(8),intent(in),optional        :: tol
    tol_=1d-12 ; if(present(tol)) tol_=tol
    test=.false.
    if(all(abs(a-b)<tol_)) test=.true.
    call assert_msg(fname,test)
  end subroutine assert_z1

  subroutine assert_ch1(a,b,fname)
    character(len=*),intent(in),dimension(:) :: a,b
    character(len=*),intent(in)              :: fname
    test=.false.
    if(all(a==b)) test=.true.
    call assert_msg(fname,test)
  end subroutine assert_ch1

  subroutine assert_b1(a,b,fname)
    logical,intent(in),dimension(:) :: a,b
    character(len=*),intent(in)     :: fname
    test=.false.
    if(all(a.eqv.b)) test=.true.
    call assert_msg(fname,test)
  end subroutine assert_b1


  !#############################################################


  subroutine assert_i2(a,b,fname)
    integer,intent(in),dimension(:,:) :: a,b
    character(len=*),intent(in)       :: fname
    test=.false.
    if(all(a==b)) test=.true.
    call assert_msg(fname,test)
  end subroutine assert_i2

  subroutine assert_d2(a,b,fname,tol)
    real(8),intent(in),dimension(:,:) :: a,b
    character(len=*),intent(in)     :: fname
    real(8),intent(in),optional     :: tol
    tol_=1d-12 ; if(present(tol)) tol_=tol
    test=.false.
    if(all(abs(a-b)<tol_)) test=.true.
    call assert_msg(fname,test)
  end subroutine assert_d2

  subroutine assert_z2(a,b,fname,tol)
    complex(8),intent(in),dimension(:,:) :: a,b
    character(len=*),intent(in)        :: fname
    real(8),intent(in),optional        :: tol
    tol_=1d-12 ; if(present(tol)) tol_=tol
    test=.false.
    if(all(abs(a-b)<tol_)) test=.true.
    call assert_msg(fname,test)
  end subroutine assert_z2

  subroutine assert_ch2(a,b,fname)
    character(len=*),intent(in),dimension(:,:) :: a,b
    character(len=*),intent(in)              :: fname
    test=.false.
    if(all(a==b)) test=.true.
    call assert_msg(fname,test)
  end subroutine assert_ch2

  subroutine assert_b2(a,b,fname)
    logical,intent(in),dimension(:,:) :: a,b
    character(len=*),intent(in)     :: fname
    test=.false.
    if(all(a.eqv.b)) test=.true.
    call assert_msg(fname,test)
  end subroutine assert_b2


  !#############################################################


  subroutine assert_i3(a,b,fname)
    integer,intent(in),dimension(:,:,:) :: a,b
    character(len=*),intent(in)       :: fname
    test=.false.
    if(all(a==b)) test=.true.
    call assert_msg(fname,test)
  end subroutine assert_i3

  subroutine assert_d3(a,b,fname,tol)
    real(8),intent(in),dimension(:,:,:) :: a,b
    character(len=*),intent(in)     :: fname
    real(8),intent(in),optional     :: tol
    tol_=1d-12 ; if(present(tol)) tol_=tol
    test=.false.
    if(all(abs(a-b)<tol_)) test=.true.
    call assert_msg(fname,test)
  end subroutine assert_d3

  subroutine assert_z3(a,b,fname,tol)
    complex(8),intent(in),dimension(:,:,:) :: a,b
    character(len=*),intent(in)        :: fname
    real(8),intent(in),optional        :: tol
    tol_=1d-12 ; if(present(tol)) tol_=tol
    test=.false.
    if(all(abs(a-b)<tol_)) test=.true.
    call assert_msg(fname,test)
  end subroutine assert_z3

  subroutine assert_ch3(a,b,fname)
    character(len=*),intent(in),dimension(:,:,:) :: a,b
    character(len=*),intent(in)              :: fname
    test=.false.
    if(all(a==b)) test=.true.
    call assert_msg(fname,test)
  end subroutine assert_ch3

  subroutine assert_b3(a,b,fname)
    logical,intent(in),dimension(:,:,:) :: a,b
    character(len=*),intent(in)     :: fname
    test=.false.
    if(all(a.eqv.b)) test=.true.
    call assert_msg(fname,test)
  end subroutine assert_b3


  !#############################################################


  subroutine assert_i4(a,b,fname)
    integer,intent(in),dimension(:,:,:,:) :: a,b
    character(len=*),intent(in)       :: fname
    test=.false.
    if(all(a==b)) test=.true.
    call assert_msg(fname,test)
  end subroutine assert_i4

  subroutine assert_d4(a,b,fname,tol)
    real(8),intent(in),dimension(:,:,:,:) :: a,b
    character(len=*),intent(in)     :: fname
    real(8),intent(in),optional     :: tol
    tol_=1d-12 ; if(present(tol)) tol_=tol
    test=.false.
    if(all(abs(a-b)<tol_)) test=.true.
    call assert_msg(fname,test)
  end subroutine assert_d4

  subroutine assert_z4(a,b,fname,tol)
    complex(8),intent(in),dimension(:,:,:,:) :: a,b
    character(len=*),intent(in)        :: fname
    real(8),intent(in),optional        :: tol
    tol_=1d-12 ; if(present(tol)) tol_=tol
    test=.false.
    if(all(abs(a-b)<tol_)) test=.true.
    call assert_msg(fname,test)
  end subroutine assert_z4

  subroutine assert_ch4(a,b,fname)
    character(len=*),intent(in),dimension(:,:,:,:) :: a,b
    character(len=*),intent(in)              :: fname
    test=.false.
    if(all(a==b)) test=.true.
    call assert_msg(fname,test)
  end subroutine assert_ch4

  subroutine assert_b4(a,b,fname)
    logical,intent(in),dimension(:,:,:,:) :: a,b
    character(len=*),intent(in)     :: fname
    test=.false.
    if(all(a.eqv.b)) test=.true.
    call assert_msg(fname,test)
  end subroutine assert_b4


  !#############################################################



  subroutine assert_i5(a,b,fname)
    integer,intent(in),dimension(:,:,:,:,:) :: a,b
    character(len=*),intent(in)       :: fname
    test=.false.
    if(all(a==b)) test=.true.
    call assert_msg(fname,test)
  end subroutine assert_i5

  subroutine assert_d5(a,b,fname,tol)
    real(8),intent(in),dimension(:,:,:,:,:) :: a,b
    character(len=*),intent(in)     :: fname
    real(8),intent(in),optional     :: tol
    tol_=1d-12 ; if(present(tol)) tol_=tol
    test=.false.
    if(all(abs(a-b)<tol_)) test=.true.
    call assert_msg(fname,test)
  end subroutine assert_d5

  subroutine assert_z5(a,b,fname,tol)
    complex(8),intent(in),dimension(:,:,:,:,:) :: a,b
    character(len=*),intent(in)        :: fname
    real(8),intent(in),optional        :: tol
    tol_=1d-12 ; if(present(tol)) tol_=tol
    test=.false.
    if(all(abs(a-b)<tol_)) test=.true.
    call assert_msg(fname,test)
  end subroutine assert_z5

  subroutine assert_ch5(a,b,fname)
    character(len=*),intent(in),dimension(:,:,:,:,:) :: a,b
    character(len=*),intent(in)              :: fname
    test=.false.
    if(all(a==b)) test=.true.
    call assert_msg(fname,test)
  end subroutine assert_ch5

  subroutine assert_b5(a,b,fname)
    logical,intent(in),dimension(:,:,:,:,:) :: a,b
    character(len=*),intent(in)     :: fname
    test=.false.
    if(all(a.eqv.b)) test=.true.
    call assert_msg(fname,test)
  end subroutine assert_b5


  !#############################################################



  subroutine assert_i6(a,b,fname)
    integer,intent(in),dimension(:,:,:,:,:,:) :: a,b
    character(len=*),intent(in)       :: fname
    test=.false.
    if(all(a==b)) test=.true.
    call assert_msg(fname,test)
  end subroutine assert_i6

  subroutine assert_d6(a,b,fname,tol)
    real(8),intent(in),dimension(:,:,:,:,:,:) :: a,b
    character(len=*),intent(in)     :: fname
    real(8),intent(in),optional     :: tol
    tol_=1d-12 ; if(present(tol)) tol_=tol
    test=.false.
    if(all(abs(a-b)<tol_)) test=.true.
    call assert_msg(fname,test)
  end subroutine assert_d6

  subroutine assert_z6(a,b,fname,tol)
    complex(8),intent(in),dimension(:,:,:,:,:,:) :: a,b
    character(len=*),intent(in)        :: fname
    real(8),intent(in),optional        :: tol
    tol_=1d-12 ; if(present(tol)) tol_=tol
    test=.false.
    if(all(abs(a-b)<tol_)) test=.true.
    call assert_msg(fname,test)
  end subroutine assert_z6

  subroutine assert_ch6(a,b,fname)
    character(len=*),intent(in),dimension(:,:,:,:,:,:) :: a,b
    character(len=*),intent(in)              :: fname
    test=.false.
    if(all(a==b)) test=.true.
    call assert_msg(fname,test)
  end subroutine assert_ch6

  subroutine assert_b6(a,b,fname)
    logical,intent(in),dimension(:,:,:,:,:,:) :: a,b
    character(len=*),intent(in)     :: fname
    test=.false.
    if(all(a.eqv.b)) test=.true.
    call assert_msg(fname,test)
  end subroutine assert_b6


  !#############################################################

  subroutine assert_i7(a,b,fname)
    integer,intent(in),dimension(:,:,:,:,:,:,:) :: a,b
    character(len=*),intent(in)       :: fname
    test=.false.
    if(all(a==b)) test=.true.
    call assert_msg(fname,test)
  end subroutine assert_i7

  subroutine assert_d7(a,b,fname,tol)
    real(8),intent(in),dimension(:,:,:,:,:,:,:) :: a,b
    character(len=*),intent(in)     :: fname
    real(8),intent(in),optional     :: tol
    tol_=1d-12 ; if(present(tol)) tol_=tol
    test=.false.
    if(all(abs(a-b)<tol_)) test=.true.
    call assert_msg(fname,test)
  end subroutine assert_d7

  subroutine assert_z7(a,b,fname,tol)
    complex(8),intent(in),dimension(:,:,:,:,:,:,:) :: a,b
    character(len=*),intent(in)        :: fname
    real(8),intent(in),optional        :: tol
    tol_=1d-12 ; if(present(tol)) tol_=tol
    test=.false.
    if(all(abs(a-b)<tol_)) test=.true.
    call assert_msg(fname,test)
  end subroutine assert_z7

  subroutine assert_ch7(a,b,fname)
    character(len=*),intent(in),dimension(:,:,:,:,:,:,:) :: a,b
    character(len=*),intent(in)              :: fname
    test=.false.
    if(all(a==b)) test=.true.
    call assert_msg(fname,test)
  end subroutine assert_ch7

  subroutine assert_b7(a,b,fname)
    logical,intent(in),dimension(:,:,:,:,:,:,:) :: a,b
    character(len=*),intent(in)     :: fname
    test=.false.
    if(all(a.eqv.b)) test=.true.
    call assert_msg(fname,test)
  end subroutine assert_b7


  !#############################################################

  subroutine assert_msg(fname,info)
    character(len=*),intent(in) :: fname
    logical,intent(in)          :: info
    if(info) then  
       write(*,"(A50,A)")str(fname),": exit status is "//bold_green("POSITIVE")
    else
       write(*,"(A50,A)")str(fname),": exit status is "//bold_red("NEGATIVE")
       error stop 2
    endif
  end subroutine assert_msg


end module asserting
