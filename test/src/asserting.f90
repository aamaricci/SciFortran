module asserting
  implicit none

  public :: assert_arr_d, assert_d, assert_i, assert_c
contains


  subroutine assert_i(a,b,info)
    integer, intent(in) :: a,b
    logical, intent(out) :: info

    info=.false.
    if(a==b) info=.true.
  end subroutine assert_i

  subroutine assert_d(a,b,info,prec)
    real(8), intent(in) :: a,b
    real(8), intent(in), optional :: prec
    logical, intent(out) :: info
    real(8) :: prec_=1.d-5
    if(present(prec)) prec_=prec

    info=.false.
    if(abs(a-b)<prec_) info=.true.
  end subroutine assert_d

  

  subroutine assert_c(a,b,info)
    character(len=*), intent(in) :: a,b
    logical, intent(out) :: info

    info=.false.
    if(a==b) info=.true.
  end subroutine assert_c
  
  subroutine assert_arr_d(a,b,info,prec)
    real(8), intent(in), dimension(:) :: a,b
    real(8), intent(in), optional :: prec
    logical, intent(out) :: info
    real(8) :: prec_=1.d-5
    if(present(prec)) prec_=prec
    
    info=.false.
    if(all(abs(a-b)<10.d0**(-14))) info=.true.
  end subroutine assert_arr_d
  
end module asserting
