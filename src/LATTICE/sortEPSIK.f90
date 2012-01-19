!+-----------------------------------------------------------------+
!PROGRAM  : 
!TYPE     : Subroutine
!PURPOSE  :   
!Copied from TOOLS.f90 to sort the epsik array and leave this 
!module independent
!+-----------------------------------------------------------------+
subroutine sort(array,order2,no_touch_array)
  implicit none
  real(8),dimension(:)                    :: array
  real(8),dimension(size(array))          :: backup
  integer,dimension(size(array))          :: order
  integer,dimension(size(array)),optional :: order2
  integer                                 :: i
  logical,optional                        :: no_touch_array
  do i=1,size(order)
     order(i)=i
  enddo
  call qsort_sort_( array, order, 1, size(array) )
  if(.not.present(no_touch_array))then
     do i=1,size(order)
        backup(i)=array(order(i))
     enddo
     array=backup
  endif
  if(present(order2)) order2=order
contains
  recursive subroutine qsort_sort_( array, order, left, right )
    implicit none
    real(8), dimension(:)         :: array
    integer, dimension(:)         :: order
    integer                       :: left
    integer                       :: right
    integer                       :: i
    integer                       :: last
    if ( left .ge. right ) return
    call qsort_swap_( order, left, qsort_rand_(left,right) )
    last = left
    do i = left+1, right
       if ( compare_(array(order(i)), array(order(left)) ) .lt. 0 ) then
          last = last + 1
          call qsort_swap_( order, last, i )
       endif
    enddo
    call qsort_swap_( order, left, last )
    call qsort_sort_( array, order, left, last-1 )
    call qsort_sort_( array, order, last+1, right )
  end subroutine qsort_sort_
  !---------------------------------------------!
  !---------------------------------------------!
  !---------------------------------------------!
  subroutine qsort_swap_( order, first, second )
    implicit none
    integer, dimension(:)         :: order
    integer                       :: first, second
    integer                       :: tmp
    tmp           = order(first)
    order(first)  = order(second)
    order(second) = tmp
  end subroutine qsort_swap_
  !---------------------------------------------!
  !---------------------------------------------!
  !---------------------------------------------!
  integer function qsort_rand_( lower, upper )
    implicit none
    integer                       :: lower, upper
    real(4)                       :: r
    r=drand_()
    qsort_rand_ =  lower + nint(r * (upper-lower))
  end function qsort_rand_
  !---------------------------------------------!
  !---------------------------------------------!
  !---------------------------------------------!
  function drand_()
    implicit none
    real(8) :: drand_
    real(4) :: r
    call random_number(r)
    drand_=dble(r)
  end function drand_
  !---------------------------------------------!
  !---------------------------------------------!
  !---------------------------------------------!
  function compare_(f,g)
    implicit none
    real(8) :: f,g
    integer :: compare_
    if(f<g) then
       compare_=-1
    else
       compare_=1
    endif
  end function compare_
end subroutine sort
!***************************************************************
!***************************************************************
!***************************************************************
