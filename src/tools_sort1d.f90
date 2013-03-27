

!+-----------------------------------------------------------------+
!PROGRAM  : 
!TYPE     : Subroutine
!PURPOSE  :   
!+-----------------------------------------------------------------+
subroutine i_uniq(AIN,AOUT,MASK)
  !USE M_UNISTA
  integer,dimension(:),intent(INOUT)           :: AIN
  integer,dimension(:),allocatable,intent(OUT) :: AOUT
  integer                                      :: NDIM
  logical,dimension(:),allocatable,intent(OUT)&
       ,optional                               :: MASK
  if(present(MASK))then
     allocate(MASK(size(AIN)))
     call i_unista(AIN,ndim,MASK)
  else
     call i_unista(AIN,ndim)
  endif
  allocate(AOUT(NDIM))
  AOUT(1:NDIM)=AIN(1:NDIM)
end subroutine i_uniq
subroutine d_uniq(AIN,AOUT,MASK)
  !USE M_UNISTA
  real(8),dimension(:),intent(INOUT)           :: AIN
  real(8),dimension(:),allocatable,intent(OUT) :: AOUT
  integer                                      :: NDIM
  logical,dimension(:),allocatable,intent(OUT)&
       ,optional                               :: MASK
  if(present(MASK))then
     allocate(MASK(size(AIN)))
     call d_unista(AIN,ndim,MASK)
  else
     call d_unista(AIN,ndim)
  endif
  allocate(AOUT(NDIM))
  AOUT(1:NDIM)=AIN(1:NDIM)
end subroutine d_uniq
!*******************************************************************
!*******************************************************************
!*******************************************************************



!+-----------------------------------------------------------------+
!PROGRAM  : 
!TYPE     : Subroutine
!PURPOSE  :   
!+-----------------------------------------------------------------+
subroutine reshuffle(array,order) 
  real(8),dimension(:)                    :: array
  integer,dimension(size(array))          :: order
  real(8),dimension(size(array))          :: dummy
  integer                                 :: i,Lk
  Lk=size(array)
  forall(i=1:Lk)dummy(i)=array(order(i))
  array=dummy
end subroutine reshuffle
!*******************************************************************
!*******************************************************************
!*******************************************************************




!+-----------------------------------------------------------------+
!PURPOSE  :   
!+-----------------------------------------------------------------+
subroutine sort(ARR,M)
  implicit none
  integer :: i,j,M
  real(8) :: a
  real(8),dimension(M) :: ARR
  do j=2, M
     a=ARR(j)
     do i=j-1,1,-1
        if (ARR(i)<=a) goto 10
        ARR(i+1)=ARR(i)
     enddo
     i=0
10   ARR(i+1)=a
  enddo
  return
end subroutine sort
!*******************************************************************
!*******************************************************************
!*******************************************************************


!+------------------------------------------------------------------+
!PURPOSE  : Sort an array, gives the new ordering of the label.
!+------------------------------------------------------------------+
subroutine sort_array(array,order2,no_touch_array)
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
  call qsort_sort( array, order, 1, size(array) )
  if(.not.present(no_touch_array))then
     do i=1,size(order)
        backup(i)=array(order(i))
     enddo
     array=backup
  endif
  if(present(order2)) order2=order
contains
  recursive subroutine qsort_sort( array, order, left, right )
    implicit none
    real(8), dimension(:)                 :: array
    integer, dimension(:)                 :: order
    integer                               :: left
    integer                               :: right
    integer                               :: i
    integer                               :: last
    if ( left .ge. right ) return
    call qsort_swap( order, left, qsort_rand(left,right) )
    last = left
    do i = left+1, right
       if ( compare(array(order(i)), array(order(left)) ) .lt. 0 ) then
          last = last + 1
          call qsort_swap( order, last, i )
       endif
    enddo
    call qsort_swap( order, left, last )
    call qsort_sort( array, order, left, last-1 )
    call qsort_sort( array, order, last+1, right )
  end subroutine qsort_sort
  !---------------------------------------------!
  subroutine qsort_swap( order, first, second )
    implicit none
    integer, dimension(:)                 :: order
    integer                               :: first, second
    integer                               :: tmp
    tmp           = order(first)
    order(first)  = order(second)
    order(second) = tmp
  end subroutine qsort_swap
  !---------------------------------------------!
  integer function qsort_rand( lower, upper )
    implicit none
    integer                               :: lower, upper
    real(4)                               :: r
    r=drand()
    qsort_rand =  lower + nint(r * (upper-lower))
  end function qsort_rand
  !---------------------------------------------!
  function drand()
    implicit none
    real(8)                               :: drand
    real(4)                               :: r
    call random_number(r)
    drand=dble(r)
  end function drand
  !---------------------------------------------!
  function compare(f,g)
    implicit none
    real(8)                               :: f,g
    integer                               :: compare
    if(f<g) then
       compare=-1
    else
       compare=1
    endif
  end function compare
end subroutine sort_array
!*******************************************************************
!*******************************************************************
!*******************************************************************
