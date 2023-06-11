function histogram_allocate(n) result(h)
  integer,intent(in) :: n
  type(histogram)    :: h
  if(n<=0)then
     print*,"histogram length must be positive integer. n=",n
     stop
  endif
  allocate(h%range(0:n));h%range=0.d0
  allocate(h%bin(0:n))    ;h%bin=0.d0
  h%n=n
end function histogram_allocate


subroutine histogram_deallocate(h)
  type(histogram)    :: h
  deallocate(h%range)
  deallocate(h%bin)
  h%n=0
end subroutine  histogram_deallocate


subroutine histogram_reset(h)
  type(histogram)    :: h
  h%range=0.d0
  h%bin=0.d0
end subroutine histogram_reset


subroutine histogram_set_range_uniform(h,xmin,xmax)
  type(histogram),intent(inout) :: h
  real(8),intent(in)            :: xmin,xmax
  integer                       :: i,n
  real(8)                       :: f1,f2
  if(xmin>=xmax)then
     print*,"histogram%range: xmin must be less than xmax:",xmin,xmax
     stop
  endif
  n=h%n
  do i=0,n
     f1= real(n-i,8)/real(n,8)
     f2= real(i,8)/real(n,8)
     h%range(i) = f1*xmin + f2*xmax
  enddo
  h%bin=0.d0
end subroutine histogram_set_range_uniform


subroutine histogram_accumulate(h,x,w)
  type(histogram),intent(inout) :: h
  real(8),intent(in)            :: x,w
  integer                       :: i,index
  index=0
  call find_index(h%n,h%range,x,index)
  if(index>=h%n)then
     print*,"index lies outside valid range of 0 .. n - 1"
     stop
  endif
  h%bin(index)=h%bin(index)+w
end subroutine histogram_accumulate


subroutine find_index(n,range,x,index)
  integer,intent(in)                :: n
  real(8),dimension(0:n),intent(in) :: range
  real(8),intent(in)                :: x
  integer,intent(out)               :: index
  integer                           :: i,upper,lower,mid
  if((x<range(0)) .OR. (x>range(n)))then
     print*,"X out of range!"
     return
  endif
  upper=n
  lower=0
  do while((upper-lower>1))
     mid = (upper+lower)/2    !int(dble(upper + lower)/2.d0)
     if( x >= range(mid))then
        lower=mid
     else
        upper=mid
     endif
  enddo
  index=lower
  if(x<range(lower) .OR. x>range(lower+1))then
     print*,"error: x not found within range!"
     stop
  endif
end subroutine find_index


subroutine histogram_get_range(h,index,lower,upper)
  type(histogram),intent(in) :: h
  integer,intent(in)         :: index
  real(8),intent(out)        :: lower,upper
  if(index>=h%n)then
     print*,"error: *i lies outside valid range=0...n-1"
     stop
  endif
  lower=h%range(index)
  upper=h%range(index+1)
end subroutine histogram_get_range


function histogram_get_value(h,index) result(value)
  type(histogram),intent(in) :: h
  integer,intent(in)         :: index
  real(8)                    :: value
  if(index>=h%n)then
     print*,"error: *index lies outside valid range=0...n-1"
     stop
  endif
  value=h%bin(index)
end function histogram_get_value


subroutine histogram_print(h,unit)
  type(histogram),intent(in) :: h
  integer,intent(in)         :: unit
  integer                    :: i,n
  real(8)                    :: lower,upper,bin_value
  n=h%n
  call histogram_get_range(h,0,lower,upper)
  write(unit,"(2F12.7)")lower,0.d0
  do i=0,n-1
     call histogram_get_range(h,i,lower,upper)
     bin_value = histogram_get_value(h,i)
     write(unit,"(2F12.7)")lower,bin_value
     write(unit,"(2F12.7)")upper,bin_value
  enddo
  write(unit,"(2F12.7)")upper,0.d0
end subroutine histogram_print

