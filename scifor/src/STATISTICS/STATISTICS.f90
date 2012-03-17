MODULE STATISTICS
  implicit none
  private

  type,public :: histogram
     integer                      :: n
     real(8),dimension(:),pointer :: range
     real(8),dimension(:),pointer  :: bin
  end type histogram

  public :: histogram_allocate
  public :: histogram_set_range_uniform
  public :: histogram_accumulate
  public :: histogram_get_range
  public :: histogram_get_value
  public :: histogram_print


  public :: get_moments
  public :: get_mean,get_sd,get_var,get_skew,get_curt
  public :: get_covariance


contains


  subroutine get_moments(data,ave,sdev,var,skew,curt)
    real(8), intent(out),optional     :: ave,sdev,var,skew,curt
    real(8)                           :: ave_,sdev_,var_,skew_,curt_
    real(8), dimension(:), intent(in) :: data
    integer                           :: n
    real(8)                           :: ep,adev
    real(8), dimension(size(data))    :: p,s
    n=size(data)
    if (n <= 1) then
       print*,'data.size must be at least 2'
       stop
    endif
    ave_=sum(data(:))/n
    s(:)=data(:)-ave_
    ep=sum(s(:))
    adev=sum(abs(s(:)))/n
    p(:)=s(:)*s(:)
    var_=sum(p(:))
    p(:)=p(:)*s(:)
    skew_=sum(p(:))
    p(:)=p(:)*s(:)
    curt_=sum(p(:))
    var_=(var_-ep**2/n)/(n-1)
    sdev_=sqrt(var_)
    if (var_ /= 0.0) then
       skew_=skew_/(n*sdev_**3)
       curt_=curt_/(n*var_**2)-3.0d0
    else
       skew_=0.d0
       curt_=0.d0
    endif
    if(present(ave))ave=ave_
    if(present(sdev))sdev=sdev_
    if(present(var))var=var_
    if(present(skew))skew=skew_
    if(present(curt))curt=curt_
  end subroutine get_moments


  function get_mean(data) result(mean)
    real(8), dimension(:), intent(in) :: data
    real(8)                           :: mean
    call get_moments(data,ave=mean)
  end function get_mean


  function get_sd(data) result(sd)
    real(8), dimension(:), intent(in) :: data
    real(8)                           :: sd
    call get_moments(data,sdev=sd)
  end function get_sd


  function get_var(data) result(var)
    real(8), dimension(:), intent(in) :: data
    real(8)                           :: var
    call get_moments(data,var=var)
  end function get_var


  function get_skew(data) result(skew)
    real(8), dimension(:), intent(in) :: data
    real(8)                           :: skew
    call get_moments(data,skew=skew)
  end function get_skew


  function get_curt(data) result(curt)
    real(8), dimension(:), intent(in) :: data
    real(8)                           :: curt
    call get_moments(data,curt=curt)
  end function get_curt


  function get_covariance(data,mean) result(covariance)
    real(8),dimension(:,:),intent(in)            :: data
    real(8),dimension(size(data,1)),intent(in)   :: mean
    real(8),dimension(size(data,1),size(data,1)) :: covariance
    real(8),dimension(size(data,2))              :: Xi,Xj
    integer                                      :: i,j,ncol,L
    ncol=size(data,1)
    L   =size(data,2)
    covariance=0.d0
    do i=1,ncol
       Xi(:)=data(i,:)-mean(i)
       do j=1,ncol
          Xj(:)=data(j,:)-mean(j)
          covariance(i,j) = sum(Xi(:)*Xj(:))/real(L-1,8)
       enddo
    enddo
  end function get_covariance



  function histogram_allocate(n) result(h)
    integer,intent(in) :: n
    type(histogram)    :: h
    if(n<=0)then
       print*,"histogram length must be positive integer. n=",n
       stop
    endif
    allocate(h%range(0:n));h%range=0.d0
    allocate(h%bin(n))    ;h%bin=0.d0
    h%n=n
  end function histogram_allocate


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
    if((x<range(0)) .OR. (x>=range(n)))then
       print*,"X out of range!"
       return
    endif
    upper=n
    lower=0
    do while((upper-lower>1))
       mid = int(dble(upper + lower)/2.d0)
       if( x >= range(mid))then
          lower=mid
       else
          upper=mid
       endif
    enddo
    index=lower
    if(x<range(lower) .OR. x>=range(lower+1))then
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
    do i=0,n-1
       call histogram_get_range(h,i,lower,upper)
       bin_value = histogram_get_value(h,i)
       write(unit,"(2F12.7)")lower,bin_value
       write(unit,"(2F12.7)")upper,bin_value
    enddo
  end subroutine histogram_print


END MODULE STATISTICS

