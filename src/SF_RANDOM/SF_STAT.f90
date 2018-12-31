MODULE SF_STAT
  USE SF_ARRAYS, only: linspace
  USE SF_INTEGRATE, only: simps
  USE SF_IOTOOLS, only: free_unit
  implicit none
  private

  type,public :: histogram
     integer                       :: n=0
     real(8),dimension(:),pointer  :: range
     real(8),dimension(:),pointer  :: bin
  end type histogram

  type,public :: pdf_kernel
     real(8)                          :: xmin
     real(8)                          :: xmax
     real(8)                          :: dx
     integer                          :: N=0
     real(8),dimension(:),allocatable :: x
     real(8),dimension(:),allocatable :: pdf
     real(8)                          :: sigma
     logical                          :: status=.false.
     logical                          :: variance=.false.
  end type pdf_kernel



  interface pdf_accumulate
     module procedure :: pdf_accumulate_s
     module procedure :: pdf_accumulate_v
  end interface pdf_accumulate
  
  interface pdf_print
     module procedure :: pdf_print_unit
     module procedure :: pdf_print_pfile
  end interface pdf_print

  public :: get_moments
  public :: get_mean,get_sd,get_var,get_skew,get_curt
  public :: get_covariance

  public :: pdf_allocate
  public :: pdf_deallocate
  public :: pdf_set_range
  public :: pdf_reset
  public :: pdf_set_sigma
  public :: pdf_get_range
  public :: pdf_get_pdf
  public :: pdf_accumulate
  public :: pdf_normalize
  public :: pdf_print
  public :: pdf_mean
  public :: pdf_var
  public :: pdf_sdev
  public :: pdf_skew
  public :: pdf_curt
  public :: pdf_moment


  public :: histogram_allocate
  public :: histogram_deallocate
  public :: histogram_set_range_uniform
  public :: histogram_accumulate
  public :: histogram_get_range
  public :: histogram_get_value
  public :: histogram_print
  public :: histogram_reset


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




  !### PDF
  subroutine pdf_allocate(self,a,b,N)
    type(pdf_kernel) :: self
    real(8)          :: a
    real(8)          :: b
    integer          :: N
    if(self%status)stop "PDF_ALLOCATE: PDF already allocated"
    self%xmin = a
    self%xmax = b
    self%N    = N
    allocate(self%x(N))
    allocate(self%pdf(N))
    self%status=.true.
    self%x = linspace(a,b,N,mesh=self%dx)
    self%pdf = 0d0
    self%sigma = 0d0
  end subroutine pdf_allocate

  subroutine pdf_deallocate(self)
    type(pdf_kernel) :: self
    if(.not.self%status)stop "PDF_DEALLOCATE: PDF not allocated"
    self%xmin=0d0
    self%xmax=0d0
    self%dx  =0d0
    self%sigma = 0d0
    self%N=0
    deallocate(self%x)
    deallocate(self%pdf)
    self%status=.false.
  end subroutine pdf_deallocate

  subroutine pdf_set_range(self,a,b)
    type(pdf_kernel) :: self
    real(8)          :: a
    real(8)          :: b
    if(.not.self%status)stop "PDF_SET_RANGE: PDF not allocated"
    self%xmin=a
    self%xmax=b
    self%x = linspace(a,b,self%N,mesh=self%dx)
  end subroutine pdf_set_range

  subroutine pdf_get_range(self,range)
    type(pdf_kernel)                 :: self
    real(8),dimension(:),allocatable :: range
    if(allocated(range))deallocate(range)
    allocate(range(self%N))
    range = self%x
  end subroutine pdf_get_range

  subroutine pdf_get_pdf(self,pdf)
    type(pdf_kernel)                 :: self
    real(8),dimension(:),allocatable :: pdf
    if(allocated(pdf))deallocate(pdf)
    allocate(pdf(self%N))
    pdf = self%pdf
  end subroutine pdf_get_pdf

  subroutine pdf_reset(self)
    type(pdf_kernel) :: self
    if(.not.self%status)stop "PDF_RESET: PDF not allocated"
    self%pdf=0d0
    self%sigma = 0d0
  end subroutine pdf_reset

  subroutine pdf_set_sigma(self,data)
    type(pdf_kernel)     :: self
    real(8),dimension(:) :: data
    real(8)              :: h
    if(.not.self%status)stop "PDF_SET_SIGMA: PDF not allocated"
    h = get_sd(data)
    h = (4d0/3d0/size(data))**(1d0/5d0)*h !Silverman's rule of thumb.
    self%sigma = h
    self%variance=.true.
  end subroutine pdf_set_sigma

  subroutine pdf_accumulate_s(self,data,w,sigma)
    type(pdf_kernel) :: self
    real(8)          :: data
    real(8)          :: w
    real(8),optional :: sigma
    if(.not.self%status)stop "PDF_ACCUMULATE: PDF not allocated"
    if(self%variance)then
       self%pdf   = self%pdf + gaussian_kernel(self%x,data,self%sigma)*w
    else
       if(present(sigma))then
          self%pdf   = self%pdf + gaussian_kernel(self%x,data,sigma)*w
       else
          stop "PDF_ACCUMULATE: PDF sigma not set or passed"
       endif
    endif
  end subroutine pdf_accumulate_s

  subroutine pdf_accumulate_v(self,data,sigma)
    type(pdf_kernel) :: self
    real(8)          :: data(:)
    real(8),optional :: sigma
    integer          :: i
    if(.not.self%status)stop "PDF_ACCUMULATE: PDF not allocated"
    if(self%variance)then
       do i=1,size(data)
          self%pdf   = self%pdf + gaussian_kernel(self%x,data(i),self%sigma)/size(data)
       enddo
    else
       if(present(sigma))then
          do i=1,size(data)
             self%pdf   = self%pdf + gaussian_kernel(self%x,data(i),sigma)/size(data)
          enddo
       else
          stop "PDF_ACCUMULATE: PDF sigma not set or passed"
       endif
    endif
  end subroutine pdf_accumulate_v

  subroutine pdf_normalize(self)
    type(pdf_kernel) :: self
    if(.not.self%status)stop "PDF_NORMALIZE: PDF not allocated"
    self%pdf = self%pdf/sum(self%pdf)/self%dx
  end subroutine pdf_normalize

  subroutine pdf_print_unit(self,unit)
    type(pdf_kernel) :: self    
    integer,optional :: unit
    integer          :: unit_
    integer          :: i
    if(.not.self%status)stop "PDF_WRITE: PDF not allocated"
    unit_ = 6; if(present(unit))unit_=unit
    call pdf_normalize(self)
    do i=1,self%N
       write(unit,*)self%x(i),self%pdf(i)
    enddo
  end subroutine pdf_print_unit

  subroutine pdf_print_pfile(self,pfile)
    type(pdf_kernel) :: self    
    character(len=*) :: pfile
    integer          :: i,unit
    if(.not.self%status)stop "PDF_WRITE: PDF not allocated"
    call pdf_normalize(self)
    open(free_unit(unit),file=trim(pfile))
    do i=1,self%N
       write(unit,*)self%x(i),self%pdf(i)
    enddo
    close(unit)
  end subroutine pdf_print_pfile

  function pdf_mean(self) result(mean)
    type(pdf_kernel) :: self
    real(8)          :: mean    
    mean = simps(self%x*self%pdf,self%xmin,self%xmax)
  end function pdf_mean

  function pdf_var(self) result(var)
    type(pdf_kernel) :: self
    real(8)          :: var,mean
    mean = pdf_mean(self)
    var  = simps((self%x-mean)**2*self%pdf,self%xmin,self%xmax)
  end function pdf_var

  function pdf_sdev(self) result(sdev)
    type(pdf_kernel) :: self
    real(8)          :: var,sdev
    var  = pdf_var(self)
    sdev = sqrt(var)
  end function pdf_sdev

  function pdf_moment(self,n,mu) result(mom)
    type(pdf_kernel) :: self
    integer          :: n
    real(8),optional :: mu
    real(8)          :: mom,mu_
    mu_ = 0d0;if(present(mu))mu_=mu
    mom = simps((self%x-mu_)**n*self%pdf,self%xmin,self%xmax)
  end function pdf_moment

  function pdf_skew(self) result(skew)
    type(pdf_kernel) :: self
    real(8)          :: mean,sdev,skew
    mean = pdf_mean(self)
    sdev = pdf_sdev(self)
    skew = pdf_moment(self,3)
    skew = skew - 3*mean*sdev**2 - mean**3
    skew = skew/sdev**3
  end function pdf_skew

  function pdf_curt(self) result(curt)
    type(pdf_kernel) :: self
    real(8)          :: mean,var,curt
    mean = pdf_mean(self)
    var  = pdf_var(self)
    curt = pdf_moment(self,4,mean)
    curt = curt/var**2
  end function pdf_curt

  elemental function gaussian_kernel(x,mean,sigma)
    real(8),intent(in) :: x,mean,sigma
    real(8)            :: gaussian_kernel
    real(8),parameter  :: pi2=2d0*acos(-1d0)
    gaussian_kernel=exp(-0.5d0*((x-mean)/sigma)**2)/sqrt(pi2)/sigma
  end function gaussian_kernel










  !### HISTOGRAMS
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


END MODULE SF_STAT

