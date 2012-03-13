MODULE STATISTICS
  implicit none
  private

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

END MODULE STATISTICS

