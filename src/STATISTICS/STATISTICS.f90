MODULE STATISTICS
  implicit none
  private

  public :: get_moments
  public :: get_covariance

contains

  

  subroutine get_moments(data,ave,sdev,var,skew,curt)
    REAL(8), INTENT(OUT)              :: ave,sdev,var,skew,curt
    REAL(8), DIMENSION(:), INTENT(IN) :: data
    INTEGER                           :: n
    REAL(8)                           :: ep,adev
    REAL(8), DIMENSION(size(data))    :: p,s
    n=size(data)
    if (n <= 1) call abort('data.size must be at least 2')
    ave=sum(data(:))/n
    s(:)=data(:)-ave
    ep=sum(s(:))
    adev=sum(abs(s(:)))/n
    p(:)=s(:)*s(:)
    var=sum(p(:))
    p(:)=p(:)*s(:)
    skew=sum(p(:))
    p(:)=p(:)*s(:)
    curt=sum(p(:))
    var=(var-ep**2/n)/(n-1)
    sdev=sqrt(var)
    if (var /= 0.0) then
       skew=skew/(n*sdev**3)
       curt=curt/(n*var**2)-3.0d0
    else
       skew=0.d0
       curt=0.d0
    endif
  end subroutine get_moments

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

