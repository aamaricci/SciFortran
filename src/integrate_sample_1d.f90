!+-----------------------------------------------------------------+
!PURPOSE:
! Trapezoidal rule for 1d data function integration between a and b.
! a,b default are 0,1
! + _ab: given a,b and f(:) integrate f(:)
! + _dh: given dh=(b-a)/L-1/2 integrate f(:)
! + _nonlin: integrate f(:) using given x(:)
!+-----------------------------------------------------------------+
function d_trapz_sample(f,dh,a,b,x) result(sum)
  real(8),dimension(:)          :: f
  real(8),optional              :: dh,a,b
  real(8),dimension(:),optional :: x
  real(8)                       :: h,a_,b_
  real(8)                       :: sum
  a_=0d0;if(present(a))a_=a
  b_=1d0;if(present(b))b_=b
  h = (b_-a_)/(size(f)-1)
  if(present(dh))h=dh
  if(.not.present(x))then
     sum = d_trapz_dh_sample(h,f)
  else
     sum = d_trapz_nonlin_sample(x,f)
  endif
end function d_trapz_sample
function c_trapz_sample(f,dh,a,b,x) result(sum)
  complex(8),dimension(:)       :: f
  real(8),optional              :: dh,a,b
  real(8),dimension(:),optional :: x
  real(8)                       :: h,a_,b_
  complex(8)                    :: sum
  a_=0d0;if(present(a))a_=a
  b_=1d0;if(present(b))b_=b
  h = (b_-a_)/(size(f)-1)
  if(present(dh))h=dh
  if(.not.present(x))then
     sum = c_trapz_dh_sample(h,f)
  else
     sum = c_trapz_nonlin_sample(x,f)
  endif
end function c_trapz_sample




function d_trapz_dh_sample(dh,f) result(sum)
  real(8) :: dh
  real(8) :: f(:)
  real(8) :: sum
  integer :: i,L
  L=size(f)
  sum=0d0
  do i=1,L-1
     sum = sum+(f(i+1)+f(i))*dh/2.d0
  enddo
end function d_trapz_dh_sample
function c_trapz_dh_sample(dh,f) result(sum)
  real(8)    :: dh
  complex(8) :: f(:)
  complex(8) :: sum
  integer    :: i,L
  L=size(f)
  sum=0.d0
  do i=1,L-1
     sum = sum+(f(i+1)+f(i))*dh/2.d0
  enddo
end function c_trapz_dh_sample

function d_trapz_ab_sample(a,b,f) result(sum)
  real(8) :: a,b,dh
  real(8) :: f(:)
  real(8) :: sum
  integer :: i,L
  L=size(f)
  dh=(b-a)/(L-1)/2.d0
  sum=0d0
  do i=1,L-1
     sum = sum+(f(i+1)+f(i))*dh
  enddo
end function d_trapz_ab_sample
function c_trapz_ab_sample(a,b,f) result(sum)
  real(8)    :: a,b,dh
  complex(8) :: f(:)
  complex(8) :: sum
  integer    :: i,L
  L=size(f)
  dh=(b-a)/(L-1)/2.d0
  sum=0d0
  do i=1,L-1
     sum = sum+(f(i+1)+f(i))*dh
  enddo
end function c_trapz_ab_sample


function d_trapz_nonlin_sample(x,f) result(sum)
  real(8) :: a,b,dh
  real(8) :: f(:),x(size(f))
  real(8) :: sum
  integer :: i,L
  L=size(f)
  a=minval(x)
  b=maxval(x)
  sum=0.d0
  do i=1,L-1
     dh  = (x(i+1)-x(i))/2.d0
     sum = sum + (f(i+1)+f(i))*dh
  enddo
end function d_trapz_nonlin_sample
function c_trapz_nonlin_sample(x,f) result(sum)
  real(8)    :: a,b,dh
  complex(8) :: f(:)
  real(8)    :: x(size(f))
  complex(8) :: sum
  integer    :: i,L
  L=size(f)
  a=minval(x)
  b=maxval(x)
  sum=0.d0
  do i=1,L-1
     dh  = (x(i+1)-x(i))/2.d0
     sum = sum + (f(i+1)+f(i))*dh
  enddo
end function c_trapz_nonlin_sample







!+-----------------------------------------------------------------+
!PURPOSE: Simpson rule for data function integration
! + int_simps: integrate f(:) using quadrature weights
! + _ab: given a,b and f(:) integrate f(:)
! + _dh: fiven dh=(b-a)/L-1/2 integrate f(:)
! + _nonlin: integrate f(:) using given x(:)
!+-----------------------------------------------------------------+
function d_simps_sample(f,dh,a,b,x) result(int_value)
  real(8),dimension(:)                :: f
  real(8),optional                    :: dh
  real(8),optional                    :: a,b
  real(8),dimension(size(f)),optional :: x
  real(8),dimension(size(f))          :: wt
  real(8)                             :: a_,b_,dh_
  real(8)                             :: int_value
  a_=0d0;if(present(a))a_=a
  b_=1d0;if(present(b))b_=b
  dh_=(b_-a_)/(size(f)-1);if(present(dh))dh_=dh
  if(.not.present(x))then
     call get_quadrature_weights(wt)
     int_value = sum(f(:)*wt(:))*dh_
  else
     int_value = d_simpson_nonlin_sample(x,f)
  endif
end function  d_simps_sample
function c_simps_sample(f,dh,a,b,x) result(int_value)
  complex(8),dimension(:)             :: f
  real(8),optional                    :: dh
  real(8),optional                    :: a,b
  real(8),dimension(size(f)),optional :: x
  real(8),dimension(size(f))          :: wt
  real(8)                             :: a_,b_,dh_
  complex(8)                          :: int_value
  a_=0d0;if(present(a))a_=a
  b_=1d0;if(present(b))b_=b
  dh_=(b_-a_)/(size(f)-1);if(present(dh))dh_=dh
  if(.not.present(x))then
     call get_quadrature_weights(wt)
     int_value = sum(f(:)*wt(:))*dh_
  else
     int_value = c_simpson_nonlin_sample(x,f)
  endif
end function  c_simps_sample


function d_simpson_dh_sample(dh,f) result(sum)
  integer :: n
  real(8) :: f(:)
  real(8) :: dh,sum,sum1,sum2,int1,int2
  integer :: i,p,m,mm,mmm
  N=size(f)
  if(N==1)then
     sum=0.d0
     return
  endif
  sum1=0.d0
  sum2=0.d0
  sum =0.d0
  int1=0.d0
  int2=0.d0
  if(mod(n-1,2)==0)then                !if n-1 is even:
     do i=2,N-1,2
        sum1 = sum1 + f(i)
     enddo
     do i=3,N-2,2
        sum2 = sum2 + f(i)
     enddo
     sum = (f(1) + 4.d0*sum1 + 2.d0*sum2 + f(n))*dh/3.d0
  else                        !if n-1 is odd, use Simpson's for N-3 slices + 3/8rule for the last
     if (N>=6) then
        do i=2,N-4,2
           sum1 = sum1 + f(i)
        enddo
        do i=3,N-5,2
           sum2 = sum2 + f(i)
        enddo
        int1 = (f(1) + 4.d0*sum1 + 2.d0*sum2 + f(n-3))*dh/3.d0
     endif
     int2 = (f(n-3)+3.d0*f(n-2)+3.d0*f(n-1)+f(n))*dh*3.d0/8.d0
     sum  = int1 + int2
  end if
end function d_simpson_dh_sample
function c_simpson_dh_sample(dh,f) result(sum)
  integer              :: n
  complex(8)           :: f(:)
  real(8)              :: dh
  complex(8)           :: sum,sum1,sum2,int1,int2
  integer              :: i,p,m
  complex(8),parameter :: zero=cmplx(0.d0,0.d0,8)
  N=size(f)
  if(N==1)then
     sum=zero
     return
  endif
  sum1=zero
  sum2=zero
  sum =zero
  int1=zero
  int2=zero
  if(mod(n-1,2)==0)then                !if n-1 is even:
     do i=2,N-1,2
        sum1 = sum1 + f(i)
     enddo
     do i=3,N-2,2
        sum2 = sum2 + f(i)
     enddo
     sum = (f(1) + 4.d0*sum1 + 2.d0*sum2 + f(n))*dh/3.d0
  else                        !if n-1 is odd, use Simpson's for N-3 slices + 3/8rule for the last
     if (N>=6) then
        do i=2,N-4,2
           sum1 = sum1 + f(i)
        enddo
        do i=3,N-5,2
           sum2 = sum2 + f(i)
        enddo
        int1 = (f(1) + 4.d0*sum1 + 2.d0*sum2 + f(n-3))*dh/3.d0
     endif
     int2 = (f(n-3)+3.d0*f(n-2)+3.d0*f(n-1)+f(n))*dh*3.d0/8.d0
     sum  = int1 + int2
  end if
end function c_simpson_dh_sample


function d_simpson_ab_sample(a,b,f) result(sum)
  real(8) :: dh,a,b
  real(8) :: f(:)
  real(8) :: sum
  integer :: L
  L=size(f)
  dh=(b-a)/real(L-1,8)
  sum = d_simpson_dh_sample(dh,f)
end function d_simpson_ab_sample
function c_simpson_ab_sample(a,b,f) result(sum)
  real(8)    :: dh,a,b
  complex(8) :: f(:)
  complex(8) :: sum
  integer    :: L
  L=size(f)
  dh=(b-a)/real(L-1,8)
  sum = c_simpson_dh_sample(dh,f)
end function c_simpson_ab_sample


function d_simpson_nonlin_sample(x,f) result(sum)
  real(8) :: f(:),x(size(f)),dx(size(f)+1)
  real(8) :: sum,sum1,sum2,sum3
  real(8) :: a,b,dh
  integer :: i,n,m
  n=size(f)
  m=n-1
  dx=0.d0
  forall(i=1:m)dx(i)=x(i+1)-x(i)
  sum1=0.d0
  sum2=0.d0
  sum3=0.d0
  i=0
  do while(i<n)
     !Simpson's 3/8 rule
     if((dx(i)==dx(i+1)).AND.(dx(i)==dx(i+2)))then
        sum1=sum1+(3.d0*dx(i)*(f(i)+&
             3.d0*(f(i+1)+f(i+2))+f(i+3)))/8.d0
        i=i+3
        !Simpson's 1/3 rule
     elseif(dx(i)==dx(i+1))then
        sum2=sum2+(2.d0*dx(i)*(f(i)+&
             4.d0*f(i+1)+f(i+2)))/6.d0
        i=i+2
        !trapezoidal rule
     elseif(dx(i)/=dx(i+1)) then
        sum3=sum3+dx(i)*(f(i)+f(i+1))/2.d0
        i = i + 1
     endif
  enddo
  sum = sum1+sum2+sum3
end function d_simpson_nonlin_sample
function c_simpson_nonlin_sample(x,f) result(sum)
  complex(8) :: f(:),sum
  real(8)    :: x(size(f)),rsum,isum
  rsum=d_simpson_nonlin_sample(x,dreal(f))
  isum=d_simpson_nonlin_sample(x,dimag(f))
  sum  = dcmplx(rsum,isum)
end function c_simpson_nonlin_sample








