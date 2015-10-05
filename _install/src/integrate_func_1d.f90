!+-----------------------------------------------------------------+
!PURPOSE:
! Trapezoidal rule for 1d function integration between a and b.
!+-----------------------------------------------------------------+
function d_trapz_func(f,a,b,N,x) result(int)
  interface
     function f(x)
       real(8) :: x
       real(8) :: f
     end function f
  end interface
  real(8),optional                 :: a,b
  real(8),dimension(:),optional    :: x
  integer,optional                 :: N
  real(8)                          :: dh,a_,b_
  integer                          :: L,i
  real(8),dimension(:),allocatable :: xx
  real(8)                          :: int
  L  = 100;if(present(N))L =N;if(present(x))L=size(x)
  a_ = 0d0;if(present(a))a_=a
  b_ = 1d0;if(present(b))b_=b
  !
  int=0.d0
  if(.not.present(x))then
     allocate(xx(L))
     xx = linspace(a_,b_,L,mesh=dh)
     do i=1,L-1
        int = int+( f(xx(i+1)) + f(xx(i)) )
     enddo
     int = int*dh/2d0
     deallocate(xx)
  else
     do i=1,L-1
        dh  = (x(i+1)-x(i))/2d0
        int = int + ( f(x(i+1)) + f(x(i)) )*dh
     enddo
  endif
end function d_trapz_func
function c_trapz_func(f,a,b,N,x) result(int)
  interface
     function f(x)
       real(8)    :: x
       complex(8) :: f
     end function f
  end interface
  real(8),optional                 :: a,b
  real(8),dimension(:),optional    :: x
  integer,optional                 :: N
  real(8)                          :: dh,a_,b_
  integer                          :: L,i
  real(8),dimension(:),allocatable :: xx
  complex(8)                       :: int
  L  = 100;if(present(N))L =N;if(present(x))L=size(x)
  a_ = 0d0;if(present(a))a_=a
  b_ = 1d0;if(present(b))b_=b
  !
  int=0.d0
  if(.not.present(x))then
     allocate(xx(L))
     xx = linspace(a_,b_,L,mesh=dh)
     do i=1,L-1
        int = int+( f(xx(i+1)) + f(xx(i)) )
     enddo
     int = int*dh/2d0
     deallocate(xx)
  else
     do i=1,L-1
        dh  = (x(i+1)-x(i))/2d0
        int = int + ( f(x(i+1)) + f(x(i)) )*dh
     enddo
  endif
end function c_trapz_func



!+-----------------------------------------------------------------+
!PURPOSE:
! Simpson's rule for 1d function integration between a and b.
!+-----------------------------------------------------------------+
function d_simps_func(f,a,b,N,x) result(int)
  interface
     function f(x)
       real(8) :: x
       real(8) :: f
     end function f
  end interface
  real(8),optional                 :: a,b
  real(8),dimension(:),optional    :: x
  integer,optional                 :: N
  real(8)                          :: dh,a_,b_
  integer                          :: L,M,i
  real(8),dimension(:),allocatable :: xx,wt,dx
  real(8)                          :: int,int1,int2,int3
  L  = 100;if(present(N))L =N;if(present(x))L=size(x)
  a_ = 0d0;if(present(a))a_=a
  b_ = 1d0;if(present(b))b_=b
  !
  int=0.d0
  if(.not.present(x))then
     allocate(xx(L),wt(L))
     xx = linspace(a_,b_,L,mesh=dh)
     call get_quadrature_weights(wt)
     do i=1,L
        int = int + f(xx(i))*wt(i)
     enddo
     int = int*dh
     deallocate(xx,wt)
  else
     allocate(dx(L+1))
     M=L-1
     dx=0d0
     forall(i=1:M)dx(i)=x(i+1)-x(i)
     int1=0d0;int2=0d0;int3=0d0
     i=0
     do while(i < L)
        if( (dx(i)==dx(i+1)) .AND. (dx(i)==dx(i+2)) )then !Simpson's 3/8 rule
           int1 = int1 + ( 3*dx(i)*( f(x(i)) + &
                3*( f(x(i+1)) + f(x(i+2)) ) + f(x(i+3)) ))/8
           i=i+3
        elseif(dx(i)==dx(i+1))then !Simpson's 1/3 rule
           int2 = int2+( 2*dx(i)*( f(x(i)) + &
                4*f(x(i+1)) + f(x(i+2))) )/6
           i=i+2
        elseif(dx(i)/=dx(i+1)) then !trapezoidal rule
           int3 = int3 + dx(i)*( f(x(i))+f(x(i+1)) )/2.d0
           i = i + 1
        endif
     enddo
     int = int1+int2+int3
  endif
end function d_simps_func
function c_simps_func(f,a,b,N,x) result(int)
  interface
     function f(x)
       real(8)    :: x
       complex(8) :: f
     end function f
  end interface
  real(8),optional                 :: a,b
  real(8),dimension(:),optional    :: x
  integer,optional                 :: N
  real(8)                          :: dh,a_,b_
  integer                          :: L,M,i
  real(8),dimension(:),allocatable :: xx,wt,dx
  complex(8)                       :: int,int1,int2,int3
  L  = 100;if(present(N))L =N;if(present(x))L=size(x)
  a_ = 0d0;if(present(a))a_=a
  b_ = 1d0;if(present(b))b_=b
  !
  int=0.d0
  if(.not.present(x))then
     allocate(xx(L),wt(L))
     xx = linspace(a_,b_,L,mesh=dh)
     call get_quadrature_weights(wt)
     do i=1,L
        int = int + f(xx(i))*wt(i)
     enddo
     int = int*dh
     deallocate(xx,wt)
  else
     allocate(dx(L+1))
     M=L-1
     dx=0d0
     forall(i=1:M)dx(i)=x(i+1)-x(i)
     int1=0d0;int2=0d0;int3=0d0
     i=0
     do while(i < L)
        if( (dx(i)==dx(i+1)) .AND. (dx(i)==dx(i+2)) )then !Simpson's 3/8 rule
           int1 = int1 + ( 3*dx(i)*( f(x(i)) + &
                3*( f(x(i+1)) + f(x(i+2)) ) + f(x(i+3)) ))/8
           i=i+3
        elseif(dx(i)==dx(i+1))then !Simpson's 1/3 rule
           int2 = int2+( 2*dx(i)*( f(x(i)) + &
                4*f(x(i+1)) + f(x(i+2))) )/6
           i=i+2
        elseif(dx(i)/=dx(i+1)) then !trapezoidal rule
           int3 = int3 + dx(i)*( f(x(i))+f(x(i+1)) )/2.d0
           i = i + 1
        endif
     enddo
     int = int1+int2+int3
  endif
end function c_simps_func

