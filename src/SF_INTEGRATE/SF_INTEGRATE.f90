! ABOUT QUADPACK:
! How to decide what routine to use, if your integration region is finite:
!
! If you can factor the integrand as F(X)=W(X)*G(X), where G is smooth on [A,B] and W(X)=COS(OMEGA*X) or SIN(OMEGA*X) then use QAWO.
! Otherwise, if you can factor F(X)=W(X)*G(X) where G is smooth and W(X)=(X-A)**ALFA * (B-X)**BETA * (LOG(X-A))**L * (LOG(B-X))**K with K, L = 0 or 1, and ALFA, BETA greater than -1, then use QAWS.
!
! Otherwise, if you can factor F(X)=W(X)*G(X) where G is smooth and W(X)=1/(X-C) for some constant C, use QAWC.
! Otherwise, if you do not care too much about possible inefficient use of computer time, and do not want to further analyze the problem, use QAGS.
! Otherwise, if the integrand is smooth, use QNG or QAG.
! Otherwise, if there are discontinuities or singularities of the integrand or of its derivative, and you know where they are, split the integration range at these points and analyze each subinterval. You can also use QAGP, which is to be provided with the x-locations of the singularities or discontinuities.
! Otherwise, if the integrand has end point singularities, use QAGS.
! Otherwise, if the integrand has an oscillatory behavior of nonspecific type, and no singularities, use QAG with KEY=6.
! Otherwise, use QAGS.
!
! Routines for an infinite region:
! If the integrand decays rapidly to zero, truncate the interval and use the finite interval decision tree.
! Otherwise, if the integrand oscillates over the entire infinite range, and the integral is a Fourier transform, use QAWF.
! Or, if the integrand oscillates over the entire infinite range, but is not a Fourier transform, then sum the successive positive and negative contributions by integrating between the zeroes of the integrand. Apply convergence acceleration with QELG.
! Otherwise, if you are not constrained by computer time, and do not wish to analyze the problem further, use QAGI.
! Otherwise, if the integrand has a non-smooth behavior in the range, and you know where it occurs, split off these regions and use the appropriate finite range routines to integrate over them. Then begin this tree again to handle the remainder of the region.
! Otherwise, truncation of the interval, or application of a suitable transformation for reducing the problem to a finite range may be possible. And you may also call QAGI.

module SF_INTEGRATE
   USE GAUSS_QUADRATURE
   implicit none
   private


   complex(8),parameter :: zero=(0.d0,0.d0)
   complex(8),parameter :: xi=(0.d0,1.d0)
   complex(8),parameter :: one=(1.d0,0.d0)
   real(8),parameter    :: pi    = 3.14159265358979323846264338327950288419716939937510d0


   !TRAPEZIODAL RULE:
   interface trapz
      !SAMPLE:
      module procedure d_trapz_ab_sample
      module procedure c_trapz_ab_sample
      module procedure d_trapz_dh_sample
      module procedure c_trapz_dh_sample
      module procedure d_trapz_nonlin_sample
      module procedure c_trapz_nonlin_sample
      !FUNCTION:
      module procedure d_trapz_ab_func
      module procedure c_trapz_ab_func
      module procedure d_trapz_nonlin_func
      module procedure c_trapz_nonlin_func
   end interface trapz
   !
   interface trapz2d
      module procedure d_trapz2d_func
      module procedure c_trapz2d_func
      module procedure d_trapz2d_func_recursive
      module procedure c_trapz2d_func_recursive
      module procedure d_trapz2d_sample
      module procedure c_trapz2d_sample
   end interface trapz2d





   !SIMPSON'S RULE
   interface simps
      !SAMPLE:
      module procedure d_simpson_dh_sample
      module procedure c_simpson_dh_sample
      module procedure d_simpson_ab_sample
      module procedure c_simpson_ab_sample
      module procedure d_simpson_nonlin_sample
      module procedure c_simpson_nonlin_sample
      !FUNCTION:
      module procedure d_simps_ab_func
      module procedure c_simps_ab_func
      module procedure d_simps_nonlin_func
      module procedure c_simps_nonlin_func
   end interface simps
   !
   interface simps2d
      module procedure d_simps2d_func
      module procedure c_simps2d_func
      module procedure d_simps2d_func_recursive
      module procedure c_simps2d_func_recursive
      module procedure d_simps2d_sample
      module procedure c_simps2d_sample
   end interface simps2d
   !


   !QUADPACK global interface
   interface quad
      module procedure :: quad_func
      module procedure :: quad_sample
   end interface quad



   !1D Adaptive QUADPACK
   public :: quad
   !nD Adapative GAUSS RULE 6-14
   public :: gauss_quad
   public :: integrate

   !1D simple
   public :: trapz
   public :: simps
   !2D simple
   public :: trapz2d
   public :: simps2d

   !KRAMERS-KRONIG:
   public :: kronig

   !AUX:
   public :: get_quadrature_weights



   !<TODO
   ! add cubature methods for arbitrary domains
   ! add Montecarlo base 1d/2d/3d integrals
   !>TODO



contains



!+-----------------------------------------------------------------------------+!
!PURPOSE:
! evaluate 1D and 2D integrals using trapz (2nd order) and simps (4th order) rule.
!+-----------------------------------------------------------------------------+!

   !+-----------------------------------------------------------------+
   !PURPOSE:
   ! Trapezoidal rule for 1d function integration between a and b.
   !+-----------------------------------------------------------------+
   function d_trapz_ab_func(f,a,b,N) result(int)
      interface
         function f(x)
            real(8) :: x
            real(8) :: f
         end function f
      end interface
      real(8),optional                 :: a,b
      integer,optional                 :: N
      real(8)                          :: dh,a_,b_
      integer                          :: L,i
      real(8),dimension(:),allocatable :: xx
      real(8)                          :: int
      L  = 100;if(present(N))L =N
      a_ = 0d0;if(present(a))a_=a
      b_ = 1d0;if(present(b))b_=b
      !
      int=0d0
      allocate(xx(L))
      xx = linspace(a_,b_,L,mesh=dh)
      do i=1,L-1
         int = int+( f(xx(i+1)) + f(xx(i)) )
      enddo
      int = int*dh/2d0
      deallocate(xx)
   end function d_trapz_ab_func
   function c_trapz_ab_func(f,a,b,N) result(int)
      interface
         function f(x)
            real(8)    :: x
            complex(8) :: f
         end function f
      end interface
      real(8),optional                 :: a,b
      integer,optional                 :: N
      real(8)                          :: dh,a_,b_
      integer                          :: L,i
      real(8),dimension(:),allocatable :: xx
      complex(8)                       :: int
      L  = 100;if(present(N))L =N
      a_ = 0d0;if(present(a))a_=a
      b_ = 1d0;if(present(b))b_=b
      !
      int=0.d0
      allocate(xx(L))
      xx = linspace(a_,b_,L,mesh=dh)
      do i=1,L-1
         int = int+( f(xx(i+1)) + f(xx(i)) )
      enddo
      int = int*dh/2d0
      deallocate(xx)
   end function c_trapz_ab_func


   function d_trapz_nonlin_func(f,x) result(int)
      interface
         function f(x)
            real(8) :: x
            real(8) :: f
         end function f
      end interface
      real(8),dimension(:) :: x
      real(8)              :: dh
      integer              :: L,i
      real(8)              :: int
      L  = size(x)
      !
      int=0.d0
      do i=1,L-1
         dh  = (x(i+1)-x(i))/2d0
         int = int + ( f(x(i+1)) + f(x(i)) )*dh
      enddo
   end function d_trapz_nonlin_func
   function c_trapz_nonlin_func(f,x) result(int)
      interface
         function f(x)
            real(8)    :: x
            complex(8) :: f
         end function f
      end interface
      real(8),dimension(:) :: x
      real(8)              :: dh
      integer              :: L,i
      complex(8)           :: int
      L  = size(x)
      !
      int=0.d0
      do i=1,L-1
         dh  = (x(i+1)-x(i))/2d0
         int = int + ( f(x(i+1)) + f(x(i)) )*dh
      enddo
   end function c_trapz_nonlin_func


   !+-----------------------------------------------------------------+
   !PURPOSE:
   ! Simpson's rule for 1d function integration between a and b.
   !+-----------------------------------------------------------------+
   function d_simps_ab_func(f,a,b,N) result(int)
      interface
         function f(x)
            real(8) :: x
            real(8) :: f
         end function f
      end interface
      real(8),optional                 :: a,b
      integer,optional                 :: N
      real(8)                          :: dh,a_,b_
      integer                          :: L,M,i
      real(8),dimension(:),allocatable :: xx,wt,dx
      real(8)                          :: int,int1,int2,int3
      L  = 100;if(present(N))L =N
      a_ = 0d0;if(present(a))a_=a
      b_ = 1d0;if(present(b))b_=b
      !
      int=0.d0
      allocate(xx(L),wt(L))
      xx = linspace(a_,b_,L,mesh=dh)
      call get_quadrature_weights(wt)
      do i=1,L
         int = int + f(xx(i))*wt(i)
      enddo
      int = int*dh
      deallocate(xx,wt)
   end function d_simps_ab_func
   function c_simps_ab_func(f,a,b,N) result(int)
      interface
         function f(x)
            real(8)    :: x
            complex(8) :: f
         end function f
      end interface
      real(8),optional                 :: a,b
      integer,optional                 :: N
      real(8)                          :: dh,a_,b_
      integer                          :: L,M,i
      real(8),dimension(:),allocatable :: xx,wt,dx
      complex(8)                       :: int,int1,int2,int3
      L  = 100;if(present(N))L =N
      a_ = 0d0;if(present(a))a_=a
      b_ = 1d0;if(present(b))b_=b
      !
      int=0.d0
      allocate(xx(L),wt(L))
      xx = linspace(a_,b_,L,mesh=dh)
      call get_quadrature_weights(wt)
      do i=1,L
         int = int + f(xx(i))*wt(i)
      enddo
      int = int*dh
      deallocate(xx,wt)
   end function c_simps_ab_func

   function d_simps_nonlin_func(f,x) result(int)
      interface
         function f(x)
            real(8) :: x
            real(8) :: f
         end function f
      end interface
      real(8),dimension(:)             :: x
      real(8)                          :: dh
      integer                          :: L,M,i
      real(8),dimension(:),allocatable :: xx,wt,dx
      real(8)                          :: int,int1,int2,int3
      L  = size(x)
      !
      int=0.d0
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
   end function d_simps_nonlin_func
   function c_simps_nonlin_func(f,x) result(int)
      interface
         function f(x)
            real(8)    :: x
            complex(8) :: f
         end function f
      end interface
      real(8),dimension(:)             :: x
      real(8)                          :: dh
      integer                          :: L,M,i
      real(8),dimension(:),allocatable :: xx,wt,dx
      complex(8)                       :: int,int1,int2,int3
      L  = size(x)
      !
      int=0.d0
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
   end function c_simps_nonlin_func

   !+-----------------------------------------------------------------+
   !PURPOSE:
   ! Trapezoidal rule for 1d data function integration between a and b
   ! or with respect to a given dh
   !
   ! + _ab: given a,b and f(:) integrate f(:)
   ! + _dh: given dh=(b-a)/L-1/2 integrate f(:)
   ! + _nonlin: integrate f(:) using given x(:)
   !+-----------------------------------------------------------------+
   function d_trapz_ab_sample(f,a,b) result(sum)
      real(8) :: f(:)
      real(8) :: a,b,dh
      real(8) :: sum
      integer :: i,L
      L=size(f)
      dh=(b-a)/(L-1)/2d0
      sum=0d0
      do i=1,L-1
         sum = sum+(f(i+1)+f(i))*dh
      enddo
   end function d_trapz_ab_sample
   function c_trapz_ab_sample(f,a,b) result(sum)
      complex(8) :: f(:)
      real(8)    :: a,b,dh
      complex(8) :: sum
      integer    :: i,L
      L=size(f)
      dh=(b-a)/(L-1)/2.d0
      sum=0d0
      do i=1,L-1
         sum = sum+(f(i+1)+f(i))*dh
      enddo
   end function c_trapz_ab_sample

   function d_trapz_dh_sample(f,dh) result(sum)
      real(8) :: f(:)
      real(8) :: dh
      real(8) :: sum
      integer :: i,L
      L=size(f)
      sum=0d0
      do i=1,L-1
         sum = sum+(f(i+1)+f(i))*dh/2d0
      enddo
   end function d_trapz_dh_sample
   function c_trapz_dh_sample(f,dh) result(sum)
      complex(8) :: f(:)
      real(8)    :: dh
      complex(8) :: sum
      integer    :: i,L
      L=size(f)
      sum=0.d0
      do i=1,L-1
         sum = sum+(f(i+1)+f(i))*dh/2.d0
      enddo
   end function c_trapz_dh_sample

   function d_trapz_nonlin_sample(f,x) result(sum)
      real(8) :: f(:)
      real(8) :: x(size(f))
      real(8) :: a,b,dh
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
   function c_trapz_nonlin_sample(f,x) result(sum)
      complex(8) :: f(:)
      real(8)    :: x(size(f))
      real(8)    :: a,b,dh
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
   !PURPOSE: Simpson rule for data function integration between extrema
   ! [a,b] or with respect to a given dh
   !
   ! + _ab: given a,b and f(:) integrate f(:)
   ! + _dh: fiven dh=(b-a)/L-1/2 integrate f(:)
   ! + _nonlin: integrate f(:) using given x(:)
   !+-----------------------------------------------------------------+
   function d_simpson_ab_sample(f,a,b) result(sum)
      real(8) :: f(:)
      real(8) :: dh,a,b
      real(8) :: sum
      integer :: L
      L=size(f)
      dh=(b-a)/dble(L-1)
      sum = d_simpson_dh_sample(f,dh)
   end function d_simpson_ab_sample
   function c_simpson_ab_sample(f,a,b) result(sum)
      complex(8) :: f(:)
      real(8)    :: dh,a,b
      complex(8) :: sum
      integer    :: L
      L=size(f)
      dh=(b-a)/dble(L-1)
      sum = c_simpson_dh_sample(f,dh)
   end function c_simpson_ab_sample


   function d_simpson_dh_sample(f,dh) result(sum)
      real(8) :: f(:)
      integer :: n
      real(8) :: dh,sum,sum1,sum2,int1,int2
      integer :: i,p,m,mm,mmm
      N=size(f)
      sum=0d0
      if(N==1)return
      sum1=0d0
      sum2=0d0
      int1=0d0
      int2=0d0
      if(mod(n-1,2)==0)then                !if n-1 is even:
         do i=2,N-1,2
            sum1 = sum1 + f(i)
         enddo
         do i=3,N-2,2
            sum2 = sum2 + f(i)
         enddo
         sum = (f(1) + 4d0*sum1 + 2d0*sum2 + f(n))*dh/3d0
      else                        !if n-1 is odd, use Simpson's for N-3 slices + 3/8rule for the last
         if (N>=6) then
            do i=2,N-4,2
               sum1 = sum1 + f(i)
            enddo
            do i=3,N-5,2
               sum2 = sum2 + f(i)
            enddo
            int1 = (f(1) + 4d0*sum1 + 2d0*sum2 + f(n-3))*dh/3d0
         endif
         int2 = (f(n-3)+3d0*f(n-2)+3d0*f(n-1)+f(n))*dh*3d0/8d0
         sum  = int1 + int2
      end if
   end function d_simpson_dh_sample
   function c_simpson_dh_sample(f,dh) result(sum)
      complex(8)           :: f(:)
      integer              :: n
      real(8)              :: dh
      complex(8)           :: sum,sum1,sum2,int1,int2
      integer              :: i,p,m
      complex(8),parameter :: zero=cmplx(0d0,0d0,8)
      N=size(f)
      sum=zero
      if(N==1)return
      sum1=zero
      sum2=zero
      int1=zero
      int2=zero
      if(mod(n-1,2)==0)then                !if n-1 is even:
         do i=2,N-1,2
            sum1 = sum1 + f(i)
         enddo
         do i=3,N-2,2
            sum2 = sum2 + f(i)
         enddo
         sum = (f(1) + 4d0*sum1 + 2d0*sum2 + f(n))*dh/3d0
      else                        !if n-1 is odd, use Simpson's for N-3 slices + 3/8rule for the last
         if (N>=6) then
            do i=2,N-4,2
               sum1 = sum1 + f(i)
            enddo
            do i=3,N-5,2
               sum2 = sum2 + f(i)
            enddo
            int1 = (f(1) + 4d0*sum1 + 2d0*sum2 + f(n-3))*dh/3d0
         endif
         int2 = (f(n-3)+3d0*f(n-2)+3d0*f(n-1)+f(n))*dh*3d0/8d0
         sum  = int1 + int2
      end if
   end function c_simpson_dh_sample


   function d_simpson_nonlin_sample(f,x) result(sum)
      real(8) :: f(:)
      real(8) :: x(size(f))
      real(8) :: dx(size(f)+1)
      real(8) :: sum,sum1,sum2,sum3
      real(8) :: a,b,dh
      integer :: i,n,m
      n=size(f)
      m=n-1
      dx=0d0
      forall(i=1:m)dx(i)=x(i+1)-x(i)
      sum1=0d0
      sum2=0d0
      sum3=0d0
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
   function c_simpson_nonlin_sample(f,x) result(sum)
      complex(8) :: f(:)
      complex(8) :: sum
      real(8)    :: x(size(f)),rsum,isum
      rsum=d_simpson_nonlin_sample(dreal(f),x)
      isum=d_simpson_nonlin_sample(dimag(f),x)
      sum  = dcmplx(rsum,isum)
   end function c_simpson_nonlin_sample


   !PURPOSE:  working procedures:
   function d_trapz2d_func(func,xrange,yrange,Nx,Ny) result(int)
      interface
         function func(x)
            real(8),dimension(:) :: x
            real(8)              :: func
         end function func
      end interface
      real(8),dimension(2) :: xrange,yrange
      integer              :: Nx,Ny,i,j
      real(8)              :: xx(Nx),yy(Ny)
      real(8)              :: hx,hy
      real(8)              :: int
      hx=xrange(2)-xrange(1)
      hx=hx/Nx
      hy=yrange(2)-yrange(1)
      hy=hy/Ny
      int=&
         func([xrange(1),yrange(1)])+&
         func([xrange(1),yrange(2)])+&
         func([xrange(2),yrange(1)])+&
         func([xrange(2),yrange(2)])
      xx=linspace(xrange(1),xrange(2),Nx)
      yy=linspace(yrange(1),yrange(2),Ny)
      do i=2,Nx
         do j=2,Ny
            int=int+4d0*func([xx(i),yy(j)])
         enddo
      enddo
      do j=2,Ny
         int=int+2d0*( func([xrange(1),yy(j)]) + func([xrange(2),yy(j)]) )
      enddo
      do i=2,Nx
         int=int+2d0*( func([xx(i),yrange(1)]) + func([xx(i),yrange(2)]) )
      enddo
      int=int*hx*hy/4d0
   end function d_trapz2d_func

   function c_trapz2d_func(func,xrange,yrange,Nx,Ny) result(int)
      interface
         function func(x)
            real(8),dimension(:) :: x
            complex(8)           :: func
         end function func
      end interface
      real(8),dimension(2) :: xrange,yrange
      integer              :: Nx,Ny,i,j
      real(8)              :: xx(Nx),yy(Ny)
      real(8)              :: hx,hy
      complex(8)           :: int
      hx=xrange(2)-xrange(1)
      hx=hx/Nx
      hy=yrange(2)-yrange(1)
      hy=hy/Ny
      int=&
         func([xrange(1),yrange(1)])+&
         func([xrange(1),yrange(2)])+&
         func([xrange(2),yrange(1)])+&
         func([xrange(2),yrange(2)])
      xx=linspace(xrange(1),xrange(2),Nx)
      yy=linspace(yrange(1),yrange(2),Ny)
      do i=2,Nx
         do j=2,Ny
            int=int+4d0*func([xx(i),yy(j)])
         enddo
      enddo
      do j=2,Ny
         int=int+2d0*( func([xrange(1),yy(j)]) + func([xrange(2),yy(j)]) )
      enddo
      do i=2,Nx
         int=int+2d0*( func([xx(i),yrange(1)]) + func([xx(i),yrange(2)]) )
      enddo
      int=int*hx*hy/4d0
   end function c_trapz2d_func

   !RECURSIVE VERSIONS:
   function d_trapz2d_func_recursive(func,xrange,yrange,N0,iterative,threshold) result(int)
      interface
         function func(x)
            real(8),dimension(:) :: x
            real(8)              :: func
         end function func
      end interface
      real(8),dimension(2)      :: xrange,yrange
      integer                   :: N,icount
      real(8)                   :: int,eps,int0
      integer,optional          :: N0
      integer                   :: N0_
      logical,optional          :: iterative
      logical                   :: iterative_
      real(8),optional          :: threshold
      real(8)                   :: threshold_
      iterative_=.false.;if(present(iterative))iterative_=iterative
      N0_=51;if(present(N0))N0_=N0
      threshold_=1d0;if(iterative_)threshold_=5.d-3
      if(present(threshold))threshold_=threshold
      N=N0_
      eps=1d0
      icount=1
      int=d_trapz2d_func(func,xrange,yrange,N,N)
      do while (eps>threshold_)
         icount=icount+1
         int0=int
         N=2*N-10
         int=d_trapz2d_func(func,xrange,yrange,N,N)
         eps=abs(int-int0)/abs(int)
      enddo
   end function d_trapz2d_func_recursive

   function c_trapz2d_func_recursive(func,xrange,yrange,N0,iterative,threshold) result(int)
      interface
         function func(x)
            real(8),dimension(:) :: x
            complex(8)           :: func
         end function func
      end interface
      real(8),dimension(2)      :: xrange,yrange
      integer                   :: N,icount
      real(8)                   :: eps,int0
      integer,optional          :: N0
      integer                   :: N0_
      logical,optional          :: iterative
      logical                   :: iterative_
      real(8),optional          :: threshold
      real(8)                   :: threshold_
      complex(8)                :: int
      iterative_=.false.;if(present(iterative))iterative_=iterative
      N0_=51;if(present(N0))N0_=N0
      threshold_=1d0;if(iterative_)threshold_=5.d-3
      if(present(threshold))threshold_=threshold
      N=N0_
      eps=1d0
      icount=1
      int=c_trapz2d_func(func,xrange,yrange,N,N)
      do while (eps>threshold_)
         icount=icount+1
         int0=int
         N=2*N-10
         int=c_trapz2d_func(func,xrange,yrange,N,N)
         eps=abs(int-int0)/abs(int)
      enddo
   end function c_trapz2d_func_recursive

   function d_simps2d_func(func,xrange,yrange,Nx,Ny) result(int)
      interface
         function func(x)
            real(8),dimension(:) :: x
            real(8)              :: func
         end function func
      end interface
      real(8),dimension(2) :: xrange,yrange
      integer :: Nx,Ny,i,j
      real(8) :: xx(2*Nx),yy(2*Ny)
      real(8) :: hx,hy
      real(8) :: int
      hx=xrange(2)-xrange(1)
      hx=hx/Nx/2
      hy=yrange(2)-yrange(1)
      hy=hy/Ny/2
      xx=linspace(xrange(1),xrange(2),2*Nx)
      yy=linspace(yrange(1),yrange(2),2*Ny)
      !
      int=&
         func([xrange(1),yrange(1)])+&
         func([xrange(1),yrange(2)])+&
         func([xrange(2),yrange(1)])+&
         func([xrange(2),yrange(2)])
      !
      do j=1,Ny
         int=int+4d0*(func([xrange(1),yy(2*j-1)])+func([xrange(2),yy(2*j-1)]))
      enddo
      do j=1,Ny-1
         int=int+2d0*(func([xrange(1),yy(2*j)])+func([xrange(2),yy(2*j)]))
      enddo
      !
      do i=1,Nx
         int=int+4d0*(func([xx(2*i-1),yrange(1)])+func([xx(2*i-1),yrange(2)]))
      enddo
      do i=1,Nx-1
         int=int+2d0*(func([xx(2*i),yrange(1)])+func([xx(2*i),yrange(2)]))
      enddo
      !
      do j=1,Ny
         do i=1,Nx
            int=int+16d0*func([xx(2*i-1),yy(2*j-1)])
         enddo
      enddo
      do j=1,Ny-1
         do i=1,Nx
            int=int+8d0*func([xx(2*i-1),yy(2*j)])
         enddo
      enddo
      !
      do j=1,Ny
         do i=1,Nx-1
            int=int+8d0*func([xx(2*i),yy(2*j-1)])
         enddo
      enddo
      do j=1,Ny-1
         do i=1,Nx-1
            int=int+4d0*func([xx(2*i),yy(2*j)])
         enddo
      enddo
      int=int*hx*hy/9d0
   end function d_simps2d_func

   function c_simps2d_func(func,xrange,yrange,Nx,Ny) result(int)
      interface
         function func(x)
            real(8),dimension(:) :: x
            complex(8)           :: func
         end function func
      end interface
      real(8),dimension(2) :: xrange,yrange
      integer              :: Nx,Ny,i,j
      real(8)              :: xx(2*Nx),yy(2*Ny)
      real(8)              :: hx,hy
      complex(8)           :: int
      hx=xrange(2)-xrange(1)
      hx=hx/Nx/2
      hy=yrange(2)-yrange(1)
      hy=hy/Ny/2
      xx=linspace(xrange(1),xrange(2),2*Nx)
      yy=linspace(yrange(1),yrange(2),2*Ny)
      !
      int=&
         func([xrange(1),yrange(1)])+&
         func([xrange(1),yrange(2)])+&
         func([xrange(2),yrange(1)])+&
         func([xrange(2),yrange(2)])
      !
      do j=1,Ny
         int=int+4d0*(func([xrange(1),yy(2*j-1)])+func([xrange(2),yy(2*j-1)]))
      enddo
      do j=1,Ny-1
         int=int+2d0*(func([xrange(1),yy(2*j)])+func([xrange(2),yy(2*j)]))
      enddo
      !
      do i=1,Nx
         int=int+4d0*(func([xx(2*i-1),yrange(1)])+func([xx(2*i-1),yrange(2)]))
      enddo
      do i=1,Nx-1
         int=int+2d0*(func([xx(2*i),yrange(1)])+func([xx(2*i),yrange(2)]))
      enddo
      !
      do j=1,Ny
         do i=1,Nx
            int=int+16d0*func([xx(2*i-1),yy(2*j-1)])
         enddo
      enddo
      do j=1,Ny-1
         do i=1,Nx
            int=int+8d0*func([xx(2*i-1),yy(2*j)])
         enddo
      enddo
      !
      do j=1,Ny
         do i=1,Nx-1
            int=int+8d0*func([xx(2*i),yy(2*j-1)])
         enddo
      enddo
      do j=1,Ny-1
         do i=1,Nx-1
            int=int+4d0*func([xx(2*i),yy(2*j)])
         enddo
      enddo
      int=int*hx*hy/9d0
   end function c_simps2d_func

   !RECURSIVE VERSIONS:
   function d_simps2d_func_recursive(func,xrange,yrange,N0,iterative,threshold) result(int)
      interface
         function func(x)
            real(8),dimension(:) :: x
            real(8)              :: func
         end function func
      end interface
      real(8),dimension(2)      :: xrange,yrange
      integer                   :: N,icount
      real(8)                   :: eps,int0
      integer,optional          :: N0
      integer                   :: N0_
      logical,optional          :: iterative
      logical                   :: iterative_
      real(8),optional          :: threshold
      real(8)                   :: threshold_
      real(8)                   :: int
      iterative_=.false.;if(present(iterative))iterative_=iterative
      N0_=51;if(present(N0))N0_=N0
      threshold_=1d0;if(iterative_)threshold_=1.d-3
      if(present(threshold))threshold_=threshold
      N=N0_
      eps=1d0
      icount=1
      int=d_simps2d_func(func,xrange,yrange,N,N)
      do while (eps>threshold_)
         icount=icount+1
         int0=int
         N=2*N-10
         int=d_simps2d_func(func,xrange,yrange,N,N)
         eps=abs(int-int0)/abs(int)
      enddo
   end function d_simps2d_func_recursive

   function c_simps2d_func_recursive(func,xrange,yrange,N0,iterative,threshold) result(int)
      interface
         function func(x)
            real(8),dimension(:) :: x
            complex(8)           :: func
         end function func
      end interface
      real(8),dimension(2)      :: xrange,yrange
      integer                   :: N,icount
      real(8)                   :: eps,int0
      integer,optional          :: N0
      integer                   :: N0_
      logical,optional          :: iterative
      logical                   :: iterative_
      real(8),optional          :: threshold
      real(8)                   :: threshold_
      complex(8)                :: int
      iterative_=.false.;if(present(iterative))iterative_=iterative
      N0_=51;if(present(N0))N0_=N0
      threshold_=1d0;if(iterative_)threshold_=1.d-3
      if(present(threshold))threshold_=threshold
      N=N0_
      eps=1d0
      icount=1
      int=c_simps2d_func(func,xrange,yrange,N,N)
      do while (eps>threshold_)
         icount=icount+1
         int0=int
         N=2*N-10
         int=c_simps2d_func(func,xrange,yrange,N,N)
         eps=abs(int-int0)/abs(int)
      enddo
   end function c_simps2d_func_recursive

   !+-----------------------------------------------------------------+
   !PURPOSE: Trapezoidal rule for 2d data function integration
   ! + _std: integrate f(:) between a and b. a,b default are 0,1
   !+-----------------------------------------------------------------+
   function d_trapz2d_sample(func,dhx,dhy,xrange,yrange) result(int)
      real(8),dimension(:,:)        :: func
      real(8),optional              :: dhx,dhy
      real(8),dimension(2),optional :: xrange,yrange
      real(8)                       :: hx,hy
      real(8),dimension(2)          :: xrange_,yrange_
      integer                       :: Nx,Ny,i,j
      real(8)                       :: int
      Nx=size(func,1)
      Ny=size(func,2)
      xrange_=[0d0,1d0];if(present(xrange))xrange_=xrange
      yrange_=[0d0,1d0];if(present(yrange))yrange_=yrange
      hx=xrange_(2)-xrange_(1);hx=hx/Nx
      if(present(dhx))hx=dhx
      hy=yrange_(2)-yrange_(1);hy=hy/Ny
      if(present(dhy))hy=dhy
      !
      int = func(1,1)+func(1,Ny)+func(Nx,1)+func(Nx,Ny)
      do i=2,Nx
         do j=2,Ny
            int=int+4d0*func(i,j)
         enddo
      enddo
      do j=2,Ny
         int=int+2d0*( func(1,j) + func(Nx,j) )
      enddo
      do i=2,Nx
         int=int+2d0*( func(i,1) + func(i,Ny) )
      enddo
      int=int*hx*hy/4d0
   end function d_trapz2d_sample

   function c_trapz2d_sample(func,dhx,dhy,xrange,yrange) result(int)
      complex(8),dimension(:,:)     :: func
      real(8),optional              :: dhx,dhy
      real(8),dimension(2),optional :: xrange,yrange
      real(8)                       :: hx,hy
      real(8),dimension(2)          :: xrange_,yrange_
      integer                       :: Nx,Ny,i,j
      complex(8)                    :: int
      Nx=size(func,1)
      Ny=size(func,2)
      xrange_=[0d0,1d0];if(present(xrange))xrange_=xrange
      yrange_=[0d0,1d0];if(present(yrange))yrange_=yrange
      hx=xrange_(2)-xrange_(1);hx=hx/Nx
      if(present(dhx))hx=dhx
      hy=yrange_(2)-yrange_(1);hy=hy/Ny
      if(present(dhy))hy=dhy
      !
      int = func(1,1)+func(1,Ny)+func(Nx,1)+func(Nx,Ny)
      do i=2,Nx
         do j=2,Ny
            int=int+4d0*func(i,j)
         enddo
      enddo
      do j=2,Ny
         int=int+2d0*( func(1,j) + func(Nx,j) )
      enddo
      do i=2,Nx
         int=int+2d0*( func(i,1) + func(i,Ny) )
      enddo
      int=int*hx*hy/4d0
   end function c_trapz2d_sample

   !+-----------------------------------------------------------------+
   !PURPOSE: Simpson rule for 2d data function integration
   !+-----------------------------------------------------------------+
   function d_simps2d_sample(func,dhx,dhy,xrange,yrange) result(int)
      real(8),dimension(:,:)        :: func
      real(8),optional              :: dhx,dhy
      real(8),dimension(2),optional :: xrange,yrange
      real(8)                       :: hx,hy
      real(8),dimension(2)          :: xrange_,yrange_
      integer                       :: Nx,Ny,i,j
      real(8)                       :: int
      Nx=size(func,1)
      Ny=size(func,2)
      xrange_=[0d0,1d0];if(present(xrange))xrange_=xrange
      yrange_=[0d0,1d0];if(present(yrange))yrange_=yrange
      hx=xrange_(2)-xrange_(1);hx=hx/Nx/2
      if(present(dhx))hx=dhx
      hy=yrange_(2)-yrange_(1);hy=hy/Ny/2
      if(present(dhy))hy=dhy
      !
      int = func(1,1)+func(1,Ny)+func(Nx,1)+func(Nx,Ny)
      !
      do j=1,Ny
         int=int+4d0*( func(1,2*j-1) + func(Nx,2*j-1) )
      enddo
      do j=1,Ny-1
         int=int+2d0*( func(1,2*j) + func(Nx,2*j) )
      enddo
      !
      do i=1,Nx
         int=int+4d0*( func(2*i-1,1) + func(2*i-1,Ny) )
      enddo
      do i=1,Nx-1
         int=int+2d0*( func(2*i,1) + func(2*i,Ny) )
      enddo
      !
      do j=1,Ny
         do i=1,Nx
            int=int+16d0*func(2*i-1,2*j-1)
         enddo
      enddo
      do j=1,Ny-1
         do i=1,Nx
            int=int+8d0*func(2*i-1,2*j)
         enddo
      enddo
      !
      do j=1,Ny
         do i=1,Nx-1
            int=int+8d0*func(2*i,2*j-1)
         enddo
      enddo
      do j=1,Ny-1
         do i=1,Nx-1
            int=int+4d0*func(2*i,2*j)
         enddo
      enddo
      int=int*hx*hy/9d0
   end function d_simps2d_sample

   function c_simps2d_sample(func,dhx,dhy,xrange,yrange) result(int)
      complex(8),dimension(:,:)     :: func
      real(8),optional              :: dhx,dhy
      real(8),dimension(2),optional :: xrange,yrange
      real(8)                       :: hx,hy
      real(8),dimension(2)          :: xrange_,yrange_
      integer                       :: Nx,Ny,i,j
      complex(8)                    :: int
      Nx=size(func,1)
      Ny=size(func,2)
      xrange_=[0d0,1d0];if(present(xrange))xrange_=xrange
      yrange_=[0d0,1d0];if(present(yrange))yrange_=yrange
      hx=xrange_(2)-xrange_(1);hx=hx/Nx/2
      if(present(dhx))hx=dhx
      hy=yrange_(2)-yrange_(1);hy=hy/Ny/2
      if(present(dhy))hy=dhy
      !
      int = func(1,1)+func(1,Ny)+func(Nx,1)+func(Nx,Ny)
      !
      do j=1,Ny
         int=int+4d0*( func(1,2*j-1) + func(Nx,2*j-1) )
      enddo
      do j=1,Ny-1
         int=int+2d0*( func(1,2*j) + func(Nx,2*j) )
      enddo
      !
      do i=1,Nx
         int=int+4d0*( func(2*i-1,1) + func(2*i-1,Ny) )
      enddo
      do i=1,Nx-1
         int=int+2d0*( func(2*i,1) + func(2*i,Ny) )
      enddo
      !
      do j=1,Ny
         do i=1,Nx
            int=int+16d0*func(2*i-1,2*j-1)
         enddo
      enddo
      do j=1,Ny-1
         do i=1,Nx
            int=int+8d0*func(2*i-1,2*j)
         enddo
      enddo
      !
      do j=1,Ny
         do i=1,Nx-1
            int=int+8d0*func(2*i,2*j-1)
         enddo
      enddo
      do j=1,Ny-1
         do i=1,Nx-1
            int=int+4d0*func(2*i,2*j)
         enddo
      enddo
      int=int*hx*hy/9d0
   end function c_simps2d_sample





   !+-----------------------------------------------------------------------------+!
   !PURPOSE:
   !generic interface to QUADPACK routines:
   !non-automatic integrator:
   !estimate the integral on [a,b] using 15, 21,..., 61
   !point rule and return an error estimate.
   ! public :: QNG
   ! public :: QK15
   ! public :: QK21
   ! public :: QK31
   ! public :: QK41
   ! public :: QK51
   ! public :: QK61
   !automatic integrator:
   !without weights:
   ! public :: QAG  !globally adaptive integrator
   ! public :: QAGI !integration over infinite intervals
   ! public :: QAGS !globally adaptive interval subdivisio with epsilon extrapolation
   ! public :: QAGP !serves the same purposes as qags, but user-supplied information about singularities
   !with weights:
   ! public :: QAWO !integration of cos(omega*x)*f(x) sin(omega*x)*f(x) over a finite interval (a,b).
   ! public :: QAWF !Fourier cosine or fourier sine transform of f(x)
   ! public :: QAWC !computes the cauchy principal value of f(x)/(x-c) over a finite interval (a,b)
   ! public :: QAWS !integrates w(x)*f(x) over (a,b) with a < b finite,
   ! !              ! and   w(x) = ((x-a)**alfa)*((b-x)**beta)*v(x)
   ! !              ! where v(x) = 1 or log(x-a) or log(b-x)
   ! !              !              or log(x-a)*log(b-x)
   ! !              ! and   -1 < alfa, -1 < beta.
   !+-----------------------------------------------------------------------------+!
   subroutine quad_func(func,a,b,epsabs,epsrel,&
      key,&
      inf,&
      singular_endpoint,&
      singular_points,&
      cpole,&
      alfa,beta,&
      omega,&
      weight_func,&
      verbose,&
      strict,&
      result)
      interface
         function func(x)
            real(8) :: x
            real(8) :: func
         end function func
      end interface
      !optional variables
      real(8)                       :: a
      real(8),optional              :: b
      real(8),optional              :: epsabs
      real(8),optional              :: epsrel
      integer,optional              :: key               !order of GK rule in QAG
      integer,optional              :: inf               !infinite integration limit in QAGI:
      !                                                  ! -1(-inf:a);1(a:inf);2(-inf:inf)
      logical,optional              :: singular_endpoint !T if singular_endpoint exists (QAGS)
      real(8),dimension(:),optional :: singular_points   !location of singular points in QAGP
      real(8),optional              :: cpole             !location of the pole in QAWC f(x)/x-cpole
      real(8),optional              :: alfa,beta         !parameters for QAWS
      real(8),optional              :: omega             !frequency of the weight funcion in QAWF,QAWO
      integer,optional              :: weight_func       !if(QAWF,QAWO)then
      !                                                  !  weight_func=1,2 -> cos(omega*x),sin(omega*x)
      !                                                  !if(QAWS)then
      !                                                  !  weight_func = 1  (x-a)**alfa*(b-x)**beta
      !                                                  !  weight_func = 2  (x-a)**alfa*(b-x)**beta*log(x-a)
      !                                                  !  weight_func = 3  (x-a)**alfa*(b-x)**beta*log(b-x)
      !                                                  !  weight_func = 4  (x-a)**alfa*(b-x)**beta*log(x-a)*log(b-x)
      logical,optional              :: verbose
      logical,optional              :: strict
      real(8)                       :: result
      !actual default variables
      real(8)                       :: epsabs_
      real(8)                       :: epsrel_
      logical                       :: verbose_
      logical                       :: strict_
      real(8)                       :: abserr
      integer                       :: Neval
      integer                       :: Ier,i
      integer                       :: Nsingular_points
      character(len=20)             :: routine_name
      !
      !
      !
      epsabs_ = 1d-12 ; if(present(epsabs))epsabs_=epsabs
      epsrel_ = 1d-6  ; if(present(epsrel))epsrel_=epsrel
      verbose_=.false.; if(present(verbose))verbose_=verbose
      strict_ =.true. ; if(present(strict))strict_=strict
      !
      if(present(key).AND.(present(alfa).OR.present(beta)))stop "ERROR in quad: key & alfa,beta"
      if(present(key).AND.present(cpole))stop "ERROR in quad: key & cpole"
      if(present(key).AND.present(singular_points))stop "ERROR in quad: key & singular_points"
      if(present(key).AND.present(singular_endpoint))stop "ERROR in quad: key & singular_endpoint"
      if(present(key).AND.present(omega))stop "ERROR in quad: key & omega"
      if(present(key).AND.present(inf))stop "ERROR in quad: key & inf"
      if(present(singular_endpoint).AND.(present(alfa).OR.present(beta)))stop "ERROR in quad: singular_endpoint & alfa,beta"
      if(present(singular_endpoint).AND.present(cpole))stop "ERROR in quad: singular_endpoint & cpole"
      if(present(singular_endpoint).AND.present(singular_points))stop "ERROR in quad: singular_endpoint & singular_points"
      if(present(singular_endpoint).AND.present(omega))stop "ERROR in quad: singular_endpoint & omega"
      if(present(singular_endpoint).AND.present(inf))stop "ERROR in quad: singular_endpoint & inf"
      if(present(cpole).AND.(present(alfa).OR.present(beta)))stop "ERROR in quad: cpole & alfa,beta"
      if(present(cpole).AND.present(singular_points))stop "ERROR in quad: cpole & singular_points"
      if(present(cpole).AND.present(omega))stop "ERROR in quad: cpole & omega"
      if(present(cpole).AND.present(inf))stop "ERROR in quad: cpole & inf"
      if(present(omega).AND.(present(alfa).OR.present(beta)))stop "ERROR in quad: omega & alfa,beta"
      if(present(omega).AND.present(singular_points))stop "ERROR in quad: omega & singular_points"
      if(present(omega).AND.present(inf))stop "ERROR in quad: omega & inf"
      if(present(singular_points).AND.(present(alfa).OR.present(beta)))stop "ERROR in quad: singular_points & alfa,beta"
      if(present(singular_points).AND.present(inf))stop "ERROR in quad: singular_points & inf"
      if(present(inf).AND.(present(alfa).OR.present(beta)))stop "ERROR in quad: inf & alfa,beta"
      !
      if(present(b))then
         !QNG>QAGS:  no need
         !QAGS: need singular_endpoint=T
         !QAG : need key=1,6
         !QAGP: need singular_points
         !QAWC: need cpole
         !QAWO: need omega + weight_func=1,2
         !QAWS: need alfa + beta + weight_func=1,2,3,4
         routine_name='QNG'
         if(present(singular_endpoint))routine_name='QAGS'
         if(present(key))routine_name='QAG'
         if(present(singular_points))routine_name='QAGP'
         if(present(cpole))routine_name='QAWC'
         if(present(omega))routine_name='QAWO'
         if(present(alfa).AND.present(beta))routine_name='QAWS'
      else
         !QAGI: need inf
         !QAWF: need omega + weight_func=1,2
         if(present(inf))then
            routine_name='QAGI'
         elseif(present(omega))then
            routine_name='QAWF'
         else
            stop 'ERROR in quad: b not present but neither inf (QAGI) nor omega (QAWF) are given. stop.'
         endif
      endif
      !
      if(verbose_)write(*,"(A,A)")"QUAD: selected procedure =", routine_name
      !
      result=0d0
      !
      select case(routine_name)
       case ('QNG')
         ! call QNG(func,a,b,epsabs_,epsrel_,result,abserr,neval,ier)
         ! call QAG(func,a,b,epsabs_,epsrel_,1,result,abserr,neval,ier)
         call QAGS(func,a,b,epsabs_,epsrel_,result,abserr,neval,ier)
         if(verbose_)then
            write(*,'(A,2F14.6)')'Endpoints (a,b)                             =', a,b
            write(*,'(A,F14.6)') 'Estimated integral                          =', result
            write(*,'(A,F14.6)') 'Estimated integral error                    =', abserr
            write(*,'(A,I8)')    'Error return code IER                       =', ier
         else
            if(ier>0)then
               write(*,'(A,I8)') 'Error return code IER =', ier
               if(strict_)stop
            endif
         endif
         !
         !
       case ('QAGS')
         call QAGS(func,a,b,epsabs_,epsrel_,result,abserr,neval,ier)
         if(verbose_)then
            write(*,'(A,2F14.6)')'(Singular) endpoints (a*,b*)                =', a,b
            write(*,'(A,F14.6)') 'Estimated integral                          =', result
            write(*,'(A,F14.6)') 'Estimated integral error                    =', abserr
            write(*,'(A,I8)')    'Error return code IER                       =', ier
         else
            if(ier>2)then
               write(*,'(A,I8)') 'Error return code IER =', ier
               if(strict_)stop
            endif
         endif
         !
         !
       case ('QAG')
         if(key<1.OR.key>6)stop "ERROR in quad: use QAG, key !=[1:6]"
         call QAG(func,a,b,epsabs_,epsrel_,key,result,abserr,neval,ier)
         if(verbose_)then
            write(*,'(A,2F14.6)')'Endpoints (a,b)                             =', a,b
            write(*,'(A,I8)')    'Approximation or                            =', key
            write(*,'(A,F14.6)') 'Estimated integral                          =', result
            write(*,'(A,F14.6)') 'Estimated integral error                    =', abserr
            write(*,'(A,I8)')    'Error return code IER                       =', ier
         else
            if(ier>2)then
               write(*,'(A,I8)') 'Error return code IER =', ier
               if(strict_)stop
            endif
         endif
         !
         !
       case ('QAGP')
         Nsingular_points=size(singular_points)
         call QAGP(func,a,b,Nsingular_points,singular_points,epsabs_,epsrel_,result,abserr,neval,ier)
         if(verbose_)then
            write(*,'(A,2F14.6)')'Endpoints (a,b)                             =', a,b
            write(*,'(A,I8)')'N_singular_points                           =', Nsingular_points
            do i=1,Nsingular_points
               write(*,"(A,I6,F14.6)")'I_singular_point                            =',i,singular_points(i)
            enddo
            write(*,'(A,F14.6)') 'Estimated integral                          =', result
            write(*,'(A,F14.6)') 'Estimated integral error                    =', abserr
            write(*,'(A,I8)')    'Error return code IER                       =', ier
         else
            if(ier>2)then
               write(*,'(A,I8)') 'Error return code IER =', ier
               if(strict_)stop
            endif
         endif
         !
         !
       case ('QAWC')
         call QAWC(func,a,b,cpole,epsabs_,epsrel_,result,abserr,neval,ier)
         if(verbose_)then
            write(*,'(A,2F14.6)')'Endpoints (a,b)                             =', a,b
            write(*,'(A,F14.6)')'C_pole                                      =', cpole
            write(*,'(A,F14.6)') 'Estimated integral                          =', result
            write(*,'(A,F14.6)') 'Estimated integral error                    =', abserr
            write(*,'(A,I8)')    'Error return code IER                       =', ier
         else
            if(ier>2)then
               write(*,'(A,I8)') 'Error return code IER =', ier
               if(strict_)stop
            endif
         endif
         !
         !
       case ('QAWO')
         if(weight_func<1.OR.weight_func>2)stop "ERROR in quad: use QAWO, weight_func !=[1,2]"
         if(verbose_)then
            if(weight_func==1)then
               write(*,'(A)')'W(x) function is cos(omega*x)'
            elseif(weight_func==2)then
               write(*,'(A)')'W(x) function is sin(omega*x)'
            endif
         endif
         call QAWO(func,a,b,omega,weight_func,epsabs_,epsrel_,result,abserr,neval,ier)
         if(verbose_)then
            write(*,'(A,2F14.6)')'Endpoints (a,b)                             =', a
            write(*,'(A,F14.6)') 'Omega                                       =', omega
            write(*,'(A,F14.6)') 'Estimated integral                          =', result
            write(*,'(A,F14.6)') 'Estimated integral error                    =', abserr
            write(*,'(A,I8)')    'Error return code IER                       =', ier
         else
            if(ier>2)then
               write(*,'(A,I8)') 'Error return code IER =', ier
               if(strict_)stop
            endif
         endif
         !
         !
       case ('QAWS')
         if(weight_func<1.OR.weight_func>4)stop "ERROR in quad: use QAWS, weight_func !=[1,..,4]"
         if(verbose_)then
            if(weight_func==1)then
               write(*,'(A)')'W(x) function is (x-a)**alfa*(b-x)**beta         '
            elseif(weight_func==2)then
               write(*,'(A)')'W(x) function is (x-a)**alfa*(b-x)**beta*log(x-a)'
            elseif(weight_func==3)then
               write(*,'(A)')'W(x) function is (x-a)**alfa*(b-x)**beta*log(b-x) '
            elseif(weight_func==4)then
               write(*,'(A)')'W(x) function is (x-a)**alfa*(b-x)**beta*log(x-a)*log(b-x)'
            endif
         endif
         call QAWS(func,a,b,alfa,beta,weight_func,epsabs_,epsrel_,result,abserr,neval,ier)
         if(verbose_)then
            write(*,'(A,2F14.6)')'Endpoints (a,b)                             =', a,b
            write(*,'(A,2F14.6)')'Parameters (alfa,beta)                      =', alfa,beta
            write(*,'(A,F14.6)') 'Estimated integral                          =', result
            write(*,'(A,F14.6)') 'Estimated integral error                    =', abserr
            write(*,'(A,I8)')    'Error return code IER                       =', ier
         else
            if(ier>2)then
               write(*,'(A,I8)') 'Error return code IER =', ier
               if(strict_)stop
            endif
         endif
         !
         !
       case ('QAGI')
         if(inf<-1.OR.inf>2)stop "ERROR in quad: use QAGI, inf !=[-1,1,2]"
         call QAGI(func,a,inf,epsabs_,epsrel_,result,abserr,neval,ier)
         if(verbose_)then
            if(inf==-1)then
               write(*,'(A,F14.6)')'Integration interval is (-infty:a)',a
            elseif(inf==1)then
               write(*,'(A,F14.6)')'Integration interval is (a:infty) ',a
            elseif(inf==2)then
               write(*,'(A)')'Integration interval is (-infty:infty)'
            endif
            write(*,'(A,F14.6)') 'Estimated integral                          =', result
            write(*,'(A,F14.6)') 'Estimated integral error                    =', abserr
            write(*,'(A,I8)')    'Error return code IER                       =', ier
         else
            if(ier>2)then
               write(*,'(A,I8)') 'Error return code IER =', ier
               if(strict_)stop
            endif
         endif
         !
         !
       case ('QAWF')
         if(weight_func<1.OR.weight_func>2)stop "ERROR in quad: use QAWF, weight_func !=[1,2]"
         if(verbose_)then
            if(weight_func==1)then
               write(*,'(A)')'W(x) function is cos(omega*x)'
            elseif(weight_func==2)then
               write(*,'(A)')'W(x) function is sin(omega*x)'
            endif
         endif
         call QAWF(func,a,omega,weight_func,epsabs_,result,abserr,neval,ier)
         if(verbose_)then
            write(*,'(A,2F14.6)')'Endpoints (a,infty)                         =', a
            write(*,'(A,F14.6)') 'Omega                                       =', omega
            write(*,'(A,F14.6)') 'Estimated integral                          =', result
            write(*,'(A,F14.6)') 'Estimated integral error                    =', abserr
            write(*,'(A,I8)')    'Error return code IER                       =', ier
         else
            if(ier>2)then
               write(*,'(A,I8)') 'Error return code IER =', ier
               if(strict_)stop
            endif
         endif
         !
         !
      end select
   end subroutine quad_func

   subroutine quad_sample(fsample,a,b,&
      epsabs,epsrel,&
      key,&
      singular_endpoint,&
      singular_points,&
      cpole,&
      alfa,beta,&
      omega,&
      weight_func,&
      verbose,&
      Ninterp,&
      result)
      real(8),dimension(:)             :: fsample
      real(8)                          :: a
      real(8)                          :: b
      !optional variables
      real(8),optional                 :: epsabs
      real(8),optional                 :: epsrel
      integer,optional                 :: key               !order of GK rule in QAG
      logical,optional                 :: singular_endpoint !T if singular_endpoint exists (QAGS)
      real(8),dimension(:),optional    :: singular_points   !location of singular points in QAGP
      real(8),optional                 :: cpole             !location of the pole in QAWC f(x)/x-cpole
      real(8),optional                 :: alfa,beta         !parameters for QAWS
      real(8),optional                 :: omega             !frequency of the weight funcion in QAWO
      integer,optional                 :: weight_func       !if(QAWF)then
      !                                                  !  weight_func=1,2 -> cos(omega*x),sin(omega*x)
      !                                                  !if(QAWS)then
      !                                                  !  weight_func = 1  (x-a)**alfa*(b-x)**beta
      !                                                  !  weight_func = 2  (x-a)**alfa*(b-x)**beta*log(x-a)
      !                                                  !  weight_func = 3  (x-a)**alfa*(b-x)**beta*log(b-x)
      !                                                  !  weight_func = 4  (x-a)**alfa*(b-x)**beta*log(x-a)*log(b-x)
      logical,optional                 :: verbose
      integer,optional                 :: Ninterp
      real(8)                          :: result
      !actual default variables
      real(8)                          :: epsabs_
      real(8)                          :: epsrel_
      logical                          :: verbose_
      real(8)                          :: abserr
      real(8),dimension(size(fsample)) :: xsample
      integer                          :: Neval
      integer                          :: Ier,i
      integer                          :: Nsingular_points
      integer                          :: Ninterp_,Lin
      character(len=20)                :: routine_name
      !
      type finter_type
         real(8),allocatable :: X(:)
         real(8),allocatable :: F(:)
         integer             :: Imin=0,Imax=0,N=0
         logical             :: status=.false.
      end type finter_type
      type(finter_type)                :: Finterp
      !
      !
      epsabs_ = 1d-12 ; if(present(epsabs))epsabs_=epsabs
      epsrel_ = 1d-6  ; if(present(epsrel))epsrel_=epsrel
      verbose_=.false.; if(present(verbose))verbose_=verbose
      Ninterp_=3      ; if(present(Ninterp))Ninterp_=Ninterp
      !
      !Build the interpolating function to be integrated:
      Lin = size(fsample)
      xsample = linspace(a,b,Lin)
      call init_finter(Finterp,xsample,fsample,Ninterp_)
      !
      if(present(key).AND.(present(alfa).OR.present(beta)))stop "ERROR in quad: key & alfa,beta"
      if(present(key).AND.present(cpole))stop "ERROR in quad: key & cpole"
      if(present(key).AND.present(singular_points))stop "ERROR in quad: key & singular_points"
      if(present(key).AND.present(singular_endpoint))stop "ERROR in quad: key & singular_endpoint"
      if(present(key).AND.present(omega))stop "ERROR in quad: key & omega"
      if(present(singular_endpoint).AND.(present(alfa).OR.present(beta)))stop "ERROR in quad: singular_endpoint & alfa,beta"
      if(present(singular_endpoint).AND.present(cpole))stop "ERROR in quad: singular_endpoint & cpole"
      if(present(singular_endpoint).AND.present(singular_points))stop "ERROR in quad: singular_endpoint & singular_points"
      if(present(singular_endpoint).AND.present(omega))stop "ERROR in quad: singular_endpoint & omega"
      if(present(cpole).AND.(present(alfa).OR.present(beta)))stop "ERROR in quad: cpole & alfa,beta"
      if(present(cpole).AND.present(singular_points))stop "ERROR in quad: cpole & singular_points"
      if(present(cpole).AND.present(omega))stop "ERROR in quad: cpole & omega"
      if(present(omega).AND.(present(alfa).OR.present(beta)))stop "ERROR in quad: omega & alfa,beta"
      if(present(omega).AND.present(singular_points))stop "ERROR in quad: omega & singular_points"
      if(present(singular_points).AND.(present(alfa).OR.present(beta)))stop "ERROR in quad: singular_points & alfa,beta"
      !
      !QNG:  no need
      !QAGS: need singular_endpoint=T
      !QAG : need key=1,6
      !QAGP: need singular_points
      !QAWC: need cpole
      !QAWO: need omega + weight_func=1,2
      !QAWS: need alfa + beta + weight_func=1,2,3,4
      routine_name='QNG'
      if(present(singular_endpoint))routine_name='QAGS'
      if(present(key))routine_name='QAG'
      if(present(singular_points))routine_name='QAGP'
      if(present(cpole))routine_name='QAWC'
      if(present(omega))routine_name='QAWO'
      if(present(alfa).AND.present(beta))routine_name='QAWS'
      !
      if(verbose_)write(*,"(A,A)")"QUAD: selected procedure =", routine_name
      !
      select case(routine_name)
       case ('QNG')
         call QNG(func,a,b,epsabs_,epsrel_,result,abserr,neval,ier)
         if(verbose_)then
            write(*,'(A,2F14.6)')'Endpoints (a,b)                             =', a,b
            write(*,'(A,F14.6)') 'Estimated integral                          =', result
            write(*,'(A,F14.6)') 'Estimated integral error                    =', abserr
            write(*,'(A,I8)')    'Error return code IER                       =', ier
         else
            if(ier/=0)then
               write(*,'(A,I8)') 'Error return code IER =', ier
               stop
            endif
         endif
         !
         !
       case ('QAGS')
         call QAGS(func,a,b,epsabs_,epsrel_,result,abserr,neval,ier)
         if(verbose_)then
            write(*,'(A,2F14.6)')'(Singular) endpoints (a*,b*)                =', a,b
            write(*,'(A,F14.6)') 'Estimated integral                          =', result
            write(*,'(A,F14.6)') 'Estimated integral error                    =', abserr
            write(*,'(A,I8)')    'Error return code IER                       =', ier
         else
            if(ier/=0)then
               write(*,'(A,I8)') 'Error return code IER =', ier
               stop
            endif
         endif
         !
         !
       case ('QAG')
         if(key<1.OR.key>6)stop "ERROR in quad: use QAG, key !=[1:6]"
         call QAG(func,a,b,epsabs_,epsrel_,key,result,abserr,neval,ier)
         if(verbose_)then
            write(*,'(A,2F14.6)')'Endpoints (a,b)                             =', a,b
            write(*,'(A,I8)')    'Approximation or                            =', key
            write(*,'(A,F14.6)') 'Estimated integral                          =', result
            write(*,'(A,F14.6)') 'Estimated integral error                    =', abserr
            write(*,'(A,I8)')    'Error return code IER                       =', ier
         else
            if(ier/=0)then
               write(*,'(A,I8)') 'Error return code IER =', ier
               stop
            endif
         endif
         !
         !
       case ('QAGP')
         Nsingular_points=size(singular_points)
         call QAGP(func,a,b,Nsingular_points,singular_points,epsabs_,epsrel_,result,abserr,neval,ier)
         if(verbose_)then
            write(*,'(A,2F14.6)')'Endpoints (a,b)                             =', a,b
            write(*,'(A,I8)')'N_singular_points                           =', Nsingular_points
            do i=1,Nsingular_points
               write(*,"(A,I6,F14.6)")'I_singular_point                            =',i,singular_points(i)
            enddo
            write(*,'(A,F14.6)') 'Estimated integral                          =', result
            write(*,'(A,F14.6)') 'Estimated integral error                    =', abserr
            write(*,'(A,I8)')    'Error return code IER                       =', ier
         else
            if(ier/=0)then
               write(*,'(A,I8)') 'Error return code IER =', ier
               stop
            endif
         endif
         !
         !
       case ('QAWC')
         call QAWC(func,a,b,cpole,epsabs_,epsrel_,result,abserr,neval,ier)
         if(verbose_)then
            write(*,'(A,2F14.6)')'Endpoints (a,b)                             =', a,b
            write(*,'(A,F14.6)')'C_pole                                      =', cpole
            write(*,'(A,F14.6)') 'Estimated integral                          =', result
            write(*,'(A,F14.6)') 'Estimated integral error                    =', abserr
            write(*,'(A,I8)')    'Error return code IER                       =', ier
         else
            if(ier/=0)then
               write(*,'(A,I8)') 'Error return code IER =', ier
               stop
            endif
         endif
         !
         !
       case ('QAWO')
         if(weight_func<1.OR.weight_func>2)stop "ERROR in quad: use QAWO, weight_func !=[1,2]"
         if(verbose_)then
            if(weight_func==1)then
               write(*,'(A)')'W(x) function is cos(omega*x)'
            elseif(weight_func==2)then
               write(*,'(A)')'W(x) function is sin(omega*x)'
            endif
         endif
         call QAWO(func,a,b,omega,weight_func,epsabs_,epsrel_,result,abserr,neval,ier)
         if(verbose_)then
            write(*,'(A,2F14.6)')'Endpoints (a,b)                             =', a,b
            write(*,'(A,F14.6)') 'Omega                                       =', omega
            write(*,'(A,F14.6)') 'Estimated integral                          =', result
            write(*,'(A,F14.6)') 'Estimated integral error                    =', abserr
            write(*,'(A,I8)')    'Error return code IER                       =', ier
         else
            if(ier/=0)then
               write(*,'(A,I8)') 'Error return code IER =', ier
               stop
            endif
         endif
         !
         !
       case ('QAWS')
         if(weight_func<1.OR.weight_func>4)stop "ERROR in quad: use QAWS, weight_func !=[1,..,4]"
         if(verbose_)then
            if(weight_func==1)then
               write(*,'(A)')'W(x) function is (x-a)**alfa*(b-x)**beta         '
            elseif(weight_func==2)then
               write(*,'(A)')'W(x) function is (x-a)**alfa*(b-x)**beta*log(x-a)'
            elseif(weight_func==3)then
               write(*,'(A)')'W(x) function is (x-a)**alfa*(b-x)**beta*log(b-x) '
            elseif(weight_func==4)then
               write(*,'(A)')'W(x) function is (x-a)**alfa*(b-x)**beta*log(x-a)*log(b-x)'
            endif
         endif
         call QAWS(func,a,b,alfa,beta,weight_func,epsabs_,epsrel_,result,abserr,neval,ier)
         if(verbose_)then
            write(*,'(A,2F14.6)')'Endpoints (a,b)                             =', a,b
            write(*,'(A,2F14.6)')'Parameters (alfa,beta)                      =', alfa,beta
            write(*,'(A,F14.6)') 'Estimated integral                          =', result
            write(*,'(A,F14.6)') 'Estimated integral error                    =', abserr
            write(*,'(A,I8)')    'Error return code IER                       =', ier
         else
            if(ier/=0)then
               write(*,'(A,I8)') 'Error return code IER =', ier
               stop
            endif
         endif
         !
         !
      end select
      !
      call delete_finter(finterp)
      !
      !
   contains
      !
      !
      subroutine init_finter(func,Xin,Fin,N)
         type(finter_type) :: func
         real(8)           :: xin(:)
         real(8)           :: fin(size(xin))
         integer           :: N,Lin
         if(func%status)deallocate(func%x,func%f)
         Lin=size(xin)
         allocate(func%x(Lin),func%f(Lin))
         func%X    = Xin
         func%F    = Fin
         func%Imax = Lin
         func%Imin = 1
         func%N    = N
         func%status=.true.
      end subroutine init_finter
      !
      subroutine delete_finter(func)
         type(finter_type) :: func
         if(allocated(func%x))deallocate(func%x)
         if(allocated(func%f))deallocate(func%f)
         func%imax=0
         func%imin=0
         func%N   =0
         func%status=.false.
      end subroutine delete_finter
      !
      function func(x)
         real(8)           :: x
         real(8)           :: func
         real(8)           :: dy
         integer           :: j,k,k0,k1
         integer           :: n
         N=finterp%N    !order of polynomial interpolation
         func=0d0
         j=locate(finterp%X(finterp%Imin:finterp%Imax),x)
         !
         k=max(j-(N-1)/2,1)
         k0=k
         if(k0 < finterp%Imin)k0=finterp%Imin
         k1=k+N+1
         if(k1 > finterp%Imax)then
            k1=finterp%Imax
            k0=k1-N-1
         endif
         call polint(finterp%X(k0:k1),finterp%F(k0:k1),x,func,dy)
      end function func
   end subroutine quad_sample




   !+-----------------------------------------------------------------+
   !PURPOSE  : Perform a fast Kramers-K\"onig integration:
   !+-----------------------------------------------------------------+
   function kronig(fi,wr,M) result(fr)
      integer :: i,j,M
      real(8),dimension(M) :: fi,wr,fr
      real(8),dimension(M) :: logo,deriv
      real(8) :: dh,sum
      dh=wr(2)-wr(1)
      logo=0.d0
      do i=2,M-1
         logo(i) = log( (wr(M)-wr(i))/(wr(i)-wr(1)) )
      enddo
      deriv(1)= (fi(2)-fi(1))/dh
      deriv(M)= (fi(M)-fi(M-1))/dh
      do i=2,M-1
         deriv(i) = (fi(i+1)-fi(i-1))/(2*dh)
      enddo
      fr=0.d0
      do i=1,M
         sum=0.d0
         do j=1,M
            if(i/=j)then
               sum=sum+(fi(j)-fi(i))*dh/(wr(j)-wr(i))
            else
               sum=sum+deriv(i)*dh
            endif
         enddo
         fr(i) = (sum + fi(i)*logo(i))/pi
      enddo
      return
   end function kronig















   !*******************************************************************
   !*******************************************************************
   !*******************************************************************
   !*******************************************************************
   !*******************************************************************
   !*******************************************************************



   !+-----------------------------------------------------------------+
   !PURPOSE: obtain quadrature weights for higher order integration (2,4)
   !+-----------------------------------------------------------------+
   subroutine get_quadrature_weights(wt,nrk)
      real(8),dimension(:) :: wt
      integer,optional     :: nrk
      integer              :: nrk_
      integer              :: N
      nrk_=4;if(present(nrk))nrk_=nrk
      N=size(wt)
      if(nrk_==4)then
         select case(n)           !n>=3
          case (1)
            wt = 1.d0
          case (2)
            wt = 0.5d0
          case (3)                 !simpson's rule
            wt(1)=1.d0/3.d0
            wt(2)=4.d0/3.d0
            wt(3)=1.d0/3.d0
          case(4)                  !simpson's 3/8 rule
            wt(1)=3.d0/8.d0
            wt(2)=9.d0/8.d0
            wt(3)=9.d0/8.d0
            wt(4)=3.d0/8.d0
          case(5)                  !Simpson's rule (E,O n)
            wt(1)=1.d0/3.d0
            wt(2)=4.d0/3.d0
            wt(3)=2.d0/3.d0
            wt(4)=4.d0/3.d0
            wt(5)=1.d0/3.d0
          case default            !Simpson's rule n>=6
            if(mod(n-1,2)==0)then
               wt(1)=1.d0/3.d0
               wt(n)=1.d0/3.d0
               wt(2:n-1:2)=4.d0/3.d0
               wt(3:n-2:2)=2.d0/3.d0
            else
               wt(1)=1.d0/3.d0
               wt(2:n-4:2)=4.d0/3.d0
               wt(3:n-5:2)=2.d0/3.d0
               wt(n-3)=17.d0/24.d0
               wt(n-2)=9.d0/8.d0
               wt(n-1)=9.d0/8.d0
               wt(n)=3.d0/8.d0
            endif
            ! case default             !Simpson's rule n>=6
            !    wt(1)=3.d0/8.d0
            !    wt(2)=7.d0/6.d0
            !    wt(3)=23.d0/24.d0
            !    wt(4:n-3)=1.d0
            !    wt(n-2)=23.d0/24.d0
            !    wt(n-1)=7.d0/6.d0
            !    wt(n)=3.d0/8.d0
         end select
      elseif(nrk_==2)then
         wt(1) = 0.5d0
         wt(2:n-1)=1.d0
         wt(n) = 0.5d0
      else
         stop "error in +get_quadrature_weights: nrk != 2,4"
      end if
   end subroutine get_quadrature_weights





   function linspace(start,stop,num,istart,iend,mesh) result(array)
      integer          :: num,i
      real(8)          :: start,stop,step,array(num)
      logical,optional :: istart,iend
      logical          :: startpoint_,endpoint_
      real(8),optional :: mesh
      if(num<0)stop "linspace: N<0, abort."
      startpoint_=.true.;if(present(istart))startpoint_=istart
      endpoint_=.true.;if(present(iend))endpoint_=iend
      if(startpoint_.AND.endpoint_)then
         if(num<2)stop "linspace: N<2 with both start and end points"
         step = (stop-start)/real(num-1,8)
         forall(i=1:num)array(i)=start + real(i-1,8)*step
      elseif(startpoint_.AND.(.not.endpoint_))then
         step = (stop-start)/real(num,8)
         forall(i=1:num)array(i)=start + real(i-1,8)*step
      elseif(.not.startpoint_.AND.endpoint_)then
         step = (stop-start)/real(num,8)
         forall(i=1:num)array(i)=start + real(i,8)*step
      else
         step = (stop-start)/real(num+1,8)
         forall(i=1:num)array(i)=start + real(i,8)*step
      endif
      if(present(mesh))mesh=step
   end function linspace
















   subroutine polint(xa,ya,x,y,dy)
      real(8), dimension(:), intent(in) :: xa,ya
      real(8), intent(in)          :: x
      real(8), intent(out)         :: y,dy
      integer                      :: m,n,ns
      real(8), dimension(size(xa)) :: c,d,den,ho
      n=assert_eq2(size(xa),size(ya),'polint')
      c=ya
      d=ya
      ho=xa-x
      ns=iminloc(abs(x-xa))
      y=ya(ns)
      ns=ns-1
      do m=1,n-1
         den(1:n-m)=ho(1:n-m)-ho(1+m:n)
         if (any(den(1:n-m) == 0.0))then
            print*,'polint: calculation failure'
            stop
         endif
         den(1:n-m)=(c(2:n-m+1)-d(1:n-m))/den(1:n-m)
         d(1:n-m)=ho(1+m:n)*den(1:n-m)
         c(1:n-m)=ho(1:n-m)*den(1:n-m)
         if (2*ns < n-m) then
            dy=c(ns+1)
         else
            dy=d(ns)
            ns=ns-1
         end if
         y=y+dy
      end do
   end subroutine polint
   function locate(xx,x)
      real(8), dimension(:), intent(in) :: xx
      real(8), intent(in) :: x
      integer :: locate
      integer :: n,jl,jm,ju
      logical :: ascnd
      n=size(xx)
      ascnd = (xx(n) >= xx(1))
      jl=0
      ju=n+1
      do
         if (ju-jl <= 1) exit
         jm=(ju+jl)/2
         if (ascnd .eqv. (x >= xx(jm))) then
            jl=jm
         else
            ju=jm
         end if
      end do
      if (x == xx(1)) then
         locate=1
      else if (x == xx(n)) then
         locate=n-1
      else
         locate=jl
      end if
   end function locate
   function iminloc(arr)
      real(8), dimension(:), intent(in) :: arr
      integer, dimension(1) :: imin
      integer :: iminloc
      imin=minloc(arr(:))
      iminloc=imin(1)
   end function iminloc
   function assert_eq2(n1,n2,string)
      character(len=*), intent(in) :: string
      integer, intent(in) :: n1,n2
      integer :: assert_eq2
      if (n1 == n2) then
         assert_eq2=n1
      else
         write (*,*) 'nrerror: an assert_eq failed with this tag:', &
            string
         stop 'program terminated by assert_eq2'
      end if
   end function assert_eq2

end module SF_INTEGRATE






















! !+-----------------------------------------------------------------+
! !PURPOSE  : A slower (but possibly more accurate) Kramers-Kr\"onig
! ! integral, using local interpolation w/ polynomial of order-5
! !+-----------------------------------------------------------------+
! function kramers_kronig(fi,x,L) result(fr)
!   integer                 :: i,L
!   real(8),dimension(L)    :: fi,fr,x
!   real(8)                 :: dx
!   integer, parameter :: LIMIT = 500
!   integer            :: IERSUM,NEVAL,IER,LAST
!   real(8)            :: EPSABS,EPSREL,A,B,C,ABSERR
!   real(8)            :: ALIST(LIMIT),BLIST(LIMIT),ELIST(LIMIT),RLIST(LIMIT)
!   integer            :: IORD(limit)
!   if(allocated(finterX))deallocate(finterX)
!   if(allocated(finterF))deallocate(finterF)
!   allocate(finterX(1:L),finterF(1:L))
!   finterX    = x
!   finterF    = fi/pi
!   finterImin = 1
!   finterImax = L
!   finterN    = 5
!   EPSABS = 0.0d0
!   EPSREL = 1.d-12
!   IERSUM = 0
!   A      = x(1)             !minval(x)
!   B      = x(L)             !maxval(x)
!   do i=1,L
!      C = x(i)
!      CALL QAWCE(kkfinter,A,B,C,EPSABS,EPSREL,LIMIT,fr(i),ABSERR,NEVAL,&
!           IER,alist,blist,rlist,elist,iord,last)
!      IERSUM=IERSUM+IER
!   enddo
!   !Some regularization at the borders: (thanks Jan Tomczak)
!   dx=x(2)-x(1)
!   fr(L) =(fr(L-2) - fr(L-1))*dx  + fr(L-1)
!   fr(1)=(fr(1+1)- fr(1+2))*dx + fr(1+1)
! end function kramers_kronig
! !
! function kkfinter(x) result(finter)
!   real(8) :: finter
!   real(8) :: x,y,dy
!   integer :: itmp,k
!   integer :: n
!   n=finterN    !order of polynomial interpolation
!   finter=0.d0
!   itmp=locate(finterX(FinterImin:finterImax),x)
!   k=max(itmp-(N-1)/2,1)
!   if (k < finterImin)k=finterImin
!   if(k+n+1 <= finterImax)then
!      call polint(finterX(k:k+n+1),finterF(k:k+n+1),x,y,dy)
!   else
!      call polint(finterX(k:finterImax),finterF(k:finterImax),x,y,dy)
!   endif
!   finter=y
!   return
! end function kkfinter
