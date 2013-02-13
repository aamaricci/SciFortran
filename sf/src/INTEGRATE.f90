  include "integrate_d_quadpack.f90"
  include "integrate_finter_mod.f90"
  !###############################################################
  ! PROGRAM  : INTEGRATE
  ! TYPE     : Module
  ! PURPOSE  : A set of routines to perform specific integrals
  !###############################################################
  module INTEGRATE
    USE INTEGRATE_FINTER_MOD
    implicit none
    private

    complex(8),parameter :: one= (1.d0,0.d0)
    complex(8),parameter :: zero=(0.d0,0.d0)
    complex(8),parameter :: xi = (0.d0,1.d0)
    real(8),parameter    :: pi=3.14159265358979323846264338327950288419716939937510D0

    interface trapz
       module procedure &
            d_trapz_ab,c_trapz_ab,&
            d_trapz_dh,c_trapz_dh,&
            d_trapz_nonlin,c_trapz_nonlin
    end interface trapz

    interface simps
       module procedure &
            d_simpson_ab,d_simpson_dh,&
            c_simpson_ab,c_simpson_dh,&
            d_simpson_nonlin,c_simpson_nonlin
    end interface simps

    public :: kramers_kronig
    public :: kronig
    public :: trapz
    public :: simps


  contains

    !+-----------------------------------------------------------------+
    !PURPOSE: Trapezoidal rule for data function integration
    !+-----------------------------------------------------------------+
    function d_trapz_ab(a,b,f) result(sum)
      real(8) :: a,b,dh
      real(8) :: f(:)
      real(8) :: sum
      integer :: i,L
      L=size(f)
      dh=(b-a)/real(L-1,8)/2.d0
      sum=0.d0
      do i=1,L-1
         sum = sum+(f(i+1)+f(i))*dh
      enddo
    end function d_trapz_ab

    function d_trapz_dh(dh,f) result(sum)
      real(8) :: dh
      real(8) :: f(:)
      real(8) :: sum
      integer :: i,L
      L=size(f)
      sum=0.d0
      do i=1,L-1
         sum = sum+(f(i+1)+f(i))*dh/2.d0
      enddo
    end function d_trapz_dh

    function c_trapz_ab(a,b,f) result(sum)
      real(8)    :: a,b,dh
      complex(8) :: f(:)
      complex(8) :: sum
      integer    :: i,L
      L=size(f)
      dh=(b-a)/real(L-1,8)/2.d0
      sum=0.d0
      do i=1,L-1
         sum = sum+(f(i+1)+f(i))*dh
      enddo
    end function c_trapz_ab

    function c_trapz_dh(dh,f) result(sum)
      real(8)    :: dh
      complex(8) :: f(:)
      complex(8) :: sum
      integer    :: i,L
      L=size(f)
      sum=0.d0
      do i=1,L-1
         sum = sum+(f(i+1)+f(i))*dh/2.d0
      enddo
    end function c_trapz_dh

    function d_trapz_nonlin(x,f) result(sum)
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
    end function d_trapz_nonlin

    function c_trapz_nonlin(x,f) result(sum)
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
    end function c_trapz_nonlin


    !+-----------------------------------------------------------------+
    !PURPOSE: Simpson rule for data function integration
    !+-----------------------------------------------------------------+
    function d_simpson_ab(a,b,f) result(sum)
      integer :: n
      real(8) :: a,b
      real(8) :: f(:)
      real(8) :: dh,sum,sum1,sum2,int1,int2
      integer :: i,p,m,mm,mmm
      if(b==a)then
         sum=0.d0
         return
      endif
      n=size(f)
      m=n-1
      dh=(b-a)/real(m,8)
      p=(m/2)*2
      sum1=0.d0
      sum2=0.d0
      sum=0.d0
      if(p==m)then
         do i=2,m,2
            sum1 = sum1 + f(i)
         enddo
         do i=3,m-1,2
            sum2 = sum2 + f(i)
         enddo
         sum = (f(1) + 4.d0*sum1 + 2.d0*sum2 + f(n))/3.d0
         sum = sum*dh
      else
         if (m>=5) then
            mm=m-3
            do i=2,m-3,2
               sum1 = sum1 + f(i)
            enddo
            mmm=mm-1
            do i=3,m-4,2
               sum2 = sum2 + f(i)
            enddo
            int1 = (f(1) + 4.d0*sum1 + 2.d0*sum2 + f(n-3))/(3.d0*(m-3))
            int1 = ((b-3.d0*dh)-a)*int1
         endif
         int2 = 3.d0*dh*(f(n-3)+3.d0*(f(n-2)+f(n-1))+f(n))/8.d0
         sum  = int1 + int2
      end if
    end function d_simpson_ab

    function c_simpson_ab(a,b,f) result(sum)
      real(8)    :: a,b
      complex(8) :: f(:)
      complex(8) :: sum
      real(8)    :: rsum,isum
      rsum = d_simpson_ab(a,b,dreal(f))
      isum = d_simpson_ab(a,b,dimag(f))
      sum  = cmplx(rsum,isum,8)
    end function c_simpson_ab

    function d_simpson_dh(dh,f) result(sum)
      real(8) :: dh,a,b
      real(8) :: f(:)
      real(8) :: sum
      a=0.d0 ; b=size(f)*dh-dh
      sum = d_simpson_ab(a,b,f)
    end function d_simpson_dh

    function c_simpson_dh(dh,f) result(sum)
      real(8)    :: dh,a,b
      complex(8) :: f(:)
      complex(8) :: sum
      a=0.d0 ; b=size(f)*dh-dh
      sum = c_simpson_ab(a,b,f)
    end function c_simpson_dh

    function d_simpson_nonlin(x,f) result(sum)
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
    end function d_simpson_nonlin

    function c_simpson_nonlin(x,f) result(sum)
      complex(8) :: f(:),sum
      real(8)    :: x(size(f)),rsum,isum
      rsum=d_simpson_nonlin(x,dreal(f))
      isum=d_simpson_nonlin(x,dimag(f))
      sum  = cmplx(rsum,isum,8)
    end function c_simpson_nonlin






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

    !+-----------------------------------------------------------------+
    !PURPOSE  : A slower (but possibly more accurate) Kramers-Kr\"onig 
    ! integral, using local interpolation w/ polynomial of order-5
    !+-----------------------------------------------------------------+
    function kramers_kronig(fi,x,L) result(fr)
      integer                 :: i,L
      real(8),dimension(L)    :: fi,fr,x
      real(8)                 :: dx
      integer, parameter :: LIMIT = 500
      integer            :: IERSUM,NEVAL,IER,LAST
      real(8)            :: EPSABS,EPSREL,A,B,C,ABSERR
      real(8)            :: ALIST(LIMIT),BLIST(LIMIT),ELIST(LIMIT),RLIST(LIMIT)
      integer            :: IORD(limit)
      if(allocated(finterX))deallocate(finterX)
      if(allocated(finterF))deallocate(finterF)
      allocate(finterX(1:L),finterF(1:L))
      finterX    = x
      finterF    = fi/pi
      finterImin = 1
      finterImax = L
      finterN    = 5
      EPSABS = 0.0d0
      EPSREL = 1.d-12
      IERSUM = 0
      A      = x(1)             !minval(x)      
      B      = x(L)             !maxval(x)
      do i=1,L
         C = x(i)
         CALL QAWCE(finter,A,B,C,EPSABS,EPSREL,LIMIT,fr(i),ABSERR,NEVAL,&
              IER,alist,blist,rlist,elist,iord,last)
         IERSUM=IERSUM+IER
      enddo
      !Some regularization at the borders: (thanks Jan Tomczak)
      dx=x(2)-x(1)
      fr(L) =(fr(L-2) - fr(L-1))*dx  + fr(L-1)
      fr(1)=(fr(1+1)- fr(1+2))*dx + fr(1+1)
    end function kramers_kronig




    !*******************************************************************
    !*******************************************************************
    !*******************************************************************





  end module INTEGRATE
