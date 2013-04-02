!###############################################################
! PROGRAM  : INTEGRATE
! PURPOSE  : A set of routines to perform specific integrals
!###############################################################
module INTEGRATE
  implicit none
  private

  complex(8),parameter :: one= (1.d0,0.d0)
  complex(8),parameter :: zero=(0.d0,0.d0)
  complex(8),parameter :: xi = (0.d0,1.d0)
  real(8),parameter    :: pi=3.14159265358979323846264338327950288419716939937510D0

  !used in the old version of kramers-kronig
  real(8), allocatable :: finterX(:)  !vector with frequencies
  real(8), allocatable :: finterF(:) !corresponding vector of functional values
  integer              :: finterImin,finterImax,finterN

  type,public :: d_finter
     integer                      :: Imin,Imax,Iorder
     real(8),dimension(:),pointer :: X,F
     logical                      :: status=.false.
  end type d_finter
  type,public :: c_finter
     integer                         :: Imin,Imax,Iorder
     real(8),dimension(:),pointer    :: X
     complex(8),dimension(:),pointer :: F
     logical                         :: status=.false.
  end type c_finter

  interface init_finter
     module procedure d_init_finter,c_init_finter
  end interface init_finter

  interface kill_finter
     module procedure d_kill_finter,c_kill_finter
  end interface kill_finter

  interface set_finter
     module procedure d_set_finter,c_set_finter
  end interface set_finter

  interface finter_func
     module procedure d_finter_func,c_finter_func
  end interface finter_func

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

  public :: init_finter
  public :: kill_finter
  public :: set_finter
  public :: finter_func


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
       CALL QAWCE(kkfinter,A,B,C,EPSABS,EPSREL,LIMIT,fr(i),ABSERR,NEVAL,&
            IER,alist,blist,rlist,elist,iord,last)
       IERSUM=IERSUM+IER
    enddo
    !Some regularization at the borders: (thanks Jan Tomczak)
    dx=x(2)-x(1)
    fr(L) =(fr(L-2) - fr(L-1))*dx  + fr(L-1)
    fr(1)=(fr(1+1)- fr(1+2))*dx + fr(1+1)
  end function kramers_kronig
  !
  function kkfinter(x) result(finter)
    real(8) :: finter
    real(8) :: x,y,dy
    integer :: itmp,k
    integer :: n
    n=finterN    !order of polynomial interpolation
    finter=0.d0
    itmp=locate(finterX(FinterImin:finterImax),x)
    k=max(itmp-(N-1)/2,1)
    if (k < finterImin)k=finterImin
    if(k+n+1 <= finterImax)then
       call polint(finterX(k:k+n+1),finterF(k:k+n+1),x,y,dy)
    else
       call polint(finterX(k:finterImax),finterF(k:finterImax),x,y,dy)
    endif
    finter=y
    return
  end function kkfinter




  !*******************************************************************
  !*******************************************************************
  !*******************************************************************

  subroutine d_init_finter(f,size)
    type(d_finter)  :: f
    integer       :: size
    if(f%status)call kill_finter(f)
    allocate(f%x(size),f%f(size))
    f%status=.true.
  end subroutine d_init_finter
  subroutine c_init_finter(f,size)
    type(c_finter)  :: f
    integer         :: size
    if(f%status)call kill_finter(f)
    allocate(f%x(size),f%f(size))
    f%status=.true.
  end subroutine c_init_finter


  subroutine d_kill_finter(f)
    type(d_finter) :: f
    if(f%status)then
       deallocate(f%X,f%F)
    else
       print*,"kill_finter: nothing to be done."
    endif
  end subroutine d_kill_finter
  subroutine c_kill_finter(f)
    type(c_finter) :: f
    if(f%status)then
       deallocate(f%X,f%F)
    else
       print*,"kill_finter: nothing to be done."
    endif
  end subroutine c_kill_finter


  subroutine d_set_finter(f,fsize,fmin,fmax,forder,fx,fy)
    type(d_finter)           :: f
    integer                  :: fsize,fmin,fmax,forder
    real(8),dimension(fsize) :: fx,fy
    f%imin=fmin
    f%imax=fmax
    f%iorder=forder
    f%x(1:fsize)=fx(1:fsize)
    f%f(1:fsize)=fy(1:fsize)
  end subroutine d_set_finter
  subroutine c_set_finter(f,fsize,fmin,fmax,forder,fx,fy)
    type(c_finter)              :: f
    integer                     :: fsize,fmin,fmax,forder
    real(8),dimension(fsize)    :: fx
    complex(8),dimension(fsize) :: fy
    f%imin=fmin
    f%imax=fmax
    f%iorder=forder
    f%x(1:fsize)=fx(1:fsize)
    f%f(1:fsize)=fy(1:fsize)
  end subroutine c_set_finter


  function d_finter_func(x,f) result(finter)
    real(8)       :: x
    type(d_finter):: f
    real(8)       :: finter
    real(8)       :: y,dy
    integer       :: itmp,k
    integer       :: n
    n=f%iorder    !order of polynomial interpolation
    finter=0.d0
    itmp=locate(f%x(f%imin:f%imax),x)
    k=max(itmp-(N-1)/2,1)
    if (k < f%imin)k=f%imin
    if(k+n+1 <= f%imax)then
       call polint(f%x(k:k+n+1),f%f(k:k+n+1),x,y,dy)
    else
       call polint(f%x(k:f%imax),f%f(k:f%imax),x,y,dy)
    endif
    finter=y
  end function d_finter_func
  function c_finter_func(x,f) result(finter)
    real(8)        :: x
    type(c_finter) :: f
    complex(8)     :: finter
    real(8)        :: rey,imy,dy
    integer        :: itmp,k
    integer        :: n
    n=f%iorder    !order of polynomial interpolation    
    itmp=locate(f%x(f%imin:f%imax),x)
    k=max(itmp-(N-1)/2,1)
    if (k < f%imin)k=f%imin
    if(k+n+1 <= f%imax)then
       call polint(f%x(k:k+n+1),dreal(f%f(k:k+n+1)),x,rey,dy)
       call polint(f%x(k:k+n+1),dimag(f%f(k:k+n+1)),x,imy,dy)
    else
       call polint(f%x(k:f%imax),dreal(f%f(k:f%imax)),x,rey,dy)
       call polint(f%x(k:f%imax),dimag(f%f(k:f%imax)),x,imy,dy)
    endif
    finter=cmplx(rey,imy,8)
  end function c_finter_func


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



end module INTEGRATE
