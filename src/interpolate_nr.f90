module INTERPOLATE_NR
  implicit none
  private

  public :: locate
  public :: polint
  public :: polin2
contains



  function locate(xx,x)
    REAL(8), DIMENSION(:), INTENT(IN) :: xx
    REAL(8), INTENT(IN) :: x
    INTEGER :: locate
    INTEGER :: n,jl,jm,ju
    LOGICAL :: ascnd
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
  END FUNCTION locate




  SUBROUTINE polint(xa,ya,x,y,dy)
    REAL(8), DIMENSION(:), INTENT(IN) :: xa,ya
    REAL(8), INTENT(IN)          :: x
    REAL(8), INTENT(OUT)         :: y,dy
    INTEGER                      :: m,n,ns
    REAL(8), DIMENSION(size(xa)) :: c,d,den,ho
    n=assert_eq(size(xa),size(ya),'polint')
    c=ya
    d=ya
    ho=xa-x
    ns=iminloc(abs(x-xa))
    y=ya(ns)
    ns=ns-1
    do m=1,n-1
       den(1:n-m)=ho(1:n-m)-ho(1+m:n)
       if (any(den(1:n-m) == 0.d0))then
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
  END SUBROUTINE polint


  subroutine polin2(x1a,x2a,ya,x1,x2,y,dy)
    real(8), dimension(:), intent(in)   :: x1a,x2a
    real(8), dimension(:,:), intent(in) :: ya
    real(8), intent(in)                 :: x1,x2
    real(8), intent(out)                :: y,dy
    integer                             :: j,m,ndum
    real(8), dimension(size(x1a))       :: ymtmp
    real(8), dimension(size(x2a))       :: yntmp
    m = size(x1a);if(m/=size(ya,1))stop "POLINT: wrong dimensions m"
    ndum=size(x2a);if(ndum/=size(ya,2))stop "POLINT: wrong dimensions ndum"
    do j=1,m
       yntmp=ya(j,:)
       call polint(x2a,yntmp,x2,ymtmp(j),dy)
    end do
    call polint(x1a,ymtmp,x1,y,dy)
  end subroutine polin2


  function iminloc(arr)
    real(8), dimension(:), intent(in) :: arr
    integer, dimension(1) :: imin
    integer :: iminloc
    imin=minloc(arr(:))
    iminloc=imin(1)
  end function iminloc

  function assert_eq(n1,n2,string)
    character(len=*), intent(in) :: string
    integer, intent(in) :: n1,n2
    integer :: assert_eq
    if (n1 == n2) then
       assert_eq=n1
    else
       write (*,*) 'nrerror: an assert_eq failed with this tag:', &
            string
       stop 'program terminated by assert_eq'
    end if
  end function assert_eq

END MODULE INTERPOLATE_NR
