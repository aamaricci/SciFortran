! now in a separate library with separate MODULE interface.
! Here we provide only standard functions and 
! will add more object oriented interface to QUADPACK
! which will require to have this lib installed at compilation time.
! public :: qc25c
! public :: qc25o
! public :: qc25s
! public :: qcheb
! public :: qextr
! public :: qfour
! public :: qk15
! public :: qk15i
! public :: qk15w
! public :: qk21
! public :: qk31
! public :: qk41
! public :: qk51
! public :: qk61
! public :: qmomo
! public :: qsort
! public :: qag
! public :: qage
! public :: qagi
! public :: qagp
! public :: qags
! public :: qawc
! public :: qawce
! public :: qawf
! public :: qawfe
! public :: qawo
! public :: qaws
! public :: qawse
! public :: qng

module SF_INTEGRATE
  USE SF_CONSTANTS, only: one,zero,xi,pi
  USE SF_ARRAYS, only: linspace
  implicit none
  private

  !TRAPEZIOSAL RULE:
  interface trapz
     module procedure d_trapz_func
     module procedure c_trapz_func
     ! module procedure d_trapz_sample
     ! module procedure c_trapz_sample
     module procedure d_trapz_dh_sample
     module procedure c_trapz_dh_sample
     module procedure d_trapz_ab_sample
     module procedure c_trapz_ab_sample
     module procedure d_trapz_nonlin_sample
     module procedure c_trapz_nonlin_sample
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
  !
  public :: trapz
  public :: trapz2d



  !SIMPSON'S RULE
  interface simps
     module procedure d_simps_func
     module procedure c_simps_func
     ! module procedure d_simps_sample
     ! module procedure c_simps_sample
     module procedure d_simpson_dh_sample
     module procedure c_simpson_dh_sample
     module procedure d_simpson_ab_sample
     module procedure c_simpson_ab_sample
     module procedure d_simpson_nonlin_sample
     module procedure c_simpson_nonlin_sample
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
  public :: simps
  public :: simps2d



  !KRAMERS-KRONIG:
  public :: kronig


  !<TODO
  ! add routines for the 3d case
  ! add interaface to QUADPACK
  ! add Montecarlo base 1d/2d/3d integrals
  !>TODO



contains


  !+-----------------------------------------------------------------------------+!
  !PURPOSE:
  ! evaluate 1D and 2D integrals using trapz (2nd order) and simps (4th order) rule.
  !+-----------------------------------------------------------------------------+!
  include "integrate_func_1d.f90"
  include "integrate_sample_1d.f90"
  include "integrate_func_2d.f90"
  include "integrate_sample_2d.f90"





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
  ! subroutine polint(xa,ya,x,y,dy)
  !   real(8), dimension(:), intent(in) :: xa,ya
  !   real(8), intent(in)          :: x
  !   real(8), intent(out)         :: y,dy
  !   integer                      :: m,n,ns
  !   real(8), dimension(size(xa)) :: c,d,den,ho
  !   n=assert_eq2(size(xa),size(ya),'polint')
  !   c=ya
  !   d=ya
  !   ho=xa-x
  !   ns=iminloc(abs(x-xa))
  !   y=ya(ns)
  !   ns=ns-1
  !   do m=1,n-1
  !      den(1:n-m)=ho(1:n-m)-ho(1+m:n)
  !      if (any(den(1:n-m) == 0.0))then
  !         print*,'polint: calculation failure'
  !         stop
  !      endif
  !      den(1:n-m)=(c(2:n-m+1)-d(1:n-m))/den(1:n-m)
  !      d(1:n-m)=ho(1+m:n)*den(1:n-m)
  !      c(1:n-m)=ho(1:n-m)*den(1:n-m)
  !      if (2*ns < n-m) then
  !         dy=c(ns+1)
  !      else
  !         dy=d(ns)
  !         ns=ns-1
  !      end if
  !      y=y+dy
  !   end do
  ! end subroutine polint
  ! function locate(xx,x)
  !   real(8), dimension(:), intent(in) :: xx
  !   real(8), intent(in) :: x
  !   integer :: locate
  !   integer :: n,jl,jm,ju
  !   logical :: ascnd
  !   n=size(xx)
  !   ascnd = (xx(n) >= xx(1))
  !   jl=0
  !   ju=n+1
  !   do
  !      if (ju-jl <= 1) exit
  !      jm=(ju+jl)/2
  !      if (ascnd .eqv. (x >= xx(jm))) then
  !         jl=jm
  !      else
  !         ju=jm
  !      end if
  !   end do
  !   if (x == xx(1)) then
  !      locate=1
  !   else if (x == xx(n)) then
  !      locate=n-1
  !   else
  !      locate=jl
  !   end if
  ! end function locate
  ! function iminloc(arr)
  !   real(8), dimension(:), intent(in) :: arr
  !   integer, dimension(1) :: imin
  !   integer :: iminloc
  !   imin=minloc(arr(:))
  !   iminloc=imin(1)
  ! end function iminloc
  ! function assert_eq2(n1,n2,string)
  !   character(len=*), intent(in) :: string
  !   integer, intent(in) :: n1,n2
  !   integer :: assert_eq2
  !   if (n1 == n2) then
  !      assert_eq2=n1
  !   else
  !      write (*,*) 'nrerror: an assert_eq failed with this tag:', &
  !           string
  !      stop 'program terminated by assert_eq2'
  !   end if
  ! end function assert_eq2



end module SF_INTEGRATE
