  include "d_quadpack.f90"
  include "MOD_FINTER.f90"

  !###############################################################
  ! PROGRAM  : INTEGRATE
  ! TYPE     : Module
  ! PURPOSE  : A set of routines to perform specific integrals
  !###############################################################
  module INTEGRATE
    USE MOD_FINTER
    implicit none
    private

    complex(8),parameter :: one= (1.d0,0.d0)
    complex(8),parameter :: xi = (0.d0,1.d0)
    real(8),parameter    :: pi=3.14159265358979323846264338327950288419716939937510D0

    public :: kramers_kronig
    public :: kronig

  contains

    !+-----------------------------------------------------------------+
    !PROGRAM  : 
    !TYPE     : Subroutine
    !PURPOSE  :   
    !+-----------------------------------------------------------------+
    function kronig(fi,wr,M) result(fr)
      real(8),dimension(M) :: fi,wr,fr
      real(8),dimension(M) :: logo,deriv
      real(8) :: dh,sum
      integer :: i,j,M

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


    !+-----------------------------------------------------------------+
    !PROGRAM  : 
    !TYPE     : Subroutine
    !PURPOSE  :   
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
      A      = x(1)            !minval(x)      
      B      = x(L)             !maxval(x)
      do i=1,L
         C = x(i)
         CALL QAWCE(finter,A,B,C,EPSABS,EPSREL,LIMIT,fr(i),ABSERR,NEVAL,&
              IER,alist,blist,rlist,elist,iord,last)
         IERSUM=IERSUM+IER
      enddo
      !Some regularization at the borders: (thank Jan Tomczak)
      dx=x(2)-x(1)
      fr(L) =(fr(L-2) - fr(L-1))*dx  + fr(L-1)
      fr(1)=(fr(1+1)- fr(1+2))*dx + fr(1+1)
    end function kramers_kronig

    ! function kramers_kronig(fi,x,L) result(fr)
    !   integer                 :: i,L
    !   real(8),dimension(-L:L) :: fi,fr,x
    !   real(8)                 :: dx
    !   !QAWC/E parameters:
    !   integer, parameter :: LIMIT = 500
    !   integer            :: IERSUM,NEVAL,IER,LAST
    !   real(8)            :: EPSABS,EPSREL,A,B,C,ABSERR
    !   real(8)            :: ALIST(LIMIT),BLIST(LIMIT),ELIST(LIMIT),RLIST(LIMIT)
    !   integer            :: IORD(limit)
    !   if(allocated(finterX))deallocate(finterX)
    !   if(allocated(finterF))deallocate(finterF)
    !   allocate(finterX(-L:L),finterF(-L:L))
    !   finterX    = x
    !   finterF    = fi/pi
    !   finterImax = L
    !   finterImin =-L
    !   finterN    = 5
    !   EPSABS = 0.0d0
    !   EPSREL = 5.0d-9
    !   IERSUM = 0
    !   A      = x(-L)-1.d-3
    !   B      = x(L)+1.d-3
    !   do i=-L,L
    !      C = x(i)
    !      CALL QAWCE(finter,A,B,C,EPSABS,EPSREL,LIMIT,fr(i),ABSERR,NEVAL,&
    !           IER,alist,blist,rlist,elist,iord,last)
    !      IERSUM=IERSUM+IER
    !   enddo
    !   !Some regularization at the borders: (thanks Jan Tomczak)
    !   dx=x(1)-x(0)
    !   fr(L) =(fr(L-2) - fr(L-1))*dx  + fr(L-1)
    !   fr(-L)=(fr(-L+1)- fr(-L+2))*dx + fr(-L+1)
    ! end function kramers_kronig

    !*******************************************************************
    !*******************************************************************
    !*******************************************************************





  end module INTEGRATE
