  include "spline_nr_mod.f90"
  include "spline_finter_mod.f90"
  include "spline_interp.f90"
  include "spline_cubspl_routines.f90"
  !###############################################################
  !     PROGRAM  : SPLINE
  !     TYPE     : Module
  !     PURPOSE  : provide an interface to interpolation routines
  ! Now I don't have time to add more stuff, that can be  found 
  ! in splinF directory as copied from internet.
  !###############################################################
  module SPLINE
    USE SPLINE_NR_MOD
    USE SPLINE_FINTER_MOD
    implicit none
    private

    interface poly_spline
       module procedure d_fspline,c_fspline
    end interface poly_spline

    interface cubic_spline
       module procedure d_cub_interp,c_cub_interp
    end interface cubic_spline

    interface linear_spline
       module procedure d_linear_spline,c_linear_spline
    end interface linear_spline

    interface extract
       module procedure extract_gtau
    end interface extract

    public :: poly_spline
    public :: cubic_spline
    public :: linear_spline
    public :: extract,extract_gtau 
    public :: interp_gtau


  contains
    !+-------------------------------------------------------------------+
    !PROGRAM  : LINEAR_SPLINE
    !TYPE     : Subroutine
    !PURPOSE  : Linear interpolation of given data
    !+-------------------------------------------------------------------+
    subroutine d_linear_spline(Fin,Xin,Fout,Xout)
      integer                       :: i,Lin,Lout
      real(8),dimension(:)          :: Fin
      real(8),dimension(:)          :: Fout
      real(8),dimension(size(Fin))  :: Xin
      real(8),dimension(size(Fout)) :: Xout
      integer                       :: k
      real(8)                       :: x,y,y0,y1,x0,x1
      logical                       :: identical=.true.
      Lin =size(Fin) ; Lout=size(Fout)
      if(Lin==Lout)then
         do i=1,Lin
            identical=identical.AND.(Xin(i)==Xout(i))
         enddo
         if(identical)then
            Fout=Fin
            write(*,"(A)")"msg: SPLINE/linear_spline: Fout=Fin. Exit."
            return
         endif
      endif
      call interp_linear(1,Lin,Xin,Fin,Lout,Xout,Fout)
    end subroutine d_linear_spline

    subroutine c_linear_spline(Fin,Xin,Fout,Xout)
      integer                      :: i,Lin,Lout
      complex(8),dimension(:)      :: Fin
      complex(8),dimension(:)      :: Fout
      real(8),dimension(size(Fin)) :: Xin
      real(8),dimension(size(Fout)):: Xout,dummyR,dummyI
      Lin =size(Fin) ; Lout=size(Fout)
      call d_linear_spline(real(Fin,8),Xin,dummyR,Xout)
      call d_linear_spline(aimag(Fin),Xin,dummyI,Xout)
      Fout=cmplx(dummyR,dummyI,8)
    end subroutine c_linear_spline

    !*******************************************************************
    !*******************************************************************
    !*******************************************************************



    !+-------------------------------------------------------------------+
    !PROGRAM  : FSPLINE
    !TYPE     : Subroutine
    !PURPOSE  : Generic SPLINE routine with polinomial approximation
    ! still in the form of -L:L...
    !+-------------------------------------------------------------------+
    subroutine d_fspline(Fin,Xin,Fout,Xout,N)
      integer                       :: i,Lin,Lout
      integer,optional              :: N
      real(8),dimension(:)          :: Fin
      real(8),dimension(:)          :: Fout
      real(8),dimension(size(Fin))  :: Xin
      real(8),dimension(size(Fout)) :: Xout
      logical                       :: identical=.true.
      Lin =size(Fin) ; Lout=size(Fout)
      if(Lin==Lout)then
         do i=1,Lin
            identical=identical.AND.(Xin(i)==Xout(i))
         enddo
         if(identical)then
            Fout=Fin
            write(*,"(A)")"msg: SPLINE/poly_spline: Fout=Fin. Exit."
            return
         endif
      endif
      if(allocated(finterX))deallocate(finterX)
      if(allocated(finterF))deallocate(finterF)
      allocate(finterX(Lin),finterF(Lin))
      finterX    = Xin
      finterF    = Fin
      finterImax = Lin
      finterImin = 1
      finterN    = 5 ; if(present(N))finterN = N
      do i=1,Lout
         Fout(i)=finter(Xout(i))
      enddo
    end subroutine d_fspline

    subroutine c_fspline(Fin,Xin,Fout,Xout,N)
      integer                      :: i,Lin,Lout
      integer,optional             :: N
      complex(8),dimension(:)      :: Fin
      complex(8),dimension(:)      :: Fout
      real(8),dimension(size(Fin)) :: Xin
      real(8),dimension(size(Fout)):: Xout,dummyR,dummyI
      Lin =size(Fin) ; Lout=size(Fout)
      if(.not.present(N))then
         call d_fspline(real(Fin,8),Xin,dummyR,Xout)
         call d_fspline(aimag(Fin),Xin,dummyI,Xout)
      else
         call d_fspline(real(Fin,8),Xin,dummyR,Xout,N)
         call d_fspline(aimag(Fin),Xin,dummyI,Xout,N)
      endif
      Fout=cmplx(dummyR,dummyI,8)
    end subroutine c_fspline



    !*******************************************************************
    !*******************************************************************
    !*******************************************************************




    !+-------------------------------------------------------------------+
    !PURPOSE  : 
    !+-------------------------------------------------------------------+
    subroutine d_cub_interp(Fin,Xin,Fout,Xout)
      integer             :: i, j
      integer             :: Lin, Lout
      real(8)             :: Fin(:),   Xin(1:size(Fin))
      real(8)             :: Fout(:),  Xout(1:size(Fout))
      real(8)             :: xa(1:size(Fin)), ya(4,1:size(Fin))
      real(8)             :: x, y
      real(8),external    :: ppvalu
      logical             :: identical=.true.
      Lin = size(Fin) ; Lout= size(Fout)
      if(Lin==Lout)then
         do i=1,Lin
            identical=identical.AND.(Xin(i)==Xout(i))
         enddo
         if(identical)then
            Fout=Fin
            write(*,"(A)")"msg: SPLINE/cubic_spline: Fout=Fin. Exit."
            return
         endif
      endif
      xa(:)  = Xin(:)
      ya(1,:)= Fin(:)
      call CUBSPL(xa,ya,Lin,0,0)
      do i=1,Lout
         x = Xout(i)
         Fout(i) = PPVALU(xa,ya,Lin-1,4,x,0)
      enddo
      if(Xin(Lin) >= Xout(Lout))then
         Fout(Lout)=Fin(Lin)
      else
         Fout(Lout) = Fout(Lout-2) + &
              (Xout(Lout)-Xout(Lout-2))/(Xout(Lout-1)-Xout(Lout-2))*(Fout(Lout-1)-Fout(Lout-2))
      endif
    end subroutine d_cub_interp
    !+-------------------------------------------------------------------+

    !+-------------------------------------------------------------------+
    subroutine c_cub_interp(Fin,Xin,Fout,Xout)
      integer             :: i, j
      integer             :: Lin, Lout
      complex(8)          :: Fin(:),Fout(:)
      real(8),dimension(size(Fin)) :: reFin,imFin,Xin
      real(8),dimension(size(Fout)) :: reFout,imFout,Xout
      Lin = size(Fin)   ; Lout= size(Fout)
      reFin=real(Fin,8) ; imFin=aimag(Fin)
      call d_cub_interp(reFin,Xin,reFout,Xout)
      call d_cub_interp(imFin,Xin,imFout,Xout)
      Fout=cmplx(reFout,imFout,8)
    end subroutine c_cub_interp



    !*******************************************************************
    !*******************************************************************
    !*******************************************************************



    !+-------------------------------------------------------------------+
    !PURPOSE  :  
    !+-------------------------------------------------------------------+
    subroutine interp_Gtau(FG1,FG2,L1,L2)
      integer             :: L1, L2
      real(8)             :: FG1(0:L1), FG2(0:L2)
      real(8)             :: FctL1(-L1:L1), FctL2(-L2:L2)
      real(8),allocatable :: xa(:), ya(:,:)!, y2(:)
      integer             :: L11, L12, L13
      real(8)             :: x, y
      integer             :: i, j
      real(8),external    :: ppvalu

      FctL1(0:L1)=FG1 ; forall(i=1:L1)FctL1(-i)=-FctL1(L1-i)
      L11 = L1 + 1
      L12 = L1 + 2
      L13 = L1 + 3
      allocate(xa(L11),ya(4,L11))
      do i=1, L11
         xa(i)=dble(i-1)/dble(L1)
         ya(1,i)=FctL1(i-1)
      enddo
      call CUBSPL(xa,ya,L11,0,0)
      do i=1, L2
         x=dble(i)/dble(L2)
         FctL2(i)=PPVALU(xa,ya,L1,4,x,0)
      enddo
      FctL2(0)=FctL1(0)
      ! do i=1, L2
      !    FctL2(-i)=-FctL2(L2-i)
      ! enddo
      FG2=FctL2(0:L2)
    end subroutine interp_Gtau
    !*******************************************************************
    !*******************************************************************
    !*******************************************************************




    !+-------------------------------------------------------------------+
    !PURPOSE  : Sample a given function G(tau) over Nfak < N points.
    !COMMENTS : Incoming function is expected to have N+1 points (tau=0,beta)
    !this is modification with respect to the std extract routine used in
    !HFqmc.
    ! g0 has N+1 points
    ! g00 will have Nfak+1 (tau_fak=0,beta)
    !+-------------------------------------------------------------------+
    subroutine extract_gtau(g0,g00)
      real(8),dimension(0:)  :: g0 !0:L
      real(8),dimension(0:) :: g00 !0:Lfak
      integer :: N,Nfak
      integer :: i,ip
      real(8) :: p,mismatch
      N=size(g0)-1
      Nfak=size(g00)-1
      g00(0)=g0(0)
      if(g0(0) > 0.d0) g00(Nfak)=1.d0-g0(0)
      if(g0(0) < 0.d0) g00(Nfak)=-(g0(0)+1.d0)
      mismatch=dble(N)/dble(Nfak)
      do i=1,Nfak-1
         p=dble(i)*mismatch
         ip=int(p)
         g00(i)=g0(ip)
      enddo
    end subroutine extract_gtau



    !*******************************************************************
    !*******************************************************************
    !*******************************************************************


    ! !-----------------------    
    ! !-----------------------    
    ! !-----------------------    
    ! subroutine Spline(x,y,n,yp1,ypn,y2)
    !   implicit none
    !   !integer, parameter :: nmax = 4096*32
    !   integer      :: n
    !   real(8)      :: x(:), y(size(x)), y2(size(x)), yp1, ypn
    !   integer      :: i, k
    !   real(8)      :: p, qn, sig, un
    !   real(8),allocatable :: u(:)
    !   n=size(x);allocate(u(n))
    !   if (yp1 .gt. 0.99e30) then
    !      y2(1) = 0.0d0
    !      u(1)  = 0.0d0
    !   else
    !      y2(1) = -0.5d0
    !      u(1)  = (3.0d0/(x(2)-x(1))) * ((y(2)-y(1))/(x(2)-x(1)) - yp1)
    !   endif
    !   do i=2, n-1
    !      sig   = (x(i)-x(i-1))/(x(i+1)-x(i-1))
    !      p     = sig * y2(i-1) + 2.0d0
    !      y2(i) = (sig-1.0d0)/p
    !      u(i)  = (6.0d0 * ((y(i+1)-y(i))/(x(i+1)-x(i)) - (y(i)-y(i-1))/(x(i)-x(i-1)))/(x(i+1)-x(i-1)) - sig*u(i-1)) /p
    !   enddo
    !   if (ypn .gt. 0.99e30) then
    !      qn = 0.0d0
    !      un = 0.0d0
    !   else
    !      qn = 0.5d0
    !      un = (3.0d0/(x(n)-x(n-1))) * (ypn - (y(n)-y(n-1))/(x(n)-x(n-1)))
    !   endif
    !   y2(n) = (un - qn*u(n-1))/(qn*y2(n-1) + 1.0d0)
    !   do k=n-1,1,-1
    !      y2(k) = y2(k)*y2(k+1) + u(k)
    !   enddo
    ! end subroutine Spline
    ! !-----------------------    
    ! !-----------------------    
    ! !-----------------------    
    ! subroutine Splint(xa,ya,y2a,n,x,y)
    !   implicit none
    !   integer      :: n
    !   real(kind=8) :: xa(n), ya(n), y2a(n), x, y
    !   integer      :: k, khi, klo
    !   real(kind=8) :: a, b, h
    !   klo=1
    !   khi=n
    !   do while(khi-klo .gt. 1)
    !      k = (khi+klo)/2
    !      if(xa(k) .gt. x)then
    !         khi = k
    !      else
    !         klo = k
    !      endif
    !   enddo
    !   h = xa(khi)-xa(klo)
    !   if (h .eq. 0.0d0) then
    !      write(*,*) "ERROR: Subroutine 'Splint'"
    !      write(*,*) "       Bad xa input"
    !      stop
    !   end if
    !   a = (xa(khi)-x)/h
    !   b = (x-xa(klo))/h
    !   y = a*ya(klo) + b*ya(khi) + ((a**3-a)*y2a(klo) + (b**3-b)*y2a(khi))*(h**2)/6.0d0
    ! end subroutine Splint


  end module SPLINE
