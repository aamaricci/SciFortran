  include "SPLINE_NR_MOD.f90"
  include "SPLINE_FINTER_MOD.f90"
  include "interp.f90"
  module SPLINE
    USE SPLINE_NR_MOD
    USE SPLINE_FINTER_MOD
    !###############################################################
    !     PROGRAM  : SPLINE
    !     TYPE     : Module
    !     PURPOSE  : provide an interface to interpolation routines
    ! Now I don't have time to add more stuff, that can be  found 
    ! in splinF directory as copied from internet.
    !     AUTHORS  : 
    !###############################################################
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

      integer :: k
      real(8) :: x,y,y0,y1,x0,x1
      Lin =size(Fin) ; Lout=size(Fout)
      if(Lin==Lout)then
         Fout=Fin
         return
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
      Lin =size(Fin) ; Lout=size(Fout)
      if(Lin==Lout)then
         Fout=Fin
         return
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
      Lin = size(Fin) ; Lout= size(Fout)
      if(Lin==Lout)then
         Fout=Fin
         return
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


    include "gtau_funx.f90"


    include "cubspl_routines.f90"


  end module SPLINE
