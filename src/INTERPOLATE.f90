  !###############################################################
  ! PURPOSE  : provide an interface to interpolation routines
  ! Now I don't have time to add more stuff, that can be  found 
  ! in splinF directory as copied from internet.
  !###############################################################
  include "interpolate_nr.f90"
  include "interpolate_pack.f90"
  include "interpolate_pppack.f90"
  module INTERPOLATE
    USE INTERPOLATE_NR
    implicit none
    private

    type finter_type
       real(8),allocatable :: X(:)
       real(8),allocatable :: F(:)
       integer             :: Imin=0,Imax=0,N=0
       logical             :: status=.false.
    end type finter_type

    type finter2d_type
       real(8),allocatable :: X(:)  !vector with frequencies
       real(8),allocatable :: Y(:)  !vector with frequencies
       real(8),allocatable :: F(:,:) !corresponding vector of functional values
       integer             :: N=0
       integer             :: Imin=0,Imax=0,Jmin=0,Jmax=0
       logical             :: status=.false.
    end type finter2d_type

    interface linear_spline
       module procedure &
            d_linear_spline_s,&
            d_linear_spline_v,&
            c_linear_spline_s,&
            c_linear_spline_v
    end interface linear_spline


    interface poly_spline
       module procedure &
            d_poly_spline_s,&
            d_poly_spline_v,&
            c_poly_spline_s,&
            c_poly_spline_v
    end interface poly_spline

    interface cubic_spline
       module procedure  &
            d_cub_interp_s,&
            d_cub_interp_v,&
            c_cub_interp_s,&
            c_cub_interp_v
    end interface cubic_spline


    interface linear_spline2d
       module procedure &
            d_linear_spline_2d_s,&
            d_linear_spline_2d_v,&
            c_linear_spline_2d_s,&
            c_linear_spline_2d_v
    end interface linear_spline2d


    interface poly_spline2d
       module procedure &
            d_poly_spline_2d_s,&
            d_poly_spline_2d_v,&
            c_poly_spline_2d_s,&
            c_poly_spline_2d_v
    end interface poly_spline2d


    interface extract
       module procedure extract_gtau
    end interface extract

    !AVAILABLE FROM NR
    public :: locate !binary search:
    public :: polint
    public :: polin2

    !AVAILABLE FROM INTERPOLATE
    !1-dimension
    public :: finter_type
    public :: init_finter
    public :: delete_finter
    public :: finter
    public :: poly_spline
    public :: cubic_spline
    public :: linear_spline

    !2-dimension
    public :: finter2d_type
    public :: init_finter2d
    public :: delete_finter2d
    public :: finter2d
    public :: linear_spline2d
    public :: poly_spline2d

    !LEGACY CODE:
    !used in some determinant qmc algorithm. to be removed
    public :: extract,extract_gtau 
    public :: interp_gtau


  contains

    !*******************************************************************
    ! 1-DIMENSIONAL SPLINES:
    !*******************************************************************
    !+-------------------------------------------------------------------+
    !PURPOSE  : Linear interpolation of given data
    !+-------------------------------------------------------------------+
    subroutine d_linear_spline_s(Xin,Fin,Xout,Fout)
      integer                      :: i,Lin,Lout
      real(8),dimension(:)         :: Xin
      real(8),dimension(size(Xin)) :: Fin
      real(8)                      :: Xout
      real(8)                      :: Fout
      real(8),dimension(1)         :: xout_,fout_
      integer                      :: k
      real(8)                      :: x,y,y0,y1,x0,x1
      logical                      :: identical=.true.
      Lin =size(Fin)
      xout_(1) = xout
      call interp_linear(1,Lin,Xin,Fin,1,xout_,fout_)
      Fout = fout_(1)
    end subroutine d_linear_spline_s
    !+-------------------------------------------------------------------+
    subroutine d_linear_spline_v(Xin,Fin,Xout,Fout)
      integer                       :: i,Lin,Lout
      real(8),dimension(:)          :: Xin
      real(8),dimension(size(Xin))  :: Fin
      real(8),dimension(:)          :: Xout
      real(8),dimension(size(Xout)) :: Fout
      integer                       :: k
      real(8)                       :: x,y,y0,y1,x0,x1
      logical                       :: identical=.true.
      Lin =size(Fin) ; Lout=size(Fout)
      if(Lin==Lout)call test_grid_equality_d(xin,xout,fin,fout)
      call interp_linear(1,Lin,Xin,Fin,Lout,Xout,Fout)
    end subroutine d_linear_spline_v
    !+-------------------------------------------------------------------+
    subroutine c_linear_spline_s(Xin,Fin,Xout,Fout)
      integer                         :: i,Lin
      real(8),dimension(:)            :: Xin
      complex(8),dimension(size(Xin)) :: Fin
      real(8)                         :: Xout
      complex(8)                      :: Fout
      real(8)                         :: ref_,imf_
      Lin = size(Fin)
      call d_linear_spline_s(Xin,dreal(Fin),xout,ref_)
      call d_linear_spline_s(Xin,dimag(Fin),xout,imf_)
      Fout=cmplx(ref_,imf_,8)
    end subroutine c_linear_spline_s
    !+-------------------------------------------------------------------+
    subroutine c_linear_spline_v(Xin,Fin,Xout,Fout)
      integer                          :: i,Lin,Lout
      real(8),dimension(:)             :: Xin
      complex(8),dimension(size(Xin))  :: Fin
      real(8),dimension(:)             :: Xout
      complex(8),dimension(size(Xout)) :: Fout
      real(8),dimension(size(Xout))    :: dummyR,dummyI
      logical                          :: identical=.true.
      Lin =size(Xin) ; Lout=size(Xout)
      if(Lin==Lout)call test_grid_equality_c(xin,xout,fin,fout)
      call d_linear_spline_v(Xin,dreal(Fin),Xout,dummyR)
      call d_linear_spline_v(Xin,dimag(Fin),Xout,dummyI)
      Fout=cmplx(dummyR,dummyI,8)
    end subroutine c_linear_spline_v



    !*******************************************************************
    !*******************************************************************
    !*******************************************************************



    !+-------------------------------------------------------------------+
    !PURPOSE  : Generic SPLINE routine with polinomial approximation
    ! uses a local function to pass exactly N points to the polynomial
    ! interpolation routines, fixing the order of the approximation
    !HINT: N>20 is not reccomended.
    !+-------------------------------------------------------------------+
    subroutine d_poly_spline_s(Xin,Fin,Xout,Fout,N)
      integer                      :: i,Lin,Lout,N_
      integer,optional             :: N
      real(8),dimension(:)         :: Xin
      real(8),dimension(size(Xin)) :: Fin
      real(8)                      :: Xout
      real(8)                      :: Fout
      type(finter_type)                :: pfinter
      Lin =size(Fin)
      N_    = 5 ; if(present(N))N_ = N
      call init_finter(pfinter,Xin,Fin,N_)
      Fout = finter(pfinter,Xout)
      call delete_finter(pfinter)
    end subroutine d_poly_spline_s
    !+-------------------------------------------------------------------+
    subroutine d_poly_spline_v(Xin,Fin,Xout,Fout,N)
      integer                       :: i,Lin,Lout,N_
      integer,optional              :: N
      real(8),dimension(:)          :: Xin
      real(8),dimension(size(Xin))  :: Fin
      real(8),dimension(:)          :: Xout
      real(8),dimension(size(Xout)) :: Fout
      type(finter_type)                 :: pfinter
      Lin =size(Fin) ; Lout=size(Fout)
      if(Lin==Lout)call test_grid_equality_d(xin,xout,fin,fout)
      N_    = 5 ; if(present(N))N_ = N
      call init_finter(pfinter,Xin,Fin,N_)
      do i=1,Lout
         Fout(i)=finter(pfinter,Xout(i))
      enddo
      call delete_finter(pfinter)
    end subroutine d_poly_spline_v
    !+-------------------------------------------------------------------+
    subroutine c_poly_spline_s(Xin,Fin,Xout,Fout,N)
      integer                         :: i,Lin,Lout
      integer,optional                :: N
      real(8),dimension(:)            :: Xin
      complex(8),dimension(size(Xin)) :: Fin
      real(8)                         :: Xout
      complex(8)                      :: Fout
      real(8)                         :: dummyR,dummyI
      Lin =size(Fin)
      if(.not.present(N))then
         call d_poly_spline_s(Xin,dreal(Fin),Xout,dummyR)
         call d_poly_spline_s(Xin,dimag(Fin),Xout,dummyI)
      else
         call d_poly_spline_s(Xin,dreal(Fin),Xout,dummyR,N)
         call d_poly_spline_s(Xin,dimag(Fin),Xout,dummyI,N)
      endif
      Fout=cmplx(dummyR,dummyI,8)
    end subroutine c_poly_spline_s
    !+-------------------------------------------------------------------+
    subroutine c_poly_spline_v(Xin,Fin,Xout,Fout,N)
      integer                          :: i,Lin,Lout
      integer,optional                 :: N
      real(8),dimension(:)             :: Xin
      complex(8),dimension(size(Xin))  :: Fin
      real(8),dimension(:)             :: Xout
      complex(8),dimension(size(Xout)) :: Fout
      real(8),dimension(size(Fout))    :: dummyR,dummyI
      Lin =size(Fin) ; Lout=size(Fout)
      if(Lin==Lout)call test_grid_equality_c(xin,xout,fin,fout)
      if(.not.present(N))then
         call d_poly_spline_v(Xin,dreal(Fin),Xout,dummyR)
         call d_poly_spline_v(Xin,dimag(Fin),Xout,dummyI)
      else
         call d_poly_spline_v(Xin,dreal(Fin),Xout,dummyR,N)
         call d_poly_spline_v(Xin,dimag(Fin),Xout,dummyI,N)
      endif
      Fout=cmplx(dummyR,dummyI,8)
    end subroutine c_poly_spline_v



    !*******************************************************************
    !*******************************************************************
    !*******************************************************************



    !+-------------------------------------------------------------------+
    !PURPOSE  : Interface to the DeBoer Cubic Spline routine.
    ! this routine works great with imaginary time GF
    !+-------------------------------------------------------------------+
    subroutine d_cub_interp_s(Xin,Fin,Xout,Fout)
      integer                      :: i, j
      integer                      :: Lin, Lout
      real(8),dimension(:)         :: Xin
      real(8),dimension(size(Xin)) :: Fin
      real(8)                      :: Xout
      real(8)                      :: Fout
      real(8)                      :: xa(1:size(Fin)), ya(4,1:size(Fin))
      real(8)                      :: x, y
      real(8),external             :: ppvalu
      Lin = size(Fin) 
      xa(:)  = Xin(:)
      ya(1,:)= Fin(:)
      call CUBSPL(xa,ya,Lin,0,0)
      x = Xout
      Fout = PPVALU(xa,ya,Lin-1,4,x,0)
      if(Xin(Lin) >= Xout)Fout=Fin(Lin)
    end subroutine d_cub_interp_s
    !+-------------------------------------------------------------------+
    subroutine d_cub_interp_v(Xin,Fin,Xout,Fout)
      integer                       :: i, j
      integer                       :: Lin, Lout
      real(8),dimension(:)          :: Xin
      real(8),dimension(size(Xin))  :: Fin
      real(8),dimension(:)          :: Xout
      real(8),dimension(size(Xout)) :: Fout
      real(8)                       :: xa(1:size(Fin)), ya(4,1:size(Fin))
      real(8)                       :: x, y
      real(8),external              :: ppvalu
      Lin = size(Fin) ; Lout= size(Fout)
      if(Lin==Lout)call test_grid_equality_d(xin,xout,fin,fout)
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
    end subroutine d_cub_interp_v
    !+-------------------------------------------------------------------+
    subroutine c_cub_interp_s(Xin,Fin,Xout,Fout)
      integer                         :: i, j
      integer                         :: Lin, Lout
      real(8),dimension(:)            :: Xin
      complex(8),dimension(size(Xin)) :: Fin
      real(8)                         :: Xout
      complex(8)                      :: Fout
      real(8)                         :: reFout,imFout
      Lin = size(Fin)
      call d_cub_interp_s(Xin,dreal(Fin),Xout,reFout)
      call d_cub_interp_s(Xin,dimag(Fin),Xout,imFout)
      Fout=cmplx(reFout,imFout,8)
    end subroutine c_cub_interp_s
    !+-------------------------------------------------------------------+
    subroutine c_cub_interp_v(Xin,Fin,Xout,Fout)
      integer                          :: i, j
      integer                          :: Lin, Lout
      real(8),dimension(:)             :: Xin
      complex(8),dimension(size(Xin))  :: Fin
      real(8),dimension(:)             :: Xout
      complex(8),dimension(size(Xout)) :: Fout
      real(8),dimension(size(Xout))    :: reFout,imFout
      Lin = size(Fin)   ; Lout= size(Fout)
      if(Lin==Lout)call test_grid_equality_c(xin,xout,fin,fout)
      call d_cub_interp_v(Xin,dreal(Fin),Xout,reFout)
      call d_cub_interp_v(Xin,dimag(Fin),Xout,imFout)
      Fout=cmplx(reFout,imFout,8)
    end subroutine c_cub_interp_v






    !*******************************************************************
    ! 2-DIMENSIONAL SPLINES:
    !*******************************************************************
    !+-------------------------------------------------------------------+
    !PURPOSE  : This function uses bilinear interpolation to estimate 
    ! the value of a function f at point (x,y).
    ! f is assumed to be sampled on a regular grid, with the grid x values specified
    ! by x_array and the grid y values specified by y_array
    ! Reference: http://en.wikipedia.org/wiki/Bilinear_interpolation
    !+-------------------------------------------------------------------+
    subroutine d_linear_spline_2d_s(xin,yin,fin,xout,yout,fout)
      real(8),dimension(:)                   :: xin
      real(8),dimension(:)                   :: yin
      real(8),dimension(size(xin),size(yin)) :: fin
      integer                                :: Lxin,Lyin
      real(8)                                :: xout
      real(8)                                :: yout
      real(8)                                :: fout
      fout = bilinear_interpolate(xin,yin,fin,xout,yout)
    end subroutine d_linear_spline_2d_s
    !+-------------------------------------------------------------------+!
    subroutine d_linear_spline_2d_v(xin,yin,fin,xout,yout,fout)
      real(8),dimension(:)                     :: xin
      real(8),dimension(:)                     :: yin
      real(8),dimension(size(xin),size(yin))   :: fin
      integer                                  :: Lxin,Lyin
      real(8),dimension(:)                     :: xout
      real(8),dimension(:)                     :: yout
      real(8),dimension(size(xout),size(yout)) :: fout
      integer                                  :: Lxout,Lyout
      integer                                  :: ix,iy
      real(8)                                  :: x,y
      Lxin = size(xin) ; Lyin = size(yin)
      Lxout= size(xout); Lyout=size(yout)
      do ix=1,Lxout
         x = xout(ix)
         do iy=1,Lyout
            y = yout(iy)
            fout(ix,iy) = bilinear_interpolate(xin,yin,fin,x,y)
         enddo
      enddo
    end subroutine d_linear_spline_2d_v
    !+-------------------------------------------------------------------+
    subroutine c_linear_spline_2d_s(xin,yin,fin,xout,yout,fout)
      real(8),dimension(:)                      :: xin
      real(8),dimension(:)                      :: yin
      complex(8),dimension(size(xin),size(yin)) :: fin
      integer                                   :: Lxin,Lyin
      real(8)                                   :: xout
      real(8)                                   :: yout
      complex(8)                                :: fout
      real(8)                                   :: ref,imf
      ref = bilinear_interpolate(xin,yin,dreal(fin),xout,yout)
      imf = bilinear_interpolate(xin,yin,dimag(fin),xout,yout)
      fout=cmplx(ref,imf,8)
    end subroutine c_linear_spline_2d_s
    !+-------------------------------------------------------------------+
    subroutine c_linear_spline_2d_v(xin,yin,fin,xout,yout,fout)
      real(8),dimension(:)                     :: xin
      real(8),dimension(:)                     :: yin
      complex(8),dimension(size(xin),size(yin))   :: fin
      integer                                  :: Lxin,Lyin
      real(8),dimension(:)                     :: xout
      real(8),dimension(:)                     :: yout
      complex(8),dimension(size(xout),size(yout)) :: fout
      integer                                  :: Lxout,Lyout
      integer                                  :: ix,iy
      real(8)                                  :: x,y,ref,imf
      Lxin = size(xin) ; Lyin = size(yin)
      Lxout= size(xout); Lyout=size(yout)
      do ix=1,Lxout
         x = xout(ix)
         do iy=1,Lyout
            y = yout(iy)
            ref = bilinear_interpolate(xin,yin,dreal(fin),x,y)
            imf = bilinear_interpolate(xin,yin,dimag(fin),x,y)
            fout(ix,iy) = cmplx(ref,imf,8)
         enddo
      enddo
    end subroutine c_linear_spline_2d_v




    !*******************************************************************
    !*******************************************************************
    !*******************************************************************




    !+-------------------------------------------------------------------+
    !PURPOSE  : Polynomial interpolation of 2D surface using 
    ! local plaquette (order N**2).
    !+-------------------------------------------------------------------+
    subroutine d_poly_spline_2d_s(xin,yin,fin,xout,yout,fout,N)
      real(8),dimension(:)                   :: xin
      real(8),dimension(:)                   :: yin
      real(8),dimension(size(xin),size(yin)) :: fin
      integer                                :: Lxin,Lyin
      integer,optional                       :: N
      integer                                :: N_
      real(8)                                :: xout
      real(8)                                :: yout
      real(8)                                :: fout
      type(finter2d_type)                    :: pfinter
      Lxin = size(xin)
      N_    = 5 ; if(present(N))N_ = N
      call init_finter2d(pfinter,Xin,Yin,Fin,N_)
      fout = finter2d(pfinter,xout,yout)
      call delete_finter2d(pfinter)
    end subroutine d_poly_spline_2d_s
    !+-------------------------------------------------------------------+
    subroutine d_poly_spline_2d_v(xin,yin,fin,xout,yout,fout,N)
      real(8),dimension(:)                     :: xin
      real(8),dimension(:)                     :: yin
      real(8),dimension(size(xin),size(yin))   :: fin
      integer                                  :: Lxin,Lyin
      integer,optional                         :: N
      integer                                  :: N_
      real(8),dimension(:)                     :: xout
      real(8),dimension(:)                     :: yout
      real(8),dimension(size(xout),size(yout)) :: fout
      integer                                  :: Lxout,Lyout
      integer                                  :: ix,iy
      real(8)                                  :: x,y
      type(finter2d_type)                      :: pfinter
      Lxin = size(xin) ; Lyin = size(yin)
      Lxout= size(xout); Lyout=size(yout)
      if(Lxin==Lxout.AND.Lyin==Lyout)call test_grid_equality_d(xin,xout,fin,fout)
      N_    = 5 ; if(present(N))N_ = N
      call init_finter2d(pfinter,Xin,Yin,Fin,N_)
      do ix=1,Lxout
         x = xout(ix)
         do iy=1,Lyout
            y = yout(iy)
            fout(ix,iy) = finter2d(pfinter,x,y)
         enddo
      enddo
      call delete_finter2d(pfinter)
    end subroutine d_poly_spline_2d_v
    !+-------------------------------------------------------------------+
    subroutine c_poly_spline_2d_s(xin,yin,fin,xout,yout,fout,N)
      real(8),dimension(:)                      :: xin
      real(8),dimension(:)                      :: yin
      complex(8),dimension(size(xin),size(yin)) :: fin
      integer                                   :: Lxin,Lyin
      integer,optional                          :: N
      integer                                   :: N_
      real(8)                                   :: xout
      real(8)                                   :: yout
      complex(8)                                :: fout
      real(8)                                   :: ref,imf
      if(.not.present(N))then
         call d_poly_spline_2d_s(Xin,Yin,dreal(Fin),Xout,Yout,ref)
         call d_poly_spline_2d_s(Xin,Yin,dimag(Fin),Xout,Yout,imf)
      else
         call d_poly_spline_2d_s(Xin,Yin,dreal(Fin),Xout,Yout,ref,N)
         call d_poly_spline_2d_s(Xin,Yin,dimag(Fin),Xout,Yout,imf,N)
      endif
      Fout=cmplx(ref,imf,8)
    end subroutine c_poly_spline_2d_s
    !+-------------------------------------------------------------------+
    subroutine c_poly_spline_2d_v(xin,yin,fin,xout,yout,fout,N)
      real(8),dimension(:)                        :: xin
      real(8),dimension(:)                        :: yin
      complex(8),dimension(size(xin),size(yin))   :: fin
      integer                                     :: Lxin,Lyin
      integer,optional                            :: N
      integer                                     :: N_
      real(8),dimension(:)                        :: xout
      real(8),dimension(:)                        :: yout
      complex(8),dimension(size(xout),size(yout)) :: fout
      real(8),dimension(size(xout),size(yout))    :: ref,imf
      integer                                     :: Lxout,Lyout
      Lxin = size(xin) ; Lyin = size(yin)
      Lxout= size(xout); Lyout=size(yout)
      if(Lxin==Lxout.AND.Lyin==Lyout)call test_grid_equality_c(xin,xout,fin,fout)
      if(.not.present(N))then
         call d_poly_spline_2d_v(Xin,Yin,dreal(Fin),Xout,Yout,ref)
         call d_poly_spline_2d_v(Xin,Yin,dimag(Fin),Xout,Yout,imf)
      else
         call d_poly_spline_2d_v(Xin,Yin,dreal(Fin),Xout,Yout,ref,N)
         call d_poly_spline_2d_v(Xin,Yin,dimag(Fin),Xout,Yout,imf,N)
      endif
      Fout=cmplx(ref,imf,8)
    end subroutine c_poly_spline_2d_v





    !*******************************************************************
    !*******************************************************************
    !*******************************************************************






    ! subroutine d_poly_fspline_2d_v(xin,yin,fin,xout,yout,fout,N)
    !   real(8),dimension(:)                     :: xin
    !   real(8),dimension(:)                     :: yin
    !   real(8),dimension(size(xin),size(yin))   :: fin
    !   integer                                  :: Lxin,Lyin
    !   real(8),dimension(:)                     :: xout
    !   real(8),dimension(:)                     :: yout
    !   real(8),dimension(size(xout),size(yout)) :: fout
    !   integer                                  :: Lxout,Lyout
    !   integer :: N
    !   integer                                  :: ix,iy
    !   real(8),dimension(:,:),allocatable       :: fxtmp
    !   Lxin = size(xin) ; Lyin = size(yin)
    !   Lxout= size(xout); Lyout=size(yout)
    !   allocate(fxtmp(Lxout,Lyin))
    !   do iy=1,Lyin              !loop sulle righe griglia-vecchia
    !      call poly_spline(xin,fin(:,iy),xout,fxtmp(:,iy),N) !spline sulle colonne, dell'ix-esima riga
    !   enddo
    !   do ix=1,Lxout
    !      call poly_spline(yin,fxtmp(ix,:),yout,fout(ix,:),N)
    !   enddo
    ! end subroutine  d_poly_fspline_2d_v
    ! !+-------------------------------------------------------------------+
    ! subroutine c_poly_fspline_2d_v(xin,yin,fin,xout,yout,fout,N)
    !   real(8),dimension(:)                        :: xin
    !   real(8),dimension(:)                        :: yin
    !   complex(8),dimension(size(xin),size(yin))   :: fin
    !   integer :: N
    !   integer                                     :: Lxin,Lyin
    !   real(8),dimension(:)                        :: xout
    !   real(8),dimension(:)                        :: yout
    !   complex(8),dimension(size(xout),size(yout)) :: fout
    !   integer                                     :: Lxout,Lyout
    !   real(8),dimension(:,:),allocatable          :: reF,imF
    !   Lxin = size(xin) ; Lyin = size(yin)
    !   Lxout= size(xout); Lyout=size(yout)
    !   allocate(reF(Lxout,Lyout),imF(Lxout,Lyout))
    !   call d_poly_spline_2d_v(xin,yin,dreal(fin),xout,yout,reF,N)
    !   call d_poly_spline_2d_v(xin,yin,dimag(fin),xout,yout,imF,N)
    !   fout=cmplx(reF,imF,8)
    !   deallocate(reF,imF)
    ! end subroutine  c_poly_fspline_2d_v






    !*******************************************************************
    !*******************************************************************
    !*******************************************************************






    !*******************************************************************
    ! computational ROUTINES:
    !*******************************************************************
    subroutine init_finter(func,Xin,Fin,N)
      type(finter_type) :: func
      real(8)       :: xin(:)
      real(8)       :: fin(size(xin))
      integer       :: N,Lin
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
    !+-------------------------------------------------------------------+
    subroutine init_finter2d(func,Xin,Yin,Fin,N)
      type(finter2d_type) :: func
      real(8)       :: xin(:),yin(:)
      real(8)       :: fin(size(xin),size(yin))
      integer       :: N,Lx,Ly
      if(func%status)deallocate(func%x,func%y,func%f)
      Lx=size(xin) ; Ly=size(yin)
      allocate(func%x(Lx),func%y(Ly),func%f(Lx,Ly))
      func%X    = Xin
      func%Y    = Yin
      func%F    = Fin
      func%Imin = 1
      func%Jmin = 1
      func%Imax = Lx
      func%Jmax = Ly
      func%N    = N
      func%status=.true.
    end subroutine init_finter2d


    !*******************************************************************
    !*******************************************************************
    !*******************************************************************


    subroutine delete_finter(func)
      type(finter_type) :: func
      if(allocated(func%x))deallocate(func%x)
      if(allocated(func%f))deallocate(func%f)
      func%imax=0
      func%imin=0
      func%N   =0
      func%status=.false.
    end subroutine delete_finter
    !+-------------------------------------------------------------------+
    subroutine delete_finter2d(func)
      type(finter2d_type) :: func
      if(allocated(func%x))deallocate(func%x)
      if(allocated(func%y))deallocate(func%y)
      if(allocated(func%f))deallocate(func%f)
      func%imin=0
      func%jmin=0
      func%imax=0
      func%jmax=0
      func%N   =0
      func%status=.false.
    end subroutine delete_finter2d


    !*******************************************************************
    !*******************************************************************
    !*******************************************************************


    function finter(func,x)
      real(8)           :: x
      type(finter_type) :: func
      real(8)           :: finter
      real(8)           :: y,dy
      integer           :: j,k,k0,k1
      integer           :: n
      N=func%N    !order of polynomial interpolation
      finter=0.d0
      j=locate(func%X(func%Imin:func%Imax),x)
      !k = min(max(j-(N-1)/2,1),func%Imax+1-N)
      k=max(j-(N-1)/2,1)
      k0=k
      if(k0 < func%Imin)k0=func%Imin
      k1=k+N+1
      if(k1 > func%Imax)then
         k1=func%Imax
         k0=k1-N-1
      endif
      call polint(func%X(k0:k1),func%F(k0:k1),x,y,dy)
      !call polint(func%X(k:k+n),func%F(k:k+n),x,y,dy)
      finter=y
      return
    end function finter



    !*******************************************************************
    !*******************************************************************
    !*******************************************************************



    function finter2d(func,x,y)
      real(8)         :: x,y
      type(finter2d_type) :: func
      real(8)         :: finter2d
      real(8)         :: f,df
      integer         :: itmp,jtmp,kx,ky,k0x,k0y,k1x,k1y
      integer         :: n
      N=func%N    !order of polynomial interpolation
      finter2d=0.d0
      itmp=locate(func%X(func%Imin:func%Imax),x)
      jtmp=locate(func%Y(func%Jmin:func%Jmax),y)
      kx=max(itmp-(N-1)/2,1)
      ky=max(jtmp-(N-1)/2,1)
      k0x = kx ; if(k0x < func%Imin)k0x=func%Imin
      k0y = ky ; if(k0y < func%Jmin)k0y=func%Jmin
      k1x = kx+N+1
      if(k1x > func%Imax)then         
         k1x=func%Imax
         k0x=k1x-N-1
      endif
      k1y = ky+N+1
      if(k1y > func%Jmax)then
         k1y=func%Jmax
         k0y=k1y-N-1
      endif
      call polin2(func%X(k0x:k1x),func%Y(k0y:k1y),func%F(k0x:k1x,k0y:k1y),x,y,f,df)
      finter2d=f
      return
    end function finter2d


    !*******************************************************************
    !*******************************************************************
    !*******************************************************************


    !+-------------------------------------------------------------------+
    !PURPOSE  :  test equality of grids and functions
    !+-------------------------------------------------------------------+
    subroutine test_grid_equality_d(xin,xout,Fin,Fout)
      real(8),dimension(:)         :: xin,xout
      real(8),dimension(size(xin)) :: fin,fout
      logical                      :: equal=.true.
      integer                      :: i,L
      L=size(xin)
      equal=.true.
      do i=1,L
         equal=equal.AND.(Xin(i)==Xout(i))
      enddo
      if(equal)then
         Fout=Fin
         write(*,"(A)")"msg: test_grid_equlity: Fout=Fin. Exit."
         return
      endif
    end subroutine test_grid_equality_d
    !
    subroutine test_grid_equality_c(xin,xout,Fin,Fout)
      real(8),dimension(:)            :: xin,xout
      complex(8),dimension(size(xin)) :: fin,fout
      logical                         :: equal=.true.
      integer                         :: i,L
      L=size(xin)
      equal=.true.
      do i=1,L
         equal=equal.AND.(Xin(i)==Xout(i))
      enddo
      if(equal)then
         Fout=Fin
         write(*,"(A)")"msg: test_grid_equlity: Fout=Fin. Exit."
         return
      endif
    end subroutine test_grid_equality_c


    !*******************************************************************
    !*******************************************************************
    !*******************************************************************


    function bilinear_interpolate(xgrid,ygrid,fgrid,x,y,delta) result(finterpolate)
      real(8), dimension(:),intent(in)                        :: xgrid
      real(8), dimension(:),intent(in)                        :: ygrid
      real(8), dimension(size(xgrid),size(ygrid)), intent(in) :: fgrid
      real(8), intent(in)                                     :: x,y
      real(8), intent(in), optional                           :: delta
      real(8)                                                 :: denom, x1, x2, y1, y2
      real(8)                                                 :: Q11,Q12,Q21,Q22
      real(8) :: finterpolate
      integer                                                 :: i,j
      !Get the indices of the bottom-leftmost point in the grid near the givem (x,y)
      i = locate(xgrid, x)!binary_search(xgrid, x)
      j = locate(ygrid, y)!binary_search(ygrid, y)
      x1 = xgrid(i)
      x2 = xgrid(i+1)
      y1 = ygrid(j)
      y2 = ygrid(j+1)
      denom = (x2-x1)*(y2-y1)
      Q11 = fgrid(i  ,j  )
      Q21 = fgrid(i+1,j  ) 
      Q22 = fgrid(i+1,j+1)
      Q12 = fgrid(i  ,j+1)    
      finterpolate = (&
           Q11*(x2-x)*(y2-y) + Q21*(x-x1)*(y2-y) + &
           Q12*(x2-x)*(y-y1) + Q22*(x-x1)*(y-y1) )/denom
    end function bilinear_interpolate








    !*******************************************************************
    ! LEGACY CODE USED IN THE DETERMINANTAL QMC. TO BE REMOVED
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



    ! !+-------------------------------------------------------------------+
    ! !PURPOSE  : Interface to NR polint. This is an alternative to the 
    ! ! poly_spline above. This one try to use all the available points
    ! ! in the in-grids
    ! !+-------------------------------------------------------------------+
    ! subroutine poly_spline_d(xin,fin,xout,fout)
    !   real(8),dimension(:)          :: xin
    !   real(8),dimension(size(xin))  :: fin
    !   real(8)                       :: xout
    !   real(8)                       :: fout
    !   integer                       :: i
    !   real(8)                       :: x,y,dy
    !   call polint(xin,fin,xout,fout,dy)
    ! end subroutine poly_spline_d
    ! !
    ! subroutine poly_spline_dv(xin,fin,xout,fout)
    !   real(8),dimension(:)          :: xin
    !   real(8),dimension(size(xin))  :: fin
    !   integer                       :: Lin
    !   real(8),dimension(:)          :: xout
    !   real(8),dimension(size(xout)) :: fout
    !   integer                       :: Lout
    !   integer                       :: i,j
    !   real(8)                       :: x,y,dy
    !   Lout= size(xout)
    !   do i=1,Lout
    !      x = xout(i)
    !      call polint(xin,fin,x,y,dy)
    !      fout(i)=y
    !   enddo
    ! end subroutine poly_spline_dv
    ! !
    ! subroutine poly_spline_c(xin,fin,xout,fout)
    !   real(8),dimension(:)            :: xin
    !   complex(8),dimension(size(xin)) :: fin
    !   real(8)                         :: xout
    !   complex(8)                      :: fout
    !   integer                         :: i,j
    !   real(8)                         :: x,rey,imy,dy
    !   call polint(xin,dreal(fin),xout,rey,dy)
    !   call polint(xin,dimag(fin),xout,imy,dy)
    !   fout=cmplx(rey,imy,8)
    ! end subroutine poly_spline_c
    ! !
    ! subroutine poly_spline_cv(xin,fin,xout,fout)
    !   real(8),dimension(:)             :: xin
    !   complex(8),dimension(size(xin))  :: fin
    !   real(8),dimension(:)             :: xout
    !   complex(8),dimension(size(xout)) :: fout
    !   integer                          :: Lout
    !   integer                          :: i,j
    !   real(8)                          :: x,rey,imy,dy
    !   Lout= size(xout)
    !   do i=1,Lout
    !      x = xout(i)
    !      call polint(xin,dreal(fin),x,rey,dy)
    !      call polint(xin,dimag(fin),x,imy,dy)
    !      fout(i)=cmplx(rey,imy,8)
    !   enddo
    ! end subroutine poly_spline_cv

  end module INTERPOLATE
