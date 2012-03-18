  !###############################################################
  ! PROGRAM  : TOOLS
  ! TYPE     : Module
  ! PURPOSE  : A "zibaldone" of useful functions/routines
  !###############################################################  
  include "uniinv.f90"
  include "unista.f90"
  module TOOLS
    USE COMMON_VARS
    USE CHRONOBAR
    implicit none
    private

    !LOOP:
    public :: start_loop
    public :: end_loop

    !GRIDS:
    public :: linspace
    public :: logspace
    public :: arange

    !SORT,UNIQE,SHUFFLE:
    public :: sort, sort_array
    public :: uniq,reshuffle
    public :: shiftFW,shiftBW

    !FUNCTIONS:
    public :: heaviside
    public :: step
    public :: fermi
    public :: sgn
    public :: dens_hyperc

    !BETHE:
    public :: gfbethe,gfbether,bethe_lattice,dens_bethe
    public :: wfun         !complex error function (Faddeeva function)
    public :: zerf   

    !SORT 2D:
    public :: find2Dmesh

    !OTHER:
    public :: get_matsubara_gf_from_dos
    public :: check_convergence,check_convergence_local
    public :: get_free_dos
    public :: sum_overk_zeta

    public :: mix

    interface shiftFW
       module procedure shiftM_fw_C,shiftA_fw_C,shiftM_fw_R,shiftA_fw_R
    end interface shiftFW

    interface shiftBW
       module procedure shiftM_bw_C,shiftA_bw_C,shiftM_bw_R,shiftA_bw_R
    end interface shiftBW

    interface uniq
       module procedure i_uniq,d_uniq
    end interface uniq

    interface sgn
       module procedure i_sgn,d_sgn
    end interface sgn

    interface check_convergence
       module procedure dv_check_convergence,zv_check_convergence,&
            dm_check_convergence,zm_check_convergence
    end interface check_convergence

    interface check_convergence_local
       module procedure dv_check_convergence_local,zv_check_convergence_local,&
            dm_check_convergence_local,zm_check_convergence_local
    end interface check_convergence_local


  contains

    include "mix.f90"


    subroutine start_loop(loop,max,name,unit,id)
      integer                   :: loop
      integer,optional          :: max,unit,id
      character(len=*),optional :: name
      character(len=16)         :: loop_name
      integer                   :: unit_,id_
      loop_name="main-loop";if(present(name))loop_name=name
      unit_    =6          ;if(present(unit))unit_=unit
      id_      =0          ;if(present(id))id_=id
      if(mpiID==id_)then
         write(unit_,*)
         if(.not.present(max))then
            write(unit_,"(A,I5)")bold("-----"//trim(adjustl(trim(loop_name)))),loop,bold("-----")
         else
            write(unit_,"(A,I5,A,I5,A)")bold("-----"//trim(adjustl(trim(loop_name)))),loop,&
                 bold(" (max:"),max,bold(")-----")
         endif
         call start_timer
      endif
    end subroutine start_loop


    !******************************************************************
    !******************************************************************
    !******************************************************************

    subroutine end_loop(unit,id)
      integer,optional :: unit,id
      integer          :: unit_,id_
      unit_=6 ; if(present(unit))unit_=unit
      id_  =0 ; if(present(id))id_=id
      if(mpiID==id_)then
         write(unit_,"(A)")bold("=====================================")
         call stop_timer
         write(unit_,*)
         write(unit_,*)
      endif
    end subroutine end_loop

    !******************************************************************
    !******************************************************************
    !******************************************************************


    !###################################################################
    ! GRIDS:
    !###################################################################
    include "grids.f90"


    !###################################################################
    ! FUNCTIONS:
    !###################################################################
    include "functions.f90"


    !###################################################################
    ! BETHE:
    !###################################################################
    include "bethe.f90"


    !###################################################################
    ! CONVERGENCE:
    !###################################################################
    include "convergence.f90"
    include "convergence_.f90"

    !###################################################################
    ! SORTING 1D:
    !###################################################################
    include "sort1d.f90"
    include "shifts.f90"


    !###################################################################
    ! USEFUL ROUTINES:
    !###################################################################
    !+-------------------------------------------------------------------+
    !PURPOSE  : Evaluate sum over k-points (very simple method)
    !+-------------------------------------------------------------------+
    function sum_overk_zeta(zeta,ek,wk) result(fg)
      complex(8)           :: zeta,fg
      real(8),dimension(:) :: ek,wk
      fg=sum(wk(:)/(zeta-ek(:)))
    end function sum_overk_zeta


    !+-----------------------------------------------------------------+
    !PURPOSE  : 
    !+-----------------------------------------------------------------+
    subroutine get_free_dos(ek,wk,dos,file,store,wmin,wmax,eps)
      real(8),dimension(:)          :: ek,wk
      real(8),dimension(:),optional :: dos
      character(len=*),optional     :: file
      character(len=32)             :: file_     
      logical,optional              :: store
      logical                       :: store_
      real(8),optional              :: wmin,wmax,eps
      real(8)                       :: wini,wfin,eta
      integer,parameter             :: M=2048
      integer                       :: i,ik,Lk,L
      real(8)                       :: w,dew
      complex(8)                    :: gf,iw
      file_="DOSfree.lattice" ; if(present(file))file_=file
      L=M           ; if(present(dos))L=size(dos)
      store_=.true. ; if(present(store))store_=store
      wini = -10.d0 ; if(present(wmin))wini=wmin
      wfin =  10.d0 ; if(present(wmax))wfin=wmax
      eta  = 0.05d0 ; if(present(eps))eta=eps
      Lk =size(ek)  ; dew=abs(wfin-wini)/real(L,8)
      open(70,file=trim(adjustl(trim(file_))))
      do i=1,L
         w  = wini + dble(i-1)*dew;iw=cmplx(w,eta)
         gf = sum(wk/(iw-ek))
         if(present(dos))dos(i)=-aimag(gf)/pi
         write(70,*)w,-aimag(gf)/pi
      enddo
      close(70)
      if(store_)then
         call system("if [ ! -d LATTICEinfo ]; then mkdir LATTICEinfo; fi")
         call system("mv "//trim(adjustl(trim(file_)))//" LATTICEinfo/ 2>/dev/null")
      endif
    end subroutine get_free_dos
    ! subroutine get_DOS(epsk_,file_,wt_,wini,wfin,eps)
    !   implicit none
    !   character(len=*)              :: file_
    !   integer,parameter             :: M=2048
    !   integer                       :: i,ik,Lk_
    !   real(8),optional              :: wini,wfin,eps
    !   real(8)                       :: w,peso,wini_,wfin_,dew_,eps_
    !   complex(8)                    :: gf,iw
    !   real(8),dimension(:)          :: epsk_
    !   real(8),dimension(:),optional :: wt_
    !   wini_ =-20.d0;if(present(wini))wini_=wini
    !   wfin_ =20.d0;if(present(wfin))wfin_=wfin
    !   eps_  =0.1d0;if(present(eps))eps_=eps
    !   Lk_   =size(epsk_)
    !   dew_=abs(wfin_-wini_)/dble(2*M)
    !   open(70,file=trim(file_))
    !   do i=1,2*M
    !      w=wini_ + dble(i-1)*dew_;iw=cmplx(w,eps_)
    !      gf=(0.d0,0.d0)
    !      do ik=1,Lk_
    !         peso=1.d0/Lk_
    !         if(present(wt_))peso=wt_(ik)
    !         gf=gf+peso/(iw-epsk_(ik))
    !      enddo
    !      write(70,*)w,-aimag(gf)/acos(-1.d0)
    !   enddo
    !   close(70)
    ! end subroutine get_DOS




    !+-----------------------------------------------------------------+
    !PURPOSE  : calculate the Matsubara GF given the spectral density  
    !G^M(\iw)= int_\RR d\e A(\e)/iw+\e=-\iw\int_\RR d\e A(\e)/w^2+\e^2 
    !+-----------------------------------------------------------------+
    subroutine get_matsubara_gf_from_dos(win,gin,gout,beta)
      implicit none
      integer :: i,j
      real(8) :: w,wm,beta,A,wini,dw
      real(8) :: gmats_re,gmats_im
      real(8),dimension(:)    :: win
      complex(8),dimension(:) :: gin
      complex(8),dimension(:) :: gout
      wini=minval(win)
      dw=abs(win(2)-win(1)) !Assume constant step in x-grid
      do j=1,size(gout)
         wm=pi/beta*dble(2*j-1)
         gmats_re=(0.d0,0.d0);gmats_im=(0.d0,0.d0)
         do i=1,size(gin)          
            w=wini+dble(i)*dw
            A=aimag(gin(i))/pi !DOS
            gmats_im=gmats_im + dw*A*wm/(wm**2+w**2)
            gmats_re=gmats_re + dw*A*w/(wm**2+w**2)
         enddo
         gout(j)=cmplx(gmats_re,gmats_im)
      enddo
    end subroutine get_matsubara_gf_from_dos


    !*******************************************************************
    !*******************************************************************
    !*******************************************************************




    !###################################################################
    ! SORTING 2D:
    !###################################################################
    !+------------------------------------------------------------------+
    !PURPOSE  : 
    !+------------------------------------------------------------------+
    subroutine find2Dmesh(gridX,gridY,xin,iout)
      implicit none
      integer :: i,j,k
      real(8) :: dist
      integer,dimension(2) :: iout
      real(8),dimension(2) :: xin
      real(8),dimension(:) :: gridX,gridY
      !Find the index of the nearest point along X-axis
      iout(1)=fastsearchreal(xin(1),gridX(:))
      !Find the index of the nearest point along Y-axis
      iout(2)=fastsearchreal(xin(2),gridY(:))

      !Further checks on X:
      k=min(iout(1)+1,size(gridX))
      dist=abs(xin(1)-gridX(iout(1)) )
      if(abs(xin(1)-gridX(k)) < dist ) iout(1)=k
      k=max(iout(1)-1,1)
      if(abs(xin(1)-gridX(k)) < dist ) iout(1)=k

      !Further checks on Y:
      k=min(iout(2)+1,size(gridY))
      dist=abs(xin(2)-gridY(iout(2)))
      if(abs(xin(2)-gridY(k)) < dist ) iout(2)=k
      k=max(iout(2)-1,1)
      if(abs(xin(2)-gridY(1)) < dist ) iout(2)=k
      return
    end subroutine find2Dmesh
    !---------------------------------------------!
    function fastsearchreal(xx,tab)
      implicit none
      !bissection non recursive
      !LES BIN COMMENCE A GAUCHE : gauche_bin_i=tab(xi)
      integer :: i1,i2,is,siz
      real(8)  :: xx,tab(:)
      integer :: fastsearchreal
      siz=size(tab)
      is=siz/2
      i1=1
      i2=siz
      fastsearchreal=0
      if(tab(1)>xx)then
         fastsearchreal=siz
         !write(*,*) 'element pas dans le tableau'
         return
      endif
      if(tab(siz)<=xx)then
         fastsearchreal=siz
         return
      endif
      do
         !*********************************
         if(tab(is)<=xx) then
            i1=is
            is=i1+max(1,(i2-i1)/2)
            goto 28
         endif
         !*********************************
         if(tab(is)>xx)then
            i2=is
            is=i1+(i2-i1)/2
            goto 28
         endif
         !**********************************
28       continue
         if(is==siz.and.tab(is)<=xx)then
            fastsearchreal=is
            return
         endif
         if(is+1<=siz)then
            if(tab(is)<=xx.and.tab(is+1)>xx)then
               fastsearchreal=is
               return
            endif
         endif
      enddo
    end function fastsearchreal



    !*******************************************************************
    !*******************************************************************
    !*******************************************************************






  end module TOOLS
