!###############################################################
! PROGRAM  : TOOLS
! TYPE     : Module
! PURPOSE  : A "zibaldone" of useful routines
!###############################################################  
module TOOLS
  USE COMMON_VARS
  USE TIMER
  implicit none
  private

  !LOOP:
  public :: start_loop
  public :: end_loop

  !GRIDS:
  public :: linspace
  public :: logspace
  public :: arange
  public :: powspace
  public :: upmspace
  public :: upminterval

  !SORT,UNIQE,SHUFFLE:
  public :: sort,sort_array
  public :: uniq,reshuffle
  public :: shiftFW,shiftBW

  !DERIVATIVE:
  public :: deriv

  !BETHE:
  public :: gfbethe
  public :: gfbether
  public :: bethe_lattice
  public :: dens_bethe
  public :: dens_hyperc

  !SORT 2D:
  public :: find2Dmesh

  !OTHER:
  public :: get_matsubara_gf_from_dos
  public :: check_convergence,check_convergence_scalar,check_convergence_local
  public :: get_free_dos
  public :: sum_overk_zeta

  interface shiftFW
     module procedure shiftM_fw_C,shiftA_fw_C,shiftM_fw_R,shiftA_fw_R
  end interface shiftFW

  interface shiftBW
     module procedure shiftM_bw_C,shiftA_bw_C,shiftM_bw_R,shiftA_bw_R
  end interface shiftBW

  interface uniinv
     module procedure d_uniinv, r_uniinv, i_uniinv
  end interface uniinv
  public :: uniinv

  interface nearless
     module procedure D_nearless, R_nearless, I_nearless
  end interface nearless

  interface unista
     module procedure d_unista, r_unista, i_unista
  end interface unista
  public  :: unista

  interface uniq
     module procedure i_uniq,d_uniq
  end interface uniq

  interface check_convergence_scalar
     module procedure &
          i0_check_convergence_scalar,&
          i1_check_convergence_scalar,&
          i2_check_convergence_scalar,&
          d0_check_convergence_scalar,&
          d1_check_convergence_scalar,&
          d2_check_convergence_scalar,&
          z0_check_convergence_scalar,&
          z1_check_convergence_scalar,&
          z2_check_convergence_scalar
  end interface check_convergence_scalar

  interface check_convergence
     module procedure &
          i0_check_convergence_function,&
          i1_check_convergence_function,&
          i2_check_convergence_function,&
          d0_check_convergence_function,&
          d1_check_convergence_function,&
          d2_check_convergence_function,&
          z0_check_convergence_function,&
          z1_check_convergence_function,&
          z2_check_convergence_function
  end interface check_convergence

  interface check_convergence_local
     module procedure &
          i0_check_convergence_function_local,&
          i1_check_convergence_function_local,&
          i2_check_convergence_function_local,&
          d0_check_convergence_function_local,&
          d1_check_convergence_function_local,&
          d2_check_convergence_function_local,&
          z0_check_convergence_function_local,&
          z1_check_convergence_function_local,&
          z2_check_convergence_function_local
  end interface check_convergence_local





contains


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
  include "tools_grids.f90"


  !###################################################################
  ! FUNCTIONS:
  !###################################################################
  function deriv(f,dh) result(df)
    real(8),dimension(:),intent(in) :: f
    real(8),intent(in)              :: dh
    real(8),dimension(size(f))      :: df
    integer                         :: i,L
    L=size(f)
    df(1)= (f(2)-f(1))/dh
    do i=2,L-1
       df(i) = (f(i+1)-f(i-1))/(2.d0*dh)
    enddo
    df(L)= (f(L)-f(L-1))/dh
  end function deriv



  !+-------------------------------------------------------------------+
  !PURPOSE  : calculate the non-interacting dos for HYPERCUBIC lattice 
  !+-------------------------------------------------------------------+
  pure function dens_hyperc(x,t1)
    real(8),optional,intent(in) :: t1
    real(8),intent(in)          :: x
    REAL(8):: dens_hyperc,t1_,pi2,sqrt2
    pi2=2.d0*acos(-1.d0)
    sqrt2=sqrt(2.d0)
    t1_=sqrt2 ; if(present(t1))t1_=t1
    dens_hyperc = (1/(t1_*sqrt(pi2)))*exp(-(x**2)/(2.d0*t1_**2))
    return
  end function dens_hyperc


  !###################################################################
  ! BETHE:
  !###################################################################
  include "tools_bethe.f90"


  !###################################################################
  ! CONVERGENCE:
  !###################################################################
  include "tools_convergence/convergence_check_scalar.f90"
  include "tools_convergence/convergence_check_function1d.f90"
  include "tools_convergence/convergence_check_function1d_local.f90"

  !###################################################################
  ! SORTING 1D:
  !###################################################################
  include "tools_sort1d.f90"
  include "tools_shifts.f90"


  !###################################################################
  ! USEFUL ROUTINES:
  !###################################################################
  !+-------------------------------------------------------------------+
  !PURPOSE  : Evaluate sum over k-points (very simple method)
  !+-------------------------------------------------------------------+
  function sum_overk_zeta(zeta,ek,wk) result(fg)
    complex(8)                    :: zeta,fg
    real(8),dimension(:)          :: ek
    real(8),dimension(:),optional :: wk
    real(8),dimension(size(ek))   :: wk_
    wk_ = 1.d0/size(ek);if(present(wk))wk_=wk
    fg=sum(wk(:)/(zeta-ek(:)))
  end function sum_overk_zeta


  !+-----------------------------------------------------------------+
  !PURPOSE  : 
  !+-----------------------------------------------------------------+
  subroutine get_free_dos(ek,wk,dos,file,store,wmin,wmax,eps)
    real(8),dimension(:)          :: ek,wk
    real(8),dimension(:),optional :: dos
    character(len=*),optional     :: file
    logical,optional              :: store
    real(8),optional              :: wmin,wmax,eps
    character(len=32)             :: file_     
    logical                       :: store_
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
    eta  = 0.01d0 ; if(present(eps))eta=eps
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




  !+-----------------------------------------------------------------+
  !PURPOSE  : calculate the Matsubara GF given the spectral density  
  !G^M(\iw)= int_\RR d\e A(\e)/iw+\e=-\iw\int_\RR d\e A(\e)/w^2+\e^2 
  !+-----------------------------------------------------------------+
  subroutine get_matsubara_gf_from_dos(win,gin,gout,beta)
    implicit none
    integer :: i,j,Nin,Nout
    real(8) :: w,wm,beta,A,wini,dw
    real(8) :: gmats_re,gmats_im
    real(8),dimension(:)    :: win
    complex(8),dimension(:) :: gin
    complex(8),dimension(:) :: gout
    complex(8),dimension(:),allocatable :: dummy_out
    wini=minval(win)
    dw=abs(win(2)-win(1)) !Assume constant step in x-grid
    Nin=size(gin)
    Nout=size(gout)
    allocate(dummy_out(4*Nout))

    do j=1,Nout
       wm=pi/beta*real(2*j-1,8)
       gmats_re=zero;gmats_im=zero
       do i=1,Nin
          w=wini+dble(i)*dw
          A=aimag(gin(i))/pi !DOS
          gmats_im=gmats_im + dw*A*wm/(wm**2+w**2)
          gmats_re=gmats_re + dw*A*w/(wm**2+w**2)
       enddo
       dummy_out(j)=cmplx(gmats_re,gmats_im)
    enddo
    gout = dummy_out(1:Nout)
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
28     continue
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







  subroutine d_uniinv (xdont, igoest)
    real (kind=8), dimension (:), intent (in) :: xdont
    integer, dimension (:), intent (out)      :: igoest
    real (kind=8) :: xtst, xdona, xdonb
    integer, dimension (size(igoest)) :: jwrkt, irngt
    integer :: lmtna, lmtnc, irng, irng1, irng2, nuni
    integer :: nval, iind, iwrkd, iwrk, iwrkf, jinda, iinda, iindb
    !
    nval = min (size(xdont), size(igoest))
    !
    select case (nval)
    case (:0)
       return
    case (1)
       igoest (1) = 1
       return
    case default
       continue
    end select
    !
    !  fill-in the index array, creating ordered couples
    !
    do iind = 2, nval, 2
       if (xdont(iind-1) < xdont(iind)) then
          irngt (iind-1) = iind - 1
          irngt (iind) = iind
       else
          irngt (iind-1) = iind
          irngt (iind) = iind - 1
       end if
    end do
    if (modulo (nval, 2) /= 0) then
       irngt (nval) = nval
    end if
    !
    !  we will now have ordered subsets a - b - a - b - ...
    !  and merge a and b couples into     c   -   c   - ...
    !
    lmtna = 2
    lmtnc = 4
    !
    !  first iteration. the length of the ordered subsets goes from 2 to 4
    !
    do
       if (nval <= 4) exit
       !
       !   loop on merges of a and b into c
       !
       do iwrkd = 0, nval - 1, 4
          if ((iwrkd+4) > nval) then
             if ((iwrkd+2) >= nval) exit
             !
             !   1 2 3
             !
             if (xdont(irngt(iwrkd+2)) <= xdont(irngt(iwrkd+3))) exit
             !
             !   1 3 2
             !
             if (xdont(irngt(iwrkd+1)) <= xdont(irngt(iwrkd+3))) then
                irng2 = irngt (iwrkd+2)
                irngt (iwrkd+2) = irngt (iwrkd+3)
                irngt (iwrkd+3) = irng2
                !
                !   3 1 2
                !
             else
                irng1 = irngt (iwrkd+1)
                irngt (iwrkd+1) = irngt (iwrkd+3)
                irngt (iwrkd+3) = irngt (iwrkd+2)
                irngt (iwrkd+2) = irng1
             end if
             exit
          end if
          !
          !   1 2 3 4
          !
          if (xdont(irngt(iwrkd+2)) <= xdont(irngt(iwrkd+3))) cycle
          !
          !   1 3 x x
          !
          if (xdont(irngt(iwrkd+1)) <= xdont(irngt(iwrkd+3))) then
             irng2 = irngt (iwrkd+2)
             irngt (iwrkd+2) = irngt (iwrkd+3)
             if (xdont(irng2) <= xdont(irngt(iwrkd+4))) then
                !   1 3 2 4
                irngt (iwrkd+3) = irng2
             else
                !   1 3 4 2
                irngt (iwrkd+3) = irngt (iwrkd+4)
                irngt (iwrkd+4) = irng2
             end if
             !
             !   3 x x x
             !
          else
             irng1 = irngt (iwrkd+1)
             irng2 = irngt (iwrkd+2)
             irngt (iwrkd+1) = irngt (iwrkd+3)
             if (xdont(irng1) <= xdont(irngt(iwrkd+4))) then
                irngt (iwrkd+2) = irng1
                if (xdont(irng2) <= xdont(irngt(iwrkd+4))) then
                   !   3 1 2 4
                   irngt (iwrkd+3) = irng2
                else
                   !   3 1 4 2
                   irngt (iwrkd+3) = irngt (iwrkd+4)
                   irngt (iwrkd+4) = irng2
                end if
             else
                !   3 4 1 2
                irngt (iwrkd+2) = irngt (iwrkd+4)
                irngt (iwrkd+3) = irng1
                irngt (iwrkd+4) = irng2
             end if
          end if
       end do
       !
       !  the cs become as and bs
       !
       lmtna = 4
       exit
    end do
    !
    !  iteration loop. each time, the length of the ordered subsets
    !  is doubled.
    !
    do
       if (2*lmtna >= nval) exit
       iwrkf = 0
       lmtnc = 2 * lmtnc
       !
       !   loop on merges of a and b into c
       !
       do
          iwrk = iwrkf
          iwrkd = iwrkf + 1
          jinda = iwrkf + lmtna
          iwrkf = iwrkf + lmtnc
          if (iwrkf >= nval) then
             if (jinda >= nval) exit
             iwrkf = nval
          end if
          iinda = 1
          iindb = jinda + 1
          !
          !  one steps in the c subset, that we create in the final rank array
          !
          !  make a copy of the rank array for the iteration
          !
          jwrkt (1:lmtna) = irngt (iwrkd:jinda)
          xdona = xdont (jwrkt(iinda))
          xdonb = xdont (irngt(iindb))
          !
          do
             iwrk = iwrk + 1
             !
             !  we still have unprocessed values in both a and b
             !
             if (xdona > xdonb) then
                irngt (iwrk) = irngt (iindb)
                iindb = iindb + 1
                if (iindb > iwrkf) then
                   !  only a still with unprocessed values
                   irngt (iwrk+1:iwrkf) = jwrkt (iinda:lmtna)
                   exit
                end if
                xdonb = xdont (irngt(iindb))
             else
                irngt (iwrk) = jwrkt (iinda)
                iinda = iinda + 1
                if (iinda > lmtna) exit! only b still with unprocessed values
                xdona = xdont (jwrkt(iinda))
             end if
             !
          end do
       end do
       !
       !  the cs become as and bs
       !
       lmtna = 2 * lmtna
    end do
    !
    !   last merge of a and b into c, with removal of duplicates.
    !
    iinda = 1
    iindb = lmtna + 1
    nuni = 0
    !
    !  one steps in the c subset, that we create in the final rank array
    !
    jwrkt (1:lmtna) = irngt (1:lmtna)
    if (iindb <= nval) then
       xtst = nearless (min(xdont(jwrkt(1)), xdont(irngt(iindb))))
    else
       xtst = nearless (xdont(jwrkt(1)))
    endif
    do iwrk = 1, nval
       !
       !  we still have unprocessed values in both a and b
       !
       if (iinda <= lmtna) then
          if (iindb <= nval) then
             if (xdont(jwrkt(iinda)) > xdont(irngt(iindb))) then
                irng = irngt (iindb)
                iindb = iindb + 1
             else
                irng = jwrkt (iinda)
                iinda = iinda + 1
             end if
          else
             !
             !  only a still with unprocessed values
             !
             irng = jwrkt (iinda)
             iinda = iinda + 1
          end if
       else
          !
          !  only b still with unprocessed values
          !
          irng = irngt (iwrk)
       end if
       if (xdont(irng) > xtst) then
          xtst = xdont (irng)
          nuni = nuni + 1
       end if
       igoest (irng) = nuni
       !
    end do
    !
    return
    !
  end subroutine d_uniinv



  !*******************************************************************
  !*******************************************************************
  !*******************************************************************




  subroutine r_uniinv (xdont, igoest)
    real, dimension (:), intent (in) :: xdont
    integer, dimension (:), intent (out) :: igoest
    real    :: xtst, xdona, xdonb
    integer, dimension (size(igoest)) :: jwrkt, irngt
    integer :: lmtna, lmtnc, irng, irng1, irng2, nuni
    integer :: nval, iind, iwrkd, iwrk, iwrkf, jinda, iinda, iindb
    !
    nval = min (size(xdont), size(igoest))
    !
    select case (nval)
    case (:0)
       return
    case (1)
       igoest (1) = 1
       return
    case default
       continue
    end select
    !
    !  fill-in the index array, creating ordered couples
    !
    do iind = 2, nval, 2
       if (xdont(iind-1) < xdont(iind)) then
          irngt (iind-1) = iind - 1
          irngt (iind) = iind
       else
          irngt (iind-1) = iind
          irngt (iind) = iind - 1
       end if
    end do
    if (modulo (nval, 2) /= 0) then
       irngt (nval) = nval
    end if
    !
    !  we will now have ordered subsets a - b - a - b - ...
    !  and merge a and b couples into     c   -   c   - ...
    !
    lmtna = 2
    lmtnc = 4
    !
    !  first iteration. the length of the ordered subsets goes from 2 to 4
    !
    do
       if (nval <= 4) exit
       !
       !   loop on merges of a and b into c
       !
       do iwrkd = 0, nval - 1, 4
          if ((iwrkd+4) > nval) then
             if ((iwrkd+2) >= nval) exit
             !
             !   1 2 3
             !
             if (xdont(irngt(iwrkd+2)) <= xdont(irngt(iwrkd+3))) exit
             !
             !   1 3 2
             !
             if (xdont(irngt(iwrkd+1)) <= xdont(irngt(iwrkd+3))) then
                irng2 = irngt (iwrkd+2)
                irngt (iwrkd+2) = irngt (iwrkd+3)
                irngt (iwrkd+3) = irng2
                !
                !   3 1 2
                !
             else
                irng1 = irngt (iwrkd+1)
                irngt (iwrkd+1) = irngt (iwrkd+3)
                irngt (iwrkd+3) = irngt (iwrkd+2)
                irngt (iwrkd+2) = irng1
             end if
             exit
          end if
          !
          !   1 2 3 4
          !
          if (xdont(irngt(iwrkd+2)) <= xdont(irngt(iwrkd+3))) cycle
          !
          !   1 3 x x
          !
          if (xdont(irngt(iwrkd+1)) <= xdont(irngt(iwrkd+3))) then
             irng2 = irngt (iwrkd+2)
             irngt (iwrkd+2) = irngt (iwrkd+3)
             if (xdont(irng2) <= xdont(irngt(iwrkd+4))) then
                !   1 3 2 4
                irngt (iwrkd+3) = irng2
             else
                !   1 3 4 2
                irngt (iwrkd+3) = irngt (iwrkd+4)
                irngt (iwrkd+4) = irng2
             end if
             !
             !   3 x x x
             !
          else
             irng1 = irngt (iwrkd+1)
             irng2 = irngt (iwrkd+2)
             irngt (iwrkd+1) = irngt (iwrkd+3)
             if (xdont(irng1) <= xdont(irngt(iwrkd+4))) then
                irngt (iwrkd+2) = irng1
                if (xdont(irng2) <= xdont(irngt(iwrkd+4))) then
                   !   3 1 2 4
                   irngt (iwrkd+3) = irng2
                else
                   !   3 1 4 2
                   irngt (iwrkd+3) = irngt (iwrkd+4)
                   irngt (iwrkd+4) = irng2
                end if
             else
                !   3 4 1 2
                irngt (iwrkd+2) = irngt (iwrkd+4)
                irngt (iwrkd+3) = irng1
                irngt (iwrkd+4) = irng2
             end if
          end if
       end do
       !
       !  the cs become as and bs
       !
       lmtna = 4
       exit
    end do
    !
    !  iteration loop. each time, the length of the ordered subsets
    !  is doubled.
    !
    do
       if (2*lmtna >= nval) exit
       iwrkf = 0
       lmtnc = 2 * lmtnc
       !
       !   loop on merges of a and b into c
       !
       do
          iwrk = iwrkf
          iwrkd = iwrkf + 1
          jinda = iwrkf + lmtna
          iwrkf = iwrkf + lmtnc
          if (iwrkf >= nval) then
             if (jinda >= nval) exit
             iwrkf = nval
          end if
          iinda = 1
          iindb = jinda + 1
          !
          !  one steps in the c subset, that we create in the final rank array
          !
          !  make a copy of the rank array for the iteration
          !
          jwrkt (1:lmtna) = irngt (iwrkd:jinda)
          xdona = xdont (jwrkt(iinda))
          xdonb = xdont (irngt(iindb))
          !
          do
             iwrk = iwrk + 1
             !
             !  we still have unprocessed values in both a and b
             !
             if (xdona > xdonb) then
                irngt (iwrk) = irngt (iindb)
                iindb = iindb + 1
                if (iindb > iwrkf) then
                   !  only a still with unprocessed values
                   irngt (iwrk+1:iwrkf) = jwrkt (iinda:lmtna)
                   exit
                end if
                xdonb = xdont (irngt(iindb))
             else
                irngt (iwrk) = jwrkt (iinda)
                iinda = iinda + 1
                if (iinda > lmtna) exit! only b still with unprocessed values
                xdona = xdont (jwrkt(iinda))
             end if
             !
          end do
       end do
       !
       !  the cs become as and bs
       !
       lmtna = 2 * lmtna
    end do
    !
    !   last merge of a and b into c, with removal of duplicates.
    !
    iinda = 1
    iindb = lmtna + 1
    nuni = 0
    !
    !  one steps in the c subset, that we create in the final rank array
    !
    jwrkt (1:lmtna) = irngt (1:lmtna)
    if (iindb <= nval) then
       xtst = nearless (min(xdont(jwrkt(1)), xdont(irngt(iindb))))
    else
       xtst = nearless (xdont(jwrkt(1)))
    endif
    do iwrk = 1, nval
       !
       !  we still have unprocessed values in both a and b
       !
       if (iinda <= lmtna) then
          if (iindb <= nval) then
             if (xdont(jwrkt(iinda)) > xdont(irngt(iindb))) then
                irng = irngt (iindb)
                iindb = iindb + 1
             else
                irng = jwrkt (iinda)
                iinda = iinda + 1
             end if
          else
             !
             !  only a still with unprocessed values
             !
             irng = jwrkt (iinda)
             iinda = iinda + 1
          end if
       else
          !
          !  only b still with unprocessed values
          !
          irng = irngt (iwrk)
       end if
       if (xdont(irng) > xtst) then
          xtst = xdont (irng)
          nuni = nuni + 1
       end if
       igoest (irng) = nuni
       !
    end do
    !
    return
    !
  end subroutine r_uniinv
  subroutine i_uniinv (xdont, igoest)
    ! __________________________________________________________
    !   uniinv = merge-sort inverse ranking of an array, with removal of
    !   duplicate entries.
    !   the routine is similar to pure merge-sort ranking, but on
    !   the last pass, it sets indices in igoest to the rank
    !   of the value in the ordered set with duplicates removed.
    !   for performance reasons, the first 2 passes are taken
    !   out of the standard loop, and use dedicated coding.
    ! __________________________________________________________
    ! __________________________________________________________
    integer, dimension (:), intent (in)  :: xdont
    integer, dimension (:), intent (out) :: igoest
    ! __________________________________________________________
    integer :: xtst, xdona, xdonb
    !
    ! __________________________________________________________
    integer, dimension (size(igoest)) :: jwrkt, irngt
    integer :: lmtna, lmtnc, irng, irng1, irng2, nuni
    integer :: nval, iind, iwrkd, iwrk, iwrkf, jinda, iinda, iindb
    !
    nval = min (size(xdont), size(igoest))
    !
    select case (nval)
    case (:0)
       return
    case (1)
       igoest (1) = 1
       return
    case default
       continue
    end select
    !
    !  fill-in the index array, creating ordered couples
    !
    do iind = 2, nval, 2
       if (xdont(iind-1) < xdont(iind)) then
          irngt (iind-1) = iind - 1
          irngt (iind) = iind
       else
          irngt (iind-1) = iind
          irngt (iind) = iind - 1
       end if
    end do
    if (modulo (nval, 2) /= 0) then
       irngt (nval) = nval
    end if
    !
    !  we will now have ordered subsets a - b - a - b - ...
    !  and merge a and b couples into     c   -   c   - ...
    !
    lmtna = 2
    lmtnc = 4
    !
    !  first iteration. the length of the ordered subsets goes from 2 to 4
    !
    do
       if (nval <= 4) exit
       !
       !   loop on merges of a and b into c
       !
       do iwrkd = 0, nval - 1, 4
          if ((iwrkd+4) > nval) then
             if ((iwrkd+2) >= nval) exit
             !
             !   1 2 3
             !
             if (xdont(irngt(iwrkd+2)) <= xdont(irngt(iwrkd+3))) exit
             !
             !   1 3 2
             !
             if (xdont(irngt(iwrkd+1)) <= xdont(irngt(iwrkd+3))) then
                irng2 = irngt (iwrkd+2)
                irngt (iwrkd+2) = irngt (iwrkd+3)
                irngt (iwrkd+3) = irng2
                !
                !   3 1 2
                !
             else
                irng1 = irngt (iwrkd+1)
                irngt (iwrkd+1) = irngt (iwrkd+3)
                irngt (iwrkd+3) = irngt (iwrkd+2)
                irngt (iwrkd+2) = irng1
             end if
             exit
          end if
          !
          !   1 2 3 4
          !
          if (xdont(irngt(iwrkd+2)) <= xdont(irngt(iwrkd+3))) cycle
          !
          !   1 3 x x
          !
          if (xdont(irngt(iwrkd+1)) <= xdont(irngt(iwrkd+3))) then
             irng2 = irngt (iwrkd+2)
             irngt (iwrkd+2) = irngt (iwrkd+3)
             if (xdont(irng2) <= xdont(irngt(iwrkd+4))) then
                !   1 3 2 4
                irngt (iwrkd+3) = irng2
             else
                !   1 3 4 2
                irngt (iwrkd+3) = irngt (iwrkd+4)
                irngt (iwrkd+4) = irng2
             end if
             !
             !   3 x x x
             !
          else
             irng1 = irngt (iwrkd+1)
             irng2 = irngt (iwrkd+2)
             irngt (iwrkd+1) = irngt (iwrkd+3)
             if (xdont(irng1) <= xdont(irngt(iwrkd+4))) then
                irngt (iwrkd+2) = irng1
                if (xdont(irng2) <= xdont(irngt(iwrkd+4))) then
                   !   3 1 2 4
                   irngt (iwrkd+3) = irng2
                else
                   !   3 1 4 2
                   irngt (iwrkd+3) = irngt (iwrkd+4)
                   irngt (iwrkd+4) = irng2
                end if
             else
                !   3 4 1 2
                irngt (iwrkd+2) = irngt (iwrkd+4)
                irngt (iwrkd+3) = irng1
                irngt (iwrkd+4) = irng2
             end if
          end if
       end do
       !
       !  the cs become as and bs
       !
       lmtna = 4
       exit
    end do
    !
    !  iteration loop. each time, the length of the ordered subsets
    !  is doubled.
    !
    do
       if (2*lmtna >= nval) exit
       iwrkf = 0
       lmtnc = 2 * lmtnc
       !
       !   loop on merges of a and b into c
       !
       do
          iwrk = iwrkf
          iwrkd = iwrkf + 1
          jinda = iwrkf + lmtna
          iwrkf = iwrkf + lmtnc
          if (iwrkf >= nval) then
             if (jinda >= nval) exit
             iwrkf = nval
          end if
          iinda = 1
          iindb = jinda + 1
          !
          !  one steps in the c subset, that we create in the final rank array
          !
          !  make a copy of the rank array for the iteration
          !
          jwrkt (1:lmtna) = irngt (iwrkd:jinda)
          xdona = xdont (jwrkt(iinda))
          xdonb = xdont (irngt(iindb))
          !
          do
             iwrk = iwrk + 1
             !
             !  we still have unprocessed values in both a and b
             !
             if (xdona > xdonb) then
                irngt (iwrk) = irngt (iindb)
                iindb = iindb + 1
                if (iindb > iwrkf) then
                   !  only a still with unprocessed values
                   irngt (iwrk+1:iwrkf) = jwrkt (iinda:lmtna)
                   exit
                end if
                xdonb = xdont (irngt(iindb))
             else
                irngt (iwrk) = jwrkt (iinda)
                iinda = iinda + 1
                if (iinda > lmtna) exit! only b still with unprocessed values
                xdona = xdont (jwrkt(iinda))
             end if
             !
          end do
       end do
       !
       !  the cs become as and bs
       !
       lmtna = 2 * lmtna
    end do
    !
    !   last merge of a and b into c, with removal of duplicates.
    !
    iinda = 1
    iindb = lmtna + 1
    nuni = 0
    !
    !  one steps in the c subset, that we create in the final rank array
    !
    jwrkt (1:lmtna) = irngt (1:lmtna)
    if (iindb <= nval) then
       xtst = nearless (min(xdont(jwrkt(1)), xdont(irngt(iindb))))
    else
       xtst = nearless (xdont(jwrkt(1)))
    endif
    do iwrk = 1, nval
       !
       !  we still have unprocessed values in both a and b
       !
       if (iinda <= lmtna) then
          if (iindb <= nval) then
             if (xdont(jwrkt(iinda)) > xdont(irngt(iindb))) then
                irng = irngt (iindb)
                iindb = iindb + 1
             else
                irng = jwrkt (iinda)
                iinda = iinda + 1
             end if
          else
             !
             !  only a still with unprocessed values
             !
             irng = jwrkt (iinda)
             iinda = iinda + 1
          end if
       else
          !
          !  only b still with unprocessed values
          !
          irng = irngt (iwrk)
       end if
       if (xdont(irng) > xtst) then
          xtst = xdont (irng)
          nuni = nuni + 1
       end if
       igoest (irng) = nuni
       !
    end do
    !
    return
    !
  end subroutine i_uniinv




  !*******************************************************************
  !*******************************************************************
  !*******************************************************************




  function d_nearless (xval) result (d_nl)
    !  nearest value less than given value
    real (kind=8), intent (in) :: xval
    real (kind=8) :: d_nl
    d_nl = nearest (xval, -1.0d0)
    return
  end function d_nearless
  function r_nearless (xval) result (r_nl)
    !  nearest value less than given value
    real, intent (in) :: xval
    real :: r_nl
    r_nl = nearest (xval, -1.0)
    return
  end function r_nearless
  function i_nearless (xval) result (i_nl)
    !  nearest value less than given value
    integer, intent (in) :: xval
    integer :: i_nl
    i_nl = xval - 1
    return
  end function i_nearless




  !*******************************************************************
  !*******************************************************************
  !*******************************************************************



  subroutine d_unista (xdont, nuni, mask)
    real(kind=8), dimension (:), intent (inout) :: xdont
    integer, intent (out)                       :: nuni
    integer, dimension (size(xdont)) :: iwrkt
    logical, dimension (size(xdont)) :: ifmptyt
    logical, dimension (size(xdont)),optional :: mask
    integer :: icrs
    call uniinv (xdont, iwrkt)

    ifmptyt = .true.
    nuni = 0
    do icrs = 1, size(xdont)
       if(present(mask))mask(icrs)=ifmptyt(iwrkt(icrs))
       if (ifmptyt(iwrkt(icrs))) then
          ifmptyt(iwrkt(icrs)) = .false. !modify the mask to that next acces is false
          nuni = nuni + 1
          xdont (nuni) = xdont (icrs)
       end if
    end do
    return
  end subroutine d_unista

  subroutine r_unista (xdont, nuni, mask)
    real, dimension (:), intent (inout) :: xdont
    integer, intent (out)               :: nuni
    integer, dimension (size(xdont)) :: iwrkt
    logical, dimension (size(xdont)) :: ifmptyt
    logical, dimension (size(xdont)),optional :: mask
    integer :: icrs
    call uniinv (xdont, iwrkt)
    ifmptyt = .true.
    nuni = 0
    do icrs = 1, size(xdont)
       if(present(mask))mask(icrs)=ifmptyt(iwrkt(icrs))
       if (ifmptyt(iwrkt(icrs))) then
          ifmptyt(iwrkt(icrs)) = .false.
          nuni = nuni + 1
          xdont (nuni) = xdont (icrs)
       end if
    end do
    return
  end subroutine r_unista

  subroutine i_unista (xdont, nuni, mask)
    integer, dimension (:), intent (inout)  :: xdont
    integer, intent (out) :: nuni
    integer, dimension (size(xdont)) :: iwrkt
    logical, dimension (size(xdont)) :: ifmptyt
    logical, dimension (size(xdont)),optional :: mask
    integer :: icrs
    call uniinv (xdont, iwrkt)
    ifmptyt = .true.
    nuni = 0
    do icrs = 1, size(xdont)
       if(present(mask))mask(icrs)=ifmptyt(iwrkt(icrs))
       if (ifmptyt(iwrkt(icrs))) then
          ifmptyt(iwrkt(icrs)) = .false.
          nuni = nuni + 1
          xdont (nuni) = xdont (icrs)
       end if
    end do
    return
  end subroutine i_unista



end module TOOLS
