module DMFT_MISC
  USE SF_TIMER
  USE SF_CONSTANTS, only: pi,zero
  implicit none
  private

  !LOOP:
  public :: start_loop
  public :: end_loop

  !SORT 2D:
  public :: find2Dmesh

  !OTHER:
  ! public :: get_local_density
  public :: order_of_magnitude
  public :: get_density_from_matsubara_gf
  public :: get_matsubara_gf_from_dos
  public :: get_free_dos
  public :: sum_overk_zeta
  public :: finalize_run


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
    write(unit_,*)
    if(.not.present(max))then
       write(unit_,"(A,I5)")"-----"//trim(adjustl(trim(loop_name))),loop,"-----"
    else
       write(unit_,"(A,I5,A,I5,A)")"-----"//trim(adjustl(trim(loop_name))),loop,&
            " (max:",max,")-----"
    endif
    call start_timer
  end subroutine start_loop


  subroutine end_loop(unit,id)
    integer,optional :: unit,id
    integer          :: unit_,id_
    unit_=6 ; if(present(unit))unit_=unit
    id_  =0 ; if(present(id))id_=id
    write(unit_,"(A)")"====================================="
    call stop_timer
    write(unit_,*)
    write(unit_,*)
  end subroutine end_loop


  subroutine finalize_run(iter,error,filename)
    integer,intent(in)                   :: iter
    real(8),intent(in)                   :: error
    character(len=*),intent(in),optional :: filename
    character(len=100)                   :: filename_
    integer                              :: unit
    filename_="job_done.out";if(present(filename))filename_=filename
    unit=897
    open(unit,file=trim(filename),status="new")
    write(unit,*)iter,error
    close(unit)
  end subroutine finalize_run



  function order_of_magnitude(x) result(norder)
    integer :: norder
    real(8) :: x
    norder = floor(log10(abs(x)))
  end function order_of_magnitude



  ! function get_local_density(giw,beta) result(n)
  !   complex(8),dimension(:) :: giw
  !   real(8)                 :: gtau(0:size(giw))
  !   real(8)                 :: beta,n
  !   call fftgf_iw2tau(giw,gtau,beta)
  !   n = -2.d0*gtau(size(giw))
  ! end function get_local_density



  function get_density_from_matsubara_gf(giw,beta) result(n)
    complex(8),dimension(:) :: giw
    complex(8)              :: tail
    real(8)                 :: sum
    real(8)                 :: n,wmax,beta,mues,At,w
    integer                 :: i,Liw
    Liw=size(giw)
    wmax = pi/beta*real(2*Liw-1,8)
    mues =-dreal(giw(Liw))*wmax**2
    sum=0.d0
    do i=1,Liw
       w=pi/beta*real(2*i-1,8)
       tail=-cmplx(mues,w,8)/(mues**2+w**2)
       sum=sum + dreal(giw(i)-tail)
    enddo
    At = -1.d0/(1.d0 + exp(-beta*mues))
    if((mues*beta) >  30.d0)At = -1.d0
    if((mues*beta) < -30.d0)At = -exp(mues*beta)
    n=sum*2.d0/beta+At+1.d0
  end function get_density_from_matsubara_gf


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
    integer :: i1,i2,is,siz
    real(8) :: xx,tab(:)
    integer :: fastsearchreal
    siz=size(tab)
    is=siz/2
    i1=1
    i2=siz
    fastsearchreal=0
    if(tab(1)>xx)then
       fastsearchreal=siz
       return
    endif
    if(tab(siz)<=xx)then
       fastsearchreal=siz
       return
    endif
    do
       if(tab(is)<=xx) then
          i1=is
          is=i1+max(1,(i2-i1)/2)
          goto 28
       endif
       if(tab(is)>xx)then
          i2=is
          is=i1+(i2-i1)/2
          goto 28
       endif
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





end module DMFT_MISC
