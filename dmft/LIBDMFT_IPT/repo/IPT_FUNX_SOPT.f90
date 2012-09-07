!###############################################################
!     PROGRAM  : IPT_FUNX_SOPT
!     TYPE     : Module
!     PURPOSE  : Contains common routines for IPT calculations
! using Second Order PErturbation theory in DMFT approximation
!     AUTHORS  : Adriano Amaricci
!###############################################################
module IPT_FUNX_SOPT
  USE IPT_VARS_GLOBAL
  implicit none
  private
  real(8),dimension(:),allocatable    :: A0m,A0p,P1,P2
  integer,allocatable,dimension(:,:)  :: iy_m_ix
  integer,save                        :: loop=1
  complex(8),dimension(:),allocatable :: fg,fg0,sigma
  logical                             :: iprint=.true.

  public :: solve_ipt_sopt

contains
  !+-------------------------------------------------------------------+
  !PROGRAM  : SOLVE_IPT_SOPT
  !TYPE     : function
  !PURPOSE  : interface function for the impurity solver (SOPT) 
  !+-------------------------------------------------------------------+
  function solve_ipt_sopt(fg0_,label_,iprint_) result(sigma_)
    complex(8),dimension(-L:L) :: fg0_,sigma_
    character(len=*),optional  :: label_
    logical,optional           :: iprint_
    if(present(label_))label=trim(adjustl(trim(label_)))
    if(present(iprint_))iprint=iprint_
    if(loop==1)then
       if(.not.allocated(fg))allocate(fg(-L:L))
       if(.not.allocated(fg0))allocate(fg0(-L:L))
       if(.not.allocated(sigma))allocate(sigma(-L:L))
       call get_frequency_index       
    endif
    fg0=fg0_
    call simpurity
    sigma_=sigma
    if(iprint)call print_out()
    loop=loop+1
    return
  end function solve_ipt_sopt
  !******************************************************************
  !******************************************************************
  !******************************************************************



  !+-------------------------------------------------------------------+
  !PROGRAM  : SIMPURITY
  !TYPE     : subroutine
  !PURPOSE  : Evaluate the 2^nd-order perturbation theory self-energy
  !+-------------------------------------------------------------------+
  subroutine Simpurity
    call getAs
    call getPolarization
    call Sopt
  end subroutine Simpurity
  !******************************************************************
  !******************************************************************
  !******************************************************************


  !+-------------------------------------------------------------------+
  !PROGRAM  : CREATE_FREQUENCY_INDEX
  !TYPE     : subroutine
  !PURPOSE  : Create an array of the indices y-x for a faster evaluation
  !+-------------------------------------------------------------------+
  subroutine get_frequency_index()
    integer :: ix,iy,iz
    if(.not.allocated(iy_m_ix))allocate(iy_m_ix(-L:L,-L:L))
    iy_m_ix=0
    do ix=-L,L
       do iy=-L,L
          iz = iy - ix 
          if(iz<-L .OR. iz> L) iz=-L-10 !out of range-> if(iz>-L)
          iy_m_ix(iy,ix)=iz
       enddo
    enddo
    if(.not.allocated(A0m))allocate(A0m(-L:L))
    if(.not.allocated(A0p))allocate(A0p(-L:L))
    if(.not.allocated(P1)) allocate(P1(-L:L))
    if(.not.allocated(P2)) allocate(P2(-L:L))
  end subroutine get_frequency_index
  !******************************************************************
  !******************************************************************
  !******************************************************************





  !+-------------------------------------------------------------------+
  !PROGRAM  : SIMPURITY
  !TYPE     : function
  !PURPOSE  : Evaluate the 2^nd-order perturbation theory self-energy
  !+-------------------------------------------------------------------+
  subroutine getAs
    real(8) :: dos(-L:L)
    dos(:) =-aimag(fg0(:))/pi
    A0p(:) = dos(:)*fermi(wr,beta)
    A0m(:) = dos(:)*(1.d0-fermi(wr,beta))
  end subroutine getAs
  !******************************************************************
  !******************************************************************
  !******************************************************************



  !+-------------------------------------------------------------------+
  !PROGRAM  : SIMPURITY
  !TYPE     : function
  !PURPOSE  : Evaluate the 2^nd-order perturbation theory self-energy
  !+-------------------------------------------------------------------+
  subroutine getPolarization
    integer :: ix,iy,iz
    real(8) :: sum1,sum2
    do ix=-L,L
       sum1=zero
       sum2=zero
       do iy=-L,L
          iz= iy_m_ix(iy,ix)
          if(iz>-L)then
             sum1=sum1 + A0m(iy)*A0p(iz)
             sum2=sum2 + A0p(iy)*A0m(iz)
          endif
       enddo
       P1(ix)=sum1*fmesh
       P2(ix)=sum2*fmesh
    enddo
  end subroutine getPolarization
  !******************************************************************
  !******************************************************************
  !******************************************************************




  !+-------------------------------------------------------------------+
  !PROGRAM  : SIMPURITY
  !TYPE     : function
  !PURPOSE  : Evaluate the 2^nd-order perturbation theory self-energy
  !+-------------------------------------------------------------------+
  subroutine Sopt
    integer :: ix,iy,iz
    real(8) :: sum1,sum2
    real(8),dimension(-L:L) :: reS,imS
    do ix=-L,L
       sum1=zero
       sum2=zero
       do iy=-L,L
          iz= iy_m_ix(iy,ix)
          if(iz>-L)then
             sum1=sum1+A0p(iy)*P1(iz)*fmesh
             sum2=sum2+A0m(iy)*P2(iz)*fmesh
          end if
       enddo
       imS(ix)=-(U**2)*(sum1+sum2)*pi
    enddo
    reS = kronig(imS,wr,size(ImS))
    sigma = reS + xi*imS
  end subroutine Sopt
  !******************************************************************
  !******************************************************************
  !******************************************************************




  !+-------------------------------------------------------------------+
  !PROGRAM  : PRINT_OUT
  !TYPE     : function
  !PURPOSE  : Print out results
  !+-------------------------------------------------------------------+
  subroutine print_out()
    real(8) :: nimp
    fg = one/(one/fg0 - sigma)
    nimp = -2.d0*sum(aimag(fg(-L:L))*fermi(wr(-L:L),beta))*fmesh/pi
    call splot("nVSiloop.ipt"//trim(adjustl(trim(label))),loop,nimp,append=TT)
    call splot("DOS.ipt"//trim(adjustl(trim(label))),wr,-aimag(fg)/pi,append=printf)
    call splot("Sigma_realw.ipt"//trim(adjustl(trim(label))),wr,sigma,append=printf)
    call splot("G_realw.ipt"//trim(adjustl(trim(label))),wr,fg,append=printf)
    call splot("G0_realw.ipt"//trim(adjustl(trim(label))),wr,fg0,append=printf)
    return
  end subroutine print_out
  !******************************************************************
  !******************************************************************
  !******************************************************************




  ! !+-------------------------------------------------------------------+
  ! !PROGRAM  : 
  ! !TYPE     : Subroutine
  ! !PURPOSE  : 
  ! !COMMENT  : 
  ! !+-------------------------------------------------------------------+
  ! subroutine get_energy_sopt(logic,totE)
  !   real(8),optional   :: totE
  !   logical            :: logic
  !   real(8)            :: nimp,docc,Ekin,Epot,Stot,nk(Lk),nkk(Nx,Nx),w
  !   integer            :: Lm,ix,iy,ik
  !   complex(8),allocatable :: gf(:),sf(:),gff(:)
  !   real(8),allocatable    :: gtau(:)

  !   Lm=int(L*beta/pi);if(beta<1)Lm=L;if(Lm<L)Lm=L
  !   allocate(gf(Lm),sf(Lm),gtau(0:L),gff(L))

  !   call getGmats(wr,sigma,sf,beta)
  !   call getGmats(wr,fg,gf,beta)

  !   call fftgf_iw2tau(gf,gtau,beta)
  !   nimp=-2.d0*gtau(L)

  !   !Energy && k-dispersed quantities
  !   Ekin=0.d0
  !   Stot=0.d0
  !   do ik=1,Lk
  !      gff=zero
  !      do i=1,L
  !         w=pi/beta*dble(2*i-1)
  !         gff(i)=one/(xi*w - epsik(ik) - sf(i))
  !      enddo
  !      call fftgf_iw2tau(gff,gtau,beta)
  !      nk(ik)=-gtau(L) 
  !      Ekin=Ekin + wt(ik)*nk(ik)*epsik(ik)
  !      Stot=Stot - wt(ik)*nk(ik)*log(nk(ik))
  !   enddo

  !   Epot=dot_product(conjg(sf(:)),gf(:))/beta/2.0
  !   docc=0.5*nimp  - 0.25d0
  !   if(u/=0.0)docc = Epot/U + 0.5*nimp - 0.25d0
  !   if(present(totE))totE=Ekin+Epot
  !   if(logic)then
  !      call splot("nkVSepsik.ipt"//trim(adjustl(trim(label))),epsik(1:Lk),nk(1:Lk))
  !      call splot("EtotVS"//trim(extension),xout,Ekin+Epot,append=TT)
  !      call splot("StotVS"//trim(extension),xout,Stot,append=TT)
  !      call splot("EkinVS"//trim(extension),xout,Ekin,append=TT)
  !      call splot("EpotVS"//trim(extension),xout,Epot,append=TT)
  !      call splot("doccVS"//trim(extension),xout,docc,append=TT)
  !   end if
  !   deallocate(gf,sf,gtau,gff)
  ! end subroutine get_energy_sopt
  ! !*******************************************************************
  ! !*******************************************************************
  ! !*******************************************************************

end module IPT_FUNX_SOPT
