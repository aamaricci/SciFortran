!###############################################################
!     PROGRAM  : MPT_FUNX_SOPT
!     TYPE     : Module
!     PURPOSE  : Contains common routines
!     AUTHORS  : Adriano Amaricci
!###############################################################
module MPT_DISORDER_FUNX_SOPT
  USE IPT_VARS_GLOBAL
  implicit none
  private

  real(8),dimension(:),allocatable      :: A0m,A0p,P1,P2
  integer,allocatable,dimension(:,:)    :: iy_m_ix
  complex(8),allocatable,dimension(:,:) :: sold
  integer                               :: loop=1,Isite,Nsite
  complex(8),dimension(:,:),allocatable :: fg0,sigma
  real(8)                               :: n,n0,xmu0

  save sold

  public :: solve_mpt_disorder_sopt

contains

  !+-------------------------------------------------------------------+
  !PROGRAM  : SOLVE_MPT_finiteT
  !TYPE     : function
  !PURPOSE  : 
  !+-------------------------------------------------------------------+
  function solve_mpt_disorder_sopt(is,ns,fg0_,n_,n0_,xmu0_,peso) result(sigma_)
    integer                    :: is,ns
    complex(8),dimension(-L:L) :: fg0_,sigma_
    real(8)                    :: peso
    real(8)                    :: n_,n0_,xmu0_

    Isite=is;Nsite=ns
    if(loop==1) then
       if(.not.allocated(fg0))allocate(fg0(Nsite,-L:L))
       if(.not.allocated(sigma))allocate(sigma(Nsite,-L:L))
       if(.not.allocated(sold))allocate(sold(Nsite,-L:L))
       call get_frequency_index
       sold=sigma
    endif
    fg0(Isite,:)=fg0_ ; n=n_ ; n0=n0_ ; xmu0=xmu0_
    call simpurity
    sigma(Isite,:) = peso*sigma(Isite,:) + (1.d0-peso)*sold(Isite,:)
    sold(Isite,:)=sigma(Isite,:)
    sigma_=sigma(Isite,:)
    loop=loop+1
  end function solve_mpt_disorder_sopt
  !******************************************************************
  !******************************************************************
  !******************************************************************



  !+-------------------------------------------------------------------+
  !PROGRAM  : SIMPURITY
  !TYPE     : function
  !PURPOSE  : Evaluate the 2^nd-order perturbation theory self-energy
  !+-------------------------------------------------------------------+
  subroutine Simpurity
    real(8) :: A,B
    call getAs
    call getPolarization
    call Sopt
    A=(n*(1.d0-n))/(n0*(1.d0-n0))
    B=((xmu0-xmu) + U*(1.d0-2.d0*n))/(n0*(1.d0-n0)*U**2)
    sigma(Isite,:) = U*(n-0.5d0) + A*sigma(Isite,:)/(1.d0-B*sigma(Isite,:))
  end subroutine Simpurity
  !******************************************************************
  !******************************************************************
  !******************************************************************



  !+-------------------------------------------------------------------+
  !PROGRAM  : CREATE_FREQUENCY_INDEX
  !TYPE     : function
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
    integer :: i
    real(8) :: w,dos
    do i=-L,L       
       w=wr(i); dos=-dimag(fg0(isite,i))/pi !G0 here is the Isite^th one
       A0p(i) = dos*fermi(w,beta)
       A0m(i) = dos*(1.d0-fermi(w,beta))
    enddo
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
    reS = kronig(imS,wr,size(wr))
    sigma(Isite,:) = reS + xi*imS
  end subroutine Sopt
  !******************************************************************
  !******************************************************************
  !******************************************************************











end module MPT_DISORDER_FUNX_SOPT


