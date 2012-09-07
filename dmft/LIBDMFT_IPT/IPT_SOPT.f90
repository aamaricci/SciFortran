!###############################################################
!     PROGRAM  : IPT_FUNX_SOPT
!     TYPE     : Module
!     PURPOSE  : Contains common routines for IPT calculations
! using Second Order PErturbation theory in DMFT approximation
!     AUTHORS  : Adriano Amaricci
!###############################################################
module IPT_SOPT
  USE IPT_VARS_GLOBAL
  implicit none
  private
  real(8),dimension(:),allocatable    :: A0m,A0p,P1,P2
  integer,allocatable,dimension(:,:)  :: iy_m_ix
  integer,save                        :: loop=1
  complex(8),dimension(:),allocatable :: fg0,sigma
  real(8),dimension(:),allocatable    :: wr
  real(8)                             :: n,n0,xmu0,mesh
  integer                             :: MM

  public :: solve_ipt_sopt
  public :: solve_mpt_sopt

contains

  !+-------------------------------------------------------------------+
  !PURPOSE  : interface function for the impurity solver (SOPT) 
  !+-------------------------------------------------------------------+
  function solve_ipt_sopt(fg0_,wr_) result(sigma_)
    complex(8),dimension(:)          :: fg0_
    complex(8),dimension(size(fg0_)) :: sigma_
    real(8),dimension(size(fg0_))    :: wr_
    MM=size(fg0_)
    if(loop==1)then
       if(.not.allocated(wr))allocate(wr(MM))
       if(.not.allocated(fg0))allocate(fg0(MM))
       if(.not.allocated(sigma))allocate(sigma(MM))
       call get_frequency_index       
    endif
    fg0=fg0_; wr=wr_ ; mesh=abs(wr(2)-wr(1))
    call getAs
    call getPolarization
    call Sopt
    sigma_=sigma
    loop=loop+1
  end function solve_ipt_sopt



  !******************************************************************
  !******************************************************************
  !******************************************************************



  !+-------------------------------------------------------------------+
  !PURPOSE  : interface function for the impurity solver (SOPT) 
  !+-------------------------------------------------------------------+
  function solve_mpt_sopt(fg0_,wr_,n_,n0_,xmu0_) result(sigma_)
    complex(8),dimension(:)          :: fg0_
    complex(8),dimension(size(fg0_)) :: sigma_
    real(8),dimension(size(fg0_))    :: wr_
    real(8)                          :: A,B,n_,n0_,xmu0_
    MM=size(fg0_)
    if(loop==1) then
       if(.not.allocated(fg0))allocate(fg0(MM))
       if(.not.allocated(sigma))allocate(sigma(MM))
       if(.not.allocated(wr))allocate(wr(MM))
       call get_frequency_index
    endif
    fg0=fg0_; wr=wr_ ; mesh=abs(wr(2)-wr(1))
    n=n_    ; n0=n0_ ; xmu0=xmu0_
    A=0.d0  ; B=0.d0
    call getAs
    call getPolarization
    call Sopt
    call get_A
    call get_B
    sigma_ = U*(n-0.5d0) + A*sigma/(1.d0-B*sigma)
    loop=loop+1
  contains
    subroutine get_A
      real(8) :: A1,A2
      A1=n*(1.d0-n)
      A2=n0*(1.d0-n0)
      A=A1/A2
    end subroutine get_A
    subroutine get_B
      real(8) :: B1,B2
      B1 = (xmu0-xmu) + U*(1.d0-2.d0*n)
      B2 = n0*(1.d0-n0)*U**2
      B=B1/B2
    end subroutine get_B
  end function solve_mpt_sopt



  !******************************************************************
  !******************************************************************
  !******************************************************************





  !+-------------------------------------------------------------------+
  !PURPOSE  : Create an array of the indices y-x for a faster evaluation
  !+-------------------------------------------------------------------+
  subroutine get_frequency_index()
    integer :: ix,iy,iz
    if(.not.allocated(iy_m_ix))allocate(iy_m_ix(MM,MM))
    iy_m_ix=0
    do ix=1,MM
       do iy=1,MM
          iz = iy - ix + MM/2 
          if(iz<1 .OR. iz>MM) iz=-1 !out of range-> if(iz>-L)
          iy_m_ix(iy,ix)=iz
       enddo
    enddo
    if(.not.allocated(A0m))allocate(A0m(MM))
    if(.not.allocated(A0p))allocate(A0p(MM))
    if(.not.allocated(P1)) allocate(P1(MM))
    if(.not.allocated(P2)) allocate(P2(MM))
  end subroutine get_frequency_index




  !******************************************************************
  !******************************************************************
  !******************************************************************



  subroutine getAs
    real(8) :: dos(MM)
    dos(:) =-aimag(fg0(:))/pi
    A0p(:) = dos(:)*fermi(wr(:),beta)
    A0m(:) = dos(:)*(1.d0-fermi(wr(:),beta))
  end subroutine getAs



  !******************************************************************
  !******************************************************************
  !******************************************************************



  subroutine getPolarization
    integer :: ix,iy,iz    
    P1=zero
    P2=zero
    do ix=1,MM
       do iy=1,MM
          iz= iy_m_ix(iy,ix)
          if(iz>0)then
             P1(ix)=P1(ix) + A0p(iy)*A0m(iz)*mesh
             P2(ix)=P2(ix) + A0m(iy)*A0p(iz)*mesh
          endif
       enddo
    enddo
  end subroutine getPolarization




  !******************************************************************
  !******************************************************************
  !******************************************************************




  !+-------------------------------------------------------------------+
  !PURPOSE  : Solve 2^nd order perturbation theory
  !+-------------------------------------------------------------------+
  subroutine Sopt
    integer :: ix,iy,iz
    real(8) :: sum1,sum2
    real(8),dimension(MM) :: reS,imS
    do ix=1,MM
       sum1=zero
       sum2=zero
       do iy=1,MM
          iz= iy_m_ix(iy,ix)
          if(iz>0)then
             sum1=sum1+A0p(MM-iz+1)*P1(iy)*mesh
             sum2=sum2+A0m(MM-iz+1)*P2(iy)*mesh
          end if
       enddo
       imS(ix)=-(U**2)*(sum1+sum2)*pi
    enddo
    reS = kronig(imS,wr,size(ImS))
    sigma = reS + xi*imS
  end subroutine Sopt


end module IPT_SOPT
