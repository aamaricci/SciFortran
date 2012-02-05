!###############################################################
! PROGRAM  : FUNX_IPT
! TYPE     : Module
! PURPOSE  : Contains common routines
! AUTHORS  : Adriano Amaricci
! COMMENTS : little bit dirty... but seems to work
!###############################################################
module IPT_SC_SOPT
  USE IPT_VARS_GLOBAL
  implicit none
  private

  real(8),dimension(:),allocatable      :: A0p11,A0m11,A0p22,A0m22,B0p,B0m,C0p,C0m
  real(8),dimension(:),allocatable      :: P1,P2,Q1,Q2,R1,R2,T1,T2
  real(8),dimension(:),allocatable      :: wr
  integer,allocatable,dimension(:,:)    :: iy_m_ix
  integer,save                          :: loop=1
  complex(8),allocatable,dimension(:,:) :: fg0,sigma
  real(8)                               :: n,n0,delta,delta0
  integer                               :: MM
  public                                :: solve_ipt_sc_sopt
  public                                :: solve_mpt_sc_sopt


contains


  !+-------------------------------------------------------------------+
  !PURPOSE  : 
  !+-------------------------------------------------------------------+
  function solve_ipt_sc_sopt(fg0_,wr_,delta_) result(sigma_)
    complex(8),dimension(:,:)                       :: fg0_
    complex(8),dimension(size(fg0_,1),size(fg0_,2)) :: sigma_
    real(8),dimension(size(fg0_,2))                 :: wr_
    real(8)                    :: delta_
    MM=size(fg0_,2)
    if(loop==1)then
       if(.not.allocated(fg0))allocate(fg0(2,MM))
       if(.not.allocated(sigma))allocate(sigma(2,MM))
       if(.not.allocated(wr))allocate(wr(MM))
       call get_frequency_index
    endif
    fg0=fg0_ ; delta=delta_ ; wr=wr_ ; fmesh=abs(wr(2)-wr(1))
    call simpurity 
    sigma_(1,:) = sigma(1,:)
    sigma_(2,:) = -delta + sigma(2,:)
    loop=loop+1
  end function solve_ipt_sc_sopt




  !******************************************************************
  !******************************************************************
  !******************************************************************




  !+-------------------------------------------------------------------+
  !PURPOSE  : 
  !+-------------------------------------------------------------------+
  function solve_mpt_sc_sopt(fg0_,wr_,n_,n0_,delta_,delta0_) result(sigma_)
    complex(8),dimension(:,:)                       :: fg0_
    complex(8),dimension(size(fg0_,1),size(fg0_,2)) :: sigma_
    real(8),dimension(size(fg0_,2))                 :: wr_
    real(8)                                         :: n_,n0_,delta_,delta0_
    real(8)                                         :: A,B
    MM=size(fg0_,2)
    if(loop==1) then
       if(.not.allocated(fg0))allocate(fg0(2,MM))
       if(.not.allocated(sigma))allocate(sigma(2,MM))
       if(.not.allocated(wr))allocate(wr(MM))
       call get_frequency_index  
    endif
    fg0=fg0_  ; wr=wr_ ; n=n_ ; n0=n0_ ; delta=delta_ ; delta0=delta0_
    fmesh=abs(wr(2)-wr(1))
    call simpurity
    A = U**2*n*(1.d0-n)-delta**2
    B = U**2*n0*(1.d0-n0)-delta0**2
    sigma_(1,:) =  u*(n-0.5d0)   + sigma(1,:)*A/B
    sigma_(2,:) = -delta         + sigma(2,:)*A/B
    loop=loop+1
  end function solve_mpt_sc_sopt




  !******************************************************************
  !******************************************************************
  !******************************************************************




  !+-------------------------------------------------------------------+
  !PURPOSE  : Create an array of the indices y-x for a faster evaluation
  !+-------------------------------------------------------------------+
  subroutine get_frequency_index()
    integer :: ix,iy,iz
    allocate(iy_m_ix(MM,MM))
    iy_m_ix=0
    do ix=1,MM
       do iy=1,MM
          iz = iy - ix + MM/2
          if(iz<1 .OR. iz>MM) iz=-MM-10
          iy_m_ix(iy,ix)=iz
       enddo
    enddo
    if(.not.allocated(A0p11))allocate(A0p11(MM))
    if(.not.allocated(A0m11))allocate(A0m11(MM))
    if(.not.allocated(A0p22))allocate(A0p22(MM))
    if(.not.allocated(A0m22))allocate(A0m22(MM))
    if(.not.allocated(B0p))allocate(B0p(MM))
    if(.not.allocated(B0m))allocate(B0m(MM))
    if(.not.allocated(C0p))allocate(C0p(MM))
    if(.not.allocated(C0m))allocate(C0m(MM))
    if(.not.allocated(P1))allocate(P1(MM))
    if(.not.allocated(P2))allocate(P2(MM))
    if(.not.allocated(Q1))allocate(Q1(MM))
    if(.not.allocated(Q2))allocate(Q2(MM))
    if(.not.allocated(R1))allocate(R1(MM))
    if(.not.allocated(R2))allocate(R2(MM))
    if(.not.allocated(T1))allocate(T1(MM))
    if(.not.allocated(T2))allocate(T2(MM))
  end subroutine get_frequency_index
  !******************************************************************
  !******************************************************************
  !******************************************************************





  !+-------------------------------------------------------------------+
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
  !PURPOSE  : 
  !+-------------------------------------------------------------------+
  subroutine getAs
    integer :: i
    real(8) :: w,dos11(MM),dos22(MM),dosF1(MM),dosF2(MM)
    do i=1,MM
       dos11(i) = -aimag(fg0(1,i))/pi
       dos22(i) = -aimag(-conjg(fg0(1,MM-i+1)))/pi
       dosF1(i) = -aimag(fg0(2,i))/pi
       dosF2(i) = -aimag(fg0(2,MM-i+1))/pi
    enddo
    A0p11(:) = dos11*fermi(wr,beta)
    A0m11(:) = dos11*(1.d0-fermi(wr,beta))
    A0p22(:) = dos22*fermi(wr,beta)
    A0m22(:) = dos22*(1.d0-fermi(wr,beta))
    B0p(:) = dosF1*fermi(wr,beta)
    B0m(:) = dosF1*(1.d0-fermi(wr,beta))
    C0p(:) = dosF2*fermi(wr,beta)
    C0m(:) = dosF2*(1.d0-fermi(wr,beta))
  end subroutine getAs




  !******************************************************************
  !******************************************************************
  !******************************************************************




  !+-------------------------------------------------------------------+
  !PURPOSE  :
  !+-------------------------------------------------------------------+
  subroutine getPolarization
    integer :: ix,iy,iz
    P1=zero
    P2=zero
    Q1=zero
    Q2=zero
    R1=zero
    R2=zero
    T1=zero
    T2=zero
    do ix=1,MM
       do iy=1,MM
          iz= iy_m_ix(iy,ix)
          if(iz>-L)then
             P1(ix)=P1(ix) + A0p11(iy)*A0m22(iz)*fmesh
             P2(ix)=P2(ix) + A0m11(iy)*A0p22(iz)*fmesh
             Q1(ix)=Q1(ix) + C0p(iy)*A0m22(iz)*fmesh !C0p(iy)*A0m22(iz)*fmesh
             Q2(ix)=Q2(ix) + C0m(iy)*A0p22(iz)*fmesh !C0m(iy)*A0p22(iz)*fmesh
             R1(ix)=R1(ix) + B0p(iy)*B0m(iz)*fmesh
             R2(ix)=R2(ix) + B0m(iy)*B0p(iz)*fmesh
             T1(ix)=T1(ix) + A0p11(iy)*B0m(iz)*fmesh
             T2(ix)=T2(ix) + A0m11(iy)*B0p(iz)*fmesh
          endif
       enddo
    enddo
  end subroutine getPolarization




  !******************************************************************
  !******************************************************************
  !******************************************************************




  !+-------------------------------------------------------------------+
  !PURPOSE  : 
  !+-------------------------------------------------------------------+
  subroutine Sopt
    integer :: ix,iy,iz
    real(8) :: sumP1,sumP2
    real(8) :: sumQ1,sumQ2
    real(8) :: sumR1,sumR2
    real(8) :: sumT1,sumT2
    real(8),dimension(MM) :: reSig,imSig
    real(8),dimension(MM) :: reS,imS
    do ix=1,MM
       sumP1=zero
       sumP2=zero
       sumQ1=zero
       sumQ2=zero
       sumR1=zero
       sumR2=zero
       sumT1=zero
       sumT2=zero
       do iy=1,MM
          iz= iy_m_ix(iy,ix)
          if(iz>-L)then
             sumP1=sumP1 + A0m22(iy)*P1(iz)*fmesh 
             sumP2=sumP2 + A0p22(iy)*P2(iz)*fmesh
             sumQ1=sumQ1 + B0m(iy)*Q1(iz)*fmesh
             sumQ2=sumQ2 + B0p(iy)*Q2(iz)*fmesh
             sumR1=sumR1 + C0m(iy)*R1(iz)*fmesh!C0m(iy)*R1(iz)*fmesh
             sumR2=sumR2 + C0p(iy)*R2(iz)*fmesh!C0p(iy)*R2(iz)*fmesh
             sumT1=sumT1 + A0m22(iy)*T1(iz)*fmesh
             sumT2=sumT2 + A0p22(iy)*T2(iz)*fmesh
          end if
       enddo
       imSig(ix)= -(U**2)*(sumP1 + sumP2 - sumQ1 - sumQ2)*pi
       imS(ix)  = -(U**2)*(sumR1 + sumR2 - sumT1 - sumT2)*pi
    enddo
    reSig = kronig(imSig,wr,size(imSig))
    reS   = kronig(imS,wr,size(imS))
    sigma(1,:) = reSig + xi*imSig
    sigma(2,:) = reS   + xi*imS
  end subroutine Sopt



  !******************************************************************
  !******************************************************************
  !******************************************************************






end module IPT_SC_SOPT
