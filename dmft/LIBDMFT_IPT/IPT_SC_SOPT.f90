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
  real(8)                               :: n,n0,delta,delta0,mesh
  integer                               :: Lm
  public                                :: solve_ipt_sc_sopt
  public                                :: solve_mpt_sc_sopt


contains


  !+-------------------------------------------------------------------+
  !PURPOSE  : 
  !+-------------------------------------------------------------------+
  function solve_ipt_sc_sopt(fg0_,wr_,delta_,L) result(sigma_)
    integer                   :: L
    complex(8),dimension(2,L) :: fg0_
    complex(8),dimension(2,L) :: sigma_
    real(8),dimension(L)      :: wr_
    real(8)                   :: delta_
    Lm=L
    if(loop==1)then
       if(.not.allocated(fg0))allocate(fg0(2,Lm))
       if(.not.allocated(sigma))allocate(sigma(2,Lm))
       if(.not.allocated(wr))allocate(wr(Lm))
       call get_frequency_index
    endif
    fg0=fg0_ ; delta=delta_ ; wr=wr_ ; mesh=abs(wr(2)-wr(1))
    call simpurity 
    sigma_(1,:) =          sigma(1,:)
    sigma_(2,:) = -delta + sigma(2,:)
    loop=loop+1
  end function solve_ipt_sc_sopt




  !******************************************************************
  !******************************************************************
  !******************************************************************




  !+-------------------------------------------------------------------+
  !PURPOSE  : 
  !+-------------------------------------------------------------------+
  function solve_mpt_sc_sopt(fg0_,wr_,n_,n0_,delta_,delta0_,L) result(sigma_)
    integer                   :: L
    complex(8),dimension(2,L) :: fg0_
    complex(8),dimension(2,L) :: sigma_
    real(8),dimension(L)      :: wr_
    real(8)                   :: n_,n0_,delta_,delta0_
    real(8)                   :: A,B
    Lm=L
    if(loop==1) then
       if(.not.allocated(fg0))allocate(fg0(2,Lm))
       if(.not.allocated(sigma))allocate(sigma(2,Lm))
       if(.not.allocated(wr))allocate(wr(Lm))
       call get_frequency_index  
    endif
    fg0=fg0_  ; wr=wr_ 
    n=n_ ; n0=n0_ ; delta=delta_ ; delta0=delta0_
    mesh=abs(wr(2)-wr(1))
    call simpurity
    A = U**2*n*(1.d0-n)-delta**2
    B = U**2*n0*(1.d0-n0)-delta0**2
    sigma_(1,:) = -u*(n-0.5d0)   + sigma(1,:)*A/B
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
    allocate(iy_m_ix(Lm,Lm))
    iy_m_ix=0
    do ix=1,Lm
       do iy=1,Lm
          iz = iy - ix + Lm/2
          ! if(iz<-Lm .OR. iz>Lm) iz=-1 !out of range, the same as the old if(iz>-L)
          iy_m_ix(iy,ix)=iz
       enddo
    enddo
    if(.not.allocated(A0p11))allocate(A0p11(Lm))
    if(.not.allocated(A0m11))allocate(A0m11(Lm))
    if(.not.allocated(A0p22))allocate(A0p22(Lm))
    if(.not.allocated(A0m22))allocate(A0m22(Lm))
    if(.not.allocated(B0p))allocate(B0p(Lm))
    if(.not.allocated(B0m))allocate(B0m(Lm))
    if(.not.allocated(C0p))allocate(C0p(Lm))
    if(.not.allocated(C0m))allocate(C0m(Lm))
    if(.not.allocated(P1))allocate(P1(Lm))
    if(.not.allocated(P2))allocate(P2(Lm))
    if(.not.allocated(Q1))allocate(Q1(Lm))
    if(.not.allocated(Q2))allocate(Q2(Lm))
    if(.not.allocated(R1))allocate(R1(Lm))
    if(.not.allocated(R2))allocate(R2(Lm))
    if(.not.allocated(T1))allocate(T1(Lm))
    if(.not.allocated(T2))allocate(T2(Lm))
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
    real(8) :: w,dos11(Lm),dos22(Lm),dosF1(Lm),dosF2(Lm)
    do i=1,Lm
       dos11(i) = -aimag(fg0(1,i))/pi
       dos22(i) = -aimag(-conjg(fg0(1,Lm+1-i)))/pi
       dosF1(i) = -aimag(fg0(2,i))/pi
       dosF2(i) = -aimag(fg0(2,i))/pi
       !OLD
       ! dos11(i) = -aimag(fg0(1,i))/pi
       ! dos22(i) = -aimag(-conjg(fg0(1,Lm+1-i)))/pi
       ! dosF1(i) = -aimag(fg0(2,i))/pi
       ! dosF2(i) = -aimag(fg0(2,i))/pi
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
    do ix=1,Lm
       do iy=1,Lm
          iz= iy_m_ix(iy,ix)
          ! posso eliminare il check di prima perche i punti fuori  dai 
          ! boundaries non verranno calcolati comunque
          if((iz>=1).and.(iz<=Lm)) then
             P2(ix)=P2(ix) + A0p11(iy)*A0m22(iz)*mesh
             P1(ix)=P1(ix) + A0m11(iy)*A0p22(iz)*mesh
             Q2(ix)=Q2(ix) + C0p(iy)*A0m22(iz)*mesh     
             Q1(ix)=Q1(ix) + C0m(iy)*A0p22(iz)*mesh     
             R2(ix)=R2(ix) + B0p(iy)*B0m(iz)*mesh
             R1(ix)=R1(ix) + B0m(iy)*B0p(iz)*mesh
             T2(ix)=T2(ix) + A0p11(iy)*B0m(iz)*mesh
             T1(ix)=T1(ix) + A0m11(iy)*B0p(iz)*mesh
          endif

          !OLDER
          !if(iz>0)then
          ! P2(ix)=P2(ix) + A0p11(iy)*A0m22(iz)*mesh
          ! P1(ix)=P1(ix) + A0m11(iy)*A0p22(iz)*mesh
          ! Q2(ix)=Q2(ix) + C0p(iy)*A0m22(iz)*mesh
          ! Q1(ix)=Q1(ix) + C0m(iy)*A0p22(iz)*mesh
          ! R2(ix)=R2(ix) + B0p(iy)*B0m(iz)*mesh
          ! R1(ix)=R1(ix) + B0m(iy)*B0p(iz)*mesh
          ! T2(ix)=T2(ix) + A0p11(iy)*B0m(iz)*mesh
          ! T1(ix)=T1(ix) + A0m11(iy)*B0p(iz)*mesh

          !VERY OLDER VERSION:
          ! P1(ix)=P1(ix) + A0p11(iy)*A0m22(iz)*mesh
          ! P2(ix)=P2(ix) + A0m11(iy)*A0p22(iz)*mesh
          ! Q1(ix)=Q1(ix) + C0p(iy)*A0m22(iz)*mesh !C0p(iy)*A0m22(iz)*mesh
          ! Q2(ix)=Q2(ix) + C0m(iy)*A0p22(iz)*mesh !C0m(iy)*A0p22(iz)*mesh
          ! R1(ix)=R1(ix) + B0p(iy)*B0m(iz)*mesh
          ! R2(ix)=R2(ix) + B0m(iy)*B0p(iz)*mesh
          ! T1(ix)=T1(ix) + A0p11(iy)*B0m(iz)*mesh
          ! T2(ix)=T2(ix) + A0m11(iy)*B0p(iz)*mesh
          !endif
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
    real(8),dimension(Lm) :: reSig,imSig
    real(8),dimension(Lm) :: reS,imS
    do ix=1,Lm
       sumP1=zero
       sumP2=zero
       sumQ1=zero
       sumQ2=zero
       sumR1=zero
       sumR2=zero
       sumT1=zero
       sumT2=zero
       do iy=1,Lm
          iz= iy_m_ix(iy,ix)

          if((iz>=1).and.(iz<=Lm)) then         ! in questo modo il contributo fuori da dove le polarizzazioni sono definite e' eliminato 
             sumP1=sumP1 + A0m22(Lm+1-iz)*P1(iy)*mesh ! Adriano mi ha spiegato il senso di com'era prima ma  
             sumP2=sumP2 + A0p22(Lm+1-iz)*P2(iy)*mesh ! questo mi sembra MORALMENTE 
             sumQ1=sumQ1 + B0m(Lm+1-iz)*Q1(iy)*mesh   ! piu' chiaro per i posteri per cui... 
             sumQ2=sumQ2 + B0p(Lm+1-iz)*Q2(iy)*mesh
             sumR1=sumR1 + C0m(Lm+1-iz)*R1(iy)*mesh!C0m(iy)*R1(iz)*mesh   !! DOUBLE CHECK
             sumR2=sumR2 + C0p(Lm+1-iz)*R2(iy)*mesh!C0p(iy)*R2(iz)*mesh
             sumT1=sumT1 + A0m22(Lm+1-iz)*T1(iy)*mesh
             sumT2=sumT2 + A0p22(Lm+1-iz)*T2(iy)*mesh
          endif

          !OLDER
          ! if(iz>0)then
          !    sumP1=sumP1 + A0m22(Lm-iz+1)*P1(iy)*mesh 
          !    sumP2=sumP2 + A0p22(Lm-iz+1)*P2(iy)*mesh
          !    sumQ1=sumQ1 + B0m(Lm-iz+1)*Q1(iy)*mesh
          !    sumQ2=sumQ2 + B0p(Lm-iz+1)*Q2(iy)*mesh
          !    sumR1=sumR1 + C0m(Lm-iz+1)*R1(iy)*mesh
          !    sumR2=sumR2 + C0p(Lm-iz+1)*R2(iy)*mesh
          !    sumT1=sumT1 + A0m22(Lm-iz+1)*T1(iy)*mesh
          !    sumT2=sumT2 + A0p22(Lm-iz+1)*T2(iy)*mesh

          !    !OLDER++ VERSION:
          !    ! sumP1=sumP1 + A0m22(iy)*P1(iz)*mesh 
          !    ! sumP2=sumP2 + A0p22(iy)*P2(iz)*mesh
          !    ! sumQ1=sumQ1 + B0m(iy)*Q1(iz)*mesh
          !    ! sumQ2=sumQ2 + B0p(iy)*Q2(iz)*mesh
          !    ! sumR1=sumR1 + C0m(iy)*R1(iz)*mesh
          !    ! sumR2=sumR2 + C0p(iy)*R2(iz)*mesh
          !    ! sumT1=sumT1 + A0m22(iy)*T1(iz)*mesh
          !    ! sumT2=sumT2 + A0p22(iy)*T2(iz)*mesh
          ! end if


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
