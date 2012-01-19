!###############################################################
! PROGRAM  : FUNX_IPT
! TYPE     : Module
! PURPOSE  : Contains common routines
! AUTHORS  : Adriano Amaricci
! COMMENTS : little bit dirty... but seems to work
!###############################################################
module IPT_SC_FUNX_SOPT
  USE IPT_VARS_GLOBAL
  implicit none
  private

  real(8),dimension(:),allocatable      :: A0p11,A0m11,A0p22,A0m22,B0p,B0m,C0p,C0m
  real(8),dimension(:),allocatable      :: P1,P2,Q1,Q2,R1,R2,T1,T2
  integer,allocatable,dimension(:,:)    :: iy_m_ix
  integer,save                          :: loop=1
  complex(8),allocatable,dimension(:,:) :: fg,fg0,sigma
  logical                               :: iprint=.true.
  real(8)                               :: delta
  public :: solve_ipt_sc_sopt

contains

  !+-------------------------------------------------------------------+
  !PROGRAM  : 
  !TYPE     : subroutine
  !PURPOSE  : 
  !+-------------------------------------------------------------------+
  function solve_ipt_sc_sopt(fg0_,fg_,delta_,label_,iprint_) result(sigma_)
    complex(8)                :: fg0_(2,-L:L),sigma_(2,-L:L),fg_(2,-L:L)
    real(8)                   :: delta_
    character(len=*),optional :: label_
    logical,optional          :: iprint_
    if(present(label_))label=trim(adjustl(trim(label_)))
    if(present(iprint_))iprint=iprint_
    if(loop==1)then
       if(.not.allocated(fg))allocate(fg(2,-L:L))
       if(.not.allocated(fg0))allocate(fg0(2,-L:L))
       if(.not.allocated(sigma))allocate(sigma(2,-L:L))
       call get_frequency_index
    endif
    fg0=fg0_ ; fg=fg_ ; delta=delta_
    call simpurity 
    sigma_=sigma
    if(iprint)call print_out()
    loop=loop+1
  end function solve_ipt_sc_sopt
  !******************************************************************
  !******************************************************************
  !******************************************************************









  !+-------------------------------------------------------------------+
  !PROGRAM  : SIMPURITY
  !TYPE     : function
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
  !PROGRAM  : SIMPURITY
  !TYPE     : function
  !PURPOSE  : Evaluate the 2^nd-order perturbation theory self-energy
  !+-------------------------------------------------------------------+
  subroutine getAs
    integer :: i
    real(8) :: w,dos11(-L:L),dos22(-L:L),dosF1(-L:L),dosF2(-L:L)
    complex(8) :: det
    complex(8) :: calG011(-L:L),calF012(-L:L)
    complex(8) :: calF021(-L:L),calG022(-L:L)

    !set the HFB corrected Weiss Field \hatG0^-1 from \calG0^-1 elements
    !\hatG0^-1 = \calG0^-1 - \Sigma_HFB--> -|U|(n-1/2)*\tau_0 - \Delta*\tau_x

    do i=-L,L
       det       =  fg0(1,i)*conjg(fg0(1,-i)) + fg0(2,i)**2!conjg(fg0(2,i))
       calG011(i)=  conjg(fg0(1,-i))/det
       calF012(i)=  fg0(2,i)/det
       calF021(i)=  conjg(fg0(2,i))/det
       calG022(i)= -fg0(1,i)/det
    end do
    call splot("calG011.ipt",wr,calG011(:),append=TT)
    call splot("calF012.ipt",wr,calF012(:),append=TT)
    call splot("calF021.ipt",wr,calF021(:),append=TT)
    call splot("calG022.ipt",wr,calG022(:),append=TT)

    do i=-L,L
       dos11(i) = -aimag(calG011(i))/pi
       dosF1(i) = -aimag(calF012(i))/pi
       dosF2(i) = -aimag(calF021(i))/pi
       dos22(i) = -aimag(calG022(i))/pi
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
  !PROGRAM  : CREATE_FREQUENCY_INDEX
  !TYPE     : function
  !PURPOSE  : Create an array of the indices y-x for a faster evaluation
  !+-------------------------------------------------------------------+
  subroutine get_frequency_index()
    integer :: ix,iy,iz
    allocate(iy_m_ix(-L:L,-L:L))
    iy_m_ix=0
    do ix=-L,L
       do iy=-L,L
          iz = iy - ix !+ L
          if(iz<-L .OR. iz> L) iz=-L-10
          iy_m_ix(iy,ix)=iz
       enddo
    enddo
    if(.not.allocated(A0p11))allocate(A0p11(-L:L))
    if(.not.allocated(A0m11))allocate(A0m11(-L:L))
    if(.not.allocated(A0p22))allocate(A0p22(-L:L))
    if(.not.allocated(A0m22))allocate(A0m22(-L:L))
    if(.not.allocated(B0p))allocate(B0p(-L:L))
    if(.not.allocated(B0m))allocate(B0m(-L:L))
    if(.not.allocated(C0p))allocate(C0p(-L:L))
    if(.not.allocated(C0m))allocate(C0m(-L:L))
    if(.not.allocated(P1))allocate(P1(-L:L))
    if(.not.allocated(P2))allocate(P2(-L:L))
    if(.not.allocated(Q1))allocate(Q1(-L:L))
    if(.not.allocated(Q2))allocate(Q2(-L:L))
    if(.not.allocated(R1))allocate(R1(-L:L))
    if(.not.allocated(R2))allocate(R2(-L:L))
    if(.not.allocated(T1))allocate(T1(-L:L))
    if(.not.allocated(T2))allocate(T2(-L:L))
  end subroutine get_frequency_index
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
    P1=zero
    P2=zero
    Q1=zero
    Q2=zero
    R1=zero
    R2=zero
    T1=zero
    T2=zero
    do ix=-L,L
       do iy=-L,L
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
  !PROGRAM  : SIMPURITY
  !TYPE     : function
  !PURPOSE  : Evaluate the 2^nd-order perturbation theory self-energy
  !+-------------------------------------------------------------------+
  subroutine Sopt
    integer :: ix,iy,iz
    real(8) :: sumP1,sumP2
    real(8) :: sumQ1,sumQ2
    real(8) :: sumR1,sumR2
    real(8) :: sumT1,sumT2
    real(8),dimension(-L:L) :: reSig,imSig
    real(8),dimension(-L:L) :: reS,imS
    do ix=-L,L
       sumP1=zero
       sumP2=zero
       sumQ1=zero
       sumQ2=zero
       sumR1=zero
       sumR2=zero
       sumT1=zero
       sumT2=zero
       do iy=-L,L
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
    sigma(2,:) = reS   + xi*imS   - delta
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
    real(8)    :: nimp
    nimp = -2.d0*sum(aimag(fg(1,-L:L))*fermi(wr(-L:L),beta))*fmesh/pi
    call splot("nVSiloop.ipt"//trim(adjustl(trim(label))),iloop,nimp,append=TT)
    call splot("deltaVSiloop.ipt"//trim(adjustl(trim(label))),iloop,delta,append=TT)

    call splot("Gloc_realw.ipt"//trim(adjustl(trim(label))),wr,fg(1,:),append=printf)
    call splot("Floc_realw.ipt"//trim(adjustl(trim(label))),wr,fg(2,:),append=printf)
    call splot("DOS.ipt"//trim(adjustl(trim(label))),wr,-aimag(fg(1,:))/pi,append=printf)

    call splot("Sigma_realw.ipt"//trim(adjustl(trim(label))),wr,sigma(1,:),append=printf)
    call splot("Self_realw.ipt"//trim(adjustl(trim(label))),wr,sigma(2,:),append=printf)

    call splot("calG0_realw.ipt"//trim(adjustl(trim(label))),wr,fg0(1,:),append=printf)
    call splot("calF0_realw.ipt"//trim(adjustl(trim(label))),wr,fg0(2,:),append=printf)
    return
  end subroutine print_out
  !******************************************************************
  !******************************************************************
  !******************************************************************

end module IPT_SC_FUNX_SOPT
