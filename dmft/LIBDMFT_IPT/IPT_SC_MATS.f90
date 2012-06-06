!###############################################################
!     PURPOSE  : Contains common routines
!     AUTHORS  : Adriano Amaricci
!###############################################################
module IPT_SC_MATS
  USE IPT_VARS_GLOBAL
  implicit none
  private

  integer                               :: i,LM
  integer,save                          :: loop=1
  !occupations:
  real(8)                               :: n,n0
  !order parameters (Real and Complex)
  real(8)                               :: rdelta,rdelta0
  complex(8)                            :: cdelta,cdelta0
  !I/O arrays to be shared between routines:
  complex(8),dimension(:,:),allocatable :: fg0,sigma

  !Common:
  complex(8),dimension(:),allocatable   :: calG11,calG22,calF
  real(8),dimension(:),allocatable      :: calG11t,calG22t
  real(8),dimension(:),allocatable      :: Sigmat
  !Real order parameter:
  real(8),dimension(:),allocatable      :: RcalFt
  real(8),dimension(:),allocatable      :: Rselft
  !Complex order parameter:
  complex(8),dimension(:),allocatable   :: CcalFt
  complex(8),dimension(:),allocatable   :: Cselft


  interface solve_ipt_sc_matsubara
     module procedure solve_ipt_sc_matsubara_r,solve_ipt_sc_matsubara_c
  end interface solve_ipt_sc_matsubara

  interface solve_mpt_sc_matsubara
     module procedure solve_mpt_sc_matsubara_r,solve_mpt_sc_matsubara_c
  end interface solve_mpt_sc_matsubara

  public                                :: solve_ipt_sc_matsubara
  public                                :: solve_mpt_sc_matsubara

contains

  !+-------------------------------------------------------------------+
  !PURPOSE  : 
  !+-------------------------------------------------------------------+
  function solve_ipt_sc_matsubara_r(fg0_,delta_) result(sigma_)
    complex(8),dimension(:,:)                       :: fg0_
    complex(8),dimension(size(fg0_,1),size(fg0_,2)) :: sigma_
    real(8)                                         :: delta_
    LM=size(fg0_,2)
    if(size(fg0_,1)/=2)call error("Error in solve_ipt_sc_matsubara: size(input,1)!= 2")
    if(loop==1)then
       if(.not.allocated(fg0))allocate(fg0(2,LM))
       if(.not.allocated(sigma))allocate(sigma(2,LM))
       !
       if(.not.allocated(calG11))allocate(calG11(LM))
       if(.not.allocated(calG22))allocate(calG22(LM))
       if(.not.allocated(calF))allocate(calF(LM))
       !
       if(.not.allocated(calG11t))allocate(calG11t(0:LM))
       if(.not.allocated(calG22t))allocate(calG22t(0:LM))
       if(.not.allocated(RcalFt))allocate(RcalFt(0:LM))
       !
       if(.not.allocated(sigmat))allocate(sigmat(0:LM))
       if(.not.allocated(Rselft))allocate(Rselft(0:LM))
    endif
    fg0 = fg0_
    rdelta = delta_
    call simpurity_r
    sigma_(1,:)=sigma(1,:)
    sigma_(2,:)=sigma(2,:) - rdelta
    loop=loop+1
  end function solve_ipt_sc_matsubara_r

  function solve_ipt_sc_matsubara_c(fg0_,delta_) result(sigma_)
    complex(8),dimension(:,:)                       :: fg0_
    complex(8),dimension(size(fg0_,1),size(fg0_,2)) :: sigma_
    complex(8)                                      :: delta_
    LM=size(fg0_,2)
    if(size(fg0_,1)/=2)call error("Error in solve_ipt_sc_matsubara: size(input,1)!= 2")
    if(loop==1)then
       if(.not.allocated(fg0))allocate(fg0(2,LM))
       if(.not.allocated(sigma))allocate(sigma(2,LM))
       !
       if(.not.allocated(calG11))allocate(calG11(LM))
       if(.not.allocated(calG22))allocate(calG22(LM))
       if(.not.allocated(calF))allocate(calF(LM))
       !
       if(.not.allocated(calG11t))allocate(calG11t(0:LM))
       if(.not.allocated(calG22t))allocate(calG22t(0:LM))
       if(.not.allocated(CcalFt))allocate(CcalFt(0:LM))
       !
       if(.not.allocated(sigmat))allocate(sigmat(0:LM))
       if(.not.allocated(Cselft))allocate(Cselft(0:LM))
    endif
    fg0 = fg0_
    Cdelta = delta_
    call simpurity_c
    sigma_(1,:)=sigma(1,:)
    sigma_(2,:)=sigma(2,:) - cdelta
    loop=loop+1
  end function solve_ipt_sc_matsubara_c

  function solve_mpt_sc_matsubara_r(fg0_,n_,n0_,delta_,delta0_) result(sigma_)
    complex(8),dimension(:,:)                       :: fg0_
    complex(8),dimension(size(fg0_,1),size(fg0_,2)) :: sigma_
    real(8)                                         :: n_,n0_
    real(8)                                         :: delta_,delta0_
    real(8)                                         :: A,B
    LM=size(fg0_,2)
    if(size(fg0_,1)/=2)call error("Error in solve_ipt_sc_matsubara: size(input,1) != 2")
    if(loop==1)then
       if(.not.allocated(fg0))allocate(fg0(2,LM))
       if(.not.allocated(sigma))allocate(sigma(2,LM))
       !
       if(.not.allocated(calG11))allocate(calG11(LM))
       if(.not.allocated(calG22))allocate(calG22(LM))
       if(.not.allocated(calF))allocate(calF(LM))
       !
       if(.not.allocated(calG11t))allocate(calG11t(0:LM))
       if(.not.allocated(calG22t))allocate(calG22t(0:LM))
       if(.not.allocated(RcalFt))allocate(RcalFt(0:LM))
       !
       if(.not.allocated(sigmat))allocate(sigmat(0:LM))
       if(.not.allocated(Rselft))allocate(Rselft(0:LM))
    endif
    fg0 = fg0_
    n     = n_      ; n0    = n0_
    rdelta = delta_ ; rdelta0= delta0_
    call simpurity_r
    A=U**2*n*(1.d0-n)-rdelta**2
    B=U**2*n0*(1.d0-n0)-rdelta0**2
    sigma_(1,:) =-u*(n-0.5d0) + sigma(1,:)*A/B
    sigma_(2,:) =-rdelta       + sigma(2,:)*A/B
    loop=loop+1
  end function solve_mpt_sc_matsubara_r

  function solve_mpt_sc_matsubara_c(fg0_,n_,n0_,delta_,delta0_) result(sigma_)
    complex(8),dimension(:,:)                       :: fg0_
    complex(8),dimension(size(fg0_,1),size(fg0_,2)) :: sigma_
    real(8)                                         :: n_,n0_
    complex(8)                                      :: delta_,delta0_
    real(8)                                         :: A,B
    LM=size(fg0_,2)
    if(size(fg0_,1)/=2)call error("Error in solve_ipt_sc_matsubara: size(input,1) != 2")
    if(loop==1)then
       if(.not.allocated(fg0))allocate(fg0(2,LM))
       if(.not.allocated(sigma))allocate(sigma(2,LM))
       !
       if(.not.allocated(calG11))allocate(calG11(LM))
       if(.not.allocated(calG22))allocate(calG22(LM))
       if(.not.allocated(calF))allocate(calF(LM))
       !
       if(.not.allocated(calG11t))allocate(calG11t(0:LM))
       if(.not.allocated(calG22t))allocate(calG22t(0:LM))
       if(.not.allocated(CcalFt))allocate(CcalFt(0:LM))
       !
       if(.not.allocated(sigmat))allocate(sigmat(0:LM))
       if(.not.allocated(Cselft))allocate(Cselft(0:LM))
    endif
    fg0 = fg0_
    n     = n_      ; n0    = n0_
    Cdelta = delta_ ; Cdelta0= delta0_
    call simpurity_c
    !This is not obvious, just a guess now but I need to check!!
    A=U**2*n*(1.d0-n)-abs(Cdelta)**2
    B=U**2*n0*(1.d0-n0)-abs(Cdelta0)**2
    sigma_(1,:) =-u*(n-0.5d0) + sigma(1,:)*A/B
    sigma_(2,:) =-Cdelta       + sigma(2,:)*A/B
    loop=loop+1
  end function solve_mpt_sc_matsubara_c




  !*******************************************************************
  !*******************************************************************
  !*******************************************************************



  !+-------------------------------------------------------------------+
  !PURPOSE  : Evaluate the 2^nd-order perturbation theory self-energy
  !+-------------------------------------------------------------------+
  subroutine Simpurity_r
    calG11 =  fg0(1,:)
    calG22 = -conjg(fg0(1,:))
    calF   =  fg0(2,:)
    call fftgf_iw2tau(calG11,calG11t,beta)
    call fftgf_iw2tau(calG22,calG22t,beta)
    call fftgf_iw2tau(calF,RcalFt(0:),beta,notail=.true.)
    forall(i=0:LM)
       sigmat(i)=  U**2*(calG11t(i)*calG22t(i) - RcalFt(i)**2)*calG22t(LM-i)
       Rselft(i)= -U**2*(RcalFt(i)**2 - calG11t(i)*calG22t(i))*RcalFt(i)
    end forall
    call fftgf_tau2iw(sigmat,sigma(1,:),beta)
    call fftgf_tau2iw(Rselft,sigma(2,:),beta)
  end subroutine Simpurity_r

  subroutine Simpurity_C
    !Get the HF-corrected Weiss-Fields:
    calG11 =  fg0(1,:)
    calG22 = -conjg(fg0(1,:))
    calF   =  fg0(2,:)
    call fftgf_iw2tau(calG11,calG11t,beta)
    call fftgf_iw2tau(calG22,calG22t,beta)
    call fftff_iw2tau(calF,CcalFt(0:),beta)
    !Get the 2nd-order Sigma:
    forall(i=0:LM)
       sigmat(i)=  U**2*(calG11t(i)*calG22t(i) - abs(CcalFt(i))**2)*calG22t(LM-i)
       Cselft(i) = -U**2*(abs(CcalFt(i))**2 - calG11t(i)*calG22t(i))*CcalFt(i) !ACTHUNG HERE: time inversion simmetry
    end forall
    call fftgf_tau2iw(sigmat(0:),sigma(1,:),beta)
    call fftff_tau2iw(Cselft(0:),sigma(2,:),beta)
    ! open(11,file="Sigma_tau.ipt",access="append")
    ! open(12,file="Self_tau.ipt",access="append")
    ! do i=0,LM
    !    write(11,*)dble(i)*beta/dble(LM),sigmat(i)
    !    write(12,*)dble(i)*beta/dble(LM),real(selft(i)),aimag(selft(i))
    ! enddo
    ! do i=11,12
    !    write(i,*)""
    !    close(i)
    ! enddo
  end subroutine Simpurity_C




  !******************************************************************
  !******************************************************************
  !******************************************************************




end module IPT_SC_MATS
