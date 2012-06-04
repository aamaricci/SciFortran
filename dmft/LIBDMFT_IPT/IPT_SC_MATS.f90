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
  real(8)                               :: n,n0
  complex(8)                            :: delta,delta0
  complex(8),dimension(:,:),allocatable :: fg0,sigma
  complex(8),dimension(:),allocatable   :: calG11,calG22,calF12,calF21
  real(8),dimension(:),allocatable      :: calG11t,calG22t
  complex(8),dimension(:),allocatable   :: sigmat,selft,calF12t,calF21t

  public                                :: solve_ipt_sc_matsubara
  public                                :: solve_mpt_sc_matsubara

contains


  !+-------------------------------------------------------------------+
  !PURPOSE  : 
  !+-------------------------------------------------------------------+
  function solve_ipt_sc_matsubara(fg0_,delta_) result(sigma_)
    complex(8),dimension(:,:)                       :: fg0_
    complex(8),dimension(size(fg0_,1),size(fg0_,2)) :: sigma_
    complex(8)                                      :: delta_

    LM=size(fg0_,2)
    if(size(fg0_,1)/=2)call error("Error in solve_ipt_sc_matsubara: size(input,1)!= 2")
    if(loop==1)then
       if(.not.allocated(fg0))allocate(fg0(2,LM))
       if(.not.allocated(sigma))allocate(sigma(2,LM))
       if(.not.allocated(calG11))allocate(calG11(LM))
       if(.not.allocated(calG22))allocate(calG22(LM))
       if(.not.allocated(calF12))allocate(calF12(LM))
       if(.not.allocated(calF21))allocate(calF21(LM))
       !
       if(.not.allocated(sigmat))allocate(sigmat(0:LM))
       if(.not.allocated(selft))allocate(selft(0:LM))
       if(.not.allocated(calG11t))allocate(calG11t(0:LM))
       if(.not.allocated(calG22t))allocate(calG22t(0:LM))
       if(.not.allocated(calF12t))allocate(calF12t(0:LM))
       if(.not.allocated(calF21t))allocate(calF21t(0:LM))
    endif
    fg0 = fg0_
    delta = delta_
    call simpurity
    sigma_(1,:)=sigma(1,:)
    sigma_(2,:)=sigma(2,:) - delta
    loop=loop+1
  end function solve_ipt_sc_matsubara



  !*******************************************************************
  !*******************************************************************
  !*******************************************************************



  !+-------------------------------------------------------------------+
  !PURPOSE  : 
  !+-------------------------------------------------------------------+
  function solve_mpt_sc_matsubara(fg0_,n_,n0_,delta_,delta0_) result(sigma_)
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
       if(.not.allocated(calG11))allocate(calG11(LM))
       if(.not.allocated(calG22))allocate(calG22(LM))
       if(.not.allocated(calF12))allocate(calF12(LM))
       if(.not.allocated(calF21))allocate(calF21(LM))
       !
       if(.not.allocated(sigmat))allocate(sigmat(0:LM))
       if(.not.allocated(selft))allocate(selft(0:LM))
       if(.not.allocated(calG11t))allocate(calG11t(0:LM))
       if(.not.allocated(calG22t))allocate(calG22t(0:LM))
       if(.not.allocated(calF12t))allocate(calF12t(0:LM))
       if(.not.allocated(calF21t))allocate(calF21t(0:LM))
    endif
    fg0 = fg0_
    n     = n_     ; n0    = n0_
    delta = delta_ ; delta0= delta0_
    call simpurity
    A=U**2*n*(1.d0-n)-delta**2
    B=U**2*n0*(1.d0-n0)-delta0**2
    sigma_(1,:) =-u*(n-0.5d0) + sigma(1,:)*A/B
    sigma_(2,:) =-delta       + sigma(2,:)*A/B
    loop=loop+1
  end function solve_mpt_sc_matsubara




  !*******************************************************************
  !*******************************************************************
  !*******************************************************************


  !+-------------------------------------------------------------------+
  !PURPOSE  : Evaluate the 2^nd-order perturbation theory self-energy
  !+-------------------------------------------------------------------+
  subroutine simpurity
    real(8),dimension(0:LM):: dummy12,dummy21
    !Get the HF-corrected Weiss-Fields:
    calG11 =  fg0(1,:)
    calG22 = -conjg(fg0(1,:))
    calF12 =  fg0(2,:)
    calF21 =  conjg(fg0(2,:))

    !FFT to imaginary time:
    call fftgf_iw2tau(calG11,calG11t,beta)
    call fftgf_iw2tau(calG22,calG22t,beta)
    !
    call fftgf_iw2tau(calF12,dummy12,beta,notail=.true.)
    call fftgf_iw2tau(calF21,dummy21,beta,notail=.true.)

    call fftff_iw2tau(calF12,calF12t,beta)
    call fftff_iw2tau(calF21,calF21t,beta)

    open(11,file="calF12_tau.ipt",access="append")    
    open(12,file="calF21_tau.ipt",access="append")    
    open(13,file="DummyF12_tau.ipt")    
    open(14,file="DummyF21_tau.ipt")    
    do i=0,LM
       write(11,*)dble(i)*beta/dble(LM),real(calF12t(i)),aimag(calF12(i))
       write(12,*)dble(i)*beta/dble(LM),real(calF21t(i)),aimag(calF21(i))
       write(13,*)dble(i)*beta/dble(LM),dummy12(i)
       write(14,*)dble(i)*beta/dble(LM),dummy21(i)
    enddo
    do i=11,14
       write(i,*)""
       close(i)
    enddo


    !THIS MAKES THINGS WORKING---> there is a problem in the calculation of Sigma(tau) if using a C-array, rather than a 
    !real one. WHY?!?!?!
    calF12t=dummy12
    calF21t=dummy21


    !Get the 2nd-order Sigma:
    forall(i=0:LM)
       sigmat(i)=  U**2*(calG11t(i)*calG22t(i) - calF12t(i)*calF21t(i))*calG22t(LM-i)
       selft(i) = -U**2*(calF12t(i)*calF21t(i) - calG11t(i)*calG22t(i))*calF12t(i) !ACTHUNG HERE: time inversion simmetry
    end forall

    call fftff_tau2iw(sigmat(0:),sigma(1,:),beta)
    call fftff_tau2iw(selft(0:),sigma(2,:),beta)


    open(11,file="Sigma_tau.ipt",access="append")
    open(12,file="Self_tau.ipt",access="append")
    do i=0,LM
       write(11,*)dble(i)*beta/dble(LM),real(sigmat(i)),aimag(sigmat(i))
       write(12,*)dble(i)*beta/dble(LM),real(selft(i)),aimag(selft(i))
    enddo
    do i=11,12
       write(i,*)""
       close(i)
    enddo

  end subroutine Simpurity


  !******************************************************************
  !******************************************************************
  !******************************************************************




end module IPT_SC_MATS
