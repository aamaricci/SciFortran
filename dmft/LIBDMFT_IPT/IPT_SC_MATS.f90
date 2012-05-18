!###############################################################
!     PURPOSE  : Contains common routines
!     AUTHORS  : Adriano Amaricci
!###############################################################
module IPT_SC_MATS
  USE IPT_VARS_GLOBAL
  implicit none
  private

  integer                         :: i,LM
  integer,save                    :: loop=1
  type(matsubara_gf),dimension(2) :: fg0,sigma
  type(matsubara_gf)              :: calG11,calG22,calF
  real(8)                         :: n,n0,delta,delta0
  public                          :: solve_ipt_sc_matsubara
  public                          :: solve_mpt_sc_matsubara

contains


  !+-------------------------------------------------------------------+
  !PURPOSE  : 
  !+-------------------------------------------------------------------+
  function solve_ipt_sc_matsubara(fg0_,delta_) result(sigma_)
    complex(8),dimension(:,:)                       :: fg0_
    complex(8),dimension(size(fg0_,1),size(fg0_,2)) :: sigma_
    real(8)                                         :: delta_
    LM=size(fg0_,2)
    if(size(fg0_,1)/=2)call error("Error in solve_ipt_sc_matsubara: size(input,1) != 2")
    if(loop==1)then
       if(.not.fg0(1)%status)call allocate_gf(fg0(1),LM)
       if(.not.fg0(2)%status)call allocate_gf(fg0(2),LM)
       if(.not.sigma(1)%status)call allocate_gf(sigma(1),LM)
       if(.not.sigma(2)%status)call allocate_gf(sigma(2),LM)
       if(.not.calG11%status)call allocate_gf(calG11,LM)
       if(.not.calG22%status)call allocate_gf(calG22,LM)
       if(.not.calF%status)call allocate_gf(calF,LM)
    endif
    fg0(1)%iw = fg0_(1,:); fg0(2)%iw = fg0_(2,:) ; delta = delta_
    call simpurity
    sigma_(1,:)=sigma(1)%iw
    sigma_(2,:)=sigma(2)%iw - delta
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
    real(8)                                         :: n_,n0_,delta_,delta0_,A,B
    LM=size(fg0_,2)
    if(size(fg0_,1)/=2)call error("Error in solve_ipt_sc_matsubara: size(input,1) != 2")
    if(loop==1)then
       if(.not.fg0(1)%status)call allocate_gf(fg0(1),LM)
       if(.not.fg0(2)%status)call allocate_gf(fg0(2),LM)
       if(.not.sigma(1)%status)call allocate_gf(sigma(1),LM)
       if(.not.sigma(2)%status)call allocate_gf(sigma(2),LM)
       if(.not.calG11%status)call allocate_gf(calG11,LM)
       if(.not.calG22%status)call allocate_gf(calG22,LM)
       if(.not.calF%status)call allocate_gf(calF,LM)
    endif
    fg0(1)%iw = fg0_(1,:) ; fg0(2)%iw = fg0_(2,:)
    n     = n_     ; n0    = n0_
    delta = delta_ ; delta0= delta0_
    call simpurity
    A=U**2*n*(1.d0-n)-delta**2
    B=U**2*n0*(1.d0-n0)-delta0**2
    sigma_(1,:) =-u*(n-0.5d0) + sigma(1)%iw*A/B
    sigma_(2,:) =-delta       + sigma(2)%iw*A/B
    loop=loop+1
  end function solve_mpt_sc_matsubara




  !*******************************************************************
  !*******************************************************************
  !*******************************************************************


  !+-------------------------------------------------------------------+
  !PURPOSE  : Evaluate the 2^nd-order perturbation theory self-energy
  !+-------------------------------------------------------------------+
  subroutine simpurity
    calG11%iw =  fg0(1)%iw
    calG22%iw = -conjg(fg0(1)%iw)
    calF%iw   =  fg0(2)%iw
    call fftgf_iw2tau(calG11%iw,calG11%tau,beta)
    call fftgf_iw2tau(calG22%iw,calG22%tau,beta)
    call fftgf_iw2tau(calF%iw,calF%tau,beta,notail=.true.)

    forall(i=0:LM)
       sigma(1)%tau(i)=  U**2*(calG11%tau(i)*calG22%tau(i) - calF%tau(i)**2)*calG22%tau(LM-i)
       sigma(2)%tau(i)= -U**2*(calF%tau(i)**2 - calG11%tau(i)*calG22%tau(i))*calF%tau(i)
    end forall


    call fftgf_tau2iw(sigma(1)%tau,sigma(1)%iw,beta)
    call fftgf_tau2iw(sigma(2)%tau,sigma(2)%iw,beta) !,normal=.false.) keep in mind


    ! open(11,file="Sigma_tau.ipt",access="append")
    ! open(12,file="Self_tau.ipt",access="append")
    ! ! open(13,file="Sigma2_iw.ipt",access="append")
    ! ! open(14,file="Self2_iw.ipt",access="append")
    ! do i=0,LM
    !    write(11,*)dble(i)*beta/dble(LM),sigma(1)%tau(i)
    !    write(12,*)dble(i)*beta/dble(LM),sigma(2)%tau(i)
    ! enddo
    ! ! do i=1,LM
    ! !    write(13,*)pi/beta*dble(2*i-1),aimag(sigma(1)%iw(i)),real(sigma(1)%iw(i))
    ! !    write(14,*)pi/beta*dble(2*i-1),aimag(sigma(2)%iw(i)),real(sigma(2)%iw(i))
    ! ! enddo
    ! do i=11,12
    !    write(i,*)""
    !    close(i)
    ! enddo

  end subroutine Simpurity


  !******************************************************************
  !******************************************************************
  !******************************************************************




end module IPT_SC_MATS
