!###############################################################
!     PROGRAM  : FUNX_MATS
!     TYPE     : Module
!     PURPOSE  : Contains common routines
!     AUTHORS  : Adriano Amaricci
!###############################################################
module IPT_AF_MATS
  USE IPT_VARS_GLOBAL
  implicit none
  private

  integer            :: i,ispin
  integer,save       :: loop=1 
  type(matsubara_gf) :: fg0(2),sigma(2)
  real(8)            :: n(2),n0(2),xmu0(2)
  public             :: solve_mpt_af_matsubara

contains


  !+-------------------------------------------------------------------+
  !PURPOSE  : 
  !+-------------------------------------------------------------------+
  function solve_mpt_af_matsubara(fg0_,n_,n0_,xmu0_) result(sigma_)
    complex(8),dimension(2,L)    :: fg0_,sigma_
    real(8),dimension(2)         :: n_,n0_,xmu0_
    real(8)                      :: A(2),B(2)

    call allocate_gf(fg0,L)
    call allocate_gf(sigma,L)

    n=n_ ; n0=n0_ ; xmu0=xmu0_
    do ispin=1,2
       fg0(ispin)%iw=fg0_(ispin,:) 
       call fftgf_iw2tau(fg0(ispin)%iw,fg0(ispin)%tau,beta)
    enddo

    forall(ispin=1:2,i=0:L)sigma(ispin)%tau(i)=&
         U**2*fg0(ispin)%tau(i)*fg0(3-ispin)%tau(i)*fg0(3-ispin)%tau(L-i)
    do ispin=1,2
       call fftgf_tau2iw(sigma(ispin)%tau,sigma(ispin)%iw,beta)
       call get_A(ispin)
       call get_B(ispin)
       sigma_(ispin,:) = U*(n(3-ispin)-0.5d0) + &
            A(ispin)*sigma(ispin)%iw/(1.d0-B(ispin)*sigma(ispin)%iw)
    enddo
    call deallocate_gf(fg0)
    call deallocate_gf(sigma)
    loop=loop+1
  contains
    subroutine get_A(idx)
      integer :: idx
      real(8) :: A1,A2
      A1=n(3-idx)*(1.d0-n(3-idx))
      A2=n0(3-idx)*(1.d0-n0(3-idx))
      A(idx)=A1/A2
    end subroutine get_A
    subroutine get_B(idx)
      integer :: idx
      real(8) :: B1,B2
      B1 = (xmu0(idx)-xmu) + U*(1.d0-2.d0*n(3-idx))
      B2 = n0(3-idx)*(1.d0-n0(3-idx))*U**2
      B(idx)=B1/B2
    end subroutine get_B
  end function solve_mpt_af_matsubara



end module IPT_AF_MATS
