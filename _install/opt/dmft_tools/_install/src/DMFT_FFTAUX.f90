MODULE DMFT_FFTAUX
  USE DMFT_FFTGF
  implicit none
  private

  public :: fft_get_density
  public :: tail_coeff_sigma
  public :: tail_coeff_glat
  public :: tail_coeff_g0
  public :: tail_coeff_gk
  public :: tail_coeff_delta

contains

  !+-------------------------------------------------------------------+
  ! PURPOSE  :
  !+-------------------------------------------------------------------+
  function fft_get_density(giw,beta,C) result(dens)
    complex(8),dimension(:)          :: giw
    real(8),dimension(:),allocatable :: gtau
    real(8),dimension(0:3),optional  :: C
    real(8)                          :: beta
    real(8)                          :: dens
    integer                          :: Liw
    Liw=size(giw)
    if(present(C))then
       allocate(gtau(Liw))
       call fft_gf_iw2tau(giw,gtau,beta,C)
    else
       allocate(gtau(0:Liw))
       call fft_gf_iw2tau(giw,gtau(0:),beta)
    endif
    dens = -gtau(Liw)
  end function fft_get_density


  !+-------------------------------------------------------------------+
  ! PURPOSE  : return the tail coefficients given the necessary
  ! informations. This is now limited to the one-band case
  ! - Sigma: S0,S1
  ! - G_lat, Gk, G0: C0,C1,C2,C3
  ! - \Delta: C0,C1
  !+-------------------------------------------------------------------+
  function tail_coeff_sigma(uloc,dens,hfmode) result(S)
    real(8)                :: uloc
    real(8)                :: dens
    logical,optional       :: hfmode
    logical                :: hfmode_
    real(8),dimension(0:1) :: S
    hfmode_=.true.;if(present(hfmode))hfmode_=hfmode
    S(0)=Uloc*dens
    if(hfmode_)S(0)=Uloc*(dens-0.5d0)
    S(1)=Uloc*Uloc*dens*(1d0-dens)
  end function tail_coeff_sigma

  !hloc = <h(k)>
  function tail_coeff_glat(uloc,dens,xmu,hloc,hfmode) result(C)
    real(8)                :: uloc
    real(8)                :: dens,xmu
    real(8)                :: hloc
    logical,optional       :: hfmode
    logical                :: hfmode_
    real(8),dimension(0:1) :: S
    real(8),dimension(0:3) :: C
    hfmode_=.true.;if(present(hfmode))hfmode_=hfmode
    S(0:1) = tail_coeff_sigma(uloc,dens,hfmode_)
    C(0)=0d0
    C(1)=1d0
    C(2)=xmu-hloc-S(0)
    C(3)=S(1)+(xmu-hloc-S(0))*(xmu-hloc-S(0))
  end function tail_coeff_glat

  !ek = h(k)
  function tail_coeff_gk(uloc,dens,xmu,ek,hfmode) result(C)
    real(8)                :: uloc
    real(8)                :: dens,xmu
    real(8)                :: ek
    logical,optional       :: hfmode
    logical                :: hfmode_
    real(8),dimension(0:1) :: S
    real(8),dimension(0:3) :: C
    hfmode_=.true.;if(present(hfmode))hfmode_=hfmode
    S(0:1) = tail_coeff_sigma(uloc,dens,hfmode_)
    C(0)=0d0
    C(1)=1d0
    C(2)=xmu-ek-S(0)
    C(3)=S(1)+(xmu-ek-S(0))*(xmu-ek-S(0))
  end function tail_coeff_gk

  !h1 = <h(k)>
  !h2 = <h(k)**2>
  !hDC = h_{double counting}
  function tail_coeff_g0(xmu,h1,h2,hDC) result(C)
    real(8)                :: dens,xmu
    real(8)                :: h1,h2,hDC
    real(8)                :: mues
    real(8),dimension(0:3) :: C
    mues = xmu-hDC
    C(0) = 0d0
    C(1) = 1d0
    C(2) = mues - h1
    C(3) = mues**2 - 2d0*mues*h1 + h2
  end function tail_coeff_g0

  function tail_coeff_delta(xmu,h1,h2,hDC) result(C)
    real(8)                :: dens,xmu
    real(8)                :: h1,h2,hDC
    real(8)                :: mues
    real(8),dimension(0:1) :: C
    mues = xmu-hDC
    C(0) = 0d0
    C(1) = h1**2 - h2
  end function tail_coeff_delta


END MODULE DMFT_FFTAUX
