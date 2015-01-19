module SF_SPECIAL
  implicit none
  private

  !FUNCTIONS:
  public :: heaviside
  interface step
     module procedure step_x,step_ij
  end interface step
  public :: step

  public :: fermi
  interface sgn
     module procedure i_sgn,d_sgn
  end interface sgn
  public :: sgn

  !ERROR FUNCS
  public :: wfun         !complex error function (Faddeeva function)
  public :: zerf   

  !BETHE:
  public :: gfbethe
  public :: gfbether
  public :: bethe_lattice
  public :: dens_bethe
  public :: bethe_guess_g0

  !HYPERCUBIC/GAUSSIAN DENS:
  public :: dens_hyperc

  !2D-SQUARE LATTICE ANALYTIC DOS:
  public :: dens_2dsquare

  !SPECIAL FUNCTIONS:
  public :: airya           ! AIRYA computes Airy functions and their derivatives.
  public :: airyb           ! AIRYB computes Airy functions and their derivatives.
  public :: airyzo          ! AIRYZO computes the first NT zeros of Ai(x) and Ai'(x).
  public :: ajyik           ! AJYIK computes Bessel functions Jv(x), Yv(x), Iv(x), Kv(x).
  public :: aswfa           ! ASWFA: prolate and oblate spheroidal angular functions of the first kind.
  public :: aswfb           ! ASWFB: prolate and oblate spheroidal angular functions of the first kind.
  public :: bernoa          ! BERNOA computes the Bernoulli number Bn.
  public :: bernob          ! BERNOB computes the Bernoulli number Bn.
  public :: betaf           ! BETA computes the Beta function B(p,q).
  public :: bjndd           ! BJNDD computes Bessel functions Jn(x) and first and second derivatives.
  public :: cbk             ! CBK computes coefficients for oblate radial functions with small argument.
  public :: cchg            ! CCHG computes the confluent hypergeometric function.
  public :: cerf            ! CERF computes the error function and derivative for a complex argument.
  public :: cerror          ! CERROR computes the error function for a complex argument.
  public :: cerzo           ! CERZO evaluates the complex zeros of the error function.
  public :: cfc             ! CFC computes the complex Fresnel integral C(z) and C'(z).
  public :: cfs             ! CFS computes the complex Fresnel integral S(z) and S'(z).
  public :: cgama           ! CGAMA computes the Gamma function for complex argument.
  public :: ch12n           ! CH12N computes Hankel functions of first and second kinds, complex argument.
  public :: chgm            ! CHGM computes the confluent hypergeometric function M(a,b,x).
  public :: chgu            ! CHGU computes the confluent hypergeometric function U(a,b,x).
  public :: chgubi          ! CHGUBI: confluent hypergeometric function with integer argument B.
  public :: chguit          ! CHGUIT computes the hypergeometric function using Gauss-Legendre integration.
  public :: chgul           ! CHGUL: confluent hypergeometric function U(a,b,x) for large argument X.
  public :: chgus           ! CHGUS: confluent hypergeometric function U(a,b,x) for small argument X.
  public :: cik01           ! CIK01: modified Bessel I0(z), I1(z), K0(z) and K1(z) for complex argument.
  public :: ciklv           ! CIKLV: modified Bessel functions Iv(z), Kv(z), complex argument, large order.
  public :: cikna           ! CIKNA: modified Bessel functions In(z), Kn(z), derivatives, complex argument.
  public :: ciknb           ! CIKNB computes complex modified Bessel functions In(z) and Kn(z).
  public :: cikva           ! CIKVA: modified Bessel functions Iv(z), Kv(z), arbitrary order, complex.
  public :: cikvb           ! CIKVB: modified Bessel functions,Iv(z), Kv(z), arbitrary order, complex.
  public :: cisia           ! CISIA computes cosine Ci(x) and sine integrals Si(x).
  public :: cisib           ! CISIB computes cosine and sine integrals.
  public :: cjk             ! CJK: asymptotic expansion coefficients for Bessel functions of large order.
  public :: cjy01           ! CJY01: complexBessel functions, derivatives, J0(z), J1(z), Y0(z), Y1(z).
  public :: cjylv           ! CJYLV: Bessel functions Jv(z), Yv(z) of complex argument and large order v.
  public :: cjyna           ! CJYNA: Bessel functions and derivatives, Jn(z) and Yn(z) of complex argument.
  public :: cjynb           ! CJYNB: Bessel functions, derivatives, Jn(z) and Yn(z) of complex argument.
  public :: cjyva           ! CJYVA: Bessel functions and derivatives, Jv(z) and Yv(z) of complex argument.
  public :: cjyvb           ! CJYVB: Bessel functions and derivatives, Jv(z) and Yv(z) of complex argument.
  public :: clpmn           ! CLPMN: associated Legendre functions and derivatives for complex argument.
  public :: clpn            ! CLPN computes Legendre functions and derivatives for complex argument.
  public :: clqmn           ! CLQMN: associated Legendre functions and derivatives for complex argument.
  public :: clqn            ! CLQN: Legendre function Qn(z) and derivative Wn'(z) for complex argument.
  public :: comelp          ! COMELP computes complete elliptic integrals K(k) and E(k).
  public :: cpbdn           ! CPBDN: parabolic cylinder function Dn(z) and Dn'(z) for complex argument.
  public :: cpdla           ! CPDLA computes complex parabolic cylinder function Dn(z) for large argument.
  public :: cpdsa           ! CPDSA computes complex parabolic cylinder function Dn(z) for small argument.
  public :: cpsi            ! CPSI computes the psi function for a complex argument.
  public :: csphik          ! CSPHIK: complex modified spherical Bessel functions and derivatives.
  public :: csphjy          ! CSPHJY: spherical Bessel functions jn(z) and yn(z) for complex argument.
  public :: cv0             ! CV0 computes the initial characteristic value of Mathieu functions.
  public :: cva1            ! CVA1 computes a sequence of characteristic values of Mathieu functions.
  public :: cva2            ! CVA2 computes a specific characteristic value of Mathieu functions.
  public :: cvf             ! CVF computes F for the characteristic equation of Mathieu functions.
  public :: cvql            ! CVQL computes the characteristic value of Mathieu functions for q <= 3*m.
  public :: cvqm            ! CVQM computes the characteristic value of Mathieu functions for q <= m*m.
  public :: cy01            ! CY01 computes complex Bessel functions Y0(z) and Y1(z) and derivatives.
  public :: cyzo            ! CYZO computes zeros of complex Bessel functions Y0(z) and Y1(z) and Y1'(z).
  public :: dvla            ! DVLA computes parabolic cylinder functions Dv(x) for large argument.
  public :: dvsa            ! DVSA computes parabolic cylinder functions Dv(x) for small argument.
  public :: e1xa            ! E1XA computes the exponential integral E1(x).
  public :: e1xb            ! E1XB computes the exponential integral E1(x).
  public :: e1z             ! E1Z computes the complex exponential integral E1(z).
  public :: eix             ! EIX computes the exponential integral Ei(x).
  public :: elit            ! ELIT: complete and incomplete elliptic integrals F(k,phi) and E(k,phi).
  public :: elit3           ! ELIT3 computes the elliptic integral of the third kind.
  public :: enxa            ! ENXA computes the exponential integral En(x).
  public :: enxb            ! ENXB computes the exponential integral En(x).
  public :: werror          ! ERROR evaluates the error function.
  public :: eulera          ! EULERA computes the Euler number En.
  public :: eulerb          ! EULERB computes the Euler number En.
  public :: fcoef           ! FCOEF: expansion coefficients for Mathieu and modified Mathieu functions.
  public :: fcs             ! FCS computes Fresnel integrals C(x) and S(x).
  public :: fcszo           ! FCSZO computes complex zeros of Fresnel integrals C(x) or S(x).
  public :: ffk             ! FFK computes modified Fresnel integrals F+/-(x) and K+/-(x).
  public :: gaih            ! GAIH computes the GammaH function.
  public :: gam0            ! GAM0 computes the Gamma function for the LAMV function.
  public :: gammaf          ! GAMMA evaluates the Gamma function.
  public :: gmn             ! GMN computes quantities for oblate radial functions with small argument.
  public :: herzo           ! HERZO computes the zeros the Hermite polynomial Hn(x).
  public :: hygfx           ! HYGFX evaluates the hypergeometric function F(A,B,C,X).
  public :: hygfz           ! HYGFZ computes the hypergeometric function F(a,b,c,x) for complex argument.
  public :: ik01a           ! IK01A compute Bessel function I0(x), I1(x), K0(x), and K1(x).
  public :: ik01b           ! IK01B: Bessel functions I0(x), I1(x), K0(x), and K1(x) and derivatives.
  public :: ikna            ! IKNA compute Bessel function In(x) and Kn(x), and derivatives.
  public :: iknb            ! IKNB compute Bessel function In(x) and Kn(x).
  public :: ikv             ! IKV compute modified Bessel function Iv(x) and Kv(x) and their derivatives.
  public :: incob           ! INCOB computes the incomplete beta function Ix(a,b).
  public :: incog           ! INCOG computes the incomplete gamma function r(a,x), ,(a,x), P(a,x).
  public :: itairy          ! ITAIRY computes the integrals of Airy functions.
  public :: itika           ! ITIKA computes the integral of the modified Bessel functions I0(t) and K0(t).
  public :: itikb           ! ITIKB computes the integral of the Bessel functions I0(t) and K0(t).
  public :: itjya           ! ITJYA computes integrals of Bessel functions J0(t) and Y0(t).
  public :: itjyb           ! ITJYB computes integrals of Bessel functions J0(t) and Y0(t).
  public :: itsh0           ! ITSH0 integrates the Struve function H0(t) from 0 to x.
  public :: itsl0           ! ITSL0 integrates the Struve function L0(t) from 0 to x.
  public :: itth0           ! ITTH0 integrates H0(t)/t from x to oo.
  public :: ittika          ! ITTIKA integrates (I0(t)-1)/t from 0 to x, K0(t)/t from x to infinity.
  public :: ittikb          ! ITTIKB integrates (I0(t)-1)/t from 0 to x, K0(t)/t from x to infinity.
  public :: ittjya          ! ITTJYA integrates (1-J0(t))/t from 0 to x, and Y0(t)/t from x to infinity.
  public :: ittjyb          ! ITTJYB integrates (1-J0(t))/t from 0 to x, and Y0(t)/t from x to infinity.
  public :: jdzo            ! JDZO computes the zeros of Bessel functions Jn(x) and Jn'(x).
  public :: jelp            ! JELP computes Jacobian elliptic functions SN(u), CN(u), DN(u).
  public :: jy01a           ! JY01A computes Bessel functions J0(x), J1(x), Y0(x), Y1(x) and derivatives.
  public :: jy01b           ! JY01B computes Bessel functions J0(x), J1(x), Y0(x), Y1(x) and derivatives.
  public :: jyna            ! JYNA computes Bessel functions Jn(x) and Yn(x) and derivatives.
  public :: jynb            ! JYNB computes Bessel functions Jn(x) and Yn(x) and derivatives.
  public :: jyndd           ! JYNDD: Bessel functions Jn(x) and Yn(x), first and second derivatives.
  public :: jyv             ! JYV computes Bessel functions Jv(x) and Yv(x) and their derivatives.
  public :: jyzo            ! JYZO computes the zeros of Bessel functions Jn(x), Yn(x) and derivatives.
  public :: klvna           ! KLVNA: Kelvin functions ber(x), bei(x), ker(x), and kei(x), and derivatives.
  public :: klvnb           ! KLVNB: Kelvin functions ber(x), bei(x), ker(x), and kei(x), and derivatives.
  public :: klvnzo          ! KLVNZO computes zeros of the Kelvin functions.
  public :: kmn             ! KMN: expansion coefficients of prolate or oblate spheroidal functions.
  public :: lagzo           ! LAGZO computes zeros of the Laguerre polynomial, and integration weights.
  public :: lamn            ! LAMN computes lambda functions and derivatives.
  public :: lamv            ! LAMV computes lambda functions and derivatives of arbitrary order.
  public :: legzo           ! LEGZO computes the zeros of Legendre polynomials, and integration weights.
  public :: lgama           ! LGAMA computes the gamma function or its logarithm.
  public :: lpmn            ! LPMN computes associated Legendre functions Pmn(X) and derivatives P'mn(x).
  public :: lpmns           ! LPMNS computes associated Legendre functions Pmn(X) and derivatives P'mn(x).
  public :: lpmv            ! LPMV computes associated Legendre functions Pmv(X) with arbitrary degree.
  public :: lpn             ! LPN computes Legendre polynomials Pn(x) and derivatives Pn'(x).
  public :: lpni            ! LPNI computes Legendre polynomials Pn(x), derivatives, and integrals.
  public :: lqmn            ! LQMN computes associated Legendre functions Qmn(x) and derivatives.
  public :: lqmns           ! LQMNS computes associated Legendre functions Qmn(x) and derivatives Qmn'(x).
  public :: lqna            ! LQNA computes Legendre function Qn(x) and derivatives Qn'(x).
  public :: lqnb            ! LQNB computes Legendre function Qn(x) and derivatives Qn'(x).
  public :: mtu0            ! MTU0 computes Mathieu functions CEM(x,q) and SEM(x,q) and derivatives.
  public :: mtu12           ! MTU12 computes modified Mathieu functions of the first and second kind.
  public :: othpl           ! OTHPL computes orthogonal polynomials Tn(x), Un(x), Ln(x) or Hn(x).
  public :: pbdv            ! PBDV computes parabolic cylinder functions Dv(x) and derivatives.
  public :: pbvv            ! PBVV computes parabolic cylinder functions Vv(x) and their derivatives.
  public :: pbwa            ! PBWA computes parabolic cylinder functions W(a,x) and derivatives.
  public :: psi             ! PSI computes the PSI function.
  public :: qstar           ! QSTAR computes Q*mn(-ic) for oblate radial functions with a small argument.
  public :: rctj            ! RCTJ computes Riccati-Bessel function of the first kind, and derivatives.
  public :: rcty            ! RCTY computes Riccati-Bessel function of the second kind, and derivatives.
  public :: refine          ! REFINE refines an estimate of the characteristic value of Mathieu functions.
  public :: rmn1            ! RMN1 computes prolate and oblate spheroidal functions of the first kind.
  public :: rmn2l           ! RMN2L: prolate and oblate spheroidal functions, second kind, large CX.
  public :: rmn2so          ! RMN2SO: oblate radial functions of the second kind with small argument.
  public :: rmn2sp          ! RMN2SP: prolate, oblate spheroidal radial functions, kind 2, small argument.
  public :: rswfo           ! RSWFO computes prolate spheroidal radial function of first and second kinds.
  public :: rswfp           ! RSWFP computes prolate spheroidal radial function of first and second kinds.
  public :: scka            ! SCKA: expansion coefficients for prolate and oblate spheroidal functions.
  public :: sckb            ! SCKB: expansion coefficients for prolate and oblate spheroidal functions.
  public :: sdmn            ! SDMN: expansion coefficients for prolate and oblate spheroidal functions.
  public :: segv            ! SEGV computes the characteristic values of spheroidal wave functions.
  public :: sphi            ! SPHI computes spherical Bessel functions in(x) and their derivatives in'(x).
  public :: sphj            ! SPHJ computes spherical Bessel functions jn(x) and their derivatives.
  public :: sphk            ! SPHK computes modified spherical Bessel functions kn(x) and derivatives.
  public :: sphy            ! SPHY computes spherical Bessel functions yn(x) and their derivatives.
  public :: stvh0           ! STVH0 computes the Struve function H0(x).
  public :: stvh1           ! STVH1 computes the Struve function H1(x).
  public :: stvhv           ! STVHV computes the Struve function Hv(x) with arbitrary order v.
  public :: stvl0           ! STVL0 computes the modified Struve function L0(x).
  public :: stvl1           ! STVL1 computes the modified Struve function L1(x).
  public :: stvlv           ! STVLV computes the modified Struve function Lv(x) with arbitary order.
  public :: vvla            ! VVLA computes parabolic cylinder function Vv(x) for large arguments.
  public :: vvsa            ! VVSA computes parabolic cylinder function V(nu,x) for small arguments.  
  public :: msta1           ! MSTA1 determines a backward recurrence starting point for Jn(x).
  public :: msta2           ! MSTA2 determines a backward recurrence starting point for Jn(x).
  public :: envj            ! ENVJ is a utility function used by MSTA1 and MSTA2.


contains


  !+------------------------------------------------------------------+
  !PURPOSE  : calculate the Heaviside  function
  !+------------------------------------------------------------------+
  elemental function heaviside(x)
    real(8),intent(in) :: x
    real(8)            :: heaviside
    if(x < 0.d0) then
       heaviside = 0.0d0
    elseif(x==0.d0)then
       heaviside = 0.50d0
    else
       heaviside = 1.0d0
    endif
  end function heaviside


  !+------------------------------------------------------------------+
  !PURPOSE  : calculate step function
  !+------------------------------------------------------------------+
  pure function step_x(x,origin) result(step)
    real(8),intent(in)          :: x
    logical,optional,intent(in) :: origin
    real(8)                     :: step
    logical                     :: w0
    step=0.d0
    w0=.true.;if(present(origin))w0=origin
    select case(w0)
    case (.true.)
       if(x>=0.d0)step=1.d0
    case (.false.)
       if(x>0.d0)step=1.d0
    end select
  end function step_x


  !+------------------------------------------------------------------+
  !PURPOSE  : calculate step function
  !+------------------------------------------------------------------+
  pure function step_ij(i,j,origin)
    integer,intent(in)          :: i,j
    logical,optional,intent(in) :: origin
    real(8)                     :: step_ij
    logical                     :: w0
    step_ij=0.d0
    w0=.true.;if(present(origin))w0=origin
    select case(w0)
    case (.true.)
       if(i.ge.j)step_ij=1.d0
    case (.false.)
       if(i.gt.j)step_ij=1.d0
    end select
  end function step_ij

  !*******************************************************************
  !*******************************************************************
  !*******************************************************************



  !+-------------------------------------------------------------------+
  !PURPOSE  : calculate the Fermi-Dirac distribution
  !+-------------------------------------------------------------------+
  elemental function fermi(x,beta)
    real(8),intent(in) :: x, beta 
    real(8)            :: fermi
    if(x*beta > 100.d0)then
       fermi=0.d0
       return
    endif
    fermi = 1.d0/(1.d0+exp(beta*x))
  end function fermi


  !*******************************************************************
  !*******************************************************************
  !*******************************************************************



  !+-------------------------------------------------------------------+
  !PURPOSE:  evaluate the sign of a given number (I,R)
  !+-------------------------------------------------------------------+
  pure function i_sgn(x) result(sgn)
    integer,intent(in) :: x
    integer            :: sgn
    sgn=x/abs(x)
  end function i_sgn
  pure function d_sgn(x) result(sgn)
    real(8),intent(in) :: x
    real(8)            :: sgn
    sgn=x/abs(x)
  end function d_sgn


  !*******************************************************************
  !*******************************************************************
  !*******************************************************************





  !+------------------------------------------------------------------+
  !PURPOSE  : Evaluate the Complex Error Functions (Faddeeva function)
  ! w(x)=exp(-x^2)erfc(-ix)
  !+------------------------------------------------------------------+
  function wfun(z)
    complex(8):: z,wfun
    real(8)   :: x,y,u,v
    logical   :: flag
    x=real(z,8)
    y=aimag(z)
    call wofz(x,y,u,v,flag)
    wfun=cmplx(u,v)
  contains
    include "functions_wofz.f90"
  end function wfun


  !*******************************************************************
  !*******************************************************************
  !*******************************************************************

  !Double precision complex argument Error function
  include "functions_zerf.f90"


  !###################################################################
  ! BETHE:
  !###################################################################
  include "functions_bethe.f90"



  !+-------------------------------------------------------------------+
  !PURPOSE  : calculate the non-interacting dos for HYPERCUBIC lattice 
  !+-------------------------------------------------------------------+
  pure function dens_hyperc(x,t1)
    real(8),optional,intent(in) :: t1
    real(8),intent(in)          :: x
    REAL(8):: dens_hyperc,t1_,pi2,sqrt2
    pi2=2.d0*acos(-1.d0)
    sqrt2=sqrt(2.d0)
    t1_=sqrt2 ; if(present(t1))t1_=t1
    dens_hyperc = (1/(t1_*sqrt(pi2)))*exp(-(x**2)/(2.d0*t1_**2))
    return
  end function dens_hyperc


  !+-------------------------------------------------------------------+
  !PURPOSE  : calculate the non-interacting dos for HYPERCUBIC lattice 
  !+-------------------------------------------------------------------+
  function dens_2dsquare(x,ts) result(dos)
    real(8),intent(in)          :: x
    real(8),intent(in),optional :: ts
    real(8)                     :: wband,y,kint,eint,dos,pi
    pi=acos(-1.d0)
    wband=4.d0;if(present(ts))wband=4.d0*ts
    dos=0.d0
    if(abs(x)<=wband)then
       y=0.5d0*(x/wband)**2-1.d0
       call comelp(y,kint,eint)
       dos=2.d0/wband/pi**2*kint!*heaviside(wband-abs(y))
    endif
  end function dens_2dsquare


  include "special_functions.f90"

END MODULE SF_SPECIAL
