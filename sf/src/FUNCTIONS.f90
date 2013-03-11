  !###############################################################
  ! PROGRAM  : FUNCTIONS
  ! PURPOSE  : give access to standard functions not implicit in F
  !###############################################################  
  ! AIRYA computes Airy functions and their derivatives.
  ! AIRYB computes Airy functions and their derivatives.
  ! AIRYZO computes the first NT zeros of Ai(x) and Ai'(x).
  ! AJYIK computes Bessel functions Jv(x), Yv(x), Iv(x), Kv(x).
  ! ASWFA: prolate and oblate spheroidal angular functions of the first kind.
  ! ASWFB: prolate and oblate spheroidal angular functions of the first kind.
  ! BERNOA computes the Bernoulli number Bn.
  ! BERNOB computes the Bernoulli number Bn.
  ! BETA computes the Beta function B(p,q).
  ! BJNDD computes Bessel functions Jn(x) and first and second derivatives.
  ! CBK computes coefficients for oblate radial functions with small argument.
  ! CCHG computes the confluent hypergeometric function.
  ! CERF computes the error function and derivative for a complex argument.
  ! CERROR computes the error function for a complex argument.
  ! CERZO evaluates the complex zeros of the error function.
  ! CFC computes the complex Fresnel integral C(z) and C'(z).
  ! CFS computes the complex Fresnel integral S(z) and S'(z).
  ! CGAMA computes the Gamma function for complex argument.
  ! CH12N computes Hankel functions of first and second kinds, complex argument.
  ! CHGM computes the confluent hypergeometric function M(a,b,x).
  ! CHGU computes the confluent hypergeometric function U(a,b,x).
  ! CHGUBI: confluent hypergeometric function with integer argument B.
  ! CHGUIT computes the hypergeometric function using Gauss-Legendre integration.
  ! CHGUL: confluent hypergeometric function U(a,b,x) for large argument X.
  ! CHGUS: confluent hypergeometric function U(a,b,x) for small argument X.
  ! CIK01: modified Bessel I0(z), I1(z), K0(z) and K1(z) for complex argument.
  ! CIKLV: modified Bessel functions Iv(z), Kv(z), complex argument, large order.
  ! CIKNA: modified Bessel functions In(z), Kn(z), derivatives, complex argument.
  ! CIKNB computes complex modified Bessel functions In(z) and Kn(z).
  ! CIKVA: modified Bessel functions Iv(z), Kv(z), arbitrary order, complex.
  ! CIKVB: modified Bessel functions,Iv(z), Kv(z), arbitrary order, complex.
  ! CISIA computes cosine Ci(x) and sine integrals Si(x).
  ! CISIB computes cosine and sine integrals.
  ! CJK: asymptotic expansion coefficients for Bessel functions of large order.
  ! CJY01: complexBessel functions, derivatives, J0(z), J1(z), Y0(z), Y1(z).
  ! CJYLV: Bessel functions Jv(z), Yv(z) of complex argument and large order v.
  ! CJYNA: Bessel functions and derivatives, Jn(z) and Yn(z) of complex argument.
  ! CJYNB: Bessel functions, derivatives, Jn(z) and Yn(z) of complex argument.
  ! CJYVA: Bessel functions and derivatives, Jv(z) and Yv(z) of complex argument.
  ! CJYVB: Bessel functions and derivatives, Jv(z) and Yv(z) of complex argument.
  ! CLPMN: associated Legendre functions and derivatives for complex argument.
  ! CLPN computes Legendre functions and derivatives for complex argument.
  ! CLQMN: associated Legendre functions and derivatives for complex argument.
  ! CLQN: Legendre function Qn(z) and derivative Wn'(z) for complex argument.
  ! COMELP computes complete elliptic integrals K(k) and E(k).
  ! CPBDN: parabolic cylinder function Dn(z) and Dn'(z) for complex argument.
  ! CPDLA computes complex parabolic cylinder function Dn(z) for large argument.
  ! CPDSA computes complex parabolic cylinder function Dn(z) for small argument.
  ! CPSI computes the psi function for a complex argument.
  ! CSPHIK: complex modified spherical Bessel functions and derivatives.
  ! CSPHJY: spherical Bessel functions jn(z) and yn(z) for complex argument.
  ! CV0 computes the initial characteristic value of Mathieu functions.
  ! CVA1 computes a sequence of characteristic values of Mathieu functions.
  ! CVA2 computes a specific characteristic value of Mathieu functions.
  ! CVF computes F for the characteristic equation of Mathieu functions.
  ! CVQL computes the characteristic value of Mathieu functions for q <= 3*m.
  ! CVQM computes the characteristic value of Mathieu functions for q <= m*m.
  ! CY01 computes complex Bessel functions Y0(z) and Y1(z) and derivatives.
  ! CYZO computes zeros of complex Bessel functions Y0(z) and Y1(z) and Y1'(z).
  ! DVLA computes parabolic cylinder functions Dv(x) for large argument.
  ! DVSA computes parabolic cylinder functions Dv(x) for small argument.
  ! E1XA computes the exponential integral E1(x).
  ! E1XB computes the exponential integral E1(x).
  ! E1Z computes the complex exponential integral E1(z).
  ! EIX computes the exponential integral Ei(x).
  ! ELIT: complete and incomplete elliptic integrals F(k,phi) and E(k,phi).
  ! ELIT3 computes the elliptic integral of the third kind.
  ! ENXA computes the exponential integral En(x).
  ! ENXB computes the exponential integral En(x).
  ! ERROR evaluates the error function.
  ! EULERA computes the Euler number En.
  ! EULERB computes the Euler number En.
  ! FCOEF: expansion coefficients for Mathieu and modified Mathieu functions.
  ! FCS computes Fresnel integrals C(x) and S(x).
  ! FCSZO computes complex zeros of Fresnel integrals C(x) or S(x).
  ! FFK computes modified Fresnel integrals F+/-(x) and K+/-(x).
  ! GAIH computes the GammaH function.
  ! GAM0 computes the Gamma function for the LAMV function.
  ! GAMMA evaluates the Gamma function.
  ! GMN computes quantities for oblate radial functions with small argument.
  ! HERZO computes the zeros the Hermite polynomial Hn(x).
  ! HYGFX evaluates the hypergeometric function F(A,B,C,X).
  ! HYGFZ computes the hypergeometric function F(a,b,c,x) for complex argument.
  ! IK01A compute Bessel function I0(x), I1(x), K0(x), and K1(x).
  ! IK01B: Bessel functions I0(x), I1(x), K0(x), and K1(x) and derivatives.
  ! IKNA compute Bessel function In(x) and Kn(x), and derivatives.
  ! IKNB compute Bessel function In(x) and Kn(x).
  ! IKV compute modified Bessel function Iv(x) and Kv(x) and their derivatives.
  ! INCOB computes the incomplete beta function Ix(a,b).
  ! INCOG computes the incomplete gamma function r(a,x), ,(a,x), P(a,x).
  ! ITAIRY computes the integrals of Airy functions.
  ! ITIKA computes the integral of the modified Bessel functions I0(t) and K0(t).
  ! ITIKB computes the integral of the Bessel functions I0(t) and K0(t).
  ! ITJYA computes integrals of Bessel functions J0(t) and Y0(t).
  ! ITJYB computes integrals of Bessel functions J0(t) and Y0(t).
  ! ITSH0 integrates the Struve function H0(t) from 0 to x.
  ! ITSL0 integrates the Struve function L0(t) from 0 to x.
  ! ITTH0 integrates H0(t)/t from x to oo.
  ! ITTIKA integrates (I0(t)-1)/t from 0 to x, K0(t)/t from x to infinity.
  ! ITTIKB integrates (I0(t)-1)/t from 0 to x, K0(t)/t from x to infinity.
  ! ITTJYA integrates (1-J0(t))/t from 0 to x, and Y0(t)/t from x to infinity.
  ! ITTJYB integrates (1-J0(t))/t from 0 to x, and Y0(t)/t from x to infinity.
  ! JDZO computes the zeros of Bessel functions Jn(x) and Jn'(x).
  ! JELP computes Jacobian elliptic functions SN(u), CN(u), DN(u).
  ! JY01A computes Bessel functions J0(x), J1(x), Y0(x), Y1(x) and derivatives.
  ! JY01B computes Bessel functions J0(x), J1(x), Y0(x), Y1(x) and derivatives.
  ! JYNA computes Bessel functions Jn(x) and Yn(x) and derivatives.
  ! JYNB computes Bessel functions Jn(x) and Yn(x) and derivatives.
  ! JYNDD: Bessel functions Jn(x) and Yn(x), first and second derivatives.
  ! JYV computes Bessel functions Jv(x) and Yv(x) and their derivatives.
  ! JYZO computes the zeros of Bessel functions Jn(x), Yn(x) and derivatives.
  ! KLVNA: Kelvin functions ber(x), bei(x), ker(x), and kei(x), and derivatives.
  ! KLVNB: Kelvin functions ber(x), bei(x), ker(x), and kei(x), and derivatives.
  ! KLVNZO computes zeros of the Kelvin functions.
  ! KMN: expansion coefficients of prolate or oblate spheroidal functions.
  ! LAGZO computes zeros of the Laguerre polynomial, and integration weights.
  ! LAMN computes lambda functions and derivatives.
  ! LAMV computes lambda functions and derivatives of arbitrary order.
  ! LEGZO computes the zeros of Legendre polynomials, and integration weights.
  ! LGAMA computes the gamma function or its logarithm.
  ! LPMN computes associated Legendre functions Pmn(X) and derivatives P'mn(x).
  ! LPMNS computes associated Legendre functions Pmn(X) and derivatives P'mn(x).
  ! LPMV computes associated Legendre functions Pmv(X) with arbitrary degree.
  ! LPN computes Legendre polynomials Pn(x) and derivatives Pn'(x).
  ! LPNI computes Legendre polynomials Pn(x), derivatives, and integrals.
  ! LQMN computes associated Legendre functions Qmn(x) and derivatives.
  ! LQMNS computes associated Legendre functions Qmn(x) and derivatives Qmn'(x).
  ! LQNA computes Legendre function Qn(x) and derivatives Qn'(x).
  ! LQNB computes Legendre function Qn(x) and derivatives Qn'(x).
  ! MTU0 computes Mathieu functions CEM(x,q) and SEM(x,q) and derivatives.
  ! MTU12 computes modified Mathieu functions of the first and second kind.
  ! OTHPL computes orthogonal polynomials Tn(x), Un(x), Ln(x) or Hn(x).
  ! PBDV computes parabolic cylinder functions Dv(x) and derivatives.
  ! PBVV computes parabolic cylinder functions Vv(x) and their derivatives.
  ! PBWA computes parabolic cylinder functions W(a,x) and derivatives.
  ! PSI computes the PSI function.
  ! QSTAR computes Q*mn(-ic) for oblate radial functions with a small argument.
  ! RCTJ computes Riccati-Bessel function of the first kind, and derivatives.
  ! RCTY computes Riccati-Bessel function of the second kind, and derivatives.
  ! REFINE refines an estimate of the characteristic value of Mathieu functions.
  ! RMN1 computes prolate and oblate spheroidal functions of the first kind.
  ! RMN2L: prolate and oblate spheroidal functions, second kind, large CX.
  ! RMN2SO: oblate radial functions of the second kind with small argument.
  ! RMN2SP: prolate, oblate spheroidal radial functions, kind 2, small argument.
  ! RSWFO computes prolate spheroidal radial function of first and second kinds.
  ! RSWFP computes prolate spheroidal radial function of first and second kinds.
  ! SCKA: expansion coefficients for prolate and oblate spheroidal functions.
  ! SCKB: expansion coefficients for prolate and oblate spheroidal functions.
  ! SDMN: expansion coefficients for prolate and oblate spheroidal functions.
  ! SEGV computes the characteristic values of spheroidal wave functions.
  ! SPHI computes spherical Bessel functions in(x) and their derivatives in'(x).
  ! SPHJ computes spherical Bessel functions jn(x) and their derivatives.
  ! SPHK computes modified spherical Bessel functions kn(x) and derivatives.
  ! SPHY computes spherical Bessel functions yn(x) and their derivatives.
  ! STVH0 computes the Struve function H0(x).
  ! STVH1 computes the Struve function H1(x).
  ! STVHV computes the Struve function Hv(x) with arbitrary order v.
  ! STVL0 computes the modified Struve function L0(x).
  ! STVL1 computes the modified Struve function L1(x).
  ! STVLV computes the modified Struve function Lv(x) with arbitary order.
  ! TIMESTAMP prints the current YMDHMS date as a time stamp.
  ! VVLA computes parabolic cylinder function Vv(x) for large arguments.
  ! VVSA computes parabolic cylinder function V(nu,x) for small arguments.  
  ! MSTA1 determines a backward recurrence starting point for Jn(x).
  ! MSTA2 determines a backward recurrence starting point for Jn(x).
  ! ENVJ is a utility function used by MSTA1 and MSTA2.
  include "functions_special_funcs.f90"
  !
  module FUNCTIONS
    USE COMMON_VARS
    implicit none
    private

    !FUNCTIONS:
    public :: heaviside
    public :: step
    public :: fermi
    interface sgn
       module procedure i_sgn,d_sgn
    end interface sgn
    public :: sgn

    !ERROR FUNCS
    public :: wfun         !complex error function (Faddeeva function)
    public :: zerf   


  contains


    !+------------------------------------------------------------------+
    !PURPOSE  : calculate the Heaviside  function
    !+------------------------------------------------------------------+
    pure function heaviside(x)
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
    pure function step(x,origin)
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
    end function step



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

  END module FUNCTIONS
