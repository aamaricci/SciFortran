!###############################################################
!PROGRAM  : RANDOM
!TYPE     : Module
!PURPOSE  : Module for Random Number generators
!###############################################################
module RANDOM
  implicit none
  private
  integer           :: i,j,k,D
  real(8),parameter :: pi=3.14159265358979d0
  real(8),parameter :: pi2=6.28318530717959d0
  real(8),parameter :: sqrt2 = 1.41421356237309504880169d0
  real(8),parameter :: sqrt3 = 1.73205080756887729352745d0
  real(8),parameter :: sqrt6 = 2.44948974278317809819728d0

  !adam Miller parameters:
  integer,parameter :: dp = SELECTED_REAL_KIND(12, 60)
  real    :: zero = 0.0, half = 0.5, one = 1.0, two = 2.0, vsmall = TINY(1.0), vlarge = HUGE(1.0)

  public :: nrand
  public :: init_random_number
  public :: random_order
  !
  public :: random_normal
  public :: random_gamma
  public :: random_gamma1
  public :: random_gamma2
  public :: random_chisq
  public :: random_exponential
  public :: random_Weibull
  public :: random_beta
  public :: random_inv_gauss
  public :: random_Poisson
  public :: random_binomial1
  public :: bin_prob
  public :: lngamma
  public :: random_binomial2
  public :: random_neg_binomial
  public :: random_von_Mises
  public :: random_Cauchy


contains


  !+-----------------------------------------------------------------+
  !PURPOSE  :   
  !+-----------------------------------------------------------------+
  real(8) function nrand(dseed)
    integer                 :: dseed
    integer,parameter       :: IM1=2147483563, IM2=2147483399, IMM1=IM1-1, IA1=40014, &
         & IA2=40692, IQ1=53668, IQ2=52774, IR1=12211, IR2=3791,  &
         & NTAB=32, NDIV=1+IMM1/NTAB
    real(kind=8), parameter :: AM=1.0d0/IM1, EPS=1.2e-7, RNMX=1.-EPS
    integer                 :: dseed2, j, k, iv(NTAB), iy
    save iv, iy, dseed2
    data dseed2/123456789/, iv/NTAB*0/, iy/0/
    if(dseed .le. 0) then
       dseed = max(-dseed,1)
       dseed2 = dseed
       do j=NTAB+8, 1, -1
          k = dseed/IQ1
          dseed = IA1*(dseed-k*IQ1)-k*IR1
          if(dseed .lt. 0) dseed = dseed+IM1
          if(j .le. NTAB) iv(j) = dseed
       enddo
       iy=iv(1)
    endif
    k = dseed/IQ1
    dseed = IA1*(dseed-k*IQ1)-k*IR1
    if(dseed .lt. 0) dseed = dseed+IM1
    k = dseed2/IQ2
    dseed2 = IA2*(dseed2-k*IQ2)-k*IR2
    if(dseed2 .lt. 0) dseed2 = dseed2+IM2
    j = 1+iy/NDIV
    iy = iv(j)-dseed2
    iv(j) = dseed
    if(iy .lt. 1) iy = iy+IMM1
    nrand = min(AM*iy,RNMX)
  end function nrand






  !+-----------------------------------------------------------------+
  !PURPOSE  : initialize the seed for the INTRINSIC RNG 
  !+-----------------------------------------------------------------+
  subroutine init_random_number(shift)
    integer,optional                 :: shift
    integer                          :: i, n, clock
    integer,dimension(:),allocatable :: seed
    call RANDOM_SEED(size = n)
    allocate(seed(n))
    call SYSTEM_CLOCK(COUNT=clock)
    seed = clock + 37 * (/ (i - 1, i = 1, n) /)
    if(present(shift))seed=seed+shift
    call RANDOM_SEED(PUT = seed)
    deallocate(seed)
  end subroutine init_random_number





  !+-----------------------------------------------------------------+
  !PURPOSE  :   
  !+-----------------------------------------------------------------+
  SUBROUTINE random_order(order, n)
    !     generate a random ordering of the integers 1 ... n.
    integer, intent(in)  :: n
    integer, intent(out) :: order(n)
    !     local variables
    integer :: i, j, k
    real(8) :: wk
    do i = 1, n
       order(i) = i
    end do
    !     starting at the end, swap the current last indicator with one
    !     randomly chosen from those preceeding it.
    do i = n, 2, -1
       call random_number(wk)
       j = 1 + i * wk
       if (j < i) then
          k = order(i)
          order(i) = order(j)
          order(j) = k
       end if
    end do
    return
  end subroutine random_order



  !+-----------------------------------------------------------------+
  !purpose  : rng library, not exported!! 
  !+-----------------------------------------------------------------+
  ! A module for random number generation from the following distributions:
  !
  !     Distribution                    Function/subroutine name
  !
  !     Normal (Gaussian)               random_normal
  !     Gamma                           random_gamma
  !     Chi-squared                     random_chisq
  !     Exponential                     random_exponential
  !     Weibull                         random_Weibull
  !     Beta                            random_beta
  !     t                               random_t
  !     Multivariate normal             random_mvnorm
  !     Generalized inverse Gaussian    random_inv_gauss
  !     Poisson                         random_Poisson
  !     Binomial                        random_binomial1   *
  !                                     random_binomial2   *
  !     Negative binomial               random_neg_binomial
  !     von Mises                       random_von_Mises
  !     Cauchy                          random_Cauchy
  !
  !  Generate a random ordering of the integers 1 .. N
  !                                     random_order
  !     Initialize (seed) the uniform random number generator for ANY compiler
  !                                     seed_random_number
  !     Lognormal - see note below.
  !  ** Two functions are provided for the binomial distribution.
  !  If the parameter values remain constant, it is recommended that the
  !  first function is used (random_binomial1).   If one or both of the
  !  parameters change, use the second function (random_binomial2).
  ! The compilers own random number generator, SUBROUTINE RANDOM_NUMBER(r),
  ! is used to provide a source of uniformly distributed random numbers.
  ! N.B. At this stage, only one random number is generated at each call to
  !      one of the functions above.
  ! The module uses the following functions which are included here:
  ! bin_prob to calculate a single binomial probability
  ! lngamma  to calculate the logarithm to base e of the gamma function
  ! Some of the code is adapted from Dagpunar's book:
  !     Dagpunar, J. 'Principles of random variate generation'
  !     Clarendon Press, Oxford, 1988.   ISBN 0-19-852202-9
  !
  ! In most of Dagpunar's routines, there is a test to see whether the value
  ! of one or two floating-point parameters has changed since the last call.
  ! These tests have been replaced by using a logical variable FIRST.
  ! This should be set to .TRUE. on the first call using new values of the
  ! parameters, and .FALSE. if the parameter values are the same as for the
  ! previous call.
  ! Lognormal distribution
  ! If X has a lognormal distribution, then log(X) is normally distributed.
  ! Here the logarithm is the natural logarithm, that is to base e, sometimes
  ! denoted as ln.  To generate random variates from this distribution, generate
  ! a random deviate from the normal distribution with mean and variance equal
  ! to the mean and variance of the logarithms of X, then take its exponential.
  ! Relationship between the mean & variance of log(X) and the mean & variance
  ! of X, when X has a lognormal distribution.
  ! Let m = mean of log(X), and s^2 = variance of log(X)
  ! Then
  ! mean of X     = exp(m + 0.5s^2)
  ! variance of X = (mean(X))^2.[exp(s^2) - 1]
  ! In the reverse direction (rarely used)
  ! variance of log(X) = log[1 + var(X)/(mean(X))^2]
  ! mean of log(X)     = log(mean(X) - 0.5var(log(X))
  ! N.B. The above formulae relate to population parameters; they will only be
  !      approximate if applied to sample values.
  ! Version 1.13, 2 October 2000
  ! Changes from version 1.01
  ! 1. The random_order, random_Poisson & random_binomial routines have been
  !    replaced with more efficient routines.
  ! 2. A routine, seed_random_number, has been added to seed the uniform random
  !    number generator.   This requires input of the required number of seeds
  !    for the particular compiler from a specified I/O unit such as a keyboard.
  ! 3. Made compatible with Lahey's ELF90.
  ! 4. Marsaglia & Tsang algorithm used for random_gamma when shape parameter > 1.
  ! 5. INTENT for array f corrected in random_mvnorm.

  include "random_routines.f90"


end module random
