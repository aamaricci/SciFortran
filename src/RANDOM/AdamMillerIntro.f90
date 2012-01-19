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
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
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
  
  !     Author: Alan Miller
  !             CSIRO Division of Mathematical & Information Sciences
  !             Private Bag 10, Clayton South MDC
  !             Clayton 3169, Victoria, Australia
  !     Phone: (+61) 3 9545-8016      Fax: (+61) 3 9545-8080
  !     e-mail: amiller @ bigpond.net.au
  
