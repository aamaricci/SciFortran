module SF_RANDOM
   use RNG_LEGACY
   implicit none
   private

   !MT parameters:
   integer,parameter :: defaultsd = 4357 !Default seed
   integer,parameter :: N=624, N1=N+1    !Period parameters
   integer,save      :: mt(0:n-1)        !the array for the state vector
   integer,save      :: mti = n1         !mti==N+1 means mt[N] is not initialized

   integer           :: i,j,k,D
   real(8),parameter :: pi    = 3.14159265358979d0
   real(8),parameter :: pi2   = 6.28318530717959d0
   real(8),parameter :: sqrt2 = 1.41421356237309504880169d0
   real(8),parameter :: sqrt3 = 1.73205080756887729352745d0
   real(8),parameter :: sqrt6 = 2.44948974278317809819728d0


   !MT interface:
   interface mersenne_init
      module procedure :: init_genrand
   end interface mersenne_init

   interface mt_init
      module procedure :: init_genrand
   end interface mt_init

   interface mersenne
      module procedure :: grnd
   end interface mersenne

   interface mt_random
      module procedure :: d_grnd_1
      module procedure :: d_grnd_2
      module procedure :: d_grnd_3
      module procedure :: d_grnd_4
      module procedure :: d_grnd_5
      module procedure :: d_grnd_6
      module procedure :: d_grnd_7
      !
      module procedure :: c_grnd_1
      module procedure :: c_grnd_2
      module procedure :: c_grnd_3
      module procedure :: c_grnd_4
      module procedure :: c_grnd_5
      module procedure :: c_grnd_6
      module procedure :: c_grnd_7
   end interface mt_random

   interface mt_uniform
      module procedure :: igrnd
      module procedure :: dgrnd_uniform
   end interface mt_uniform

   interface mt_normal
      module procedure :: gaussrnd
      module procedure :: normalrnd
   end interface mt_normal

   interface mt_exponential
      module procedure :: exponentialrnd
   end interface mt_exponential

   interface mt_gamma
      module procedure :: gammarnd
   end interface mt_gamma

   interface mt_chi_square
      module procedure :: chi_squarernd
   end interface mt_chi_square

   interface mt_inverse_gamma
      module procedure :: inverse_gammarnd
   end interface mt_inverse_gamma

   interface mt_weibull
      module procedure :: weibullrnd
   end interface mt_weibull

   interface mt_cauchy
      module procedure :: cauchyrnd
   end interface mt_cauchy

   interface mt_student_t
      module procedure :: student_trnd
   end interface mt_student_t

   interface mt_laplace
      module procedure :: laplacernd
   end interface mt_laplace

   interface mt_log_normal
      module procedure :: log_normalrnd
   end interface mt_log_normal

   interface mt_beta
      module procedure :: betarnd
   end interface mt_beta

   ! Overload procedures for saving and getting mt state
   interface mt_save
      module procedure :: mtsavef
      module procedure :: mtsaveu
   end interface mt_save

   interface mt_get
      module procedure :: mtgetf
      module procedure :: mtgetu
   end interface mt_get



   public :: mersenne
   public :: mersenne_init
   !
   public :: mt_random
   public :: mt_uniform
   public :: mt_normal
   public :: mt_exponential
   public :: mt_gamma
   public :: mt_chi_square
   public :: mt_inverse_gamma
   public :: mt_weibull
   public :: mt_cauchy
   public :: mt_student_t
   public :: mt_laplace
   public :: mt_log_normal
   public :: mt_beta
   !
   public :: mt_save
   public :: mt_get
   public :: mt_init
   !
   public :: random_number_init
   public :: random_number_seed
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
   !
   public :: nrand
   public :: random_order


contains




   !+-----------------------------------------------------------------+
   !PURPOSE  : INTRINSIC RNG initialization
   !+-----------------------------------------------------------------+
   subroutine random_number_init(shift)
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
   end subroutine random_number_init





   !------------------------------------------------------!
   !  Get the RNG seed from /dev/urandom device.          !
   !                                                      !
   !  In order to get positive seed the most              !
   !  significant bit in the number read from the         !
   !  device is cleared (by anding it with LMASK).        !
   !                                                      !
   !  NOTE: Routine uses the default integer type.        !
   !                                                      !
   !  If the device can not be opened or read routine     !
   !  falls back to calculating seed from the current     !
   !  time.                                               !
   !                                                      !
   !  Note that stream i/o is used which is a Fortran     !
   !  2003 feature.                                       !
   !                                                      !
   !  Input parameters:                                   !
   !    info : integer, if /=0 print info to stdout       !
   !    file : integer, 0: use /dev/urandom               !
   !                    1: use /dev/random                !
   !                                                      !
   !  Both parameters are optional, so that the simplest  !
   !  way to call the function is 'getseed()'.            !
   !                                                      !
   !  Generating a large amount of random numbers using   !
   !  /dev/random may be slow because quality of random   !
   !  bits from this device is guaranteed and system may  !
   !  have to wait while enough 'entropy' is collected    !
   !  from network traffic, keyboard etc.                 !
   !                                                      !
   !  A.Kuronen, antti.kuronen@helsinki.fi, 2008-2014     !
   !------------------------------------------------------!
   integer function random_number_seed(info,file)
      integer,optional,intent(in) :: info,file
      integer                     :: t(8),rn,is
      integer,parameter           :: LMASK=huge(rn) ! = 0111...111
      integer,parameter           :: LUN=676769
      character (len=80)          :: rdev0='/dev/urandom',rdev1='/dev/random',rdev
      logical                     :: openok,readok,printinfo
      openok=.true.
      readok=.true.
      if (present(file)) then
         if (file==0) then
            rdev=rdev0
         else
            rdev=rdev1
         end if
      else
         rdev=rdev0
      end if
      if (present(info)) then
         printinfo=(info/=0)
      else
         printinfo=.false.
      end if
      open(LUN,file=rdev,form='unformatted',access='stream',action='read',iostat=is)
      if (is/=0) then
         openok=.false.
         print *,'open',is
      else
         read(LUN,iostat=is) rn
         if (is/=0) then
            readok=.false.
         end if
      end if
      if (openok) close(LUN)
      if (openok.and.readok) then
         rn=iand(rn,LMASK) ! Make it positive, i.e. zero the leftmost bit
         if (printinfo) write(6,'(a,a,a,i0)') 'Seed from ',trim(rdev),': ',rn
      else
         call date_and_time(values=t)
         rn=t(7)+60*(t(6)+60*(t(5)+24*(t(3)-1+31*(t(2)-1+12*t(1)))))+t(8)
         if (printinfo) write(6,'(a,i12)') 'Seed from time:',rn
      end if
      random_number_seed=rn
      return
   end function random_number_seed






   !+-----------------------------------------------------------------+
   !PURPOSE  : MERSENNE TWISTER RNG
   !+-----------------------------------------------------------------+
   !---------------------------------------------------------------!
! Initialization subroutine                                     !
!---------------------------------------------------------------!
   subroutine sgrnd(seed)
      ! Setting initial seeds to mt[N] using the generator Line 25 of Table 1 in
      ! [KNUTH 1981, The Art of Computer Programming Vol. 2 (2nd Ed.), pp102]
      integer,intent(in) :: seed
      mt(0) = iand(seed,-1)
      do mti=1,N-1
         mt(mti) = iand(69069 * mt(mti-1),-1)
      enddo
      return
   end subroutine sgrnd

   subroutine init_genrand(seed)
      ! This initialization is based upon the multiplier given on p.106 of the
      ! 3rd edition of Knuth, The Art of Computer Programming Vol. 2.
      ! This version assumes that integer overflow does NOT cause a crash.
      integer,intent(in) :: seed !
      integer             :: latest
      mt(0) = seed
      latest = seed
      do mti = 1,n-1
         latest = IEOR( latest,ISHFT( latest,-30 ) )
         latest = latest * 1812433253 + mti
         mt(mti) = latest
      end do
      return
   end subroutine init_genrand





   !---------------------------------------------------------------!
   ! Random number generator: [0,1[                                !
   !---------------------------------------------------------------!
   function grnd()
      real(8)            :: grnd
      integer,parameter :: M = 397, MATA  = -1727483681 ! constant vector a
      integer,parameter :: LMASK =  2147483647          ! least significant r bits
      integer,parameter :: UMASK = -LMASK - 1           ! most significant w-r bits
      integer,parameter :: TMASKB= -1658038656, TMASKC= -272236544
      integer,save       :: mag01(0:1)=[0,MATA] ! mag01(x) = x * MATA for x=0,1
      integer            :: kk,y
      ! These statement functions have been replaced with separate functions
      ! tshftu(y) = ISHFT(y,-11)
      ! tshfts(y) = ISHFT(y,7)
      ! tshftt(y) = ISHFT(y,15)
      ! tshftl(y) = ISHFT(y,-18)
      if (mti>=N) then               ! generate N words at one time
         if (mti==N+1) then          ! if sgrnd() has not been called,
            call sgrnd( defaultsd ) ! a default initial seed is used
         endif
         do kk=0,N-M-1
            y=ior(iand(mt(kk),UMASK),iand(mt(kk+1),LMASK))
            mt(kk)=ieor(ieor(mt(kk+M),ishft(y,-1)),mag01(iand(y,1)))
         enddo
         do kk=N-M,N-2
            y=ior(iand(mt(kk),UMASK),iand(mt(kk+1),LMASK))
            mt(kk)=ieor(ieor(mt(kk+(M-N)),ishft(y,-1)),mag01(iand(y,1)))
         enddo
         y=ior(iand(mt(N-1),UMASK),iand(mt(0),LMASK))
         mt(N-1)=ieor(ieor(mt(M-1),ishft(y,-1)),mag01(iand(y,1)))
         mti = 0
      endif
      y=mt(mti)
      mti = mti + 1
      y=ieor(y,TSHFTU(y))
      y=ieor(y,iand(TSHFTS(y),TMASKB))
      y=ieor(y,iand(TSHFTT(y),TMASKC))
      y=ieor(y,TSHFTL(y))
      if (y<0) then
         grnd=(dble(y)+2d0**32)/(2d0**32-1d0)
      else
         grnd=dble(y)/(2d0**32-1d0)
      endif
      return
   contains
      integer function TSHFTU(y)
         integer,intent(in) :: y
         TSHFTU=ishft(y,-11)
         return
      end function TSHFTU
      integer function TSHFTS(y)
         integer,intent(in) :: y
         TSHFTS=ishft(y,7)
         return
      end function TSHFTS
      integer function TSHFTT(y)
         integer,intent(in) :: y
         TSHFTT=ishft(y,15)
         return
      end function TSHFTT
      integer function TSHFTL(y)
         integer,intent(in) :: y
         TSHFTL=ishft(y,-18)
         return
      end function TSHFTL
   end function grnd





   subroutine d_grnd_1(A)
      real(8),dimension(:) :: A
      integer              :: i
      do i=1,size(A)
         A(i) = mersenne()
      enddo
   end subroutine d_grnd_1

   subroutine d_grnd_2(A)
      real(8),dimension(:,:) :: A
      integer              :: i1,i2
      do i1=1,size(A,1)
         do i2=2,size(A,2)
            A(i1,i2) = mersenne()
         enddo
      enddo
   end subroutine d_grnd_2

   subroutine d_grnd_3(A)
      real(8),dimension(:,:,:) :: A
      integer                  :: i1,i2,i3
      do i1=1,size(A,1)
         do i2=2,size(A,2)
            do i3=2,size(A,3)
               A(i1,i2,i3) = mersenne()
            enddo
         enddo
      enddo
   end subroutine d_grnd_3

   subroutine d_grnd_4(A)
      real(8),dimension(:,:,:,:) :: A
      integer                  :: i1,i2,i3,i4
      do i1=1,size(A,1)
         do i2=2,size(A,2)
            do i3=2,size(A,3)
               do i4=1,size(A,4)
                  A(i1,i2,i3,i4) = mersenne()
               enddo
            enddo
         enddo
      enddo
   end subroutine d_grnd_4

   subroutine d_grnd_5(A)
      real(8),dimension(:,:,:,:,:) :: A
      integer                  :: i1,i2,i3,i4,i5
      do i1=1,size(A,1)
         do i2=2,size(A,2)
            do i3=2,size(A,3)
               do i4=1,size(A,4)
                  do i5=1,size(A,5)
                     A(i1,i2,i3,i4,i5) = mersenne()
                  enddo
               enddo
            enddo
         enddo
      enddo
   end subroutine d_grnd_5

   subroutine d_grnd_6(A)
      real(8),dimension(:,:,:,:,:,:) :: A
      integer                  :: i1,i2,i3,i4,i5,i6
      do i1=1,size(A,1)
         do i2=2,size(A,2)
            do i3=2,size(A,3)
               do i4=1,size(A,4)
                  do i5=1,size(A,5)
                     do i6=1,size(A,6)
                        A(i1,i2,i3,i4,i5,i6) = mersenne()
                     enddo
                  enddo
               enddo
            enddo
         enddo
      enddo
   end subroutine d_grnd_6


   subroutine d_grnd_7(A)
      real(8),dimension(:,:,:,:,:,:,:) :: A
      integer                  :: i1,i2,i3,i4,i5,i6,i7
      do i1=1,size(A,1)
         do i2=2,size(A,2)
            do i3=2,size(A,3)
               do i4=1,size(A,4)
                  do i5=1,size(A,5)
                     do i6=1,size(A,6)
                        do i7=1,size(A,7)
                           A(i1,i2,i3,i4,i5,i6,i7) = mersenne()
                        enddo
                     enddo
                  enddo
               enddo
            enddo
         enddo
      enddo
   end subroutine d_grnd_7










   subroutine c_grnd_1(A)
      complex(8),dimension(:) :: A
      integer                 :: i
      do i=1,size(A)
         A(i) = dcmplx(mersenne(),mersenne())
      enddo
   end subroutine c_grnd_1

   subroutine c_grnd_2(A)
      complex(8),dimension(:,:) :: A
      integer                   :: i1,i2
      do i1=1,size(A,1)
         do i2=2,size(A,2)
            A(i1,i2) = dcmplx(mersenne(),mersenne())
         enddo
      enddo
   end subroutine c_grnd_2

   subroutine c_grnd_3(A)
      complex(8),dimension(:,:,:) :: A
      integer                     :: i1,i2,i3
      do i1=1,size(A,1)
         do i2=2,size(A,2)
            do i3=2,size(A,3)
               A(i1,i2,i3) = dcmplx(mersenne(),mersenne())
            enddo
         enddo
      enddo
   end subroutine c_grnd_3

   subroutine c_grnd_4(A)
      complex(8),dimension(:,:,:,:) :: A
      integer                       :: i1,i2,i3,i4
      do i1=1,size(A,1)
         do i2=2,size(A,2)
            do i3=2,size(A,3)
               do i4=1,size(A,4)
                  A(i1,i2,i3,i4) = dcmplx(mersenne(),mersenne())
               enddo
            enddo
         enddo
      enddo
   end subroutine c_grnd_4

   subroutine c_grnd_5(A)
      complex(8),dimension(:,:,:,:,:) :: A
      integer                         :: i1,i2,i3,i4,i5
      do i1=1,size(A,1)
         do i2=2,size(A,2)
            do i3=2,size(A,3)
               do i4=1,size(A,4)
                  do i5=1,size(A,5)
                     A(i1,i2,i3,i4,i5) = dcmplx(mersenne(),mersenne())
                  enddo
               enddo
            enddo
         enddo
      enddo
   end subroutine c_grnd_5

   subroutine c_grnd_6(A)
      complex(8),dimension(:,:,:,:,:,:) :: A
      integer                           :: i1,i2,i3,i4,i5,i6
      do i1=1,size(A,1)
         do i2=2,size(A,2)
            do i3=2,size(A,3)
               do i4=1,size(A,4)
                  do i5=1,size(A,5)
                     do i6=1,size(A,6)
                        A(i1,i2,i3,i4,i5,i6) = dcmplx(mersenne(),mersenne())
                     enddo
                  enddo
               enddo
            enddo
         enddo
      enddo
   end subroutine c_grnd_6


   subroutine c_grnd_7(A)
      complex(8),dimension(:,:,:,:,:,:,:) :: A
      integer                             :: i1,i2,i3,i4,i5,i6,i7
      do i1=1,size(A,1)
         do i2=2,size(A,2)
            do i3=2,size(A,3)
               do i4=1,size(A,4)
                  do i5=1,size(A,5)
                     do i6=1,size(A,6)
                        do i7=1,size(A,7)
                           A(i1,i2,i3,i4,i5,i6,i7) = dcmplx(mersenne(),mersenne())
                        enddo
                     enddo
                  enddo
               enddo
            enddo
         enddo
      enddo
   end subroutine c_grnd_7








   !---------------------------------------------------------------!
   ! Integer random number generator:                              !
   ! return a random integer between in [l,h]                      !
   ! (boundaries l and h included)                                 !
   !                                                               !
   ! A.Kuronen, 2009                                               !
   !---------------------------------------------------------------!
   integer function igrnd(l,h)
      integer,intent(in) :: l,h
      real(8)            :: u,r
      u=grnd()
      r=(h-l+1)*u+l
      igrnd=int(r)
      return
   end function igrnd


   function dgrnd_uniform(a,b) result(c)
      real(8) :: a,b,c
      c = a + (b-a)*mersenne()
   end function dgrnd_uniform






   !---------------------------------------------------------------!
   ! Random numbers with normal (Gaussian) distribution.           !
   ! Mean is 0 and standard deviation is 1                         !
   ! See W.H.Press et al., Numerical Recipes 1st ed., page 203     !
   !                                                               !
   ! A.Kuronen, 2009                                               !
   !---------------------------------------------------------------!
   function gaussrnd()
      real(8)       :: gaussrnd
      real(8)       :: fac,v1,v2,r
      real(8),save :: gset
      integer,save :: iset=0
      if (iset==0) then ! Create a new RN
         r=100.0
         do while (r>1.0)
            v1 = 2.0*grnd()-1.0
            v2 = 2.0*grnd()-1.0
            r = v1*v1+v2*v2
         end do
         fac = sqrt(-2.0*log(r)/r)
         gset = v1*fac
         gaussrnd = v2*fac
         iset = 1
      else ! Use the 2nd NR from the previous call
         gaussrnd = gset
         iset = 0
      endif
      return
   end function gaussrnd


   ! Random Sample from normal (Gaussian) distribution
   !
   function normalrnd(mean,stdev) result(c)
      real(8)           :: mean,stdev,c,r,theta
      real(8),parameter :: pi=acos(-1d0)
      if(stdev <= 0d0) then
         Write(*,*) "NORMALRND: Standard Deviation must be +ve"
      else
         r     = (-2d0*log(mersenne()))**0.5d0
         theta = 2d0*pi*mersenne()
         c     = mean+stdev*r*sin(theta)
      end if
   end function normalrnd




   !  Random smaple from an exponential distribution
   !
   function exponentialrnd(mean) result(c)
      real(8) :: mean,c
      if (mean <= 0d0) then
         write(*,*) "EXPONENTIALRND: mean must be positive"
      else
         c=-mean*log(mersenne())
      end if
   end function exponentialrnd



   !  Return a random sample from a gamma distribution
   !
   recursive function gammarnd(shape,scale) result(ans)
      real(8) ::  shape,scale,u,w,d,c,x,xsq,g,v,ans
      if (shape <= 0d0) then
         write(*,*) "GAMMARND: Shape parameter must be positive"
      end if
      if (scale <= 0d0) then
         write(*,*) "GAMMARND: Scale parameter must be positive"
      end if
      !    ## Implementation based on "A Simple Method for Generating Gamma Variables"
      !    ## by George Marsaglia and Wai Wan Tsang.
      !    ## ACM Transactions on Mathematical Software
      !    ## Vol 26, No 3, September 2000, pages 363-372.
      if (shape >= 1d0) then
         d = shape - 1d0/3d0
         c = 1d0/(9d0*d)**0.5d0
         do while (.true.)
            x = normalrnd(0d0,1d0)
            v = 1.d0 + c*x
            do while (v <= 0d0)
               x = normalrnd(0d0,1d0)
               v = 1d0 + c*x
            end do
            v = v*v*v
            u = mersenne()
            xsq = x*x
            if ((u < 1d0 -.0331d0*xsq*xsq) .or. (log(u) < 0.5d0*xsq + d*(1d0 - v + log(v))) )then
               ans=scale*d*v
               return
            end if
         end do
      else
         g = gammarnd(shape+1d0,1d0)
         w = mersenne()
         ans=scale*g*(w)**(1d0/shape)
         return
      end if
   end function gammarnd



   ! ## return a random sample from a chi square distribution
   ! ## with the specified degrees of freedom
   !
   function chi_squarernd(dof) result(ans)
      real(8) :: ans,dof
      ans=gammarnd(0.5d0,2d0*dof)
   end function chi_squarernd



   ! ## return a random sample from an inverse gamma random variable
   !
   function inverse_gammarnd(shape,scale) result(ans)
      real(8) :: shape,scale,ans
      ! If X is gamma(shape, scale) then
      ! 1/Y is inverse gamma(shape, 1/scale)
      ans= 1d0/gammarnd(shape,1d0/scale)
   end function inverse_gammarnd



   !## return a sample from a Weibull distribution
   !
   function weibullrnd(shape,scale) result(ans)
      real(8) :: shape,scale,temp,ans
      if (shape <= 0d0) then
         write(*,*) "WEIBULLRND: Shape parameter must be positive"
      end if
      if (scale <= 0d0) then
         write(*,*) "WEIBULLRND: Scale parameter must be positive"
      end if
      ans= scale * (-log(mersenne()))**(1.d0/shape)
   end function weibullrnd

   !## return a random sample from a Cauchy distribution
   !
   function cauchyrnd(median,scale) result(ans)
      real(8) :: ans,median,scale,p
      if (scale <= 0d0) then
         write(*,*) "CAUCHYRND: Scale parameter must be positive"
      end if
      ans = median + scale*tan(PI*(mersenne() - 0.5))
   end function cauchyrnd

   !
   !## return a random sample from a Student t distribution
   !
   function student_trnd(dof) result(ans)
      real(8) :: ans,dof,y1,y2
      if (dof <= 0.d0) then
         write(*,*) "STUDENT_TRND: Degrees of freedom must be positive"
      end if
      !
      ! ## See Seminumerical Algorithms by Knuth
      y1 = normalrnd(0d0,1d0)
      y2 = chi_squarernd(dof)
      ans= y1 / (y2 / dof)**0.5d0
      !
   end function student_trnd

   !
   !## return a random sample from a Laplace distribution
   !## The Laplace distribution is also known as the double exponential distribution.
   !
   function laplacernd(mean,scale)  result(ans)
      real(8) :: ans,mean,scale,u
      if (scale <= 0d0) then
         write(*,*) "LAPLACERND: Scale parameter must be positive"
      end if
      u = mersenne()
      if (u < 0.5d0) then
         ans = mean + scale*log(2d0*u)
      else
         ans = mean - scale*log(2*(1-u))
      end if

   end function laplacernd

   !
   ! ## return a random sample from a log-normal distribution
   !
   function log_normalrnd(mu,sigma) result(ans)
      real(8) :: ans,mu,sigma
      ans= exp(normalrnd(mu,sigma))
   end function log_normalrnd

   !
   ! ## return a random sample from a beta distribution
   !
   function betarnd(a,b) result(ans)
      real(8) :: a,b,ans,u,v
      if ((a <= 0d0) .or. (b <= 0d0)) then
         write(*,*) "BETARND: Beta parameters must be positive"
      end if
      !    ## There are more efficient methods for generating beta samples.
      !    ## However such methods are a little more efficient and much more complicated.
      !    ## For an explanation of why the following method works, see
      !    ## http://www.johndcook.com/distribution_chart.html#gamma_beta
      u = gammarnd(a,1d0)
      v = gammarnd(b,1d0)
      ans = u / (u + v)
   end function betarnd



   !---------------------------------------------------------------!
   ! State saving subroutines.                                     !
   !                                                               !
   !  Usage:  call mtsave( file_name, format_character )           !
   !     or   call mtsave( unit_number, format_character )         !
   !  where   format_character = 'u' or 'U' will save in           !
   !          unformatted form, otherwise state information will   !
   !          be written in formatted form.                        !
   !---------------------------------------------------------------!
   subroutine mtsavef( fname,forma )
      !NOTE: This subroutine APPENDS to the end of the file "fname".
      character(*),intent(in) :: fname
      character,intent(in)    :: forma
      select case (forma)
       case('u','U')
         open(unit=10,file=trim(fname),status='UNKNOWN',form='UNFORMATTED',position='APPEND')
         write(10) mti
         write(10) mt
       case default
         open(unit=10,file=trim(fname),status='UNKNOWN',form='FORMATTED',position='APPEND')
         write(10,*) mti
         write(10,*) mt
      end select
      close(10)
      return
   end subroutine mtsavef

   subroutine mtsaveu(unum,forma)
      integer,intent(in)    :: unum
      character,intent(in)  :: forma
      select case (forma)
       case('u','U')
         write(unum) mti
         write(unum) mt
       case default
         write(unum,*) mti
         write(unum,*) mt
      end select
      return
   end subroutine mtsaveu



   !---------------------------------------------------------------!
   ! State getting subroutines.                                    !
   !                                                               !
   !  Usage:  call mtget( file_name, format_character )            !
   !     or   call mtget( unit_number, format_character )          !
   !  where   format_character = 'u' or 'U' will read in           !
   !          unformatted form, otherwise state information will   !
   !          be read in formatted form.                           !
   !---------------------------------------------------------------!
   subroutine mtgetf(fname,forma)
      character(*),intent(in) :: fname
      character,intent(in)    :: forma
      select case (forma)
       case('u','U')
         open(unit=10,file=trim(fname),status='OLD',form='UNFORMATTED')
         read(10) mti
         read(10) mt
       case default
         open(unit=10,file=trim(fname),status='OLD',form='FORMATTED')
         read(10,*) mti
         read(10,*) mt
      end select
      close(10)
      return
   end subroutine mtgetf

   subroutine mtgetu(unum,forma)
      integer,intent(in)    :: unum
      character,intent(in)  :: forma
      select case (forma)
       case('u','U')
         read(unum) mti
         read(unum) mt
       case default
         read(unum,*) mti
         read(unum,*) mt
      end select
      return
   end subroutine mtgetu







   !+-----------------------------------------------------------------+
   !PURPOSE:  Numerical Recipes
   !+-----------------------------------------------------------------+
   real(8) function nrand(dseed_)
      integer,optional        :: dseed_
      integer                 :: dseed
      integer,parameter       :: IM1=2147483563, IM2=2147483399, IMM1=IM1-1, IA1=40014, &
      & IA2=40692, IQ1=53668, IQ2=52774, IR1=12211, IR2=3791,  &
      & NTAB=32, NDIV=1+IMM1/NTAB
      real(kind=8), parameter :: AM=1.0d0/IM1, EPS=1.2e-7, RNMX=1.-EPS
      integer                 :: dseed2, j, k, iv(NTAB), iy
      save iv, iy, dseed2
      data dseed2/123456789/, iv/NTAB*0/, iy/0/
      dseed = 123456 ;if(present(dseed_))dseed=dseed_
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


end module SF_RANDOM


