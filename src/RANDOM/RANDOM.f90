  include 'mkl_vsl.fi' 
  !###############################################################
  !PROGRAM  : RANDOM
  !TYPE     : Module
  !PURPOSE  : Module for Random Number generators
  !###############################################################
  module RANDOM
    USE MKL_VSL_TYPE
    USE MKL_VSL
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

    interface rand
       module procedure is_random,iv_random,im_random,&
            ds_random,dv_random,dm_random,&
            cs_random,cv_random,cm_random
    end interface rand

    public :: nrand
    public :: irand,drand,crand
    public :: init_random_number
    public :: random_order
    public :: rand

  contains

    function irand()
      integer :: irand
      real(8) :: r
      call random_number(r)
      irand=nint(r*10.d0)
    end function irand
    function drand()
      real(8) :: drand
      call random_number(drand)
    end function drand
    function crand()
      complex(8) :: crand
      real(8) :: re,im
      call random_number(re)
      call random_number(im)
      crand=dcmplx(re,im)
    end function crand




    !*******************************************************************
    !*******************************************************************
    !*******************************************************************


    !+-----------------------------------------------------------------+
    !PURPOSE  :  interface to multi-dimensional RNG
    !+-----------------------------------------------------------------+
    subroutine is_random(M,dim)
      integer,intent(inout) :: M
      integer,optional      :: dim
      D=1;if(present(dim))D=dim
      M=nint(drand()*D)
    end subroutine is_random
    subroutine iv_random(M,dim)
      integer,dimension(:),intent(inout) :: M
      integer,optional                   :: dim
      D=1;if(present(dim))D=dim
      do i=1,size(M)
         M(i)=nint(drand()*D)
      end do
    end subroutine iv_random
    subroutine im_random(M,dim)
      integer,dimension(:,:),intent(inout) :: M
      integer,optional                     :: dim
      D=1;if(present(dim))D=dim
      do i=1,size(M,1)
         do j=1,size(M,2)
            M(i,j)=nint(drand()*D)
         end do
      end do
    end subroutine im_random
    subroutine ds_random(M,dim)
      real(8),intent(inout) :: M
      integer,optional      :: dim
      D=1;if(present(dim))D=dim
      M=drand()*D
    end subroutine ds_random
    subroutine dv_random(M,dim)
      real(8),dimension(:),intent(inout) :: M
      integer,optional                   :: dim
      D=1;if(present(dim))D=dim
      do i=1,size(M)
         M(i)=drand()*D
      end do
    end subroutine dv_random
    subroutine dm_random(M,dim)
      real(8),dimension(:,:),intent(inout) :: M
      integer,optional                     :: dim
      D=1;if(present(dim))D=dim
      do i=1,size(M,1)
         do j=1,size(M,2)
            M(i,j)=drand()*D
         end do
      end do
    end subroutine dm_random
    subroutine cs_random(M,dim)
      complex(8),intent(inout) :: M
      integer,optional      :: dim
      D=1;if(present(dim))D=dim
      M=cmplx(drand(),drand())*D
    end subroutine cs_random
    subroutine cv_random(M,dim)
      complex(8),dimension(:),intent(inout) :: M
      integer,optional                   :: dim
      D=1;if(present(dim))D=dim
      do i=1,size(M)
         M(i)=cmplx(drand(),drand())*D
      end do
    end subroutine cv_random
    subroutine cm_random(M,dim)
      complex(8),dimension(:,:),intent(inout) :: M
      integer,optional                     :: dim
      D=1;if(present(dim))D=dim
      do i=1,size(M,1)
         do j=1,size(M,2)
            M(i,j)=cmplx(drand(),drand())*D
         end do
      end do
    end subroutine cm_random


    !***************************************************************
    !***************************************************************
    !***************************************************************


    !+-----------------------------------------------------------------+
    !PURPOSE  :   
    !+-----------------------------------------------------------------+
    real(8) function nrand(dseed)
      implicit none
      integer :: dseed
      integer,      parameter :: IM1=2147483563, IM2=2147483399, IMM1=IM1-1, IA1=40014, &
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




    !*******************************************************************
    !*******************************************************************
    !*******************************************************************









    !+-----------------------------------------------------------------+
    !PURPOSE  : initialize the seed for the INTRINSIC RNG 
    !+-----------------------------------------------------------------+
    subroutine init_random_number(shift)
      integer,optional :: shift
      integer :: i, n, clock
      integer,dimension(:),allocatable :: seed
      call RANDOM_SEED(size = n)
      allocate(seed(n))
      call SYSTEM_CLOCK(COUNT=clock)
      seed = clock + 37 * (/ (i - 1, i = 1, n) /)
      if(present(shift))seed=seed+shift
      call RANDOM_SEED(PUT = seed)
      deallocate(seed)
    end subroutine init_random_number




    !*******************************************************************
    !*******************************************************************
    !*******************************************************************




    !+-----------------------------------------------------------------+
    !PURPOSE  :   
    !+-----------------------------------------------------------------+
    function mkl_random_number_uniform()
      real(4)                 :: mkl_random_number_uniform
      integer                 :: brng_
      real(4)                 :: a_,b_
      real(4),dimension(1)    :: r
      integer                 :: method,seed,errcode
      TYPE (VSL_STREAM_STATE) :: stream
      brng_=VSL_BRNG_MT19937 !;if(present(brng))brng_=brng
      a_=0.0 !; if(present(a))a_=real(a,4)
      b_=1.0 !; if(present(b))b_=real(b,4)
      method=0
      seed=1654
      errcode=vslnewstream( stream, brng_,  seed )       !Initialize
      errcode=vsrnguniform( method, stream, 1, r, a_, b_ )!Call RNG Uniform
      mkl_random_number_uniform=r(1)
      errcode=vsldeletestream( stream )                 !Close
    end function mkl_random_number_uniform




    !*******************************************************************
    !*******************************************************************
    !*******************************************************************




    !+-----------------------------------------------------------------+
    !PURPOSE  :   
    !+-----------------------------------------------------------------+
    SUBROUTINE random_order(order, n)
      !     Generate a random ordering of the integers 1 ... n.
      INTEGER, INTENT(IN)  :: n
      INTEGER, INTENT(OUT) :: order(n)
      !     Local variables
      INTEGER :: i, j, k
      REAL(8) :: wk
      DO i = 1, n
         order(i) = i
      END DO
      !     Starting at the end, swap the current last indicator with one
      !     randomly chosen from those preceeding it.
      DO i = n, 2, -1
         CALL RANDOM_NUMBER(wk)
         j = 1 + i * wk
         IF (j < i) THEN
            k = order(i)
            order(i) = order(j)
            order(j) = k
         END IF
      END DO
      RETURN
    END SUBROUTINE random_order



    !*******************************************************************
    !*******************************************************************
    !*******************************************************************





    !+-----------------------------------------------------------------+
    !PURPOSE  : RNG library, not exported!! 
    !+-----------------------------------------------------------------+
    !include "Routines.f90"



    !*******************************************************************
    !*******************************************************************
    !*******************************************************************


  end module RANDOM
