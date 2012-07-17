      subroutine cfft_1d_forward(func,L)
      implicit none
      integer i,L
      complex*16 func(L)
      real*8 wsave(5*L)
      call zffti(L,wsave)
      call zfftf(L,func,wsave)
      end 


      subroutine cfft_1d_backward(func,L)
      implicit none
      integer i,L
      complex*16 func(L)
      real*8 wsave(5*L)
      call zffti(L,wsave)
      call zfftb(L,func,wsave)
      end 


      subroutine cfft_1d_shift(fin,fout,L)
      integer i,L
      complex*16 fin(2*L),fout(-L:L),dout(-L:L)
      do i=1,2*L
         dout(i-L-1)=fin(i)     ![1,2*L]---> [-L,L-1]
      enddo
      do i=-L,-1
         fout(i+L)=dout(i)      !g[0,M-1]<--- x[-M,-1]
      enddo
      do i=0,L-1
         fout(i-L)=dout(i)      !g[-L,-1]<--- x[0,L-1]         
      enddo
      end


      subroutine swap_fftrt2rw(func_in,Nsize)
      integer i,Nsize,Nhalf
      complex*16 func_in(Nsize),dummy(Nsize)
      Nhalf=Nsize/2
      do i=1,Nsize
         dummy(i)=func_in(i)
      enddo
      do i=1,Nhalf
         func_in(i)=dummy(Nhalf+i)
         func_in(Nhalf+i)=dummy(i)
      enddo
      end


!     +-------------------------------------------------------------------+
!     PROGRAM  : CFFT_IW2IT
!     TYPE     : Subroutine
!     PURPOSE  : Evaluate the FFT of a given function from Matsubara frequencies
!     to imaginary time. 
!     "gw"= the GF to FFT. dimension  NMAX (fg(NMAX))
!     "gt"= output function. dimension 0:n (funct(0:N))
!     "beta"= inverse temperature (directly passed)
!     COMMENT  : The routine implements all the necessary manipulations required 
!     by the FFT in this formalism: tail subtraction and reshaping. 
!     Output is only for tau\in[0,beta]
!     if g(-tau) is required this has to be implemented in the calling code
!     using the transformation: g(-tau)=-g(beta-tau) for tau>0
!     +-------------------------------------------------------------------+
      subroutine fft_iw2tau(gw,gt,beta,L)
      implicit none
      integer i,L
      complex*16 gw(L),tmpGw(2*L),tmpGt(-L:L)
      real*8 gt(0:L)
      real*8 pi,wmax,beta,mues,tau,dtau,At,w
      complex*16 :: tail
      pi=acos(-1.d0)
      dtau=beta/dble(L)
      wmax=pi/beta*dble(2*L-1)
      mues=-dreal(gw(L))*wmax**2
      do i=1,2*L
         tmpGw(i)=(0.d0,0.d0)
      enddo
      do i=1,L
         w=pi/beta*dble(2*i-1)
         tail=-dcmplx(mues,w)/(mues**2+w**2)
         tmpGw(2*i)= gw(i)-tail
      enddo
      call cfft_1d_forward(tmpGw,2*L)
      call cfft_1d_shift(tmpGw,tmpGt,L)
      do i=0,L-1
         tau=dble(i)*dtau
         if(mues > 0.d0)then
            if((mues*beta) > 30.d0)then
               At = -exp(-mues*tau)
            else
               At = -exp(-mues*tau)/(1.d0 + exp(-beta*mues))
            endif
         else
            if((mues*beta) < -30.d0)then
               At = -exp(mues*(beta-tau))
            else
               At = -exp(-mues*tau)/(1.d0 + exp(-beta*mues))
            endif
         endif
         gt(i)=dreal(tmpGt(i))*2.d0/beta+At
      enddo
      gt(L)=-(gt(0)+1.d0)       !fix the end point:
      end subroutine fft_iw2tau
