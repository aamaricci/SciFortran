!+-------------------------------------------------------------------+
!PURPOSE  : Evaluate the FFT of a given function from Matsubara frequencies
! to imaginary time. 
!+-------------------------------------------------------------------+
subroutine fft_iw2tau(gw,gt,beta,L)
  implicit none
  complex(8),dimension(L)    :: gw
  complex(8),dimension(2*L)  :: tmpGw
  real(8),dimension(0:L)     :: gt
  complex(8),dimension(-L:L) :: tmpGt
  integer    :: i,L
  real(8)    :: pi,beta
  complex(8) :: tail,zero
  pi=acos(-1.d0) ; zero=cmplx(0.d0,0.d0)
  do i=1,L
     tmpGw(2*i)  = gw(i)
     tmpGw(2*i-1)= zero
  enddo
  call four1(tmpGw,2*L,-1)
  forall(i=1:2*L)tmpGt(i-L-1)=2.d0*dble(tmpGw(i))/beta
  gt(0:L-1) = real(tmpGt(0:L-1),8)
  gt(L)=-gt(0) !fix the end point:
end subroutine fft_iw2tau




!*******************************************************************
!*******************************************************************
!*******************************************************************




!+-------------------------------------------------------------------+
!PURPOSE  :  
!+-------------------------------------------------------------------+
subroutine fft_tau2iw(gt,gw,beta,L)
  implicit none
  real(8)                :: gt(0:L),xgt(-L:L)
  complex(8)             :: gw(2*L),xgw(2*L)
  real(8)                :: beta
  integer                :: i,L,n,M
  !Get eXtended G(tau) 0:L --> -L:L
  xgt(0:L)=gt(0:L) ; forall(i=1:L)xgt(-i)=-xgt(L-i)
  !FFT to Matsubara freq.
  forall(i=1:2*L)xgw(i)=(1.d0,0.d0)*xgt(i-L-1)
  call four1(xgw,2*L,1)
  !Exit:
  forall(i=1:L)gw(2*i)=xgw(2*i)*beta/dble(L)/2.d0
end subroutine fft_tau2iw




!*******************************************************************
!*******************************************************************
!*******************************************************************



!+-------------------------------------------------------------------+
!PROGRAM  : FOUR1
!TYPE     : Subroutine
!PURPOSE  : adapted from Num. Rec. 
!+-------------------------------------------------------------------+
subroutine four1(Fct,nn,isign)
  implicit none
  integer    :: nn, isign
  complex(8) :: Fct(nn)
  integer    :: i, istep, j, m, mmax, n
  real(8)    :: Fcttmp(2*nn)
  real(8)    :: tempi, tempr
  real(8)    :: theta, wi, wpi, wpr, wr, wtemp
  do i=1, nn
     Fcttmp(2*i-1) = real(Fct(i))
     Fcttmp(2*i)   = aimag(Fct(i))
  enddo
  n = 2*nn
  j = 1
  do i=1, n, 2
     if(j .gt. i)then
        tempr = Fcttmp(j)
        tempi = Fcttmp(j+1)
        Fcttmp(j)   = Fcttmp(i)
        Fcttmp(j+1) = Fcttmp(i+1)
        Fcttmp(i)   = tempr
        Fcttmp(i+1) = tempi
     endif
     m = n/2
     do while((m .ge. 2) .and. (j .gt. m))
        j = j-m
        m = m/2
     enddo
     j = j+m
  enddo
  mmax = 2
  do while(n .gt. mmax)
     istep = 2*mmax
     theta = 6.28318530717959d0/(isign*mmax)
     wpr = -2.d0*dsin(0.5d0*theta)**2
     wpi = dsin(theta)
     wr = 1.d0
     wi = 0.d0
     do m=1, mmax, 2
        do i=m, n, istep
           j = i+mmax
           tempr = sngl(wr)*Fcttmp(j)-sngl(wi)*Fcttmp(j+1)
           tempi = sngl(wr)*Fcttmp(j+1)+sngl(wi)*Fcttmp(j)
           Fcttmp(j)   = Fcttmp(i)-tempr
           Fcttmp(j+1) = Fcttmp(i+1)-tempi
           Fcttmp(i)   = Fcttmp(i)+tempr
           Fcttmp(i+1) = Fcttmp(i+1)+tempi
        enddo
        wtemp = wr
        wr = wr*wpr-wi*wpi+wr
        wi = wi*wpr+wtemp*wpi+wi
     enddo
     mmax = istep
  enddo
  do i=1, nn
     Fct(i) = cmplx(Fcttmp(2*i-1),Fcttmp(2*i))
  enddo
end subroutine four1
