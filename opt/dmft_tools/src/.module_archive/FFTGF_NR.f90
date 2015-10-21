module FFTGF_NR
  USE CONSTANTS, only: pi,one,zero,xi
  USE ARRAYS, only: linspace
  USE INTERPOLATE, only: linear_spline,cubic_spline

  implicit none 
  private


  !GREEN'S FUNCTION FFT:
  !- - - - - - - - - - - - - - - - - - - - - - - - - - - - -!
  !Real freq.      <==> real time
  public :: fftgf_rw2rt
  public :: fftgf_rt2rw
  !- - - - - - - - - - - - - - - - - - - - - - - - - - - - -!
  !Matsubara freq. <==> imaginary time
  public :: fftgf_iw2tau
  public :: fftgf_tau2iw
  !- - - - - - - - - - - - - - - - - - - - - - - - - - - - -!
  public :: fftff_iw2tau , fftff_tau2iw
  public :: fftff_iw2tau_



contains



  !*********************************************************************
  !               FOURIER TRANSFORM OF GREEN'S FUNCTIONS:
  !*********************************************************************  
  !+-------------------------------------------------------------------+
  !PURPOSE  : Evaluate the FFT of a given Green's function from 
  !real frequencies to real time.
  !+-------------------------------------------------------------------+
  subroutine fftgf_rw2rt(func_in,func_out)
    complex(8),dimension(:),intent(inout)        :: func_in
    complex(8),dimension(size(func_in)),optional :: func_out
    complex(8),dimension(size(func_in))          :: ftmp
    integer                                      :: i,M
    M=size(func_in)
    call nr_radix2_test(M)
    ftmp = func_in
    call cfour1(ftmp,M,-1)
    call fftex(ftmp)
    if(present(func_out))then
       func_out = fftshift(ftmp)
    else
       func_in = fftshift(ftmp)
    endif
  end subroutine fftgf_rw2rt


  !+-------------------------------------------------------------------+
  !PURPOSE  : Evaluate the FFT of a given Green's function from 
  !real time to real frequency.
  !+-------------------------------------------------------------------+
  subroutine fftgf_rt2rw(func_in,func_out)
    complex(8),dimension(:),intent(inout)        :: func_in
    complex(8),dimension(size(func_in)),optional :: func_out
    complex(8),dimension(size(func_in))          :: ftmp
    integer                                      :: i,M
    M=size(func_in)
    call nr_radix2_test(M)
    ftmp=func_in
    call cfour1(ftmp,M,1)
    call fftex(ftmp)
    if(present(func_out))then
       func_out = ifftshift(ftmp)/size(ftmp)
    else
       func_in  = ifftshift(ftmp)/size(ftmp)
    endif
  end subroutine fftgf_rt2rw












  !+-------------------------------------------------------------------+
  !PURPOSE  : Evaluate the FFT of a given function from Matsubara frequencies
  ! to imaginary time. 
  !COMMENT  : The routine implements all the necessary manipulations required 
  !by the FFT in this formalism: tail subtraction and reshaping. 
  !Output is only for tau\in[0,beta]
  !if g(-tau) is required this has to be implemented in the calling code
  !using the transformation: g(-tau)=-g(beta-tau) for tau>0
  !+-------------------------------------------------------------------+
  subroutine fftgf_iw2tau(gw,gt,beta,notail,nofix)
    complex(8),dimension(:)     :: gw
    real(8),dimension(:)        :: gt
    real(8)                     :: beta
    logical,optional            :: notail,nofix
    logical                     :: notail_,nofix_
    integer                     :: Liw,Ltau
    notail_=.false. ; if(present(notail))notail_=notail
    nofix_ =.false. ; if(present(nofix))nofix_=nofix
    Liw=size(gw)
    Ltau=size(gt)-1
    call fftgf_iw2tau_(gw,Liw,gt,Ltau,beta,notail_,nofix_)
  end subroutine fftgf_iw2tau
  subroutine fftgf_iw2tau_(gw,N,gt,L,beta,notail_,nofix_)    
    complex(8),dimension(N)       :: gw
    integer                       :: N
    real(8),dimension(L+1)        :: gt
    integer                       :: L
    real(8)                       :: beta
    logical                       :: notail_,nofix_
    complex(8),dimension(2*L)     :: tmpGw
    real(8),dimension(2*L)        :: tmpGt
    integer                       :: i
    real(8)                       :: wmax,mues,fmom
    complex(8)                    :: tail
    real(8)                       :: tau,dtau,At,w
    if(N<L)stop "fftgf_iw2tau: Liw must be > Ltau"
    dtau=beta/dble(L+1)
    wmax=pi/beta*dble(2*L-1)
    mues=-dreal(gw(L-1))*wmax**2
    fmom=-dimag(gw(L-1))*wmax
    tmpGw=dcmplx(0.d0,0.d0)
    select case(notail_)
    case default
       do i=1,L-1
          w=pi/beta*dble(2*i-1)
          tail=-fmom*dcmplx(mues,w)/(mues**2+w**2)
          tmpGw(2*i) =(gw(i)-tail)
       enddo
       call fft(tmpGw)
       tmpGt = -dreal(tmpGw)*2.d0/beta       
       do i=1,L
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
          gt(i)=tmpGt(i)+fmom*At
       enddo
       if(.not.nofix_)gt(L+1)=-(gt(1)+1.d0)
    case (.true.)
       if(N/=L)stop "fftgf_iw2tau: notail requires Liw==Ltau"
       tmpGw=dcmplx(0.d0,0.d0)
       forall(i=1:L-1)tmpGw(2*i)=gw(i)
       call fft(tmpGw)
       tmpGt = dreal(tmpGw)*size(tmpGw)*2.d0/beta
       gt(1:L) = tmpGt(1:L)
       if(.not.nofix_)gt(L)=-gt(1)
    end select
  end subroutine fftgf_iw2tau_


  subroutine fftff_iw2tau_(gw,gt,beta)
    integer                             :: i,n,L,itau,M
    complex(8),dimension(:)             :: gw
    real(8),dimension(size(gw))         :: wm
    complex(8),dimension(0:)            :: gt
    real(8)                             :: wmax,beta,mues,tau,dtau,At,w
    n=size(gw) ; L=size(gt)-1
    if(L/=N)then
       print*,"error in fftgf_iw2tau: call w/ notail and L/=N"
       stop
    endif
    dtau=beta/real(L,8) 
    gt = (0.d0,0.d0)
    forall(i=1:n)wm(i)=pi/beta*real(2*i-1,8)
    forall(i=0:L)gt(i)=sum(cos(real(i,8)*dtau*wm(:))*gw(:))
    gt=gt*2.d0/beta
  end subroutine fftff_iw2tau_

  subroutine fftff_iw2tau(gw,gt,beta)
    integer                             :: i,n,L,itau,M
    complex(8),dimension(:)             :: gw
    real(8),dimension(2*size(gw))       :: reF,imF,tau_long
    complex(8),dimension(0:size(gw))    :: gt
    real(8)                             :: beta,dtau
    real(8),dimension(0:size(gw))       :: tau_short
    N=size(gw) ; L=size(gt)-1
    if(L/=N)then
       print*,"error in fftff_iw2tau: L/=N"
       stop
    endif
    N=2*L
    reF=0.d0 ; imF=0.d0
    forall(i=1:L)
       reF(2*i) = real(gw(i),8)
       imF(2*i) = aimag(gw(i))
    end forall
    call cosft2(reF,2*L,-1)
    call cosft2(imF,2*L,-1)
    tau_long = linspace(0.d0,beta,2*L)
    tau_short= linspace(0.d0,beta,L+1)
    call linear_spline(tau_long,dcmplx(reF,imF),tau_short,gt)
    gt(L)=-gt(0)
    gt=gt/beta*2.d0
  end subroutine fftff_iw2tau




  !+-------------------------------------------------------------------+
  !PURPOSE  : 
  !+-------------------------------------------------------------------+
  subroutine fftgf_tau2iw(gt,gw,beta)
    real(8)    :: gt(0:)
    complex(8) :: gw(:)
    integer    :: Liw,Ltau
    real(8)    :: beta
    Liw = size(gw)
    Ltau= size(gt)-1
    call fftgf_tau2iw_(gt(0:),Ltau,gw,Liw,beta)
  end subroutine fftgf_tau2iw
  subroutine fftgf_tau2iw_(gt,L,gw,N,beta)
    real(8)                :: gt(0:L)
    integer                :: L
    complex(8)             :: gw(N),one=(1.d0,0.d0)
    integer                :: N
    real(8)                :: beta
    integer                :: i,M
    real(8),allocatable    :: xgt(:)
    complex(8),allocatable :: Igw(:)
    real(8),allocatable    :: Igt(:)
    M=32*L
    allocate(xgt(2*L),Igt(0:2*M),Igw(2*M))
    forall(i=0:L)
       xgt(L+i)= gt(i)
       xgt(i)  =-gt(L-i)                      !Valid for fermionic GF (no bosonic)
    end forall
    call interp_gtau(xgt(L+1:2*L),Igt(M+1:2*M),L,M) !interpolate to treat long tails
    forall(i=0:M)Igt(i)=-Igt(2*M-i)           !Valid for fermionic GF (no bosonic)
    Igw=one*Igt
    call ifft(Igw)
    Igw = Igw/size(Igw)*beta
    forall(i=1:L)gw(i)=Igw(2*i)
    deallocate(Xgt,Igt,Igw)
  contains
    subroutine interp_gtau(Fin,Fout,Lin,Lout)
      integer :: Lin,Lout
      real(8) :: Fin(Lin),Xin(Lin)
      real(8) :: Fout(Lout),Xout(Lout)
      Xin = linspace(0.d0,beta,Lin)
      Xout= linspace(0.d0,beta,Lout)
      call cubic_spline(Xin,Fin,Xout,Fout)
    end subroutine interp_gtau
  end subroutine fftgf_tau2iw_


  subroutine fftff_tau2iw(gt,gw,beta)
    complex(8)             :: gt(0:)
    complex(8)             :: gw(:)
    real(8)                :: beta
    integer                :: i,L,n,M
    real(8),allocatable    :: reGt(:),imGt(:)
    complex(8),allocatable :: Igw(:),Igt(:)
    L=size(gt)-1    ; N=size(gw)
    M=32*L
    allocate(Igt(-M:M),Igw(2*M))
    allocate(reGt(0:M),imGt(0:M))
    call interp_gtau(dreal(gt(0:L)),reGt(0:M),L,M)
    call interp_gtau(dimag(gt(0:L)),imGt(0:M),L,M)
    Igt(0:M)=cmplx(reGt(0:M),imGt(0:M),8)
    !
    forall(i=1:M)Igt(-i)=-Igt(M-i) !Valid for every fermionic GF (bosonic case not here)
    call fftgf_rt2rw(Igt,Igw)
    Igw=Igw*beta/real(M,8)/2.d0
    forall(i=1:n)gw(i)=Igw(2*i)
    deallocate(Igt,Igw)
  contains
    subroutine interp_gtau(Fin,Fout,Lin,Lout)
      integer :: Lin,Lout
      real(8) :: Fin(Lin),Xin(Lin)
      real(8) :: Fout(Lout),Xout(Lout)
      Xin = linspace(0.d0,beta,Lin)
      Xout= linspace(0.d0,beta,Lout)
      call cubic_spline(Xin,Fin,Xout,Fout)
    end subroutine interp_gtau
  end subroutine fftff_tau2iw




  !*******************************************************************
  !*******************************************************************
  !*******************************************************************




  !+-------------------------------------------------------------------+
  !PROGRAM  : FOUR1, COSFT2, REALFT
  !TYPE     : Subroutines
  !PURPOSE  : adapted from Num. Rec. 
  !+-------------------------------------------------------------------+
  subroutine cosft2(y,n,isign)
    integer :: i,isign,n
    real(8) :: y(n)
    real(8) ::  sum,sum1,y1,y2,ytemp
    real(8) ::  theta,wi,wi1,wpi,wpr,wr,wr1,wtemp
    theta=0.5d0*pi/n
    wr=1.0d0
    wi=0.0d0
    wr1=cos(theta)
    wi1=sin(theta)
    wpr=-2.0d0*wi1**2
    wpi=sin(2.d0*theta)
    if(isign.eq.1)then
       do i=1,n/2
          y1=0.5*(y(i)+y(n-i+1))
          y2=wi1*(y(i)-y(n-i+1))
          y(i)=y1+y2
          y(n-i+1)=y1-y2
          wtemp=wr1
          wr1=wr1*wpr-wi1*wpi+wr1
          wi1=wi1*wpr+wtemp*wpi+wi1
       enddo
       call realft(y,n,1)
       do i=3,n,2
          wtemp=wr
          wr=wr*wpr-wi*wpi+wr
          wi=wi*wpr+wtemp*wpi+wi
          y1=y(i)*wr-y(i+1)*wi
          y2=y(i+1)*wr+y(i)*wi
          y(i)=y1
          y(i+1)=y2
       enddo
       sum=0.5*y(2)
       do i=n,2,-2
          sum1=sum
          sum=sum+y(i)
          y(i)=sum1
       enddo
    else if(isign.eq.-1)then
       ytemp=y(n)
       do i=n,4,-2
          y(i)=y(i-2)-y(i)
       enddo
       y(2)=2.0*ytemp
       do i=3,n,2
          wtemp=wr
          wr=wr*wpr-wi*wpi+wr
          wi=wi*wpr+wtemp*wpi+wi
          y1=y(i)*wr+y(i+1)*wi
          y2=y(i+1)*wr-y(i)*wi
          y(i)=y1
          y(i+1)=y2
       enddo
       call realft(y,n,-1)
       do i=1,n/2
          y1=y(i)+y(n-i+1)
          y2=(0.5/wi1)*(y(i)-y(n-i+1))
          y(i)=0.5*(y1+y2)
          y(n-i+1)=0.5*(y1-y2)
          wtemp=wr1
          wr1=wr1*wpr-wi1*wpi+wr1
          wi1=wi1*wpr+wtemp*wpi+wi1
       enddo
    endif
  end subroutine cosft2

  subroutine realft(data,n,isign)
    integer :: isign,n
    real(8) :: data(n)
    integer ::  i,i1,i2,i3,i4,n2p3
    real(8) ::  c1,c2,h1i,h1r,h2i,h2r,wis,wrs
    real(8) ::  theta,wi,wpi,wpr,wr,wtemp
    theta=3.141592653589793d0/dble(n/2)
    c1=0.5
    if (isign.eq.1) then
       c2=-0.5
       call rfour1(data,n/2,+1)
    else
       c2=0.5
       theta=-theta
    endif
    wpr=-2.0d0*sin(0.5d0*theta)**2
    wpi=sin(theta)
    wr=1.0d0+wpr
    wi=wpi
    n2p3=n+3
    do i=2,n/4
       i1=2*i-1
       i2=i1+1
       i3=n2p3-i2
       i4=i3+1
       wrs=sngl(wr)
       wis=sngl(wi)
       h1r=c1*(data(i1)+data(i3))
       h1i=c1*(data(i2)-data(i4))
       h2r=-c2*(data(i2)+data(i4))
       h2i=c2*(data(i1)-data(i3))
       data(i1)=h1r+wrs*h2r-wis*h2i
       data(i2)=h1i+wrs*h2i+wis*h2r
       data(i3)=h1r-wrs*h2r+wis*h2i
       data(i4)=-h1i+wrs*h2i+wis*h2r
       wtemp=wr
       wr=wr*wpr-wi*wpi+wr
       wi=wi*wpr+wtemp*wpi+wi
    enddo
    if (isign.eq.1) then
       h1r=data(1)
       data(1)=h1r+data(2)
       data(2)=h1r-data(2)
    else
       h1r=data(1)
       data(1)=c1*(h1r+data(2))
       data(2)=c1*(h1r-data(2))
       call rfour1(data,n/2,-1)
    endif
    return
  end subroutine realft

  !This routine is modified version w/ respect to original NumRec!
  !Below is the original one renamed rfour1 (real input twice the size)
  subroutine cfour1(Fct,nn,isign)
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
       Fct(i) = cmplx(Fcttmp(2*i-1),Fcttmp(2*i),8)
    enddo
  end subroutine cfour1

  subroutine rfour1(data,nn,isign)
    integer :: isign,nn
    real(8) :: data(2*nn)
    integer :: i,istep,j,m,mmax,n
    real(8) ::  tempi,tempr
    real(8) ::  theta,wi,wpi,wpr,wr,wtemp
    n=2*nn
    j=1
    do i=1,n,2
       if(j.gt.i)then
          tempr=data(j)
          tempi=data(j+1)
          data(j)=data(i)
          data(j+1)=data(i+1)
          data(i)=tempr
          data(i+1)=tempi
       endif
       m=nn
1      if ((m.ge.2).and.(j.gt.m)) then
          j=j-m
          m=m/2
          goto 1
       endif
       j=j+m
    enddo
    mmax=2
2   if (n.gt.mmax) then
       istep=2*mmax
       theta=6.28318530717959d0/(isign*mmax)
       wpr=-2.d0*sin(0.5d0*theta)**2
       wpi=sin(theta)
       wr=1.d0
       wi=0.d0
       do m=1,mmax,2
          do i=m,n,istep
             j=i+mmax
             tempr=sngl(wr)*data(j)-sngl(wi)*data(j+1)
             tempi=sngl(wr)*data(j+1)+sngl(wi)*data(j)
             data(j)=data(i)-tempr
             data(j+1)=data(i+1)-tempi
             data(i)=data(i)+tempr
             data(i+1)=data(i+1)+tempi
          enddo
          wtemp=wr
          wr=wr*wpr-wi*wpi+wr
          wi=wi*wpr+wtemp*wpi+wi
       enddo
       mmax=istep
       goto 2
    endif
    return
  end subroutine rfour1

  !*******************************************************************
  !*******************************************************************
  !*******************************************************************


end module FFTGF_NR
