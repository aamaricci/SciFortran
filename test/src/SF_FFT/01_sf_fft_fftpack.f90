program testFFT
  USE SCIFOR
  USE ASSERTING
  implicit none



  write(*,*)"Testing FFT:"
  call test_notable_functions()   
  write(*,*)""


contains


  subroutine test_notable_functions()
    integer,parameter :: N=4096
    real(8)           :: Ft(N),Fw(N),r(N),rfoo(N),w(N),t(N),dw,dt,wmax,tmax,A
    complex(8)        :: c(N),cfoo(N)
    integer           :: i,j
    dt=0.1d0
    dw=pi/dt/dble(N/2)
    tmax=dt*dble(N/2)
    wmax=dw*dble(N/2)
    t = linspace(-tmax,tmax-dt,N,mesh=dt)
    w = linspace(-wmax,wmax-dw,N,mesh=dw)
    print*,"Test1: Rectangular function"
    A = 10.d0
    do i=1,N
       if(abs(t(i))<=A)Ft(i)=1.d0/A/2.d0
       Fw(i)=sin(A*w(i))/w(i)/A
       if(w(i)==0.d0)Fw(i)=1.0
    enddo
    !Ft and Fw are exact
    ! call splot("Rect_f_w.dat",w,Fw)
    ! call splot("Rect_f_t.dat",t,Ft)

    C = one*Ft
    R = Ft
    call tfft(c);c=c*dt         !FFT complex
    call tfft(r);r=r*dt         !FFT dble
    ! call splot("Rect_FT(f_t).dat",w,abs(r-Fw))
    call assert(dreal(c),Fw,"Rectangular function A=10: cmplx F(w)=FT(F(t))",tol=1d-2)
    call assert(r,Fw,"Rectangular function A=10: dble F(w)=FT(F(t))",tol=1d-2)

    C = one*Fw
    R = Fw
    call itfft(c);c=c/dt
    call itfft(r);r=r/dt
    ! call splot("Rect_iFT(f_w).dat",t,abs(r-Ft))
    call assert(dreal(c),Ft,"Rectangular function A=10: cmplx F(t)=iFT(F(w))",tol=0.03d0)
    call assert(r,Ft,"Rectangular function A=10: dble F(t)=iFT(F(w))",tol=0.03d0)



    !TEST2:
    print*,"Test2: Gaussian function"
    A = 1.d0
    do i=1,N
       Ft(i)= exp(-A*t(i)**2)
       Fw(i)= exp(-w(i)**2/A/4.d0)*sqrt(pi/A)
    enddo
    !Ft and Fw are exact
    ! call splot("Gauss_f_w.dat",w,Fw)
    ! call splot("Gauss_f_t.dat",t,Ft)

    C = one*Ft
    R = Ft
    call tfft(c);c=c*dt         !FFT complex
    call tfft(r);r=r*dt         !FFT dble
    ! call splot("Gauss_FT(f_t).dat",w,abs(r-Fw))
    call assert(dreal(c),Fw,"Gaussian function A=1: cmplx F(w)=FT(F(t))",tol=1d-8)
    call assert(r,Fw,"Gaussian function A=1: dble F(w)=FT(F(t))",tol=1d-8)

    C = one*Fw
    R = Fw
    call itfft(c);c=c/dt
    call itfft(r);r=r/dt
    ! call splot("Gauss_iFT(f_w).dat",t,abs(r-Ft))
    call assert(dreal(c),Ft,"Gaussian function A=1: cmplx F(t)=iFT(F(w))",tol=1d-8)
    call assert(r,Ft,"Gaussian function A=1: dble F(t)=iFT(F(w))",tol=1d-8)


    print*,"Test3: Lorentzian function"
    A = 1.d0
    c = 0.d0 
    r = 0.d0
    do i=1,N
       Ft(i)= exp(-A*abs(t(i)))
       Fw(i)= 2.d0*A/(w(i)**2 + A**2)
    enddo
    !Ft and Fw are exact
    ! call splot("Lorentzian_f_w.dat",w,Fw)
    ! call splot("Lorentzian_f_t.dat",t,Ft)

    C = one*Ft
    R = Ft
    call tfft(c);c=c*dt         !FFT complex
    call tfft(r);r=r*dt         !FFT dble
    call splot("Lorentzian_FT(f_t).dat",w,abs(r-Fw))
    call assert(dreal(c),Fw,"Lorentzian function A=1: cmplx F(w)=FT(F(t))",tol=0.03d0)
    call assert(r,Fw,"Lorentzian function A=1: dble F(w)=FT(F(t))",tol=0.03d0)

    C = one*Fw
    R = Fw
    call itfft(c);c=c/dt
    call itfft(r);r=r/dt
    call splot("Lorentzian_iFT(f_w).dat",t,abs(r-Ft))
    call assert(dreal(c),Ft,"Lorentzian function A=1: cmplx F(t)=iFT(F(w))",tol=0.03d0)
    call assert(r,Ft,"Lorentzian function A=1: dble F(t)=iFT(F(w))",tol=0.03d0)


  end subroutine test_notable_functions



end program testFFT
