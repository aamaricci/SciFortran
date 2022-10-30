program interpALL
  USE SCIFOR
  USE ASSERTING
  implicit none

  integer,parameter               :: Lin=50, Lout=2*Lin-1, Np=10
  real(8),parameter               :: dt=0.1d0
  real(8),dimension(Lin)          :: xin,yin
  real(8),dimension(Lin)          :: xfin
  complex(8),dimension(Lin)       :: cxfin
  real(8),dimension(Lin,Lin)      :: fin
  complex(8),dimension(Lin,Lin)   :: cfin
  real(8),dimension(Lout)         :: xout,yout
  real(8),dimension(Lout)         :: xfout
  complex(8),dimension(Lout)      :: cxfout

  real(8),dimension(Lout,Lout)    :: fout
  complex(8),dimension(Lout,Lout) :: cfout
  real(8)                         :: xmax

  xmax=dt*real(Lin-1,8)

  xin = linspace(0.d0,xmax,Lin)
  yin = linspace(0.d0,xmax,Lin)
  !Exact functions
  xfin = sin(xin)*cos(xin)  
  cxfin= exp(-xi*0.5d0*xin)
  !Poly Interp arrays:
  xout = linspace(0.d0,xmax,Lout)
  yout = linspace(0.d0,xmax,Lout)
  print*,"1D interpolate:"
  
  call linear_spline(xin,xfin,xout,xfout)
  call linear_spline(xin,cxfin,xout,cxfout)
  call assert(xfout,sin(xout)*cos(xout),"Sin(x)Cos(x) Linear interpolation N=10:",tol=1d-2)
  call assert(cxfout,exp(-xi*xout*0.5d0),"Exp(-i.x/2) Linear interpolation N=10:",tol=1d-2)

  call poly_spline(xin,xfin,xout,xfout,Np)
  call poly_spline(xin,cxfin,xout,cxfout,Np)
  call assert(xfout,sin(xout)*cos(xout),"Sin(x)Cos(x) Polynomial interpolation N=10:",tol=1d-10)
  call assert(cxfout,exp(-xi*xout*0.5d0),"Exp(-i.x/2) Polynomial interpolation N=10:",tol=1d-10)

  call cubic_spline(xin,xfin,xout,xfout)
  call cubic_spline(xin,cxfin,xout,cxfout)
  call assert(xfout,sin(xout)*cos(xout),"Sin(x)Cos(x) Cubic interpolation N=10:",tol=1d-5)
  call assert(cxfout,exp(-xi*xout*0.5d0),"Exp(-i.x/2) Cubic interpolation N=10:",tol=1d-5)


  print*,"2D interpolate:"
  fin = outerprod(sin(xin),cos(yin))
  call linear_spline(xin,yin,fin,xout,yout,fout)
  call assert(fout,outerprod(sin(xout),cos(yout)),"Outprod(Sin(x)Cos(y)) Linear interpolation 2D:",tol=0.025d0)

  cfin= outerprod(exp(-xi*0.5d0*xin),exp(xi*0.5d0*yin))
  call poly_spline(xin,yin,fin,xout,yout,fout,Np)
  call assert(fout,outerprod(sin(xout),cos(yout)),"Oprod(Sin(x),Cos(y)) Poly interp 2D N=10:",tol=1d-10)
  call poly_spline(xin,yin,cfin,xout,yout,cfout,Np)
  call assert(cfout,outerprod(exp(-xi*0.5d0*xout),exp(xi*0.5d0*yout)),"Oprod(Exp(-i*x/2),exp(iy/2)) Poly interp 2D N=10:",tol=1d-10)



end program interpALL
