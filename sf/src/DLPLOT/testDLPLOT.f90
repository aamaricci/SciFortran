!This is a (rather simple) example about how to use the DLPLOT module: 
program testDLPLOT
  USE DLPLOT
  IMPLICIT NONE
  INTEGER :: NDIG
  INTEGER :: I,J,M
  INTEGER,PARAMETER        :: L=30
  INTEGER, PARAMETER       :: N=100
  REAL(8)    :: ZMAT(N,N),YT(N,N)
  REAL(8)    :: PI,FPI,STEP,X(N),FUNC(N),Y(N),Y1(N),Y2(N),Y3(L,N)
  REAL(8)    :: VX(N,N),VY(N,N)
  REAL(8)    :: VTX(L,N,N),VTY(L,N,N),ZMATT(N,N,L)
  real(4)    :: xmin,xmax,dex
  real(4)    :: ymin,ymax,dey


  PI=acos(-1.d0)
  FPI=PI/180.d0
  STEP = 360.d0/dble(N-1)

  !1-dim array && 2-dim array-notime && vector field:
  DO I=1,N
     Y1(I)=I
     X(I)=(I-1)*STEP
     FUNC(I)=COS(X(I)*FPI)
     Y1(I)=FUNC(I)*EXP(-0.01*X(I))
     Y2(I)=COS(X(I)*FPI)*SIN(X(I)*FPI)
     DO J=1,N
        Y(J)=(J-1)*STEP
        YT(I,J)=Y(J)-2.*DBLE(I)
        ZMAT(I,J)=COS(X(I)*FPI) + COS(Y(J)*FPI)
        VX(I,J)=COS(X(I)*FPI)
        VY(I,J)=SIN(Y(J)*FPI)
     ENDDO
  ENDDO

  !1-dim array + time dependence:
  DO I=1,L
     DO J=1,N
        Y3(I,J)=SIN(X(J)*FPI + DBLE(I))
     ENDDO
  ENDDO

  !2-dim arrays + time dependence:
  DO M=1,L
     DO I=1,N
        X(I)=(I-1)*STEP
        DO J=1,N
           Y(J)=(J-1)*STEP
           ZMATT(I,J,M)=COS(X(I)*FPI+ 0.1*DBLE(M)) + COS(Y(J)*FPI+ 0.1*DBLE(M))
           VTX(M,I,J)=COS(X(I)*FPI + 0.1*DBLE(M))
           VTY(M,I,J)=SIN(Y(J)*FPI + 0.1*DBLE(M)) 
        ENDDO
     ENDDO
  ENDDO

  ! !1-dim array: PLOT
  call dplot("plotFUNCS",X,FUNC,"c(x)",Y1,"$c(x)*e(-x)$",Y2,"$c(x)*s(x)$",wlp="wlp") 

  ! !1-dim + animation:
  call dplot_movie("plotANI",X,Y3)

  ! !2-dim: plot3D
  ! call dplot_3d("plotZMAT","X","Y","Z",X,Y,ZMAT)
  ! !change the point of view of the "observer": s1,s2 
  ! call dplot_3d("plotZMAT_newVIEW","X","Y","Z",X,Y,ZMAT,theta=250.d0,phi=40.d0,rho=5.6d0)
  ! !produce a rotating figure in addition to the static plots for better visualization
  call dplot_3D("plotZMAT","X","Y","Z",X,Y,ZMAT,imovie=.true.)

  ! !2-dim + time animatino:
  call dplot_3D_intensity_animated("plot3Dintensity","X","Y","Z",X*FPI,Y*FPI,ZMATT)

  ! !2-dim + time animatino:
  call dplot_3D_surface_animated("plot3Dsurface","X","Y","Z",X*FPI,Y*FPI,ZMATT)

  ! !Vector Field:
  call dplot_vector_field("plotVF",X,Y,VX,VY,XLABEL="X",YLABEL="Y")

  ! !Vector Field + time animation:
  ! call dplot_vector_field("plotVF",X,Y,VTX,VTY,XLABEL="X",YLABEL="Y")


END program testDLPLOT

