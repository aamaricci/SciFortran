!+-----------------------------------------------------------------+
!PROGRAM  : PLOT_DISLIN3D
!TYPE     : Subroutine
!PURPOSE  : 
!+-----------------------------------------------------------------+
subroutine R_plot_dislin3D(pname,xname,yname,zname,X,Y,GF,theta,phi,rho,imovie)
  character(len=*)        :: pname
  character(len=*)        :: xname,yname,zname
  real(8),dimension(:)    :: X,Y
  real(8),dimension(:,:)  :: GF
  logical,optional        :: imovie
  real(8),optional        :: theta,phi,rho    !cf. WIKIPEDIA spherical coordinates
  real(4)                 :: theta_,phi_,rho_
  theta_ = 210.0;if(present(theta))theta_=real(theta,4)
  phi_   = 30.0;if(present(phi))phi_=real(phi,4)
  rho_   = 5.5 ;if(present(rho))rho_=real(rho,4)
  call system("rm -rf "//trim(pname))
  call system("mkdir "//trim(pname))
  write(*,"(A,A)")"print: ",adjustl(trim(pname))
  call plot3Dintensity("3dPlot_"//adjustl(trim(pname))//".png",xname,yname,zname,X,Y,GF)
  call plot3Dsurface("3dSurface_"//adjustl(trim(pname))//".png",xname,yname,zname,X,Y,GF,theta_,phi_,rho_)
  if(present(imovie))then!if imovie;Animated Gif
     write(*,"(A)")"get a movie too:        "
     call plot3Dsurface_rotating("3dSurface_"//adjustl(trim(pname)),xname,yname,zname,X,Y,GF)
  endif
  call system("mv *"//trim(pname)//".png *.gif "//trim(pname)//"/ 2>/dev/null")
end subroutine R_plot_dislin3D
!-------------------
!-------------------
!-------------------
subroutine C_plot_dislin3D(pname,xname,yname,zname,X,Y,GF,theta,phi,rho,imovie)
  character(len=*)         :: pname
  character(len=*)         :: xname,yname,zname
  real(8),dimension(:)     :: X,Y
  complex(8),dimension(:,:):: GF
  logical,optional         :: imovie
  real(8),optional        :: theta,phi,rho    !cf. WIKIPEDIA spherical coordinates
  real(4)                 :: theta_,phi_,rho_
  theta_ = 210.0;if(present(theta))theta_=real(theta,4)
  phi_   = 30.0;if(present(phi))phi_=real(phi,4)
  rho_   = 5.5 ;if(present(rho))rho_=real(rho,4)
  call system("rm -rf "//trim(pname))
  call system("mkdir "//trim(pname))
  write(*,"(A,A)")"print: ",adjustl(trim(pname))
  call plot3Dintensity("3dPlot_Im"//adjustl(trim(pname))//".png",xname,yname,zname,X,Y,aimag(GF))
  call plot3Dintensity("3dPlot_Re"//adjustl(trim(pname))//".png",xname,yname,zname,X,Y,real(GF,8))
  call plot3Dsurface("3dSurface_Im"//adjustl(trim(pname))//".png",xname,yname,zname,X,Y,aimag(GF),&
       theta_,phi_,rho_)
  call plot3Dsurface("3dSurface_Re"//adjustl(trim(pname))//".png",xname,yname,zname,X,Y,real(GF,8),&
       theta_,phi_,rho_)
  if(present(imovie))then
     write(*,"(A)")"get a movie too:       "
     call plot3Dsurface_rotating("3dSurface_Im"//adjustl(trim(pname)),xname,yname,zname,X,Y,aimag(GF))
     call plot3Dsurface_rotating("3dSurface_Re"//adjustl(trim(pname)),xname,yname,zname,X,Y,real(GF,8))
  endif
  call system("mv *"//trim(pname)//".png *.gif "//trim(pname)//"/ 2>/dev/null")
  call system("mv *"//trim(pname)//"*_rotating  "//trim(pname)//"/ 2>/dev/null")
end subroutine C_plot_dislin3D
!********************************************************************
!********************************************************************
!********************************************************************








!+-----------------------------------------------------------------+
!PROGRAM  : PLOT3DINTENSITY
!TYPE     : Subroutine
!PURPOSE  : 
!+-----------------------------------------------------------------+
subroutine plot3Dintensity(pname,xname,yname,zname,X,Y,GF)
  character(len=*)      :: pname,xname,yname,zname
  integer               :: Nx,Ny
  real(4)               :: xmin,xmax,dex
  real(4)               :: ymin,ymax,dey
  real(4)               :: zmin,zmax,dez
  real(8),dimension(:)  :: X,Y
  real(8),dimension(:,:):: GF

  !Get dimension of the function to plot:
  Nx=size(GF,1)  ; Ny=size(GF,2)

  !Get min&max values of the X,Y arrays and set the step:
  xmin=minval(X) ; xmax=maxval(X)
  ymin=minval(Y) ; ymax=maxval(Y)
  dex=abs(xmax-xmin)/5.
  dey=abs(ymax-ymin)/5.

  !Set the steps for the functions to plot:
  zmin=minval((GF));zmax=maxval((GF))
  GF=GF+abs(zmin) !shift to 0
  if(maxval(abs(GF))>1.d-12) GF=GF/maxval(abs(GF)) !normalization
  zmin=minval((GF));zmax=maxval((GF))
  if(zmin==zmax)zmax=zmin+0.01d0
  dez=max(abs(zmin),abs(zmax))/2.

  CALL METAFL('PNG')         !DEFINES THE METAFILE FORMAT (\ie THE OUTPUT TYPE)
  CALL WINSIZ(1920,1280)     !16:9, SET THE DIMENSION/RESOLUTION OF THE IMAGE
  CALL SCRMOD('REVERSE')     !SET BACKGROUND TO WHITE
  CALL IMGFMT('RGB')         !SET PALETTE TO R.G.B. def:8palette
  CALL SETFIL(trim(pname))   !SET THE FILE NAME
  CALL FILMOD('DELETE')      !ACTION OF A PRE-EXISENT FILE def:COUNT
  !CALL UNITS('POINTS')       !CHANGE UNITS TO POINTS. def:CM
  CALL SETPAG('DA4L')        !SET PAGE SIZE: DA4L/P DIN A4, USAL/P US LETTER 
  CALL ERRDEV('APPEND')     !SET BEHAVIOR ON ERROR MESSAGE
  CALL ERRFIL('dislin.out') !SET THE FILE NAME FOR ERROR AND OUTPUT MESSAGES

  !INIT DISLIN:
  CALL DISINI
  CALL ERRMOD('PROTOCOL','FILE')!SET MOD ON ERROR MESSAGE

  !SET SOME PAGE PROPERTIES: 
  CALL PAGORG('BOTTOM')     !SETS THE ORIGIN AT "BOTTOM"-LEFT
  !CALL PAGHDR("","",3,0)    !PLOTS A HEADER WITH INFO
  !CALL PAGERA               !PLOTS THE BORDER

  !SET SOME AXIS PROPERTIES:
  !CALL AXSPOS(2970,2100)     !SET THE AXIS-POSITION. don't need
  !CALL AXSLEN(2100,1500)     !SET THE LENGTH OF THE AXIS. weird units
  CALL TICKS(5,"XYZ")        !SET THE NUMBER OF TICKS per AXIS
  CALL TICLEN(45,15)         !LENGTH OF MAJ/MIN TICKS
  CALL LABDIG(3,'Z')         !SET NUMBER OF DIGITS: -2=AUTO, -1=INT, n=#DIGITS
  CALL LABDIG(0,'XY')        !SET NUMBER OF DIGITS: -2=AUTO, -1=INT, n=#DIGITS
  CALL FRAME(3)              !SET THE THICKNESS OF THE FRAME

  !COLOR:
  CALL SETVLT("RAIN")        !SET COLOR PALETTE: (R)RAIN,(R)TEMP,(R)GREY,VGA

  !CHARACTER:
  CALL HEIGHT(50)            !SET CHARACHTER HEIGHT
  CALL HWFONT                !SET TO HARDWARE/SYSTEM FONT
  CALL TEXMOD('ON')          !SET LaTeX MODE

  !TITLE:
  CALL TITLIN(trim(pname),1) !SET THE TITLE AND ITS LINE

  !AXIS:
  CALL NAME(trim(xname),'X')
  CALL NAME(trim(yname),'Y')
  CALL NAME(trim(zname)//" ${\small\rm (normalized)}$",'Z')

  !PLOT:
  CALL INTAX
  CALL AUTRES(Nx,Ny)
  CALL GRAF3(xmin,xmax,xmin,dex,ymin,ymax,ymin,dey,0.,1.01,0.,1.)
  CALL CRVMAT(real(GF,4),Nx,Ny,20,20)

  !DISLIN FINALIZE:
  CALL TITLE
  CALL DISFIN
end subroutine plot3Dintensity
!********************************************************************
!********************************************************************
!********************************************************************





!+-----------------------------------------------------------------+
!PROGRAM  : PLOT3DSURFACE
!TYPE     : Subroutine
!PURPOSE  : 
!+-----------------------------------------------------------------+
subroutine plot3Dsurface(pname,xname,yname,zname,X,Y,GF,theta,phi,rho)
  character(len=*)      :: pname,xname,yname,zname
  integer               :: Nx,Ny
  real(4)               :: xmin,xmax,dex
  real(4)               :: ymin,ymax,dey
  real(4)               :: zmin,zmax,dez
  real(4)               :: theta,phi,rho
  real(8),dimension(:)  :: X,Y
  real(8),dimension(:,:):: GF
  integer               :: sg1,sg2

  Nx=size(GF,1) ; Ny=size(GF,2)
  xmin=minval(X);xmax=maxval(X)
  if(xmin==xmax)then
     xmin=xmin-0.25
     xmax=xmax+0.25
  endif
  dex=abs(xmax-xmin)/5.

  ymin=minval(Y);ymax=maxval(Y)
  if(ymin==ymax)then
     ymin=ymin-0.25
     ymax=ymax+0.25
  endif
  dey=abs(ymax-ymin)/5.

  zmin=minval((GF));zmax=maxval((GF))
  if(zmin==zmax)zmax=zmin+0.01d0
  dez=max(abs(zmin),abs(zmax))    
  !dez=abs(zmax-zmin)
  dez=dez/5.

  CALL METAFL('PNG')
  CALL SETFIL(trim(pname))
  CALL SCRMOD('REVERSE')
  CALL WINSIZ(1920,1280)

  CALL ERRDEV('APPEND')
  CALL ERRFIL('dislin.out')
  CALL DISINI
  CALL ERRMOD('PROTOCOL','FILE')

  CALL PAGERA
  CALL HWFONT  
  CALL HEIGHT(50)            !SET CHARACHTER HEIGHT
  CALL TITLIN(trim(pname),2)

  CALL TEXMOD('ON')
  CALL NAME(trim(xname),'X')
  CALL NAME(trim(yname),'Y')
  CALL NAME(trim(zname),'Z')
  CALL LABDIG(6,'Z')
  CALL LABDIG(0,'XY')

  !CALL VIEW3D(sg1*3.75,sg2*3.75,2.,'ABS')
  CALL VIEW3D(theta,phi,rho,'ANGLE')
  CALL SETVLT('RAIN') !TEMP
  CALL GRAF3D(xmin,xmax,xmin,dex,ymin,ymax,ymin,dey,zmin,zmax,zmin,dez)
  CALL SHDMOD("SMOOTH","SURFACE")
  CALL SURSHD(real(X,4),Nx,real(Y,4),Ny,real(GF,4))
  CALL SURMAT(real(GF,4),Nx,Ny,1,1)
  CALL TITLE
  CALL DISFIN
  return
end subroutine plot3Dsurface
!********************************************************************
!********************************************************************
!********************************************************************




!+-----------------------------------------------------------------+
!PROGRAM  : PLOT3DSURFACE_rotating
!TYPE     : Subroutine
!PURPOSE  : realizes an animated .gif out of the .png surface plots
!+-----------------------------------------------------------------+
subroutine plot3Dsurface_rotating(pname,xname,yname,zname,X,Y,GF)
  character(len=*) :: xname,yname,zname,pname
  character(len=7) :: char_temp
  character(len=32):: xsuff,ysuff
  character(len=128):: cmd
  integer,parameter :: Nangle=30
  integer :: Nx,Ny,ixstep,iystep,i
  real(4) :: xstep,ystep,angle
  real(4) :: xmin,xmax,dex
  real(4) :: ymin,ymax,dey
  real(4) :: zmin,zmax,dez
  real(8),dimension(:) :: X,Y
  real(8),dimension(:,:) :: GF

  Nx=size(GF,1) ;Ny=size(GF,2)
  xmin=minval(X);xmax=maxval(X)
  if(xmin==xmax)then
     xmin=xmin-0.25
     xmax=xmax+0.25
  endif
  dex=abs(xmax-xmin)/5.

  ymin=minval(Y);ymax=maxval(Y)
  if(ymin==ymax)then
     ymin=ymin-0.25
     ymax=ymax+0.25
  endif
  dey=abs(ymax-ymin)/5.

  zmin=minval((GF));zmax=maxval((GF))
  if(zmin==zmax)then
     zmin=zmin-0.25
     zmax=zmax+0.25
  endif
  dez=max(abs(zmin),abs(zmax))    
  !dez=abs(zmax-zmin)
  dez=dez/5.

  do i=1,Nangle
     angle=360./dble(Nangle)*dble(i)
     write(char_temp,'(f4.0)')angle
     !print*,char_temp
     CALL METAFL('PNG')
     CALL SETFIL(trim(adjustl(trim(char_temp)//trim(pname)//".png")))
     CALL SCRMOD('REVERSE')
     CALL WINSIZ(800,600)

     CALL ERRDEV('APPEND')
     CALL ERRFIL('dislin.out')
     CALL DISINI
     CALL ERRMOD('PROTOCOL','FILE')
     CALL HEIGHT(50)            !SET CHARACHTER HEIGHT
     CALL TEXMOD('ON')

     CALL PAGERA
     CALL HWFONT  
     CALL TITLIN(trim(pname)//" angle:"//trim(adjustl(char_temp)),2)

     CALL NAME(trim(xname),'X')
     CALL NAME(trim(yname),'Y')
     CALL NAME(trim(zname),'Z')
     CALL LABDIG(6,'Z')
     CALL LABDIG(0,'XY')


     CALL VIEW3D(angle,30.,5.5,'ANGLE')
     CALL SETVLT('RAIN') !TEMP
     CALL GRAF3D(xmin,xmax,xmin,dex,  &
          ymin,ymax,ymin,dey,  &
          zmin,zmax,zmin,dez)
     CALL SHDMOD("SMOOTH","SURFACE")
     CALL SURSHD(real(X,4),Nx,real(Y,4),Ny,real(GF,4))
     CALL SURMAT(real(GF,4),Nx,Ny,1,1)

     CALL TITLE
     CALL DISFIN
  enddo

  call system("if [ ! -d "//trim(pname)//"_rotating ]; then mkdir "//trim(pname)//"_rotating ; fi")
  call system("mv *."//trim(pname)//".png "//trim(pname)//"_rotating/ 2>/dev/null")    
  call system("echo '#!/bin/bash' > make_gif.sh")
  call system("echo 'echo 'Converting to .gif'' >> make_gif.sh")
  call system("echo 'convert -delay 50 `ls *.*.png|sort -n` -loop 0 "//trim(pname)//".gif' >> make_gif.sh")
  call system("echo 'if [ ! -d "//trim(pname)//".gif ]; then rm -f *.*.png ; fi' >> make_gif.sh")
  call system("chmod u+x make_gif.sh")
  call system("mv make_gif.sh "//trim(pname)//"_rotating/ 2>/dev/null")

end subroutine plot3Dsurface_rotating
!********************************************************************
!********************************************************************
!********************************************************************
